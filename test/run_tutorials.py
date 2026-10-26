#!/usr/bin/env python3
"""Execute documented Julia notebooks without Jupyter or extra Python packages."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time


SMOKE = ("examples/lin_elast_2D_static.ipynb", "examples/plane-beam-structure.ipynb")
EXCLUDED = {
    "examples/lin_elast_2D_transient.ipynb": "Transient analysis is outside this test suite's scope.",
    "examples/space-time-FEM.ipynb": "Space-time dynamics is outside this test suite's scope.",
}


def discover(root):
    pattern = re.compile(
        r"https://github\.com/perebalazs/LowLevelFEM\.jl/blob/main/(examples/[^)\s]+\.ipynb)"
    )
    paths = set()
    for page in sorted((root / "docs/src/tutorials").glob("*.md")):
        paths.update(pattern.findall(page.read_text(encoding="utf-8")))
    for name in paths:
        path = (root / name).resolve()
        if not path.is_relative_to(root / "examples") or not path.is_file():
            raise ValueError(f"Invalid documented notebook: {name}")
    return sorted(paths)


def select(catalog, requested, all_notebooks):
    if all_notebooks and requested:
        raise ValueError("Use notebook names or --all, not both.")
    selected = [p for p in catalog if p not in EXCLUDED] if all_notebooks else list(SMOKE)
    if requested:
        selected = []
        for name in requested:
            matches = [p for p in catalog if p == name or Path(p).name == name]
            if len(matches) != 1:
                raise ValueError(f"Expected one documented notebook for {name!r}; use --list.")
            if matches[0] in EXCLUDED:
                raise ValueError(f"{matches[0]}: {EXCLUDED[matches[0]]}")
            if matches[0] not in selected:
                selected.append(matches[0])
    if any(p not in catalog for p in selected):
        raise ValueError("A default smoke notebook is no longer linked by the documentation.")
    return selected


def prepare(root, notebook, work):
    source = root / notebook
    data = json.loads(source.read_text(encoding="utf-8"))
    if data.get("nbformat") != 4:
        raise ValueError("Only nbformat 4 notebooks are supported.")
    language = data.get("metadata", {}).get("kernelspec", {}).get("language", "julia")
    if language.lower() != "julia":
        raise ValueError(f"Expected a Julia notebook, got {language!r}.")
    # Copy neighboring input files, so geometry paths resolve and notebook
    # outputs cannot overwrite the examples. Notebook contents stay unchanged.
    work.mkdir()
    for asset in source.parent.iterdir():
        if asset.is_file() and not asset.name.startswith(".") and asset.suffix != ".ipynb":
            shutil.copy2(asset, work / asset.name)
    cells = []
    for index, cell in enumerate(data["cells"], 1):
        if cell.get("cell_type") != "code":
            continue
        text = cell.get("source", "")
        text = "".join(text) if isinstance(text, list) else text
        if not text.strip():
            continue
        # Preserve original notebook cell numbers in filenames and error logs.
        name = f"notebook_cell_{index:04d}.jl"
        (work / name).write_text(text, encoding="utf-8")
        cells.append(name)
    if not cells:
        raise ValueError("Notebook has no nonempty code cells.")
    return cells


def stop_process(proc):
    if proc.poll() is not None:
        return
    os.killpg(proc.pid, signal.SIGTERM)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()


def execute(root, scripts, notebook, output, julia, threads, timeout):
    started = time.monotonic()
    work = output / Path(notebook).stem
    log = output / f"{Path(notebook).stem}.log"
    result = {"notebook": notebook, "status": "failed", "log": str(log), "workdir": str(work)}
    try:
        cells = prepare(root, notebook, work)
        result["code_cells"] = len(cells)
        command = [julia, "--startup-file=no", "--history-file=no", f"--project={root}",
                   f"--threads={threads}", str(scripts / "tutorial_worker.jl"),
                   str(work), notebook, str(scripts / "tutorial_checks.jl")]
        env = dict(os.environ, GKSwstype="100")
        with log.open("w", encoding="utf-8") as stream:
            proc = subprocess.Popen(command, cwd=work, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, start_new_session=True)
            try:
                result["exit_code"] = proc.wait(timeout=timeout)
                result["status"] = "passed" if proc.returncode == 0 else "failed"
            except subprocess.TimeoutExpired:
                stop_process(proc)
                result.update(status="timeout", error=f"Exceeded {timeout:g} seconds.")
            except KeyboardInterrupt:
                stop_process(proc)
                raise
        result["headless_calls"] = [line.strip() for line in log.read_text(encoding="utf-8").splitlines()
                                    if line.startswith("HEADLESS_SKIP ")]
    except (OSError, ValueError, KeyError) as error:
        result["error"] = str(error)
    result["seconds"] = round(time.monotonic() - started, 2)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebooks", nargs="*", help="Documented paths or unique filenames; default: two smoke notebooks.")
    parser.add_argument("--all", action="store_true", help="Run every documented, non-transient tutorial (can be expensive).")
    parser.add_argument("--list", action="store_true", help="List scope without executing notebooks.")
    parser.add_argument("--julia", default="julia", help="Existing Julia executable; no packages are installed.")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=600, help="Seconds per notebook, including compilation.")
    parser.add_argument("--output-dir", type=Path, help="New report directory; default: a fresh directory under /tmp.")
    args = parser.parse_args()
    if args.threads < 1 or args.timeout <= 0:
        parser.error("--threads and --timeout must be positive.")
    scripts = Path(__file__).resolve().parent
    root = scripts.parent
    try:
        catalog = discover(root)
        selected = select(catalog, args.notebooks, args.all)
    except ValueError as error:
        parser.error(str(error))
    if args.list:
        for name in catalog:
            status = f"excluded: {EXCLUDED[name]}" if name in EXCLUDED else (
                "default smoke + numerical checks" if name in SMOKE else "available; not in default run")
            print(f"{name}: {status}")
        return 0
    julia = shutil.which(args.julia)
    if julia is None:
        parser.error(f"Julia executable not found: {args.julia}")
    if args.output_dir:
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=False)
    else:
        output = Path(tempfile.mkdtemp(prefix="lowlevelfem-tutorials-"))
    print(f"Reports and isolated notebook workdirs: {output}", flush=True)
    results = []
    for name in selected:
        print(f"RUN {name}", flush=True)
        result = execute(root, scripts, name, output, julia, args.threads, args.timeout)
        results.append(result)
        print(f"{result['status'].upper()} {name} ({result['seconds']:.1f}s): {result['log']}", flush=True)
        if result["status"] != "passed":
            if "error" in result:
                print(result["error"], flush=True)
            elif Path(result["log"]).exists():
                print("\n".join(Path(result["log"]).read_text(encoding="utf-8").splitlines()[-35:]), flush=True)
    report = {"root": str(root), "julia": julia, "threads": args.threads,
              "results": results, "excluded": EXCLUDED,
              "not_selected": [p for p in catalog if p not in selected and p not in EXCLUDED]}
    report_path = output / "report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    passed = sum(r["status"] == "passed" for r in results)
    print(f"Tutorial notebooks: {passed}/{len(results)} passed. Report: {report_path}", flush=True)
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

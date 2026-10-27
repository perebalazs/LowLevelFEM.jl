#!/usr/bin/env python3
"""Run explicitly marked, self-contained documentation examples without installing packages."""
import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile

def examples(folder):
    for page in sorted(folder.rglob('*.md')):
        for code in re.findall(r'^```julia\s*\n(.*?)^```\s*$',page.read_text(),re.M|re.S):
            first = code.splitlines()[0] if code.splitlines() else ''
            marker = re.fullmatch(r'# docs-check: ([\w-]+)(; optional)?',first)
            if marker:
                yield page,marker[1],bool(marker[2]),code

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--julia',default='julia')
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--timeout',type=float,default=180)
    parser.add_argument('--with-optional',action='store_true',help='also run examples needing already installed optional packages')
    args = parser.parse_args()
    if args.threads<1 or args.timeout<=0:
        parser.error('threads and timeout must be positive')
    root=Path(__file__).resolve().parents[1]
    count=failed=skipped=0
    with tempfile.TemporaryDirectory(prefix='llfem-doc-examples-') as folder:
        for page,name,optional,code in examples(root/'docs/src'):
            if optional and not args.with_optional:
                print('SKIP',name,'(optional packages; use --with-optional)',flush=True)
                skipped+=1
                continue
            script=Path(folder)/(name+'.jl')
            script.write_text(code)
            count+=1
            try:
                result=subprocess.run([args.julia,'--startup-file=no','--compiled-modules=existing',
                    '--project='+str(root),'--threads='+str(args.threads),str(script)],
                    cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                    timeout=args.timeout,check=False)
                ok=result.returncode==0
                if not ok: print(result.stdout,flush=True)
            except subprocess.TimeoutExpired:
                print('TIMEOUT',name,flush=True)
                ok=False
            print(('PASS' if ok else 'FAIL'),name,'—',page.relative_to(root),flush=True)
            failed+=not ok
    if not count:
        print('ERROR: no examples executed',file=sys.stderr)
        return 1
    print(f'Documentation examples: {count-failed}/{count} passed; {skipped} optional skipped.')
    return int(bool(failed))

if __name__=='__main__':
    sys.exit(main())

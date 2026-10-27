#!/usr/bin/env python3
"""Check loaded public Julia bindings using the existing project environment."""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

KNOWN_UNDEFINED = {'getEigenVectors', 'nodalForce!', 'nodalAcceleration!', 'probe_field'}
JULIA_CHECK = r'''
using LowLevelFEM
for name in sort(names(LowLevelFEM); by=string)
    name === :LowLevelFEM && continue
    status = !isdefined(LowLevelFEM,name) ? "UNDEFINED" :
             Base.Docs.hasdoc(LowLevelFEM,name) ? "DOCUMENTED" : "UNDOCUMENTED"
    println("LLFEM_DOC\t",status,"\t",name)
end
'''

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict', action='store_true', help='also fail for known undefined exports')
    parser.add_argument('--julia', default='julia', help='Julia executable (default: julia)')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix='llfem-doc-audit-') as folder:
        script = Path(folder)/'audit.jl'
        script.write_text(JULIA_CHECK)
        result = subprocess.run([args.julia, '--startup-file=no', '--compiled-modules=existing',
                                 '--project='+str(root), str(script)], cwd=root,
                                text=True, stdout=subprocess.PIPE, check=False)
    if result.returncode:
        return result.returncode
    counts = {'DOCUMENTED':0,'UNDOCUMENTED':0,'UNDEFINED':0}
    failed = False
    for line in result.stdout.splitlines():
        if not line.startswith('LLFEM_DOC\t'):
            print(line)
            continue
        _,status,name = line.split('\t')
        counts[status] += 1
        if status == 'UNDOCUMENTED':
            print('ERROR: missing docstring:',name)
            failed = True
        elif status == 'UNDEFINED':
            known = name in KNOWN_UNDEFINED
            print(('WARNING' if known and not args.strict else 'ERROR')+': exported but undefined:',name)
            failed |= args.strict or not known
    print('Public documentation: {DOCUMENTED} documented, {UNDOCUMENTED} missing, '
          '{UNDEFINED} undefined exports.'.format(**counts))
    if not sum(counts.values()):
        print('ERROR: no public API records received', file=sys.stderr)
        return 1
    return int(failed)

if __name__ == '__main__':
    sys.exit(main())

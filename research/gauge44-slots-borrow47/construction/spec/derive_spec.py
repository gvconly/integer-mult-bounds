"""Re-derive data/SPEC.json: discovery -> unpin -> freeze (the research/gmx-t74-borrow43-handoff procedure), at T = 8.

Base SPEC = research/gauge44-slots-arcflow's data/SPEC.json plus the T = 8 geometry of arcflow's parent558 (the literals
of arcflow's code/assemble.py and invoice: parent live/physical/idle 11704/12218/16128, 504 parent helpers, payload 512),
the outer refinement (endpoint 24, grid 18, margin 28 digits, 12 bootstrap steps) and this package's PLAN/OWNERS hashes.
A copy of the package runs verify_global.py --discover; new/unpin_spec.py (main) removes the arcflow values that the
independent recomputation disagrees with; freeze_spec.py writes the recomputed ones. Not run by verify_borrow47.py (the
strict replay checks every frozen value); usage: derive_spec.py MAIN LOCAL_DIR WORK OUT_SPEC.json
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import hashlib, json, shutil, subprocess, sys
from pathlib import Path
main, local, work, out = (Path(x).resolve() for x in sys.argv[1:5])
P = Path(__file__).resolve().parents[2]
assert not work.exists(); work.mkdir(parents=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
pkg = work / 'pkg'; shutil.copytree(P, pkg, ignore=shutil.ignore_patterns('__pycache__'))
spec = json.loads((main / 'research/gauge44-slots-arcflow/data/SPEC.json').read_text())
spec['geometry'] = dict(charged_prefix_bits=512, dimension=73, pairs=4, parent_helpers=504, parent_idle=16128, parent_live=11704,
                        parent_physical=12218, replicas=8, stages=9)
spec['outer_refinement'] = dict(grid_digits=18, margin_digits=28, bootstrap_steps=12, endpoint_digits=24)
spec['data'] = dict(owners_sha256=sha(pkg / 'data/OWNERS.json.gz'), plan_sha256=sha(pkg / 'data/PLAN.json'))
(pkg / 'data/SPEC.json').write_text(json.dumps(spec, indent=2) + '\n')
run = lambda *a: subprocess.run([sys.executable, '-B', *map(str, a)], check=True)
run(pkg / 'verify_global.py', '--main', main, '--local', local, '--work', work / 'discovery', '--discover')
run(main / 'research/gmx-t74-borrow43-handoff/new/unpin_spec.py', pkg / 'data/SPEC.json', work / 'discovery')
run(pkg / 'freeze_spec.py', '--spec', pkg / 'data/SPEC.json', '--discovery', work / 'discovery')
shutil.copy2(pkg / 'data/SPEC.json', out); print('SPEC', out, sha(out))

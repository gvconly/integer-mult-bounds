"""Derive data/PLAN.json from research/gauge44-slots-arcflow's PLAN, the native borrowed provider and the owner map.

Changed fields only: the gauge44 provider's widths (|final|-|start| of the native canonical word), its actual paid profile,
artifact and semantic hashes and orientation program hashes (exactly the provider-derived fields that
research/gmx-t74-borrow43-handoff/new/rebind.py rebinds), every call's local_program_sha256 and per-call widths, the global
profile (72 x local paid profile + the idle remainders of OWNERS) and the nomination metadata. Schedule, windows, dispatch,
ports, centres, sides and pairs are arcflow's, unchanged.
usage: make_plan.py MAIN LOCAL_DIR OWNERS.json.gz OUT_PLAN.json
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import collections, gzip, hashlib, json, sys
from pathlib import Path

main, local, owners, out = (Path(x).resolve() for x in sys.argv[1:5])
def rd(p): return json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
arc_path = main / 'research/gauge44-slots-arcflow/data/PLAN.json'
plan = rd(arc_path)
sys.path.insert(0, str(main / 'research/gmx-t74-borrow43-handoff/code'))
from generated_identity import artifact_identities
desc = plan['provider_registry']['gauge44']; word = rd(local / 'CANONICAL-WORD.json.gz'); checks = rd(local / 'LOCAL-CHECKS.json')
assert word['ports'] == desc['ports'] and word['centers'] == desc['centers'] and (word['v'], word['h']) == (112, 9)
widths = [len(word['frames'][word['final'][r]]) - len(word['frames'][word['start'][r]]) for r in range(224, 224 + word['R'])]
desc.update(widths=widths,
            actual_paid_profile={str(k): v for k, v in sorted({int(k): v for k, v in checks['forward']['paid_profile'].items()}.items())},
            artifact_sha256={n: sha(local / n) for n in desc['artifact_sha256']},
            orientation_program_sha256={o: checks[o]['sha256'] for o in ('backward', 'forward')},
            artifact_semantic_sha256=artifact_identities(local))
for c in plan['per_call_widths']: c['widths'] = widths
for c in plan['provider_calls']: c['local_program_sha256'] = checks[c['orientation']]['sha256']
roles = rd(owners); used = collections.Counter()
for z in roles: used[z['owner_replica'], z['child'], z['port']] += z['width']
idle = collections.Counter()
for key, (cap, times) in plan['windows'].items():
    o, ch = map(int, key.split(','))
    for port in range(112):
        if cap - used[o, ch, port]: idle[cap - used[o, ch, port]] += 1
prof = collections.Counter()
for r, n in desc['actual_paid_profile'].items(): prof[int(r)] += len(plan['provider_calls']) * n
prof.update(idle)
assert plan['M'] * plan['live'] - sum(r * n for r, n in prof.items()) == plan['deficit'] == 2560, 'rank-mass identity'
plan['profile'] = {str(k): v for k, v in sorted(prof.items())}
plan.update(status='FROZEN_EXPLICIT_GAUGE44_WHOLE_OWNER_NOMINATION_BORROW47',
            solver_status='window-LP support MILP (HiGHS, optimal on the support and on the remainder<=3 neighbourhood), '
                          'exact per-pattern time-incidence MILP, Konig colouring (construction/owners/t8pack.py)',
            gap=None, source_plan_sha256=sha(arc_path), assignment_sha256=sha(owners), source_owners_sha256=sha(owners))
out.write_text(json.dumps(plan, indent=2) + '\n')
print('PLAN', out, 'R', len(widths), 'idle', dict(sorted(idle.items())), 'profile', plan['profile'])

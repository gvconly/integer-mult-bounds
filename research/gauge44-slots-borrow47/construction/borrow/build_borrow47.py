"""Build and audit the 47-alias BORROW-WORD on the selected slot word.

Runs the unchanged BORROW-WORD builder, auditor, negative controls and guards of research/gmx-t74-borrow43-handoff
(new/borrow/) from the pinned main checkout:
  build_word.py   selected word + block certificate + WITNESS-f47.json -> BORROW-WORD (content pinned in BORROW-PINS.json)
  audit_local.py  --base selected word: PASS_INDEPENDENT_EXACT_BORROWING_WORD
  borrow_controls.py   every mutation rejected for its expected reason
  guard_cut_centres / guard_helper_indices (PR2478 / PR2508 guards) with their negative controls
and checks that the word's residual widths are the per-helper widths of data/PLAN.json.
usage: build_borrow47.py MAIN WORK   (WORK must not exist); writes WORK/BORROW-WORD-f47.json.gz and WORK/BUILD.json
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import contextlib, gzip, hashlib, io, json, subprocess, sys
from pathlib import Path

main, work = (Path(x).resolve() for x in sys.argv[1:3])
assert __debug__ and not work.exists()
work.mkdir(parents=True)
HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]
B = main / 'research/gmx-t74-borrow43-handoff/new/borrow'
SRC = main / 'research/anticipatory-dirty-borrowing/inputs'
pins = json.loads((HERE / 'BORROW-PINS.json').read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for rel, h in pins['main_sources'].items(): assert sha(main / rel) == h, 'source pin: ' + rel
assert sha(HERE / 'WITNESS-f47.json') == pins['witness_sha256']
selected = SRC / 'selected-complex-word.json.gz'
# The selected word is the canonical word bound by research/gauge44-slots-arcflow's PLAN (the #1094/#1175 word).
arc = json.loads((main / 'research/gauge44-slots-arcflow/data/PLAN.json').read_text())
assert sha(selected) == arc['provider_registry']['gauge44']['artifact_sha256']['CANONICAL-WORD.json.gz']
word_path = work / 'BORROW-WORD-f47.json.gz'
r = subprocess.run([sys.executable, '-B', str(B / 'build_word.py'), str(selected), str(SRC / 'selected-block-certificate.json'),
                    str(HERE / 'WITNESS-f47.json'), str(word_path), 'selected1094'], cwd=B, capture_output=True, text=True)
assert r.returncode == 0, r.stdout + r.stderr
content = hashlib.sha256(gzip.decompress(word_path.read_bytes())).hexdigest()
assert content == pins['word_content_sha256'], ('word content', content)
r = subprocess.run([sys.executable, '-B', str(B / 'audit_local.py'), '--word', str(word_path), '--base', str(selected),
                    '--output', str(work / 'AUDIT-BORROW-WORD.json')], capture_output=True, text=True)
assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
audit = json.loads((work / 'AUDIT-BORROW-WORD.json').read_text())
assert audit['status'] == 'PASS_INDEPENDENT_EXACT_BORROWING_WORD' and audit['content_sha256'] == content
assert (audit['removed_helpers'], audit['warps'], audit['zero_buffers']) == (47, 47, 73)
r = subprocess.run([sys.executable, '-B', str(B / 'borrow_controls.py'), str(word_path), str(selected), str(work / 'BORROW-CONTROLS.json')],
                   capture_output=True, text=True)
assert r.returncode == 0 and 'ALL CONTROLS OK' in r.stdout, r.stdout[-3000:] + r.stderr[-3000:]
controls = json.loads((work / 'BORROW-CONTROLS.json').read_text())
assert controls['all_controls_rejected_for_expected_reason'] is True
sys.path.insert(0, str(B))
import guard_cut_centres as gcc, guard_helper_indices as ghi
w = json.loads(gzip.decompress(word_path.read_bytes()))
bad_cut = json.loads(json.dumps(w))
bad_cut['cut'][0] = ['add', bad_cut['cut'][0][1], next(x for x in range(224, 704) if x not in set(bad_cut['centers'])), 1, 1] + list(bad_cut['cut'][0][5:])
bad_idx = json.loads(json.dumps(w)); bad_idx['helper_final']['-1'] = list(next(iter(bad_idx['helper_final'].values())))
guards = dict(cut_centre_violations=len(gcc.violations(w)), cut_centre_negative_control_rejected=len(gcc.violations(bad_cut)) > 0,
              helper_index_violations=len(ghi.violations(w)), helper_index_negative_control_rejected=len(ghi.violations(bad_idx)) > 0)
assert guards == dict(cut_centre_violations=0, cut_centre_negative_control_rejected=True, helper_index_violations=0, helper_index_negative_control_rejected=True), guards
# residual widths of the remaining helpers, in order, are the provider widths of data/PLAN.json
rem = [x for x in range(224, 704) if x not in set(w['removed'])]
se = {int(k): v for k, v in w['start_extra'].items()}; hf = {int(k): v for k, v in w['helper_final'].items()}
widths = [len(hf.get(x, range(9))) - len(se.get(x, [])) for x in rem]
plan = json.loads((PKG / 'data/PLAN.json').read_text())
assert widths == plan['provider_registry']['gauge44']['widths']
receipt = dict(status='PASS_BORROW47_WORD_BUILT_AUDITED_AND_CONTROLLED', word_content_sha256=content, removed=w['removed'],
               local_profile=audit['local_profile'], local_mass=audit['local_mass'], width_inventory=audit['width_inventory'],
               zero_buffers=audit['zero_buffers'], warps=audit['warps'], borrow_controls=len(controls['controls']), guards=guards,
               widths_bound_to_plan=True)
(work / 'BUILD.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(receipt['status'], content)

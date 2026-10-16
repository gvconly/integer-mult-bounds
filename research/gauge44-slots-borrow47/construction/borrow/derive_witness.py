"""Derive the 47-alias borrowing witness from PR2226's 80-alias witness on the selected slot word.

The subset is fixed by a rule, not by a list: keep every PR2226 alias that is not a low alias and whose
lender is a helper register (borrow_from >= 2v). That drops the 32 low aliases (cold lenders) and the one
source-lender alias (register 653, lender 73). PR2226's own replay.py and operators.py are then run on the
subset; the only change to replay.py is the removal of its whole-witness count assertion
(len(pairs)==80 and 32 low), applied as an exact one-line text substitution of the pinned file.
Each kept pair gains retirement_basis = the literal frame at its retirement_frame in the selected word.
usage: derive_witness.py MAIN WORK OUT.json   (WORK must not exist)
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import gzip, json, subprocess, sys
from pathlib import Path

main, work, out = (Path(x).resolve() for x in sys.argv[1:4])
assert __debug__ and not work.exists()
work.mkdir(parents=True)
src = main / 'research/anticipatory-dirty-borrowing'
word_path = src / 'inputs/selected-complex-word.json.gz'
blocks_path = src / 'inputs/selected-block-certificate.json'
witness = json.loads((src / 'WITNESS.json').read_text())
word = json.loads(gzip.decompress(word_path.read_bytes()))
v = word['v']
pairs = witness['pairs']
assert len(pairs) == 80 and sum(p['low'] for p in pairs) == 32
keep = [i for i, p in enumerate(pairs) if not p['low'] and p['borrow_from'] >= 2 * v]
dropped_nonlow = [pairs[i]['register'] for i in range(80) if not pairs[i]['low'] and i not in keep]
assert len(keep) == 47 and dropped_nonlow == [653], (len(keep), dropped_nonlow)
idx = {old + 1: new + 1 for new, old in enumerate(keep)}   # operator numbering of the subset
corrections = [dict(c, operator_word=[(1 if o > 0 else -1) * idx[abs(o)] for o in c['operator_word']])
               for c in witness['initial_corrections'] if all(abs(o) in idx for o in c['operator_word'])]
inventory = dict(pairs=[pairs[i] for i in keep], early_tail_events=witness['early_tail_events'], initial_corrections=corrections)
(work / 'inventory.json').write_text(json.dumps(inventory))
line = "    assert len(pairs)==80 and sum(x['low'] for x in pairs)==32\n"
text = (src / 'replay.py').read_text()
assert text.count(line) == 1
(work / 'replay_subset.py').write_text(text.replace(line, '    pass\n'))
env_run = lambda args: subprocess.run([sys.executable, '-B', *map(str, args)], check=True, cwd=work)
env_run([work / 'replay_subset.py', '--word', word_path, '--blocks', blocks_path, '--inventory', work / 'inventory.json',
         '--recover-cold', '--output', work / 'replay.json'])
env_run([src / 'operators.py', '--replay', work / 'replay.json', '--output', work / 'operators.json'])
replay = json.loads((work / 'replay.json').read_text())
gaussian = json.loads((work / 'operators.json').read_text())
assert replay['status'] == 'PASS_EXACT_FREE_OPERATOR_BORROWING_REPLAY' and gaussian['status'] == 'PASS_EXACT_GAUSSIAN_BORROW_HANDOFFS'
assert replay['removed_helper_registers'] == 47 and gaussian['operators'] == 47
out_pairs = [dict(p, retirement_basis=word['frames'][p['retirement_frame']]) for p in inventory['pairs']]
result = dict(credit='#2226 anticipatory dirty borrowing (Cinna Davis / 1337isnot1337), the 47 full aliases on the selected slot word '
                     '(subset chosen by G. V. Chudnovsky with Anthropic Claude assistance)',
              variant='full47', pairs=out_pairs, initial_corrections=corrections, early_tail_events=inventory['early_tail_events'],
              replay=replay, gaussian=gaussian)
out.write_text(json.dumps(result))
print('PASS_DERIVED_FULL47_WITNESS', len(out_pairs), 'aliases; local profile', replay['local_profile'], 'mass', replay['local_mass'])

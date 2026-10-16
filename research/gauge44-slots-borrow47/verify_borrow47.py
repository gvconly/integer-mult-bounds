"""Complete verification of research/gauge44-slots-borrow47 from pinned main sources.

Stages (each a subprocess with assertions enabled; logs and STAGES.json in --work):
  manifest        exact package file inventory and hashes
  witness         construction/borrow/derive_witness.py: PR2226's replay/operators on the rule-selected 47 aliases;
                  must equal the shipped WITNESS-f47.json except timing fields
  word            construction/borrow/build_borrow47.py: build, audit (--base), negative controls, guards, width binding
  native          research/gmx-t74-borrow43-handoff/new/run_native_borrow.py on the built word (or --native-from reuse,
                  accepted only if its FULL-WORD is the built word and all five steps passed)
  plan            construction/owners/make_plan.py regenerates data/PLAN.json from arcflow's PLAN, the native provider
                  and data/OWNERS.json.gz (compared on every field; receipt hashes compared for the three byte-stable artifacts)
  recount         construction/owners/recount_owners.py (standard library) on data/PLAN.json and data/OWNERS.json.gz
  strict          verify_global.py: bind, emit, audit with mutation controls, invoice, assembly, admission and bit controls
  independent     research/gauge44-slots-arcflow/independent/price_profile.py on data/PLAN.json
usage: python -B research/gauge44-slots-borrow47/verify_borrow47.py --main <repo root> --work <new absolute dir> [--native-from DIR]
Prepared by G. V. Chudnovsky with Anthropic Claude assistance. Apache-2.0."""
import argparse, gzip, hashlib, json, os, subprocess, sys, time
from pathlib import Path

assert __debug__, 'assertions required'
P = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--main', type=Path, required=True); ap.add_argument('--work', type=Path, required=True)
ap.add_argument('--native-from', type=Path, help='completed run_native_borrow.py directory of the same word (optional reuse)')
a = ap.parse_args(); main = a.main.resolve(); work = a.work.resolve()
assert not work.exists() and not work.is_relative_to(P); work.mkdir(parents=True)
def rd(p): return json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1', PYTHONINTMAXSTRDIGITS='0')
for k in ('PYTHONOPTIMIZE', 'WITNESS_DISCOVER', 'WITNESS_RECORD', 'WITNESS_SPEC'): env.pop(k, None)
stages = []
def run(name, script, *args, cwd=None):
    t = time.monotonic(); log = work / (name + '.log')
    with log.open('w') as f:
        r = subprocess.run([sys.executable, '-B', str(script), *map(str, args)], env=env, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    stages.append(dict(name=name, exit_code=r.returncode, seconds=round(time.monotonic() - t, 3), log_sha256=sha(log)))
    (work / 'STAGES.json').write_text(json.dumps(stages, indent=2) + '\n')
    assert r.returncode == 0, (name, log.read_text()[-6000:])
    print('PASS', name, stages[-1]['seconds'], flush=True)
run('manifest', P / 'verify_manifest.py')
run('witness', P / 'construction/borrow/derive_witness.py', main, work / 'witness-work', work / 'WITNESS-f47.json')
def untimed(x):
    if isinstance(x, dict): return {k: untimed(v) for k, v in x.items() if k != 'seconds'}
    if isinstance(x, list): return [untimed(v) for v in x]
    return x
assert untimed(rd(work / 'WITNESS-f47.json')) == untimed(rd(P / 'construction/borrow/WITNESS-f47.json')), 'derived witness differs'
run('word', P / 'construction/borrow/build_borrow47.py', main, work / 'word')
word = work / 'word/BORROW-WORD-f47.json.gz'
native_script = main / 'research/gmx-t74-borrow43-handoff/new/run_native_borrow.py'
if a.native_from:
    nat = a.native_from.resolve()
    st = rd(nat / 'STEPS.json')
    assert [x['step'] for x in st] == ['native', 'ledgers', 'local-replay', 'wrappers', 'phases'] and all(x['exit_code'] == 0 for x in st)
    stages.append(dict(name='native', reused=str(nat), exit_code=0, seconds=0.0))
else:
    nat = work / 'native'; run('native', native_script, word, nat)
assert gzip.decompress((nat / 'local/FULL-WORD.json.gz').read_bytes()) == gzip.decompress(word.read_bytes()), 'native directory is not bound to the built word'
local = nat / 'local'
run('plan', P / 'construction/owners/make_plan.py', main, local, P / 'data/OWNERS.json.gz', work / 'PLAN.json')
STABLE = ('LOCAL-CHECKS.json', 'CLIFFORD-RECIPES.json', 'RESIDUAL-ADAPTERS.json')
def bind_view(plan):
    plan = json.loads(json.dumps(plan))
    for d in plan['provider_registry'].values(): d['artifact_sha256'] = {n: h for n, h in d['artifact_sha256'].items() if n in STABLE}
    return plan
assert bind_view(rd(work / 'PLAN.json')) == bind_view(rd(P / 'data/PLAN.json')), 'regenerated PLAN differs from data/PLAN.json'
run('recount', P / 'construction/owners/recount_owners.py', P / 'data/PLAN.json', P / 'data/OWNERS.json.gz')
run('strict', P / 'verify_global.py', '--main', main, '--local', local, '--work', work / 'strict')
run('independent', main / 'research/gauge44-slots-arcflow/independent/price_profile.py', '--main', main, '--profile', P / 'data/PLAN.json',
    '--output', work / 'INDEPENDENT-PRICE.json', cwd=main)
spec = rd(P / 'data/SPEC.json'); replay = rd(work / 'strict/REPLAY.json'); ind = rd(work / 'INDEPENDENT-PRICE.json'); claim = rd(P / 'RESULT.json')
assert replay['status'] == 'PASS_GAUGE44_STRICT_GLOBAL_REPLAY' and replay['kappa_decimal'] == spec['expected']['kappa_decimal'] == claim['kappa_decimal']
assert ind['status'] == 'PASS_FROZEN_PROFILE_TWO_EXACT_ENGINE_AND_OUTER47_CROSSCHECK' and ind['global_profile'] == rd(P / 'data/PLAN.json')['profile']
build = rd(work / 'word/BUILD.json')
result = dict(status='PASS_BORROW47_COMPLETE_T8_VERIFICATION', kappa=replay['kappa'], kappa_decimal=replay['kappa_decimal'],
              complex_lower=replay['complex_lower'], bit_lower=replay['bit_lower'], word_content_sha256=build['word_content_sha256'],
              native_reuse=bool(a.native_from), independent_price_status=ind['status'], stages=stages)
(work / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
print(result['status'], result['kappa_decimal'], flush=True)

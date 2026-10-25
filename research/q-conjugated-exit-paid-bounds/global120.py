"""GLOBAL120 Gamma termination: paid-call lower bounds for the q-conjugated exit and the nominated retirement cut.

Inputs (live, pinned in SOURCES.json):
  research/coded-native-caller/NORMALIZERS.json              112 port normalizers and the 8 mode columns M
  research/native-rank8-helper-provider/PROVIDER.json        8 mode-owner normalizers (input_K literal words)
  research/global120-optimistic-budget/science/inputs/10-BALANCED-J-NATIVE-v1.json   definition of J
  research/global120-optimistic-budget/science/inputs/04-CONTROLLER.json            the 162 physical events
  research/global120-optimistic-budget/evidence/global120-Cstar-v1.json             the 916-job allowance
  research/global120-optimistic-budget/science/inputs/05-INVOICE-NEUTRAL-LITERAL.json paid jobs per controller event

Model and lemmas: PROOF.md.  Every bound below is a rank of class-coordinate rows; a paid call of any rank adds at most
one vector to that span (PROOF.md Lemma A, the output-coset form of the coded-active-query-bound / encoded-exit-obstructions
counting)."""
from fractions import Fraction as Fr
from common import (load, par, canon, phase_class, ZERO, inv_f2_rows, column0, f2_apply, inverse, matmul, rank,
                    kappa_of, sym_monomial_inv, sym_xrot)

NORM = 'research/coded-native-caller/NORMALIZERS.json'
PROV = 'research/native-rank8-helper-provider/PROVIDER.json'
JDEF = 'research/global120-optimistic-budget/science/inputs/10-BALANCED-J-NATIVE-v1.json'
CTRL = 'research/global120-optimistic-budget/science/inputs/04-CONTROLLER.json'
CSTAR = 'research/global120-optimistic-budget/evidence/global120-Cstar-v1.json'
INVOICE = 'research/global120-optimistic-budget/science/inputs/05-INVOICE-NEUTRAL-LITERAL.json'
NP, NM = 112, 8
N = NP + NM


def owners():
    """Centre functionals r_a = row 0 of G_a^-1 (chart w = G u + z) and, as a convention cross-check, G_a e0."""
    nz, pv = load(NORM), load(PROV)
    rows = [rec['G_rows'] for rec in nz['ports']]
    for mb in pv['mode_boundaries']:
        word = mb['input_K']['literal_word']
        P = [x for x in word if x[0] == 'P']
        assert len(P) == 1
        rows.append(P[0][1])
    assert len(rows) == N
    r = [inv_f2_rows(g)[0] for g in rows]
    c = [column0(g) for g in rows]
    for g, ra in zip(rows, r):      # r_a . (G u) = u_0 for every u: r_a really reads the physical centre bit
        assert all(par(ra & f2_apply(g, 1 << t)) == int(t == 0) for t in range(9))
    return rows, r, c, nz['mode_columns']


def code(M):
    """J = [[C, -4W], [3Z, 13 I8]] in the normalized chart: W = M/2, Z = M^T/2, C = I - WZ (BALANCED-J-NATIVE-v1)."""
    defs = load(JDEF)['source_operator_definitions']
    assert defs['J'] == [['C', '-4W'], ['3Z', '13I8']] and defs['C'] == 'I-WZ' and defs['ZW'] == '14I8'
    M = [[Fr(x) for x in col] for col in M]
    assert all(sum(M[k][i] * M[l][i] for i in range(NP)) == (56 if k == l else 0) for k in range(NM) for l in range(NM))
    J = [[Fr(0)] * N for _ in range(N)]
    for a in range(NP):
        for b in range(NP):
            J[a][b] = Fr(int(a == b)) - sum(M[k][a] * M[k][b] for k in range(NM)) / 4
        for k in range(NM):
            J[a][NP + k] = -4 * M[k][a] / 2
            J[NP + k][a] = 3 * M[k][a] / 2
    for k in range(NM):
        J[NP + k][NP + k] = Fr(13)
    return J


class Classes:
    def __init__(self, r):
        self.r = r
        self.ids = {}
        self.memo = {}
        self.q = [self.cid(phase_class((ra, -1))) for ra in r]      # class of q_a^-1 (sign immaterial, see PROOF.md)

    def cid(self, t):
        return None if t == ZERO else self.ids.setdefault(t, len(self.ids))

    def pair(self, a, b):
        """Class of q_a^-1 q_b (entry (a, b) of q^-1 X q)."""
        key = (a, b)
        if key not in self.memo:
            self.memo[key] = self.cid(phase_class((self.r[a], -1), (self.r[b], 1)))
        return self.memo[key]


def conj_rows(cl, X, bank, rows_wanted, cols_wanted=None):
    """Non-free class rows of q^-1 X q applied to input bank `bank`."""
    out = []
    for a in rows_wanted:
        v = {}
        for b in (cols_wanted or range(len(X))):
            if X[a][b] != 0:
                c = cl.pair(a, b)
                if c is not None:
                    v[(c, bank, b)] = v.get((c, bank, b), 0) + X[a][b]
        out.append(v)
    return out


def query_rows(cl, X, banks, a):
    """Row a of q^-1 X (s), s = sum of the given input banks: one class (q_a^-1) per output owner."""
    return {(cl.q[a], bk, b): X[a][b] for b in range(len(X)) for bk in banks if X[a][b] != 0}


def boundary(cl, J, Ji, P, decode, decode_rows=range(N)):
    """Target (5) of the GLOBAL120 chart with B = J: Y0 <- barB^-1 Y0 + q^-1 P B^-1 s, Y1 <- Y1 + q^-1 B P B^-1 s."""
    PJi = matmul(P, Ji)
    JPJi = matmul(J, PJi)
    dec = dict(zip(decode_rows, conj_rows(cl, Ji, 'Y0', decode_rows))) if decode else {}
    rows = []
    for a in range(N):
        y0 = dict(dec.get(a, {}))
        y0.update(query_rows(cl, PJi, ('X0', 'X1'), a))
        rows.append(y0)
        rows.append(query_rows(cl, JPJi, ('X0', 'X1'), a))
    return rank(rows)


def controller_retirements():
    """For contacts k = 0..14 (QUERY_SHEAR(stage, k) reads X_k = bank 2k and adds q^-1 X_k to Y_k = bank 2k+1).

    Checked: LEFT, the return tail and the signed swaps are free invertible scalar programs; after contact k no event
    before the tail moves Y_k's data into another bank (no ADD reads bank 2k+1, no swap involves it) and no free event
    writes it.  Paid calls that act on Y_k after contact k are read per event from the neutral literal invoice: a
    LOCAL or IDLE on pair k pays its listed main jobs on both banks (half on Y_k), and GLOBAL120 adds eight mode jobs
    on Y_k per IDLE ("Both mode banks are included", 01-GLOBAL120-CONDITIONAL-DELTA-v1.md)."""
    ctrl = load(CTRL)
    E = ctrl['physical_events']
    left = ctrl['boundaries']['left_program']
    assert all((e[0] == 'ADD' and e[1] != e[2] and e[3] in ('1', '-1')) or e[0] == 'SWAP' or
               (e[0] == 'SCALE' and e[2] == '-1') for e in left)
    tail = next(i for i, e in enumerate(E) if i > 107 and e[0] == 'ADD')
    assert all(e[0] in ('ADD', 'SIGNED_SWAP') for e in E[tail:])
    assert all(e[0] == 'SIGNED_SWAP' or (e[1] != e[2] and int(e[3]) in (1, -1)) for e in E[tail:])
    inv = load(INVOICE)['chronology']
    assert [r['index'] for r in inv] == list(range(len(E))) and all(r['event'] == e for r, e in zip(inv, E))
    out = []
    for k in range(15):
        qi = [i for i, e in enumerate(E) if e[0] == 'QUERY_SHEAR' and e[2] == k]
        assert len(qi) == 1
        qi = qi[0]
        y = 2 * k + 1
        sig = [[x if i < 3 else int(x) for i, x in enumerate(e)] for e in (E[qi - 2], E[qi - 1], E[qi + 1], E[qi + 2])]
        assert sig == [['ADD', 2 * k, 2 * k + 2, 1], ['ADD', 2 * k + 3, 2 * k + 1, -1],
                       ['ADD', 2 * k, 2 * k + 2, -1], ['ADD', 2 * k + 3, 2 * k + 1, 1]]
        assert sum(inv[qi]['paid'].values()) == NP           # the retained port-query bank of this contact
        after = range(qi + 3, tail)
        reads = [i for i in after if (E[i][0] == 'ADD' and E[i][2] == y) or (E[i][0] == 'SIGNED_SWAP' and y in E[i][1:3])]
        writes = [i for i in after if E[i][0] == 'ADD' and E[i][1] == y]
        on_pair = [i for i in after if (E[i][0] in ('LOCAL', 'QUERY_SHEAR') and E[i][2] == k) or (E[i][0] == 'IDLE' and E[i][1] == k)]
        assert all(E[i][0] == 'IDLE' for i in on_pair)
        calls_on_Yk = sum(sum(inv[i]['paid'].values()) // 2 + NM for i in on_pair)
        out.append(dict(contact=k, query_index=qi, data_leaves_Yk_before_tail=bool(reads), free_writes_to_Yk=len(writes),
                        later_idles_on_pair=[i for i in on_pair], scheduled_paid_calls_on_Yk_after_contact=calls_on_Yk))
    return tail, out


def idle_structure():
    """IDLEs supply only stages a pair does not run itself.  For every IDLE(p, H): H is disjoint from pair p's own FULL
    stages, and its paid jobs are 224 = 2*112 calls of rank 8|H| (pure private-axis stage Gaussians, no central axis).
    For every contact k the last FULL on pair k before it is stage 2k+1 (the chart's B = J1 at contact 0)."""
    E = load(CTRL)['physical_events']
    inv = load(INVOICE)['chronology']
    full = {}
    for e in E:
        if e[0] == 'LOCAL':
            full.setdefault(e[2], set()).add(e[1])
    idles = [(i, e) for i, e in enumerate(E) if e[0] == 'IDLE']
    for i, e in idles:
        pair, H = e[1], set(e[2])
        assert not (H & full.get(pair, set()))
        assert inv[i]['paid'] == {str(8 * len(H)): 2 * NP}
    last = []
    for k in range(15):
        qi = next(i for i, e in enumerate(E) if e[0] == 'QUERY_SHEAR' and e[2] == k)
        st = [e[1] for e in E[:qi] if e[0] == 'LOCAL' and e[2] == k]
        last.append(st[-1])
    assert last == [2 * k + 1 for k in range(15)]
    return len(idles), {p: sorted(v) for p, v in full.items()}


def kappa_census(rows, Ji):
    """Lemma H premise: every off-diagonal exit word q^-1 K_a^-1 K_b q (J^-1_ab != 0) has kappa = 2, read with either
    linear part (G as given, or its inverse transpose); diagonal words have kappa 0.  Q is the central X-rotation."""
    n = 9
    Q = sym_xrot(1)
    def transpose_inv(G):
        Gi = inv_f2_rows(G, n)
        return [sum(((Gi[b] >> i) & 1) << b for b in range(n)) for i in range(n)]
    out = {}
    for conv, Gs in (('G', rows), ('G_inverse_transpose', [transpose_inv(G) for G in rows])):
        hist = {}
        for a in range(N):
            for b in range(N):
                if Ji[a][b] == 0:
                    continue
                # m = K_a^-1 K_b has linear part G_a^-1 G_b; S_m^-1 = S_{K_b}^-1 S_{K_a}
                Ga = Gs[a]
                Gai = inv_f2_rows(Ga, n)
                fa_inv = sym_monomial_inv(Gai, n)            # inverse map of K_a^-1, i.e. S_{K_a}
                fb_inv = sym_monomial_inv(Gs[b], n)
                k = kappa_of([Q, fa_inv, fb_inv, Q], n)
                key = ('diag' if a == b else 'off', k)
                hist[key] = hist.get(key, 0) + 1
        out[conv] = {'%s_kappa_%d' % kk: v for kk, v in sorted(hist.items())}
    return out


def allowance():
    q = [x for x in load(CSTAR)['quotas'] if x['rank'] == 1][0]
    assert q['additional_paid_tariff']['lower'] == '1' and q['additional_paid_alone']['exact_integer']
    return q['additional_paid_alone']['guaranteed_maximum'], load(CSTAR)['C']


def main(verbose=True):
    rows, r, c, M = owners()
    res = {}
    res['distinct_centre_functionals'] = len(set(r))
    res['distinct_G_e0'] = len(set(c))
    # chart-convention robustness: every reading of "the centre line in owner a's chart" separates all 120 owners
    res['distinct_row0_G'] = len({g[0] for g in rows})
    res['distinct_G_inverse_e0'] = len({column0(inv_f2_rows(g)) for g in rows})
    cl = Classes(r)
    res['distinct_nonfree_q_classes'] = len({x for x in cl.q if x is not None})
    # sign independence of classes: -Q = Q mod 2*linear because Q mod 2 is linear (PROOF.md)
    assert all(phase_class((r[a], 1), (r[b], -1)) == phase_class((r[a], -1), (r[b], 1)) for a, b in ((0, 5), (17, 113), (112, 119)))
    J = code(M)
    Ji = inverse(J)
    I = [[Fr(int(a == b)) for b in range(N)] for a in range(N)]
    assert matmul(J, Ji) == I
    res['J_inverse_offdiagonal_support'] = sum(1 for a in range(N) for b in range(N) if a != b and Ji[a][b] != 0)
    pair_classes = {cl.pair(a, b) for a in range(N) for b in range(N) if a != b and Ji[a][b] != 0}
    res['distinct_pair_classes_on_support'] = len(pair_classes)
    res['decode_barJinv_120'] = rank(conj_rows(cl, Ji, 'Y', range(N)))
    A = [[Fr(int(a == b)) - sum(Fr(M[k][a] * M[k][b]) for k in range(NM)) / 4 for b in range(NP)] for a in range(NP)]
    res['decode_old112_calibration'] = rank(conj_rows(cl, inverse(A), 'Y', range(NP)))
    PP = [[Fr(int(a == b and a < NP)) for b in range(N)] for a in range(N)]
    res['cut5_CP_nominated'] = boundary(cl, J, Ji, PP, True)
    res['cut5_CP_ports_only_decode'] = boundary(cl, J, Ji, PP, True, range(NP))
    res['cut5_CP_Y0_already_decoded'] = boundary(cl, J, Ji, PP, False)
    res['cut5_CI_nominated'] = boundary(cl, J, Ji, I, True)
    res['cut5_CI_Y0_already_decoded'] = boundary(cl, J, Ji, I, False)
    unc = []
    for a in range(N):                      # uncoded C_P contact: one response copied to both receivers
        v = {(cl.q[a], bk, a): Fr(1) for bk in ('X0', 'X1')} if a < NP else {}
        unc += [v, dict(v)]
    res['calibration_uncoded_CP_contact'] = rank(unc)
    tail, ret = controller_retirements()
    res['tail_start_index'] = tail
    res['retirements'] = ret
    assert not any(x['data_leaves_Yk_before_tail'] or x['free_writes_to_Yk'] for x in ret)
    need = res['decode_barJinv_120']
    res['retirements_without_scheduled_paid_call_on_Yk'] = sum(1 for x in ret if x['scheduled_paid_calls_on_Yk_after_contact'] == 0)
    res['retirements_whose_IDLE_matches_the_count'] = sum(1 for x in ret if x['scheduled_paid_calls_on_Yk_after_contact'] >= need)
    res['forced_new_or_rerouted_calls'] = sum(max(0, need - x['scheduled_paid_calls_on_Yk_after_contact']) for x in ret)
    res['kappa_census_exit_words'] = kappa_census(rows, Ji)
    n_idle, full = idle_structure()
    res['idles_disjoint_from_own_FULL_stages'] = n_idle
    res['last_FULL_before_contact_k_is_stage_2k_plus_1'] = True
    # PROOF.md Lemma E (in-sector): scheduled IDLE calls act only on other stages' private axes and factor out, so every
    # retirement needs >= 120 non-IDLE paid calls on Y_k after its contact; none are scheduled for any k.
    res['in_sector_required_calls_15_retirements'] = 15 * need
    allow, C = allowance()
    res['allowance_rank1_jobs'] = allow
    res['C'] = C
    res['forced_share_of_allowance'] = str(Fr(res['forced_new_or_rerouted_calls'], allow))
    res['cut5_extra_over_120_retained_calls'] = min(res['cut5_CP_nominated'], res['cut5_CP_Y0_already_decoded'],
                                                    res['cut5_CI_nominated'], res['cut5_CI_Y0_already_decoded']) - 120
    if verbose:
        for k, v in res.items():
            if k != 'retirements':
                print('%-52s %s' % (k, v))
        for x in ret:
            print('   ', x)
    return res


if __name__ == '__main__':
    main()

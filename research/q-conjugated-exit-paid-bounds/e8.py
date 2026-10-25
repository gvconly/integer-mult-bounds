"""Odd-E8 negative exit (PR5464 source, PR5508 target): every paid post-decoder needs one paid call per owner.

Target (PR5508, odd-e8-negative-exit-algebra/PROOF.md): J = Q^-1 E1^-1 E2 E3^-1 Q on the emitted negative bank of
v = 112^3 owners, E_s = L_s^-1 A L_s, A^-1 = K_odd + C_odd in normalized owner coordinates.  The three-stage private
witness lemma (PROOF.md, Lemma G) gives every output owner d = (a,b,c) its own non-free class coordinate, through
the path j -> d -> d -> d with j = (a',b,c).  Its coefficient is (E^-1)_dd * E_dd * (E^-1)_dj
= 1 * E_dd * (K + P)_aa' / 2, so it suffices that each base row of K + P has an off-diagonal entry and E_dd != 0.

Inputs (live): research/coded-native-caller/NORMALIZERS.json (byte-identical to the odd-E8 source copy),
research/odd-e8-negative-exit-algebra/source/K-TEMPLATE.json, research/odd-e8-negative-exit-algebra/RESULT.json."""
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
from itertools import product
from common import load, par

NORM = 'research/coded-native-caller/NORMALIZERS.json'
KT = 'research/odd-e8-negative-exit-algebra/source/K-TEMPLATE.json'
RES = 'research/odd-e8-negative-exit-algebra/RESULT.json'
A_SEL = Fr(605457647, 200000000000)          # selected kappa used by the source screen


def base():
    nz, kt = load(NORM), load(KT)
    M = [[Fr(x) for x in col] for col in nz['mode_columns']]        # M[k][a]
    n = len(M[0])
    assert n == 112 and len(M) == 8
    P = [[sum(M[k][a] * M[k][b] for k in range(8)) / 4 for b in range(n)] for a in range(n)]
    assert all(P[a][a] == 1 for a in range(n))
    assert all(sum(M[k][i] * M[l][i] for i in range(n)) == (56 if k == l else 0) for k in range(8) for l in range(8))
    K = [[Fr(0)] * n for _ in range(n)]
    cubes = kt['literal_cube_partition']
    assert len(cubes) == 14 and sorted(p for c in cubes for p in c['ports']) == list(range(n))
    for ci, c in enumerate(cubes):
        g = c['ports']
        blk = [[Fr(int(i == j)) - sum(M[k][g[i]] * M[k][g[j]] for k in range(8)) / 4 for j in range(8)] for i in range(8)]
        assert blk == [[Fr(x) for x in row] for row in kt['base_K_coefficients'][ci]]   # source block, same port order
        for i in range(8):
            for j in range(8):
                K[g[i]][g[j]] = blk[i][j]
    KK = [[sum(K[a][c] * K[c][b] for c in range(n) if K[a][c] and K[c][b]) for b in range(n)] for a in range(n)]
    assert KK == [[Fr(int(a == b)) for b in range(n)] for a in range(n)]
    KP = [[sum(K[a][c] * P[c][b] for c in range(n) if K[a][c]) for b in range(n)] for a in range(n)]
    assert KP == [[-x for x in row] for row in P]
    assert all(K[a][a] == 0 for a in range(n))
    return nz, kt, M, P, K, n


def labels(nz, kt):
    """Base labels distinct, norm one, not all-ones; all product labels lie in U25 (so 25-bit labels are distinct)."""
    q = [p['q'] for p in nz['ports']]
    assert len(set(q)) == 112 and all(par(x & x) == 1 and x != 511 for x in q)
    orth = kt['orthonormal_basis_raw27']
    assert len(orth) == 25 and all(par(u & w) == int(i == j) for i, u in enumerate(orth) for j, w in enumerate(orth))

    # exact membership test: reduce by an echelon basis of U25
    ech = {}
    for b in orth:
        for hb in sorted(ech, reverse=True):
            if b >> hb & 1:
                b ^= ech[hb]
        if b:
            ech[b.bit_length() - 1] = b

    def member(x):
        for hb in sorted(ech, reverse=True):
            if x >> hb & 1:
                x ^= ech[hb]
        return x == 0

    assert len(ech) == 25
    raw = lambda a, b, c: q[a] | (q[b] << 9) | (q[c] << 18)
    gens = [raw(0, 0, 0)] + [raw(a, 0, 0) ^ raw(0, 0, 0) for a in range(112)] + \
           [raw(0, b, 0) ^ raw(0, 0, 0) for b in range(112)] + [raw(0, 0, c) ^ raw(0, 0, 0) for c in range(112)]
    assert all(member(x) for x in gens)          # hence every raw product label lies in U25
    return len(set(q))


def e_diagonal():
    """E = K_odd - C_odd (I - C_odd)^-1; its diagonal from the sector projectors (diag of P/14 is 1/14)."""
    p = Fr(1, 14)
    lam2, lam3 = Fr(196, 2), Fr(3 * 196 - 2744, 2)
    d2, d3 = 3 * p * p * (1 - p), p ** 3
    return -(lam2 / (1 - lam2) * d2 + lam3 / (1 - lam3) * d3), lam2, lam3


def margin(ratio):
    """(2/3) 73^t - (8^t + 9^t + 15^t + 16^t) - ratio 25^t, t = 1 - a: the sources' per-port per-local-call screen."""
    getcontext().prec = 50
    t = 1 - Decimal(A_SEL.numerator) / Decimal(A_SEL.denominator)
    pw = lambda r: (Decimal(r).ln() * t).exp()
    rat = Decimal(ratio.numerator) / Decimal(ratio.denominator)
    return Decimal(2) / 3 * pw(73) - (pw(8) + pw(9) + pw(15) + pw(16)) - rat * pw(25)


def main(verbose=True):
    nz, kt, M, P, K, n = base()
    res = {'base_labels_distinct': labels(nz, kt)}
    KP = [[K[a][b] + P[a][b] for b in range(n)] for a in range(n)]
    wit = [sum(1 for b in range(n) if b != a and KP[a][b] != 0) for a in range(n)]
    res['min_offdiagonal_support_of_K_plus_P_row'] = min(wit)
    Edd, l2, l3 = e_diagonal()
    src = load(RES)
    assert Fr(src['E_diagonal']) == Edd and Fr(src['E_inverse_diagonal']) == 1
    res['E_diagonal'] = str(Edd)
    res['C_odd_nonzero_eigenvalues'] = [str(l2), str(l3)]
    # the PR5508 example pair is one of our witnesses: d=(1,0,0), j=(0,0,0) has (E^-1)_dj = (K+P)_10/2
    assert KP[1][0] / 2 == Fr(src['E_inverse_cross_entry'])
    v = n ** 3
    res['owners_v'] = v
    res['owners_with_private_witness'] = v if min(wit) > 0 and Edd != 0 else None
    res['paid_postdecoder_lower_bound_calls'] = res['owners_with_private_witness']
    res['naive_postdecoder_calls'] = 2 * v
    # per port per local call: the three-call word emits one negative bank of v owners; screens divide by 3v
    assert Fr(61440, 3 * v) == Fr(5, 343)
    m = margin(Fr(5, 343))
    res['source_margin_per_port_per_local_call'] = str(m)[:22]
    res['decoder_floor_per_port_per_local_call'] = '1/3'
    res['floor_over_margin'] = str((Decimal(1) / 3 / m))[:8]
    m120 = margin(Fr(43, 3375))
    res['full120_declared_margin_per_port_per_local_call'] = str(m120)[:22]
    res['full120_floor_over_margin_if_premises_hold'] = str((Decimal(1) / 3 / m120))[:8]
    if verbose:
        for k, x in res.items():
            print('%-52s %s' % (k, x))
    return res


if __name__ == '__main__':
    main()

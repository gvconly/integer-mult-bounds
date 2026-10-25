"""Shared exact helpers: repository-relative input reading, F2/Z4 phase classes, rational rank and inverse.

Phase classes.  On a nine-bit chart w, the rank-one central Gaussian of an owner whose centre functional is r acts as
i^(s[r.w]), with [.] the 0/1 parity.  Two phase words lie in the same free (monomial) coset exactly when their
difference is 2*linear + constant mod 4 (a Pauli sign and a global phase).  `phase_class` returns that canonical
representative as a 512-tuple; ZERO is the free class.
"""
import json
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
W9 = range(512)


def load(rel):
    return json.loads((ROOT / rel).read_text())


def par(x):
    return bin(x).count('1') & 1


def canon(Q):
    """Representative of a Z4-valued function on F2^9 modulo constants and 2*(linear functionals)."""
    q0 = Q[0]
    L = sum(1 << t for t in range(9) if ((Q[1 << t] - q0) % 4) in (2, 3))
    return tuple((Q[w] - q0 - 2 * par(L & w)) % 4 for w in W9)


ZERO = tuple([0] * 512)


def gauss(r, s):
    """Z4 exponent table of i^(s[r.w])."""
    return [(s * par(r & w)) % 4 for w in W9]


def phase_class(*terms):
    """Class of the product of Gaussians given as (r, s) pairs."""
    tot = [0] * 512
    for r, s in terms:
        g = gauss(r, s)
        tot = [(a + b) % 4 for a, b in zip(tot, g)]
    return canon(tot)


def inv_f2_rows(rows, n=9):
    """Rows of the inverse of the n x n F2 matrix whose row b is the bitmask rows[b]."""
    A = [[(rows[i] >> j) & 1 for j in range(n)] + [int(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c])
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c]:
                A[r] = [x ^ y for x, y in zip(A[r], A[c])]
    return [sum(A[i][n + j] << j for j in range(n)) for i in range(n)]


def f2_apply(rows, x):
    return sum(par(r & x) << b for b, r in enumerate(rows))


def column0(rows):
    """G e0 for the matrix with the given rows (bit b of the result is entry (b, 0))."""
    return sum((r & 1) << b for b, r in enumerate(rows))


def inverse(A):
    n = len(A)
    X = [list(row) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if X[r][c] != 0)
        X[c], X[p] = X[p], X[c]
        pv = X[c][c]
        X[c] = [x / pv for x in X[c]]
        for r in range(n):
            if r != c and X[r][c] != 0:
                f = X[r][c]
                X[r] = [x - f * y for x, y in zip(X[r], X[c])]
    return [row[n:] for row in X]


def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum(x * y for x, y in zip(r, c) if x and y) for c in Bt] for r in A]


def rank(rows):
    """Exact rank over Q of sparse rows given as {coordinate: Fraction}."""
    piv = {}
    rk = 0
    for r in rows:
        r = {k: v for k, v in r.items() if v != 0}
        while r:
            k = min(r)
            if k in piv:
                f = r[k] / piv[k][k]
                for kk, vv in piv[k].items():
                    nv = r.get(kk, 0) - f * vv
                    if nv:
                        r[kk] = nv
                    else:
                        r.pop(kk, None)
            else:
                piv[k] = r
                rk += 1
                break
    return rk


# ---- Gaussian-rank budget (PROOF.md Lemma H): Lagrangian x-rank of Clifford words, over F2, n binary axes.
# A Clifford word U acts on Pauli labels (x|z) by a symplectic map S_U; free monomials are exactly the Cliffords whose
# S preserves L_Z = {(0|z)}.  kappa(U) = dim of the x-projection of S_U^-1 L_Z; it is a right-coset invariant, is
# unchanged by free factors on either side, and kappa(gamma U) <= kappa(U) + kappa(gamma).

def f2_rank(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def sym_monomial_inv(G_rows, n):
    """S_m^-1 for the address monomial x -> G x (Pauli/diagonal-Pauli phases only): (x|z) -> (G^-1 x | G^T z)."""
    Gi = inv_f2_rows(G_rows, n)
    cols = [sum(((G_rows[b] >> i) & 1) << b for b in range(n)) for i in range(n)]   # column i of G
    def f(x, z):
        return f2_apply(Gi, x), sum(par(cols[i] & z) << i for i in range(n))       # (G^T z)_i = (col_i of G) . z
    return f


def sym_xrot(axes_mask):
    """S for a product of single-axis X-rotations (positive C Gaussians) on the given axes: (x|z) -> (x + z|S, z)."""
    def f(x, z):
        return x ^ (z & axes_mask), z
    return f


def kappa_of(maps, n):
    """kappa(U) for U = U_1 U_2 ... U_k (U_k acts first on data).  S_U = S_{U_1} ... S_{U_k}, so applying S_U^-1 to
    L_Z applies S_{U_1}^-1 first; `maps` lists [S_{U_1}^-1, ..., S_{U_k}^-1] in that application order."""
    xs = []
    for i in range(n):
        x, z = 0, 1 << i
        for f in maps:
            x, z = f(x, z)
        xs.append(x)
    return f2_rank(xs)

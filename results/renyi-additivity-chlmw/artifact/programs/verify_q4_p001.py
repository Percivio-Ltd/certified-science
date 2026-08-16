#!/usr/bin/env python3
"""Independent analyst replay of the Q4 p=1/100 Renyi-additivity candidate.

Clean-room implementation. Inputs consumed:
  1. The frozen candidate bytes renyi-p001-candidate.json
     (sha256 980b1af4df2f4984d8b8fce7822a67073cbd3bb3b45e8426c97c08a51a9de419):
     delta, p, and the two rational decomposable witnesses Q_A, Q_B.
  2. The printed subspace bases R and S of arXiv:0712.3628v3, transcribed by
     the analyst directly from the paper's LaTeX source (counter-rank5.tex,
     sha256 a358179291e577f26f4086164630b693775b9aa7058008dc589082086816f448,
     lines 453-517 of the section "An explicit construction in small
     dimension").
  3. The declared conventions of the terminal dossier (progress/answer.json):
     unnormalized-Choi trace-preservation convention Tr_B(J) = I_4, channels
     normalized by inverse input marginals, natural log,
     S_p(rho) = ln(Tr rho^p)/(1-p), p = 1/100, floor decomposition
     J - delta*I = P + Q^{Gamma_B}.

No participant-authored program was read or executed. All arithmetic is
exact: rationals (fractions.Fraction) and the number field
Q(i, sqrt2, sqrt3, sqrt5) held as 16 rational components. PSD certification
is exact Hermitian LDL over that field with exact sign decisions -- a
different proof path from the participant's rational-approximant shifted-LDL
with Choi-entry operator-error allowance (no error allowance is needed here
because the Choi matrices are represented exactly). 100th roots are enclosed
by outward-rounded 192-bit dyadic binary search on integers.

Mathematical chain being verified (analyst's own derivation):
  (i)   J Hermitian, Tr_B J = I_4, J PSD  =>  N(X) = Tr_A[J (X^T x I)] is
        CPTP, and <psi|N(phi phi*)|psi> = <conj(phi) x psi|J|conj(phi) x psi>.
  (ii)  J - delta*I = P + Q^{Gamma_B} with P, Q PSD  =>  for every unit
        product vector u x v,  <u x v|J|u x v> - delta =
        <u x v|P|u x v> + <u x conj(v)|Q|u x conj(v)> >= 0,
        hence every channel output N(rho) >= delta*I (mixed inputs by
        convexity).
  (iii) For 0<p<1 the map lam -> sum_i lam_i^p is symmetric concave, hence
        Schur-concave; every point of {lam_i >= delta, sum lam = 1} is
        majorized by (1-2delta, delta, delta); therefore
        Tr[N(rho)^p] >= (1-2delta)^p + 2 delta^p =: A  for every input rho,
        so S_p^min(N) >= ln(A)/(1-p) for each channel.
  (iv)  For the accepted normalized entangled input, the exact joint output
        sigma is rational with charpoly t (46t-5)^2 (161t-18)^2 (322t-51)^2
        (322t-39)^2 (monic), i.e. spectrum {0, 5/46 x2, 18/161 x2,
        51/322 x2, 39/322 x2}; hence
        S_p^min(N_A x N_B) <= S_p(sigma) = ln(B)/(1-p), B = Tr sigma^p.
  (v)   A_lower^2 > B_upper with outward rounding  =>  B < A^2
        =>  S_p^min(N_A x N_B) < S_p^min(N_A) + S_p^min(N_B).
"""

import json
import sys
import math
import hashlib
from fractions import Fraction as Fr

# ----------------------------------------------------------------------------
# Exact arithmetic in Q(i, sqrt2, sqrt3, sqrt5): a + i*b with a, b in the
# real field Q(sqrt2, sqrt3, sqrt5) held as 8 rational components on the
# basis (1, s2, s3, s5, s6, s10, s15, s30).
# ----------------------------------------------------------------------------

RADICAND = (1, 2, 3, 5, 6, 10, 15, 30)
_Z8 = tuple(Fr(0) for _ in range(8))

def _mul8(x, y):
    out = [Fr(0)] * 8
    for i in range(8):
        xi = x[i]
        if xi == 0:
            continue
        ri = RADICAND[i]
        for j in range(8):
            yj = y[j]
            if yj == 0:
                continue
            r = ri * RADICAND[j]
            s, f = 1, r
            for pr in (2, 3, 5):
                p2 = pr * pr
                while f % p2 == 0:
                    f //= p2
                    s *= pr
            out[RADICAND.index(f)] += xi * yj * s
    return tuple(out)

def sign_q(a):
    return (a > 0) - (a < 0)

def sign2(a, b, r=2):
    """Exact sign of a + b*sqrt(r), r squarefree > 1."""
    if b == 0:
        return sign_q(a)
    if a == 0:
        return sign_q(b)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    d = a * a - r * b * b
    if d == 0:
        return 0
    return sign_q(a) * sign_q(d)

def _sign_q23(y):
    """Sign of y0 + y1 s2 + y2 s3 + y3 s6 (components over (1,s2,s3,s6))."""
    g = (y[0], y[1])            # g0 + g1 s2
    h = (y[2], y[3])            # (h0 + h1 s2) s3
    sg = sign2(*g)
    sh = sign2(*h)
    if sh == 0:
        return sg
    if sg == 0:
        return sh
    if sg == sh:
        return sg
    g2 = (g[0] * g[0] + 2 * g[1] * g[1], 2 * g[0] * g[1])
    h2 = (h[0] * h[0] + 2 * h[1] * h[1], 2 * h[0] * h[1])
    d = (g2[0] - 3 * h2[0], g2[1] - 3 * h2[1])
    sd = sign2(*d)
    if sd == 0:
        return 0
    return sg if sd > 0 else sh

def sign8(x):
    """Exact sign of sum_k x_k sqrt(RADICAND_k)."""
    u = (x[0], x[1], x[2], x[4])            # Q(s2,s3) part
    v = (x[3], x[5], x[6], x[7])            # coefficient of s5 over (1,s2,s3,s6)
    su = _sign_q23(u)
    sv = _sign_q23(v)
    if sv == 0:
        return su
    if su == 0:
        return sv
    if su == sv:
        return su
    u8 = (u[0], u[1], u[2], Fr(0), u[3], Fr(0), Fr(0), Fr(0))
    v8 = (v[0], v[1], v[2], Fr(0), v[3], Fr(0), Fr(0), Fr(0))
    u2 = _mul8(u8, u8)
    v2 = _mul8(v8, v8)
    d = tuple(a - 5 * b for a, b in zip(u2, v2))
    assert d[3] == 0 and d[5] == 0 and d[6] == 0 and d[7] == 0
    sd = _sign_q23((d[0], d[1], d[2], d[4]))
    if sd == 0:
        return 0
    return su if sd > 0 else sv

class F:
    """Element of Q(i, sqrt2, sqrt3, sqrt5)."""
    __slots__ = ("a", "b")

    def __init__(self, a=_Z8, b=_Z8):
        self.a = a
        self.b = b

    @staticmethod
    def from_q(q):
        z = list(_Z8)
        z[0] = Fr(q)
        return F(tuple(z))

    @staticmethod
    def from_qi(re, im):
        za, zb = list(_Z8), list(_Z8)
        za[0], zb[0] = Fr(re), Fr(im)
        return F(tuple(za), tuple(zb))

    @staticmethod
    def omega(power=1):
        """omega = -1/2 + i*sqrt3/2; power in {1, 2} (omega^2 = conj)."""
        za, zb = list(_Z8), list(_Z8)
        za[0] = Fr(-1, 2)
        zb[2] = Fr(1, 2) if power == 1 else Fr(-1, 2)
        return F(tuple(za), tuple(zb))

    @staticmethod
    def sqrt_rat(x):
        """Exact sqrt of a positive rational, when it lies in the field."""
        num, den = x.numerator, x.denominator
        nd = num * den
        s, f = 1, nd
        for pr in (2, 3, 5):
            p2 = pr * pr
            while f % p2 == 0:
                f //= p2
                s *= pr
        if f not in RADICAND:
            raise AssertionError(f"sqrt({x}) not in field")
        z = list(_Z8)
        z[RADICAND.index(f)] = Fr(s, den)
        return F(tuple(z))

    def __add__(self, o):
        return F(tuple(x + y for x, y in zip(self.a, o.a)),
                 tuple(x + y for x, y in zip(self.b, o.b)))

    def __sub__(self, o):
        return F(tuple(x - y for x, y in zip(self.a, o.a)),
                 tuple(x - y for x, y in zip(self.b, o.b)))

    def __neg__(self):
        return F(tuple(-x for x in self.a), tuple(-x for x in self.b))

    def __mul__(self, o):
        aa = _mul8(self.a, o.a)
        bb = _mul8(self.b, o.b)
        ab = _mul8(self.a, o.b)
        ba = _mul8(self.b, o.a)
        return F(tuple(x - y for x, y in zip(aa, bb)),
                 tuple(x + y for x, y in zip(ab, ba)))

    def conj(self):
        return F(self.a, tuple(-x for x in self.b))

    def is_zero(self):
        return all(v == 0 for v in self.a) and all(v == 0 for v in self.b)

    def is_real(self):
        return all(v == 0 for v in self.b)

    def is_rational(self):
        return self.is_real() and all(v == 0 for v in self.a[1:])

    def rational(self):
        assert self.is_rational()
        return self.a[0]

    def real_sign(self):
        assert self.is_real(), "sign of non-real element requested"
        return sign8(self.a)

    def inv(self):
        num = self.conj()
        d8 = tuple(x + y for x, y in
                   zip(_mul8(self.a, self.a), _mul8(self.b, self.b)))
        real_num = tuple(Fr(1) if i == 0 else Fr(0) for i in range(8))
        d = d8
        for kill in (5, 3, 2):
            dc = tuple((-v if RADICAND[i] % kill == 0 else v)
                       for i, v in enumerate(d))
            real_num = _mul8(real_num, dc)
            d = _mul8(d, dc)
        assert all(d[i] == 0 for i in range(1, 8)), "inversion chain failed"
        q = d[0]
        assert q != 0, "attempted inversion of zero"
        real_num = tuple(x / q for x in real_num)
        return F(_mul8(num.a, real_num), _mul8(num.b, real_num))

    def __eq__(self, o):
        return self.a == o.a and self.b == o.b

    def to_complex(self):
        rads = [math.sqrt(r) for r in RADICAND]
        re = sum(float(c) * r for c, r in zip(self.a, rads))
        im = sum(float(c) * r for c, r in zip(self.b, rads))
        return complex(re, im)

F_ZERO = F()
F_ONE = F.from_q(1)

# ----------------------------------------------------------------------------
# Exact Hermitian LDL positive-definiteness certification
# ----------------------------------------------------------------------------

def ldl_is_pd(M, n):
    """True iff the n x n Hermitian F-matrix M is positive definite,
    by exact LDL^dagger with exact pivot signs. Early-exits on the first
    nonpositive pivot."""
    L = [[F_ZERO] * n for _ in range(n)]
    D = [F_ZERO] * n
    for j in range(n):
        acc = M[j][j]
        for k in range(j):
            acc = acc - L[j][k] * L[j][k].conj() * D[k]
        if not acc.is_real():
            return False, f"pivot {j} not real (Hermiticity broken)"
        if acc.real_sign() <= 0:
            return False, f"pivot {j} nonpositive"
        D[j] = acc
        dinv = acc.inv()
        for i in range(j + 1, n):
            acc2 = M[i][j]
            for k in range(j):
                acc2 = acc2 - L[i][k] * L[j][k].conj() * D[k]
            L[i][j] = acc2 * dinv
    return True, "all pivots positive"

def is_hermitian(M, n):
    return all((M[i][j] - M[j][i].conj()).is_zero()
               for i in range(n) for j in range(i, n))

def shift_diag(M, q):
    out = [row[:] for row in M]
    qq = F.from_q(q)
    for i in range(len(M)):
        out[i][i] = out[i][i] + qq
    return out

# ----------------------------------------------------------------------------
# Printed subspace bases (analyst transcription of arXiv:0712.3628v3)
# Vectors of C^4 (x) C^3 as 3x4 matrices: |v> = sum_{b,a} M[b][a] |a>_A |b>_B.
# ----------------------------------------------------------------------------

R_MATS = [
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
    [[1, 0, 0, 0], [0, "w", 0, 0], [0, 0, "w2", 0]],
    [[0, 1, 0, 0], [0, 0, "w2", 0], [0, 0, 0, "w"]],
    [[0, 0, 1, 0], [0, 0, 0, -1], [0, 0, 0, 0]],
    [[0, 0, 0, 1], [0, 0, 0, 0], [-1, 0, 0, 0]],
]
S_MATS = [
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, -1, 0, 0]],
    [[1, 0, 0, 0], [0, "w2", 0, 0], [0, 0, "w", 0]],
    [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    [[0, 1, 0, 0], [0, 0, "w", 0], [0, 0, 0, "w2"]],
    [[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]],
    [[0, 0, 0, 1], [0, 0, 0, 0], [1, 0, 0, 0]],
]

def entry_to_F(e, conjugate):
    if e == "w":
        return F.omega(2 if conjugate else 1)
    if e == "w2":
        return F.omega(1 if conjugate else 2)
    return F.from_q(e)

def make_vectors(mats, conjugate, index_order):
    vecs = []
    for Mt in mats:
        v = [F_ZERO] * 12
        for b in range(3):
            for a in range(4):
                fe = entry_to_F(Mt[b][a], conjugate)
                if fe.is_zero():
                    continue
                idx = 3 * a + b if index_order == "A-major" else 4 * b + a
                v[idx] = fe
        vecs.append(v)
    return vecs

def gram(vs):
    n = len(vs)
    return [[sum((vs[i][t].conj() * vs[j][t] for t in range(12)), F_ZERO)
             for j in range(n)] for i in range(n)]

def projector(vecs):
    P = [[F_ZERO] * 12 for _ in range(12)]
    for v in vecs:
        n2 = sum((v[t].conj() * v[t] for t in range(12)), F_ZERO)
        inv_n2 = n2.inv()
        nz = [t for t in range(12) if not v[t].is_zero()]
        for i in nz:
            for j in nz:
                P[i][j] = P[i][j] + v[i] * v[j].conj() * inv_n2
    return P

def mat_mul(X, Y, n):
    return [[sum((X[i][t] * Y[t][j] for t in range(n)), F_ZERO)
             for j in range(n)] for i in range(n)]

def partial_trace_B(M, index_order):
    W = [[F_ZERO] * 4 for _ in range(4)]
    for a in range(4):
        for a2 in range(4):
            acc = F_ZERO
            for b in range(3):
                i = 3 * a + b if index_order == "A-major" else 4 * b + a
                j = 3 * a2 + b if index_order == "A-major" else 4 * b + a2
                acc = acc + M[i][j]
            W[a][a2] = acc
    return W

def build_choi(P, index_order):
    """J = (W^{-1/2} x I) P (W^{-1/2} x I), W = Tr_B P (rational diagonal)."""
    Wm = partial_trace_B(P, index_order)
    w = []
    for a in range(4):
        for a2 in range(4):
            if a == a2:
                assert Wm[a][a2].is_rational(), "marginal not rational"
            else:
                assert Wm[a][a2].is_zero(), "marginal not diagonal"
        w.append(Wm[a][a].rational())
    J = [[F_ZERO] * 12 for _ in range(12)]
    inv_sqrt = {}
    for a in range(4):
        for a2 in range(4):
            key = (a, a2)
            inv_sqrt[key] = F.sqrt_rat(w[a] * w[a2]).inv()
    for a in range(4):
        for b in range(3):
            i = 3 * a + b if index_order == "A-major" else 4 * b + a
            for a2 in range(4):
                for b2 in range(3):
                    j = 3 * a2 + b2 if index_order == "A-major" else 4 * b2 + a2
                    if not P[i][j].is_zero():
                        J[i][j] = P[i][j] * inv_sqrt[(a, a2)]
    return J, w

def check_TP(J, index_order):
    W = partial_trace_B(J, index_order)
    return all((W[a][a2] - (F_ONE if a == a2 else F_ZERO)).is_zero()
               for a in range(4) for a2 in range(4))

def partial_transpose_B(M, index_order):
    def idx(a, b):
        return 3 * a + b if index_order == "A-major" else 4 * b + a
    out = [[F_ZERO] * 12 for _ in range(12)]
    for a in range(4):
        for b in range(3):
            for a2 in range(4):
                for b2 in range(3):
                    out[idx(a, b)][idx(a2, b2)] = M[idx(a, b2)][idx(a2, b)]
    return out

# ----------------------------------------------------------------------------
# Candidate parsing
# ----------------------------------------------------------------------------

def parse_candidate(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    sha = hashlib.sha256(raw).hexdigest()
    d = json.loads(raw)
    assert d["schema_version"] == "renyi-block-positive-candidate-v1"
    wits = {}
    for w in d["witnesses"]:
        Q = [[(Fr(int(nr), int(dr)), Fr(int(ni), int(di)))
              for (nr, dr), (ni, di) in row] for row in w["Q"]]
        assert len(Q) == 12 and all(len(r) == 12 for r in Q)
        wits[w["name"]] = {"Q": Q, "mu": Fr(w["mu"]), "nu": Fr(w["nu"]),
                           "delta": Fr(w["delta"]),
                           "J_err": Fr(w["J_operator_error_bound"])}
    return d, Fr(d["delta"]), Fr(d["p"]), wits, sha

def q_hermitian(Q):
    return all(Q[i][j][0] == Q[j][i][0] and Q[i][j][1] == -Q[j][i][1]
               for i in range(12) for j in range(i, 12))

def q_to_F(Q):
    return [[F.from_qi(re, im) for re, im in row] for row in Q]

# ----------------------------------------------------------------------------
# Joint output
# ----------------------------------------------------------------------------

def joint_output_sigma(J1, J2, w1, w2, index_order, trial):
    """sigma = (N1 x N2)(|Psi><Psi|), N(X) = Tr_A[J (X^T x I)].
    trial 'weighted': Psi propto sum_j sqrt(w1_j w2_j) |jj>;
    trial 'maxent':   Psi propto sum_j |jj>."""
    if trial == "weighted":
        tprod = [w1[a] * w2[a] for a in range(4)]
    else:
        tprod = [Fr(1)] * 4
    norm2 = sum(tprod)
    inv_norm2 = F.from_q(Fr(1) / norm2)
    def idx(a, b):
        return 3 * a + b if index_order == "A-major" else 4 * b + a
    amp = {}
    for j in range(4):
        for k in range(4):
            amp[(j, k)] = F.sqrt_rat(tprod[j] * tprod[k])
    sig = [[F_ZERO] * 9 for _ in range(9)]
    for b1 in range(3):
        for b2 in range(3):
            for b1p in range(3):
                for b2p in range(3):
                    acc = F_ZERO
                    for j in range(4):
                        for k in range(4):
                            e1 = J1[idx(j, b1)][idx(k, b1p)]
                            if e1.is_zero():
                                continue
                            e2 = J2[idx(j, b2)][idx(k, b2p)]
                            if e2.is_zero():
                                continue
                            acc = acc + e1 * e2 * amp[(j, k)]
                    sig[3 * b1 + b2][3 * b1p + b2p] = acc * inv_norm2
    return sig

def charpoly_frac(M9):
    """det(tI - M) for a 9x9 Fraction matrix by Faddeev-LeVerrier.
    Returns monic coefficient list c[0..9], p(t) = sum c[k] t^k."""
    n = 9
    B = [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    a = []
    for k in range(1, n + 1):
        A = [[sum(M9[i][t] * B[t][j] for t in range(n)) for j in range(n)]
             for i in range(n)]
        ak = sum(A[i][i] for i in range(n)) / k
        a.append(ak)
        B = [[A[i][j] - (ak if i == j else 0) for j in range(n)]
             for i in range(n)]
    out = [Fr(0)] * (n + 1)
    out[n] = Fr(1)
    for k in range(1, n + 1):
        out[n - k] = -a[k - 1]
    return out

def poly_mul(x, y):
    out = [Fr(0)] * (len(x) + len(y) - 1)
    for i, xv in enumerate(x):
        if xv == 0:
            continue
        for j, yv in enumerate(y):
            out[i + j] += xv * yv
    return out

def spectrum_charpoly(spectrum):
    """Monic charpoly with the given (eigenvalue, multiplicity) list."""
    poly = [Fr(1)]
    for ev, mult in spectrum:
        num, den = ev.numerator, ev.denominator
        for _ in range(mult):
            poly = poly_mul(poly, [Fr(-num), Fr(den)])
    lead = poly[-1]
    return [c / lead for c in poly]

CLAIMED_SPECTRUM = [
    (Fr(0), 1),
    (Fr(5, 46), 2),
    (Fr(18, 161), 2),
    (Fr(51, 322), 2),
    (Fr(39, 322), 2),
]

# ----------------------------------------------------------------------------
# Outward-rounded dyadic root enclosures, logs, and the strict comparison
# ----------------------------------------------------------------------------

def rootP_bounds(x, power, bits):
    """l, u with l <= x^(1/power) <= u, dyadic with denominator 2^bits,
    for rational 0 < x <= 1."""
    num, den = x.numerator, x.denominator
    assert 0 < x <= 1
    target = num * (1 << (power * bits))
    lo, hi = 0, 1 << bits
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if pow(mid, power) * den <= target:
            lo = mid
        else:
            hi = mid - 1
    l = Fr(lo, 1 << bits)
    u = l if pow(lo, power) * den == target else Fr(lo + 1, 1 << bits)
    return l, u

def ln_bounds(x_lo, x_hi, terms=120):
    """Rational enclosure of [ln x_lo, ln x_hi], x_lo <= x_hi, both > 0,
    via atanh series with geometric tail bound; monotone endpointwise."""
    def enclose(x):
        t = (x - 1) / (x + 1)
        t2 = t * t
        s = Fr(0)
        tp = t
        for k in range(terms):
            s += tp / (2 * k + 1)
            tp *= t2
        tail = 2 * abs(t) ** (2 * terms + 1) / ((2 * terms + 1) * (1 - t2))
        v = 2 * s
        return (v, v + tail) if t >= 0 else (v - tail, v)
    lo, _ = enclose(x_lo)
    _, hi = enclose(x_hi)
    return lo, hi

def strict_compare(delta, spectrum, p, bits):
    """Certify A^2 > B with outward rounding; returns (status, info) with
    status in {'CERTIFIED', 'ABSTAIN-PRECISION', 'FALSE'}."""
    assert p.numerator == 1
    power = p.denominator
    A_lo = A_hi = Fr(0)
    for x, mult in ((Fr(1) - 2 * delta, 1), (delta, 2)):
        l, u = rootP_bounds(x, power, bits)
        A_lo += mult * l
        A_hi += mult * u
    B_lo = B_hi = Fr(0)
    for ev, mult in spectrum:
        if ev == 0:
            continue
        l, u = rootP_bounds(ev, power, bits)
        B_lo += mult * l
        B_hi += mult * u
    gap_lo = A_lo * A_lo - B_hi
    gap_hi = A_hi * A_hi - B_lo
    info = {"power": power, "bits": bits,
            "A_enclosure": [fstr(A_lo), fstr(A_hi)],
            "B_enclosure": [fstr(B_lo), fstr(B_hi)],
            "trace_power_gap_enclosure": [fstr(gap_lo), fstr(gap_hi)]}
    if gap_lo > 0:
        lnA = ln_bounds(A_lo, A_hi)
        lnB = ln_bounds(B_lo, B_hi)
        m_lo = (2 * lnA[0] - lnB[1]) / (1 - p)
        m_hi = (2 * lnA[1] - lnB[0]) / (1 - p)
        info["entropy_margin_nats_enclosure"] = [fstr(m_lo), fstr(m_hi)]
        return "CERTIFIED", info
    if gap_hi <= 0:
        return "FALSE", info
    return "ABSTAIN-PRECISION", info

def fstr(fr, digits=25):
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    ip = fr.numerator // fr.denominator
    rem = fr - ip
    ds = []
    for _ in range(digits):
        rem *= 10
        d = rem.numerator // rem.denominator
        ds.append(str(d))
        rem -= d
    return f"{sign}{ip}.{''.join(ds)}"

# ----------------------------------------------------------------------------
# Main verification
# ----------------------------------------------------------------------------

def run(candidate_path, sigma_claim_path, report_path,
        delta_override=None, p_override=None, bits=192, mutate=None,
        spectrum=None):
    report = {"schema": "analyst-q4-p001-verification-report-v1"}
    cand, delta, p, wits, cand_sha = parse_candidate(candidate_path)
    if delta_override is not None:
        delta = Fr(delta_override)
    if p_override is not None:
        p = Fr(p_override)
    spectrum = spectrum or CLAIMED_SPECTRUM
    report["candidate_sha256"] = cand_sha
    report["delta"], report["p"] = str(delta), str(p)
    report["mutation"] = mutate

    if mutate == "corrupt-Q-entry":
        re, im = wits["A"]["Q"][0][1]
        wits["A"]["Q"][0][1] = (re + Fr(1, 10**9), im)
    elif mutate == "identical-witnesses":
        wits["B"]["Q"] = wits["A"]["Q"]
    elif mutate == "swap-witnesses":
        wits["A"], wits["B"] = wits["B"], wits["A"]

    herm = {n: q_hermitian(w["Q"]) for n, w in wits.items()}
    report["witness_hermitian"] = herm
    if not all(herm.values()):
        report["verdict"] = "REJECT: witness not Hermitian"
        emit(report, report_path)
        return report

    # Q PSD checks are convention-independent
    q_psd = {}
    for name, w in wits.items():
        QF = q_to_F(w["Q"])
        ok0, why0 = ldl_is_pd(QF, 12)
        oknu, _ = ldl_is_pd(shift_diag(QF, -w["nu"]), 12)
        q_psd[name] = {"Q_pd": ok0, "Q_minus_nu_pd": oknu, "detail": why0}
    report["witness_Q_psd"] = q_psd
    if not all(v["Q_pd"] for v in q_psd.values()):
        report["verdict"] = "REJECT: witness Q not positive definite"
        emit(report, report_path)
        return report

    # scan conventions
    scan = {}
    winners = {}
    chois_by_order = {}
    for index_order in ("A-major", "B-major"):
        vecs = {"R": make_vectors(R_MATS, False, index_order),
                "Sbar": make_vectors(S_MATS, True, index_order),
                "S": make_vectors(S_MATS, False, index_order)}
        G = gram(vecs["R"] + vecs["S"])
        ortho = all(i == j or G[i][j].is_zero()
                    for i in range(12) for j in range(12))
        norms_rat = all(G[i][i].is_rational() for i in range(12))
        scan[f"transcription-{index_order}"] = {
            "R+S pairwise orthogonal": ortho,
            "norms rational": norms_rat,
            "norms": [str(G[i][i].rational()) for i in range(12)],
        }
        assert ortho and norms_rat, "transcription self-check failed"
        chois = {}
        weights = {}
        for name in ("R", "Sbar", "S"):
            P = projector(vecs[name])
            PH = is_hermitian(P, 12)
            P2 = mat_mul(P, P, 12)
            idem = all((P2[i][j] - P[i][j]).is_zero()
                       for i in range(12) for j in range(12))
            assert PH and idem, f"projector self-check failed: {name}"
            J, w = build_choi(P, index_order)
            assert check_TP(J, index_order), f"TP failed: {name}"
            assert is_hermitian(J, 12), f"J not Hermitian: {name}"
            chois[name] = J
            weights[name] = w
        chois_by_order[index_order] = (chois, weights)
        # orthogonal supports: P_R P_S = 0
        PR = projector(vecs["R"])
        PS = projector(vecs["S"])
        prod = mat_mul(PR, PS, 12)
        scan[f"orthogonal-supports-{index_order}"] = all(
            prod[i][j].is_zero() for i in range(12) for j in range(12))
        for wname in ("A", "B"):
            QF = q_to_F(wits[wname]["Q"])
            QG = partial_transpose_B(QF, index_order)
            for sub in ("R", "Sbar", "S"):
                Pm = [[chois[sub][i][j] - QG[i][j] for j in range(12)]
                      for i in range(12)]
                Pm = shift_diag(Pm, -delta)
                okP, whyP = ldl_is_pd(Pm, 12)
                scan[f"{index_order}|{wname}|{sub}"] = {"P_pd": okP,
                                                        "detail": whyP}
                if okP:
                    winners.setdefault(index_order, {}).setdefault(
                        wname, []).append(sub)
    report["scan"] = scan

    complete = {}
    for io, m in winners.items():
        if set(m) == {"A", "B"} and all(len(v) == 1 for v in m.values()):
            complete[io] = {k: v[0] for k, v in m.items()}
    report["winning_assignments"] = {
        io: m for io, m in winners.items()}
    if not complete:
        report["verdict"] = ("REJECT: no convention certifies both channel "
                             "floors with a unique witness assignment")
        emit(report, report_path)
        return report

    # sigma: for each complete convention, find the trial input reproducing
    # the claimed exact spectrum
    accepted = None
    sigma_info = {}
    for io, assign in complete.items():
        if assign["A"] == assign["B"]:
            sigma_info[io] = "witnesses map to the same subspace"
            continue
        chois, weights = chois_by_order[io]
        sub1, sub2 = assign["A"], assign["B"]
        for trial in ("weighted", "maxent"):
            sig = joint_output_sigma(chois[sub1], chois[sub2],
                                     weights[sub1], weights[sub2], io, trial)
            rational = all(sig[i][j].is_rational()
                           for i in range(9) for j in range(9))
            hermit = is_hermitian(sig, 9)
            key = f"{io}|{trial}"
            row = {"entries_rational": rational, "hermitian": hermit}
            if rational and hermit:
                M9 = [[sig[i][j].rational() for j in range(9)]
                      for i in range(9)]
                tr = sum(M9[i][i] for i in range(9))
                row["trace"] = str(tr)
                cp = charpoly_frac(M9)
                row["charpoly_matches_claimed_spectrum"] = (
                    cp == spectrum_charpoly(spectrum))
                if tr == 1 and row["charpoly_matches_claimed_spectrum"]:
                    accepted = (io, assign, trial, M9)
            sigma_info[key] = row
            if accepted:
                break
        if accepted:
            break
    report["sigma_scan"] = sigma_info
    if accepted is None:
        report["verdict"] = ("REJECT: no convention/trial input reproduces "
                             "the claimed exact joint spectrum")
        emit(report, report_path)
        return report
    io, assign, trial, M9 = accepted
    report["accepted_convention"] = {
        "index_order": io, "witness_assignment": assign,
        "trial_input": trial,
        "note": ("transpose side of the channel map is unobservable here: "
                 "the trial inputs are real diagonal, and the floor "
                 "certificate covers all product inputs symmetrically")}

    # participant margin diagnostic: P - (mu - J_err) I still PD
    chois, weights = chois_by_order[io]
    margin = {}
    for wname, sub in assign.items():
        QF = q_to_F(wits[wname]["Q"])
        QG = partial_transpose_B(QF, io)
        m = wits[wname]["mu"] - wits[wname]["J_err"]
        Pm = [[chois[sub][i][j] - QG[i][j] for j in range(12)]
              for i in range(12)]
        Pm = shift_diag(Pm, -delta - m)
        okm, _ = ldl_is_pd(Pm, 12)
        margin[wname] = {"claimed_floor": str(m), "P_minus_floor_pd": okm}
    report["participant_margin_diagnostic"] = margin

    # diagnostic: compare with the participant's printed sigma matrix
    try:
        with open(sigma_claim_path) as fh:
            sc = json.load(fh)
        their = [[Fr(x) for x in row] for row in sc["matrix"]]
        report["sigma_equals_participant_matrix"] = (their == M9)
    except Exception as exc:
        report["sigma_equals_participant_matrix"] = f"unavailable: {exc}"

    tr_spec = sum(ev * mult for ev, mult in spectrum)
    nneg = all(ev >= 0 for ev, _ in spectrum)
    report["spectrum_checks"] = {
        "trace_one": tr_spec == 1,
        "all_nonnegative": nneg,
        "rank": sum(mult for ev, mult in spectrum if ev != 0),
    }
    if not (tr_spec == 1 and nneg):
        report["verdict"] = "REJECT: claimed spectrum is not a state spectrum"
        emit(report, report_path)
        return report

    status, cmp_info = strict_compare(delta, spectrum, p, bits)
    report["strict_comparison"] = cmp_info
    if status == "FALSE":
        report["verdict"] = ("REJECT: trace-power comparison is false at "
                             f"p={p} (B >= A^2)")
        emit(report, report_path)
        return report
    if status == "ABSTAIN-PRECISION":
        report["verdict"] = ("ABSTAIN: enclosure precision insufficient to "
                             "certify the strict comparison")
        emit(report, report_path)
        return report

    report["verdict"] = "CONFIRMED"
    report["established"] = (
        "For the exact CPTP pair built from the printed R/S subspaces of "
        "arXiv:0712.3628v3 (unnormalized-Choi convention Tr_B J = I_4, "
        "inverse-input-marginal normalization), every output of either "
        f"channel has all eigenvalues >= {delta}; the exact joint output of "
        "the accepted entangled input is rational with spectrum "
        "{0, 5/46 x2, 18/161 x2, 51/322 x2, 39/322 x2}; and "
        "A_lower^2 > B_upper certifies "
        f"S_p^min(N1 x N2) < S_p^min(N1) + S_p^min(N2) at p = {p} "
        "(natural log).")
    emit(report, report_path)
    return report

def emit(report, path):
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
        fh.write("\n")

# ----------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------

def self_test():
    a = F(tuple(Fr(x, y) for x, y in
                [(1, 3), (-2, 7), (0, 1), (1, 2), (2, 9), (0, 1), (1, 5), (0, 1)]),
          tuple(Fr(x, y) for x, y in
                [(1, 5), (0, 1), (-1, 3), (2, 1), (0, 1), (1, 8), (0, 1), (1, 7)]))
    b = F(tuple(Fr(x, y) for x, y in
                [(-1, 2), (1, 4), (1, 9), (0, 1), (0, 1), (2, 3), (0, 1), (1, 11)]),
          tuple(Fr(x, y) for x, y in
                [(0, 1), (3, 8), (1, 1), (-1, 6), (1, 2), (0, 1), (2, 7), (0, 1)]))
    fa, fb = a.to_complex(), b.to_complex()
    assert abs((a * b).to_complex() - fa * fb) < 1e-8
    assert abs((a + b).to_complex() - (fa + fb)) < 1e-10
    assert ((a * b) * (a * b).inv() - F_ONE).is_zero()
    assert (a.conj().conj() - a).is_zero()
    w = F.omega()
    assert (w * w * w - F_ONE).is_zero()
    assert (w * w - F.omega(2)).is_zero()
    for comps in [(1, 1, -1, 0, 0, 0, 0, 0), (0, -2, 1, 0, 0, 0, 0, 0),
                  (7, -5, 0, 0, 0, 0, 0, 0), (-10, 1, 2, 1, 0, 0, 0, 0),
                  (1, 0, 0, -1, 0, 0, 0, 0), (0, 3, 0, 0, 0, -1, 0, 0),
                  (2, 0, -1, 0, -1, 0, 1, 0), (-1, 1, 1, -1, 1, -1, 1, -1)]:
        x = tuple(Fr(c) for c in comps)
        fx = sum(float(c) * math.sqrt(r) for c, r in zip(x, RADICAND))
        assert sign8(x) == (fx > 0) - (fx < 0), (comps, fx)
    # LDL controls
    G5 = [[F.from_qi(Fr((i * 7 + j * 3) % 5 - 2, 3), Fr((i - j) % 3 - 1, 4))
           for j in range(5)] for i in range(5)]
    M = [[sum((G5[k][i].conj() * G5[k][j] for k in range(5)), F_ZERO)
          for j in range(5)] for i in range(5)]
    ok, _ = ldl_is_pd(shift_diag(M, Fr(1)), 5)
    assert ok
    ok2, _ = ldl_is_pd(shift_diag(M, Fr(-1000)), 5)
    assert not ok2
    # charpoly
    Mt = [[Fr(0)] * 9 for _ in range(9)]
    Mt[0][0], Mt[1][1] = Fr(1, 2), Fr(1, 3)
    want = spectrum_charpoly([(Fr(1, 2), 1), (Fr(1, 3), 1), (Fr(0), 7)])
    assert charpoly_frac(Mt) == want
    # roots
    l, u = rootP_bounds(Fr(1, 400), 100, 192)
    assert pow(l.numerator, 100) * 400 <= pow(l.denominator, 100)
    if u != l:
        assert pow(u.numerator, 100) * 400 > pow(u.denominator, 100)
    fl = float(Fr(1, 400)) ** 0.01
    assert abs(float(l) - fl) < 1e-12 and abs(float(u) - fl) < 1e-12
    lo, hi = ln_bounds(Fr(2), Fr(2))
    assert lo <= hi and abs(float(lo) - math.log(2)) < 1e-25
    # sqrt_rat coverage used by the run
    for x in (Fr(25, 9), Fr(9, 4), Fr(16, 9), Fr(5, 2), Fr(20, 9), Fr(2),
              Fr(5), Fr(400, 81), Fr(81, 16), Fr(1)):
        s = F.sqrt_rat(x)
        assert abs(s.to_complex().real - math.sqrt(float(x))) < 1e-10
    print("self-test: OK", flush=True)

if __name__ == "__main__":
    self_test()
    if len(sys.argv) < 3:
        print("usage: verify_q4_p001.py <bundle-dir> <report-out>")
        raise SystemExit(2)
    base, out = sys.argv[1], sys.argv[2]
    cand = f"{base}/progress/artifacts/renyi-p001-candidate.json"
    sig = f"{base}/progress/artifacts/joint_spectrum_exact.json"
    rep = run(cand, sig, out)
    print("verdict:", rep["verdict"], flush=True)

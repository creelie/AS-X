"""Theta-decomposition computation of the twined single-centred indices at P^2=2.
Series in t = q^{1/4} as flint fmpq_poly truncated at t^L."""
import flint, pickle, math
from fractions import Fraction as F
from psi import hurwitz, SHAPE, ORDER, dg
_fq = flint.fmpq_poly
def _cv(x):
    if isinstance(x, F): return flint.fmpq(x.numerator, x.denominator)
    return x
def fq(c):
    if isinstance(c, list): c = [_cv(x) for x in c]
    return _fq(c)

def trunc(p, L):
    c = p.coeffs()[:L]
    return fq(c) if c else fq([0])

def mulT(a, b, L):
    return trunc(a * b, L)

def from_q(coeffs_q, L):
    """q-series list -> t-series (t^4 = q)"""
    c = [0] * L
    for n, x in enumerate(coeffs_q):
        if 4 * n < L: c[4 * n] = x
    return fq(c)

def etaprod_inv_q(shape, Nq):
    """1/prod (1-q^{an})^{m_a} as q-list of length Nq (ints)"""
    r = [0] * Nq; r[0] = 1
    for a, m in shape.items():
        for n in range(1, Nq):
            if a * n >= Nq: break
            for _ in range(m):
                for i in range(a * n, Nq):
                    r[i] += r[i - a * n]
    return r

def sig(n, k):
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d ** k
            if d * d != n: s += (n // d) ** k
        d += 1
    return s

def Ek(k, d, Nq):
    c = {2: -24, 4: 240, 6: -504}[k]
    r = [0] * Nq; r[0] = 1
    for n in range(1, Nq):
        if n % d == 0: r[n] = c * sig(n // d, k - 1)
    return r

def phiN_q(N, Nq):
    a = Ek(2, N, Nq); b = Ek(2, 1, Nq)
    return [F(N * x - y, N - 1) for x, y in zip(a, b)]

def thetas(L):
    th0 = [0] * L; th1 = [0] * L; d0 = [0] * L; d1 = [0] * L
    n = 0
    while n * n < L:
        for s in ((n, -n) if n else (0,)):
            e = s * s
            if s % 2 == 0:
                th0[e] += 1; d0[e] += F(e, 4)
            else:
                th1[e] += 1; d1[e] += F(e, 4)
        n += 1
    return fq(th0), fq(th1), fq(d0), fq(d1)

class Engine:
    def __init__(self, Dmax):
        self.L = L = Dmax + 8
        self.Nq = Nq = L // 4 + 4
        th0, th1, d0, d1 = thetas(L)
        # 1/eta^6 = t^{-1} / (q)^6 ; keep shift separately: store inv6 = 1/(q)^6
        inv6 = from_q(etaprod_inv_q({1: 6}, Nq), L)
        E2 = from_q(Ek(2, 1, Nq), L)
        # components without the t^{-1} shift (shift applied at end)
        self.hA = [mulT(-th1, inv6, L), mulT(th0, inv6, L)]
        self.hB = [mulT(24 * d1 - mulT(E2, th1, L), inv6, L), mulT(-24 * d0 + mulT(E2, th0, L), inv6, L)]
        # H components: H0 = sum H(4n) q^n -> exponent 4n=D ; H1 = sum H(4n+3) q^{n+3/4} -> exponent D
        H0 = [0] * L; H1 = [0] * L
        for D in range(0, L):
            if D % 4 == 0: H0[D] = hurwitz(D)
            elif D % 4 == 3: H1[D] = hurwitz(D)
        self.H = [fq(H0), fq(H1)]
        self.cache = {}

    def mod_form(self, name):
        kind, d = name
        Nq = self.Nq; L = self.L
        if kind == 'E4': return from_q(Ek(4, d, Nq), L)
        if kind == 'phi': return from_q(phiN_q(d, Nq), L)
        if kind == 'P':
            a, b = d
            return mulT(from_q(phiN_q(a, Nq), L), from_q(phiN_q(b, Nq), L), L)
        raise ValueError

    def C(self, cls, Wfit):
        """returns dict D -> C_g(D) for -1 <= D < Dmax (D = 0,3 mod 4)"""
        terms, d2 = Wfit[cls]
        L = self.L
        comp = []
        for j in (0, 1):
            num = fq([0])
            for (kind, d), c in terms:
                if kind == 'B':
                    num += mulT(self.mod_form(('phi', d)), self.hB[j], L) * c
                else:
                    num += mulT(self.mod_form((kind, d)), self.hA[j], L) * c
            # A,B parts carry t^{-1}; H part does not: num_total = t^{-1} num - 2 d2 H_j
            # divide by eta_g = q * prod(...) = t^4 prod: psi = t^{-4} inv_g * (t^{-1} num - 2 d2 H)
            inv_g = from_q(etaprod_inv_q(SHAPE[cls], self.Nq), L)
            a = mulT(num, inv_g, L)          # exponent e -> actual e - 5
            b = mulT(self.H[j], inv_g, L)    # exponent e -> actual e - 4
            comp.append((a, b, d2))
        out = {}
        for D in range(-1, L - 8):
            j = 0 if D % 4 == 0 else (1 if D % 4 == 3 else None)
            if j is None: continue
            a, b, d2 = comp[j]
            ac = a.coeffs(); bc = b.coeffs()
            va = ac[D + 5] if 0 <= D + 5 < len(ac) else 0
            vb = bc[D + 4] if 0 <= D + 4 < len(bc) else 0
            out[D] = va - 2 * d2 * vb
        return out

if __name__ == '__main__':
    import sys
    Wfit = pickle.load(open('Wfit.pkl', 'rb'))
    Wfit = {k: ([((n[0], n[1]), flint.fmpq(int(c.p), int(c.q))) for n, c in v[0]], v[1]) for k, v in Wfit.items()}
    E = Engine(60)
    for cls in ['1A', '2A', '3A', '4B']:
        C = E.C(cls, Wfit)
        print(cls, [(D, (-1) ** (D + 1) * C[D]) for D in [3, 4, 7, 8, 11, 12]])

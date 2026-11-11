"""Verify W_g = f_g A + h_g B up to the Sturm bound for Gamma(4N), via Taylor moments."""
import flint, sys, math
from fractions import Fraction as Fr
from fast1d import fq, mulT, from_q, etaprod_inv_q, Ek, phiN_q, thetas, trunc
from psi import TG, SQ, SHAPE, ORDER, dg
from scan import Wfit, CLS

def index1_C(hcomps, L):
    """hcomps: (h0, h1) t-series WITH the t^{-1} already included as a shift s: we pass coefficient lists
    where entry e corresponds to t^{e-1}. Return dict D->C(D)."""
    h0, h1 = hcomps
    C = {}
    for j, h in ((0, h0), (1, h1)):
        c = h.coeffs()
        for e, x in enumerate(c):
            D = e - 1
            if x != 0:
                assert (D % 4) == (0 if j == 0 else 3), (j, D)
                C[D] = Fr(int(x.p), int(x.q))
    return C

def moments_from_C(C, nmax, numax):
    """m_nu(n) = sum_l C(4n-l^2) l^{2nu}"""
    out = [[Fr(0)] * (nmax + 1) for _ in range(numax + 1)]
    for n in range(nmax + 1):
        L = math.isqrt(4 * n + 1)
        for l in range(-L, L + 1):
            D = 4 * n - l * l
            v = C.get(D, 0)
            if v:
                for nu in range(numax + 1):
                    out[nu][n] += v * l ** (2 * nu)
    return out

def pmul(a, b, n):
    r = [Fr(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j in range(0, n + 1 - i):
                if b[j]: r[i + j] += x * b[j]
    return r

def mprod(F, G, n, numax):
    out = []
    for nu in range(numax + 1):
        s = [Fr(0)] * (n + 1)
        for i in range(nu + 1):
            t = pmul(F[i], G[nu - i], n)
            c = math.comb(2 * nu, 2 * i)
            s = [x + c * y for x, y in zip(s, t)]
        out.append(s)
    return out

def run(cls, Bq):
    nmax = 2 * Bq + 4
    L = 4 * nmax + 12
    Nq = L // 4 + 4
    th0, th1, d0, d1 = thetas(L)
    inv6 = from_q(etaprod_inv_q({1: 6}, Nq), L)
    E2 = from_q(Ek(2, 1, Nq), L)
    hA = (mulT(-th1, inv6, L), mulT(th0, inv6, L))
    hB = (mulT(24 * d1 - mulT(E2, th1, L), inv6, L), mulT(-24 * d0 + mulT(E2, th0, L), inv6, L))
    CA = index1_C(hA, L); CB = index1_C(hB, L)
    def Cg(c):
        chi, coef = TG[c]
        T = [Fr(0)] * Nq
        for d, cd in coef.items():
            p = phiN_q(d, Nq)
            T = [x + cd * y for x, y in zip(T, p)]
        Tt = from_q(T, L)
        comps = (mulT(Tt, hA[0], L) + hB[0] * flint.fmpq(chi, 12), mulT(Tt, hA[1], L) + hB[1] * flint.fmpq(chi, 12))
        return index1_C(comps, L)
    numax = 2
    MA = moments_from_C(CA, nmax, numax); MB = moments_from_C(CB, nmax, numax)
    Mg = moments_from_C(Cg(cls), nmax, numax); Mg2 = moments_from_C(Cg(SQ[cls]), nmax, numax)
    n = Bq + 2
    # X = 1/2 phi_g^2 + U2 phi_g + 1/2 D2 phi_{g^2}
    sq = mprod(Mg, Mg, n, numax)
    U2 = [[Mg[nu][2 * k] for k in range(n + 1)] for nu in range(numax + 1)]
    D2 = [[(4 ** nu) * Mg2[nu][k // 2] if k % 2 == 0 else Fr(0) for k in range(n + 1)] for nu in range(numax + 1)]
    X = [[sq[nu][k] / 2 + U2[nu][k] + D2[nu][k] / 2 for k in range(n + 1)] for nu in range(numax + 1)]
    d2 = dg(SHAPE[cls], 2)
    E4 = [Fr(x) for x in Ek(4, 1, n + 1)]
    BB = mprod(MB, MB, n, numax); AA = mprod(MA, MA, n, numax)
    E4AA = [pmul(E4, AA[nu], n) for nu in range(numax + 1)]
    Nm = [[X[nu][k] - Fr(d2, 144) * (BB[nu][k] - E4AA[nu][k]) for k in range(n + 1)] for nu in range(numax + 1)]
    assert all(x == 0 for x in Nm[0]), 'numerator does not vanish at z=0'
    m0W = [x / 2 for x in Nm[1]]
    m2A = MA[2][:n + 1]
    m1W = [(Nm[2][k] - sum(m0W[i] * m2A[k - i] for i in range(k + 1))) / 12 for k in range(n + 1)]
    h = [x / 12 for x in m0W]
    E2l = [Fr(x) for x in Ek(2, 1, n + 1)]
    f = [(m1W[k] - 2 * sum(E2l[i] * h[k - i] for i in range(k + 1))) / 2 for k in range(n + 1)]
    # fitted
    terms, _ = Wfit[cls]
    ff = [Fr(0)] * (n + 1); hf = [Fr(0)] * (n + 1)
    for (kind, d), c in terms:
        c = Fr(int(c.p), int(c.q))
        if kind == 'E4': s = [Fr(x) for x in Ek(4, d, n + 1)]; tgt = ff
        elif kind == 'P':
            a, b = d; s = pmul(phiN_q(a, n + 1), phiN_q(b, n + 1), n); tgt = ff
        else: s = phiN_q(d, n + 1); tgt = hf
        for k in range(n + 1): tgt[k] += c * s[k]
    okf = all(f[k] == ff[k] for k in range(Bq + 1)); okh = all(h[k] == hf[k] for k in range(Bq + 1))
    return okf, okh

def sturm_bound(N, k=4):
    M = 4 * N
    prod = 1.0
    for p in range(2, M + 1):
        if M % p == 0 and all(p % r for r in range(2, p)): prod *= (1 - 1 / p ** 2)
    return math.floor(k * M * M * prod / 24) + 1

if __name__ == '__main__':
    for cls in ['1A', '2A', '3A', '4B', '5A', '6A', '7A', '8A']:
        N = ORDER[cls]
        Bq = max(sturm_bound(N), 8)
        print(cls, 'Sturm bound (q-units)', Bq, run(cls, Bq), flush=True)

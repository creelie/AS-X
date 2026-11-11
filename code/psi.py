from jac import *
from fractions import Fraction as F
import math, sys

SHAPE = {'1A':{1:24},'2A':{1:8,2:8},'3A':{1:6,3:6},'4B':{1:4,2:2,4:4},'5A':{1:4,5:4},
         '6A':{1:2,2:2,3:2,6:2},'7A':{1:3,7:3},'8A':{1:2,2:1,4:1,8:2}}
TG = {'1A':(24,{}),'2A':(8,{2:F(4,3)}),'3A':(6,{3:F(3,2)}),'4B':(4,{2:F(-1,3),4:F(2)}),'5A':(4,{5:F(5,3)}),
      '6A':(2,{2:F(-1,6),3:F(-1,2),6:F(5,2)}),'7A':(3,{7:F(7,4)}),'8A':(2,{4:F(-1,2),8:F(7,3)})}
SQ = {'1A':'1A','2A':'1A','3A':'3A','4B':'2A','5A':'5A','6A':'3A','7A':'7A','8A':'4B'}
ORDER = {'1A':1,'2A':2,'3A':3,'4B':4,'5A':5,'6A':6,'7A':7,'8A':8}

def hurwitz(N):
    if N == 0: return F(-1,12)
    if N % 4 in (1,2): return F(0)
    h = F(0)
    a = 1
    while 3*a*a <= N:
        for b in range(-a+1, a+1):
            if (b*b + N) % (4*a) == 0:
                c = (b*b + N)//(4*a)
                if c < a: continue
                if c == a and b < 0: continue
                w = F(1)
                if b == 0 and a == c: w = F(1,2)
                elif b == a and a == c: w = F(1,3)
                h += w
        a += 1
    return h

def Hjac(Q):
    d = {}
    for n in range(Q):
        L = int(math.isqrt(4*n))
        for l in range(-L, L+1):
            v = hurwitz(4*n - l*l)
            if v: d[(n,l)] = v
    return d

def dg(shape, k):
    r = eta_prod({a:-m for a,m in shape.items()}, k+1)
    return r[k]

def phi_g(cls, Q, A, B):
    chi, c = TG[cls]
    T = [0]*Q
    for d, cd in c.items():
        p = phiN(d, Q)
        T = [x + cd*y for x, y in zip(T, p)]
    return add(scal(B, F(chi,12)), mul(qser(T), A, Q))

def U2(phi, Q):
    return {(n//2, l): c for (n, l), c in phi.items() if n % 2 == 0 and n//2 < Q}

def V2dil(phi, Q):
    return {(2*n, 2*l): c for (n, l), c in phi.items() if 2*n < Q}

def divA(Fs, Q):
    """exact W with W*phi_{-2,1} = Fs (truncate Q)."""
    # multiply by prod (1-q^n)^4 / ((1-q^n y)^2 (1-q^n/y)^2)
    s = mul(Fs, qser(eta_prod({1:4}, Q)), Q)
    for n in range(1, Q):
        geo_p = {(n*k, k): 1 for k in range(0, Q//n + 1) if n*k < Q}
        geo_m = {(n*k, -k): 1 for k in range(0, Q//n + 1) if n*k < Q}
        for _ in range(2):
            s = mul(s, geo_p, Q); s = mul(s, geo_m, Q)
    # divide each q order by (y-1)^2/y
    out = {}
    for n in range(Q):
        poly = {l: c for (m, l), c in s.items() if m == n}
        if not poly: continue
        lo = min(poly); hi = max(poly)
        # poly * y = (y-1)^2 * R ; coefficients list from lo+1..hi+1
        coeffs = [poly.get(l, 0) for l in range(lo, hi+1)]  # degree index i <-> y^{lo+i}, times y -> y^{lo+1+i}
        # synthetic division by (y-1) twice, highest-first
        def divlin(cs):
            # cs low->high, divide by (y-1)
            hiq = cs[::-1]
            q = []; r = 0
            for x in hiq:
                r = r + x
                q.append(r)
            rem = q.pop()
            return q[::-1], rem
        q1, r1 = divlin(coeffs)
        q2, r2 = divlin(q1)
        assert r1 == 0 and r2 == 0, ('nonexact', n, r1, r2)
        # coeffs*y^{lo} * y = (y-1)^2 * q2*y^{lo}  => W part = q2 * y^{lo}... careful: poly = sum c_i y^{lo+i}
        # poly*y = (y-1)^2 * sum q2_i y^{lo+i}  -> W = sum q2_i y^{lo+i} ... wait degrees: poly*y has deg hi+1, (y-1)^2 * y^{lo+i}
        for i, c in enumerate(q2):
            if c: out[(n, lo + 1 + i)] = c
    return out

def psiF1(cls, Q):
    A = phim21(Q+2); B = phi01(Q+2)
    Qe = Q + 2
    pg = phi_g(cls, Qe, A, B); pg2 = phi_g(SQ[cls], Qe, A, B)
    X = add(add(scal(mul(pg, pg, Qe), F(1,2)), U2(pg, Qe)), scal(V2dil(pg2, Qe), F(1,2)))
    d2 = dg(SHAPE[cls], 2)
    E4s = qser(E4(Qe))
    Fs = add(X, scal(add(mul(B,B,Qe), mul(E4s, mul(A,A,Qe), Qe), -1), F(d2,144)), -1)
    W = divA(Fs, Qe)
    num = add(W, scal(Hjac(Qe), 2*d2), -1)
    inv = eta_prod({a:-m for a,m in SHAPE[cls].items()}, Qe)
    psi = mul(num, qser(inv), Qe)
    psi = {(n-1, l): c for (n, l), c in psi.items() if n-1 < Q}
    return psi, W, d2, X

if __name__ == '__main__':
    Q = 6
    for cls in ['1A','2A','3A','4B']:
        psi, W, d2, X = psiF1(cls, Q)
        Xz = {}
        for (n,l),c in X.items(): Xz[n] = Xz.get(n,0)+c
        print(cls, 'd2=',d2, 'X(tau,0)=', [Xz.get(n,0) for n in range(4)])
        for n in range(1,4):
            print('  n',n, [(-1)**(l+1)*psi.get((n,l),0) for l in range(0,3)])

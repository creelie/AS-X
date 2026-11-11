import math
from jac import *
from fractions import Fraction

FIX = {  # order: list of (alpha numerator k, multiplicity) with alpha=k/N
 '2A': (2, [(1,8)]), '3A': (3, [(1,6)]), '4B': (4, [(1,4)]), '5A': (5, [(1,2),(2,2)]),
 '6A': (6, [(1,2)]), '7A': (7, [(1,1),(2,1),(3,1)]), '8A': (8, [(1,1),(3,1)])}

def f_alpha(c, Q):
    # (y+1/y-c)/(2-c) * prod (1 - c q^n y + q^2n y^2)(1 - c q^n/y + q^2n/y^2)/(1 - c q^n + q^2n)^2
    s = {(0,1): 1/(2-c), (0,0): -c/(2-c), (0,-1): 1/(2-c)}
    for n in range(1, Q):
        s = mul(s, {(0,0):1.0, (n,1):-c, (2*n,2):1.0}, Q)
        s = mul(s, {(0,0):1.0, (n,-1):-c, (2*n,-2):1.0}, Q)
        den = [0.0]*Q; den[0]=1.0
        if n < Q: den[n] = -c
        if 2*n < Q: den[2*n] = 1.0
        inv = sinv(den, Q)
        inv2 = smul(inv, inv, Q)
        s = mul(s, qser(inv2), Q)
    return s

def twining_fp(cls, Q):
    N, data = FIX[cls]
    tot = {}
    for k, m in data:
        c = 2*math.cos(2*math.pi*k/N)
        tot = add(tot, scal(f_alpha(c, Q), m))
    return tot

if __name__ == '__main__':
    Q = 10
    A = phim21(Q); B = phi01(Q)
    for cls in FIX:
        N = FIX[cls][0]
        g = twining_fp(cls, Q)
        chi = sum(m for _, m in FIX[cls][1])
        rest = add(g, scal(B, Fraction(chi,12)), -1)
        # rest = T * A; T = sum c_d phiN(d); fit using coefficients (n, l=1) : A has (0,1)=1
        divs = [d for d in range(2, N+1) if N % d == 0]
        basis = [mul(qser(phiN(d, Q)), A, Q) for d in divs]
        import numpy as np
        keys = sorted(set(k for b in basis for k in b) | set(rest))
        M = np.array([[float(b.get(k,0)) for b in basis] for k in keys])
        v = np.array([float(rest.get(k,0)) for k in keys])
        sol, res, *_ = np.linalg.lstsq(M, v, rcond=None)
        print(cls, chi, dict(zip(divs, [round(x,6) for x in sol])), 'resid', float(np.max(np.abs(M@sol-v))))
        print('   q^0:', {k:round(v,4) for k,v in g.items() if k[0]==0}, ' q^1 y^0..2:', [round(g.get((1,l),0),4) for l in range(3)])

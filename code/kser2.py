"""K'_j(t) = 3E4 theta_{1-j}/(q)^6 +/- 648 t H_j - Rmaj_j - t*Hsum*|H_j|  (all over 1/(q)^24).
Then d* - sum_b |d*_{g^b}| >= [t^D] t^{-5} K'_j/(q)^24."""
import pickle, flint, sys
from fast1d import *
from scan import Wfit, CLS
from math import gcd
from fractions import Fraction as Fr

def build(N, Dmax):
    E = Engine(Dmax); L = E.L; Nq = E.Nq
    th0, th1, d0, d1 = thetas(L)
    E2 = from_q(Ek(2,1,Nq), L); E4 = from_q(Ek(4,1,Nq), L)
    inv6 = from_q(etaprod_inv_q({1:6}, Nq), L)
    numA = [-th1, th0]; numB = [24*d1 - mulT(E2, th1, L), -24*d0 + mulT(E2, th0, L)]
    out = []
    for j in (0,1):
        main = 3*mulT(mulT(E4, [th1, th0][j], L), inv6, L)
        R = [0]*L; Hsum = 0
        for b in range(1, N):
            cls = CLS[N//gcd(b,N)]
            terms, d2 = Wfit[cls]
            s = fq([0])
            for (kind, d), c in terms:
                if kind == 'B': s += mulT(E.mod_form(('phi', d)), numB[j], L)*c
                else: s += mulT(E.mod_form((kind, d)), numA[j], L)*c
            for i, x in enumerate(s.coeffs()): R[i] += abs(x)
            Hsum += 2*d2
        Hj = E.H[j]
        sign = 1 if j == 0 else -1
        Hc = Hj.coeffs()
        Hpart = [0] + [sign*648*x - Hsum*abs(x) for x in Hc]
        K = main - fq(R) + trunc(fq(Hpart), L)
        c = K.coeffs()
        # q-series k_j: exponents e = 4m + eps (eps=1 for j=0, 0 for j=1)
        eps = 1 if j == 0 else 0
        k = [Fr(int(c[4*m+eps].p), int(c[4*m+eps].q)) if 4*m+eps < len(c) else 0 for m in range((L-1)//4)]
        # sanity: other exponents zero
        assert all(c[e] == 0 for e in range(len(c)) if e % 4 != eps), 'stray exponents'
        out.append((k, Hsum))
    return out

SIG={}
def sigma1(n):
    if n not in SIG: SIG[n]=sum(d for d in range(1,n+1) if n%d==0)
    return SIG[n]

if __name__ == '__main__':
    Dmax = int(sys.argv[1]); kk = int(sys.argv[2]) if len(sys.argv)>2 else 24
    for N in range(2, 9):
        for j, (k, Hsum) in enumerate(build(N, Dmax)):
            a0 = k[0]
            worst = min((k[m]*max(m,2)/(kk*abs(a0)*sigma1(m)) if a0 < 0 else 10**9, m) for m in range(1, len(k)-2))
            neg = [m for m in range(1, len(k)-2) if k[m] <= 0]
            print(CLS[N], j, 'a0=', a0, 'a1..3=', [float(x) for x in k[1:4]], 'neg m>=1:', neg[:5], 'min ratio a_m/bound', (float(worst[0]), worst[1]))

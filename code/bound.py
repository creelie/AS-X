import math
from fractions import Fraction as Fr
from scan import Wfit, CLS
from math import gcd
A0 = {}  # filled from kser2 run
def consts(N):
    cA = cB = 0; Hs = 0
    for b in range(1, N):
        terms, d2 = Wfit[CLS[N//gcd(b,N)]]
        for (kind, d), c in terms:
            c = abs(Fr(int(c.p), int(c.q)))
            if kind == 'B': cB += c
            else: cA += c
        Hs += 2*d2
    return cA, cB, Hs
def M0(N, j, a0):
    cA, cB, Hs = consts(N)
    mu = 2 if j == 0 else 1
    m = 1
    last_bad = 0
    for m in range(1, 200000):
        lhs = 3*mu*(math.comb(m+5,5) + 240*m**3)
        rhs = (2*math.sqrt(m)+2)*(337*float(cA+cB)*m**3 + 576*float(cB)*(m+0.25)*m*(1+math.log(m))) + 4*m*(648*(j==1)+Hs) + 24*abs(a0)*(1+math.log(m))
        if lhs <= rhs: last_bad = m
        elif m > 10*last_bad+1000: break
    return last_bad+1, (float(cA), float(cB), Hs)
if __name__ == '__main__':
    a0 = {(2,0):64,(2,1):0,(3,0):72,(3,1):3,(4,0):80,(4,1):6,(5,0):80,(5,1):9,(6,0):96,(6,1):12,(7,0):84,(7,1):15,(8,0):96,(8,1):18}
    for N in range(2,9):
        for j in (0,1):
            print(CLS[N], j, M0(N,j,a0[(N,j)]))

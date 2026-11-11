import math
from bound import consts
from scan import CLS
a0 = {(2,0):64,(2,1):0,(3,0):72,(3,1):3,(4,0):80,(4,1):6,(5,0):80,(5,1):9,(6,0):96,(6,1):12,(7,0):84,(7,1):15,(8,0):96,(8,1):18}
def ok(m,N,j):
    cA,cB,Hs=consts(N); cA=float(cA); cB=float(cB); mu=2 if j==0 else 1
    lhs=3*mu*m**5/120
    rhs=(2*math.sqrt(m)+2)*(289*(cA+cB)*m**3+576*cB*(m+0.25)*m*(1+math.log(m)))+4*m*(648*(j==1)+Hs)+24*a0[(N,j)]*(1+math.log(m))
    return lhs > rhs*(1+1e-9)
for N in range(2,9):
    for j in (0,1):
        lo,hi=1,10**6
        while hi-lo>1:
            mid=(lo+hi)//2
            if ok(mid,N,j): hi=mid
            else: lo=mid
        print(CLS[N],j,hi, [round(float(x),3) for x in consts(N)])

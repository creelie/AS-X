import pickle, flint, math, sys, time
from fast1d import Engine
from math import gcd
from fractions import Fraction as F
Wfit = pickle.load(open('Wfit.pkl','rb'))
Wfit = {k: ([((n[0], n[1]), flint.fmpq(int(c.p), int(c.q))) for n, c in v[0]], v[1]) for k, v in Wfit.items()}
def mu(n):
    r=1; p=2
    while p*p<=n:
        if n%p==0:
            n//=p
            if n%p==0: return 0
            r=-r
        p+=1
    return -r if n>1 else r
def ram(q,a):
    g=gcd(q,a)
    return sum(mu(q//d)*d for d in range(1,g+1) if g%d==0)
CLS = {1:'1A',2:'2A',3:'3A',4:'4B',5:'5A',6:'6A',7:'7A',8:'8A'}
if __name__=='__main__':
    Dmax = int(sys.argv[1])
    t=time.time(); E = Engine(Dmax)
    C = {cls: E.C(cls, Wfit) for cls in CLS.values()}
    print('time', time.time()-t)
    pickle.dump({k:{D:(int(v.p),int(v.q)) for D,v in c.items()} for k,c in C.items()}, open(f'C_{Dmax}.pkl','wb'))
    for N in range(2,9):
        worst=(0,None); allpos=True; allint=True; minS=None
        divs=[e for e in range(1,N+1) if N%e==0]
        for D in range(3, Dmax):
            if D%4 not in (0,3): continue
            s=(-1)**(D+1)
            dd={e: F(int(C[CLS[N//e]][D].p),int(C[CLS[N//e]][D].q))*s for e in divs}  # g^e has order N/e
            # number of b in [1,N-1] with gcd(b,N)=e is phi(N/e)
            r = sum(abs(dd[e])*sum(1 for b in range(1,N) if gcd(b,N)==e) for e in divs if e<N)/dd[N]
            if r>worst[0]: worst=(float(r),D)
            for a in range(N):
                Sa = sum(dd[e]*ram(N//e,a) for e in divs)/N
                if Sa<=0: allpos=False
                if Sa.denominator!=1: allint=False
                rel = Sa/dd[N]
                if minS is None or rel<minS[0]: minS=(float(rel),D,a)
        print(CLS[N],'worst sum|d|/d*',worst,'all S_a>0',allpos,'integral',allint,'min S_a/d*',minS)

import pickle, math
from fractions import Fraction as F
from math import gcd
from kser2 import build, sigma1
from scan import CLS, ram
C=pickle.load(open('C_600.pkl','rb'))
out={'ratio':{}, 'S':{}, 'marg':{}}
for N in range(2,9):
    divs=[e for e in range(1,N+1) if N%e==0]
    rs=[]; Ss=[]
    for D in range(3,600):
        if D%4 not in (0,3): continue
        s=(-1)**(D+1)
        dd={e: F(*C[CLS[N//e]][D])*s for e in divs}
        r=sum(abs(dd[e])*sum(1 for b in range(1,N) if gcd(b,N)==e) for e in divs if e<N)/dd[N]
        rs.append((D,float(r)))
        Ss.append((D,[float(sum(dd[e]*ram(N//e,a) for e in divs)/N/dd[N]) for a in range(N)]))
    out['ratio'][CLS[N]]=rs; out['S'][CLS[N]]=Ss
    for j,(k,Hs) in enumerate(build(N,8000)):
        a0=k[0]
        out['marg'][(CLS[N],j)]=[(m, float(k[m]*max(m,2)/(24*abs(a0)*sigma1(m)))) if a0<0 else (m,float(k[m]/(24*sigma1(m)))) for m in range(1,len(k)-2)]
pickle.dump(out,open('figdata.pkl','wb'))

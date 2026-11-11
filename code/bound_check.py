"""Corollary (large m): check mu m^5/40 > rho(m) with 40-digit arithmetic, at m0 of the table and at m = 11000.
The right side divided by m^5 is decreasing, so this proves condition (cond) for all m >= 11000."""
import mpmath as mp
from bound import consts
from bound2 import a0
from scan import CLS
mp.mp.dps = 40
M0 = {(2,0):1089,(2,1):1720,(3,0):1720,(3,1):2719,(4,0):2483,(4,1):3929,(5,0):4914,(5,1):7784,
      (6,0):4005,(6,1):6343,(7,0):6810,(7,1):10792,(8,0):4645,(8,1):7356}
def ratio(m, N, j):
    cA, cB, Hs = [mp.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else mp.mpf(x) for x in consts(N)]
    mu = 2 if j == 0 else 1; m = mp.mpf(m)
    rho = (2*mp.sqrt(m)+2)*(289*(cA+cB)*m**3 + 576*cB*(m+mp.mpf(1)/4)*m*(1+mp.log(m))) \
          + 4*m*(648*(j == 1) + Hs) + 24*a0[(N, j)]*(1+mp.log(m))
    return (mu*m**5/40)/rho
for N in range(2, 9):
    for j in (0, 1):
        print(CLS[N], j, 'ratio at m0', mp.nstr(ratio(M0[(N,j)], N, j), 8), 'at m0-1', mp.nstr(ratio(M0[(N,j)]-1, N, j), 8),
              'at 11000', mp.nstr(ratio(11000, N, j), 8))

"""Check (B^2 - E4 A^2) = 144 A (A_{2,1} - 2 H) coefficientwise, with A_{2,1} expanded in the
strip -Im tau < Im z < 0 (|q| < |y|^{-1} < 1), for q-exponents n < Q and |l| <= L - margin."""
import sys
from fractions import Fraction as F
from jac import mul, add, scal, qser, phim21, phi01, E4
from psi import hurwitz

Q = int(sys.argv[1]) if len(sys.argv) > 1 else 40
L = 3 * Q + 20
A = phim21(Q); B = phi01(Q)
lhs = add(mul(B, B, Q), mul(mul(qser(E4(Q)), A, Q), A, Q), -1)
# A_{2,1} in strip S_0
A21 = {}
for s in range(-Q, Q):
    for k in range(0, 4 * Q + 4 * L):
        if s >= 1:
            n, l = s * s + s + k * s, 2 * s + 1 + k
        else:
            n, l = s * s - s * (k + 1), 2 * s - k - 1
        if n >= Q or abs(l) > L:
            if s >= 1 and n >= Q: break
            if s <= 0 and (abs(l) > L and l < 0): break
            continue
        A21[(n, l)] = A21.get((n, l), 0) + (k + 1)
H = {}
for n in range(Q):
    for l in range(-L, L + 1):
        D = 4 * n - l * l
        if D == 0: H[(n, l)] = F(-1, 12)
        elif D > 0:
            h = hurwitz(D)
            if h: H[(n, l)] = h
rhs = scal(mul(A, add(A21, H, -2), Q), 144)
bad = [k for k in set(lhs) | set(rhs) if abs(k[1]) <= L - 2 * Q - 5 and lhs.get(k, 0) != rhs.get(k, 0)]
print('Q =', Q, 'compared |l| <=', L - 2 * Q - 5, 'mismatches:', sorted(bad)[:10])

"""Exact check of Proposition (twining genus): the fixed-point formula for phi_g equals
chi_g/12 phi_{0,1} + T_g phi_{-2,1}, coefficient by coefficient up to q^Q, computed in the
cyclotomic field Q(zeta_N) with exact rational arithmetic (sympy)."""
import sympy as sp
from fractions import Fraction as F
from jac import mul, add, scal, qser, phim21, phi01, phiN
from twin import FIX
from psi import TG

x = sp.symbols('x')

def check(cls, Q):
    N, data = FIX[cls]
    K = sp.QQ.algebraic_field(sp.exp(2 * sp.pi * sp.I / N))
    z = K.from_sympy(sp.exp(2 * sp.pi * sp.I / N))
    one = K.one
    tot = {}
    for k, mult in data:
        zk = z ** k
        c = zk + zk ** (N - 1)
        inv2c = one / (2 * one - c)
        s = {(0, 1): inv2c, (0, 0): -c * inv2c, (0, -1): inv2c}
        for n in range(1, Q):
            s = mul(s, {(0, 0): one, (n, 1): -c, (2 * n, 2): one}, Q)
            s = mul(s, {(0, 0): one, (n, -1): -c, (2 * n, -2): one}, Q)
            # 1/(1 - c q^n + q^{2n})^2 as a power series in q
            den = [K.zero] * Q; den[0] = one
            den[n] = -c
            if 2 * n < Q: den[2 * n] = one
            inv = [K.zero] * Q; inv[0] = one
            for j in range(1, Q):
                inv[j] = -sum((den[i] * inv[j - i] for i in range(1, j + 1)), K.zero)
            inv2 = [sum((inv[i] * inv[j - i] for i in range(j + 1)), K.zero) for j in range(Q)]
            s = mul(s, {(j, 0): inv2[j] for j in range(Q) if inv2[j] != K.zero}, Q)
        for key, v in s.items():
            tot[key] = tot.get(key, K.zero) + mult * v
    chi, T = TG[cls]
    A = phim21(Q); B = phi01(Q)
    Tq = [F(0)] * Q
    for d, cd in T.items():
        for i, v in enumerate(phiN(d, Q)): Tq[i] += cd * v
    rhs = add(scal(B, F(chi, 12)), mul(qser(Tq), A, Q))
    keys = set(tot) | set(rhs)
    bad = []
    for key in keys:
        a = tot.get(key, K.zero)
        r = rhs.get(key, 0)
        r = K.from_sympy(sp.Rational(r.numerator, r.denominator)) if isinstance(r, F) else K.from_sympy(sp.Integer(r))
        if a != r: bad.append(key)
    return chi, len(keys), bad

if __name__ == '__main__':
    import sys
    Q = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for cls in FIX:
        chi, n, bad = check(cls, Q)
        print(cls, 'chi', chi, 'coefficients compared', n, 'mismatches', bad)

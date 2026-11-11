from psi import *
import sympy as sp, pickle, itertools
Q = 22
A = phim21(Q+2); B = phi01(Q+2)
res = {}
for cls in ['1A','2A','3A','4B','5A','6A','7A','8A']:
    psi, W, d2, X = psiF1(cls, Q)
    N = ORDER[cls]
    divs = [d for d in range(1, N+1) if N % d == 0]
    d2s = [d for d in divs if d > 1]
    m4 = [('E4', d, E4(Q, d)) for d in divs]
    for a, b in itertools.combinations_with_replacement(d2s, 2):
        m4.append(('P', (a, b), smul(phiN(a, Q), phiN(b, Q), Q)))
    m2 = [('phi', d, phiN(d, Q)) for d in d2s]
    basis = [(name, mul(qser(s), A, Q)) for name in [(x[0], x[1]) for x in m4] for s in [dict(((y[0],y[1]),y[2]) for y in m4)[name]]]
    basis += [(('B', x[1]), mul(qser(x[2]), B, Q)) for x in m2]
    Wt = {k: v for k, v in W.items() if k[0] < Q//2}
    keys = sorted(set(Wt) | set(k for _, b in basis for k in b if k[0] < Q//2))
    M = sp.Matrix([[sp.Rational(b.get(k, 0)) for _, b in basis] for k in keys])
    v = sp.Matrix([sp.Rational(Wt.get(k, 0)) for k in keys])
    # least squares exact via pseudo: solve normal equations with rank handling
    sol, params = M.gauss_jordan_solve(v)
    sol = sol.subs({p: 0 for p in params})
    print(cls, 'd2', d2, 'free params', len(params))
    for (name, _), c in zip(basis, sol):
        if c != 0: print('   ', name, c)
    res[cls] = ([(name, c) for (name, _), c in zip(basis, sol) if c != 0], d2)
pickle.dump(res, open('Wfit.pkl', 'wb'))

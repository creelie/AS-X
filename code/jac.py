"""2D q,y series utilities with exact Fraction/int coefficients.
A series is a dict {(n, l): coeff}, truncated at q-order < Q."""
from fractions import Fraction
from collections import defaultdict
import math

def mul(a, b, Q):
    r = defaultdict(int)
    for (n1, l1), c1 in a.items():
        if c1 == 0: continue
        for (n2, l2), c2 in b.items():
            n = n1 + n2
            if n < Q:
                r[(n, l1 + l2)] += c1 * c2
    return {k: v for k, v in r.items() if v != 0}

def add(a, b, s=1):
    r = defaultdict(int, a)
    for k, v in b.items():
        r[k] += s * v
    return {k: v for k, v in r.items() if v != 0}

def scal(a, s):
    return {k: v * s for k, v in a.items() if v * s != 0}

def qser(coeffs, shift=0):
    """1-variable q-series list -> 2D series (l=0)."""
    return {(n + shift, 0): c for n, c in enumerate(coeffs) if c != 0}

# ---------- one-variable q series (lists) ----------
def smul(a, b, Q):
    r = [0] * Q
    for i, x in enumerate(a[:Q]):
        if x == 0: continue
        for j, y in enumerate(b[:Q - i]):
            r[i + j] += x * y
    return r

def sinv(a, Q):
    a0 = a[0]
    r = [0] * Q
    r[0] = Fraction(1) / a0 if not isinstance(a0, int) or a0 not in (1, -1) else a0
    for n in range(1, Q):
        s = 0
        for k in range(1, min(n, len(a) - 1) + 1):
            s += a[k] * r[n - k]
        r[n] = -s * r[0]
    return r

def eta_prod(shape, Q):
    """prod_a prod_n (1-q^{a n})^{m_a}  (no q^{1/24} prefactor). shape: dict a->m_a (m may be negative)."""
    r = [0] * Q
    r[0] = 1
    for a, m in shape.items():
        for n in range(1, Q):
            if a * n >= Q: break
            # multiply by (1-q^{an})^m
            if m > 0:
                for _ in range(m):
                    for i in range(Q - 1, a * n - 1, -1):
                        r[i] -= r[i - a * n]
            else:
                for _ in range(-m):
                    for i in range(a * n, Q):
                        r[i] += r[i - a * n]
    return r

def sigma(n, k):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)

def E2(Q, d=1):
    r = [0] * Q
    r[0] = 1
    for n in range(1, Q):
        if n % d == 0:
            r[n] = -24 * sigma(n // d, 1)
    return r

def E4(Q, d=1):
    r = [0] * Q
    r[0] = 1
    for n in range(1, Q):
        if n % d == 0:
            r[n] = 240 * sigma(n // d, 3)
    return r

def E6(Q, d=1):
    r = [0] * Q
    r[0] = 1
    for n in range(1, Q):
        if n % d == 0:
            r[n] = -504 * sigma(n // d, 5)
    return r

def phiN(N, Q):
    """(N E2(N tau) - E2(tau))/(N-1), constant term 1."""
    a = E2(Q, N); b = E2(Q, 1)
    return [Fraction(N * x - y, N - 1) for x, y in zip(a, b)]

# ---------- Jacobi forms ----------
def phim21(Q):
    """phi_{-2,1} = (y-2+1/y) prod (1-q^n y)^2 (1-q^n/y)^2 / (1-q^n)^4"""
    s = {(0, 1): 1, (0, 0): -2, (0, -1): 1}
    for n in range(1, Q):
        f = {(0, 0): 1, (n, 1): -1}
        g = {(0, 0): 1, (n, -1): -1}
        for _ in range(2):
            s = mul(s, f, Q); s = mul(s, g, Q)
    inv = eta_prod({1: -4}, Q)
    return mul(s, qser(inv), Q)

def phi01(Q):
    """phi_{0,1} = 4 sum_i theta_i(z)^2/theta_i(0)^2, i=2,3,4. Use phi01 = (3/pi^2) wp phi_{-2,1}; instead use
    phi01 = 12 phi_{-2,1} * wp/(2pi i)^2 ... simpler: phi01 = (phi_{-2,1} * ...). We compute via theta functions."""
    # theta_i(tau,z) with q^{1/8} units handled: use squares ratio, exponents in q^{1/2}
    # Work in units q^(1/8): build thetas as dicts with n in eighths.
    def th(kind, Qe):
        d = {}
        M = int(math.isqrt(8 * Qe)) + 3
        for k in range(-M, M + 1):
            if kind == 2:
                # sum q^{(k+1/2)^2/2} y^{k+1/2}: exponent (2k+1)^2/8 in q ; y^{(2k+1)/2}
                e = (2 * k + 1) ** 2  # in q^{1/8}
                d[(e, 2 * k + 1)] = d.get((e, 2 * k + 1), 0) + 1
            elif kind == 3:
                e = 4 * k * k
                d[(e, 2 * k)] = d.get((e, 2 * k), 0) + 1
            elif kind == 4:
                e = 4 * k * k
                d[(e, 2 * k)] = d.get((e, 2 * k), 0) + (-1) ** (k % 2)
        return {kk: v for kk, v in d.items() if kk[0] < 8 * Qe + 8}
    Qe = Q + 2
    tot = {}
    for kind in (2, 3, 4):
        t = th(kind, Qe)
        t2 = mul(t, t, 8 * Qe + 8)
        t0 = defaultdict(int)
        for (e, l), c in t.items():
            t0[e] += c
        # theta(0)^2 as series in q^{1/8}
        L = 8 * Qe + 8
        t0l = [0] * L
        for e, c in t0.items():
            if e < L: t0l[e] += c
        t02 = smul(t0l, t0l, L)
        # leading power
        lead = next(i for i, x in enumerate(t02) if x != 0)
        base = t02[lead:]
        inv = sinv(base, L)
        # t2 / t02: shift exponents by -lead
        r = defaultdict(int)
        for (e, l), c in t2.items():
            for j, iv in enumerate(inv):
                ee = e - lead + j
                if ee >= 8 * Q: break
                if iv != 0:
                    r[(ee, l)] += c * iv
        for (e, l), c in r.items():
            tot[(e, l)] = tot.get((e, l), 0) + 4 * c
    out = {}
    for (e, l), c in tot.items():
        if c != 0:
            assert e % 8 == 0 and l % 2 == 0, (e, l, c)
            out[(e // 8, l // 2)] = c
    return out

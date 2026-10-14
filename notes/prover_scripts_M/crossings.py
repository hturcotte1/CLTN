"""
HEURISTIC ONLY (floating point). For fixed beta, find all t>0 at which K_t = I + beta N + t J has a
nonzero purely imaginary eigenvalue, and the direction in which it crosses.
det(mu I + K_t) = phi(mu) + t omega(mu). At mu = i y: t = -phi(iy)/omega(iy) real.
Eigenvalue lambda = -mu of K. good crossing (K gains stability as t grows) <=> Re dlambda/dt > 0.
"""
import numpy as np
from numpy.polynomial import polynomial as P

def phi_omega(N, beta):
    m = len(N)
    C = np.eye(m) + beta * N
    phi = np.poly(-C)[::-1]                    # det(mu I + C), low->high
    # det(mu I + C + J) = phi + omega
    phiJ = np.poly(-(C + np.ones((m, m))))[::-1]
    omega = phiJ - phi
    return phi, omega[:m]

def crossings(N, beta):
    phi, om = phi_omega(N, beta)
    # F(y) = Im( phi(iy) * conj(omega(iy)) ) as a real polynomial in y
    m = len(N)
    def ipow_poly(c):
        # c low->high in mu, return complex coefficients in y for c(i y)
        return np.array([c[k] * (1j) ** k for k in range(len(c))])
    pc = ipow_poly(phi); oc = ipow_poly(om)
    prod = P.polymul(pc, np.conj(oc))
    imag = prod.imag
    roots = P.polyroots(imag) if np.any(np.abs(imag) > 1e-12) else np.array([])
    out = []
    for y in roots:
        if abs(y.imag) > 1e-7 or y.real <= 1e-9:
            continue
        y = y.real
        mu = 1j * y
        o = P.polyval(mu, om)
        if abs(o) < 1e-12:
            continue
        t = -P.polyval(mu, phi) / o
        if abs(t.imag) > 1e-6 * (1 + abs(t)) or t.real <= 0:
            continue
        t = t.real
        dphi = P.polyder(phi); dom = P.polyder(om)
        dmu = -o / (P.polyval(mu, dphi) + t * P.polyval(mu, dom))
        dlam = -dmu
        out.append((t, y, dlam.real))
    return out

def a_vec(N, beta):
    m = len(N)
    C = np.eye(m) + beta * N
    d = np.linalg.det(C)
    return d * np.linalg.solve(C, np.ones(m))

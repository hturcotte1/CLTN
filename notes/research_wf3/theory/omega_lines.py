import json, itertools, sys
from fractions import Fraction as F
import tlib, weighted, boundary_tools as bt, rstar, omega_weights as ow
d6s = ['&FNFL~[udo?'] + [r['d6'] for r in json.load(open('cube_n8_witnesses.json'))]
amax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
res = {}
for s in d6s:
    A = tlib.from_d6(s); H, mu0 = ow.closed_twin_quotient(A)
    m = len(H)
    for x in range(m):
        line = []
        okline = True
        for a in range(1, amax + 1):
            mu = list(mu0); mu[x] = mu0[x] + a - 1; mu = tuple(mu)
            if ow.mod_phi3(ow.phi_mu_coeffs(H, mu)) != (0, 0): okline = False; break
            B = weighted.blowup(H, mu)
            mult, others = bt.boundary_split(B)
            pos = tlib.positive_at(B, F(2) + F(1, 10**4))
            R = rstar.R_at(B)[0] if (mult == 1 and others) else None
            good = (mult == 1 and others and pos and R is not None and R > F(1, 8))
            line.append((a, mult, others, pos, None if R is None else str(R), good))
        if okline and len(line) == amax:
            ng = sum(1 for t in line if t[-1])
            print(s, 'node', x, 'goodcount', ng, line if ng >= 2 else '', flush=True)

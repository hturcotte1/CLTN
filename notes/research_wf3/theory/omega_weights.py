"""For a base graph H and integer node weights mu (closed-twin blow-up sizes), decide exactly whether omega is a zero of
Phi_mu(z) = mu^T adj(zI - D M) 1 = chi_{DM - 1 mu^T}(z) - chi_{DM}(z)   (D-picture; weighted.py, Lemma W).
Reduction mod z^2+z+1 is exact integer arithmetic."""
import sys, json, itertools
import tlib, weighted, boundary_tools as bt
from fastregion import icharpoly
def phi_mu_coeffs(H, mu):
    m = len(H)
    X = [[(0 if i == j else 1 - H[i][j]) * mu[j] for j in range(m)] for i in range(m)]
    Y = [[X[i][j] - mu[j] for j in range(m)] for i in range(m)]
    c1 = icharpoly([r[:] for r in Y]); c0 = icharpoly([r[:] for r in X])
    return [a - b for a, b in zip(c1, c0)]          # hi->lo, leading 0
def mod_phi3(p):
    """p hi->lo integer; return (u, v) with p(omega) = u + v*omega."""
    # powers of omega cycle with period 3: omega^0=1, omega^1=omega, omega^2=-1-omega
    u = v = 0; d = len(p) - 1
    for i, c in enumerate(p):
        e = (d - i) % 3
        if e == 0: u += c
        elif e == 1: v += c
        else: u -= c; v -= c
    return u, v
def closed_twin_quotient(A):
    n = len(A); cls = []
    assigned = [-1] * n
    for u in range(n):
        if assigned[u] >= 0: continue
        assigned[u] = len(cls); cl = [u]
        for v in range(u + 1, n):
            if assigned[v] < 0 and A[u][v] and A[v][u] and all(A[i][u] == A[i][v] and A[u][i] == A[v][i] for i in range(n) if i not in (u, v)):
                assigned[v] = assigned[u]; cl.append(v)
        cls.append(cl)
    reps = [c[0] for c in cls]
    H = [[A[a][b] for b in reps] for a in reps]
    return H, tuple(len(c) for c in cls)
if __name__ == '__main__':
    d6s = ['&FNFL~[udo?'] + [r['d6'] for r in json.load(open('cube_n8_witnesses.json'))]
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for s in d6s:
        A = tlib.from_d6(s); H, mu0 = closed_twin_quotient(A)
        assert mod_phi3(phi_mu_coeffs(H, mu0)) == (0, 0)
        sols = [mu for mu in itertools.product(range(1, B + 1), repeat=len(H)) if mod_phi3(phi_mu_coeffs(H, mu)) == (0, 0)]
        print(s, 'base m=%d' % len(H), 'mu0', mu0, 'omega-locus points in box:', len(sols), sols[:12], flush=True)

"""Sanity (not proof): Lemma I.2 -- bounded eigenvalues of K=C+tJ approach spec(C|perp); big one = m t + kappa + o(1).
mpmath at 60 digits; random digraphs, C = I + beta N."""
import random, mpmath as mp
mp.mp.dps = 60
random.seed(7)
def helmert(m):
    Q = mp.zeros(m, m-1)
    for k in range(1, m):
        for i in range(k): Q[i, k-1] = 1
        Q[k, k-1] = -k
        nrm = mp.sqrt(k*(k+1))
        for i in range(m): Q[i, k-1] /= nrm
    return Q
worst = []
for trial in range(60):
    m = random.randint(2, 7)
    A = [[0 if i == j else int(random.random() < 0.6) for j in range(m)] for i in range(m)]
    beta = mp.mpf(random.randint(110, 600))/100
    N = mp.matrix([[0 if a == b else 1 - A[b][a] for b in range(m)] for a in range(m)])
    C = mp.eye(m) + beta*N
    U = helmert(m); Cp = U.T*C*U
    zs = [Cp[0,0]] if m == 2 else list(mp.eig(Cp, left=False, right=False))
    e = sum(N[a, b] for a in range(m) for b in range(m)); kappa = 1 + beta*e/m
    for t in [mp.mpf(10)**6, mp.mpf(10)**10]:
        K = C + t*mp.ones(m, m)
        ev = list(mp.eig(K, left=False, right=False))
        big = max(ev, key=lambda z: abs(z)); rest = [z for z in ev if z is not big]
        assert abs(mp.im(big)) < mp.mpf(10)**-30
        # match rest to zs greedily (sanity only)
        err = 0; pool = list(zs)
        for z in rest:
            j = min(range(len(pool)), key=lambda i: abs(pool[i]-z)); err = max(err, abs(pool[j]-z)); pool.pop(j)
        worst.append((float(t), m, float(err), float(abs(big - m*t - kappa))))
for t in [1e6, 1e10]:
    ws = [w for w in worst if w[0] == t]
    print(f"t={t:g}: max |bounded eig - spec(C|perp)| = {max(w[2] for w in ws):.3e}; max |big - m t - kappa| = {max(w[3] for w in ws):.3e}")

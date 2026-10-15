"""Arms-forest H/Z design ('Theorem candidate'): H = 3-cycle c0->c1->c2->c0 (D-arcs = non-edges) + arms a (D-arc c2 -> a),
each arm a has exactly one Z child z_a (D-arc a -> z_a) and z_a has exactly one D-arc z_a -> target(a), target in
{c2, c1, c0} or a parent arm.  Claimed: with theta = e_{c2}, chi_{D - theta 1^T}(z) = z^(m-d) p(z),
p(z) = z^d + z^(d-1) + z^(d-2) + sum_a z^(d - len(a)),  len(a) = 2 depth(a) + delta(exit), delta(c2,c1,c0) = 0,1,2."""
import numpy as np
def build_forest(tau):
    """tau[j] (j >= 2) = number of arms with len j.  j = 2 -> exit c2 (depth 1); even j >= 4 -> exit c0 at depth j/2 - 1; odd j -> exit c1 at depth (j-1)/2.
    Arms at depth k >= 2 all hang under the first arm of depth k-1 of the same exit (a chain spine + bushes)."""
    nodes = ['c0', 'c1', 'c2']; arcs = [('c0', 'c1'), ('c1', 'c2'), ('c2', 'c0')]
    spine = {}   # (exit, depth) -> first arm name
    arms = []
    for j in sorted(tau):
        for i in range(tau[j]):
            if j == 2: ex, dep = 'c2', 1
            elif j % 2 == 0: ex, dep = 'c0', j // 2 - 1
            else: ex, dep = 'c1', (j - 1) // 2
            a = 'a%d_%d' % (j, i); z = 'z%d_%d' % (j, i)
            tgt = ex if dep == 1 else spine[(ex, dep - 1)]
            nodes += [a, z]; arcs += [('c2', a), (a, z), (z, tgt)]
            spine.setdefault((ex, dep), a)
            arms.append((a, z, tgt, j))
    idx = {v: i for i, v in enumerate(nodes)}; m = len(nodes)
    D = np.zeros((m, m), int)
    for u, v in arcs: D[idx[u], idx[v]] = 1
    return D, idx, arms
def p_coeffs(tau):
    d = max(tau) if tau else 2
    d = max(d, 2)
    c = [0] * (d + 1); c[0] = 1; c[1] = 1; c[2] += 1
    for j, t in tau.items(): c[j] += t
    return c
if __name__ == '__main__':
    tau = {2: 1, 3: 2, 4: 2, 5: 1, 6: 1}
    D, idx, arms = build_forest(tau); m = len(D)
    th = np.zeros(m); th[idx['c2']] = 1
    L = D - np.outer(th, np.ones(m))
    print('m', m, 'charpoly(L) =', np.round(np.poly(L), 6))
    c = p_coeffs(tau); print('p =', c, ' padded:', c + [0] * (m - len(c) + 1))
    print('roots p', np.round(np.sort_complex(np.roots(c)), 4))

"""
theorem_U.py -- exact instances of the construction in the proof of Theorem U (sup_n q_n* = infinity).

For N >= 3 and an integer W >= 1:  s = 1/(2N),  target  p_W(z) = prod_{k=1..N} ((z + s)^2 + W^2 k^2)  (exact rationals),
c_1(W) = 2Ns = 1 exactly; c_j := nearest integer to c_j(W) for j >= 2 (exact Fraction rounding).  The arms-forest core
D_{N,W} realises p = sum c_j z^(2N-j) (Lemma T).  For each N we find the least W in a range for which the exact
certificate of certify_forest.certify passes at beta* = N (so q -> N - 1), and report the core size and k0.
Floating point: none.  Run: python3 theorem_U.py [Nmax]  (writes theorem_U_instances.json)
"""
import sys, json, time
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN'); sys.path.insert(0, '/home/user/CLTN/notes/research_wf3/theory')
from certify_forest import certify, closed_form_positivity
from forest_exact import roots_right_of

def polymul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out

def target(N, W):
    s = F(1, 2 * N); p = [F(1)]
    for k in range(1, N + 1): p = polymul(p, [F(1), 2 * s, s * s + W * W * k * k])
    return p                                        # hi->lo, p[1] == 1

def nearest(x):
    fl = x.numerator // x.denominator
    return fl + 1 if (x - fl) * 2 >= 1 else fl

def rounded(N, W):
    p = target(N, W); assert p[0] == 1 and p[1] == 1
    c = [1, 1] + [nearest(x) for x in p[2:]]
    assert c[2] >= 1 and all(v >= 1 for v in c[3:])
    return c

if __name__ == '__main__':
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out = []
    for N in range(3, Nmax + 1):
        t0 = time.time(); found = None
        for W in range(1, 400):
            c = rounded(N, W)
            if roots_right_of(c, F(1, N))[0] is not True: continue
            if not closed_form_positivity(c, F(N))['ok']: continue
            rec = certify(c, F(N))
            if rec['certified']:
                found = (W, rec); break
        if found is None:
            print('N=%d: no W < 400 found' % N, flush=True); continue
        W, rec = found
        rec['N'] = N; rec['W'] = W
        out.append(rec)
        print('N=%d (beta*=%d, q->%d): W=%d  d=%d  core m=%s (%d digits)  k0=%s  %.1fs' % (N, N, N - 1, W, rec['d'], str(rec['m_core'])[:12], len(str(rec['m_core'])), rec['C3_k0'], time.time() - t0), flush=True)
        json.dump(out, open('theorem_U_instances.json', 'w'), indent=1)

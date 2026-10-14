"""Exact cross-check of filter7p's a_eval and Descartes coefficients against sympy (DomainMatrix over ZZ[beta])."""
import sys, random, subprocess, math
from fractions import Fraction as F
from sympy import symbols, Poly, ZZ
from sympy.polys.matrices import DomainMatrix
b, u = symbols('b u')
R = ZZ[b]
HERE = sys.argv[1]; NG = int(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed)

def to_d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for x in bits[k:k+6]: v = 2*v + x
        s += chr(v + 63)
    return s

def a_polys(A):
    n = len(A)
    Nm = [[1 if (i != j and A[j][i] == 0) else 0 for j in range(n)] for i in range(n)]   # N[i][j]=1 iff j-/->i (CLAUDE.md W convention)
    out = []
    for i in range(n):
        rows = []
        for r in range(n):
            row = []
            for c in range(n):
                if c == i: row.append(R(1))
                else: row.append(R.convert((1 if r == c else 0)) + R.convert(Nm[r][c]) * R.gens[0])
            rows.append(row)
        dm = DomainMatrix(rows, (n, n), R)
        out.append(Poly(R.to_sympy(dm.det()), b))
    return out

def dominated(A):
    n = len(A)
    for j in range(n):
        for k in range(n):
            if j == k: continue
            if A[j][k] == 1 and A[k][j] == 0 and all(not (A[i][j] == 1 and A[i][k] == 0) for i in range(n) if i not in (j, k)):
                return True
    return False

SAMPLES = [(1,1),(1,2),(1,4),(1,8),(2,1),(1,16),(4,1),(1,32),(8,1),(16,1),(48,1),(1,64),(3,8),(3,4)]
graphs = []
for g in range(NG):
    n = 8 if g % 10 == 9 else 7
    dens = rng.choice([0.05, 0.2, 0.4, 0.6, 0.75, 0.85, 0.92, 0.97, 1.0])
    A = [[0 if i == j else (1 if rng.random() < dens else 0) for j in range(n)] for i in range(n)]
    samp = [(2*q + p, q) for (p, q) in SAMPLES]
    samp += [(rng.randint(129, 590 if n == 7 else 191), 64) for _ in range(6)]   # q = 64, inside the proven int128-safe range (bounds.py)
    samp += [(rng.randint(1, 50), 1) for _ in range(2)] + [(rng.randint(1, 127), 64)]   # beta up to 50 at q=1; beta < 2 at q=64
    graphs.append((A, samp))
inp = ''.join(to_d6(A) + ' ' + str(len(s)) + ' ' + ' '.join(f'{p} {q}' for p, q in s) + '\n' for A, s in graphs)
res = subprocess.run([HERE + '/harness'], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
assert len(res) == NG
bad = 0; stats = {'drop': 0, 'keep': 0, 'dom': 0, 'n8': 0, 'evals': 0, 'P': 0}
for (A, samp), line in zip(graphs, res):
    n = len(A); stats['n8'] += n == 8
    parts = line.split('#')
    head = parts[0].split()
    if len(head) < 2: print("PARSE", repr(line[:200])); bad += 1; continue
    df = int(head[1].split("=")[1])
    if df != (0 if dominated(A) else 1): print('DOMFREE MISMATCH', head); bad += 1
    stats['dom'] += 1 - df
    ap = a_polys(A)
    mf = math.factorial(n)
    drop_true = False
    for i in range(n):
        seg = parts[1 + i].split(':')[1].split()
        ccoef = [int(x) for x in seg[:-1]]; np_flag = int(seg[-1].split('=')[1])
        sh = Poly(ap[i].as_expr().subs(b, 2 + u), u)
        true = [0] * (n + 1)
        for (k,), c in sh.terms(): true[k] = int(c) * mf
        if ap[i].degree() > n - 1: print('DEGREE > m-1', line[:20]); bad += 1
        if ccoef != true: print('COEF MISMATCH', head[0], i, ccoef, true); bad += 1
        nonpos = all(c <= 0 for c in true)
        if nonpos != bool(np_flag): print('NP MISMATCH'); bad += 1
        drop_true |= nonpos
    cdrop = int(parts[1 + n].split('=')[1])
    if cdrop != drop_true: print('DROP MISMATCH'); bad += 1
    stats['drop' if drop_true else 'keep'] += 1
    evs = [s.split() for s in parts[2 + n][2:].split(';') if s.strip()]
    assert len(evs) == len(samp)
    for (p, q), ev in zip(samp, evs):
        for i in range(n):
            tv = F(ap[i].eval(F(p, q))) * q**n
            assert tv.denominator == 1
            stats['evals'] += 1
            if int(ev[i]) != tv.numerator:
                # flag overflow vs genuine
                print('EVAL MISMATCH', head[0], 'p/q', p, q, 'i', i, ev[i], tv.numerator, 'in filter samples' if (p, q) in [(2*b_ + a_, b_) for a_, b_ in SAMPLES] else ''); bad += 1
            elif (int(ev[i]) > 0) != (tv > 0): bad += 1
print(stats, 'BAD', bad)

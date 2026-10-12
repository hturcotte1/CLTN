"""
test_ctln.py -- validation suite for the exact CTLN verifier.

These tests do NOT check the code against itself. Each one checks it against a
fact that is known independently from the theory. If any of these fail, nothing
the pipeline produces downstream is trustworthy.

Run:  python3 test_ctln.py
"""

import random
from fractions import Fraction as F
from itertools import combinations

from ctln import (
    analyze_support, charpoly, _charpoly_direct, hurwitz_stable,
    identity_check_IminusW, is_clique, is_target_free_clique, legal,
    q_of, scan_supports, targets,
)

PASS, FAIL = [], []


def check(name, cond, note=""):
    (PASS if cond else FAIL).append(name)
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + (f"   [{note}]" if note else ""))


def complete_bidirected(m, n=None):
    """K_m on nodes 0..m-1, padded out to n nodes with no other edges."""
    n = n or m
    A = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            if i != j:
                A[i][j] = 1
    return A


def random_digraph(n, p=0.5, rng=None):
    rng = rng or random
    return [[0 if i == j else (1 if rng.random() < p else 0) for j in range(n)]
            for i in range(n)]


# ---------------------------------------------------------------------------
print("\n=== 1. linear algebra primitives ===")

random.seed(11)
ok = True
for _ in range(40):
    m = random.randint(1, 5)
    M = [[F(random.randint(-6, 6), random.randint(1, 4)) for _ in range(m)]
         for _ in range(m)]
    if charpoly(M) != _charpoly_direct(M):
        ok = False
        break
check("charpoly matches independent interpolation on 40 random matrices", ok)

# Hurwitz against a polynomial with known roots: (lam+1)(lam+2)(lam+3)
check("Hurwitz accepts a known-stable cubic",
      hurwitz_stable([F(1), F(6), F(11), F(6)])[0] is True)
# (lam-1)(lam+2)(lam+3) = lam^3 +4lam^2 +lam -6  -> one positive root
check("Hurwitz rejects a cubic with a positive root",
      hurwitz_stable([F(1), F(4), F(1), F(-6)])[0] is False)
# (lam^2+1)(lam+1) -> purely imaginary pair, exact boundary
check("Hurwitz flags a boundary case rather than guessing",
      hurwitz_stable([F(1), F(1), F(1), F(1)])[0] is None)

# ---------------------------------------------------------------------------
print("\n=== 2. closed-form identity  I - W_s = -dI + (1+d)J - (e+d)A_s ===")

random.seed(7)
eps, delta = F(1, 10), F(9, 10)
ok = True
for _ in range(30):
    n = random.randint(2, 6)
    A = random_digraph(n, 0.5)
    size = random.randint(1, n)
    sigma = sorted(random.sample(range(n), size))
    if not identity_check_IminusW(A, sigma, eps, delta):
        ok = False
        break
check("closed form agrees with directly built matrix on 30 random cases", ok)

# ---------------------------------------------------------------------------
print("\n=== 3. cliques: the closed-form solution ===")

# For a clique of size m the fixed point is x_i = theta / (eps + m(1-eps)).
ok = True
for m in range(1, 6):
    A = complete_bidirected(m)
    rec = analyze_support(A, list(range(m)), eps, delta, theta=F(1))
    expected = F(1) / (eps + m * (1 - eps))
    got = F(rec["x"][0])
    if not (rec["stable_fixed_point"] and all(F(v) == expected for v in rec["x"])):
        ok = False
        print(f"      m={m}: expected {expected}, got {rec['x']}, "
              f"stable={rec['stable']}")
check("K_m whole-graph clique is a stable fixed point with x = theta/(eps+m(1-eps))",
      ok, "m = 1..5")

# ---------------------------------------------------------------------------
print("\n=== 4. the target rule ===")

# Clique {0,1,2} plus node 3 receiving from ALL of them -> sigma has a target,
# so the off-support condition must FAIL and it must not be a fixed point.
A = complete_bidirected(3, n=4)
for j in range(3):
    A[j][3] = 1
rec = analyze_support(A, [0, 1, 2], eps, delta)
check("clique WITH a target is not a fixed point",
      (not rec["is_fixed_point"]) and (not rec["off_support_ok"]),
      f"targets={rec['targets']}")

# Same clique, but node 3 only receives from 0 and 1 -> no target -> survives.
A2 = complete_bidirected(3, n=4)
A2[0][3] = 1
A2[1][3] = 1
rec2 = analyze_support(A2, [0, 1, 2], eps, delta)
check("clique with NO target is a stable fixed point",
      rec2["stable_fixed_point"] and rec2["target_free"],
      f"targets={rec2['targets']}")

# ---------------------------------------------------------------------------
print("\n=== 5. singletons: {i} is stable iff i is a sink ===")

random.seed(3)
ok = True
detail = ""
for _ in range(25):
    n = random.randint(2, 5)
    A = random_digraph(n, 0.45)
    for i in range(n):
        rec = analyze_support(A, [i], eps, delta)
        is_sink = all(A[i][k] == 0 for k in range(n) if k != i)
        if rec["stable_fixed_point"] != is_sink:
            ok = False
            detail = f"node {i}, sink={is_sink}, got={rec['stable_fixed_point']}"
            break
    if not ok:
        break
check("singleton support is a stable fixed point exactly when the node is a sink",
      ok, detail)

# ---------------------------------------------------------------------------
print("\n=== 6. THE REAL TEST: conjecture must hold in Geneson's proven region ===")
# He proved the target-free clique conjecture holds for q >= n - 2 - (n-3)eps/2.
# So at large q there must be ZERO violations. Any violation here means the
# verifier is broken, NOT that we found something.

random.seed(99)
violations_found = []
boundary_found = 0
graphs_tested = 0

for n in (3, 4, 5):
    eps_s = F(1, 100)
    # pick delta so that q comfortably exceeds n-2
    delta_s = F(n)  # q = n(1 - 1/100)/(1/100) = 99n, far above n-2
    assert legal(eps_s, delta_s), "parameters must be legal"
    for _ in range(40):
        A = random_digraph(n, 0.5)
        graphs_tested += 1
        out = scan_supports(A, eps_s, delta_s)
        boundary_found += len(out["boundary_cases"])
        if out["conjecture_violations"]:
            violations_found.append((n, A, out["conjecture_violations"][0]))

check("zero conjecture violations at large q, as the theorem requires",
      not violations_found,
      f"{graphs_tested} graphs, q={q_of(F(1,100), F(3))} for n=3, "
      f"{boundary_found} boundary cases")
if violations_found:
    n, A, v = violations_found[0]
    print(f"      first violation: n={n} A={A} sigma={v['sigma']} "
          f"stable={v['stable_fixed_point']} target_free={v['target_free']}")

# ---------------------------------------------------------------------------
print("\n=== 7. known reductions that prune the n=6 search ===")
# The conjecture is already proven for oriented graphs (no 2-cycles) and for
# symmetric graphs. Confirm no violations there at SMALL q, where counterexamples
# are otherwise possible -- these classes should be clean at every legal q.

random.seed(21)
eps_s, delta_s = F(1, 2), F(21, 10)   # q = 2.1*0.5/0.5 = 2.1, small
assert legal(eps_s, delta_s)
oriented_bad, symmetric_bad = [], []
for _ in range(60):
    n = 4
    A = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        r = random.random()
        if r < 0.33:
            A[i][j] = 1          # one direction only -> oriented
        elif r < 0.66:
            A[j][i] = 1
    if scan_supports(A, eps_s, delta_s)["conjecture_violations"]:
        oriented_bad.append(A)

    S = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if random.random() < 0.5:
            S[i][j] = S[j][i] = 1   # symmetric
    if scan_supports(S, eps_s, delta_s)["conjecture_violations"]:
        symmetric_bad.append(S)

check("no violations on oriented graphs (proven class)", not oriented_bad,
      f"q={q_of(eps_s, delta_s)}")
check("no violations on symmetric graphs (proven class)", not symmetric_bad,
      f"q={q_of(eps_s, delta_s)}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print(f"{len(PASS)} passed, {len(FAIL)} failed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
    raise SystemExit(1)
print("Verifier is consistent with every independently known fact tested.")
print("NOTE: the sharpest remaining test is reproducing Geneson's 6-node")
print("counterexample from arXiv:2607.21396. That paper's explicit graph is NOT")
print("encoded here -- pulling it and confirming the verifier flags it is the")
print("first task for the run, before any new search begins.")

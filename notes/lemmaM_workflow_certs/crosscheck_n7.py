"""Independent exact cross-checks for the n=7 Lemma-M counterexample (no floats in any verdict):
 (a) sympy: W built from CLAUDE.md conventions directly (not via ctln), exact charpoly of -I+W,
     Routh array -> number of sign changes = number of RHP roots; and x = (I-W)^{-1} 1 exactly.
 (b) graphs.is_domination_free
 (c) factorisation of chi_perp over Q (exact value of alpha(A|perp) where possible)
 (d) canonical digraph6 (nauty labelg) and membership in n7 survivor files.
Re-run: python3 crosscheck_n7.py"""
import sys, subprocess, gzip
sys.path.insert(0, "/home/user/CLTN")
import sympy as sp
from graphs import is_domination_free

ROWS = ["0110000", "0000111", "0001111", "1100100", "1001011", "1111101", "1111110"]
A = [[int(c) for c in r] for r in ROWS]
n = len(A)
lam = sp.Symbol("lam")


def routh_rhp(coeffs):
    """number of sign changes in the first column of the Routh array (exact). None if a zero pivot occurs."""
    c = [sp.nsimplify(x) for x in coeffs]
    r0 = c[0::2]; r1 = c[1::2]
    rows = [r0, r1 + [0] * (len(r0) - len(r1))]
    while len(rows) < len(c):
        a, b = rows[-2], rows[-1]
        if b[0] == 0:
            return None
        new = [(b[0] * a[k + 1] - a[0] * b[k + 1]) / b[0] for k in range(len(a) - 1)] + [0]
        rows.append(new)
    col = [r[0] for r in rows]
    if any(x == 0 for x in col):
        return None
    return sum(1 for u, v in zip(col, col[1:]) if (u > 0) != (v > 0))


def check(beta, t):
    eps = sp.Rational(1) / (1 + t); delta = (beta - 1) / (1 + t)
    W = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                W[i, j] = (-1 + eps) if A[j][i] == 1 else (-1 - delta)
    x = (sp.eye(n) - W).LUsolve(sp.ones(n, 1))
    cp = sp.Poly((lam * sp.eye(n) - (-sp.eye(n) + W)).det(method="berkowitz"), lam).all_coeffs()
    q = delta * (1 - eps) / eps
    return q, min(x), routh_rhp(cp)


beta = sp.Rational(100, 49)
for t in [30, 35, 45, 60, 1000]:
    q, xmin, rhp = check(beta, sp.Integer(t))
    print(f"(a) beta=100/49 t={t}: q={q} (legal {q > 1}), min x = {xmin} (>0: {xmin > 0}), RHP roots of -I+W: {rhp}")
print("(b) domination-free:", is_domination_free(A))
J = sp.ones(n, n); Am = sp.Matrix(A)
chiperp = sp.factor(sp.expand(((lam * sp.eye(n) - Am).det(method="berkowitz") -
                               (lam * sp.eye(n) - Am - J).det(method="berkowitz")) / n))
print("(c) chi_perp factored over Q:", chiperp)
d6 = "&" + chr(n + 63)
bits = [A[i][j] for i in range(n) for j in range(n)]
bits += [0] * (-len(bits) % 6)
d6 += "".join(chr(63 + int("".join(map(str, bits[k:k + 6])), 2)) for k in range(0, len(bits), 6))
can = subprocess.run(["/home/user/CLTN/n7/labelg", "-q"], input=d6 + "\n", capture_output=True, text=True).stdout.strip()
print("(d) digraph6:", d6, " canonical:", can)
hit = [l.split()[1] for l in gzip.open("/home/user/CLTN/n7/n7_survivors_all.txt.gz", "rt") if l.split()[0] == can]
print("    in n7_survivors_all (limit S/B list):", hit or "NOT PRESENT")

# Recipe: check one candidate with sympy (exact rationals).  Usage: python3 recipe_check.py '<digraph6>' eps delta
import sys; from sympy import Matrix, Rational as R, symbols, Poly, eye, ones
d6, eps, delta = sys.argv[1], R(sys.argv[2]), R(sys.argv[3]); n = ord(d6[1]) - 63
bits = ''.join(format(ord(c) - 63, '06b') for c in d6[2:]); A = Matrix(n, n, [int(bits[i*n+j]) for i in range(n) for j in range(n)])
IW = -delta*eye(n) + (1+delta)*ones(n, n) - (eps+delta)*A.T          # I - W
x = IW.LUsolve(ones(n, 1)); s = symbols('s'); c = Poly((s*eye(n) + IW).det(), s).all_coeffs()  # charpoly of -I+W
H = Matrix(n, n, lambda i, j: c[2*j-i+1] if 0 <= 2*j-i+1 <= n else 0); D = [H[:k, :k].det() for k in range(1, n+1)]
print('q =', delta*(1-eps)/eps, '| x > 0:', all(v > 0 for v in x), '| clique:', all(A[i, j] == 1 for i in range(n) for j in range(n) if i != j))
print('charpoly coeffs:', c); print('Hurwitz minors:', D); print('all minors > 0 (stable):', all(d > 0 for d in D), '| any zero (UNDECIDED):', any(d == 0 for d in D))

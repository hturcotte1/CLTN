s = open('check_lemmas.py').read()
reps = [('''def psi(Dm):
    n = Dm.shape[0]
    return sp.expand((sp.ones(1, n) * (z * sp.eye(n) - Dm).adjugate() * sp.ones(n, 1))[0])

def chi(Mm):
    n = Mm.shape[0]; return sp.expand((z * sp.eye(n) - Mm).det(method='berkowitz'))''',
'''def chi(Mm):
    return Mm.charpoly(z).as_expr()

def uadjv(X, u, v):
    """v^T adj(zI - X) u  via the matrix determinant lemma: det(zI - X + u v^T) - det(zI - X)."""
    return sp.expand(chi(X - u * v.T) - chi(X))

def psi(Dm):
    n = Dm.shape[0]
    return uadjv(Dm, sp.ones(n, 1), sp.ones(n, 1))'''),
("rhs = sp.expand(z ** (n - m) * (sp.Matrix([mu]) * (z * sp.eye(m) - D * M).adjugate() * sp.ones(m, 1))[0])",
 "rhs = sp.expand(z ** (n - m) * uadjv(D * M, sp.ones(m, 1), sp.Matrix(mu)))"),
("l2 = sp.expand((z * sp.eye(n) + Ap).det(method='berkowitz'))", "l2 = chi(-Ap)"),
("r2 = sp.expand((z - 1) ** (n - m) * (z * sp.eye(m) + AQ).det(method='berkowitz'))", "r2 = sp.expand((z - 1) ** (n - m) * chi(-AQ))"),
("l3 = sp.expand((sp.ones(1, n) * (z * sp.eye(n) + Ap).adjugate() * sp.ones(n, 1))[0])", "l3 = uadjv(-Ap, sp.ones(n, 1), sp.ones(n, 1))"),
("r3 = sp.expand((z - 1) ** (n - m) * (sp.Matrix([mu]) * (z * sp.eye(m) + AQ).adjugate() * sp.ones(m, 1))[0])", "r3 = sp.expand((z - 1) ** (n - m) * uadjv(-AQ, sp.ones(m, 1), sp.Matrix(mu)))"),
("(sp.ones(1, m) * (z * sp.eye(m) - D).adjugate() * Bv)[0]", "uadjv(D, Bv, sp.ones(m, 1))")]
for a, b in reps:
    assert a in s, a[:40]
    s = s.replace(a, b)
open('check_lemmas.py', 'w').write(s)

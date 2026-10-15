s = open('check_F_algebra.py').read()
old_start = s.index("g = 1 - 1/b\n")
old_end = s.index("print('    w_G1(b) =', w)")
new = """Mb = sp.Matrix(6, 6, lambda i, j: b*A0[j][i] + ((b - 1) if i == j else 0))   # = b (A^T + g I)
detMb = sp.expand(Mb.det(method='berkowitz'))
adj1 = [sp.expand(e) for e in (Mb.adjugate(method='berkowitz') * sp.ones(6, 1))]
assert all(sp.expand(sum(Mb[i, j]*adj1[j] for j in range(6)) - detMb) == 0 for i in range(6))
w = [sp.factor(b*e/detMb) for e in adj1]       # w = b adj(Mb) 1 / det Mb
print('    det Mb =', sp.factor(detMb)); print('    adj(Mb)1 =', [sp.factor(e) for e in adj1])
"""
s = s[:old_start] + new + s[old_end:]
s = s.replace("for kk in range(0, 9):", "for kk in range(0, 7):").replace("k = 0..8: OK", "k = 0..6: OK")
s = s.replace("    return list(Mx.LUsolve(sp.ones(n, 1)))", "    import sys; sys.path.insert(0, '/home/user/CLTN'); import ctln\n    M2 = [[Fr(A[j][i]) + (gg if i == j else 0) for j in range(n)] for i in range(n)]\n    return ctln.solve_exact(M2, [Fr(1)]*n)")
open('check_F_algebra.py', 'w').write(s)

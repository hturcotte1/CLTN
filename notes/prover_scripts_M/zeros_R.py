"""HEURISTIC. For a digraph A: zeros z0 of s(z) = 1^T (zI + A)^{-1} 1 (= eigenvalues of -A|perp that move with c)
and the first-order coefficient R(z0) = Re(1/s'(z0)): eigenvalue of cJ - A near z0 is z0 + 1/(c s'(z0)) + O(1/c^2)."""
import numpy as np
def compress(m):
    V = np.zeros((m, m - 1))
    for j in range(m - 1): V[j, j] = 1; V[m - 1, j] = -1
    U, _ = np.linalg.qr(V); return U
def zeros_R(A):
    m = len(A); U = compress(m)
    ev = np.linalg.eigvals(-(U.T @ A @ U))
    one = np.ones(m); out = []
    for z0 in ev:
        M = z0 * np.eye(m) + A
        try:
            R1 = np.linalg.solve(M, one)
        except np.linalg.LinAlgError:
            out.append((z0, None, 'pole')); continue
        s = one @ R1
        if abs(s) > 1e-7 * (1 + np.abs(R1).max()):
            out.append((z0, None, 'fixed')); continue
        sp = -(R1 @ R1) if False else -(one @ np.linalg.solve(M, R1))
        out.append((z0, (1 / sp).real, 'moving'))
    return out
if __name__ == '__main__':
    A = np.array([[0, 0, 1, 0, 0], [1, 0, 1, 0, 0], [1, 1, 0, 0, 0], [0, 1, 1, 0, 1], [0, 1, 1, 1, 0]], float)
    print(zeros_R(A))
    m = 5
    for c in [100, 1000, 10000]:
        print(c, sorted(np.linalg.eigvals(c * np.ones((m, m)) - A), key=lambda z: z.real)[:2])

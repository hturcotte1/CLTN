"""
EXACT certification that Lemma M (both M_s and M_w), WITH the positivity hypothesis a(beta) > 0, is FALSE
at m = 7.  No floating point in any decision below; floats were used only to locate the point.

For a digraph G (digraph6 given on the command line) and rational beta, t_lo < t_hi:
  (1) analyze_support at (eps, delta) for t = t_lo: legal, positive, stable  (an SFP at a legal finite point)
  (2) analyze_support at t = t_hi (same beta): legal, positive, NOT stable  -> T(beta) not upward closed (M_s false)
  (3) omega(mu; beta) = det(mu I + C + J) - det(mu I + C) is NOT Hurwitz (exact Hurwitz minors)
      <=> alpha(A|perp) >= -(1 - 1/beta)  (Lemma C'/I)  -> M_w false if the inequality is strict;
      strictness: some Hurwitz minor is strictly negative (not merely zero), and omega(0)=s != 0.
"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs

def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A

def charpoly_coeffs(M):
    # det(mu I - M) coefficients high->low, exact
    return ctln.charpoly(M)

def hurwitz_minors(c):
    # c: high->low coefficients c_0..c_n ; returns Delta_1..Delta_n exactly
    n = len(c) - 1
    def a(k): return c[k] if 0 <= k <= n else F(0)
    out = []
    for k in range(1, n + 1):
        H = [[a(2 * (i + 1) - (j + 1)) for j in range(k)] for i in range(k)]
        out.append(ctln.det_exact(H))
    return out

def main(s, beta, t_lo, t_hi):
    A = from_d6(s); m = len(A); sigma = list(range(m))
    print("graph", s, " adjacency A[i][j]=1 iff i->j:")
    for r in A: print("   ", r)
    print("domination-free:", graphs.is_domination_free(A), " complete bidirected:", graphs.is_complete_bidirected(A))
    out = {"d6": s, "adjacency": A, "beta": str(beta)}
    for lab, t in (("t_lo", t_lo), ("t_hi", t_hi)):
        eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        rec = ctln.analyze_support(A, sigma, eps, delta)
        print(f"{lab}: t={t} eps={eps} delta={delta} q={ctln.q_of(eps, delta)} legal={ctln.legal(eps, delta)} "
              f"positive={rec['positive']} stable={rec['stable']} ({rec['stability_detail']}) SFP={rec['stable_fixed_point']}")
        out[lab] = {"t": str(t), "eps": str(eps), "delta": str(delta), "legal": ctln.legal(eps, delta),
                    "positive": rec["positive"], "stable": rec["stable"], "SFP": rec["stable_fixed_point"],
                    "x": rec["x"], "charpoly": rec["charpoly"]}
    # omega at beta, exact
    N = [[(1 if (i != j and A[j][i] == 0) else 0) for j in range(m)] for i in range(m)]
    C = [[F(int(i == j)) + beta * N[i][j] for j in range(m)] for i in range(m)]
    negC = [[-C[i][j] for j in range(m)] for i in range(m)]
    negCJ = [[-C[i][j] - 1 for j in range(m)] for i in range(m)]
    phi = charpoly_coeffs(negC)       # det(mu I + C)
    phiJ = charpoly_coeffs(negCJ)     # det(mu I + C + J)
    omega = [phiJ[k + 1] - phi[k + 1] for k in range(m)]   # degree m-1, high->low
    assert phiJ[0] == phi[0] == 1 and omega[0] == m
    Hm = hurwitz_minors(omega)
    hs = ctln.hurwitz_stable(omega)
    print("omega coeffs (high->low):", [str(c) for c in omega])
    print("Hurwitz minors of omega signs:", [(1 if h > 0 else (-1 if h < 0 else 0)) for h in Hm], " hurwitz_stable:", hs)
    print("1 - 1/beta =", 1 - 1 / beta)
    out["omega"] = [str(c) for c in omega]
    out["omega_hurwitz_minor_signs"] = [(1 if h > 0 else (-1 if h < 0 else 0)) for h in Hm]
    out["omega_hurwitz_stable"] = hs[0]
    return out

if __name__ == '__main__':
    s = sys.argv[1]; beta = F(sys.argv[2]); t_lo = F(sys.argv[3]); t_hi = F(sys.argv[4])
    out = main(s, beta, t_lo, t_hi)
    if len(sys.argv) > 5:
        with open(sys.argv[5], 'w') as fh: json.dump(out, fh, indent=1)

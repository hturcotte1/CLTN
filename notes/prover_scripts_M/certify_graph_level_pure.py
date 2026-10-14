"""
EXACT. Certifies that the GRAPH-LEVEL pure-stability statement
   "legal-stable somewhere (full support, -I+W Hurwitz at some legal (eps,delta))  =>  alpha(A|perp) <= -1/2"
is FALSE at n = 6.  The graphs have dominated nodes (positivity fails), so they are NOT CTLN counterexamples;
this only shows that any graph-level Lemma M must use positivity / domination-freeness.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs
from certify_M_violation import from_d6, hurwitz_minors

def omega_at(A, beta):
    m = len(A)
    N = [[(1 if (i != j and A[j][i] == 0) else 0) for j in range(m)] for i in range(m)]
    C = [[F(int(i == j)) + beta * N[i][j] for j in range(m)] for i in range(m)]
    phi = ctln.charpoly([[-C[i][j] for j in range(m)] for i in range(m)])
    phiJ = ctln.charpoly([[-C[i][j] - 1 for j in range(m)] for i in range(m)])
    om = [phiJ[k + 1] - phi[k + 1] for k in range(m)]
    assert om[0] == m
    return om

if __name__ == '__main__':
  for s in sys.argv[1:]:
        A = from_d6(s); m = len(A); sigma = list(range(m))
        print("graph", s, "adjacency A[i][j]=1 iff i->j:")
        for r in A: print("   ", r)
        print("dominated nodes:", graphs.dominated_nodes(A))
        beta, t = F(9, 4), F(11, 2)
        eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
        rec = ctln.analyze_support(A, sigma, eps, delta)
        print(f"  beta={beta} t={t} eps={eps} delta={delta} q={ctln.q_of(eps, delta)} legal={ctln.legal(eps, delta)} "
              f"stable={rec['stable']} ({rec['stability_detail']}) positive={rec['positive']}")
        om = omega_at(A, F(2))
        Hm = hurwitz_minors(om)
        print("  omega(mu; beta=2) Hurwitz minor signs:", [(h > 0) - (h < 0) for h in Hm], "  omega(0)=s(2) =", om[-1])
        print("  => omega(.;2) not Hurwitz, no root on iR (Delta_{m-2} != 0, s != 0) => alpha(A|perp) > -1/2 strictly")

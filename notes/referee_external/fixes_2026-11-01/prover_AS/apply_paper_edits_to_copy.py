# Applies the proposed paper.tex edits to a COPY: python3 apply_paper_edits_to_copy.py <copy/of/paper.tex> <out.json>. Test build of the edited copy: pdflatex+bibtex+pdflatex x2, 33 pages, no errors.
import sys, json
p = sys.argv[1]
s = open(p).read()
E = []
E.append((r"""Assume $\nu_\infty:=-\alpha\in(1/2,1)$ and that every root $z_0$ of $\psi$ with $\Rea z_0=\nu_\infty$ is simple and not a root of
$\chi$. Let $R^*=\min\Rea\bigl(1/\hat s\,'(z_0)\bigr)$ over these roots. Then $\nu_A(c)=\nu_\infty+R^*/c+O(c^{-2})$, and the
boundary of the stable region near $(\beta_c,\infty)$ is $\beta=\beta_c+\beta_c^3R^*/t+O(t^{-2})$.
\begin{enumerate}
\item[(a)] If $R^*>0$, then for every $\beta$ slightly above $\beta_c$ the stable $t$-set is nonempty and bounded above; if
moreover $a(\beta_c)>0$, these points are SFPs.
\item[(b)] If $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\beta_c$.
\end{enumerate}""",
r"""Assume $\nu_\infty:=-\alpha\in[1/2,1)$, so $\beta_c=1/(1-\nu_\infty)\in[2,\infty)$, and that every root $z_0$ of $\psi$ with
$\Rea z_0=\nu_\infty$ is simple and not a root of $\chi$. Let $R^*=\min\Rea\bigl(1/\hat s\,'(z_0)\bigr)$ over these roots. Then
$\nu_A(c)=\nu_\infty+R^*/c+O(c^{-2})$, and near $(\beta_c,\infty)$ the curve $\nu_A(c)=\gamma$, on which $K$ loses positive
stability, is $\beta=\beta_c+\beta_c^3R^*/t+O(t^{-2})$; for $\nu_\infty>1/2$ it is the boundary of the stable region there.
\begin{enumerate}
\item[(a)] If $\nu_\infty>1/2$ and $R^*>0$, then for every $\beta$ slightly above $\beta_c$ the stable $t$-set is nonempty and
bounded above; if moreover $a(\beta_c)>0$, these points are SFPs.
\item[(b)] If $\nu_\infty>1/2$ and $R^*<0$, the stable $t$-set is upward closed for large $t$ near $\beta_c$.
\item[(c)] If $\nu_\infty=1/2$, then $\beta_c=2$ and the legality boundary $\beta=2+1/t$ is of the same order as the curve
above. If $R^*>1/8$, there are legal stable points with $t\to\infty$ and $\beta\to2^+$, and they are SFPs if $a(2)>0$
(Lemma~\ref{lem:Bmech}); if $R^*<1/8$, no legal point with $\beta$ bounded and $t$ large is stable.
\end{enumerate}"""))
E.append((r"""theorem gives an analytic branch $z(e)=z_0+e\,\chi_1(z_0)/\psi_1'(z_0)+O(e^2)$, and $\chi_1(z_0)/\psi_1'(z_0)=1/\hat s\,'(z_0)$.
\cf{notes/lemma_M_refuted.md}, \S4.""",
r"""theorem gives an analytic branch $z(e)=z_0+e\,\chi_1(z_0)/\psi_1'(z_0)+O(e^2)$, and $\chi_1(z_0)/\psi_1'(z_0)=1/\hat s\,'(z_0)$.
The other bounded eigenvalues stay near roots of $\psi$ with real part $>\nu_\infty$ and one eigenvalue tends to $+\infty$, so
$\nu_A(c)$ is the minimum of $\Rea z(1/c)$ over the critical branches. None of this uses $\nu_\infty>1/2$; that bound enters
only (a) and (b), where it makes every large $c$ legal for all $\gamma$ near $\nu_\infty$. For (c), legal stability means
$c-\sqrt{c^2-c}<\gamma<\nu_A(c)$, and $c-\sqrt{c^2-c}=1/(1+\sqrt{1-1/c})=\tfrac12+\tfrac1{8c}+O(c^{-2})$ exceeds
$\tfrac12+\tfrac1{8c}$ for every $c>1$. \cf{notes/lemma_M_refuted.md}, \S4."""))
E.append((r"""By Theorem~\ref{thm:AS}, $\nu_A(c)=\tfrac12+R/c+O(c^{-2})$, while the legality boundary is""",
r"""By (i) and (ii), $\nu_\infty=1/2$ and the roots of $\psi$ on $\Rea z=1/2$ are exactly $z_0,\bar z_0$, simple and not roots
of $\chi$; at a root of $\psi$, $1/\hat s\,'(z_0)=\chi(z_0)/\psi'(z_0)$, and conjugate roots give the same real part, so
$R^*=R$. The expansion of Theorem~\ref{thm:AS}, which holds at $\nu_\infty=1/2$ (parts (a) and (b) are not used), gives
$\nu_A(c)=\tfrac12+R/c+O(c^{-2})$, while the legality boundary is"""))
for old, new in E:
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)
open(p, 'w').write(s)
json.dump([{"old": o, "new": n} for o, n in E], open(sys.argv[2], 'w'), indent=1)
print("applied", len(E))

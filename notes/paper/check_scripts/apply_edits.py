import sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
E = []
def e(old, new): E.append((old, new))
# E1 abstract: priority-neutral wording
e(r"""and a
new ray reduction to the legality curve""", r"""and a
ray reduction to the legality curve""")
# E2 intro: verification scope
e(r"""Graphs on $n\le7$ nodes reported as counterexamples have been checked by hand and by
four independently written exact verifiers;""", r"""Graphs on $n\le7$ nodes reported as counterexamples have been checked by hand (the two new
six-node graphs also by independent adversarial verification workflows, the $124$ seven-node graphs also by
four independently written exact verifiers);""")
e(r"""all $n\le7$ counterexamples were
hand-verified and re-decided by four independently written exact verifiers;""", r"""all $n\le7$ counterexamples were
hand-verified, and the $124$ seven-node ones were re-decided by four independently written exact verifiers;""")
# E3 Lemma A.2: Lean scope
e(r"""argues instead that $x\mapsto\det(xI+K)$ is monic with no root in $[0,\infty)$.
\end{proof}""", r"""argues instead that $x\mapsto\det(xI+K)$ is monic with no root in $[0,\infty)$.
\end{proof}
The Lean development covers $\det K>0$ and the positivity clause ($x>0\iff a>0$ under stability); the off-support
clause $o_k<0$ is not formalised.""")
# E4 Lemma I.3/I.4 dependence
e(r"""statement in Section~\ref{sec:families}; it is not formalised in Lean.""", r"""statement in Section~\ref{sec:families}; it is not formalised in Lean. Parts (b) and (c) are Lemmas~I.3 and~I.4 of
the project notes; they were refereed when the $n=6$ map was audited and were not re-derived by the referees of the
family theorems.""")
# E5 k<=4 citations
e(r"""negative exactly at legal points), $k=3$ by Theorem~4.1 of \cite{Geneson2026}, $k=4$ by a result of Curto, Geneson and Morrison as
summarised in \cite{Geneson2026}, and all three cases are re-certified by our two-parameter engine on all
digraphs with at most four nodes.""", r"""negative exactly at legal points; a non-adjacent pair is excluded by Lemma~\ref{lem:twin}), $k=3$ by Theorem~4.1 of
\cite{Geneson2026}, $k=4$ is stated as known in \cite{Geneson2026} (for nondegenerate CTLNs, citing work of Curto, Geneson and
Morrison whose exact hypotheses we have not checked), and all three cases are re-certified by our two-parameter engine on all
digraphs with at most four nodes; our claims for $k\le4$ rest on that certification.""")
# E6 all five prunes are Lean-checked
e(r"""Every prune below is a theorem valid at \emph{every} legal parameter; four of the five are machine-checked.""",
  r"""Every prune below is a theorem valid at \emph{every} legal parameter, and all five are machine-checked. The twin
obstruction is not a separate stage of the $n=7$ chain.""")
# E7 CPU time attribution
e(r"""Stage 5 took $16.4$ CPU-hours ($6.6$ ms per line).""",
  r"""The single run performing stages 4 and 5 over the $8{,}909{,}123$ positivity survivors took $16.4$ CPU-hours ($6.6$ ms per line).""")
# E8 illegal first witnesses in the limit file
e(r"""\footnote{\cf{certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl}, \cf{certificates/n7_limit_tally.json}.}""",
  r"""\footnote{\cf{certificates/COUNTEREXAMPLE_CANDIDATES_n7_limit.jsonl}, \cf{certificates/n7_limit_tally.json}. For $15$ graphs the first
stored pointwise witness in the limit file is an illegal point ($q\le1$); use \cf{certificates/n7_candidates_legal_witnesses.json}.}""")
# E9 FNFL cell
e(r"""on the cell $2<\beta<27/13$, at the sample $\beta=53/26$, the stable range is""",
  r"""on its single positivity cell $2<\beta<\beta_1'$ ($\beta_1'\approx2.08302$, the quintic root of Table~\ref{tab:n6},
isolated here in $(27/13,25/12)$), at the sample $\beta=53/26$, the stable range is""")
e(r"""about $1.0414$, reached at the
corner""", r"""about $1.0414$, approached at the
corner""")
# E10 supremum wording
e(r"""and other graphs attain the maximum.""", r"""and other graphs approach the supremum.""")
# E11 section scope / Lemma I dependence
e(r"""All statements in this section are proved for every member of a family; instances on $n\ge8$ nodes are certified by
the pointwise exact verifier but have \emph{not} been hand-verified, and are reported as candidates in that sense. The
proofs use the base graphs only as explicit matrices. Every ``for all $k$'' statement rests on Lemma~\ref{lem:limit}
(transfer from the limit to finite $t$).""",
  r"""The theorems in this section are proved for every member of a family; tables and individual instances are certified
only. Instances on $n\ge8$ nodes are certified by the pointwise exact verifier but have \emph{not} been hand-verified, and
are reported as candidates in that sense. The proofs use the base graphs only as explicit matrices. Theorems~\ref{thm:famF},
\ref{thm:famE}, \ref{thm:QQ2Y} and~\ref{thm:U} rest on Lemma~\ref{lem:limit}(b),(c) (Lemmas~I.3 and~I.4 of the notes:
transfer from the limit to finite $t$, and $t$-free positivity); Theorem~\ref{thm:A} uses only part~(c).""")
# E12 E_1 vs E_0 region: not identical
e(r"""The two-parameter engine applied to $E_1$ returns exactly the stable cell of $E_0$, as Lemma~\ref{lem:W} predicts; so the whole
family shares one region.\footnote{\cf{certificates/n8_instances_region2d.json}.}""",
  r"""The two-parameter engine applied to $E_1$ returns one stable cell, over the same $\beta$-cell and with the upper $t$-end at
$\beta=53/26$ in the same isolating bracket $(53,54)$ as for $E_0$.\footnote{\cf{certificates/n8_instances_region2d.json}.} The two
regions are not identical: at $\beta=53/26$, $t=107/2$ the full support of $E_0$ is stable and that of $E_1$ is not (Hurwitz
determinant~7 negative), so the upper end of the window depends on $a$.""")
# E13 the 15 eight-node graphs: float-steered scan
e(r"""Among all one-node extensions of the $124$
seven-node graphs, exactly $15$ eight-node graphs satisfy every hypothesis of Lemma~\ref{lem:Bmech}, each with an exact legal
witness.""", r"""A float-steered scan of all one-node extensions of the
$124$ seven-node graphs, followed by exact triage of the $9{,}808$ extensions it flagged (an eigenvalue of $\perpc{A}$ near a
primitive cube root of unity), found $15$ eight-node graphs that satisfy every hypothesis of Lemma~\ref{lem:Bmech}, each with an
exact legal witness; extensions not flagged by the float prefilter were not examined exactly.""")
# E14 Theorem A extensions: referee scope
e(r"""Also proved and refereed:""", r"""Also proved in the notes (the referee spot-checked the first statement at $28$ exact legal points and
re-proved the self-join statement independently):""")
# E15 referee assumption
e(r"""were each refereed by a separate reviewer who re-derived the steps; no blocking issue was found.""",
  r"""were each refereed by a separate reviewer who re-derived the steps, taking Lemma~\ref{lem:limit}(b),(c) as given; no blocking
issue was found.""")
# E16 human verification confidence
e(r"""and all $124$ seven-node counterexamples on
15 October 2026. The latter were additionally""", r"""and all $124$ seven-node counterexamples on
15 October 2026; the project owner recorded the latter hand check as not fully confident and made sign-off conditional on
the machine re-verification below. The $124$ were additionally""")
# E17 n=8 open problem
e(r"""and one eight-node limit-stable graph, \dsix+&GLjlNZ}NzH\o+ ($41$ edges), is
already known from a random sample.""", r"""and besides the family members of Section~\ref{sec:families} one further eight-node
limit-stable graph, \dsix+&GLjlNZ}NzH\o+ ($41$ edges; certified, not hand-verified), was found in a random dense sample.""")
# E18 proper supports update
e(r"""realising a counterexample needs condition~(2): describe the admissible in-sets for the $n\le7$ counterexamples.""",
  r"""realising a counterexample needs condition~(2): describe the admissible in-sets for the $n\le7$ counterexamples. At $n=7$ all
$12{,}288$ pairs formed by a six-node counterexample and a seventh node with any in/out pattern have been certified
($7{,}744$ stable somewhere, $4{,}544$ never, $0$ escalations; \cf{certificates/n7_proper_supports_all.jsonl}; not
hand-verified); a structural description is open.""")
# E19 appendix: n=8 record source
e(r"""\cf{family_E_witnesses.json}, \cf{cert_q_910.json}, \cf{cert_q_11.json}, \cf{families_catalogue_n{6,7}.json} & \cf{family_*_proof.py}, \cf{certify_q.py},""",
  r"""\cf{family_E_witnesses.json}, \cf{top_n8.json}, \cf{cert_q_910.json}, \cf{cert_q_11.json}, \cf{families_catalogue_n{6,7}.json} & \cf{family_*_proof.py}, \cf{certify_top.py}, \cf{certify_q.py},""")
for i,(o,n) in enumerate(E):
    c = s.count(o)
    if c != 1:
        print('EDIT', i+1, 'matches', c); sys.exit(1)
    s = s.replace(o, n)
open(dst,'w').write(s)
print('applied', len(E), 'edits')

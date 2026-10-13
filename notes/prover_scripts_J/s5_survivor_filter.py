# Classify the 34,071 'survivors' by which positivity filter they pass, and locate where in the shards
# the filter changes. Exact (fastregion polynomials, Fractions).
exec(open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/proverJ/common.py').read())
import collections
surv=[int(l) for l in open('/home/user/CLTN/certificates/n6_survivors.idx') if l.strip()]
strict=0; loose_only=0; examples=[]
for g in surv:
    A=graph_of_idx(g)
    r=fastregion.region(A,list(range(6)),None,positivity_only=True)
    if r['verdict']=='positivity_possible': strict+=1
    else:
        loose_only+=1
        if len(examples)<5: examples.append(g)
print('survivors passing current (a_i>0) filter:',strict,' passing only the common-sign filter:',loose_only, examples)
# where do the records come from: per shard, the max idx of a loose-only record and the min idx of a 'never'-without-cert record that the loose filter would keep

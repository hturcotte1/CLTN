import json, collections, statistics
from fractions import Fraction as F
rs=[json.loads(l) for l in open('results.jsonl')]
idx=[r['idx'] for r in rs]; assert sorted(idx)==list(range(9133)) and len(set(idx))==9133
cw=json.load(open('/home/user/CLTN/certificates/n8_candidates_legal_witnesses.json'))['graphs']
nulls=set(i for i,g in enumerate(cw) if g['witness_stored_legal'] is None)
C=collections.Counter
T={}
T['n']=len(rs)
T['with_problems']=sum(1 for r in rs if r['problems'])
T['errors']=sum(1 for r in rs if any(p.startswith('ERROR') for p in r['problems']))
T['W_definition_check_true']=sum(1 for r in rs if r.get('W_definition_check'))
T['complete_graphs']=sum(1 for r in rs if r.get('complete'))
T['beta_cell_in_first_stable_interval']=sum(1 for r in rs if r.get('beta_in_first_stable_interval'))
T['cell_S_contains_all_large_t']=sum(1 for r in rs if r['cell']['full']['S_contains_all_large_t'])
T['cell_stabpos_S_contains_all_large_t']=sum(1 for r in rs if r['cell']['stab_pos_only']['S_contains_all_large_t'])
T['cell_t_in_S_true']=sum(1 for r in rs if r.get('cell_t_in_S') is True)
T['cell_pointwise_stable_positive_legal']=sum(1 for r in rs if r['cell_pointwise']['stable'] is True and r['cell_pointwise']['positive'] is True and r['cell_pointwise']['legal'])
T['cell_crosscheck_all_consistent']=sum(1 for r in rs if all(a['hurwitz_signs_agree'] and a['x_signs_agree'] and a['consistent'] for a in r['cell_crosscheck']))
T['cell_crosscheck_points']=sum(len(r['cell_crosscheck']) for r in rs)
T['cell_crosscheck_undecided_points']=sum(1 for r in rs for a in r['cell_crosscheck'] if a['in_S'] is None)
T['stored_witness_points']=sum(len(r['stored_pointwise']) for r in rs)
T['stored_witness_points_agree']=sum(1 for r in rs for s in r['stored_pointwise'] if s['agree'])
T['stored_witness_points_my_undecided']=sum(1 for r in rs for s in r['stored_pointwise'] if s['my_stable'] is None or s['my_positive'] is None)
T['stored_beta_analyses']=sum(len(r['stored']) for r in rs)
T['stored_beta_checks']=sum(len(v['witness_checks']) for r in rs for v in r['stored'].values())
T['stored_beta_checks_consistent']=sum(1 for r in rs for v in r['stored'].values() for c in v['witness_checks'] if c['consistent'])
T['stored_beta_checks_undecided']=sum(1 for r in rs for v in r['stored'].values() for c in v['witness_checks'] if c['in_S_stab_pos'] is None)
T['stored_beta_in_interval']=sum(1 for r in rs for v in r['stored'].values() if v['beta_in_first_stable_interval'])
T['stored_beta_unbounded']=sum(1 for r in rs for v in r['stored'].values() if v['stab_pos_set_contains_all_large_t'])
T['stored_beta_expectation_met']=sum(1 for r in rs for v in r['stored'].values() if v['expected_from_limit_interval'] is None or v['expected_from_limit_interval']==v['stab_pos_set_contains_all_large_t'])
T['stored_beta_expectation_undetermined']=sum(1 for r in rs for v in r['stored'].values() if v['expected_from_limit_interval'] is None)
# the 37 nulls
T['nulls']=len(nulls)
T['nulls_cell_unbounded']=sum(1 for r in rs if r['idx'] in nulls and r['cell']['full']['S_contains_all_large_t'])
T['nulls_stored_beta_bounded']=sum(1 for r in rs if r['idx'] in nulls and all(not v['stab_pos_set_contains_all_large_t'] for v in r['stored'].values()))
T['nulls_stored_beta_outside_interval']=sum(1 for r in rs if r['idx'] in nulls and all(not v['beta_in_first_stable_interval'] for v in r['stored'].values()))
# t* statistics at beta_cell (full set incl. legality)
tstar=[]; binding=C(); deg=C(); ncomp=C(); nroots=C()
tstar_sp=[]
for r in rs:
    full=r['cell']['full']
    ts=full['t_star']; assert ts is not None
    hi=F(ts['hi']); lo=F(ts['lo']); tstar.append((hi,r['idx'],lo))
    van=[rt for rt in full['roots'] if rt['lo']==ts['lo'] and rt['hi']==ts['hi']][0]['vanishing']
    binding[tuple(van)]+=1
    deg[len(full['factors_hi_to_lo'][ts['factor_index']])-1]+=1
    ncomp[len(full['S'])]+=1; nroots[full['n_roots']]+=1
    sp=r['cell']['stab_pos_only']
    tstar_sp.append((F(sp['t_star']['hi']) if sp['t_star'] else F(0), r['idx']))
T['t_star_binding_condition_counts']={'+'.join(k):v for k,v in binding.items()}
T['t_star_factor_degree_counts']=dict(deg)
T['S_number_of_components']=dict(ncomp)
T['n_positive_roots_counts']=dict(sorted(nroots.items()))
tstar.sort(); T['t_star_max']={'idx':tstar[-1][1],'lo':str(tstar[-1][2]),'hi':str(tstar[-1][0]),'float':float(tstar[-1][0])}
T['t_star_min']={'idx':tstar[0][1],'lo':str(tstar[0][2]),'hi':str(tstar[0][0]),'float':float(tstar[0][0])}
T['t_star_median_float']=float(tstar[len(tstar)//2][0])
T['t_star_le_100']=sum(1 for h,_,_ in tstar if h<100); T['t_star_le_1000']=sum(1 for h,_,_ in tstar if h<1000); T['t_star_le_10000']=sum(1 for h,_,_ in tstar if h<10000); T['t_star_le_100000']=sum(1 for h,_,_ in tstar if h<100000)
tstar_sp.sort(); T['stabpos_t_star_max']={'idx':tstar_sp[-1][1],'float':float(tstar_sp[-1][0])}
T['stabpos_t_star_le_100']=sum(1 for h,_ in tstar_sp if h<100)
# legality threshold vs stab/pos threshold: binding = legality?
T['binding_is_legality_only']=sum(1 for r in rs if [rt for rt in r['cell']['full']['roots'] if rt['lo']==r['cell']['full']['t_star']['lo']][0]['vanishing']==['L'])
T['seconds_total']=round(sum(r['seconds'] for r in rs),1); T['seconds_max']=max(r['seconds'] for r in rs); T['seconds_mean']=round(statistics.mean(r['seconds'] for r in rs),3)
T['degrees_check']=sum(1 for r in rs if r['cell']['degrees']=={'H1':1,'H2':2,'H3':3,'H4':4,'H5':5,'H6':6,'H7':7,'D':1,'N0':0,'N1':0,'N2':0,'N3':0,'N4':0,'N5':0,'N6':0,'N7':0,'L':1})
json.dump(T,open('tally.json','w'),indent=1)
print(json.dumps(T,indent=1))
probs=[(r['idx'],r['problems']) for r in rs if r['problems']]
print('PROBLEMS', probs[:20])
# per-graph table
with open('table.tsv','w') as f:
    f.write('idx\td6\tedges\tbeta_cell\tbeta_in_interval\tt_legal\tt_star_lo\tt_star_hi\tt_star_binding\tS_unbounded\tn_roots\tS_components\tcell_t\tcell_t_in_S\tcrosscheck_ok\tbeta_stored\tstored_pointwise_agree\tstored_unbounded\tstored_expect\tproblems\n')
    for r in rs:
        full=r['cell']['full']; ts=full['t_star']
        van=[rt for rt in full['roots'] if rt['lo']==ts['lo'] and rt['hi']==ts['hi']][0]['vanishing']
        bs=r['beta_stored'][0]; sv=r['stored'][bs]
        f.write('\t'.join(map(str,[r['idx'],r['d6'],r['edges'],r['beta_cell'],r['beta_in_first_stable_interval'],r['cell']['legal_threshold_t'],ts['lo'],ts['hi'],'+'.join(van),full['S_contains_all_large_t'],full['n_roots'],len(full['S']),
            cw[r['idx']]['witness_in_first_stable_cell']['t'],r['cell_t_in_S'],all(a['hurwitz_signs_agree'] and a['x_signs_agree'] and a['consistent'] for a in r['cell_crosscheck']),bs,all(s['agree'] for s in r['stored_pointwise']),sv['stab_pos_set_contains_all_large_t'],sv['expected_from_limit_interval'],';'.join(r['problems'])]))+'\n')

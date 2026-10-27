#!/bin/bash
# after_prefilter2.sh -- when every n = 8 job is done: close the accounting chain, run the exact limit stage on
# the S/O lines of shards 0..8 and of all jobs, and commit.  Replaces after_prefilter.sh (geng-shard scheme).
cd "$(dirname "$0")/.." || exit 1
until grep -q ALL_JOBS_DONE n8/jobs/progress.log 2>/dev/null; do sleep 600; done
if ! python3 n8/check_accounting.py > n8/check_accounting.log 2>&1; then
  echo "$(date -u +%H:%M:%S) ACCOUNTING_MISMATCH (see n8/check_accounting.log)" >> n8/jobs/progress.log; exit 1
fi
# the per-B ledgers of the two enumeration swaps (K8 batches: n8/k8decomp; orbitenum parts: n8/orbitenum) must close too
for c in n8/k8decomp/check_k8decomp.py n8/orbitenum/check_orbit.py; do
  if ! python3 "$c" > "n8/$(basename "$c" .py).log" 2>&1; then
    echo "$(date -u +%H:%M:%S) LEDGER_MISMATCH ($c, see n8/$(basename "$c" .py).log)" >> n8/jobs/progress.log; exit 1
  fi
done
python3 n8/run_n8_limit.py --parts n8/parts n8/jobs --out n8 > n8/run_n8_limit.log 2>&1
cat n8/parts/*.stats n8/jobs/*.stats | awk '{for(i=1;i<=NF;i++){if($i=="total")t+=$(i+1); if($i=="domfree")d+=$(i+1); if($i=="survivors")s+=$(i+1); if($i=="boundary")b+=$(i+1); if($i=="overflow")o+=$(i+1)}} END {print "total",t,"domfree",d,"survivors",s,"boundary",b,"overflow",o}' > n8/n8_prefilter_TOTAL.txt
git add n8/n8_limit_results.jsonl n8/n8_limit_HITS.jsonl n8/n8_limit_tally.json n8/n8_prefilter_TOTAL.txt n8/n8_accounting.json n8/parts/*.done n8/parts/*.stats n8/parts/*.txt n8/jobs/*.done n8/jobs/*.stats n8/jobs/*.txt n8/jobs/*.genlog 2>/dev/null
if python3 test_ctln.py 2>&1 | grep -q "20 passed, 0 failed"; then
  git commit -q -m "n=8 limit stage complete (local run): see n8/n8_limit_tally.json and n8/n8_limit_HITS.jsonl (candidates)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W9KhFwNapBjm91tofywDXi" && git push -q -u origin claude/compassionate-maxwell-mtoyeb
fi
echo "$(date -u +%H:%M:%S) AFTER_PREFILTER_DONE" >> n8/jobs/progress.log

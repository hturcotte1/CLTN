#!/bin/bash
# after_prefilter.sh -- when every n = 8 shard is done, run the exact limit stage on the S/O lines and commit.
cd "$(dirname "$0")/.." || exit 1
until grep -q ALL_SHARDS_DONE n8/parts/progress.log 2>/dev/null; do sleep 600; done
python3 n8/run_n8_limit.py --parts n8/parts --out n8 > n8/run_n8_limit.log 2>&1
cat n8/parts/*.stats | awk '{for(i=1;i<=NF;i++){if($i=="total")t+=$(i+1); if($i=="domfree")d+=$(i+1); if($i=="survivors")s+=$(i+1); if($i=="boundary")b+=$(i+1); if($i=="overflow")o+=$(i+1)}} END {print "total",t,"domfree",d,"survivors",s,"boundary",b,"overflow",o}' > n8/n8_prefilter_TOTAL.txt
git add n8/n8_limit_results.jsonl n8/n8_limit_HITS.jsonl n8/n8_limit_tally.json n8/n8_prefilter_TOTAL.txt n8/parts/*.done n8/parts/*.stats n8/parts/*.txt 2>/dev/null
if python3 test_ctln.py 2>&1 | grep -q "20 passed, 0 failed"; then
  git commit -q -m "n=8 limit stage complete (local run): see n8/n8_limit_tally.json and n8/n8_limit_HITS.jsonl (candidates)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W9KhFwNapBjm91tofywDXi" && git push -q -u origin claude/compassionate-maxwell-mtoyeb
fi
echo AFTER_PREFILTER_DONE >> n8/parts/progress.log

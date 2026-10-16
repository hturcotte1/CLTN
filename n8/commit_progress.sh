#!/bin/bash
# commit_progress.sh -- every 30 minutes, commit newly completed n = 8 outputs (shards 0..8 in n8/parts, jobs in
# n8/jobs: S/O lines, stats, done markers) so progress survives a container loss.  Runs test_ctln.py first.
cd "$(dirname "$0")/.." || exit 1
while true; do
  sleep 1800
  n_done=$(ls n8/jobs/job*.done 2>/dev/null | wc -l); n_jobs=$(($(wc -l < n8/jobs/jobs.tsv) - 1))
  pct=$(cat n8/parts/*.stats n8/jobs/*.stats 2>/dev/null | awk '{for(i=1;i<=NF;i++) if($i=="total") t+=$(i+1)} END {printf "%.2f", 100*t/801558826511}')
  git add n8/parts/*.done n8/parts/*.stats n8/parts/*.txt n8/jobs/*.done n8/jobs/*.stats n8/jobs/*.txt 2>/dev/null
  if ! git diff --cached --quiet; then
    if python3 test_ctln.py 2>&1 | grep -q "20 passed, 0 failed"; then
      git commit -q -m "n=8 local prefilter progress: $n_done of $n_jobs jobs done ($pct% of the digraphs)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W9KhFwNapBjm91tofywDXi" && git push -q -u origin claude/compassionate-maxwell-mtoyeb
    fi
  fi
done

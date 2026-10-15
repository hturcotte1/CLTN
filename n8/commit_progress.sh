#!/bin/bash
# commit_progress.sh -- every 30 minutes, commit newly completed n = 8 shard outputs (S/O lines, stats, done
# markers) so progress survives a container loss. Runs test_ctln.py before each commit (project rule).
cd "$(dirname "$0")/.." || exit 1
while true; do
  sleep 1800
  n_done=$(ls n8/parts/*.done 2>/dev/null | wc -l)
  git add n8/parts/*.done n8/parts/*.stats n8/parts/*.txt 2>/dev/null
  if ! git diff --cached --quiet; then
    if python3 test_ctln.py 2>&1 | grep -q "20 passed, 0 failed"; then
      git commit -q -m "n=8 local prefilter progress: $n_done of 4000 shards done

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01W9KhFwNapBjm91tofywDXi" && git push -q -u origin claude/compassionate-maxwell-mtoyeb
    fi
  fi
done

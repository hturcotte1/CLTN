#!/bin/bash
# run_n8_jobs.sh [W] -- n = 8 limit prefilter over the (graph, directg part) jobs in n8/jobs/jobs.tsv
# (built by n8/make_jobs.py; replaces the geng res/mod shards of run_n8_local.sh beyond shard 8).
# Each job:  job<j>.g6 | directg -q -e29:56 [-s p/P] | filter8  -> job<j>.stats (counts), job<j>.txt (S/O lines;
# B lines dropped: a vanishing Hurwitz minor at alpha = -1/2 means alpha >= -1/2, never limit-stable, exact).
# Resumable (job<j>.done), largest jobs first, at most W filter8 pipelines at once (counts ALL filter8 processes,
# so a second copy of this script cannot overload the machine).  A pipeline failure or a count mismatch leaves no
# .done (job is retried at the next relaunch) and is logged in progress.log.
cd "$(dirname "$0")/.." || exit 1
W=${1:-4}; J=n8/jobs; mkdir -p "$J"; touch "$J/progress.log"
run_job() { local j=$1 split=$2 expected=$3 opt=""
  [ "$split" != "-" ] && opt="-s$split"
  if nice -n 19 bash -c "set -o pipefail; n7/directg -q -e29:56 $opt < $J/job$j.g6 | n7/filter8 2> $J/job$j.stats | awk '\$3!=\"B\"' > $J/job$j.txt.tmp"; then
    if [ "$split" = "-" ]; then
      tot=$(awk '{for(i=1;i<=NF;i++) if($i=="total") print $(i+1)}' "$J/job$j.stats")
      if [ "$tot" != "$expected" ]; then
        mv "$J/job$j.txt.tmp" "$J/job$j.txt.bad"; echo "$(date -u +%H:%M:%S) MISMATCH job $j total=$tot expected=$expected" >> "$J/progress.log"; return 1; fi
    fi
    mv "$J/job$j.txt.tmp" "$J/job$j.txt" && touch "$J/job$j.done"
    echo "$(date -u +%H:%M:%S) finished job $j" >> "$J/progress.log"
  else
    echo "$(date -u +%H:%M:%S) FAILED job $j (pipeline exit)" >> "$J/progress.log"; rm -f "$J/job$j.txt.tmp"; return 1
  fi; }
tail -n +2 "$J/jobs.tsv" | while IFS=$'\t' read -r j expected split graphs; do
  [ -f "$J/job$j.done" ] && continue
  while [ "$(pgrep -xc filter8)" -ge "$W" ]; do sleep 30; done
  awk -F'\t' -v j="$j" '$1==j {n=split($4,a," "); for(i=1;i<=n;i++) print a[i]}' "$J/jobs.tsv" > "$J/job$j.g6"
  run_job "$j" "$split" "$expected" < /dev/null &
  echo "$(date -u +%H:%M:%S) started job $j (expected $expected, split $split)" >> "$J/progress.log"
  sleep 5
done
wait
n_done=$(ls "$J"/job*.done 2>/dev/null | wc -l); n_jobs=$(($(wc -l < "$J/jobs.tsv") - 1))
if [ "$n_done" = "$n_jobs" ]; then echo "$(date -u +%H:%M:%S) ALL_JOBS_DONE" >> "$J/progress.log"
else echo "$(date -u +%H:%M:%S) DRIVER_EXIT with $n_done of $n_jobs jobs done (failed jobs will be retried at relaunch)" >> "$J/progress.log"; fi

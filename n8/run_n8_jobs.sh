#!/bin/bash
# run_n8_jobs.sh [W] -- n = 8 limit prefilter over the (graph, directg part) jobs in n8/jobs/jobs.tsv
# (built by n8/make_jobs.py).  Each job:  job<j>.g6 | directg -q -e29:56 [-s p/P] | filter8  -> job<j>.stats
# (counts), job<j>.txt (S/O lines; B lines dropped: a vanishing Hurwitz minor at alpha = -1/2 means
# alpha >= -1/2, never limit-stable, exact).
# Resumable (job<j>.done), largest jobs first, at most W filter8 pipelines at once (counts ALL filter8
# processes).  A job already being run by another pipeline (its job<j>.g6 open in a directg) is skipped.
# The main loop's PID is in n8/jobs/driver.pid: stop it with n8/stop_driver.sh, which leaves the per-job
# wrappers alive so that in-flight jobs still finalize (mv + .done).  A pipeline failure or a count
# mismatch leaves no .done (retried at the next relaunch) and is logged in progress.log.
cd "$(dirname "$0")/.." || exit 1
W=${1:-4}; J=n8/jobs; mkdir -p "$J"; touch "$J/progress.log"
echo $$ > "$J/driver.pid"; trap 'rm -f "$J/driver.pid"' EXIT
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
while IFS=$'\t' read -r j expected split graphs; do
  [ -f "$J/job$j.done" ] && continue
  if pgrep -f "$J/job$j\\.g6" > /dev/null; then echo "$(date -u +%H:%M:%S) skip job $j (already running)" >> "$J/progress.log"; continue; fi
  while [ "$(pgrep -xc filter8)" -ge "$W" ]; do sleep 30; done
  [ -f "$J/job$j.done" ] && continue
  pgrep -f "$J/job$j\\.g6" > /dev/null && continue
  awk -F'\t' -v j="$j" '$1==j {n=split($4,a," "); for(i=1;i<=n;i++) print a[i]}' "$J/jobs.tsv" > "$J/job$j.g6"
  run_job "$j" "$split" "$expected" < /dev/null &
  echo "$(date -u +%H:%M:%S) started job $j (expected $expected, split $split)" >> "$J/progress.log"
  sleep 5
done < <(tail -n +2 "$J/jobs.tsv")
wait
n_done=$(ls "$J"/job*.done 2>/dev/null | wc -l); n_jobs=$(($(wc -l < "$J/jobs.tsv") - 1))
if [ "$n_done" = "$n_jobs" ]; then echo "$(date -u +%H:%M:%S) ALL_JOBS_DONE" >> "$J/progress.log"
else echo "$(date -u +%H:%M:%S) DRIVER_EXIT with $n_done of $n_jobs jobs done (unfinished jobs are retried at relaunch)" >> "$J/progress.log"; fi

#!/bin/bash
# relaunch.sh -- idempotent: start whichever of the n = 8 job driver, the progress committer and the
# after-prefilter watcher is not running.  Safe to run any time (after a container restart, from the scheduled
# check, by hand).  Prints what it did.
cd "$(dirname "$0")/.." || exit 1
start() { local pat=$1 cmd=$2 arg=$3 log=$4
  if pgrep -f "$pat" > /dev/null; then echo "running: $cmd"
  else (nohup bash "$cmd" $arg > "$log" 2>&1 &); echo "started: $cmd $arg"; fi; }
start "run_n8_""jobs.sh" n8/run_n8_jobs.sh 4 n8/run_n8_jobs.log
start "commit_""progress.sh" n8/commit_progress.sh "" n8/commit_progress.log
start "after_""prefilter2.sh" n8/after_prefilter2.sh "" n8/after_prefilter2.log
echo "filter8 pipelines: $(pgrep -xc filter8); jobs done: $(ls n8/jobs/job*.done 2>/dev/null | wc -l); $(tail -n 1 n8/jobs/progress.log 2>/dev/null)"

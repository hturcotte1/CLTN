#!/bin/bash
# relaunch.sh -- idempotent: start whichever of the n = 8 job driver, the progress committer and the
# after-prefilter watcher is not running.  Safe to run any time (after a container restart, from the scheduled
# check, by hand).  The driver is detected by its pidfile (its job wrappers share its process name).
cd "$(dirname "$0")/.." || exit 1
# the three gitignored binaries the jobs need; rebuild any that a container restore lost (plain C, no dependencies)
[ -x n7/filter8 ] || { gcc -O2 -o n7/filter8 n7/filter8.c && echo "rebuilt n7/filter8"; }
[ -x n8/k8decomp/addarcs ] || { gcc -O2 -o n8/k8decomp/addarcs n8/k8decomp/addarcs.c && echo "rebuilt n8/k8decomp/addarcs"; }
[ -x n8/orbitenum/orbitenum ] || { bash n8/orbitenum/build.sh; }
p=$(cat n8/jobs/driver.pid 2>/dev/null)
if [ -n "$p" ] && kill -0 "$p" 2>/dev/null && tr "\0" " " < "/proc/$p/cmdline" 2>/dev/null | grep -q "run_n8_jobs"; then echo "running: driver (pid $p)"
else rm -f n8/jobs/driver.pid; (nohup bash n8/run_n8_jobs.sh 4 > n8/run_n8_jobs.log 2>&1 &); echo "started: driver"; fi
start() { local pat=$1 cmd=$2 log=$3
  if pgrep -f "$pat" > /dev/null; then echo "running: $cmd"
  else (nohup bash "$cmd" > "$log" 2>&1 &); echo "started: $cmd"; fi; }
start "commit_""progress.sh" n8/commit_progress.sh n8/commit_progress.log
start "after_""prefilter2.sh" n8/after_prefilter2.sh n8/after_prefilter2.log
sleep 1; echo "filter8 pipelines: $(pgrep -xc filter8); jobs done: $(ls n8/jobs/job*.done 2>/dev/null | wc -l); $(tail -n 1 n8/jobs/progress.log 2>/dev/null)"

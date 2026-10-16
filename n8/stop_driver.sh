#!/bin/bash
# stop_driver.sh -- stop the n = 8 job driver's main loop only (in-flight jobs keep running and finalize)
cd "$(dirname "$0")/.." || exit 1
p=$(cat n8/jobs/driver.pid 2>/dev/null)
if [ -n "$p" ] && kill -0 "$p" 2>/dev/null; then kill "$p" && echo "stopped driver loop pid $p"; else echo "no driver loop running"; fi
rm -f n8/jobs/driver.pid

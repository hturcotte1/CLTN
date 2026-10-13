#!/bin/sh
# Resume the n=6 limit sweep shards after an interruption and re-arm the continuation.
cd /home/user/CLTN
S=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad
mkdir -p $S
for s in 0 1 2; do
  if [ ! -f certificates/region_n6_limit_shard${s}of3_tally.json ]; then
    nohup python3 run_region_sweep.py --n 6 --mode limit --shard $s/3 --resume >> $S/limit_shard$s.log 2>&1 &
  fi
done
if [ ! -f $S/after_limit.done ]; then
  nohup ./after_limit.sh >> $S/after_limit.log 2>&1 &
fi
echo "pipeline (re)started"

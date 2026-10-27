#!/bin/bash
# launch.sh -- (re)start the four limit-stage shards; each resumes from its own results file (idempotent)
cd /home/user/CLTN; SH=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/limit_shards
for k in 0 1 2 3; do
  if ps -eo args | grep -q "^python3 n8/run_n8_limit.py --parts $SH/in$k "; then echo "shard $k running"; continue; fi
  if [ -f $SH/out$k/n8_limit_tally.json ]; then echo "shard $k done"; continue; fi
  (nohup python3 n8/run_n8_limit.py --parts $SH/in$k --out $SH/out$k >> $SH/out$k/log 2>&1 &); echo "shard $k started"
done

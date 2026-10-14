#!/bin/bash
cd /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_prefilter_logic
PS=129,130,131,133,136,144,160,176,192
for f in "$@"; do python3 -I check.py $PS --region < $f > $f.chk 2>&1 & done
wait

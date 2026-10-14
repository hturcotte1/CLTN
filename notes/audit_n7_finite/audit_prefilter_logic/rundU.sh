#!/bin/bash
cd /tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_prefilter_logic
ls dUc/c_?? | xargs -P 4 -I{} sh -c 'python3 -I check.py 129,130,131,133,136,144,160,176,192 --region < {} > {}.chk 2>&1'

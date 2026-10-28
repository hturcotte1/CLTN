#!/bin/bash
cd /home/user/CLTN/notes/verify_n8_2/dynamics_float
echo "START $(date)"
nice -n 10 python3 run_cands.py candidates_sample.json candidates_out.jsonl > candidates.log 2>&1
echo "CANDS DONE $(date)"
nice -n 10 python3 run_never.py never_sample.json never_out.jsonl > never.log 2>&1
echo "NEVER DONE $(date)"
nice -n 10 python3 run_anchors.py > anchors.log 2>&1
echo "ANCHORS DONE $(date)"

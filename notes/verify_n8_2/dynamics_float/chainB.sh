#!/bin/bash
cd /home/user/CLTN/notes/verify_n8_2/dynamics_float
echo "START B $(date)"
CAP_FULL=6e5 CAP_LIMITED=3e5 SCALES=0.3,0.003 LONG_KINDS=perturb0.003,uniform nice -n 10 python3 run_cands.py candidates_sample_B.json candidates_out_B.jsonl > candidates_B.log 2>&1
echo "CANDS B DONE $(date)"
nice -n 10 python3 run_never.py never_sample.json never_out.jsonl > never.log 2>&1
echo "NEVER DONE $(date)"
CAP_FULL=3e6 nice -n 10 python3 run_anchors.py > anchors.log 2>&1
echo "ANCHORS DONE $(date)"

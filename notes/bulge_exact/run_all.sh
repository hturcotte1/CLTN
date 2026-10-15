#!/bin/sh
# Re-run everything (about 8 min single core): exact R*, join with the 2-D certificates, exact adjacency and beta_c cross-checks.
D=/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf5/bulge
cd $D && nice python3 rstar_exact.py /home/user/CLTN/certificates/COUNTEREXAMPLE_CANDIDATES_n7_finite.jsonl rstar_raw.json > run.log \
 && nice python3 compare.py && nice python3 adjacency.py | head -1 && nice python3 crosscheck_betac.py

# limit_sharding -- how the n = 8 exact limit stage was finished after the container restart of 2026-10-27:
# four shards of the unmodified n8/run_n8_limit.py over disjoint symlinked subsets of the job files (launch.sh),
# merged and cross-checked by merge_finalize.py (see its docstring and the ANOMALIES.md entry of 2026-10-28);
# shard<k>_tally.json are the four shard tallies, write_summary.py wrote certificates/n8_limit_summary.json.

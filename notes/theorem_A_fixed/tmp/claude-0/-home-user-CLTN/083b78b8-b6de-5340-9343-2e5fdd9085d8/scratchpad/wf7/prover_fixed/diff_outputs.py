"""diff_outputs.py -- compare regenerated JSON outputs with the backed-up originals (orig_outputs/), ignoring timing fields."""
import json, sys, os
def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k not in ('seconds',)}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o
for fn in sys.argv[1:]:
    a = strip(json.load(open(os.path.join('orig_outputs', fn)))); b = strip(json.load(open(fn)))
    print(fn, 'IDENTICAL (modulo timing)' if a == b else 'DIFFERS')

#!/bin/sh
# setup_remote.sh -- prepare a fresh Linux box (Debian/Ubuntu, gcc, python3 >= 3.10) for the
# n = 7 finite-eps pipeline. Run from the repository root. Idempotent.
#   sh n7/setup_remote.sh
set -e
python3 -m pip install --quiet sympy 2>/dev/null || python3 -m pip install --quiet --break-system-packages sympy
# nauty 2.8.8 (geng, directg, labelg) from the pynauty source distribution, as used for every n = 7 count
mkdir -p build && cd build
if [ ! -x nauty/geng ]; then
  python3 -m pip download --no-binary :all: --no-deps -q pynauty==2.8.8.1
  tar xzf pynauty-2.8.8.1.tar.gz
  cd pynauty-2.8.8.1/src/nauty2_8_8 && ./configure --quiet >/dev/null && make -s geng directg labelg showg && cd ../../..
  mkdir -p nauty && cp pynauty-2.8.8.1/src/nauty2_8_8/geng pynauty-2.8.8.1/src/nauty2_8_8/directg \
     pynauty-2.8.8.1/src/nauty2_8_8/labelg pynauty-2.8.8.1/src/nauty2_8_8/showg nauty/
fi
cd ..
gcc -O2 -o n7/filter7p n7/filter7p.c
for b in geng directg labelg showg; do ln -sf "$(pwd)/build/nauty/$b" "n7/$b"; done
# sanity: the enumerator must reproduce A000273 for n <= 6 and the filter must keep Geneson's graph
test "$(n7/geng -q 6 | n7/directg -q -u 2>&1 | grep -o '[0-9]* digraphs' | cut -d' ' -f1)" = "1540944" || { echo "A000273(6) check FAILED"; exit 1; }
n7/filter7p < n7/geneson.d6 2>/dev/null | grep -q ' P ' || { echo "Geneson check FAILED"; exit 1; }
python3 test_ctln.py | tail -1
echo "setup OK"

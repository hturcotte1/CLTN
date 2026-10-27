#!/bin/bash
# build.sh -- compile the orbit enumerator (plain C, no dependencies)
cd "$(dirname "$0")" && gcc -O2 -Wall -Wextra -o orbitenum orbitenum.c && echo built orbitenum

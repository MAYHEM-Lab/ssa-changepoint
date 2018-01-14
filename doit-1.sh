#!/bin/bash

./ssa-decomp -x sin+N_0_1.txt -l 100 -N 200 -e 16-16 > zzz000
./ssa-decomp -x sin+N_0_1.txt -l 100 -N 200 -e 33-33 > zzz111
./ssa-decomp -x sin+N_0_1.txt -l 100 -N 200 -e 66-66 > zzz222

abut zzz000 zzz111 zzz222 | awk '{print $1+$2+$3}' 

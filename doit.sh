#!/bin/bash

R=100
CNT=2

while ( test $CNT -le $R ) ; do
	./ssa-decomp -x sin+N_0_1.txt -l 100 -N 200 -e $CNT-$R > ggg
	MEAN=`stats mean < ggg`
	SD=`stats sd < ggg`
	echo $CNT-$R $MEAN $SD
	CNT=$(($CNT+1))
done


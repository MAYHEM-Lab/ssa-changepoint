#!/bin/bash

WIN=30
LAST=`wc -l $1 | awk '{print $1}'`
LAST=$(($LAST-$WIN-$WIN))

CNT=1
while ( test $CNT -lt $LAST ) ; do
	LEFT=$(($CNT+$WIN-1))
	RIGHT=$(($LEFT+$WIN-1))

	head -n $LEFT $1 | tail -n $WIN > test.left
	head -n $RIGHT $1 | tail -n $WIN > test.right
	t-test -x test.left -y test.right

	CNT=$(($CNT+1))
done


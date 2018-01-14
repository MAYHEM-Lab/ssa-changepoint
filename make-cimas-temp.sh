#!/bin/bash

# assumes that input file is CIMAS csv

BIN=/Users/rich/bin

CNT=1
SIZE=`wc -l $1 | awk '{print $1}'`

while ( test $CNT -le $SIZE ) ; do
	LINE=`head -n $CNT $1 | tail -n 1`
	DSTR=`echo $LINE | awk -F ',' '{print $1}'`
	TSTR=`echo $LINE | awk -F ',' '{print $3}'`
	TS=`./time-change.sh $DSTR $TSTR | $BIN/convert_time`
	TEMP=`echo $LINE | awk -F ',' '{print $4}'`
	echo $TS $TEMP
	CNT=$(($CNT+1))
done

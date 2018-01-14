#!/bin/bash

# assumes that input file is CIMAS csv

BIN=/Users/rich/bin

CIMAS=$1
MTEMP=$2

CNT=1
SIZE=`wc -l $CIMAS | awk '{print $1}'`
MSIZE=`wc -l $MTEMP | awk '{print $1}'`

while ( test $CNT -le $SIZE ) ; do
	LINE=`head -n $CNT $1 | tail -n 1`
	TS=`echo $LINE | awk '{print $1}'`
	TEMP=`echo $LINE | awk '{print $2}'`
	MIN=9999999999
	MCNT=1
	LASTDIFF=0
	while ( test $MCNT -le $MSIZE ) ; do
		MLINE=`head -n $MCNT $MTEMP | tail -n 1`
		MTS=`echo $MLINE | awk '{print $1}' | sed 's/\./ /' | awk '{print $1}'`
		DIFF=$(($TS - $MTS))
		if ( test $DIFF -lt 0 ) ; then
			DIFF=$(($DIFF * -1))
		fi
		if ( test $DIFF -lt $MIN ) ; then
			BESTTS=$MTS
			BESTTEMP=`echo $MLINE | awk '{print $2}'`
			MIN=$DIFF
		fi
		if ( test $LASTDIFF -eq 0 ) ; then
			LASTDIFF=$DIFF
		else
			if ( test $DIFF -gt $LASTDIFF ) ; then
				break
			fi
		fi
		MCNT=$(($MCNT+1))
	done
	echo $TS $TEMP $BESTTS $BESTTEMP
	
	CNT=$(($CNT+1))
done

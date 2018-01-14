#!/bin/bash

#HERE=/root/bin
HERE=`pwd`
BIN=/Users/rich/bin


# 2013-02-06 21:04:15

DSTR=$1
TSTR=$2

MONTH=`echo $DSTR | sed 's$/$ $g' | awk '{print $1}'`
DAY=`echo $DSTR | sed 's$/$ $g' | awk '{print $2}'`
YEAR=`echo $DSTR | sed 's$/$ $g' | awk '{print $3}'`
YEAR="20"$YEAR


HR=$(($TSTR / 100))
HMIN=$(($HR * 100)) 
MIN=$(($TSTR - $HMIN))

echo "$YEAR-$MONTH-$DAY $HR:$MIN:00"

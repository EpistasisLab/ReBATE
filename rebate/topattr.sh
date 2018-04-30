#!/bin/bash
#===============================================================================
#
#          FILE:  maketop.sh
# 
#         USAGE:  ./maketop.sh 
# 
#   DESCRIPTION:  Create subset for data based on top number of scores
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  input data file and scores results
#        AUTHOR:  Peter Robert Schmitt (@HOME), pschmitt@upenn.edu
#       COMPANY:  University of Pennsylvania
#       VERSION:  1.0
#       CREATED:  07/12/2016 08:58:37 EDT
#      REVISION:  ---
#===============================================================================

INPUT=$1
SCORES=$2
TOP=$3
CLASS=$4
H='header.txt'
A='attributes.txt'
T="top-${TOP}-$INPUT"
#touch $T
if test "$SCORES" = "" -o "$INPUT" = "" -o "$TOP" = "" -o "$CLASS" = ""
then
    echo "Syntax: $0 input-file scores-file top-value class-name"
    exit 1
fi
if test -f $SCORES -a -f $INPUT
then
    :
else
    echo "$SCORES and/or $INPUT do not exist"
    exit 2
fi

echo "creating $H"
head -1 $INPUT > $H

echo "creating $A"
a=1
for i in `cat $H`
do
    echo "$a:$i" >> $A
    echo "$a:$i"
    a=$((a + 1))
done

echo "creating columns"
while
    read line
do
    SNP=`echo $line | awk '{print $1}'`
    CNT=`echo $line | awk '{print $3}'`
    col=`grep ":$SNP$" $A | cut -f1 -d:`
    printf "snp:$NP\tscore:$CNT\tcolumn:$col\n"
    cut -f$col $INPUT > $SNP.dat
    if test $CNT = $TOP
    then
        break
    fi
done < $SCORES
echo $CLASS
col=`grep ":$CLASS$ $A | cut -f1 -d:`
cut -f$col $INPUT > $CLASS.dat

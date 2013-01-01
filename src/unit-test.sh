#!/bin/bash
n=$1
k=$2
file_name=$3
conf_file=conf-$n-$k-$file_name
#echo $k
#echo $file_name
#echo $conf_file
chunk_name=""
declare -i i=1
declare -i number=1
while [ $i -le $k ]
do
#declare -i number=$RANDOM*$k/32768
#echo $number

let "number = n-k-1+i"

chunk_name=_$number\_$file_name

echo $chunk_name

echo -e $chunk_name >> $conf_file
let "i += 1"
#echo $number
done
#echo $chunk_name

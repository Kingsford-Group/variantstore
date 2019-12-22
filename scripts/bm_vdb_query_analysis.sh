#!/bin/bash

# vcf file path:
vdb_file="$1"
len="$2"
query_len="$3"
mode="$4"
sample="$5"
echo "Query length: $query_len"

start=16050000  # for chr22
#start=10000   # for chr2
end=$(($start+$len-$query_len))
# Generate n random starting positions
pos=$(echo $(shuf -i $start-$end -n 1))


r1=''
r2=''

for i in $pos
do
  j=$(($i+$query_len))
  r1="${r1}${i}:${j},"
  r2="${r2}${i},"
done

r1=${r1::-1}
r2=${r2::-1}

#set -x
#echo "Benchmarking query type 1"
#./variantdb query -t 1 -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
#echo "Benchmarking query type 2"
#./variantdb query -t 2 -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
#echo "Benchmarking query type 4"
#./variantdb query -t 4 -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
#echo "Benchmarking query type 5"
#./variantdb query -t 5 -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
echo "Benchmarking query type 6"
./variantdb query -t 6 -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt

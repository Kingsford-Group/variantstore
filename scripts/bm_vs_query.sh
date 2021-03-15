#!/bin/bash

# vcf file path:
vs_file="$1"
num="$2"
query_len="$3"
len="$4"
type="$5"
mode="$6"
sample="$7"
echo "Query length: $query_len"

#start=16050000  # for chr22
start=10000   # for chr2
end=$(($start+$len-$query_len))
# Generate n random starting positions
pos=$(echo $(shuf -i $start-$end -n $num))


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

set -x
if (($type == 1))
then
  echo "Benchmarking query type $type"
  ./variantstore query -t $type -p $vs_file -r $r2 -m $mode -o temp.txt
else
  echo "Benchmarking query type $type"
  ./variantstore query -t $type -p $vs_file -r $r1 -m $mode -s $sample -o temp.txt
fi

#!/bin/bash
# Author: Yinjie Gao, yinjieg@andrew.cmu.edu
# THis is the query benchmark using bcftools
# equivalent to variantdb function -- get_var_in_ref

# vcf file path:
vdb_file="$1"
num="$2"
len="$4"
type="$5"
pct="$3"
query_len=$(($len/$pct))
mode="$6"
sample="$7"
echo "Query length: $query_len"

start=16050000
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
if (($type == 3))
then
  echo "Benchmarking query type $type"
  ./variantdb query -t $type -p $vdb_file -r $r2 -m $mode -o temp.txt
else
  echo "Benchmarking query type $type"
  ./variantdb query -t $type -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
fi

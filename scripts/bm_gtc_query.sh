#!/bin/bash

# vcf file path:
vdb_file="$1"
num="$2"
query_len="$3"
len="$4"
type="$5"
mode="$6"
sample="$7"
echo "Query length: $query_len"

start=16050000  # for chr22
#start=10000   # for chr2
end=$(($start+$len-$query_len))
# Generate n random starting positions
pos=$(echo $(shuf -i $start-$end -n $num))

for i in $pos
do
  j=$(($i+$query_len))
if (($type == 1))
then
  echo "Benchmarking query type $type"
  ./gtc view -r 2:$i-$j -o query.out $vdb_file
else
  echo "Benchmarking query type $type"
  ./gtc view -r 2:$i-$j -s $sample -o query.out $vdb_file
fi
done



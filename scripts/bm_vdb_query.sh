#!/bin/bash
# Author: Yinjie Gao, yinjieg@andrew.cmu.edu
# THis is the query benchmark using bcftools
# equivalent to variantdb function -- get_var_in_ref

# vcf file path:
vdb_file="$1"

num="$2"
pct="$3"
len="$4"
type="$5"

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

r3=''
a=''
b=''

if (($type == 7))
then
  vcf_file="$8"
  zcat $vcf_file | tail -n +254 |shuf -n $num > random_vcf.txt
  regions=$(echo $(less random_vcf.txt |cut -f 2))
  refs=$(echo $(less random_vcf.txt |cut -f 4))
  alts=$(echo $(less random_vcf.txt |cut -f 5))

  for i in $regions
  do
    r3="${r3}${i},"
  done
  r3=${r3::-1}

  for i in $refs
  do
    b="${b}${i},"
  done
  b=${b::-1}

  for i in $alts
  do
    a="${a}${i},"
  done
  a=${a::-1}
fi

echo $r3
echo $a
echo $b

rm random_vcf.txt


set -x
if (($type == 3))
then
  echo "Benchmarking query type $type"
  ./variantdb query -t $type -p $vdb_file -r $r2 -m $mode -o temp.txt
else
  if (($type == 7))
  then
    echo "Benchmarking query type $type"
    ./variantdb query -t $type -p $vdb_file -r $r3 -a $a -b $b -m $mode

  else
    echo "Benchmarking query type $type"
    ./variantdb query -t $type -p $vdb_file -r $r1 -m $mode -s $sample -o temp.txt
  fi
fi

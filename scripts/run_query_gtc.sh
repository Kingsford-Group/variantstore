set -x
for i in {1..2}
	do
		echo $i;
    for j in 10 100 1000
      		do
      #/usr/bin/time ./bm_gtc_query.sh gtc_v2/chr2/output $j 43185 243199373 $i 0 HG01872 >> query.txt;
      /usr/bin/time ./bm_gtc_query.sh gtc_v2/chr22/output $j 43185 51304566 $i 0 HG01872 >> query.txt;
		done
	done


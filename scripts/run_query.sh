:for i in {1..6}
	do
		echo $i;
		for j in 10 100 1000
      		do
			./bm_vdb_query.sh vdb_v4/chr22/ $j 43185 51304566 $i 0 HG01872 >> query.txt;
			#./bm_vdb_query.sh vdb_v4/chr2/ $j 43185 243199373 $i 0 HG01872 >> query.txt;
      #./bm_vdb_query.sh /dbgap/ckingsf_A018251/ppandey2/tcga-new/vdb_luad/chr22/ $j 42776 50818468 $i 0 184:TUMOR >> query_luad.out
      #./bm_vdb_query.sh /dbgap/ckingsf_A018251/ppandey2/tcga-new/vdb_luad/chr2/ $j 42776 242193529 $i 0 184:TUMOR >> query_luad.out
		done
	done


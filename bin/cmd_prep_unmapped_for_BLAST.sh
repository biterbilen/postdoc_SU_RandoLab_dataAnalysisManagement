
for i in Trl Trx; do
	for f in RM_TR_ada_vanVelthoven_2015_PullDown\$i/*txt.gz; do 
		echo $f; 
		less $f | head -n 100000 | sort | uniq -c | sort -k 1,1gr | \
			head -n 10 | awk '{ print ">"$1"\n"$2}'; 
	done  > $i\_unmapped_reads.fasta
done

less dta1.txt | grep -Pv "\tH[0-9][^0-9][0-9]" | cut -f 3-4 | sort | uniq | bedtools groupby -g 1 -c 2 -o count  | sort -k 2,2gr | head

for cl in K562 GM12878 HepG2 HeLa-S3 H1-hESC    myocyte C2C12 ES-E14; do
  less dta1.txt | grep -Pv "\tH[0-9][^0-9][0-9]" | grep -w $cl | cut -f 4 | sort | uniq > $cl.TF
done

cat myocyte.TF C2C12.TF ES-E14.TF | sort | uniq -c | sort -k 1,1gr | head 
echo
cat myocyte.TF C2C12.TF ES-E14.TF | sort | uniq -c | sort -k 1,1gr | awk '$1==3{ print $2}' 

#ln -sf ~/Projects/biter_biter_shared/bin/fastTrack4NoiseLM.R .
# Usage: sh cmds_4_NoiseLM2016.sh Brett_2014_LCMS1_QSC
# Usage: sh cmds_4_NoiseLM2016.sh Brett_2014_LCMS2_QSC
# ----
# less Ikeda_2014_RNAseq_cleaned_repCorr.txt | sed 's/_SC_/_SC:/g' | awk -F ':' '{ OFS="\t"; print $1,$2,$3 }' | awk '$1==$3{ OFS="\t"; print $1,$5}' > b
# Rscript -e 'b <- read.table("b"); t.test(b[b$V1=="X18m_NA_muscle_SC",2], y = b[b$V1=="X04m_NA_muscle_SC",2], "greater")'
# t = 3.4419, df = 42.668, p-value = 0.0006525
# Rscript -e 'b <- read.table("b"); t.test(b[b$V1=="X18m_21dE_muscle_SC",2], y = b[b$V1=="X04m_21dE_muscle_SC",2], "greater")'
# t = 3.2064, df = 26.008, p-value = 0.001773
# Rscript -e 'b <- read.table("b"); t.test(b[b$V1=="X04m_NA_muscle_SC",2], y = b[b$V1=="X04m_04dPI_Cmuscle_SC",2], "greater")'
# t = 1.6219, df = 6.8711, p-value = 0.07483
# Rscript -e 'b <- read.table("b"); t.test(b[b$V1=="X04m_NA_muscle_SC",2], y = b[b$V1=="X04m_21dE_muscle_SC",2], "greater")'
# ----
# less Brett_2015_RNAseq_cleaned_repCorr.txt | sed 's/_SC_18hGM/_SC_18hGM:/g' | awk -F ':' '{ OFS="\t"; print $1,$2,$3 }' | awk '$1==$3{ OFS="\t"; print $1,$5}' > bb2
# Rscript -e 'b <- read.table("bb2"); t.test(b[b$V1=="X22m_NA_muscle_SC_18hGM",2], y = b[b$V1=="X04m_NA_muscle_SC_18hGM",2], "greater")'
# t = 1.9032, df = 2.6787, p-value = 0.08204
# less Brett_2015_RNAseq_cleaned_repCorr.txt | sed 's/_SC_NA/_SC_NA:/g' | awk -F ':' '{ OFS="\t"; print $1,$2,$3 }' | awk '$1==$3{ OFS="\t"; print $1,$5}' > bb
# Rscript -e 'b <- read.table("bb"); t.test(b[b$V1=="X22m_03dPI_muscleTG_SC_NA",2], y = b[b$V1=="X04m_03dPI_muscleTG_SC_NA",2], "greater")'
# t = 2.1851, df = 5.4297, p-value = 0.03815
# Rscript -e 'b <- read.table("bb"); t.test(b[b$V1=="X22m_NA_muscle_SC_NA",2], y = b[b$V1=="X04m_28dPI_muscleTG_SC_NA",2], "greater")'
# t = 1.6134, df = 9.9788, p-value = 0.0689

# ----
# ----
# less Brett_2014_LCMS1_QSC_all_repCorr.txt | sed 's/_SC_/_SC:/g' | awk -F ':' '{ OFS="\t"; print $1,$2,$3 }' | awk '$1==$3{ OFS="\t"; print $1,$5}' > a
# Rscript -e 'b <- read.table("a"); t.test(b[b$V1=="X22m_NA_SC",2], y = b[b$V1=="X04m_NA_SC",2], "greater")'
# t = 1.4664, df = 38.303, p-value = 0.07535
# ----
# less Brett_2014_LCMS2_QSC_all_repCorr.txt | sed 's/_SC_/_SC:/g' | awk -F ':' '{ OFS="\t"; print $1,$2,$3 }' | awk '$1==$3{ OFS="\t"; print $1,$5}' > aa
# Rscript -e 'b <- read.table("aa"); t.test(b[b$V1=="X22m_NA_SC",2], y = b[b$V1=="X04m_NA_SC",2], "less")'
# t = -1.72, df = 129.79, p-value = 0.04391
# ----

tag=$1

less ~/Projects/biter_jbrett_metabolomics/results/2015-09-19/XC_Brett_2014_$tag/peakTable.txt.gz | \
	sed 's/_//g' | sed 's/OQ/22m_NA_SC_/g' | sed 's/YQ/04m_NA_SC_/g' | \
	awk 'BEGIN{OFS="\t"}{ if (NR==1) { print "Id",$0; } else if (NR>2) { $1=NR; print $0; } }' | \
	gzip -c > XC_${tag}_peakTable.txt.gz

tf=XC_${tag}_peakTable.txt.gz
tag2=${tag}_all
Rscript -e "source('auxggplot2.R'); test.get.exploratoryAnalysisPlots('$tag2','$tf')" &> log.explore_$tag2



exit
#---------------------------

#ln -sf ~/Projects/biter_biter_shared/bin/auxggplot2.R .
#ln -sf ~/Projects/biter_biter_shared/bin/fastTrack4NoiseLM.R .
# Usage: sh cmds_4_NoiseLM2016.sh Liu_2014_RNAseq
# Usage: sh cmds_4_NoiseLM2016.sh Ikeda_2014_RNAseq 4
# Usage: sh cmds_4_NoiseLM2016.sh Brett_2015_RNAseq 

tag=$1
indices2filters=${2:-F}
echo $tag
#for tag in Liu_2014_RNAseq Ikeda_2014_RNAseq Brett_2015_RNAseq; do
#for tag in Liu_2014_RNAseq; do
# summary
ln -sf ../QA_$tag/baseCalls.txt QA_$tag.summary
ln -sf ../FC_$tag/summary FC_$tag.summary
ln -sf ../EG_$tag/summary EG_$tag.summary
Rscript -e "source('fastTrack4NoiseLM.R'); test.combineFiles('$tag')" &> log.combine_$tag

# explore
sf=../EG_$tag/featureCounts_GMM2expressed.txt.gz
tf=EG_${tag}_featureCounts_GMM2expressed.txt.gz
ln -sf $sf $tf 

if [ $indices2filters == F ]; then
	samples2filter=($(less FC_$tag.summary | sed 's/RNAseq_/\t/g' | awk '$3<5000000{ print $2}'))
	samples2filters=$(echo ${samples2filter[@]} | sed 's/ /\|/g')
	indices2filter=($(less EG_${tag}_featureCounts_GMM2expressed.txt.gz | head -n 1 | cut -f 2- | sed 's/\t/\n/g' | awk '{ print NR"\t"$1}' | grep -P "$samples2filters" | cut -f 1)) 
	indices2filters="$(echo ${indices2filter[@]} | sed 's/ /,/g')"
	echo ${samples2filter[@]}
	#echo $samples2filters
	#echo ${indices2filter[@]}
fi
echo Filtering $indices2filters

tag1=${tag}_cleaned
Rscript -e "source('auxggplot2.R'); test.get.exploratoryAnalysisPlots('$tag1','$tf',c($indices2filters))" &> log.explore_$tag1

tag2=${tag}_all
Rscript -e "source('auxggplot2.R'); test.get.exploratoryAnalysisPlots('$tag2','$tf')" &> log.explore_$tag2
#done

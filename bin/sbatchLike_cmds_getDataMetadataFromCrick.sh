
# Usage ./sbatchLike_cmds_getDataMetadataFromCrick.sh sdir[2015/aug]

command="rsync -rqt"

sdir='2014/may'
sdir=$1

# TODO set this
#dir=/srv/gs1/projects/scg/Archive/SolexaRuns/$sdir; account="biter@franklin.stanford.edu" # from archive for 2011-2013
dir=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/$sdir; account="biter@crick.stanford.edu" # for 2013-current

if [[ $sdir =~ 2011 ]]; then
	less metadata | grep -v Alias | sed -E 's/(.*)_L._/\1\t/' | cut -f 1 | sort | uniq | \
		while read i; do $command $account:$dir/*/$i*.csv .; done
else
	less metadata | grep -v Alias | sed -E 's/(_L.)_/\1\t/' | cut -f 1 | sort | uniq | \
		while read i; do $command $account:$dir/*/$i*.csv .; done
fi

if [[ $sdir =~ 2011 ]]; then
	less metadata | grep -v Alias | cut -f 1 | \
		while read i; do $command $account:$dir/*/$i.fastq.gz rawData/.; done
else
	less metadata | grep -v Alias | cut -f 1 | \
		while read i; do $command $account:$dir/*/*/$i.fastq.gz rawData/.; done
fi

echo
grep Submitter *csv



exit

# Liu_2014
dir=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2014/aug
$command $account:$dir/140801_TENNISON_0307_AC4D2TACXX/140801_TENNISON_0307_AC4D2TACXX_L6*csv .
$command $account:$dir/140801_TENNISON_0307_AC4D2TACXX/140801_TENNISON_0307_AC4D2TACXX_L7*csv .
$command $account:$dir/140805_BRISCOE_0171_BC4CR5ACXX/140805_BRISCOE_0171_BC4CR5ACXX_L7*csv .
$command $account:$dir/140805_BRISCOE_0171_BC4CR5ACXX/140805_BRISCOE_0171_BC4CR5ACXX_L8*csv .

# Cheung_2013
#130308_MONK_0271_BC1WGLACXX_L1_results.html	Young_QSC1
#130308_MONK_0271_BC1WGLACXX_L2_results.html	Young_QSC2
#130320_TENNISON_0218_AD1RYHACXX_L7_results.html	MuSC_ASC1
#130320_TENNISON_0218_AD1RYHACXX_L8_results.html	MuSC_ASC2
dir=/srv/gs1/projects/scg/Archive/SolexaRuns/2013/mar # from archive
$command $account:$dir/130308_MONK_0271_BC1WGLACXX/L1/*csv .
$command $account:$dir/130308_MONK_0271_BC1WGLACXX/L2/*csv .
$command $account:$dir/130320_TENNISON_0218_AD1RYHACXX/L7/*csv .
$command $account:$dir/130320_TENNISON_0218_AD1RYHACXX/L8/*csv .
$command $account:$dir/130308_MONK_0271_BC1WGLACXX/L1/130308_MONK_0271_BC1WGLACXX_L1_?_pf.fastq.gz .
$command $account:$dir/130308_MONK_0271_BC1WGLACXX/L2/130308_MONK_0271_BC1WGLACXX_L2_?_pf.fastq.gz .
$command $account:$dir/130320_TENNISON_0218_AD1RYHACXX/L7/130320_TENNISON_0218_AD1RYHACXX_L7_?_pf.fastq.gz .
$command $account:$dir/130320_TENNISON_0218_AD1RYHACXX/L8/130320_TENNISON_0218_AD1RYHACXX_L8_?_pf.fastq.gz .

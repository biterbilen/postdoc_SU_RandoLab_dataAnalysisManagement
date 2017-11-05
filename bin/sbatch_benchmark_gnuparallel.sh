#!/bin/bash -l 
# sbatch command: sbatch -J PF_HrWTK9me2_Liu_2016_ChIPseq -p normal -t 6:0:0 -c 8 			-e slurm_PF_HrWTK9me2_Liu_2016_ChIPseq.out -o slurm_PF_HrWTK9me2_Liu_2016_ChIPseq.out bin/sbatch_PF_HrWTK9me2_Liu_2016_ChIPseq.sh
samples=(HrWT_young_TMX_SC_NA_H3K9me2_1.sorted.bam HrWT_young_TMX_SC_NA_H3K9me2_2.sorted.bam)
controls=(HrWT_young_TMX_SC_NA_NA_1.sorted.bam)
cat > a.txt <<- END
n=\$(basename ${samples[0]} .sorted.bam); echo \$n; macs2 callpeak -t ${samples[0]} -c ${controls[@]} -g mm -n \$n &> slurm_\$n 
n=\$(basename ${samples[1]} .sorted.bam); echo \$n; macs2 callpeak -t ${samples[1]} -c ${controls[@]} -g mm -n \$n &> slurm_\$n 
END

# TODO rm echo to run the commands
parallel -k --delay 0 --plus --tmpdir . -j0 'echo {} {.} {/}' :::: a.txt # keep order
parallel --delay 0 --plus --tmpdir . -j0 'echo {} {.} {/}' :::: a.txt
parallel --delay 0 --plus --tmpdir . --xapply echo ::: A B C ::: D E F # one argument from each source
parallel --delay 0 --plus --tmpdir . -j0 'echo {} {.} {/} {/.} {..} {/..} &> {/..}.log' ::: A/a.txt.gz
parallel --delay 0 --plus --tmpdir . -j0 --rpl '{gpdr} s:/*[^/]*/::' 'echo {/} {gpdr}' ::: /C/B/A/a.txt.gz #custom replacement string for the grandparent directory removed DIR
parallel --delay 0 --plus --tmpdir . -j0 'echo --age {1} --sex {2} --chr {3}' ::: {1..3} ::: M F ::: {1..22} X Y
parallel --delay 0 --plus --tmpdir . --results outputdir 'echo --age {1} --sex {2}' ::: {1..3} ::: M F 
parallel --delay 0 --plus --tmpdir . --results outputdir --header : 'echo --age {AGE} --sex {SEX}' ::: AGE {1..3} ::: SEX M F 
parallel --delay 0 --plus --tmpdir . --header : 'echo --sex {SEX.} --chr {CHR}' ::: SEX M.txt F.txt ::: CHR {1..22} X Y
parallel --delay 0 --plus --tmpdir . --header : 'echo --sex {SEX..} --chr {CHR}' ::: SEX M.txt.gz F.txt.gz ::: CHR {1..22} X Y # Does NOT work
parallel --delay 0 --plus --tmpdir . --header : 'sleep {};echo Jobslot {JS%} slept {} seconds' ::: JS 4 3 2 1
# No workaround to mix stderr and stdout

# TODO check --shebang-wrap option to parallelize own scripts 
# https://www.biostars.org/p/63816/
inlf() {
	local odir=$3
	mkdir -p $3; pushd $3
	echo $1 > $2
	popd
}
export -f inlf
# FTAG TAG might be reserved word
r=OUT; rm -rf $r; mkdir -p $r
parallel --delay 1 --plus --tmpdir . -j0 --header : --tag --results $r "inlf {TAG} {FN} $r/TAG/{TAG}/FN/{FN}" ::: TAG content1 content2 ::: FN file1 file2
r=OUT2; rm -rf $r; mkdir -p $r
parallel --delay 1 --plus --tmpdir . -j0 --header : --tag --results $r "inlf {FTAG} {FN} $r/CNT/{CNT}/FN/{FN}" ::: CNT content1 content2 ::: FN file1 file2
tree OUT OUT2

# Regular expressions and globbing
# pattern matches suffix prefix
array=( "${array[@]/%/_content}" )
#will append the 'content' string to each element.

array=( "${array[@]/#/prefix_}" )
#will prepend 'prefix_' string to each element

# http://www.softpanorama.org/Scripting/Shellorama/Reference/string_operations_in_shell.shtml
# replace all occurences of ,; mark double //
a=05hdPI_SC_ATR_1,NA_SC_ORC1_1,05hdPI_SC_ORC1_1vs05hdPI_SC_ATR_2,NA_SC_ORC1_2,05hdPI_SC_ORC1_2
echo ${a//,/ }


#Given:
var=/tmp/my.dir/filename.tar.gz
#We can use these expressions:
path=${var%/*}
#To get: /tmp/my.dir (like dirname)
file=${var##*/}
#To get: filename.tar.gz (like basename)
base=${file%%.*}
#To get: filename
ext=${file#*.}
#To get: tar.gz

# Additionally - parameter expansion
#https://www.gnu.org/software/bash/manual/bashref.html#Shell-Parameter-Expansion


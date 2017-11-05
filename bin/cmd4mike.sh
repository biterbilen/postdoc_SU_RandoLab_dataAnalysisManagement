  828  vi metadata 
  829  column -t metadata  | less
  830  rm -rf QA_Wold_2011_ChIPseq
  831  less C2C12_Fosl1.fastq.gz
  832  d ..
  833  ln -sf /home/biter/PI_SCRATCH/Projects/biter_biter_shared/bin/sbatch_tr.sh sbatch_tr.sh
  834  vi QA_Wold_2011_ChIPseq/
  835  cd QA_Wold_2011_ChIPseq/
  836  vi C2C12_Fosl1.fastq.gz 
  837  less C2C12_Fosl1.fastq.gz 
  838  ll QA_Wold_2011_ChIPseq/
  839  shl
  840  im
  841  pw
  842  vi /home/biter/.bashrc 
  843  vi /home/biter/.bash_aliases 
  844  head TR_Wold_2011_ChIPseq/summary 
  845  head RM_Wold_2011_ChIPseq/summary 
  846  ln -sf /home/biter/PI_SCRATCH/Projects/biter_biter_shared/bin/sbatch_fc.sh 
  847  less FC_Wold_2011_ChIPseq/slurm_FC_15938036_0.out 
  848  less featureCounts.txt.gz 
  849  vi QA_Wold_2011_ChIPseq/metadata 
  850  vi sbatch_qa.sh 
  851  rm -rf QA_Wold_2011_ChIPseq TR_Wold_2011_ChIPseq RM_Wold_2011_ChIPseq FC_Wold_2011_ChIPseq
  852  sh sbatch_qa.sh 
  853  sh sbatch_tr.sh 
  854  vi sbatch_rm.sh sqbit
  855  sh sbatch_tr.sh T
  856  sh sbatch_rm.sh 
  857  samtools view -F 0x40 C2C12_Fosl1_1.sorted.bam | less
  858  samtools view -f 0x40 C2C12_Fosl1_1.sorted.bam | less
  859  samtools view -f 0x40 C2C12_Fosl1_1.sorted.bam 
  860  java --version
  861  java -version
  862  wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar
  863  java -jar /path/to/picard.jar -h
  864  java -jar picard.jar -h
  865  vi /home/biter/Projects/biter_biter_management/bin/bb_apps.sh
  866  java -jar picard.jar 
  867  java -jar picard.jar  -h
  868  java -jar picard.jar MarkDuplicates C2C12_Fosl1_1.sorted.bam
  869  java -jar picard.jar MarkDuplicates -h
  870  java -jar picard.jar MarkDuplicates -I=C2C12_Fosl1_1.sorted.bam -O=wdup_C2C12_Fosl1_1.sorted.bam -M=log
  871  \
  872  java -jar picard.jar MarkDuplicates I=C2C12_Fosl1_1.sorted.bam O=wdup_C2C12_Fosl1_1.sorted.bam M=log
  873  sdev -N 2
  874  less sbatch_fc.sh 
  875  less FC_Wold_2011_ChIPseq/slurm_FC_15945110_0.out 
  876  cd 2017-06-26
  877  rm picard.jar 
  878  sh sbatch_rm.sh T
  879  chip
  880  less FC_Wold_2011_ChIPseq/bin/FCFC_Wold_2011_ChIPseq.sh 
  881  less FC_Wold_2011_ChIPseq/
  882  vi FC_Wold_2011_ChIPseq/bin/FCFC_Wold_2011_ChIPseq.sh 
  883  less FC_Wold_2011_ChIPseq/slurm_FC_15958076_0.out 
  884  sh sbatch_fc.sh T
  885  less FC_Wold_2011_ChIPseq/summary 
  886  ll FC_Wold_2011_ChIPseq/
  887  les slurm_FC_15958076_0.out 
  888  less slurm_FC_15958076_0.out 
  889  cd FC_Wold_2011_ChIPseq/
  890  vi sbatch_fc.sh 
  891  rm -rf FC_Wold_2011_ChIPseq
  892  ln -sf /home/biter/PI_SCRATCH/Projects/biter_biter_shared/bin/sbatch_pf.sh .
  893  vi sbatch_pf.sh 
  894  less sbatch_rm.sh 
  895  vi sbatch_tr.sh 
  896  less sbatch_qa.sh 
  897  less sbatch_tr.sh 
  898  ll PF_Wold_2011_ChIPseq/
  899  less slurm_PF_15963353_0.out 
  900  scancel 15962989
  901  rm -rf FC_Wold_2011_ChIPseq PF_Wold_2011_ChIPseq
  902  scancel 15963353
  903  sh sbatch_fc.sh 
  904  less slurm_PF_15963485_0.out 
  905  scancel 15963485
  906  less slurm_PF_15964144_0.out 
  907  less slurm_PF_15964144_0.out  | grep "predicted fragment length "
  908  less slurm_PF_15964144_0.out  | grep "predicted fragment length " | awk '{ print ${NF-1)}'
  909  less slurm_PF_15964144_0.out  | grep "predicted fragment length " | awk '{ print $(NF-1)}'
  910  less slurm_PF_15964144_0.out  | grep "predicted fragment length " | awk '{ print $(NF-1) }'
  911  grep "predicted fragment length " slurm_PF_15964144_0.out | awk '{ print $(NF-1) }'
  912  scancel 15964144
  913  sqit
  914  less slurm_PF_15964144_0.out
  915  l
  916  less slurm_PF_15965204_0.out 
  917  less slurm_PF_15965204_1.out 
  918  les bin/PFPF_Wold_2011_ChIPseq.sh 
  919  less bin/PFPF_Wold_2011_ChIPseq.sh 
  920  grep "predicted fragment length" slurm_PF_*_0.out | awk '{ print $(NF-1) }'
  921  cat bin/PFPF_Wold_2011_ChIPseq.sh 
  922  extsize=$(grep "predicted fragment length" slurm_PF_*_${o[stiv]}.out | awk '{ print $(NF-1) }');
  923  Rscript C2C12_Fosl1_1_predictd.R
  924  less C2C12_Fosl1_1_predictd.R
  925  less C2C12_Fosl1_1_predictd.R | grep altd
  926  less C2C12_Fosl1_1_predictd.R | grep ^altd | sed 's/(//g'
  927  less C2C12_Fosl1_1_predictd.R | grep ^altd | sed 's/[c()]//g'
  928  less C2C12_Fosl1_1_predictd.R | grep ^altd | sed 's/[c()]//g' | awk '{ print $2}'
  929  less C2C12_Fosl1_1_predictd.R | grep ^altd | sed 's/[c()]//g' | awk '{ print $3}'
  930  cd PF_Wold_2011_ChIPseq/
  931  less slurm_PF_15965213_0.out 
  932  less slurm_PF_15965213_1.out 
  933  less C2C12_Fosl1_1_predictd.R 
  934  Rscript C2C12_Fosl1_1_predictd.R 
  935  Rscript C2C12_NA_1_predictd.R 
  936  cd ..
  937  scancel 15965213
  938  rm -rf PF_Wold_2011_ChIPseq
  939  sh sbatch_pf.sh 
  940  sqbit
  941  pwd
  942  bin
  943  vi sbatch_rm.sh 
  944  less auxML.R | grep edgeR
  945  vi auxggplot2.R
  946  less RNAseq_DE.txt.gz 
  947  less PullDown_topKEGG.txt.gz 
  948  vi auxML.R
  949  join -1 <(less PullDown_DE.txt.gz | sort -k 8,8) -2 <(less RNAseq_DE.txt.gz | sort -k 8,8) -i 8 -j 8
  950  join <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8
  951  join <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  952  man join
  953  join -t "\t" <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  954  join -t$"\t" <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  955  join -t`$\t` <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  956  join -t`\t` <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  957  join -t$`\t` <(less PullDown_DE.txt.gz | sort -k 8,8) <(less RNAseq_DE.txt.gz | sort -k 8,8) -1 8 -2 8 | less
  958  join -t$`\t` <(less PullDown_DE.txt.gz | sort -k 7,7) <(less RNAseq_DE.txt.gz | sort -k 7,7) -1 7 -2 7 | less
  959  less PullDown_DE.txt.gz | sort -k 7,7
  960  less PullDown_DE.txt.gz | sort -k 7,7 | head
  961  less PullDown_DE.txt.gz | head
  962  less PullDown_DE.txt.gz | grep -v "^#" | sort -k 7,7 > pulldown
  963  less RNAseq_DE.txt.gz | grep -v "^#" | sort -k 7,7 > rnaseq
  964  less RNAseq_DE.txt.gz | grep -v "^#" | sort -k 1,1 > rnaseq
  965  less PullDown_DE.txt.gz | grep -v "^#" | sort -k 1,1 > pulldown
  966  join pulldown rnaseq | less
  967  join pulldown rnaseq -t $`\t` | less
  968  join pulldown rnaseq -t$`\t` | less
  969  join pulldown rnaseq -t$'\t' | less
  970  join pulldown rnaseq -t$'\t' | column -t | less
  971  join pulldown rnaseq -t$'\t' | column -t | head -n 1
  972  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $14 <0.01{ print $8 }' | less
  973  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $14 < 0.1{ print $8 }' | less
  974  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $13 < 0.1{ print $8 }' | less
  975  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $13 < 0.1{ print $1,$7 }' | less
  976  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $13 < 0.1{ OFS="\t"; print $1,$7 }' | less
  977  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $13 < 0.1{ OFS="\t"; print $1,$7 }' | wc -l
  978  join pulldown rnaseq -t$'\t' | awk '$6<0.1 && $13 < 0.1{ OFS="\t"; print $1,$7 }' | head
  979  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $1,$7 }' | wc -l
  980  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $1,$7 }' 
  981  echo Cend1 Cntfr Crispld2 Cyp1b1 Draxin Dtnb Gfra1 Has2 Htra4 Itga11 Kcnip3 Kcnj15 P4ha3 Pet100 Podn Prg4 Rnf150 Runx1 Slc41a1 Sod3 Susd2 Tmem45a Tubb4a
  982  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $1,$7 }' | sort -k 2,2
  983  echo `join pulldown rnaseq -t$'\t' | awk '$6<0.05 && $13 < 0.05{ OFS="\t"; print $7 }' | sort -k 2,2`
  984  echo `join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2`
  985  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2
  986  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 1{ OFS="\t"; print $7 }' | sort -k 2,2 | wc -l
  987  join pulldown rnaseq -t$'\t' | awk '$6<1 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2 | wc -l
  988  join pulldown rnaseq -t$'\t' | awk '$6<1 && $13 < 1{ OFS="\t"; print $7 }' | sort -k 2,2 | wc -l
  989  find ~/Projects/biter_biter_EWSR1/bin/ -name "fisher"
  990  find ~/Projects/biter_biter_EWSR1/bin/ -name "*fisher*"
  991  R --no-save /home/biter/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R 24 298 994 11966 1
  992  R --no-save --args 24 298 994 11966 1 < /home/biter/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R
  993  R --no-save 24 298 994 11966 1  /home/biter/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R
  994  Rscript /home/biter/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R 24 298 994 11966 1
  995  Rscript $HOME/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R 24 298 994 11966 1
  996  less /home/biter/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R
  997  join pulldown /scratch/PI/casco/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FF_final/tfsmirs.txt | less
  998  join pulldown <(less /scratch/PI/casco/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FF_final/tfsmirs.txt | sort -k 1,1) | less
  999  join -t $'\t' pulldown <(less /scratch/PI/casco/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FF_final/tfsmirs.txt | sort -k 1,1) ) 
 1000  join -t $'\t' pulldown <(less /scratch/PI/casco/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FF_final/tfsmirs.txt | sort -k 1,1) | less
 1001  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2 | wc -l
 1002  join pulldown rnaseq -t$'\t' | awk '$6<0.05 && $13 < 0.05{ OFS="\t"; print $7 }' | sort -k 2,2 
 1003  join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2 
 1004  echo `join pulldown rnaseq -t$'\t' | awk '$6<1 && $13 < 0.01{ OFS="\t"; print $7 }' | sort -k 2,2 `
 1005  echo `join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13 < 1{ OFS="\t"; print $7 }' | sort -k 2,2 `
 1006  less pulldown 
 1007  less pulldown | head -n 1
 1008  echo `join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $2>0 && $13 < 1{ OFS="\t"; print $7 }' | sort -k 2,2 `
 1009  echo `join pulldown rnaseq -t$'\t' | awk '$6<1 && $13<0.01 && $9>0{ OFS="\t"; print $7 }' | sort -k 2,2 `
 1010  echo `join pulldown rnaseq -t$'\t' | awk '$6<0.01 && $13<1 { OFS="\t"; print $7 }' | sort -k 2,2 `
 1011  echo `join pulldown rnaseq -t$'\t' | awk '$6<1 && $13<0.01 { OFS="\t"; print $7 }' | sort -k 2,2 `
 1012  R --no-save --args 24 298 994 11966 1 < $HOME/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R
 1013  R --no-save --args 24 298 994 11966 1 fisher.ex.t.R 
 1014  R --no-save --args 24 298 994 11966 1 < fisher.ex.t.R 
 1015  R CMD BATCH --no-save --args '24 298 994 11966 1' fisher.ex.t.R 
 1016  R CMD BATCH --no-save --args '24 298 994 11966 1' fisher.ex.t.R a.out
 1017  cat a.out
 1018  rm a.out
 1019  cp $HOME/Projects/biter_biter_EWSR1/bin/scripts/_DIS3/fisher.ex.t.R .
 1020  less fisher.ex.t.R 
 1021  tmux
 1022  ll
 1023  rm fisher.ex.t.R 
 1024  ln -sf ~/PI_HOME/Applications/bilebi00/Project_EWSR1/scripts/_EWSR1/fisher.ex.t.R 
 1025  R --no-save --args 24 298 994 11966 <  fisher.ex.t.R
 1026  R --no-save --args 24 298 994 11966 1 <  fisher.ex.t.R
 1027  history | tail -n 200 > cmd4mike.sh 

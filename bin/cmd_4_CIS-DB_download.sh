#From Bulk Downloads - Logos, E-scores, Z-scores, Probe Intensities, TF info, PWMs; Download Species Archive
#http://cisbp.ccbr.utoronto.ca/bulk.php
# human
wget http://cisbp.ccbr.utoronto.ca/tmp/Homo_sapiens_2016_06_28_11:22_pm.zip
unzip Homo_sapiens_2016_06_28_11\:22_pm.zip
rm Homo_sapiens_2016_06_28_11\:22_pm.zip
for f in pwms_all_motifs/*; do echo ">"$(basename $f .txt); cat $f | grep -v Pos | cut -f 2-; done > all_pwms.txt

exit
# mouse
wget http://cisbp.ccbr.utoronto.ca/tmp/Mus_musculus_2016_01_28_7:52_pm.zip
unzip Mus_musculus_2016_01_28_7\:52_pm.zip 
rm Mus_musculus_2016_01_28_7\:52_pm.zip 
for f in pwms_all_motifs/*; do echo ">"$(basename $f .txt); cat $f | grep -v Pos | cut -f 2-; done > all_pwms.txt
exit

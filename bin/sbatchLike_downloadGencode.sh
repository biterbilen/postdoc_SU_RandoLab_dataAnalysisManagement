#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatchLike_downloadGencode.sh
# 
#         USAGE: ./sbatchLike_downloadGencode.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 08/16/2016 21:24
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

wget -N ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz


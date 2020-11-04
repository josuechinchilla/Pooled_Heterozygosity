#!/bin/bash
#==============================================================================#
# Header
#==============================================================================#
# ISU - Animal Breeding and Genetics Group
# Author: Josue Chinchilla-Vargas
# Created (date): 11/3/2020
# Version (of the script): 1.0
# Program (if applicable): Plink
# Program Version (if applicable): 1.9
#==============================================================================#
# Description
#==============================================================================#
#Calculates pooled heterozygosity (Hp) for 1 Mb windows

#==============================================================================#
#Keywords (Please put between 3 and 5)
#==============================================================================#
#signatures of selection
#poooled heterozygosity




#==============================================================================#
# Setup
#==============================================================================#
###Input files###

#  .frq.counts (from Plink1.9 --freq counts)
#Format: CHR SNP A1 A2 C1 C2 G0


#   .map used to calculate allele counts
#Format: CHR SNP DISTANCE BP

#Script works as follows:
#./Hp.sh 	myfile.frq.counts 	myfile.map	output_name

###Output files###
#Hp.txt
#Format:  CHR POS CHR_POS COUNT_MIN COUNT_MAJ NUM NUM Hp

#==============================================================================#

frq=$1
map=$2 
out=$3


set -o errexit -o nounset -o xtrace -o pipefail


#get the lines we need from .frq and map files.
awk -F ' ' '{print $1" "$2" "$5" "$6}' $frq | tail -n +2  > cutout.jo 
awk -F ' ' '{print $4}' $map > pos.jo

# calculate position in megabases and remove decimals. 
#This produces scientific notation, don't know how to avoid it so I removed the first 2 lines of map and frq file.
awk '{$1/=1000000}1' pos.jo | awk -F[.] '{print $1}' > pos_mb.jo 

#join and organize all info we need.
paste cutout.jo pos_mb.jo  >step1.jo
awk -F ' ' '{print $1"_"$5" "$1" "$5" "$2" "$3" "$4}' step1.jo > step2.jo
awk -F ' ' '{print $1" "$5}' step1.jo | awk '!seen[$0]++' > chrpos.jo

#generate files to save results.
touch min_sum.jo > min_sum.jo
touch maj_sum.jo > maj_sum.jo

#generate windows to use in loops.
awk -F ' ' '{print $1}' step2.jo | awk '!seen[$0]++' > windows.jo

#loops for minor and major allele counts.
for f in `awk -F ' ' '{print $1}' step2.jo | awk '!seen[$0]++'`

do

awk -F ' ' -v f="$f" '{if ($1 == f) print $0;}' step2.jo > min.jo
awk '{s+=$5}END{print s}' min.jo > results_min.jo
cat results_min.jo >> min_sum.jo #add new line at the bottom of the file to keep order
done

for t in `awk -F ' ' '{print $1}' step2.jo | awk '!seen[$0]++'`

do
awk -F ' ' -v t="$t" '{if ($1 == t) print $0;}' step2.jo > maj.jo
awk '{s+=$6}END{print s}' maj.jo > results_maj.jo
cat results_maj.jo >> maj_sum.jo #add new line at the bottom of the file to keep order

done

#generate "backbone" for final output.
paste chrpos.jo windows.jo min_sum.jo maj_sum.jo >noheader.jo

#calculate numerator and denominator, add them to "backbone".
awk '{$6 = 2 * ($4 + $5)}1' noheader.jo| awk '{$7 = ($4 + $5) ** 2}1' > almostHp.jo

#calculate Hp for each window.
awk '{$8 = $6 / $7}1' almostHp.jo > Hp_noheader.jo

#add header and paste everything together.
echo "CHR POS CHR_POS COUNT_MIN COUNT_MAJ NUM NUM Hp" > Hp_header.jo 
cat Hp_header.jo Hp_noheader.jo > $out.txt

#get rid of intermediate files.
rm *.jo


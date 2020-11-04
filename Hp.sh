# Header
#==============================================================================#
# ISU - Animal Breeding and Genetics Group
# Author: Josue Chinchilla-Vargas
# Created (date): 11/3/2020
# Version (of the script): 2.0
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

#####################_Hp calculations for windows starting at .5 megabases_#####################

#get the lines we need from .frq and map files.
awk -F ' ' '{print $1" "$2" "$5" "$6}' $frq | tail -n +2  > cutout_05mb.jo 

awk -F ' ' '{print $4}' $map > pos_05mb.jo

# calculate position in megabases and remove decimals. 
#This produces scientific notation, don't know how to avoid it so I removed the first 2 lines of map and frq file.
awk '{$1/=1000000}1' pos_05mb.jo > pos_halfmb.jo 
awk '{ printf("%.1f\n", $1) }' pos_halfmb.jo > almost_pos_05mb.jo
awk '{$2 = $1 - 0.5}1' almost_pos_05mb.jo > final_pos_05mb.jo
awk -F ' ' '{print $2}' final_pos_05mb.jo > add_to.jo
#join and organize all info we need.
paste cutout_05mb.jo add_to.jo | awk '$5>=0' > almost1_05mb.jo
awk -F[.] '{print $1}' almost1_05mb.jo > almost2._05mb.jo
awk -F ' ' '{print $1" "$2" "$3" "$4" "$5".5"}' almost2._05mb.jo > step1_05mb.jo


awk -F ' ' '{print $1"_"$5" "$1" "$5" "$2" "$3" "$4}' step1_05mb.jo >step2_05mb.jo

awk -F ' ' '{print $2" "$3}' step2_05mb.jo | awk '!seen[$0]++' > chrpos_05mb.jo

#generate files to save results.
touch min_sum_05mb.jo > min_sum_05mb.jo
touch maj_sum_05mb.jo > maj_sum_05mb.jo

#generate windows to use in loops.
awk -F ' ' '{print $1}' step2_05mb.jo | awk '!seen[$0]++' > windows_05mb.jo

#loops for minor and major allele counts.
for f in `awk -F ' ' '{print $1}' step2_05mb.jo | awk '!seen[$0]++'`

do

awk -F ' ' -v f="$f" '{if ($1 == f) print $0;}' step2_05mb.jo > min_05mb.jo
awk '{s+=$5}END{print s}' min_05mb.jo > results_min_05mb.jo
cat results_min_05mb.jo >> min_sum_05mb.jo #add new line at the bottom of the file to keep order
done

for t in `awk -F ' ' '{print $1}' step2_05mb.jo | awk '!seen[$0]++'`

do
awk -F ' ' -v t="$t" '{if ($1 == t) print $0;}' step2_05mb.jo > maj_05mb.jo
awk '{s+=$6}END{print s}' maj_05mb.jo > results_maj_05mb.jo
cat results_maj_05mb.jo >> maj_sum_05mb.jo #add new line at the bottom of the file to keep order

done

#generate "backbone" for final output.
paste chrpos_05mb.jo windows_05mb.jo min_sum_05mb.jo maj_sum_05mb.jo >noheader_05mb.jo

#calculate numerator and denominator, add them to "backbone".
awk '{$6 = 2 * ($4 + $5)}1' noheader_05mb.jo| awk '{$7 = ($4 + $5) ** 2}1' > almostHp_05mb.jo

#calculate Hp for each window.
awk '{$8 = $6 / $7}1' almostHp_05mb.jo > Hp_noheader_05mb.jo

#####################_Hp calculations for windows starting at each megabase_#####################

#get the lines we need from .frq and map files.
awk -F ' ' '{print $1" "$2" "$5" "$6}' $frq | tail -n +2  > cutout.jo 
awk -F ' ' '{print $4}' $map > pos.jo

# calculate position in megabases and remove decimals. 
#This produces scientific notation, don't know how to avoid it so I removed the first 2 lines of map and frq file.
awk '{$1/=1000000}1' pos.jo | awk -F[.] '{print $1}' > pos_mb.jo 

#join and organize all info we need.
paste cutout.jo pos_mb.jo  >step1.jo
awk -F ' ' '{print $1"_"$5" "$1" "$5" "$2" "$3" "$4}' step1.jo > step2.jo
awk -F ' ' '{print $2" "$3}' step2.jo | awk '!seen[$0]++' > chrpos.jo

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

cat Hp_noheader_05mb.jo Hp_noheader.jo | sort -n -k2 > Hp_catted_noheader.jo

#add header and paste everything together.
echo "CHR POS CHR_POS COUNT_MIN COUNT_MAJ NUM NUM Hp" > Hp_header_05mb.jo 
cat Hp_header_05mb.jo Hp_final_noheader.jo > $out-Hp_utput.txt

#get rid of intermediate files.
rm *.jo



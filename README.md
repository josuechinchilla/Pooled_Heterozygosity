 # Pooled_Heterozygosity
**ISU - Animal Breeding and Genetics Group  
Author: Josue Chinchilla-Vargas  
Created (date): 11/3/2020  
Version (of the script): 1.0  
Program (if applicable): Plink  
Program Version (if applicable): 1.9**  
#==============================================================================#  
##**Description**  
 
Calculates pooled heterozygosity (Hp) for 1 Mb windows  
#==============================================================================#  
##**Setup  
######Input files**  

.frq.counts (from Plink1.9 --freq counts)  
Format: CHR SNP A1 A2 C1 C2 G0  
 
.map used to calculate allele counts  
Format: CHR SNP DISTANCE BP  
  
######**Script works as follows:**   
./Hp.sh myfile.frq.counts myfile.map  

######**Output files**  
Hp.txt  
Format:  CHR POS CHR_POS COUNT_MIN COUNT_MAJ NUM NUM Hp  
#==============================================================================#

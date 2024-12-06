#/bin/sh
dataset=$1

set -e

#load software for vcf sort
ml VCFtools tabix

if [ -d /nfs/DATA/$dataset/GWAS/Imputed/Raw/imputation_prep/ ]; 
then
    echo "Entering directory for imputation prep..."
    cd /nfs/DATA/$dataset/GWAS/Imputed/Raw/imputation_prep/
else
    echo "Directory for imputation prep for $dataset does not exist. Creating and entering now...";
    mkdir /nfs/DATA/$dataset/GWAS/Imputed/Raw/imputation_prep/;
    cd /nfs/DATA/$dataset/GWAS/Imputed/Raw/imputation_prep/;
fi

#get cleaned genotype data
cp /nfs/DATA/$dataset/GWAS/Genotyped/Cleaned/${dataset}_genotyped_final.* ./

printf "\n\nRunning the perl script that checks genotype files and builds the plink script.\n"

#run the script that checks and builds the plink script
perl /home/clarklc/Programs/HRC-1000G-check-bim.pl -b ${dataset}_genotyped_final.bim -f ${dataset}_genotyped_final.frq -r /home/clarklc/Programs/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

#edit to run with plink2 rather than plink
sed -i 's/plink/plink2/g' Run-plink.sh

printf "\n\nRunning the plink script to edit the genotype files...\n"

#run created script with all the plink commands; breaks into files for each chr
sh Run-plink.sh

printf "\n\nNow converting genotype files to dosage format...\n"

#convert to vcf and sort
for chr in $(seq 1 22); do
plink2 --bfile ${dataset}_genotyped_final-updated-chr${chr} --recode vcf --out ${dataset}_genotyped_final-updated-chr${chr};
vcf-sort ${dataset}_genotyped_final-updated-chr${chr}.vcf | bgzip -c > ${dataset}_genotyped_final-updated-chr${chr}.vcf.gz;
done

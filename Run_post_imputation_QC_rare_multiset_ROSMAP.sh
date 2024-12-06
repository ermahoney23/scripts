#!/bin/bash
cohort=$1
set_name=$2

set -e

set_name_version=$( echo ${cohort}_${set_name} )

printf "Running post-imputation file conversion and QC for the ${set_name_version} set of ${cohort}.\n"

#load the most recent version of plink to have vcf conversion capabilities
module load PLINK/2.00-alpha1

#work in the Imputed/Raw/ folder for the cohort where the imputation results live
cd /nfs/DATA/${cohort}/GWAS/Imputed/Raw/${set_name_version}/

#remove tmp files from common variant filtering
for i in $(seq 1 22); do
    #restrict to variants with R2>=0.90 (set all FIDs to be constant since the IID had _ in it and the default delimiter for the imputed ids is a _)
    plink2 --vcf chr${i}.dose.vcf.gz --const-fid --exclude-if-info "R2<0.8" --make-bed --out  chr${i}_temp_rare --threads 8
    #get list of duplicates SNPs (they have more than 1 row in the bim file) and  remove
    awk '{ print $2 }' chr${i}_temp_rare.bim | uniq -d > chr${i}_temp_rare.dups ;
    plink2 --bfile chr${i}_temp_rare --exclude chr${i}_temp_rare.dups --make-bed --out chr${i}_temp_rare_nodups --threads 8
done

#switch back to old plink
module unload PLINK/2.00-alpha1
module load PLINK/1.9b_5.2

#create merge file (deleting old one if present)
rm merge_list.txt
for i in $( seq 2 22 ); do echo "chr${i}_temp_rare_nodups.bed   chr${i}_temp_rare_nodups.bim  chr${i}_temp_rare_nodups.fam" >> merge_list.txt ;  done

#merge individual chromosome files to be one large file
plink --bfile chr1_temp_rare_nodups --merge-list merge_list.txt --make-bed --out  ${set_name_version}_rare --threads 8

#update SNP names and add sex back into fam file
plink --bfile ${set_name_version}_rare --update-sex  /nfs/DATA/${cohort}/GWAS/Genotyped/Cleaned/${set_name_version}_genotyped_final.fam 3  --update-name /home/clarklc/Programs/HRC_snps_nodups --make-bed --out  /nfs/DATA/${cohort}/GWAS/Imputed/Cleaned/QC/${set_name_version}_rare_names

#clean-up intermediate files
rm *temp*
cd /nfs/DATA/${cohort}/GWAS/Imputed/Cleaned/QC/

# SNP filters
plink --bfile ${set_name_version}_rare_names --maf 0.005 --geno 0.05 --hwe 0.000005 --make-bed --out ${set_name_version}_rare_maf005_geno05_hwe6

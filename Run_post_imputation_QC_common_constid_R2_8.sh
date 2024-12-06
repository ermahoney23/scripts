#/bin/sh
cohort_stem=$1 #as in BIOCARD or BIOCARD_allraces
imputation_results_folder=$2 # as in /nfs/DATA/BIOCARD/GWAS/Imputed/Raw/ or /nfs/DATA/BIOCARD/GWAS/Imputed/Raw/allraces/
cleaned_imputed_folder=$3 # as in /nfs/DATA/BIOCARD/GWAS/Imputed/Cleaned/QC/
final_genotype_file=$4 # as in /nfs/DATA/BIOCARD/GWAS/Genotyped/Cleaned/BIOCARD_genotyped_final and so on... this is to pull sex from

set -e

#load the most recent version of plink to have vcf conversion capabilities
module load PLINK/2.00-alpha1

#work where the imputation results live
cd $imputation_results_folder

for i in $(seq 1 22); do
    #restrict to variants with R2>=0.80
    plink2 --vcf chr${i}.dose.vcf.gz --const-fid 0 --exclude-if-info "R2<0.8" --make-bed --out chr${i}_temp

    #get list of duplicates SNPs (they have more than 1 row in the bim file) and remove
    awk '{ print $2 }' chr${i}_temp.bim | uniq -d > chr${i}_temp.dups ;
    plink2 --bfile chr${i}_temp --exclude chr${i}_temp.dups --make-bed --out chr${i}_temp_nodups
done

#switch back to old plink
module unload PLINK/2.00-alpha1
module load PLINK/1.9b_5.2

#create merge file (emptying old version if present)
echo "" | tee merge_list.txt
for i in $( seq 2 22 ); do echo "chr${i}_temp_nodups.bed  chr${i}_temp_nodups.bim  chr${i}_temp_nodups.fam" >> merge_list.txt ; done

#merge individual chromosome files to be one large file
plink --bfile chr1_temp_nodups --merge-list merge_list.txt --make-bed --out ${cohort_stem}

#update SNP names and add sex back into fam file
plink --bfile ${cohort_stem} --update-sex ${final_genotype_file} 3 --update-name /home/clarklc/Programs/HRC_snps_nodups --make-bed --out ${cleaned_imputed_folder}/${cohort_stem}_names

cd $cleaned_imputed_folder

# SNP filters
plink --bfile ${cohort_stem}_names --maf 0.01 --geno 0.05 --hwe 0.000001 --make-bed --out ${cohort_stem}_maf01_geno05_hwe6

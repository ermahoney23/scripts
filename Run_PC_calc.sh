#!/bin/sh
filestem=$1 #give whole path and file stem, ie /nfs/DATA/A4/GWAS/Imputed/QC/A4_maf....

set -e

path=$(  echo ${filestem%/*} )
filename=$( echo ${filestem##*/} )
cohort=$( echo ${filename%%_*} )

if [[ $path == *"Genotyped"* ]]; then
    datatype=Genotyped ; 
else
    datatype=Imputed ;
fi

module purge
module load GCC/5.4.0-2.26  GSL/2.1  OpenBLAS/0.2.18-LAPACK-3.6.1

#move to right directory
cd /nfs/DATA/${cohort}/GWAS/${datatype}/Cleaned/PCs/

echo "Pruning $cohort $datatype genetic data..."

#prune and set phenotypes to 1
plink2 --bfile $filestem --indep-pairwise 200 100 0.2 --allow-no-sex --out  ${filename}
plink2 --bfile $filestem --output-missing-phenotype 1 --extract  ${filename}.prune.in --make-bed --out ${filename}_pruned

echo "Pruning complete. Running smartpca for raw ${cohort}..."

printf "genotypename: ${filename}_pruned.bed
snpname: ${filename}_pruned.bim
indivname: ${filename}_pruned.fam
evecoutname: ${filename}.pca.evec
evaloutname: ${filename}.eigenvalues
altnormstyle: NO
numoutevec: 10
numoutlieriter: 0
numoutlierevec: 10
outliersigmathresh: 6
qtmode: 0" > ${filename}_pccalc.par

smartpca -p ${filename}_pccalc.par > ${filename}_pccalc.log

echo "smartpca for ${cohort} ${datatype} genetic data complete!"

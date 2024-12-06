#!/bin/sh
set_name=$1
genotype_file_loc=$2
dosage_dest=$3
predex_dest=$4
db=$5

#Set name is the short form/prefix for files associated with this that don't have an automatic prefix (as in SNP lists and perhaps log files)
#Give full file paths for genotype files (which must have a frequency file with them of the same name),
# and dosage and predicted expression file locations (these will be created during this process.
#Must have / at the ends of destination arguments
#db must be either all (which will indicate GTEx,CommonMind,Monocyte) or a comma-separated list of which db to do (no spaces)

#####################################

set -e


#print out the command line options supplied to make sure it makes sense
printf "Calculating predicted expression from these genotype files: $genotype_file_loc \nDosage files will be created in this folder: $dosage_dest \nPredicted expression files will be saved to this folder: $predex_dest"

#be in the dosage directory
cd $dosage_dest

#start with supplied genotype files
if [ ! -f ${genotype_file_loc}.frq ]; then
    printf "Frequency file not present for these genotype files. Creating a full copy of the files plus a .frq file in the dosage directory...\n"
    #get genotype file name without path (ie everything after the last /)
    genotype_file_name=${genotype_file_loc##*/}
    #copy the files to the current directory, same file stem as before
    plink2 --bfile $genotype_file_loc --freq --make-bed --out $genotype_file_name
    #update the genotype_file_loc variable with the current directory
    genotype_file_loc=$( printf "$dosage_dest/$genotype_file_name" )
fi

# Figure out whether this should run for all db (GTEx, CommonMind, Monocyte) or just one or a few
if [ $db == "all" ] ; 
then
    db=$( echo GTEx CommonMind Monocyte )
    printf "Building predicted expression for tissues from all databases (GTEx, CommonMind, and Monocyte)\n"
elif [[ $db =~ "," ]] ;
then
    db=$( echo $db | sed 's/,/ /g' )
    printf "Building predicted expression for tissues from $db \n"
else
    printf "Building predicted expression for tissues from $db only\n"
fi

#get the SNP lists for filtering
for i in $db ; do Rscript Get_SNP_lists_to_filter.R $set_name $genotype_file_loc $dosage_dest $i ; done

#filter SNPs
plink2 --bfile $genotype_file_loc --extract ${set_name}_SNPs_to_keep.txt --allow-no-sex --make-bed --out  ${set_name}_filtered



#### break here to send to slurm



# flip the strand
plink2 --bfile ${set_name}_filtered --flip ${set_name}_SNPs_flip.txt --allow-no-sex --make-bed --out ${set_name}_filtered_flipped

# force the right reference allele
plink2 --bfile ${set_name}_filtered_flipped --reference-allele ${set_name}_SNPs_recode.txt --allow-no-sex --make-bed --out ${set_name}_filtered_flipped_refallele

echo "Beginning conversion from plink to dosage files"

#convert from plink to dosage files
python /scratch/mahone1/PrediXcan/Scripts/convert_plink_to_dosage.py -b ${set_name}_filtered_flipped_refallele -o chr

## predict expression
#get the correct db file
for i in $db ;
    do
    if [ $i == "GTEx" ] ; 
    then
	for tissue in $( echo  Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood Brain_Amygdala Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Minor_Salivary_Gland ) ;
	do
	python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${dosage_dest}${set_name}_filtered_flipped_refallele.fam --weights /scratch/clarklc/PrediXcan/Databases/GTEx-V7_HapMap-2017-11-29/gtex_v7_${tissue}_imputed_europeans_tw_0.5_signif.db --output_prefix ${predex_dest}/${tissue}
	done
    elif [ $i == "CommonMind" ] ; 
	then 
	python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${genotype_file_loc}.fam --weights /scratch/clarklc/PrediXcan/Databases/CommonMind/DLPFC_newMetax.db --output_prefix ${predex_dest}/CommonMind
    elif [ $i == "Monocyte" ] ;
	then
	python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${genotype_file_loc}.fam --weights /scratch/clarklc/PrediXcan/Databases/MESA/CAU_imputed_10_peer_3_pcs_v2.db --output_prefix ${predex_dest}/Monocyte
    else 
	printf "Either database name is incorrect or script has not been updated to handle this database: $i"
    fi
done

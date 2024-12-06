#!/bin/sh
set_name=$1
genotype_file_loc=$2
dosage_dest=$3
predex_dest=$4
db=$5

#This is a continuation of what has been done before, but set up to go onto slurm since it may take several hours.
#db will be just one database name, because the slurm script will be created in a for loop of the db names


#####################################

set -e

printf "All filtering should have been done already to the genotype files to prepare them for prediction. Now creating dosage files for the $db database(s) from these final, filtered genotype files: ${dosage_dest}${set_name}_filtered_flipped_refallele\n"

#be in the dosage directory
cd $dosage_dest

#things that have changed with variables:
#genotype_file_loc is potentially updated
#db is potentially updated

echo "Beginning conversion from plink to dosage files"

#convert from plink to dosage files
python /scratch/mahone1/PrediXcan/Scripts/convert_plink_to_dosage.py -b ${set_name}_filtered_flipped_refallele -o chr -p plink

printf "Dosage files finished. Now predicting expression...\n"

## predict expression with the correct db file
if [ $db == "GTEx" ] ; 
then
    printf "Predicting expression for GTEx tissues...\n"
    for tissue in $( echo  Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood Brain_Amygdala Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Minor_Salivary_Gland ) ;
    do
	python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${dosage_dest}${set_name}_filtered_flipped_refallele.fam --weights /scratch/clarklc/PrediXcan/Databases/GTEx-V7_HapMap-2017-11-29/gtex_v7_${tissue}_imputed_europeans_tw_0.5_signif.db --output_prefix ${predex_dest}/${tissue}
    done

elif [ $db == "CommonMind" ] ; 
then 
    printf "Predicting expression in CommonMind...\n"
    python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${genotype_file_loc}.fam --weights /scratch/clarklc/PrediXcan/Databases/CommonMind/DLPFC_newMetax.db --output_prefix ${predex_dest}/CommonMind

elif [ $db == "Monocyte" ] ;
then
    printf "Predicting expression in Monocyte...\n"
    python /scratch/clarklc/PrediXcan/Scripts/PrediXcan_update.py --predict --dosages $dosage_dest --dosages_prefix chr --samples ${genotype_file_loc}.fam --weights /scratch/clarklc/PrediXcan/Databases/MESA/CAU_imputed_10_peer_3_pcs_v2.db --output_prefix ${predex_dest}/Monocyte

else 
    printf "Either database name is incorrect or script has not been updated to handle this database: $db \n"
fi

printf "Predicted expression for $db tissue(s) built for the $genotype_file_loc file set in this folder: $predex_dest.\n"

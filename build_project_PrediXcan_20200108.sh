#!/bin/sh
set_name=$1
genotype_file_loc=$2
predex_dest=$3
db=$4

#Set name is the short form/prefix for files associated with this that don't have an automatic prefix (as in SNP lists and perhaps log files)
#Give full file paths for genotype files,
# and predicted expression file location (these will be created during this process).
# This will be for most cohorts something like /scratch/mahone1/PrediXcan/COHORT_NAME/Predicted_Expression/.
#  The assumption is that this is the top directory to work from, in which to build all the predicted expression.
#Must have / at the ends of destination arguments
#db must be either all (which will indicate GTEx,CommonMind,Monocyte) or a comma-separated list of which db to do (no spaces)


#Updates

#on 2020/01/08
# updated the assumption about the dosage folder, so that whatever is supplied as the dosage folder is used as the dosage folder, rather than there being an assumption of it changing


#####################################

set -e

#print out the command line options supplied to make sure it makes sense
printf "Calculating predicted expression from these genotype files: $genotype_file_loc \n This is the directory in which predicted expression will be built and saved: $predex_dest\n"

# Figure out whether this should run for all db (GTEx, CommonMind, Monocyte) or just one or...
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


#prep the directory structure
##Predicted_Expression folder
if [ ! -d "${predex_dest}" ];
then
    printf "${predex_dest} does not exist. Creating now...\n"
    mkdir ${predex_dest}
fi
 
##db folders
for i in $db ; do
    if [ ! -d "${predex_dest}/${i}/" ];
    then
	printf "Folder for ${i} does not exist. Creating now...\n"
	mkdir ${predex_dest}/${i}/
    elif [ "$(ls -A ${predex_dest}/${i}/)" ];
    then
	printf "Careful! ${predex_dest}/${i}/ is not empty.\n"
    fi
done

##dosage folders
for i in $db ; do
    if [ ! -d "${predex_dest}/${i}/dosages/" ];
    then
        printf "Dosage folder for ${i} does not exist. Creating now...\n"
        mkdir ${predex_dest}/${i}/dosages/
    elif [ "$(ls -A ${predex_dest}/${i}/dosages/)" ];
    then
	printf "Careful! ${predex_dest}/${i}/dosages/ is not empty.\n"
    fi
done

#be in the Predicted_Expression folder for the cohort cohort
cd ${predex_dest}/

#Check if there is a frequency file with the genotype files. Also check if the genotype files are on scratch. If not (for either), create a copy with frequency file on scratch.
if [ ! -f ${genotype_file_loc}.frq ] || [[ ! $genotype_file_loc =~ "scratch" ]]; then
    printf "Frequency file not present for these genotype files or files are not on scratch (hence inacessible from a slurm job). Creating a full copy of the files plus a .frq file in the dosage directory...\n"
    #get genotype file name without path (ie everything after the last /)
    genotype_file_name=${genotype_file_loc##*/}
    #copy the files to the current directory, same file stem as before
    plink --bfile $genotype_file_loc --freq --allow-no-sex --make-bed --out $genotype_file_name
    #update the genotype_file_loc variable with the current directory
    genotype_file_loc=$( printf "$predex_dest/$genotype_file_name" )
fi


# Prepare genotype files and submit job to slurm -- for each database
for i in $db ; 
do 

    printf "Preparing for filtering of $genotype_file_loc file set for use with the $i database.\n"

    #create variable for the dosage folder for ease of reference
    dosage_dest=$( echo "${predex_dest}/${i}/dosages/" )

    #get the SNP lists for filtering
    Rscript /scratch/mahone1/scripts/Get_SNP_lists_to_filter.R $set_name $genotype_file_loc $dosage_dest $i ; 

    cd $dosage_dest

    #filter SNPs
    plink --bfile $genotype_file_loc --extract ${set_name}_SNPs_to_keep.txt --allow-no-sex --make-bed --out  ${set_name}_filtered

    # flip the strand
    plink --bfile ${set_name}_filtered --flip ${set_name}_SNPs_flip.txt --allow-no-sex --make-bed --out ${set_name}_filtered_flipped

    # force the right reference allele
    plink --bfile ${set_name}_filtered_flipped --reference-allele ${set_name}_SNPs_recode.txt --allow-no-sex --make-bed --out ${set_name}_filtered_flipped_refallele

    printf "Filtering complete for ${i}! Now submitting the script which converts final genotype files to dosage format and predicts expression to the scheduler...\n"

    #edit the slurm script to be for this database
    script_name=$( echo ${predex_dest}build_${set_name}_PrediXcan_${i}.slurm )
    db_predex_dest=$( echo ${predex_dest}/${i}/ )
    sed -e "s|set_name|$set_name|g" -e "s|genotype_file_loc|$genotype_file_loc|g" -e "s|dosage_dest|$dosage_dest|g" -e "s|predex_dest|$db_predex_dest|g" -e "s|db|$i|g" /scratch/mahone1/scripts/build_project_PrediXcan.slurm > $script_name

    sbatch $script_name

done


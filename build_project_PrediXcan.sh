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
#dosage folder should be /path/to/project/dosages/; then within that folder will be GTEx/ CommonMind/ and Monocyte/ folders.


#####################################

set -e

#print out the command line options supplied to make sure it makes sense
printf "Calculating predicted expression from these genotype files: $genotype_file_loc \nDosage files will be created in this folder: $dosage_dest \nPredicted expression files will be saved to this folder: $predex_dest\n"

#Before starting, check if the place for the predicted expression files exists. If not, create.
if [ ! -d $predex_dest ];
then
    printf "Supplied directory for final predicted expression files (${predex_dest}) does not currently exist. Creating now...\n"
    mkdir $predex_dest
elif [ "$(ls -A $predex_dest)" ]; 
then
     printf "Careful! $predex_dest is not empty."
fi

#check if top dosage directory exists
if [ ! -d $dosage_dest ];
then
    printf "Supplied directory for dosages (${dosage_dest}) does not currently exist. Creating now...\n"
    mkdir $dosage_dest
fi

#be in the top dosage directory
cd $dosage_dest

#Check if there is a frequency file with the genotype files. Also check if the genotype files are on scratch. If not (for either), create a copy with frequency file on scratch.
if [ ! -f ${genotype_file_loc}.frq ] || [[ ! $genotype_file_loc =~ "scratch" ]]; then
    printf "Frequency file not present for these genotype files or files are not on scratch (hence inacessible from a slurm job). Creating a full copy of the files plus a .frq file in the dosage directory...\n"
    #get genotype file name without path (ie everything after the last /)
    genotype_file_name=${genotype_file_loc##*/}
    #copy the files to the current directory, same file stem as before
    plink --bfile $genotype_file_loc --freq --allow-no-sex --make-bed --out $genotype_file_name
    #update the genotype_file_loc variable with the current directory
    genotype_file_loc=$( printf "$dosage_dest/$genotype_file_name" )
fi

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


for i in $db ; 
do 

printf "Preparing for filtering of $genotype_file_loc file set for use with the $i database.\n"

#reset the dosage destination, just in case multiple databases are being run; drop everything before "dosages/" in the dosage_dest variable and then add back dosages/ so that it will be ready to add on the db name for the correct file path
dosage_dest=$( echo ${dosage_dest%*dosages/*}dosages/ )

#update the dosage destination variable to have the db name
dosage_dest=$( printf "${dosage_dest}${i}/" );
#check to see if the directory does not exist and, if so, create it
if [ ! -d "$dosage_dest" ] 
then
    printf "Directory for dosages (${dosage_dest}) does not exist. Creating now...\n" 
    mkdir $dosage_dest
fi

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
sed -e "s|set_name|$set_name|g" -e "s|genotype_file_loc|$genotype_file_loc|g" -e "s|dosage_dest|$dosage_dest|g" -e "s|predex_dest|$predex_dest|g" -e "s|db|$i|g" /scratch/mahone1/scripts/build_project_PrediXcan.slurm > $script_name

sbatch $script_name

done


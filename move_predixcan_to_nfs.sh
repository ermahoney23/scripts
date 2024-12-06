#!/bin/sh

set -e
#copying over predixcan results to /nfs/, zipping dosages

cohort=$1

printf "This script zips dosage files and moves predicted expression files and the zipped dosages from /scratch to /nfs. It does remove source files from scratch, but is done using rsync, so if there is a failure, it should be recoverable.\n Current cohort is ${cohort}. Predicted expression will be copied from /scratch/mahone1/PrediXcan/${cohort}/Predicted_Expression/ to /nfs/DATA/${cohort}/PrediXcan/.\n\n"

#zip dosages, removing the original
for i in GTEx CommonMind Monocyte ; 
do
printf "Zipping ${cohort} dosgaes for the ${i} database\n"
cd /scratch/mahone1/PrediXcan/${cohort}/Predicted_Expression/${i}
tar -zcvf ${cohort}_${i}_dosages.tar.gz dosages/ --remove-files
done

#now copy all the files over to /nfs/
cd /scratch/mahone1/PrediXcan/${cohort}/
for i in GTEx CommonMind Monocyte ; 
do
printf "Copying $i from scratch to nfs\n"
rsync -avh --remove-source-files Predicted_Expression/${i}/ /nfs/DATA/${cohort}/PrediXcan/${i}/ 
done


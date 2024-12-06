#!/bin/bash
#Generates Polygenic Risk Scores
#Logan Dumitrescu, July 2019

#--bind /cog/murat/for_logan/April2019:/new_folder_1

#Set variables

        GENETICS="/new_folder_1/michigan_imputation_server/merged_results/BLSA_HRC_merged_final_noOutliers"
        CLUMP1=.01    #Significance threshold for index SNPs
        CLUMP2=.5     #Secondary significance threshold for clumped SNP
        CLUMP3=.25    #LD threshold for clumping
        CLUMP4=200    #Physical distance threshold for clumping

#################
#CSF Ab42

        SUMSTATS="/scripts/Summary_Stats/CSF_AB42_GWAS_summary_results.txt"
        OUTCOME=CSF_Ab42

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
fi

#################

#CSF Ab42 – no APOE

        SUMSTATS="/scripts/Summary_Stats/CSF_AB42_GWAS_summary_results.txt"
        OUTCOME=CSF_Ab42_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
fi

#################
#CSF Tau

        SUMSTATS="/scripts/Summary_Stats/CSF_tau_GWAS_summary_results.txt"
        OUTCOME=CSF_Tau

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
fi

#################

#CSF Tau – no APOE

        SUMSTATS="/scripts/Summary_Stats/CSF_tau_GWAS_summary_results.txt"
        OUTCOME=CSF_Tau_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
fi


#################
#CSF Ptau

        SUMSTATS="/scripts/Summary_Stats/CSF_ptau_GWAS_summary_results.txt"
        OUTCOME=CSF_Ptau

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps
fi

#################

#CSF Ptau – no APOE

        SUMSTATS="/scripts/Summary_Stats/CSF_ptau_GWAS_summary_results.txt"
        OUTCOME=CSF_Ab42_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi


#################
# Educ

        SUMSTATS="/scripts/Summary_Stats/GWAS_EA_excl23andMe.txt"
        OUTCOME=Educ

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_Educ.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################

#Educ – no APOE

        SUMSTATS="/scripts/Summary_Stats/GWAS_EA_excl23andMe.txt"
        OUTCOME=Educ_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_Educ.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################
#HIPV

        SUMSTATS="/scripts/Summary_Stats/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL_edited"
        OUTCOME=HIPV

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_HIPV.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi


#################

#HIPV – no APOE

        SUMSTATS="/scripts/Summary_Stats/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL_edited"
        OUTCOME=HIPV_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_HIPV.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################
#IGAP

        SUMSTATS="/scripts/Summary_Stats/IGAP_stage_1.txt"
        OUTCOME=IGAP

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_IGAP.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################

#IGAP – no APOE

        SUMSTATS="/scripts/Summary_Stats/IGAP_stage_1.txt"
        OUTCOME=IGAP_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_IGAP.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################
#PGC

        SUMSTATS="/scripts/Summary_Stats/AD_sumstats_Jansenetal.txt"
        OUTCOME=PGC

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_PGC.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################

#PGC – no APOE

        SUMSTATS="/scripts/Summary_Stats/AD_sumstats_Jansenetal.txt"
        OUTCOME=PGC_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_PGC.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi


#################
#UKB_LHIPV

        SUMSTATS="/scripts/Summary_Stats/ukb_roi_volume_dec12_2018_phase1and2_pheno81_allchr_withA2.csv"
        OUTCOME=UKB_LHIPV

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################

#UKB_LHIPV – no APOE

        SUMSTATS="/scripts/Summary_Stats/ukb_roi_volume_dec12_2018_phase1and2_pheno81_allchr_withA2.csv"
        OUTCOME=UKB_LHIPV_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}


if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################
#UKB_RHIPV

        SUMSTATS="/scripts/Summary_Stats/ukb_roi_volume_dec12_2018_phase1and2_pheno94_allchr_withA2.csv"
        OUTCOME=UKB_RHIPV

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

#################

#UKB_RHIPV – no APOE

        SUMSTATS="/scripts/Summary_Stats/ukb_roi_volume_dec12_2018_phase1and2_pheno94_allchr_withA2.csv"
        OUTCOME=UKB_RHIPV_noAPOE

#Run R script which 1) removes ambiguous SNPs and 2) generates list of overlapping SNPs

        Rscript /scripts/Determine_overlapping_SNPs_CSF.R $SUMSTATS $GENETICS $OUTCOME

#Extract overlapping SNPs from the genotype data

        plink19 --bfile $GENETICS --allow-no-sex --extract ${OUTCOME}_overlapping_SNPs.txt --exclude range /scripts/APOErange.txt --make-bed --out $OUTCOME

#Perform LD clumping
#NOTE: If just rerunning analyses using different clumping criteria, start script here and comment out previous steps

        plink19 --bfile $OUTCOME --allow-no-sex --clump ${OUTCOME}_summary_stats_updated.txt --clump-p1 $CLUMP1 --clump-p2 $CLUMP2 --clump-r2 $CLUMP3 --clump-kb $CLUMP4 --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}

#Pull list of SNPs from clumped file

        awk '{print $3}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.clumped > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.SNPlist

#Merge with updated summary statistics file and generate PRS input file

        Rscript /scripts/Generate_score_input_file.R $OUTCOME $CLUMP1 $CLUMP2 $CLUMP3 $CLUMP4

#Create .profile file

        plink19 --bfile $OUTCOME --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}


if [ -f  ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred ];
then

printf "\n\nSome variants had mismatched alleles for the PRS. Attempting to flip the strands and rebuild the score.\n"

#Pull out list of SNPs with allele code mismatches

        awk '{print $2}' ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.nopred > ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps

#Flip strand for list of SNPs

        plink19 --bfile $OUTCOME --allow-no-sex --flip ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.flipsnps --make-bed --out ${OUTCOME}_flipsnps

#Create .profile file

        plink19 --bfile ${OUTCOME}_flipsnps --allow-no-sex --score ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}.score --out ${OUTCOME}_LDclump_${CLUMP1}_${CLUMP2}_${CLUMP3}_${CLUMP4}_flipsnps

fi

###
#Clean up

rm *.bim
rm *.bed
rm *.fam

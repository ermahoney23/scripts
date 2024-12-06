#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --output=predex_destset_name_PrediXcan_db_%A.log
#SBATCH --mail-user=emily.mahoney@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=set_name_PrediXcan_db
#SBATCH --time=4-12:00:00
#SBATCH --account=h_vmac

sh /scratch/mahone1/scripts/build_project_PrediXcan_for_slurm.sh set_name genotype_file_loc dosage_dest predex_dest db

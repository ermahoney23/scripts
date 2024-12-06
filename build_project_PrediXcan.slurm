#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --output=predex_destset_name_PrediXcan_db_%A.log
#SBATCH --mail-user=emily.mahoney@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --job-name=set_name_PrediXcan_db
#SBATCH --time=4-12:00:00
#SBATCH --account=h_vmac

module purge
module load PLINK/1.9b_5.2 Intel/2017.4.196 Python/2.7.14 IntelMPI/2017.3.196 numpy/1.13.1-Python-2.7.14

sh /scratch/mahone1/scripts/build_project_PrediXcan_for_slurm.sh set_name genotype_file_loc dosage_dest predex_dest db

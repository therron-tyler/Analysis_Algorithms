#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 46:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=180gb
#SBATCH --job-name=rdata2rds
#SBATCH -N 1
#SBATCH -n 10


module load R/4.2.0

Rscript /home/ttm3567/63_tylert/Analysis_Algorithms/Convert_Rdata_2_RDS.R

#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=180gb
#SBATCH --job-name=10x2seurat
#SBATCH -N 1
#SBATCH -n 10


module load R/4.2.0

Rscript /home/ttm3567/63_tylert/Analysis_Algorithms/Merge_2_Seurat_Objects.R /home/ttm3567/rootdir_scratch/Day0_BMC_STIA_Macs_MERGED/outs/filtered_feature_bc_matrix /home/ttm3567/rootdir_scratch/Day7_BMC_STIA_Macs_MERGED/outs/filtered_feature_bc_matrix Day0_BMC_STIA_Macs Day7_BMC_STIA_Macs Day0_and_Day7_STIA_Macs_MERGED /home/ttm3567/63_tylert/

#!/bin/sh

#SBATCH --job-name=QuantumClone
#SBATCH --ntasks-per-node=6
#SBATCH --time=20:00:00
#SBATCH -o RunClustering.out
#SBATCH -e RunClustering.err
#SBATCH --mem=150gb

ml R/4.4.1-gfbf-2023b
module load R-bundle-CRAN/2024.06-foss-2023b
cd /path/to/quantumclonefolder
Rscript "/path/to/quantumclonefolder/clustering.R"
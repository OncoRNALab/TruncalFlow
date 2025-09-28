#!/bin/sh

#SBATCH --job-name=QuantumClone_run2
#SBATCH --ntasks-per-node=12
#SBATCH --time=20:00:00
#SBATCH -o RunClustering_all.out
#SBATCH -e RunClustering_all.err
#SBATCH --mem=200gb

ml R/4.4.1-gfbf-2023b
module load R-bundle-CRAN/2024.06-foss-2023b
cd /path/to/quantumclonecnvfolder
Rscript "/path/to/quantumclonecnvfolder/clustering_sample.R"
#!/bin/sh

#SBATCH --job-name=SciClone
#SBATCH --ntasks-per-node=6
#SBATCH --time=20:00:00
#SBATCH -o SciClone.out
#SBATCH -e SciClone.err
#SBATCH --mem=150gb

ml R/4.4.1-gfbf-2023b
module load R-bundle-CRAN/2024.06-foss-2023b 
module load R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1
cd /path/to/sciclonefolder
Rscript "/path/to/sciclonefolder/clustering.R"
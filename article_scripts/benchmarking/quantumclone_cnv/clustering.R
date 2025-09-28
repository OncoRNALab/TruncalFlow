library(ggplot2)
library(doSNOW)
library(fpc)
library(QuantumClone)

setwd('/path/to/quantumclonecnvfolder')

current_dir <- getwd()
samples <- list.dirs(current_dir, full.names = TRUE, recursive = FALSE)
print(samples)

for (sample in samples) {

  folder_name <- basename(sample)
  
  setwd(sample)
  file <- list.files(sample, pattern = "_quantumclone\\.txt$", full.names = TRUE)
  file_freec <- list.files(sample, pattern = "_freec\\.txt$", full.names = TRUE)

  if (length(file) == 1) {

    vcf <- read.table(file, header = TRUE, sep='\t')
    freec <- read.table(file_freec, header = TRUE, sep='\t')

    vcf$Chr <- as.integer(vcf$Chr)
    vcf$Start <- as.integer(vcf$Start)
    vcf$Depth <- as.numeric(vcf$Depth)
    vcf$Alt <- as.numeric(vcf$Alt)
    #vcf$Genotype <- as.factor(vcf$Genotype)
    
    vcf_list <- list(vcf)
    freec_list <- list(freec)
    
    output_dir <- sample
    
    clustering <- One_step_clustering(vcf_list, FREEC_list = freec_list, contamination = 0, 
                                  nclone_range = 2:5, clone_priors = NULL, prior_weight = NULL, 
                                  Initializations = 1, preclustering = "FLASH", simulated = FALSE, 
                                  epsilon = NULL, save_plot = TRUE, ncores = 1, restrict.to.AB = FALSE, 
                                  output_directory = NULL, model.selection = "BIC", optim = "default", 
                                  keep.all.models = FALSE, force.single.copy = FALSE)
    df_clusters <- data.frame(
    Position = seq_along(clustering$cluster),
    Number = clustering$cluster)
    write.csv(df_clusters, "clustering.csv")

    fik_df <- as.data.frame(clustering$EM.output$fik)
    write.csv(fik_df, "FIK.csv")
    
    weights_df <- data.frame(weights = clustering$EM.output$weights)
    write.csv(weights_df, "weights.csv")
    
    val_df <- data.frame(val = clustering$EM.output$val)
    write.csv(val_df, "VAL.csv")
    
    centers_df <- do.call(rbind, lapply(clustering$EM.output$centers, as.data.frame))
    centers_df <- data.frame(centers_df)
    write.csv(centers_df, "centers.csv")
    
    normalized_centers_df <- do.call(rbind, lapply(clustering$EM.output$normalized.centers, as.data.frame))
    normalized_centers_df <- data.frame(normalized_centers_df)
    write.csv(normalized_centers_df, "normalized_centers.csv")
    
    filtered_df <- do.call(rbind, lapply(clustering$filtered.data, as.data.frame))
    filtered_df <- data.frame(filtered_df)
    write.csv(filtered_df, "filtered.csv")
  }
}




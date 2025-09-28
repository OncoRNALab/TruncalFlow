library(sciClone)

setwd('/path/to/sciclonefolder')

current_dir <- getwd()
samples <- list.dirs(current_dir, full.names = TRUE, recursive = FALSE)

for (sample in samples) {

  folder_name <- basename(sample)

  if (folder_name %in% c("sampleP1", "sampleP2", "sampleP3")) {
    next
  }
  
  setwd(sample)

  file_vaf <- list.files(sample, pattern = "_sciclone\\.txt$", full.names = TRUE)
  file_loh <- list.files(sample, pattern = "_loh\\.txt$", full.names = TRUE)
  file_cn <- list.files(sample, pattern = "_copynumber\\.txt$", full.names = TRUE)

  vaf <- read.table(file_vaf, header = TRUE, sep='\t')

  vaf <- vaf[c("CHROM", "POS", "TUM1_REF_READS", "TUM1_VAR_READS", "TUM1_VAF")]

  #read in regions to exclude (commonly LOH)
  #format is 3-col bed
  regions <- read.table(file_loh, header=T)

  #read in segmented copy number data
  #4 columns - chr, start, stop, segment_mean   
  cn <- read.table(file_cn, header=T)

  sc <- sciClone(vafs=vaf,
         copyNumberCalls=cn,
         sampleNames=sample,
         regionsToExclude=regions)

  output_dir <- "results"

  dir.create(file.path(sample, output_dir))

  writeClusterTable(sc, file.path(sample, output_dir, "clusters.txt"))

}
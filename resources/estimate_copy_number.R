#!/usr/bin/Rscript

library(magrittr)
library(QDNAseq)
library(QDNAseq.hg38)
library(ACE)

#################### READ INPUTS ################################

# bam_dir <- commandLineArgs(sys.argv[0])
args <- commandArgs(trailingOnly = TRUE)
bam_dir <- args[1]
binsize <- 50
bins <- getBinAnnotations(binSize=binsize, genome="hg38")
readCounts <- binReadCounts(bins, path = bam_dir)
filename <- tools::file_path_sans_ext(bam_dir)


################## QDNASEQ PIPELINE #############################

# Segmentation uses (sqrt(x + 3/8)) to stabilise variance 
# of a Poisson distribution.
# Blacklist filter set to 100 to not restrict filtering
# and allow overlap initially.
copy_Numbers_segmented <- readCounts %>%
  applyFilters(residual=TRUE, blacklist=TRUE, chromosomes=NA) %>%
  
  # Estimate and correct read counts as a function of GC content and
  # mappability (removes further noise from e.g. FFPE)
  estimateCorrection() %>%
  correctBins() %>%
  
  # Normalise and smooth corrected read counts
  smoothOutlierBins()  %>%
  normalizeBins() %>%
  segmentBins(transformFun="sqrt")

# Generate called read count plots
# for final CNV call profiles
copyNumbersCalled <- readCounts %>%
  
  # Filter reads based on a) residuals after fitting to a model, 
  # and b) on amount  of acceptable overlap with 
  # blacklisted (i.e. problematic) regions/known artifacts
  applyFilters(residual=TRUE, blacklist=TRUE, chromosomes=NA) %>%
  
  # Estimate and correct read counts as a function of GC content and
  # mappability (removes further noise from e.g. FFPE)
  estimateCorrection() %>%
  correctBins() %>%
  
  # Normalise and smooth corrected read counts
  normalizeBins() %>%
  smoothOutlierBins() %>%
  
  # Segment bins: i.e. group similar regions together
  # in a single "line" or "window" in terms of read count value
  segmentBins(transformFun="sqrt") %>%
  normalizeSegmentedBins() %>%
  callBins()

# Save estimated copy number plots
#PNG naming needs fixing - "%s_" generates error
#Error in png(filename = "%s_called_copy_numbers.png", width = 1300, height = 800) : 
#  invalid 'filename'
png(filename = paste0(filename, ".copy_number_segmented.png"), width = 1300, height = 800)
plot(copy_Numbers_segmented)
dev.off()

png(filename = paste0(filename, ".called_copy_numbers.png"), width = 1300, height = 800)
plot(copyNumbersCalled)
dev.off()

# Save tabular data

## Summarised calls - one file per sample
QDNAseq::exportBins(copyNumbersCalled,
                    format = "seg", 
                    type = "calls",
                    file = paste("%s.summarised_calls_binsize_",binsize,".tsv", sep = ""))

## Raw calls - one file that contains *all* samples together 
### MAY NOT NEED THIS FUNCTION
QDNAseq::exportBins(copyNumbersCalled, 
                    format = "tsv", 
                    type = "calls",
                    file = paste(filename, ".all_calls_binsize_",
                                 binsize, "_merged.raw_calls_all_depths.tsv",
                                 sep = ""))

################## ACE FUNCTIONS #############################

generate_sample_model <- function(object, QDNAseqobjectsample){
  ploidy <- get_ploidy_for_sample(QDNASeqobjectsample)
  # Perform model fitting
  model1 <- singlemodel(object, QDNAseqobjectsample = QDNAseqobjectsample)
  # Select best fit by examining error lists (cellularity where relative error is lowest)
  bestfit1 <- model1$minima[tail(which(model1$rerror==min(model1$rerror)), 1)]
  besterror1 <- min(model1$rerror)
  lastfit1 <- tail(model1$minima, 1)
  lasterror1 <- tail(model1$rerror, 1)
}


generate_absolute_plot <- function(object, QDNAseqobjectsample) {
  # Use best fit for absolute copy number plot of sample
  absolute_plot <- singleplot(object, QDNAseqobjectsample = QDNAseqobjectsample, cellularity = bestfit1,
                              error = besterror1, standard = model1$standard,
                              title = paste0(filename, "- binsize", binsize, "kbp_- 2N fit 1"))
  dev.off()
}


save_absolute_plot <- function(absolute_plot){
  png(paste0(filename,"_plot.png"))
  ggsave(paste0(filename, "_absoluteCN_plot.png"), absolute_plot, limitsize = FALSE)
  dev.off()
}


generate_error_plots <- function(object, QDNAseqobjectsample) {
  error_plot <- print(model1$errorplot + ggtitle(paste0("sample ", QDNAseqobjectsample, " - errorlist")) +
                        theme(plot.title = element_text(hjust = 0.5)))
}


save_error_plot <- function(error_plot){
  png(paste0(filename,"_errorplot.png"))
  ggsave(paste0(filename, "_errorplot.png"), error_plot, limitsize = FALSE)
  dev.off()
}


generate_CN_dfs <- function(object, QDNAseqobjectsample) {
  # generate all absolute copy numbers df
  template <- objectsampletotemplate(object, index = QDNAseqobjectsample)
  write.table(template, file = paste0(filename,"absolute_copy_number",".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # generate segmented copy number df
  adjusted_segments <- getadjustedsegments(template, log = TRUE)
  write.table(adjusted_segments, file = paste0(filename,"_adjusted_segments.tsv"),
              sep =    "\t", quote = FALSE, row.names = FALSE)
}
  

save_CN_dfs <- function(template, adjusted_segments){
  write.table(template, file = paste0(filename,"absolute_copy_number",".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(adjusted_segments, file = paste0(filename,"_adjusted_segments.tsv"),
              sep =    "\t", quote = FALSE, row.names = FALSE)
}

################## ACE WORKFLOW #############################

run_ace_pipeline <- function(object, QDNAseqobjectsample){
  generate_sample_model()
  generate_absolute_plot()
  save_absolute_plot()
  generate_error_plots()
  save_error_plot()
  generate_CN_dfs()
  save_CN_dfs()
}

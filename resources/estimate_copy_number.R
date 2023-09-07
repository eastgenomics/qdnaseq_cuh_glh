#!/usr/bin/Rscript

library(magrittr)
library(QDNAseq)
library(QDNAseq.hg38)
library(ACE)
library(ggplot2)

#################### READ INPUTS ################################

# This script will use qDNASeq and ACE and generate outputs on per
# sample basis usiung hg38 reference genome. This is because uploading 
# .bam files takes alot of memory and allows separate outputs per sample 
# should any given sample on a run fail (e.g. fail QC metrics)
# Argument input follows a specific order

# args 1 = BAM directory
# args 2 = Selected binsize 

# bam_dir <- commandLineArgs(sys.argv[0])
args <- commandArgs(trailingOnly = TRUE)
bam_dir <- args[1]
binsize <- args[2]
ploidy <- args[3]
bins <- getBinAnnotations(binSize=binsize, genome="hg38")
readCounts <- binReadCounts(bins, path = bam_dir)
files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
filename <- tools::file_path_sans_ext(basename(files))
QDNAseqobjectsample <- 1

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

generate_sample_model <- function(object, QDNAseqobjectsample, filename, ploidy){
  # Perform model fitting on a single sample
  model1 <- singlemodel(object, QDNAseqobjectsample = QDNAseqobjectsample, ploidy = ploidy)
  # Select best fit by examining error lists (cellularity where relative error is lowest)
  bestfit1 <- model1$minima[tail(which(model1$rerror==min(model1$rerror)), 1)]
  besterror1 <- min(model1$rerror)
  lastfit1 <- tail(model1$minima, 1)
  lasterror1 <- tail(model1$rerror, 1)
  error_plot <- print(model1$errorplot + ggtitle(paste0("sample ", filename, " - errorlist")) +
                        theme(plot.title = element_text(hjust = 0.5)))
  #  Use best fit and best error to select best model for absolute CN plot
  png(paste0(filename,"_errorplot.png"))
  ggsave(paste0(filename, "_errorplot.png"), error_plot, limitsize = FALSE)
  dev.off()
  
  return(model1)
}

generate_absolute_plot <- function(object, QDNAseqobjectsample, bestfit1, 
                                   besterror1, model1, filename) {
  absolute_plot <- singleplot(object, QDNAseqobjectsample = QDNAseqobjectsample, cellularity = bestfit1,
                              error = besterror1, standard = model1$standard,
                              title = paste0(filename, "- binsize", binsize, "kbp_- 2N fit 1"))
  # Create absolute CN plot
  png(paste0(filename,"_plot.png"))
  ggsave(paste0(filename, "_absoluteCN_plot.png"), absolute_plot, limitsize = FALSE)
  dev.off()
}

generate_CN_dfs <- function(object, QDNAseqobjectsample, filename) {
  template <- objectsampletotemplate(object, index = QDNAseqobjectsample)
  write.table(template, file = paste0(filename,"_absolute_copy_number.tsv"),
              sep = "\t")
  adjusted_segments <- getadjustedsegments(template, log = TRUE)
  write.table(adjusted_segments, file = paste0(filename,"_adjusted_segments.tsv"),
              sep =	"\t")
}

run_ace_pipeline <- function(object, QDNAseqobjectsample, filename, ploidy){
  model1 <- generate_sample_model(object, QDNAseqobjectsample, filename, ploidy)
  bestfit1 <- model1$minima[tail(which(model1$rerror==min(model1$rerror)), 1)]
  besterror1 <- min(model1$rerror)
  generate_absolute_plot(object, QDNAseqobjectsample, bestfit1, besterror1, model1, filename)
  generate_CN_dfs(object, QDNAseqobjectsample, filename)
}

run_ace_pipeline(copyNumbersCalled, QDNAseqobjectsample, filename)

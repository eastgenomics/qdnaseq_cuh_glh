#!/usr/bin/Rscript

library(QDNAseq)
library(magrittr)

# read bin annotations
bins <- getBinAnnotations(binSize=commandArgs()[2])

# read BAMs
readCounts <- binReadCounts(bins, path=commandArgs()[1])

# segmentation step uses a random number generator,
# so for reproducibility we should set a seed value
set.seed(5)

# call copy numbers
copyNumbersCalled <- readCounts %>%

    # Filter reads based on a) residuals after fitting to a model, 
    # and b) on amount  of acceptable overlap with 
    # blacklisted (i.e. problematic) regions/known artifacts
    applyFilters(residual=TRUE, blacklist=100) %>%

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

## write plots


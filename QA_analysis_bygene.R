library(stringr)
library(dplyr)
library(tidyr)
library(bedr)
library(DT)

######################### FUNCTIONS ####################################

# function to clean qDNASeq summarised data
make_bed_from_sum_df_qDNA <- function(SWGS_sum) {
  # make loc column
  SWGS_sum$location <- paste(SWGS_sum$CHROMOSOME, SWGS_sum$START, sep=":")
  SWGS_sum$location <- paste(SWGS_sum$location, SWGS_sum$STOP, sep="-")
  # order and make bed file
  SWGS_sum <- as.data.frame(SWGS_sum$location)
  SWGS_sum <- as.data.frame(SWGS_sum[order(SWGS_sum$`SWGS_sum$location`),])
  # remove dup rows
  SWGS_sum <- as.data.frame(SWGS_sum[!duplicated(SWGS_sum), ])
  colnames(SWGS_sum) <- "location"
  SWGS_sum <- separate(data = SWGS_sum, col = location, into = c("chr", "start", "end"))
  
  return(SWGS_sum)
}

# function to clean ACE summarised data
make_bed_from_sum_df_ACE <- function(SWGS_sum) {
  # make loc column
  SWGS_sum$location <- paste(SWGS_sum$chr, SWGS_sum$start, sep=":")
  SWGS_sum$location <- paste(SWGS_sum$location, SWGS_sum$end, sep="-")
  # order and make bed file
  SWGS_sum <- as.data.frame(SWGS_sum$location)
  SWGS_sum <- as.data.frame(SWGS_sum[order(SWGS_sum$`SWGS_sum$location`),])
  # remove dup rows
  SWGS_sum <- as.data.frame(SWGS_sum[!duplicated(SWGS_sum), ])
  colnames(SWGS_sum) <- "location"
  SWGS_sum <- separate(data = SWGS_sum, col = location, into = c("chr", "start", "end"))
  
  return(SWGS_sum)
}

bed_intersect <- function(sample_bed, truth_bed, sample, mydir, subfolder){
  # both dataframes must be chr,start,end
  # truth is expected to have X or Y chr 
  # pred is not expected to have X or Y 
  
  # finds the interaction of the sample and truth bed file
  
  # convert truth data to int, do not touch chr
  truth_bed$start <- as.integer(truth_bed$start)
  truth_bed$end <- as.integer(truth_bed$end)
  
  # convert pred data to int
  sample_bed$chr <- as.character(sample_bed$chr)
  sample_bed$start <- as.integer(sample_bed$start)
  sample_bed$end <- as.integer(sample_bed$end)
  
  # sort the bed files
  truth_bed.sort <- bedr.sort.region(truth_bed, check.chr = F, verbose = F)
  sample_bed.sort <- bedr.sort.region(sample_bed, check.chr = F, verbose = F)
  
  # remove dup rows
  truth_bed.sort <- truth_bed.sort[!duplicated(truth_bed.sort), ] 
  sample_bed.sort <- sample_bed.sort[!duplicated(sample_bed.sort), ] 
  
  # intersect
  write.table(sample_bed.sort, "sample_bed.bed", col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(truth_bed.sort, "truth_bed.bed", col.names = F, row.names = F, quote = F, sep = "\t")
  command = paste('bedtools intersect -a sample_bed.bed -b truth_bed.bed -loj > ',mydir,"/", subfolder, "/", sample, '_samplevstruth_intersect.bed', sep ="")
  system(command)
  pred.int.truth <- read.delim(paste(mydir,"/", subfolder, "/", sample, '_samplevstruth_intersect.bed', sep =""), sep = "\t", header = F)
  
  # remove the tmp bed files
  system('rm sample_bed.bed  truth_bed.bed')
  return(pred.int.truth)
}

bed_missed_in_sample <- function(sample_bed, truth_bed, sample, mydir, subfolder){
  # both dataframes must be chr,start,end
  # truth is expected to have X or Y chr 
  # pred is not expected to have X or Y 
  # finds the interaction of the sample and truth bed file
  # convert truth data to int, do not touch chr
  truth_bed$start <- as.integer(truth_bed$start)
  truth_bed$end <- as.integer(truth_bed$end)
  
  # convert pred data to int
  sample_bed$chr <- as.character(sample_bed$chr)
  sample_bed$start <- as.integer(sample_bed$start)
  sample_bed$end <- as.integer(sample_bed$end)
  
  # sort the bed files
  truth_bed.sort <- bedr.sort.region(truth_bed, check.chr = F, verbose = F)
  sample_bed.sort <- bedr.sort.region(sample_bed, check.chr = F, verbose = F)
  
  # remove dup rows
  truth_bed.sort <- truth_bed.sort[!duplicated(truth_bed.sort), ] 
  sample_bed.sort <- sample_bed.sort[!duplicated(sample_bed.sort), ] 
  
  # not intersect
  write.table(sample_bed.sort, "sample_bed.bed", col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(truth_bed.sort, "truth_bed.bed", col.names = F, row.names = F, quote = F, sep = "\t")
  command = paste('bedtools intersect -b sample_bed.bed -a truth_bed.bed -loj > ', mydir,"/", subfolder, "/",sample, '_samplevstruth_missing.bed', sep ="")
  system(command)
  pred.int.truth <- read.delim(paste(mydir,"/", subfolder, "/",sample, '_samplevstruth_missing.bed', sep =""), sep = "\t", header = F) 
  
  # remove the tmp bed files
  system('rm sample_bed.bed  truth_bed.bed')
  return(pred.int.truth)
  
}

sWGS_vs_WGS <- function(Case_WGS, metadata, sample, copynumber_df, mydir, subfolder){
  ######### Big function which will:
  # - calculate TP,FP and FN
  
  ## Input definitions:
  # Case_WGS -> sWGS segment file
  # metadata -> file that contains info on sample, tissue type, case number etc
  # sample -> case number e.g CUH123
  # copynumber_df -> sWGS dataframe 
  # mydir -> root directory of project that contains all subfolders 
  # subfolders -> output from qDNASeq using a specific depth & bin
  
  # get the summaried data for samples associated to this case
  C_samples <- unique(metadata[which(metadata$CUH_number %in% sample), "samplename"])
  C_fullsamplenames <-  unique(metadata[which(metadata$CUH_number %in% sample), "fullsamplename"])
  Case_truthvspredicted_df <- as.data.frame(matrix(NA, nrow =1,ncol=3))
  colnames(Case_truthvspredicted_df) <- c("TP", "FP", "FN")
  rownames(Case_truthvspredicted_df) <- copynumber_df
  
  # if ACE then use ACE make bed else use the qDNASeq make bed function
  print(copynumber_df)
  if(isTRUE(grepl("qDNAseq", copynumber_df, fixed = TRUE))) {
    print("qDNASeq")
    tool <- "qDNASeq"
    assign(paste(copynumber_df, "_bed", sep=""), make_bed_from_sum_df_qDNA(get(copynumber_df)))
    print(paste(copynumber_df, " has CN segements = ", nrow(make_bed_from_sum_df_qDNA(get(copynumber_df))), sep = ""))
  } else {
    print("ACE")
    tool <- "ACE"
    assign(paste(copynumber_df, "_bed", sep=""), make_bed_from_sum_df_ACE(get(copynumber_df)))
    print(paste(copynumber_df, " has CN segements = ", nrow(make_bed_from_sum_df_ACE(get(copynumber_df))), sep = ""))
  }
  # compare to the truth dataset
  # contains all intersect with TP and FP
  intersect_loj <- bed_intersect(sample_bed = get(paste(copynumber_df, "_bed", sep="")), 
                                 truth_bed = Case_WGS_location_bed, sample = sample,
                                 mydir = mydir, subfolder = subfolder)
  # row = sample, column = TP
  # if there is overlap, all rows will be coordinates
  TP <- intersect_loj[which(intersect_loj[,4] > 0),]
  Case_truthvspredicted_df[copynumber_df,1] <- nrow(TP)
  # if true, write yes as to whether this yes is for the FF or the FFPE tissue
  # first get the truth location as chr:start-end
  TP$truth_location <- paste(TP$chr.b, ":", TP$start.b, "-", TP$end.b, sep = "")
  # row = sample, column = FP
  # if FP, they will exist in 4,5,6 column (truth dataset) as -1
  Case_truthvspredicted_df[copynumber_df,2] <- nrow(intersect_loj[which(intersect_loj[,4] == "."),])
  # row = sample, column = FN
  FN <- bed_missed_in_sample(sample_bed =  get(paste(copynumber_df, "_bed", sep="")),
                             truth_bed = Case_WGS_location_bed, sample = sample,
                             mydir = mydir, subfolder = subfolder)
  Case_truthvspredicted_df[copynumber_df,3] <- nrow(FN[which(FN[,4] == "."),])
  
  datatable(Case_truthvspredicted_df)
  write.table(Case_truthvspredicted_df, 
              file = paste(mydir,"/", subfolder, "/", sample, "_", tool ,"_TP_FP_FN_gene.tsv", sep =""),
              sep = "\t", col.names = T, row.names = T, quote = F)
}


################################# DATA ANALYSIS ###################################
### Data input
# read  data in
mydir = "~/Documents/Projects/cancer/sWGS/"
metadata <- read.delim(paste(mydir, "metadata.txt", sep= "/"), sep = "\t")
metadata$fullsamplename <- paste(metadata$samplename,"_",metadata$Tissue.Type, sep="")
ucsc_hgnc <- read.delim(paste(mydir, "ucsc_hgnc", sep= "/"), sep = "\t")
colnames(ucsc_hgnc)[colnames(ucsc_hgnc) == "X.chrom"] ="chrom"

### Case 1
samples <- c()
# read in truth WGS
sample="CUH182"
# read in truth data
Case_WGS <- read.csv(paste(mydir, "/truth_dataset_gene/", paste(sample, "CNVs.csv"), sep =""))
# remove rows which are n/a or blank in first column
Case_WGS <-  Case_WGS[!(is.na(Case_WGS$Gene.s) | Case_WGS$Gene.s=="" | Case_WGS$Gene.s=="n/a"), ]

## some CN have to genes(gene1/gene2) on a row and we need to seperate it
for (i in 1:nrow(Case_WGS)) {
  if (grepl("/",Case_WGS[i,'Gene.s'], fixed = TRUE)) {
    genes_vector <- str_split_1(Case_WGS[i,'Gene.s'], "/")
    # copy the rows based on the number of genes in gene_vector
    Case_WGS_repgenes <- Case_WGS[rep(i, length(genes_vector)), ]
    Case_WGS_repgenes$Gene.s <- genes_vector
    Case_WGS <- rbind(Case_WGS, Case_WGS_repgenes)
    Case_WGS <- Case_WGS[-i,]
  } 
}

# append the gene start and end to the Case_WGS file
Case_WGS$chr <- ucsc_hgnc$chrom[match(Case_WGS$Gene.s,ucsc_hgnc$symbol)]
Case_WGS$start <- ucsc_hgnc$chromStart[match(Case_WGS$Gene.s,ucsc_hgnc$symbol)]
Case_WGS$end <-ucsc_hgnc$chromEnd[match(Case_WGS$Gene.s,ucsc_hgnc$symbol)]

# make ned format
# remove chr 
Case_WGS$chr<-gsub("chr","",as.character(Case_WGS$chr))
Case_WGS_location_bed <- Case_WGS[,c("chr", "start", "end")]



# read in qDNAseq segment files
subfolder = "qDNAseq_ACE/FF_N1/"
total_files =  list.files(path = paste0(mydir,subfolder, sep= "" ), pattern= ("*.summarised_calls_binsize_50.tsv"))
sum_df  <- as.data.frame(matrix(ncol = 6, nrow = 0))
for(i in 1:length(total_files)) { #up to 4 samples
  file <- total_files[i]
  file_name <- strsplit(file, ".summarised_calls_binsize_50.tsv") #removes the .tsv taken from file name
  tsv_file <- read.csv(paste(mydir, subfolder, "/", file, sep= ""), sep="\t", header = T) #read in the file as a dataframe
  sum_df <- rbind(sum_df, tsv_file)
  file_name <- file_name[[1]] #get the dataframe from the list
  samples <- c(gsub("_S.*", "", file_name), samples)
  file_name <- gsub("_S.*", "_qDNAseq", file_name)
  assign(file_name, tsv_file) #assigns the tsv_file to the object file_name which is the dataframe
  file_name <- NULL
  tsv_file <- NULL
}
# read in ACE segment files
total_files =  list.files(path = paste0(mydir,subfolder, sep= "" ), pattern= ("*adjusted_segments.tsv"))
sum_df  <- as.data.frame(matrix(ncol = 6, nrow = 0))
for(i in 1:length(total_files)) { #up to 4 samples
  file <- total_files[i]
  file_name <- strsplit(file, "adjusted_segments.tsv") #removes the .tsv taken from file name
  tsv_file <- read.csv(paste(mydir, subfolder, "/", file, sep= ""), sep="\t", header = T) #read in the file as a dataframe
  colnames(tsv_file) <- c("chr", "start", "end", "Num_Probes", "Segment_Mean")
  tsv_file <- tsv_file[complete.cases(tsv_file), ]
  sum_df <- rbind(sum_df, tsv_file)
  file_name <- file_name[[1]] #get the dataframe from the list
  file_name <- gsub("_S.*", "_ACE", file_name)
  assign(file_name, tsv_file) #assigns the tsv_file to the object file_name which is the dataframe
  file_name <- NULL
  tsv_file <- NULL
}

## Comparison
# truth vs qDNAseq
copynumber_df <- paste(samples, "_qDNAseq", sep="")
sWGS_vs_WGS(Case_WGS, metadata,sample, copynumber_df, mydir, subfolder)
copynumber_df <- paste(samples, "_ACE", sep="")
sWGS_vs_WGS(Case_WGS, metadata, sample, copynumber_df, mydir, subfolder)























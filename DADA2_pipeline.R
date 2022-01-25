#############################################################################
#############################################################################
###########                                                       ###########
###########       DADA2 pipeline for microbiome sequencing data   ###########
###########             Author(s): Diego Montiel Gonzalez         ###########
###########                                                       ###########
###########        Erasmus MC University Medical Centre           ###########
###########               Rotterdam, The Netherlands              ###########
###########                                                       ###########
###########                     genid@erasmusmc.nl                ###########
###########                                                       ###########
#############################################################################
#############################################################################


# Verify if library is correctly installed
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
# loading libraries
pkgTest("optparse")
pkgTest("dada2")
pkgTest("DECIPHER")
pkgTest("data.table")
pkgTest("phangorn")
pkgTest("ShortRead")


filter.taxa <- function(taxa.eHOMD.filt, taxa.levels, seqtab.nochim.filt, threshold){
  taxa.levels <- as.data.frame(taxa.levels, row.names = rownames(taxa.eHOMD.filt))
  taxa.collapse <- data.frame(rownames(seqtab.nochim.filt))
  taxa.unique <- as.vector(unique(taxa.levels$taxa.levels))
  
  for(taxa in taxa.unique)
  {
    rownames.taxa <- rownames(taxa.levels)
    seq.taxa <- rownames.taxa[taxa.levels$taxa.levels == taxa]
    tmp.taxa.seq <- seqtab.nochim.filt[,seq.taxa]
    condition <- length(seq.taxa)
    if(condition >= 2)
    {
      taxa.collapse[taxa] <- as.vector(rowSums(tmp.taxa.seq))
    }else{
      taxa.collapse[taxa] <- as.vector(tmp.taxa.seq)
    }
  }
  # discard those taxa that appear only in two or 
  rownames(taxa.collapse) <- taxa.collapse$rownames.seqtab.nochim.filt.
  taxa.collapse$rownames.seqtab.nochim.filt. <- NULL
  taxa.collapse.filt <- taxa.collapse
  table(as.vector(unlist(lapply(taxa.collapse, function(x){ length(which(x!=0))}))) > 2)
  taxa.collapse.filt <- taxa.collapse[,as.vector(unlist(lapply(taxa.collapse.filt, function(x){ length(which(x!=0))}))) > 2]
  tmp.filt <- taxa.collapse.filt[,as.vector(colSums(taxa.collapse.filt))/sum(taxa.collapse.filt) < threshold]
  taxa.collapse.filt <- taxa.collapse.filt[,as.vector(colSums(taxa.collapse.filt))/sum(taxa.collapse.filt) >= threshold]
  if(length(tmp.filt) > 0){
    if(sum(tmp.filt) > 0)
    {
      taxa.collapse.filt["other"] <- tmp.filt
    }
    
  }
  return(taxa.collapse.filt)
}


option_list = list(
  
  make_option(c("-w", "--working-dir"), type="character", default=NULL, 
              help="PATH of working directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="PATH sequencing data", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output folder", 
              metavar="character"),
  make_option(c("-p", "--preffix"), type="character", default=NULL, 
              help="Preffix for the additional files", 
              metavar="character"),
  make_option(c("-r", "--reads"), type="character", default=NULL, 
              help="Paired-end [PE] or Single-end [SE]", 
              metavar="character"),
  make_option(c("-t", "--threads"), type="character", default=NULL, 
              help="Number of threads to use", 
              metavar="character"),
  make_option(c("-d", "--db"), type="character", default=NULL, 
              help="FASTA expanded Human Oral Microbiome Database (eHOMD)",  metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd(opt$w)
path <- paste0(opt$w, "/", opt$input)
path.out <- paste0(opt$w, "/", opt$output)
preffix <- opt$preffix
reads <- opt$reads
database <- opt$db
threads <- as.integer(opt$threads)

# path <- "data/PE_reduced/" # for raw fastq sequences
# path.out <- "data"
# preffix <- "PRJNA484874"
# reads <- "PE"
# threads <- TRUE
# database <- "/media/disk1/diego/Microbiome/smoking/db/eHOMD_v15.2_assignTaxonomy.fasta"
# setwd("/media/disk1/diego/genid/smoking-microbiome/")
# path <- paste0(getwd(), "/", path)


if(reads == "PE")
{
  print("PE reads")
  fnFs <- base::sort(list.files(path, pattern="1.fastq", full.names = TRUE)) # If TRUE, the directory path is prepended
  fnRs <- base::sort(list.files(path, pattern="2.fastq", full.names = TRUE)) # fnRs for reverse sequencing reads
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #Extract sample names,
  #FILTER AND TRIM
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # place filtered forward files in 
  #a subdirectory
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) #filtered reverse files
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2,
                       compress=TRUE, multithread=threads, minLen = 100, verbose = T, truncLen=c(200,150))
  # LEARN THE ERROR RATES
  errF <- learnErrors(filtFs, multithread=threads) # forward files
  errR <- learnErrors(filtRs, multithread=threads) # reverse files
  # DEREPLICATION
  derepFs <- derepFastq(filtFs, verbose=TRUE) # Forward files, combines all identical sequencing reads into "unique sequences"
  # with a corresponding "abundance"
  derepRs <- derepFastq(filtRs, verbose=TRUE) # reverse files
  names(derepFs) <- sample.names  # name de derep-class objects by sample names
  names(derepRs) <- sample.names
  # SAMPLE INFERENCE
  dadaFs <- dada(derepFs, err=errF, multithread=threads) # the error model developed earlier calculates abundance p-values
  # for each unique seqeunce
  dadaRs <- dada(derepRs, err=errR, multithread=threads)
  # MERGE PAIRED READS
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  # CONSTRUCT SEQUENCE TABLE
  seqtab <- makeSequenceTable(mergers)
}else if(reads == "SE")
{
  print("SE reads")
  fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE)) # If TRUE, the directory path is prepended
  sample.names <- sapply(strsplit(basename(fnFs), ".filt"), `[`, 1) #Extract sample names, 
  filtFs <- file.path(path, "filtered", paste0(sample.names, ".gz")) # place filtered forward files in 
  length(fnFs)
  length(filtFs)
  out <- filterAndTrim(fnFs, filtFs, maxN = 0, maxEE = 2, truncQ =2, compress = TRUE, 
                       multithread = TRUE, verbose = TRUE, minLen = 100, maxLen = 500) 
  #plotQualityProfile(filtFs[50:51]) #visualize the quality profiles of the forward reads of the two first samples
  # LEARN THE ERROR RATES
  errF <- learnErrors(filtFs, multithread = threads) # forward files
  plotErrors(errF, nominalQ = TRUE) # estimated error rates
  # DEREPLICATION
  derepFs <- derepFastq(filtFs, verbose = TRUE) # Forward files, combines all identical sequencing reads into "unique sequences"
  names(derepFs) <- sample.names  # name de derep-class objects by sample names
  # SAMPLE INFERENCE
  dadaFs <- dada(derepFs, err = errF, multithread = threads) # the error model developed earlier calculates abundance p-values
  # CONSTRUCT SEQUENCE TABLE
  seqtab <- makeSequenceTable(dadaFs)
}

# REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)
##### ASSIGN TAXONOMY ######
#eHOMD database
set.seed(100)
taxa.eHOMD <- dada2::assignTaxonomy(seqtab.nochim, database, multithread=threads)
seqtab.nochim <- unique(seqtab.nochim) # remove duplicates
taxa.eHOMD.filt <- taxa.eHOMD
#1)	Discard those samples with less than 1,000 reads
seqtab.nochim.filt <- seqtab.nochim[,rownames(taxa.eHOMD.filt)]
seqtab.nochim.filt <- seqtab.nochim.filt[rowSums(seqtab.nochim.filt) >= 1000,] # trim samples

taxa.species.names <- paste(taxa.eHOMD.filt[,1],taxa.eHOMD.filt[,2],taxa.eHOMD.filt[,3],
                          taxa.eHOMD.filt[,4],taxa.eHOMD.filt[,5],taxa.eHOMD.filt[,6],
                          taxa.eHOMD.filt[,7], sep="_")
threshold <- 1*(10^-4)
taxa.species <- filter.taxa(taxa.eHOMD.filt, taxa.species.names, seqtab.nochim.filt, threshold)
taxa.species <- taxa.species[rowSums(taxa.species) >= 1000,] # trim samples
taxa.species$other <- NULL
bools <- !(colnames(taxa.species) %like%  c("Bacteria_NA"))
taxa.species <- taxa.species[colnames(taxa.species)[bools]]
taxa.eHOMD.filt <- as.matrix.data.frame(taxa.eHOMD.filt)
rownames(taxa.eHOMD.filt) <- rownames(taxa.eHOMD)
colnames(taxa.eHOMD.filt) <- colnames(taxa.eHOMD)
taxa.eHOMD.filt <- as.data.frame(taxa.eHOMD.filt)
taxa.eHOMD.filt$species.merged <- taxa.species.names

l <- c()
for (i in rownames(taxa.species)){l <- c(l, unlist(strsplit(i, ".", fixed = TRUE))[1])}
rownames(seqtab.nochim.filt) <- l
rownames(taxa.species) <- l
preffix.sub <- paste0(path.out, "/", preffix, "_", sep = "")
write.csv(taxa.species, file = paste0(preffix.sub, "species_eHOMD.csv"))
#write.csv(seqtab.nochim.filt, file= paste0(preffix.sub, "sequencing_eHOMD.csv"))
write.csv(taxa.eHOMD.filt, file = paste0(preffix.sub, "taxas_eHMOD.csv"))
print("-- DADA2 pipelie and Assign Taxonomy finished --")

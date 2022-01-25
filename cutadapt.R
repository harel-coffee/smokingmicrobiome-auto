#############################################################################
#############################################################################
###########                                                       ###########
###########       Removes/trimms adapters from sequencing data    ###########
###########                    cutadapt v.1.15                    ###########
###########             Author(s): Liese Boonstra,                ###########
###########                        Diego Montiel Gonzalez         ###########
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
pkgTest("ShortRead")
pkgTest("Biostrings")
pkgTest("optparse")


# Build function that creates all orientations of a primer
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert DNAString object back to character vector
}

# Counts number of times primers appear in the reads, considers all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="PATH sequencing data", metavar="character"),
  make_option(c("-o", "--option"), type="character", default=NULL, 
              help="paired-end or single-end e.g. --option [PE-SE] (SE=PRJNA434300-PRJNA434312, PE=PRJNA484874)", 
              metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
path <- opt$input
option <- opt$option

# Select the path to cutadapt
cutadapt <- "cutadapt"

if (option == "SE"){
  #Primers
  # PRJNA434300-PRJNA434312 Single-end
  FWD <- "GGAGGCAGCAGTRRGGAAT" 
  REV <- "CTACCRGGGTATCTAATCC" 
  fnFs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
  # Create all orientations of the forward primer
  FWD.orients <- allOrients(FWD)
  # Create all orientations of the reverse primer
  REV.orients <- allOrients(REV)
  # Check if there are primers in the reads
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]))
  
  # Check if cutadapt is installed correctly
  path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  # Make list with names for cutadapt reads files in cutadapt/ subdirectory
  fnFs.cut <- file.path(path.cut, basename(fnFs))
  # Make reverse complementary sequence of forward primer
  FWD.RC <- dada2:::rc(FWD)
  # Make reverse complementary sequence of reverse primer
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of reads
  R1.flags <- paste("-g", FWD, "-a", REV.RC)
  
  # Run Cutadapt to remove primers
  for(i in seq_along(fnFs.cut)) {
    system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], # output files
                               fnFs[i], "--minimum-length", 100)) # input files
  }
  
  # Check if primers are removed
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))
  
}else if(option == "PE"){
  # PRJNA484874 Paired-end
  FWD = "GTGYCAGCMGCCGCGGTA"
  REV = "GGACTACHVGGGTWTCTAAT"
  fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))
  # Create all orientations of the forward primer
  FWD.orients <- allOrients(FWD)
  # Create all orientations of the reverse primer
  REV.orients <- allOrients(REV)
  # Make function that counts the number of times the primers appear in the reads,
  # while considering all possible primer orientations.
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  # Check if there are primers in the reads
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
  
  path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  # Make lists with names for cutadapt reads files in cutadapt/ subdirectory for both forward and reverse reads files
  fnFs.cut <- file.path(path.cut, basename(fnFs))
  fnRs.cut <- file.path(path.cut, basename(fnRs))
  # Make reverse complementary sequence of forward primer
  FWD.RC <- dada2:::rc(FWD)
  # Make reverse complementary sequence of reverse primer
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of forward reads
  R1.flags <- paste("-g", FWD, "-a", REV.RC)
  # Trim REV and the reverse-complement of FWD off of reverse reads
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run Cutadapt to remove primers
  for(i in seq_along(fnFs.cut)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs[i], fnRs[i], "--minimum-length", 100)) # input files
  }
  # Check if primers are removed
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
}else{
  print("-- Please check the README.md file --")
}


#############################################################################
#############################################################################
###########                                                       ###########
###########       Amplicon Sequence variants Species Consensus    ###########
###########             Author(s): Diego Montiel Gonzalez         ###########
###########                                                       ###########
###########        Erasmus MC University Medical Centre           ###########
###########               Rotterdam, The Netherlands              ###########
###########                                                       ###########
###########                     genid@erasmusmc.nl                ###########
###########                                                       ###########
#############################################################################
#############################################################################

### setwd("genid/smoking-microbiome/") # set project as working directory
#############################
### input data
#############################
filename <- "data/species_intersect.csv"
PRJNA434300.PRJNA434312.taxas <- "data/PRJNA434300-PRJNA434312_taxas_eHMOD.csv"
PRJNA484874.taxas <- "data/PRJNA484874_taxas_eHMOD.csv"
output.dir <- "data/"
setwd(".")

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
pkgTest("data.table")
pkgTest("Biostrings")
pkgTest("doParallel")
pkgTest("DECIPHER")
pkgTest("foreach")

# functions
get.consensus <- function(seqs){
  seqs <- DNAStringSet(seqs)
  #seqs <- OrientNucleotides(seqs)
  seqs.aln <- AlignSeqs(seqs)
  seqs.cons <- ConsensusSequence(seqs.aln, ignoreNonBases=TRUE, threshold=0.03)
  DNAStringSet(seqs.cons)
  return(as.character(seqs.cons))
}

process.df <- function(df){
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df$V1 <- NULL
  return(df)
}


df.taxa <- fread(filename, header = T)
df.taxa <- process.df(df.taxa)
s1.taxas <- fread(PRJNA434300.PRJNA434312.taxas, header = T)
s5.taxas <- fread(PRJNA484874.taxas, header = T)
taxas <- rbind(s1.taxas, s5.taxas)
taxas <- process.df(taxas)
taxas.OTU <- taxas

if(grepl("species", filename)){
  clusters <- as.data.frame(cbind(taxas$species.merged, rownames(taxas)))
  taxas.OTU$genus.merged <- NULL
  taxas.OTU <- taxas.OTU[taxas.OTU$species.merged %in% colnames(df.taxa), ]
  table(duplicated(taxas.OTU))
  taxas.OTU = taxas.OTU[!duplicated(taxas.OTU$species.merged),]
  print(table(taxas.OTU$species.merged == colnames(df.taxa)))
  taxas.OTU$species.merged <- NULL
}else{
  print("Error!, something failed check the format of the taxas files")
  stop()
}


colnames(clusters) <- c("cluster", "sequence")
length(intersect(colnames(df.taxa), clusters$cluster))
taxas <- intersect(colnames(df.taxa), as.vector(clusters$cluster))

numCores <- detectCores()
registerDoParallel(numCores)  # use multicore, set to the number of our cores
list.seqs.cons <- foreach(i = 1:length(taxas)) %dopar% {
  seqs <- clusters[clusters$cluster == as.vector(taxas[i]),]$sequence
  if (length(seqs) > 1){
    return(get.consensus(seqs))
  }else if(length(seqs) == 1){
    return(as.vector(seqs))
  }
}

list.seqs.cons <- matrix(unlist(list.seqs.cons), ncol = 1, byrow = TRUE)
rownames(list.seqs.cons) <- taxas
colnames(list.seqs.cons) <- "sequence"
table(colnames(df.taxa) == rownames(list.seqs.cons))

rownames(taxas.OTU) <- as.vector(list.seqs.cons[,1])

if(grepl("species", filename)){
  write.csv(list.seqs.cons, file=paste0(output.dir,"Species_OTU_seqs_consensus.csv"), quote = F)
  write.csv(taxas.OTU, file=paste0(output.dir, "Species_OTU_diversity.csv"), quote = F)
}

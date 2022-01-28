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
pkgTest("optparse")

#############################
### input data
#############################

option_list = list(
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="ASV at the Species-level, output from DADA2", 
              metavar="character"),
  make_option(c("-t", "--taxas"), type="character", default=NULL, 
              help="Input file of assigned taxas from DADA2 pipeline, at the
              Species-level", 
              metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output folder", 
              metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

file.asv <- "data/species_intersect.csv" # input ASV
file.taxas <- "data/taxas_eHMOD.csv" # taxas input from DADA2
output.dir <- "data/" # output folder
setwd(".")
file.asv <- opt$species
file.taxas <- opt$taxas
path.out <- opt$output


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


df.taxa <- fread(file.asv, header = T)
df.taxa <- process.df(df.taxa)
taxas <- fread(file.taxas, header = T)
taxas <- process.df(taxas)
taxas.species.merged <- paste(taxas[,1],taxas[,2],taxas[,3],
                              taxas[,4],taxas[,5],taxas[,6],
                              taxas[,7], sep="_")
taxas$species.merged <- taxas.species.merged
taxas.OTU <- taxas

#### how define for which taxa to be used

clusters <- as.data.frame(cbind(taxas$species.merged, rownames(taxas)))
taxas.OTU$genus.merged <- NULL
taxas.OTU <- taxas.OTU[taxas.OTU$species.merged %in% colnames(df.taxa), ]
table(duplicated(taxas.OTU))
taxas.OTU = taxas.OTU[!duplicated(taxas.OTU$species.merged),]
print(table(taxas.OTU$species.merged == colnames(df.taxa)))
taxas.OTU$species.merged <- NULL


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
write.csv(list.seqs.cons, file=paste0(output.dir,"Species_OTU_seqs_consensus.csv"), quote = F)


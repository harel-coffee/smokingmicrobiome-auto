#############################################################################
#############################################################################
###########                                                       ###########
###########       Batch effect assessment with guided PCA (gPCA)  ###########
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
##########################
## input data
##########################
PRJNA434300.PRJNA434312.species <-"data/PRJNA434300-PRJNA434312_species_eHOMD.csv"
PRJNA484874.species <-"data/PRJNA484874_species_eHOMD.csv"
metadata <- "data/metadata.csv"
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
pkgTest("gPCA")
pkgTest("data.table")
pkgTest("DataCombine")
pkgTest("ggplot2")


df.metadata <- read.csv(metadata)
df.metadata <- subset(df.metadata, smoking_status == "current" | smoking_status == "former"  | smoking_status == "never")
# Create replacements data frame
Replaces <- data.frame(from = c("1A", "1B","5"), to = c("S1", "S1", "S2"))
# Replace patterns and return the Var as a vector
df.metadata$study <- FindReplace(data = df.metadata, Var = "study", replaceData = Replaces,
                           from = "from", to = "to", vector = TRUE)

df.study.1 <- as.data.frame(fread(PRJNA434300.PRJNA434312.species, sep = ","))
rownames(df.study.1) <- df.study.1[,1]
df.study.1$V1 <- NULL
df.study.5 <- as.data.frame(fread(PRJNA484874.species, sep = ","))
rownames(df.study.5) <- df.study.5[,1]
df.study.5$V1 <- NULL
# intersect
taxas <- intersect(colnames(df.study.1), colnames(df.study.5))
df.taxa <- rbind(df.study.1[,taxas], df.study.5[,taxas])

samples <- intersect(df.metadata$SRR, rownames(df.taxa))
df.metadata <- subset(df.metadata, SRR %in% samples)
df.taxa <- df.taxa[samples, ]
df.taxa <- df.taxa[, colnames(df.taxa) != "Bacteria_NA_NA_NA_NA_NA_NA"]

###### Batch effect assessment with gPCA
batch <- as.numeric(as.factor((as.vector(df.metadata$study))))
table(batch)
#batch <- as.numeric(as.factor((as.vector(df.metadata$smoking_status))))

df.rel <- (df.taxa) / colMeans(df.taxa) # relative abundance
res.pca <- prcomp(df.rel, scale = T, center = T)
#res.pca <- prcomp(df.taxa, scale = T)
df_out <- as.data.frame(res.pca$x)
df_out$group <- sapply( strsplit(as.character(df.metadata$study), "_"), "[[", 1 )
head(df_out)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
percentage <- round(res.pca$sdev / sum(res.pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])
p

out <- gPCA.batchdetect(x = df.rel, batch = batch, center = F,  nperm = 1000)
out$delta ; out$p.val
gDist(out)
CumulativeVarPlot(out,ug="unguided",col="blue")
PCplot(out,ug="unguided",type="1v2")
PCplot(out,ug="unguided",type="comp",npcs=3)

write.csv(df.taxa, file = paste0(output.dir, "species_intersect.csv"))
write.csv(df_out, file = paste0(output.dir, "batch_effct_pca.csv"))

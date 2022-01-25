#############################################################################
#############################################################################
###########                                                       ###########
###########       Kruskal Wallist test for data type augmented    ###########
###########       comparison from results of NestCV in            ###########
###########       smoking microbiome sequencing data              ###########
###########             Author(s): Diego Montiel Gonzalez         ###########
###########                                                       ###########
###########        Erasmus MC University Medical Centre           ###########
###########               Rotterdam, The Netherlands              ###########
###########                                                       ###########
###########                     genid@erasmusmc.nl                ###########
###########                                                       ###########
#############################################################################
#############################################################################


setwd("/media/disk1/diego/genid/smoking-microbiome/data/results/")
df <- read.csv("NestedCV_results_feature_selection.csv", sep = ",")
               
               
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
pkgTest("tidyverse")
pkgTest("rstatix")
pkgTest("ggpubr")
pkgTest("data.table")


head(df)
df <- df[df$model == "SVML",]
hist(df[df$type  == "DEFAULT",]$mcc)
hist(df[df$type  == "SMOTE_over",]$mcc)
hist(df[df$type  == "SMOTE_both",]$mcc)
hist(df[df$type  == "TADA",]$mcc)

## Non-normal distribution
res.kruskal <- df %>% kruskal_test(mcc ~ type)
res.kruskal # if significant there some pairwise comparison differneces
df %>% kruskal_effsize(mcc ~ type)
#pwc <- df %>% dunn_test(mcc ~ type, p.adjust.method = "BH")
pwc <- df %>% wilcox_test(mcc ~ type, p.adjust.method = "BH")
pwc <- pwc %>% add_xy_position(x = "type")
ggboxplot(df, x = "type", y = "mcc", add = "jitter", color = "type", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07",  "#FC4E87", "#D68910", "#D68920" )) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.kruskal, detailed = TRUE), caption = get_pwc_label(pwc))
tmp <- as.data.frame(pwc)
fwrite(x = tmp, file="SVML_stats_mcc.csv")

##################################################################
res.kruskal <- df %>% kruskal_test(auc ~ type)
res.kruskal
df %>% kruskal_effsize(auc ~ type)
pwc <- df %>% wilcox_test(auc ~ type, p.adjust.method = "BH")
pwc <- pwc %>% add_xy_position(x = "type")
ggboxplot(df, x = "type", y = "auc", add = "jitter", color = "type", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07",  "#FC4E87", "#D68910", "#D68920" )) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.kruskal, detailed = TRUE), caption = get_pwc_label(pwc))
tmp <- as.data.frame(pwc)
fwrite(x = tmp, file="SVML_stats_auc.csv")

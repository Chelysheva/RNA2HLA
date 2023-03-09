#Title:
#RNA2HLA - checking the identity of RNA-seq samples based on HLA types
#Synopsys:
#Script is taking the final output of the RNA2HLA - all_samples_comparison_matrix.csv and plots a heatmap of the % of HLA similarity between each pair of samples
#Release: 1.1
#
#Author: Irina Chelysheva, 2019-2023
#Oxford Vaccine Group, University of Oxford
#Contact: irina.chelysheva@paediatrics.ox.ac.uk
#

all_samples_comparison_matrix <- read.delim("all_samples_comparison_matrix.csv", row.names=1, stringsAsFactors = F)

for (j in 1:dim(all_samples_comparison_matrix)){
  for (i in 1:dim(all_samples_comparison_matrix)){
    
    x=all_samples_comparison_matrix[i,j]
    x1=gsub("[()]", "",x)
    y=sub('\\,.*', '', x1)
    y=substring(y, 1)
    y<-as.character(y)
    all_samples_comparison_matrix[i,j]<-as.numeric(y)
  }
}

all_samples_comparison_matrix<-sapply(all_samples_comparison_matrix, as.numeric)
rownames(all_samples_comparison_matrix)<-colnames(all_samples_comparison_matrix)
all_samples_comparison_matrix<-data.matrix(all_samples_comparison_matrix,rownames.force = T)

#install.packages("gplots")
library(gplots)
#install.packages("RColorBrewer")
library(RColorBrewer)
pdf(file="all_samples_comparison_matrix_heatmap.pdf")
par(mar = rep(2, 4))
heatmap.2(all_samples_comparison_matrix,dendrogram='none', Rowv=F, Colv="Rowv", trace='none',col=brewer.pal(9,"Blues"), 
          main = "", density.info='none', key.xlab = "% HLA identity",cexRow = 0.5, cexCol = 0.5)

dev.off()

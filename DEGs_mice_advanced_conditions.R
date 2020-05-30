library(DESeq2)
library(stringr)
library(ggplot2)

setwd("~/Bioinfo/M2/Stage/breast/DATA/STAGE - DATA/count")

counts=read.csv("featureCounts.txt", sep="", head=T, skip=1, row.names = "Geneid")

colnames(counts)[12:21]  #ATTENTION!!! le premier élement n'est pas le a1 mais le a10!

conditions = c('mutant tumor 1', rep('control', 3), rep('mutant healthy', 3), 'mutant tumor 2', 'mutant tumor 3', 'mutant tumor 4')
mice = c('Advanced_1', paste('Control_',1:3, sep=''), paste('Early_',1:3, sep=''), paste('Advanced_',2:4, sep=''))

samples=cbind(mice, conditions);samples
rownames(samples)=samples[,1]
samples=as.data.frame(samples[,-1])
colnames(samples)="condition"

#Création du DESEQdataSet.
dds = DESeqDataSetFromMatrix(countData = counts[,6:15],colData = samples,design = ~ condition)

#PREFILTERING
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#DATA VISUALISATION WITH PCA
plotPCA(rlog(dds), intgroup="condition") +
  geom_point(size = 5) +
  geom_text(check_overlap = FALSE, label=mice ,vjust=0.5, colour = "#000000", size = 3.5) + 
  xlim(-50,50) + ylim(-25,40) 

######## Lehmann Genes #######
#Pour observer les gènes de la liste de Lehmann
setwd("~/Bioinfo/M2/Stage/breast/Final_all")
LehmannGenes = read.csv("genesLehmann_Mouse.csv")
LehmannGenes = as.vector(LehmannGenes[,2])

setwd("~/Bioinfo/M2/Stage/breast/DATA/STAGE - DATA/count")

LehmannGenes_final = NULL

for (i in 1:length(LehmannGenes)){
  if(LehmannGenes[i] %in% row.names(dds) == FALSE){
    print(LehmannGenes[i])
  }
  else{
    LehmannGenes_final <- c(LehmannGenes_final, LehmannGenes[i])
  }
}
#vst ou rlog, avec control 3 ou sans control 3
#DATA VISUALISATION WITH PCA
plotPCA(rlog(dds[LehmannGenes_final,]), intgroup="condition") +
  geom_point(size = 5) +
  geom_text(check_overlap = FALSE, label=mice ,vjust=0.5, colour = "#000000", size = 3.5) 

#DATA VISUALISATION WITH PCA
plotPCA(rlog(dds[,-4]), intgroup="condition") +
  geom_point(size = 5) +
  geom_text(check_overlap = FALSE, label=mice,vjust=0.5, colour = "#000000", size = 3.5) 


# dds <- estimateSizeFactors(dds)
# sizeFactors(dds)
# my_counts <- counts(dds, normalized=F)
# vst_my_counts = vst(my_counts)
# write.table(vst_my_counts, file="Mice_vst_my_counts.csv", sep="\t")

############################################################
################## Differential expression  ################
############################################################

library(cowplot)
library(mygene)
dds <- dds[,-4] #Sans le control 3
dds <- DESeq(dds)

resLFC1 <- lfcShrink(dds, coef="condition_mutant.tumor.1_vs_control", type="apeglm")
resLFC2 <- lfcShrink(dds, coef="condition_mutant.tumor.2_vs_control", type="apeglm")
resLFC3 <- lfcShrink(dds, coef="condition_mutant.tumor.3_vs_control", type="apeglm")
resLFC4 <- lfcShrink(dds, coef="condition_mutant.tumor.4_vs_control", type="apeglm")

# OBSERVATION TNBC RECEPTOR #

resLFC1[c(which(rownames(resLFC) == "Esr1"),
         which(rownames(resLFC) == "Pgr"),
         which(rownames(resLFC) == "Erbb2")),]

resLFC2[c(which(rownames(resLFC) == "Esr1"),
         which(rownames(resLFC) == "Pgr"),
         which(rownames(resLFC) == "Erbb2")),]

resLFC3[c(which(rownames(resLFC) == "Esr1"),
         which(rownames(resLFC) == "Pgr"),
         which(rownames(resLFC) == "Erbb2")),]

resLFC4[c(which(rownames(resLFC) == "Esr1"),
         which(rownames(resLFC) == "Pgr"),
         which(rownames(resLFC) == "Erbb2")),]


######################################################################################################
#JOLIE BARRRE DE PROGRESSION

progress <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}
# RECUPERATION DES GENE UP REGULE
filterUP <- function(res){
  gene_to_keep = NULL
  pval = NULL
  FClog2 = NULL
  for(i in 1:nrow(as.matrix(res))){
    if(res$padj[i] <= 0.05 && res$log2FoldChange[i] >= 0 && is.na(res$padj[i]) != T && is.na(res$log2FoldChange[i]) != T)
    {
      gene_to_keep = c(gene_to_keep, rownames(res)[i])
      pval = c(pval, res$padj[i])
      FClog2 = c(FClog2, res$log2FoldChange[i])
    }
    progress(i, max = nrow(as.matrix(res)))
  }
  return(matrix(c(gene_to_keep, pval, FClog2), ncol = 3))
}
# RECUPERATION DES GENE DOWN REGULE
filterDOWN <- function(res){
  gene_to_keep = NULL
  pval = NULL
  FClog2 = NULL
  for(i in 1:nrow(as.matrix(res))){
    if(res$padj[i] <= 0.05 && res$log2FoldChange[i] <= 0 && is.na(res$padj[i]) != T && is.na(res$log2FoldChange[i]) != T)
    {
      gene_to_keep = c(gene_to_keep, rownames(res)[i])
      pval = c(pval, res$padj[i])
      FClog2 = c(FClog2, res$log2FoldChange[i])
    }
    progress(i, max = nrow(as.matrix(res)))
  }
  return(matrix(c(gene_to_keep, pval, FClog2), ncol = 3))
}

#J'aurais pu tout faire dans une fonction, et meme faire le tri de - et des plus sur pythons mais différencié les 2 listes est plus
#simple par la suite pour l'observation des genes par les biologistes.

######### TUMOR 1 VS CONTROL ######### 

UP_genes_control_vs_tumor1 <- filterUP(resLFC1)
DOWN_genes_control_vs_tumor1 <- filterDOWN(resLFC1)

nrow(UP_genes_control_vs_tumor1)
nrow(DOWN_genes_control_vs_tumor1)

UP_genes_control_vs_tumor1 <- UP_genes_control_vs_tumor1[order(as.numeric(UP_genes_control_vs_tumor1[,3]), decreasing = TRUE),] #Tri dans l'ordre croissant des pvalue ajusté
DOWN_genes_control_vs_tumor1 <- DOWN_genes_control_vs_tumor1[order(as.numeric(DOWN_genes_control_vs_tumor1[,3]), decreasing = FALSE),] #Tri dans l'ordre croissant des pvalue ajusté

queryUP <- queryMany(UP_genes_control_vs_tumor1[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')
queryDOWN <- queryMany(DOWN_genes_control_vs_tumor1[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')

functionGenesUP = list()
functionGenesDOWN = list()

for(i in 1:nrow(UP_genes_control_vs_tumor1)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesUP[[i]] <- paste(queryUP[i, 'go.BP'][[1]]$term, collapse = " / ")
}

for(i in 1:nrow(DOWN_genes_control_vs_tumor1)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesDOWN[[i]] <- paste(queryDOWN[i, 'go.BP'][[1]]$term, collapse = " / ")
}

tableauExportUP = cbind(UP_genes_control_vs_tumor1, functionGenesUP)
tableauExportDOWN = cbind(DOWN_genes_control_vs_tumor1, functionGenesDOWN)

colnames(tableauExportUP) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")
colnames(tableauExportDOWN) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")

write.table(tableauExportUP, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/UP_genes_control_vs_tumor1.csv", sep= ",")
write.table(tableauExportDOWN, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/DOWN_genes_control_vs_tumor1.csv", sep = ",")
#Pour chaque sample on récupere les genes + ou - différencié et on leur rajoute une information biologique pour l'analyse par les biologistes

######### TUMOR 2 VS CONTROL ######### 

UP_genes_control_vs_tumor2 <- filterUP(resLFC2)
DOWN_genes_control_vs_tumor2 <- filterDOWN(resLFC2)

nrow(UP_genes_control_vs_tumor2)
nrow(DOWN_genes_control_vs_tumor2)

UP_genes_control_vs_tumor2 <- UP_genes_control_vs_tumor2[order(as.numeric(UP_genes_control_vs_tumor2[,3]), decreasing = TRUE),] #Tri dans l'ordre croissant des pvalue ajusté
DOWN_genes_control_vs_tumor2 <- DOWN_genes_control_vs_tumor2[order(as.numeric(DOWN_genes_control_vs_tumor2[,3]), decreasing = FALSE),] #Tri dans l'ordre croissant des pvalue ajusté

queryUP <- queryMany(UP_genes_control_vs_tumor2[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')
queryDOWN <- queryMany(DOWN_genes_control_vs_tumor2[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')

functionGenesUP = list()
functionGenesDOWN = list()

for(i in 1:nrow(UP_genes_control_vs_tumor2)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesUP[[i]] <- paste(queryUP[i, 'go.BP'][[1]]$term, collapse = " / ")
}

for(i in 1:nrow(DOWN_genes_control_vs_tumor2)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesDOWN[[i]] <- paste(queryDOWN[i, 'go.BP'][[1]]$term, collapse = " / ")
}

tableauExportUP = cbind(UP_genes_control_vs_tumor2, functionGenesUP)
tableauExportDOWN = cbind(DOWN_genes_control_vs_tumor2, functionGenesDOWN)

colnames(tableauExportUP) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")
colnames(tableauExportDOWN) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")

write.table(tableauExportUP, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/UP_genes_control_vs_tumor2.csv", sep= ",")
write.table(tableauExportDOWN, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/DOWN_genes_control_vs_tumor2.csv", sep = ",")

######### TUMOR 3 VS CONTROL ######### 

UP_genes_control_vs_tumor3 <- filterUP(resLFC3)
DOWN_genes_control_vs_tumor3 <- filterDOWN(resLFC3)

nrow(UP_genes_control_vs_tumor3)
nrow(DOWN_genes_control_vs_tumor3)

UP_genes_control_vs_tumor3 <- UP_genes_control_vs_tumor3[order(as.numeric(UP_genes_control_vs_tumor3[,3]), decreasing = TRUE),] #Tri dans l'ordre croissant des pvalue ajusté
DOWN_genes_control_vs_tumor3 <- DOWN_genes_control_vs_tumor3[order(as.numeric(DOWN_genes_control_vs_tumor3[,3]), decreasing = FALSE),] #Tri dans l'ordre croissant des pvalue ajusté

queryUP <- queryMany(UP_genes_control_vs_tumor3[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')
queryDOWN <- queryMany(DOWN_genes_control_vs_tumor3[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')

functionGenesUP = list()
functionGenesDOWN = list()

for(i in 1:nrow(UP_genes_control_vs_tumor3)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesUP[[i]] <- paste(queryUP[i, 'go.BP'][[1]]$term, collapse = " / ")
}

for(i in 1:nrow(DOWN_genes_control_vs_tumor3)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesDOWN[[i]] <- paste(queryDOWN[i, 'go.BP'][[1]]$term, collapse = " / ")
}

tableauExportUP = cbind(UP_genes_control_vs_tumor3, functionGenesUP)
tableauExportDOWN = cbind(DOWN_genes_control_vs_tumor3, functionGenesDOWN)

colnames(tableauExportUP) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")
colnames(tableauExportDOWN) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")

write.table(tableauExportUP, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/UP_genes_control_vs_tumor3.csv", sep= ",")
write.table(tableauExportDOWN, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/DOWN_genes_control_vs_tumor3.csv", sep = ",")



######### TUMOR 4 VS CONTROL ######### 

UP_genes_control_vs_tumor4 <- filterUP(resLFC4)
DOWN_genes_control_vs_tumor4 <- filterDOWN(resLFC4)

nrow(UP_genes_control_vs_tumor4)
nrow(DOWN_genes_control_vs_tumor4)

UP_genes_control_vs_tumor4 <- UP_genes_control_vs_tumor4[order(as.numeric(UP_genes_control_vs_tumor4[,3]), decreasing = TRUE),] #Tri dans l'ordre croissant des pvalue ajusté
DOWN_genes_control_vs_tumor4 <- DOWN_genes_control_vs_tumor4[order(as.numeric(DOWN_genes_control_vs_tumor4[,3]), decreasing = FALSE),] #Tri dans l'ordre croissant des pvalue ajusté

queryUP <- queryMany(UP_genes_control_vs_tumor4[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')
queryDOWN <- queryMany(DOWN_genes_control_vs_tumor4[,1], scopes='symbol', fields=c('entrezgene', 'go'), species='mouse')

functionGenesUP = list()
functionGenesDOWN = list()

for(i in 1:nrow(UP_genes_control_vs_tumor4)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesUP[[i]] <- paste(queryUP[i, 'go.BP'][[1]]$term, collapse = " / ")
}

for(i in 1:nrow(DOWN_genes_control_vs_tumor4)){
  # cat("\n###########", genes_control_vs_tumor[i], '###########\n')
  functionGenesDOWN[[i]] <- paste(queryDOWN[i, 'go.BP'][[1]]$term, collapse = " / ")
}

tableauExportUP = cbind(UP_genes_control_vs_tumor4, functionGenesUP)
tableauExportDOWN = cbind(DOWN_genes_control_vs_tumor4, functionGenesDOWN)

colnames(tableauExportUP) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")
colnames(tableauExportDOWN) <- c("index\tgene name", "pvalue", "log2FC", "fonctions")

write.table(tableauExportUP, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/UP_genes_control_vs_tumor4.csv", sep= ",")
write.table(tableauExportDOWN, file = "~/Bioinfo/M2/Stage/breast/DESEQ2/Sans3/DOWN_genes_control_vs_tumor4.csv", sep = ",")



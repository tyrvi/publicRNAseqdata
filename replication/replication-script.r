# Replication of results from [here](https://github.com/hussius/publicRNAseqdata)
  
library(pheatmap)
library(reshape)
library(gplots)
library(ops)
library(calibrate)
library(biomaRt)
library(sva)
library(ggplot2)

# Implement Principal Component Analysis (PCA)
plot.pca.published <- function(df,x,y,z){
  
  p <- prcomp(t(df))
  
  v <- 100*(p$sdev)^2 / sum(p$sdev^2)
  v.x <- v[x]
  v.y <- v[y]
  
  colors <- c("indianred", "dodgerblue", "forestgreen",
              "indianred", "dodgerblue",
              "indianred", "dodgerblue", "forestgreen", 
              "indianred", "dodgerblue", "forestgreen")
  
  shapes <- c(rep(15,3),rep(16,2),rep(17,3),rep(8,3))
  
  plot(p$x[,x],p$x[,y],pch=shapes,cex=1.5,col=colors,xlab=paste(paste("PC",x),round(v.x),"% of variance"),ylab=paste(paste("PC",y),round(v.y),"% of variance"),main=paste(z," Published FPKM/RPKM values \n n=13,078"))
  
}

plot.pca.reprocessed <- function(df,x,y,z){
  
  p_repr <- prcomp(t(df))
  
  v <- 100*(p_repr$sdev)^2 / sum(p_repr$sdev^2)
  v.x <- v[x]
  v.y <- v[y]
  
  colors_repr <- c("dodgerblue", "indianred", "forestgreen",
                   "dodgerblue", "indianred", "forestgreen",
                   "dodgerblue", "indianred", "forestgreen",
                   "dodgerblue", "indianred", "forestgreen",
                   "dodgerblue", "indianred")          
  
  shapes_repr <- c(rep(11,3),rep(8,3),rep(17,3),rep(15,3),rep(16,2))
  
  plot(p_repr$x[,x],p_repr$x[,y],pch=shapes_repr,col=colors_repr,xlab=paste(paste("PC",x),round(v.x),"% of variance"),ylab=paste(paste("PC",y),round(v.y),"% of variance"),main=paste(z," reprocessed  cufflinks FPKM values \n n=18,175"))
  
}

# Causes data to only include protein coding genes using ensembl
getPcs <- function(x) {
  ensembl = useMart("ensembl",dataset= "hsapiens_gene_ensembl")
  type <- getBM(attributes=c("ensembl_gene_id","gene_biotype"),filters = "ensembl_gene_id",values = as.vector(x$ENSEMBL_ID),mart = ensembl)
  pc <- subset(type[,1],type[,2]=="protein_coding")
  results <- x[match(pc,x$ENSEMBL_ID),]              
  return(results)
}

# Correlations between Principal Components (PCs) and experimental factors
print_PCA_corrs <- function(data,sampleinfo,caption="PCA correlations",include.quant=F){
  pca <- prcomp(t(data[,]))
  rot <- pca$r
  x <- pca$x
  
  if(include.quant){
    pc1 <- rep(1,8)
    names(pc1) <- c("Raw reads","Mapped reads", "Tissue", "Library prep", "Study", "Read type", "Read length","Quantification method")
    pc2 <- rep(1,8)
  }
  else{pc1 <- rep(1,7)
  names(pc1) <- c("Raw reads","Mapped reads", "Tissue", "Library prep", "Study", "Read type", "Read length")
  pc2 <- rep(1,7)
  }
  names(pc2) <- names(pc1)
  
  # Test correlations between number of seq'd reads and PCs 1-4 from prcomp
  pval.nraw.pc1 <- cor.test(x[,1], sampleinfo$NumberRaw,method="spearman")$p.value
  pval.nraw.pc2 <- cor.test(x[,2], sampleinfo$NumberRaw,method="spearman")$p.value
  pval.nraw.pc3 <- cor.test(x[,3], sampleinfo$NumberRaw,method="spearman")$p.value
  
  cat(sprintf("Number_of_rawreads~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\n", pval.nraw.pc1,pval.nraw.pc2,pval.nraw.pc3))
  
  pc1[1] <- pval.nraw.pc1
  pc2[1] <- pval.nraw.pc2
  
  pval.nmapped.pc1 <- cor.test(x[,1], sampleinfo$Numbermapped,method="spearman")$p.value
  pval.nmapped.pc2 <- cor.test(x[,2], sampleinfo$Numbermapped,method="spearman")$p.value
  pval.nmapped.pc3 <- cor.test(x[,3], sampleinfo$Numbermapped,method="spearman")$p.value
  
  cat(sprintf("Number_of_mappedreads~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\n", pval.nmapped.pc1,pval.nmapped.pc2,pval.nmapped.pc3))
  
  pc1[2] <- pval.nmapped.pc1
  pc2[2] <- pval.nmapped.pc2
  
  # For tissue, use kruskal.test which handles ordinal variables 
  pval.tissue.pc1<-kruskal.test(x[,1], sampleinfo$Tissue)$p.value
  pval.tissue.pc2<-kruskal.test(x[,2], sampleinfo$Tissue)$p.value
  pval.tissue.pc3<-kruskal.test(x[,3], sampleinfo$Tissue)$p.value
  
  cat(sprintf("Tissues~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.tissue.pc1,pval.tissue.pc2,pval.tissue.pc3))
  
  pc1[3] <- pval.tissue.pc1
  pc2[3] <- pval.tissue.pc2
  
  # Library prep 
  pval.prep.pc1<-kruskal.test(x[,1], sampleinfo$Preparation)$p.value
  pval.prep.pc2<-kruskal.test(x[,2], sampleinfo$Preparation)$p.value
  pval.prep.pc3<-kruskal.test(x[,3], sampleinfo$Preparation)$p.value
  
  cat(sprintf("LibPrep~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.prep.pc1,pval.prep.pc2,pval.prep.pc3))
  
  pc1[4] <- pval.prep.pc1
  pc2[4] <- pval.prep.pc2
  
  # Study  
  pval.study.pc1<-kruskal.test(x[,1], sampleinfo$Study)$p.value
  pval.study.pc2<-kruskal.test(x[,2], sampleinfo$Study)$p.value
  pval.study.pc3<-kruskal.test(x[,3], sampleinfo$Study)$p.value
  
  cat(sprintf("Study~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.study.pc1,pval.study.pc2,pval.study.pc3))
  
  pc1[5] <- pval.study.pc1
  pc2[5] <- pval.study.pc2
  
  # Layout
  pval.layout.pc1<-kruskal.test(x[,1], sampleinfo$Readtype)$p.value
  pval.layout.pc2<-kruskal.test(x[,2], sampleinfo$Readtype)$p.value
  pval.layout.pc3<-kruskal.test(x[,3], sampleinfo$Readtype)$p.value
  
  cat(sprintf("Study~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.layout.pc1,pval.layout.pc2,pval.layout.pc3))
  
  pc1[6] <- pval.layout.pc1
  pc2[6] <- pval.layout.pc2
  
  # Read length
  pval.readlength.pc1<-cor.test(x[,1], sampleinfo$readlength)$p.value
  pval.readlength.pc2<-cor.test(x[,2], sampleinfo$readlength)$p.value
  pval.readlength.pc3<-cor.test(x[,3], sampleinfo$readlength)$p.value
  
  cat(sprintf("ReadType~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.layout.pc1,pval.layout.pc2,pval.layout.pc3))
  
  library("ggplot2")
  
  pc1[7] <- pval.readlength.pc1
  pc2[7] <- pval.readlength.pc2
  
  # Quantification
  if(include.quant){
    pval.quant.pc1<-kruskal.test(x[,1], as.factor(sampleinfo$quantification))$p.value
    pval.quant.pc2<-kruskal.test(x[,2], as.factor(sampleinfo$quantification))$p.value
    pval.quant.pc3<-kruskal.test(x[,3], as.factor(sampleinfo$quantification))$p.value
    
    cat(sprintf("QuantificationMethod~PCAs: PCA1=\"%f\"PCA2=\"%f\"PCA3=\"%f\"\n", pval.layout.pc1,pval.layout.pc2,pval.layout.pc3))
    
    pc1[8] <- pval.quant.pc1
    pc2[8] <- pval.quant.pc2
  }
  par(mfrow=c(2,1))
  barplot(-log(pc1),las=2,main=paste(caption, "PC1"),ylab="-log(p)")
  barplot(-log(pc2),las=2,main=paste(caption, "PC2"),ylab="-log(p)")
}

# ANOVA (Analysis of Variance) for estimating influence of factors
do_anova <- function(data, sampleinfo, caption="ANOVA", include.quant=F){
  m <- melt(data)
  colnames(m) <- c("sample_ID","RPKM")
  if (include.quant){
    meta <- sampleinfo[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype","readlength","quantification")]}
  else{
    meta <- sampleinfo[,c("Study","Tissue","Preparation","NumberRaw","Numbermapped","Readtype","readlength")]}
  rownames(meta) <- colnames(data)
  tissue <- rep(meta$Tissue, each=nrow(data))
  study <- rep(meta$Study, each=nrow(data))
  prep <- rep(meta$Preparation, each=nrow(data))
  layout <- rep(meta$Readtype, each=nrow(data))
  raw <- rep(meta$NumberRaw, each=nrow(data))
  mapped <- rep(meta$Numbermapped, each=nrow(data))
  readlen <- rep(meta$readlength, each=nrow(data))
  if(include.quant){
    quant <- rep(meta$quantification, each=nrow(data))
    matrix <- data.frame(m, tissue=tissue, study=study, prep=prep, layout=layout, readlen=readlen, nraw=raw,nmapped=mapped, quant=quant)
    fit <- lm(RPKM ~ layout + readlen + prep + nraw + quant + study + tissue, data=matrix)
  }
  else{
    matrix <- data.frame(m, tissue=tissue, study=study, prep=prep, layout=layout, readlen=readlen, nraw=raw,nmapped=mapped)
    fit <- lm(RPKM ~ layout + readlen + prep + nraw + study + tissue, data=matrix)
  }
  
  a <- anova(fit)
  nfac <- length(a[,1])-1
  maxval = 100
  barplot(100*a$"Sum Sq"[1:nfac]/sum(a$"Sum Sq"[1:nfac]),names.arg=rownames(a[1:nfac,]),main=caption,ylim=c(0,maxval))
}

# Plot highest gene loading
plot.highest.loadings <- function(p, fpkm.table, caption="Highest loading"){
  par(mfrow=c(3,1))
  for(i in 1:3){
    load <- p$rotation[,i][order(p$rotation[,i])]
    extreme <- c(tail(load), head(load))
    extreme.ensg <- names(extreme)
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") #select the ensembl database
    extreme.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                             filters = "ensembl_gene_id",
                             values=extreme.ensg,
                             mart=ensembl)
    q <- extreme.symbols[,2]
    names(q) <- extreme.symbols[,1]
    fpkm <- cbind(q[extreme.ensg],fpkm.table[extreme.ensg,])
    names(fpkm)[names(fpkm) == 'q[extreme.ensg]'] <- 'Gene Symbol'
    barplot(extreme, names.arg=q[extreme.ensg],las=2,main=paste0(caption, ", PC", i))
    print(fpkm)
  }
}

# Adjust data for distributions with ggplot
distr_data <- function(x,y,z,k) {
  
  if(k==1){
    if(y=="heart"){
      data <- data.frame(x[grep("heart", names(x), value = TRUE)])
    }
    else if(y=="brain"){
      data <- data.frame(x[grep("brain", names(x), value = TRUE)])
    }
    else{
      data <- data.frame(x[grep("kidney", names(x), value = TRUE)])
    }
    colnames(data) <- "FPKM"
    data$study <- z
    return(data)
  }
  else if(k==2){
    if(y=="heart"){
      tissue_data <- data.frame(x[grep("heart", names(x), value = TRUE)])
      data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
    }
    else if(y=="brain"){
      tissue_data <- data.frame(x[grep("brain", names(x), value = TRUE)])
      data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
    }
    else{
      tissue_data <- data.frame(x[grep("kidney", names(x), value = TRUE)])
      data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
    }
    colnames(data) <- "FPKM"
    data$study <- z
    return(data)
  }
}

# Adjust data for distributions with ggplot for reprocessed data
distr_data_repr <- function(x,y,z) {
  if(y=="heart"){
    tissue_data <- data.frame(x[grep("heart", names(x), value = TRUE)])
    data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
  }
  else if(y=="brain"){
    tissue_data <- data.frame(x[grep("brain", names(x), value = TRUE)])
    data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
  }
  else{
    tissue_data <- data.frame(x[grep("kidney", names(x), value = TRUE)])
    data <- data.frame(tissue_data[grep(z, names(tissue_data), value = TRUE)])
  }
  colnames(data) <- "FPKM"
  data$study <- paste(z,"_repr.",sep ="")
  return(data)
}

# Quantile Normalization
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


# More Libaries
library(org.Hs.eg.db) # for transferring gene identifiers
library(data.table) # for collapsing transcript RPKMs
library(pheatmap) # for nicer visualization

# Data analysis begins
# All of the data was already downloaded, combined and written to the file **published_rkpms.txt**
# The data appears as the following
#
# HPA_heart HPA_brain HPA_kidney AltIso_heart AltIso_brain GTEx_heart GTEx_brain GTEx_kidney Atlas_heart Atlas_brain Atlas_kidney
# ENSG00000000003 6.7 11.8 44.4 2.51 2.41 4.15468311309814 5.41635084152222 11.2788972854614 2.748 4.999 14.768

# Get data from local files
published <- read.delim("published_rpkms.txt",sep=" ")
sampleinfo_published <- read.table("sample_info_published.txt",header=TRUE)

# Figure 1
published.nozero <- published[-which(rowMeans(published[,])<=0.01),]

# Fig. 1B
pheatmap(cor(published.nozero, method="spearman"))

# Figure 1C
plot.pca.published(published.nozero,1,2,"")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20)
legend("top",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)

# Figure 1D
do_anova(published.nozero,sampleinfo_published,"ANOVA, published data")

# Figure 2
pseudo <- 1
published.log <- log2(published.nozero + pseudo)

# Figure 2A
pheatmap(cor(published.log),method="spearman")

# Figure 2B
plot.pca.published(published.log,1,2,"log2")
legend("left",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20)
legend("top",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)

# Figure 2C
plot.pca.published(published.log,2,3,"log2")
legend("bottom",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20)
legend("topright",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)

# Figure 2D
colors <- c("indianred","dodgerblue","forestgreen",
            "indianred","dodgerblue","indianred",
            "dodgerblue","forestgreen", "indianred", 
            "dodgerblue", "forestgreen")

p.loo <- published.log[,-c(4,5)]
colors.loo <- colors[-c(4,5)]
p.loo.log <- prcomp(t(p.loo))
p.add <- published.log[,c(4,5)]
projection <- t(p.add) %*% p.loo.log$rotation
p.original.plus.new <- rbind(p.loo.log$x, projection)
col.original.plus.new <- c(colors.loo, colors[c(4,5)])
plot(p.original.plus.new[,2],p.original.plus.new[,3],pch=c(rep(15,3),rep(17,3),rep(8,3),rep(22,nrow(projection))),col=col.original.plus.new,xlab="PC2",ylab="PC3",main="log2 Published FPKM/RPKM values; AltIso projected onto existing PCs \n n=13,323",xlim=c(-150,100))
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("HPA","GTEx","Atlas","AltIso"),col="black",pch=c(15,17,8,22),ncol=2)

# Figure 2E
do_anova(published.log,sampleinfo_published,"ANOVA, published data (log)")

# Figure 3
meta <- data.frame(study=c(rep("HPA",3),rep("AltIso",2),rep("GTex",3),rep("Atlas",3)),tissue=c("Heart","Brain","Kidney","Heart","Brain","Heart","Brain","Kidney","Heart","Brain","Kidney"))
batch <- meta$study
design <- model.matrix(~1,data=meta)
combat <- ComBat(dat=published.log,batch=batch,mod=design,par.prior=TRUE)
# Original code numCovs was removed since it gave an error of 
# Error in ComBat(dat = published.log, batch = batch, mod = design, numCovs = NULL,  : 
# unused argument (numCovs = NULL)
# When removing the numCovs argument the result is the same as in the paper
# combat <- ComBat(dat=published.log,batch=batch,mod=design,numCovs=NULL,par.prior=TRUE)

# Figure 3A
pheatmap(cor(combat, method="spearman")) 

# Figure 3B
plot.pca.published(combat,1,2,"COMBAT")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20)
legend("top",legend=c("HPA","AltIso","GTEx","Atlas"),col="black",pch=c(15,16,17,8),ncol=2)

# Figure 3C
do_anova(combat, sampleinfo_published, caption="ANOVA, ComBat adj log F/RPKM")

# Figure 4
cufflinks <- read.delim("fpkm_table_tophat.txt")
sampleinfo_cufflinks <- read.delim("sample_info_reprocessed.txt")
cufflinks_pc <- getPcs(cufflinks)
cufflinks_pc_nozero <- cufflinks_pc[-which(rowMeans(cufflinks_pc[,3:16])<=0.01),]

# Figure 4A
pheatmap(cor(cufflinks_pc_nozero[,3:16], method="spearman")) 

# Figure 4B
cufflinks_fpkms <- cufflinks_pc_nozero[,3:16]
rownames(cufflinks_fpkms) <- cufflinks_pc_nozero[,1]
plot.pca.reprocessed(cufflinks_fpkms,1,2,"")
legend("bottomright",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)

# Figure 4C
do_anova(cufflinks_fpkms,sampleinfo_cufflinks,caption="ANOVA, Cufflinks FPKM")

# Figure 4D
pseudo <- 1
cufflinks_log <- log2(cufflinks_fpkms + pseudo)
pheatmap(cor(cufflinks_log) ,method="spearman")

# Figure 4E
plot.pca.reprocessed(cufflinks_log,1,2,"log2")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)

# Figure 4F
plot.pca.reprocessed(cufflinks_log,2,3,"log2")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("topleft",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)

# Figure 4G Does not work
# p.log.cufflinks <- prcomp(t(cufflinks_log))
# plot.highest.loadings(p.log.cufflinks,cufflinks_log,caption="Cufflinks log FPKM")
# par(mfrow=c(4,4))
# for (i in 1:6){
#   for(j in 1:6){
#     if (i<j){ 
#   	plot(p.log.cufflinks$x[,i],p.log.cufflinks$x[,j],pch=shapes_cufflinks,col=colors,xlab=paste("PC",i),ylab=paste("PC",j),main="log2 reprocessed FPKM values \n n=19475")
# 		}
# 	}
# }
# colors <- c("dodgerblue", "indianred", "forestgreen",
#             "dodgerblue", "indianred", "forestgreen",
#             "dodgerblue", "indianred", "forestgreen",
#             "dodgerblue", "indianred", "forestgreen",
#             "dodgerblue", "indianred") 
# p.loo <- cufflinks_log[,-c(13,14)]
# colors.loo <- colors[-c(13,14)]
# p.loo.log <- prcomp(t(p.loo))
# p.add <- cufflinks_log[,c(13,14)]
# projection <- t(p.add) %*% p.loo.log$rotation
# p.original.plus.new <- rbind(p.loo.log$x, projection)
# col.original.plus.new <- c(colors.loo, colors[c(13,14)])
# plot(p.original.plus.new[,2],p.original.plus.new[,3],pch=c(shapes_cufflinks[1:12],rep(22,nrow(projection))),col=col.original.plus.new,xlab="PC2",ylab="PC3",main="log2 Cufflinks FPKM values; AltIso projected onto existing PCs \n n=19,475,",xlim=c(-150,100))
# legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
# legend("topleft",legend=c("HPA","GTEx","Atlas","AltIso"),col="black",pch=c(11,8,17,15,22),ncol=2)
```

# Figure 4H
meta <- data.frame(study=c(rep("EoGE",3),rep("Atlas",3),rep("BodyMap",3),rep("HPA",3),rep("AltIso",2)),tissue=c("Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart","Kidney","Brain","Heart"),prep=c(rep("poly-A",3),rep("rRNA-depl",3),rep("poly-A",8)),layout=c(rep("PE",3),rep("SE",3),rep("PE",6),rep("SE",2)))
batch <- meta$study
design <- model.matrix(~1,data=meta)
combat.cufflinks <- ComBat(dat=cufflinks_log,batch=batch,mod=design,par.prior=TRUE)
pheatmap(cor(combat.cufflinks),method="spearman")

# Figure 4I
plot.pca.reprocessed(combat.cufflinks,1,2,"COMBAT")
legend("bottomleft",legend=c("Heart","Brain","Kidney"),col=c("indianred", "dodgerblue", "forestgreen"),cex=1.5,pch=20,bty="n")
legend("top",legend=c("EoGE","Atlas","BodyMap","HPA","AltIso"),col="black",pch=c(11,8,17,15,16),ncol=2)

# Figure 4J
do_anova(combat.cufflinks,sampleinfo_cufflinks,caption="ANOVA, Cufflinks ComBat adj log FPKM")


# Figure 4K-L
j <- merge(published.nozero, cufflinks_pc_nozero, by.x=0, by.y="ENSEMBL_ID")
j <- j[,-which(colnames(j)=="Gene_ID"),]
rown <- j[,1]
j <- j[,2:ncol(j)]
rownames(j) <- rown

sampleinfo_cufflinks$quantification <- "topcuff"
sampleinfo_published$quantification <- c(rep("topcuff",3),rep("custom_AltIso",2),rep("fluxcapacitor",3),rep("custom_Atlas",3))
studylabels_cuff <- c("EoGE_brain","EoEG_heart","EoEG_kidney","Atlas_brain","Atlas_heart","Atlas_kidney","BodyMap_brain","BodyMap_heart","BodyMap_kidney","HPA_brain","HPA_heart","HPA_kidney","AltIso_brain","AltIso_heart")
sampleinfo_cufflinks <- data.frame(Study_labels=studylabels_cuff, sampleinfo_cufflinks)
sampleinfo <- rbind(sampleinfo_published, sampleinfo_cufflinks)
par(mfrow=c(3,1))
# figure 4K
do_anova(j, sampleinfo, caption="Joint analysis, raw F/RPKM", include.quant=T)
# figure 4L
do_anova(combat.j, sampleinfo, caption="Joint analysis, ComBat adj log F/RPKM", include.quant=T)
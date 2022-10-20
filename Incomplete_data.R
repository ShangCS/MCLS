library(survival)
library("MASS")
source("run_data_col.R")
source("specClust.R")
source("enrichment.R")
assignInNamespace('specClust',specClust,ns='kknn')
environment(specClust) <- asNamespace('kknn')


get_emb<-function(data,complete_names,all_names,dim,dim_pca,pca_scale=F){
  print("This time is for incomplete data (For clustering).",quote=FALSE)
  mRNA_pca<-prcomp(t(data[["mRNA"]]), center=pca_scale, scale=pca_scale)
  Methy_pca<-prcomp(t(data[["Methy"]]), center=pca_scale, scale=pca_scale)
  miRNA_pca<-prcomp(t(data[["miRNA"]]), center=pca_scale, scale=pca_scale)
  
  complete_mRNA_pca <- mRNA_pca$x[1:length(complete_names),1:dim_pca]
  complete_Methy_pca <- Methy_pca$x[1:length(complete_names),1:dim_pca]
  complete_miRNA_pca <- miRNA_pca$x[1:length(complete_names),1:dim_pca]
  
  #incomplete_mRNA_pca <- mRNA_pca$x[(length(complete_names)+1):length(all_names),1:dim_pca]
  incomplete_Methy_pca <- Methy_pca$x[(length(complete_names)+1):length(all_names),1:dim_pca]
  incomplete_miRNA_pca <- miRNA_pca$x[(length(complete_names)+1):length(all_names),1:dim_pca]

  completed_data_pca<-cbind(complete_mRNA_pca,complete_Methy_pca,complete_miRNA_pca)
  complete_data_svd<-svd(completed_data_pca)
  yc<-complete_data_svd[["u"]][,1:dim]
  
  #yw_data1 <- t(yc) %*% t(ginv(complete_mRNA_pca)) %*% t(incomplete_mRNA_pca)
  yw_data1 <- t(yc) %*% t(ginv(complete_Methy_pca)) %*% t(incomplete_Methy_pca)
  yw_data2 <- t(yc) %*% t(ginv(complete_miRNA_pca)) %*% t(incomplete_miRNA_pca)
  
  yw <- t((yw_data1 + yw_data2)/2)
  
  y<-rbind(yc,yw)
  
  return(y)
}
############################Incomplete data######################################################
file1 <- "BIC" #COAD GBM KRCCC LSCC BIC Bladder Brain
file2 <- "BREAST" #COLON GLIO KIDNEY LUNG BREAST BLCA LGG

file1_name <- paste0("data/", file1, "/partial_data/", file2, "_Gene_0.1.txt")
file2_name <- paste0("data/", file1, "/", file2, "_Methy_Expression.txt")
file3_name <- paste0("data/", file1, "/", file2, "_Mirna_Expression.txt")

#file1_name <- paste0("data/", file1, "/", file2, "_Gene_Expression.txt")
#file2_name <- paste0("data/", file1, "/partial_data/", file2, "_Methy_0.7.txt")
#file3_name <- paste0("data/", file1, "/", file2, "_Mirna_Expression.txt")

#file1_name <- paste0("data/", file1, "/", file2, "_Gene_Expression.txt")
#file2_name <- paste0("data/", file1, "/", file2, "_Methy_Expression.txt")
#file3_name <- paste0("data/", file1, "/partial_data/", file2, "_Mirna_0.1.txt")

#file1_name <- paste0("data/", file1, "/partial_data/", file2, "_miRNA_0.1.txt")
#file2_name <- paste0("data/", file1, "/", file2, "_protein.txt")
#file3_name <- paste0("data/", file1, "/", file2, "_RNAseq.txt")

#file1_name <- paste0("data/", file1, "/", file2, "_miRNAseq.txt")
#file2_name <- paste0("data/", file1, "/partial_data/", file2, "_protein_0.7.txt")
#file3_name <- paste0("data/", file1, "/", file2, "_RNAseq.txt")

#file1_name <- paste0("data/", file1, "/", file2, "_miRNAseq.txt")
#file2_name <- paste0("data/", file1, "/", file2, "_protein.txt")
#file3_name <- paste0("data/", file1, "/partial_data/", file2, "_RNA_0.7.txt")

v1<- read.table(file1_name, header = T, sep = "\t", row.names = 1)
v2<- read.table(file2_name, header = T, sep = "\t", row.names = 1)
v3<- read.table(file3_name, header = T, sep = "\t", row.names = 1)

all_names <- names(v1)
complete_names <- list()
incomplete_names <- list()

for(i in seq(1, length(all_names))){
  if (v1[1,i] == 0){
    incomplete_names[as.character(i)] <- all_names[i]
  }
  else{
    complete_names[as.character(i)] <- all_names[i]
  }
}

data<-list(mRNA=v1,Methy=v2,miRNA=v3)

#survive analyze
emb <- get_emb(data,complete_names,all_names,dim=4,dim_pca=15,pca_scale=F)
cl_s <- specClust(emb, 5, nn=20)
labels <- cl_s$cluster
names(labels) <- names(v1)

survfile <- paste("data/", file1, "/", file2, "_Survival.txt", sep="")
surv <- read.table(survfile, header = T)

survresult<-survdiff(Surv(Survival, Death)~labels, data=surv)
#enrichment
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
para_num <- enrich(folder=file1, sign=T, label=labels)

print(p.val)
print(para_num)
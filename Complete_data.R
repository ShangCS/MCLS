library(survival)
library("MASS")
source("run_data_col.R")
source("specClust.R")
source("enrichment.R")
assignInNamespace('specClust',specClust,ns='kknn')
environment(specClust) <- asNamespace('kknn')



############################complete data######################################################
file1 <- "BIC" #COAD  GBM   KRCCC   LSCC BIC    Bladder Brain
file2 <- "BREAST" #COLON GLIO  KIDNEY  LUNG BREAST BLCA    LGG

file1_name <- paste0("data/", file1, "/", file2, "_Gene_Expression.txt")
file2_name <- paste0("data/", file1, "/", file2, "_Methy_Expression.txt")
file3_name <- paste0("data/", file1, "/", file2, "_Mirna_Expression.txt")

#file1_name <- paste0("data/", file1, "/", file2, "_miRNAseq.txt")
#file2_name <- paste0("data/", file1, "/", file2, "_protein.txt")
#file3_name <- paste0("data/", file1, "/", file2, "_RNAseq.txt")

v1<- read.table(file1_name, header = T, sep = "\t", row.names = 1)
v2<- read.table(file2_name, header = T, sep = "\t", row.names = 1)
v3<- read.table(file3_name, header = T, sep = "\t", row.names = 1)

data<-list(mRNA=v1,Methy=v2,miRNA=v3,clinical=NULL)

simul_result<-data_col(data,incomplete_data=F,incomplete_sample_name,remain_view=1,dim1=13,dim2=4,pca_scale=F)
emb<- simul_result$l_space

cl_s <- specClust(emb, 5, nn=20)

#survive analyze
labels <- cl_s$cluster
names(labels) <- names(v1)
survfile <- paste0("data/", file1, "/", file2, "_Survival.txt")
surv <- read.table(survfile, header = T, sep = "\t")
survresult<-survdiff(Surv(Survival, Death)~labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)

#enrichment
para_num <- enrich(folder=file1, sign = T, label=labels)

print(p.val)
print(para_num)
data_col<-function(data,incomplete_data=F,incomplete_sample_name=F,remain_view=F,dim1=2,dim2=2,pca_scale=F,seed=0){
  if (incomplete_data==F){
    print("This time is for complete data (For clustering).",quote=FALSE)
    data_col_result<-list()
    mRNA_pca<-prcomp(t(data[["mRNA"]]), center=pca_scale, scale=pca_scale)
    Methy_pca<-prcomp(t(data[["Methy"]]), center=pca_scale, scale=pca_scale)
    miRNA_pca<-prcomp(t(data[["miRNA"]]), center=pca_scale, scale=pca_scale)
    f<-list(mRNA_pca[["rotation"]],Methy_pca[["rotation"]],miRNA_pca[["rotation"]])
    
    completed_data<-cbind(mRNA_pca$x[,1:dim1],Methy_pca$x[,1:dim1],miRNA_pca$x[,1:dim1])
    data_svd<-svd(completed_data)
    U1<-data_svd[["u"]][,1:dim2]
    
    g<-list(ginv(mRNA_pca$x[,1:dim1]) %*% U1,
            ginv(Methy_pca$x[,1:dim1]) %*% U1,
            ginv(miRNA_pca$x[,1:dim1]) %*% U1)
    return(list(l_space=U1,f=f,g=g))
    
  }
  else{
    print("This time is for incomplete data.",quote=FALSE)
    data_col_result<-list()
    sample_names<-colnames(data[[1]])
    complete_sample_name<-setdiff(sample_names,incomplete_sample_name)
    completed_data<-list()
    for (i in seq(1,length(data)-1)){
      completed_data[[i]]<-subset(data[[i]], select=complete_sample_name)
    }
    mRNA_pca <- prcomp(t(completed_data[[1]]), center=pca_scale, scale=pca_scale)
    Methy_pca <- prcomp(t(completed_data[[2]]), center=pca_scale, scale=pca_scale)
    miRNA_pca <- prcomp(t(completed_data[[3]]), center=pca_scale, scale=pca_scale)
    f<-list(mRNA_pca[["rotation"]],Methy_pca[["rotation"]],miRNA_pca[["rotation"]])
    integrated_data<-cbind(mRNA_pca$x[,1:dim1],Methy_pca$x[,1:dim1],miRNA_pca$x[,1:dim1])
    data_svd<-svd(integrated_data)
    U1<-data_svd[["u"]][,1:dim2]
    
    g<-list(ginv(mRNA_pca$x[,1:dim1]) %*% U1,
            ginv(Methy_pca$x[,1:dim1]) %*% U1,
            ginv(miRNA_pca$x[,1:dim1]) %*% U1)
    
    random_selected_view<-list()
    set.seed(seed)
    for (i in incomplete_sample_name){
      random_selected_view[[i]]<-sample(c(1,2,3),size = remain_view)
    }
    print("Random views selection finished.",quote=FALSE)
    if (remain_view == 2){
      incomlete_data1<-list()
      for (i in incomplete_sample_name){
        incomlete_data1[[i]]<-data[[as.numeric(random_selected_view[[i]][1])]][,i]
      }
      incomlete_data2<-list()
      for (i in incomplete_sample_name){
        incomlete_data2[[i]]<-data[[as.numeric(random_selected_view[[i]][2])]][,i]
      }
      
      
      converted_incomplete_data<-list()
      for (i in incomplete_sample_name){
        a<-(t(incomlete_data1[[i]]) %*% 
              as.matrix(f[[as.numeric(random_selected_view[[i]][1])]]))[,1:dim1] %*% 
          as.matrix(g[[as.numeric(random_selected_view[[i]][1])]])
        b<-(t(incomlete_data2[[i]]) %*% 
              as.matrix(f[[as.numeric(random_selected_view[[i]][2])]]))[,1:dim1] %*% 
          as.matrix(g[[as.numeric(random_selected_view[[i]][2])]])
        converted_incomplete_data[[i]]<-(a+b)/2
      }
      U<-U1
      for(i in incomplete_sample_name){
        U<-rbind(U,converted_incomplete_data[[i]])
      }
      print("Data collaboration finished",quote=FALSE)
      return(list(l_space=U,f=f,g=g,views=random_selected_view))
    }
    else{
      incomlete_data1<-list()
      for (i in incomplete_sample_name){
        incomlete_data1[[i]]<-data[[as.numeric(random_selected_view[[i]][1])]][,i]
      }
      print(incomlete_data1[1])
      converted_incomplete_data<-list()
      for (i in incomplete_sample_name){
        converted_incomplete_data[[i]]<-(t(incomlete_data1[[i]]) %*% 
              as.matrix(f[[as.numeric(random_selected_view[[i]][1])]]))[,1:dim1] %*% 
              as.matrix(g[[as.numeric(random_selected_view[[i]][1])]])
      }
      U<-U1
      for(i in incomplete_sample_name){
        U<-rbind(U,converted_incomplete_data[[i]])
      }
      print("Data collaboration finished",quote=FALSE)
      return(list(l_space=U,f=f,g=g,views=random_selected_view))
    }
  }
}
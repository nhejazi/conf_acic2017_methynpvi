#Find optimal combination of outcome variables
library(r2weight)
library(SummarizedExperiment)

data<-data.frame(assays(neuroIQse)$meth450k)
meta<-data.frame(colData(neuroIQse))

#Remove samples with missing outcomes
meta_comp<-meta[(!is.na(meta$vciq_i_7y) & !is.na(meta$priq_i_7y) & !is.na(meta$wmiq_i_7y) & !is.na(meta$psiq_7y) & !is.na(meta$fsiq_i_7y) & !is.na(meta$baattss_7y) & !is.na(meta$pbdss_7y) & !is.na(meta$ftdr_7y)),]
data_comp<-data[,match(row.names(meta_comp), names(data))]

Y<-data.frame(meta_comp[,1:8])

#With all covariates... 
#out1 <- optWeight(Y = Y, X = data.frame(t(data_comp)), SL.library = c("SL.mean"))

#R2Weight based on clusters:
data_comp_island<-cbind.data.frame(clusterCpGs,data_comp)

#Slow
data_comp_isl_split<-split(data_comp, as.factor(clusterCpGs))
optWeight_res <- vector("list", length(data_comp_isl_split))

for(i in 1:length(data_comp_isl_split)){
  
  X<-data.frame(t(data_comp_isl_split[[i]]))
  optWeight_res[[i]]<-optWeight(Y = Y, X = X, SL.library = c("SL.glm","SL.mean","SL.step"))

}






















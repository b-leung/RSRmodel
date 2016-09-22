# RSRmodel
Relative Suitability Richness Model
requires mgcv library
sample code to run model. 


library(pROC)

ec<-1 #environmental columns
sp<-2:6 #biotic columns (presence/absence)
nis=6 #nis column (in prediction data set, can fill with NA)
fit_dat=read.csv('test_fit.csv') # fitting data set
predict_dat=read.csv('test_pred.csv') #prediction data set
protocolB=T #do protocolB (protocol A always done)
mem_nis=predict_dat[,nis]

predict_dat[,nis]=NA #setting all nis to NA, so that RSR fills in values, for validation. Otherwise, uses real information when available
source('RSR_protocols.R')

R<-RSR(ec,sp,nis,fit_dat,predict_dat,protocolB)

sdm_roc=roc(mem_nis,R$pred_sdm)$auc
cdm_roc=roc(mem_nis,R$pred_cdm)$auc
print(paste(sdm_roc,cdm_roc)) 

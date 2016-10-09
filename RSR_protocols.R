#NOTE: Column names must match in fitting and prediction data sets
library(mgcv, quietly=T)

#ec=environment columns
#sp=species presence/absence columns
# nis=column number for the NIS presence/absence
#fit_dat=input fitting data
#predict_dat=input prediction data
#protocolB - true/false to run protocol B
RSR<-function(ec,sp,nis,fit_dat,predict_dat,do_protocolB)
{
	rnd_avg=1000
	nis_pos=which(sp==nis)

	#get variable names, plus the formula fragment
	ec_nm=names(fit_dat)[ec]
	sp_nm=names(fit_dat)[sp]
	nis_nm=names(fit_dat)[nis]
	ec_form=paste(paste("s(",ec_nm,")",sep=''),collapse='+')
	sp_form=paste(sp_nm[-nis_pos],collapse="+")

	pres_data=as.matrix(predict_dat[,sp])
	nmiss_inv=sum(is.na(predict_dat[,nis]))
	n_sp=length(sp)
	#protocol A
	protocolA<-function()
	{
		retA<-list()
		#Fitting set
		#Step 1
		rich=apply(fit_dat[,sp],1,sum) 
		fit_rich=gam(formula(paste("rich ~ ", ec_nm)), family="poisson",data=fit_dat) #all species
		#Step 2
		sdm=list()
		#note, it assumes that there is some variation in presence/absence, otherwise will crash
		for(i in 1:n_sp)
		{
			sdm[[i]]<-gam(formula(paste("fit_dat$",sp_nm[i],"~",ec_form,sep='')), family=binomial(link = "logit"),data=fit_dat)
		}

		cdm<-gam(formula(paste("fit_dat$",nis_nm,"~",ec_form,"+",sp_form,sep='')), family=binomial(link = "logit"),data=fit_dat)

		#step 3
		missing=is.na(pres_data) # find all missing values (NA)

		#step 4
		sdm_pred=matrix(0,ncol=n_sp,nrow=nrow(pres_data))
		for(i in 1:length(sp))
		{
			sdm_pred[,i]=predict.gam(sdm[[i]],type="response",newdata=predict_dat)
		}

		#step 5a
		predict_rich=predict.gam(fit_rich,type="response",newdata=predict_dat)
		#make sure are within bounds, since are using approximate fitting
		predict_rich[predict_rich>n_sp]=n_sp
		predict_rich[predict_rich<0]=0

		#step 5-7
		cdm_pred=pres_data[,nis_pos]#replace missing values below
		#RSR protocols depend on the following functions
		monte_carlo<-function(miss,rich,pred,dat,nis_pos,full_dat,niter)
		{
			cdm_rnd_pres=dat
			#find missing locations
			miss_inv=miss[,nis_pos]
			avg_miss=rep(0,sum(miss_inv))
			tmp_dat=full_dat[miss_inv,]
			for(i2 in 1:niter)			
			{
				#fill in 0/1 for all missing data based on probabilities. 
				#Step 5
				cdm_rnd_pres[miss]=choose_sp_loc(miss,rich,pred,dat)[miss]
				#Step 6 
				tmp_dat[,sp]=cdm_rnd_pres[miss_inv,]
				tmp_gam=predict.gam(cdm,type="response",newdata=tmp_dat)
				#just in case it goes beyond bounds
				tmp_gam[tmp_gam>n_sp]=n_sp
				tmp_gam[tmp_gam<0]=0
				avg_miss[miss_inv]=avg_miss[miss_inv]+tmp_gam
			}
			#update_cdm_pred
			return(avg_miss[miss_inv]/niter)

		}

		if(nmiss_inv>0)
			cdm_pred[missing[,nis_pos]]=monte_carlo(missing,predict_rich,sdm_pred,pres_data,nis_pos,predict_dat,rnd_avg)

		retA$sdm=sdm[[nis_pos]]
		retA$cdm=cdm
		retA$pred_sdm=sdm_pred[,nis_pos]
		retA$pred_cdm=cdm_pred
		return(retA)
	}
	ret=protocolA()	
	
	protocolB<-function()
	{
		#Fitting set: protocol B
		#Step 1
		retB=list()
		nat_rich=apply(fit_dat[,sp[-nis_pos]],1,sum) 
		fit_nat_rich=gam(formula(paste("nat_rich ~ ", ec_nm,"+",nis_nm)), family="poisson",data=fit_dat) #all native species using env and the invader as predictors
		#Step 2 (step 2a, is protocol A)
		#Step 2b
		if(nmiss_inv>0)
		{#only need to do sampling if had missing predictues
			rnd_avg=1
		}
		natr=0
		tmp_dat=predict_dat
		for(i in 1:rnd_avg)
		{
			tmp_i=rbinom(nrow(pres_data),1,ret$pred_cdm)
			tmp_dat[,nis]=tmp_i
			tmp_gam=predict.gam(fit_nat_rich,type="response",newdata=tmp_dat)
			tmp_gam[tmp_gam>(n_sp-1)]=n_sp-1
			tmp_gam[tmp_gam<0]=0
			natr=natr+tmp_gam/rnd_avg
		}
		#step 2c
		sum_natr=sum(natr)

		#step 3
		tmp_dat[,nis]=0
		cf_rich=predict.gam(fit_nat_rich,type="response",newdata=tmp_dat)
		cf_rich[cf_rich>(n_sp-1)]=n_sp-1
		cf_rich[cf_rich<0]=0
		sum_cf=sum(cf_rich)
		#step 4
		pred_loss=sum_cf-sum_natr

		retB$pred_loss=pred_loss
		retB$nat_rich=sum_natr
		retB$cf_rich=sum_cf
		return(retB)

	}
	if(do_protocolB==T)
	{
		ret$B<-protocolB()
	}
	return(ret)
}

choose_sp_loc<-function(miss, rch,pred, dat)
{	#figure out empty spots predicted
	mx=ncol(miss)
# Protocol A
#Step 5a - second part
	rch<-as.integer(rch)+rbinom(length(rch),1,rch-as.integer(rch))
	#but make sure does not exceed total number possible
	rch[rch>mx]=mx
#Step 5b
	rch=rch-apply(dat, 1, function(x){return(sum(x,na.rm=TRUE))}) #sum all ones that are not NA
	fill=matrix(0,nrow=nrow(miss),ncol=ncol(miss))

#Step 5c
		#which species
	for(i in 1:nrow(miss))
	{
		if(rch[i]>0 && sum(pred[i,miss[i,]])>0)
		{	
			if(sum(miss[i,])<=rch[i]) # everything known
			{
				fill[i,which(miss[i,])]=1
			}
			else
			{
				mx1=sum(pred[i,miss[i,]]>0)
				if(rch[i]>mx1)
					rch[i]=mx1
				pos=sample(which(miss[i,]),rch[i],replace=FALSE,pred[i,miss[i,]]) #these become present
				fill[i,pos]=1 
			}
			#the others remain zero
		}
	}
	return(fill)
	
}





#======================
# utility functions 
#======================


#    preProcessing
#=====================================================================================
#=====================================================================================

preProcessing<-function(inputParameters=inputParameters){
		folder.output=inputParameters$folder.output
		directory=inputParameters$directory
		directory.cluster=inputParameters$directory.cluster
		clstr=inputParameters$clstr
		NormalizedData=inputParameters$NormalizedData
		ExperimentalDesign=inputParameters$ExperimentalDesign
		Significant=inputParameters$Significant
		Controlnormalization=inputParameters$Controlnormalization
		qfactors=inputParameters$qfactors
		
		dataAnalysis<-list()
		
		# Set output folder name
		if (is.null(folder.output)){
			folder.output<-tempname
		}
		dataAnalysis$folder.output<-folder.output
		# Create new directory
		if (length(folder.output)==0){
			folder.output<-paste("NLCAResults",format(Sys.time(),"%d%m%Y%H%M%OS"),sep="")
		}
		dir.create(folder.output)
		#Reading the input files
		#Reading the file with the normalized data
		data<-read.table(NormalizedData, header = TRUE, sep = "\t", 
			stringsAsFactors = default.stringsAsFactors(),row.names=1,check.names=FALSE)
		dataAnalysis$rawdata<-data
		#Reading the experimental design file
		expdes<-NULL	# Initializate the experimental design matrix
		for (aa in 1:length(ExperimentalDesign)){
			temp<-read.table(ExperimentalDesign[aa], header = FALSE, sep = "\t", stringsAsFactors = default.stringsAsFactors())	
			expdes<-rbind(expdes,temp)
		}
		dataAnalysis$expdes<-expdes
		#Read the significant data
		Significant<-read.table(Significant, header = TRUE, sep = "\t", stringsAsFactors = default.stringsAsFactors(),row.names=1)
		dataAnalysis$Significant<-Significant
		#Determine the number of transcription factors
		ntf<-which(colnames(data)=="repeats")-1
		tf_ids<-colnames(data)[1:ntf]
		dataAnalysis$ntf<-ntf
		dataAnalysis$tf_ids<-tf_ids
		#Determine the number of factors
		nfactor<-ncol(expdes)-2
		dataAnalysis$nfactor<-nfactor
		#Determine the number of levels for each factor
		if (nfactor==1) {
			factorlevels<-length(levels(expdes[,3]))
		} else{
			factorlevels<-apply(expdes[,3:ncol(expdes)],2,function(x) {length(unique(x))})
		}
		dataAnalysis$factorlevels<-factorlevels
		#Names of the levels for each factor
		if (nfactor==1) {
			factorlevels_ids<-levels(expdes[,3])
		} else{
			factorlevels_ids<-apply(expdes[,3:ncol(expdes)],2,function(x) {levels(as.factor(x))})
		}
		treatment_ids<-levels(expdes[,3])
		dataAnalysis$factorlevels_ids<-factorlevels_ids
		dataAnalysis$treatment_ids<-treatment_ids
		#Number of repeats
		nrepeats_TF<-sum(!is.na(data[which(as.logical((data$time==1)& (data$experiment==1)& (data$treatment==1))),1]))
		nrepeats_TA<-sum(!is.na(data[which(as.logical((data$time==1)& (data$experiment==1)& (data$treatment==1))),name_TA]))
		if (length(which(colnames(data)==name_no_dna))>0){
			nrepeats_NODNA<-sum(!is.na(data[which(
				as.logical((data$time==1)& (data$experiment==1)& (data$treatment==1))),name_no_dna]))
		} else {
			nrepeats_NODNA=0
		}
		dataAnalysis$nrepeats_TF<-nrepeats_TF
		dataAnalysis$nrepeats_TA<-nrepeats_TA
		dataAnalysis$nrepeats_NODNA<-nrepeats_NODNA		
		nrepeats<-max(nrepeats_TF,nrepeats_TA,nrepeats_NODNA)
		dataAnalysis$nrepeats<-nrepeats		
		#Number of timepoints
		ntime<-length(unique(data$time))
		dataAnalysis$ntime<-ntime	
		#Number of treatments
		ntreatment<-length(unique(data$treatment))
		dataAnalysis$ntreatment<-ntreatment		
		#Number of experiments
		nexp<-length(unique(data$experiment))
		dataAnalysis$nexp<-nexp			
		#Normalization control
		if (Controlnormalization==1){
			tempnor<-data
			colnames(tempnor)[1:ntf]<-tf_ids
			rownames(tempnor)<-rownames(data)
			for (aa in 1:ntime){
				for (bb in 1:ntreatment){
					for (dd in 1:ntf){
						for (ee in 1:nexp){
						tempnor[which(data$time==aa & data$treatment==bb & data$experiment==ee),dd]<-
							data[which(data$time==aa & data$treatment==bb & data$experiment==ee),dd]-
							mean(data[which(data$time==aa & data$treatment==which(treatment_ids==name_untreated)& data$experiment==ee),dd])
						}
					}
				}
			}
			data<-tempnor
		}
		dataAnalysis$data<-data	
		#Calculate the mean and standard deviation
		#Just by repeats
		average<-matrix(0,ntime*ntreatment*nexp,ntf)
		stdev<-matrix(0,ntime*ntreatment*nexp,ntf)
		count_cc=0
		for (aa in 1:nexp){
			for (bb in 1:ntreatment){
				for (cc in 1:ntime){
					count_cc=count_cc+1
					for (dd in 1:ntf){
						average[count_cc,dd]<-mean(data[which(as.logical((data$experiment==aa) & (data$treatment==bb) & (data$time==cc))),dd],na.rm=TRUE)
						stdev[count_cc,dd]<-sd(data[which(as.logical((data$experiment==aa) & (data$treatment==bb) & (data$time==cc))),dd],na.rm=TRUE)
					}
				}
			}
		}
		colnames(average)[1:ntf]<-tf_ids
		colnames(stdev)[1:ntf]<-tf_ids
		average<-data.frame(average,data[which(as.logical(data$repeats==1)),(ntf+1):ncol(data)])
		stdev<-data.frame(stdev,data[which(as.logical(data$repeats==1)),(ntf+1):ncol(data)])
		dataAnalysis$average<-average	
		dataAnalysis$stdev<-stdev	
		#Ignoring different experiments
		average_all<-matrix(0,ntime*ntreatment,ntf)
		stdev_all<-matrix(0,ntime*ntreatment,ntf)
		count_cc=0
		for (bb in 1:ntreatment){
			for (cc in 1:ntime){
				count_cc=count_cc+1
				for (dd in 1:ntf){
					average_all[count_cc,dd]<-mean(data[which(as.logical((data$treatment==bb)&(data$time==cc))),dd],
					na.rm=TRUE)
					stdev_all[count_cc,dd]<-sd(data[which(as.logical((data$treatment==bb)&(data$time==cc))),dd],
					na.rm=TRUE)
				}
			}
		}
		colnames(average_all)[1:ntf]<-tf_ids
		colnames(stdev_all)[1:ntf]<-tf_ids
		average_all<-data.frame(average_all,data[which(as.logical((data$repeats==1)&(data$experiment==1))),(ntf+1):ncol(data)])
		stdev_all<-data.frame(stdev_all,data[which(as.logical((data$repeats==1)&(data$experiment==1))),(ntf+1):ncol(data)])
		dataAnalysis$average_all<-average_all	
		dataAnalysis$stdev_all<-stdev_all	
		#Ignoring time
		average_time<-data[which(data$time==1),]
		stdev_time<-data[which(data$time==1),]
		for (aa in 1:nexp){
			for (bb in 1:ntreatment){
				for (cc in 1:nrepeats){
					for (dd in 1:ntf){
						average_time[which(as.logical((average_time$experiment==aa)&(average_time$treatment==bb)&(average_time$repeats==cc))),dd]<-
							mean(data[which(as.logical((data$experiment==aa)&(data$treatment==bb)&(data$repeats==cc))),dd],na.rm=TRUE)
						stdev_time[which(as.logical((stdev_time$experiment==aa)&(stdev_time$treatment==bb)&(stdev_time$repeats==cc))),dd]<-
							sd(data[which(as.logical((data$experiment==aa)&(data$treatment==bb)&(data$repeats==cc))),dd],na.rm=TRUE)
					}
				}
			}
		}
		colnames(average_time)[1:ntf]<-tf_ids
		colnames(stdev_time)[1:ntf]<-tf_ids
		dataAnalysis$average_time<-average_time	
		dataAnalysis$stdev_time<-stdev_time
		#Calculate fold change versus non-treated
		FcvsUntreated<-average_all
		for (ii in 1:ntreatment){
			for (jj in 1:ntime){
				FcvsUntreated[which(as.logical(average_all$Treatment==ii & average_all$time==jj)),]<-
					average_all[which(as.logical(average_all$Treatment==ii & average_all$time==jj)),]-
					average_all[which(as.logical(average_all$Treatment==which(treatment_ids==name_untreated) & average_all$time==jj)),]
			}
		}
		FcvsUntreated[,1:ntf]<-apply(FcvsUntreated[,1:ntf],c(1,2),function(x){if (x==0) x=1 else {sign(x)*2^abs(x)}})
		FcvsUntreated<-FcvsUntreated[-which(as.logical(average_all$Treatment==which(treatment_ids==name_untreated)| average_all$time==timevector[1])),1:ntf]
		FcvsUntreated<-t(FcvsUntreated)
		dataAnalysis$FcvsUntreated<-FcvsUntreated
		#Generation of variables due to qualitative nature
		ncolq=0 #dimensions of the new factorial matrix
		for (pp in 1:nfactor){
			if(qfactors[pp]==1){
				ncolq=ncolq+factorlevels[pp]-1
			} else {
				ncolq=ncolq+1
			}
		}
		qfactormatrix<-matrix(0,nrow(data),ncolq)
		qaveragematrix<-matrix(0,nrow(average),ncolq)
		qstdevmatrix<-matrix(0,nrow(stdev),ncolq)
		qaverage_allmatrix<-matrix(0,nrow(average_all),ncolq)
		qstdev_allmatrix<-matrix(0,nrow(stdev_all),ncolq)
		qaverage_timematrix<-matrix(0,nrow(average_time),ncolq)
		qstdev_timematrix<-matrix(0,nrow(stdev_time),ncolq)
		qfactormatrix_names<-NULL
		oo=1
		for (nn in 1:length(qfactors)){
			if(qfactors[nn]==1){
				for (mm in 1:factorlevels[nn]){
					if ((nfactor==1)&& (factorlevels_ids[mm]!=baseline[nn])){
						qfactormatrix[which(data[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qaveragematrix[which(average[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qstdevmatrix[which(stdev[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qaverage_allmatrix[which(average_all[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qstdev_allmatrix[which(stdev_all[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qaverage_timematrix[which(average_time[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qstdev_timematrix[which(stdev_time[,ntf+3+nn]==which(factorlevels_ids==factorlevels_ids[mm])),oo]<-1 
						qfactormatrix_names<-c(qfactormatrix_names,factorlevels_ids[mm])
						oo=oo+1
					}else {
						if ((nfactor>1)&& (factorlevels_ids[nn][[1]][mm]!=baseline[nn])){
							qfactormatrix[which(data[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1
							qaveragematrix[which(average[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1
							qstdevmatrix[which(stdev[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1 
							qaverage_allmatrix[which(average_all[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1
							qstdev_allmatrix[which(stdev_all[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1 
							qaverage_timematrix[which(average_time[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1
							qstdev_timematrix[which(stdev_time[,ntf+3+nn]==which(factorlevels_ids[nn][[1]]==factorlevels_ids[nn][[1]][mm])),oo]<-1 
							qfactormatrix_names<-c(qfactormatrix_names,factorlevels_ids[mm])
							oo=oo+1
							}
					} 
				}
			}else{
				qfactormatrix[,oo]<-data[,ntf+3+nn]
				qaveragematrix[,oo]<-average[,ntf+3+nn]
				qstdevmatrix[,oo]<-stdev[,ntf+3+nn]
				qaverage_allmatrix[,oo]<-average_all[,ntf+3+nn]
				qstdev_allmatrix[,oo]<-stdev_all[,ntf+3+nn]
				qaverage_timematrix[,oo]<-average_time[,ntf+3+nn]
				qstdev_timematrix[,oo]<-stdev_time[,ntf+3+nn]
				qfactormatrix_names<-c(qfactormatrix_names,colnames(data[ntf+3+nn]))
			}
		}
		colnames(qfactormatrix)<-qfactormatrix_names
		colnames(qaveragematrix)<-qfactormatrix_names
		colnames(qstdevmatrix)<-qfactormatrix_names
		colnames(qaverage_allmatrix)<-qfactormatrix_names
		colnames(qstdev_allmatrix)<-qfactormatrix_names
		colnames(qaverage_timematrix)<-qfactormatrix_names
		colnames(qstdev_timematrix)<-qfactormatrix_names
		dataAnalysis$qfactormatrix<-qfactormatrix
		dataAnalysis$qaveragematrix<-qaveragematrix	
		dataAnalysis$qstdevmatrix<-qstdevmatrix		
		dataAnalysis$qaverage_allmatrix<-qaverage_allmatrix		
		dataAnalysis$qaverage_timematrix<-qaverage_timematrix
		dataAnalysis$qstdev_timematrix<-qstdev_timematrix
		dataq<-data.frame(data[,1:(ntf+3)],experiment=data[,(ntf+nfactor+4)],qfactormatrix)
		averageq<-data.frame(average[,1:(ntf+3)],experiment=average[,(ntf+nfactor+4)],qaveragematrix)
		stdevq<-data.frame(stdev[,1:(ntf+3)],experiment=stdev[,(ntf+nfactor+4)],qstdevmatrix)
		average_allq<-data.frame(average_all[,1:(ntf+3)],experiment=average_all[,(ntf+nfactor+4)],qaverage_allmatrix)
		stdev_allq<-data.frame(stdev_all[,1:(ntf+3)],experiment=stdev_all[,(ntf+nfactor+4)],qstdev_allmatrix)
		average_timeq<-data.frame(average_time[,1:(ntf+3)],experiment=average_time[,(ntf+nfactor+4)],qaverage_timematrix)
		stdev_timeq<-data.frame(stdev_time[,1:(ntf+3)],experiment=stdev_time[,(ntf+nfactor+4)],qstdev_timematrix)
		dataAnalysis$dataq<-dataq
		dataAnalysis$averageq<-averageq
		dataAnalysis$stdevq<-stdevq
		dataAnalysis$average_allq<-average_allq
		dataAnalysis$stdev_allq<-stdev_allq
		dataAnalysis$average_timeq<-average_timeq		
		dataAnalysis$stdev_timeq<-stdev_timeq			
		if (nfactor==1){
			dataf<-data.frame(data[,1:(ntf+3)],experiment=data[,(ntf+nfactor+4)],data[,(ntf+4)])
			averagef<-data.frame(average[,1:(ntf+3)],experiment=average[,(ntf+nfactor+4)],average[,(ntf+4)])
			stdevf<-data.frame(stdev[,1:(ntf+3)],experiment=stdev[,(ntf+nfactor+4)],stdev[,(ntf+4)])
			average_allf<-data.frame(average_all[,1:(ntf+3)],experiment=average_all[,(ntf+nfactor+4)],average_all[,(ntf+4)])
			stdev_allf<-data.frame(stdev_all[,1:(ntf+3)],experiment=stdev_all[,(ntf+nfactor+4)],stdev_all[,(ntf+4)])
			average_timef<-data.frame(average_time[,1:(ntf+3)],experiment=average_time[,(ntf+nfactor+4)],average_time[,(ntf+4)])
			stdev_timef<-data.frame(stdev_time[,1:(ntf+3)],experiment=stdev_time[,(ntf+nfactor+4)],stdev_time[,(ntf+4)])
			colnames(dataf)[ncol(dataf)]<-factornames
			colnames(averagef)[ncol(averagef)]<-factornames
			colnames(stdevf)[ncol(stdevf)]<-factornames
			colnames(average_allf)[ncol(average_allf)]<-factornames
			colnames(stdev_allf)[ncol(stdev_allf)]<-factornames
			colnames(average_timef)[ncol(average_timef)]<-factornames
			colnames(stdev_timef)[ncol(stdev_timef)]<-factornames
		}else{
			dataf<-data.frame(data[,1:(ntf+3)],experiment=data[,(ntf+nfactor+4)],data[,(ntf+4):(ntf+nfactor)])
			averagef<-data.frame(average[,1:(ntf+3)],experiment=average[,(ntf+nfactor+4)],average[,(ntf+4):(ntf+nfactor)])
			stdevf<-data.frame(stdev[,1:(ntf+3)],experiment=stdev[,(ntf+nfactor+4)],stdev[,(ntf+4):(ntf+nfactor)])
			average_allf<-data.frame(average_all[,1:(ntf+3)],experiment=average_all[,(ntf+nfactor+4)],average_all[,(ntf+4):(ntf+nfactor)])
			stdev_allf<-data.frame(stdev_all[,1:(ntf+3)],experiment=stdev_all[,(ntf+nfactor+4)],stdev_all[,(ntf+4):(ntf+nfactor)])
			average_timef<-data.frame(average_time[,1:(ntf+3)],experiment=average_time[,(ntf+nfactor+4)],average_time[,(ntf+4):(ntf+nfactor)])
			stdev_timef<-data.frame(stdev_time[,1:(ntf+3)],experiment=stdev_time[,(ntf+nfactor+4)],stdev_time[,(ntf+4):(ntf+nfactor)])
			colnames(dataf)[(ncol(dataf)-nfactor+1):ncol(dataf)]<-factornames
			colnames(averagef)[(ncol(averagef)-nfactor+1):ncol(averagef)]<-factornames
			colnames(stdevf)[(ncol(stdevf)-nfactor+1):ncol(stdevf)]<-factornames
			colnames(average_allf)[(ncol(average_allf)-nfactor+1):ncol(average_allf)]<-factornames
			colnames(stdev_allf)[(ncol(stdev_allf)-nfactor+1):ncol(stdev_allf)]<-factornames
		}
		dataAnalysis$dataf<-dataf
		dataAnalysis$averagef<-averagef
		dataAnalysis$stdevf<-stdevf
		dataAnalysis$average_allf<-average_allf
		dataAnalysis$stdev_allf<-stdev_allf
		return (dataAnalysis)
	}
#=====================================================================================


#     splineinterpolation
#=====================================================================================
#=====================================================================================

splineinterpolation<-function(data,method="natural",ninter=ninter,clstr=clstr,filename=filename,timevector=timevector,
	ntf=ntf,folder.output=folder.output,lib.loc=lib.loc,tf_ids = tf_ids,ntreatment=ntreatment,ntime=ntime,
	treatment_ids=treatment_ids){

	time<-timevector
	xout<-c()
	for (ss in 2:length(timevector)){
		seqlen=(timevector[ss]-timevector[ss-1])/(ninter+1)
		tmp<-seq(timevector[ss-1],timevector[ss],by=seqlen)
		xout<-c(xout,tmp)
	}
	xout<-unique(xout)
	xout<-xout[order(xout,decreasing=FALSE)]
	nexp<-length(unique(data$experiment))
	ntreat<-length(unique(data$treatment))
	ntimeint<-length(xout)
	nrep<-length(unique(data$repeats))
	temp<-c()
	for (oo in 1:nexp){
	for (pp in 1:ntreat){
		templine<-rep(data[which(data$experiment==oo & data$treatment==pp & data$time==1),],ntimeint)
		templine<-unstack(stack(templine))
		temp<-rbind(temp,templine)
	}}
	temp<-temp[,match(colnames(data),colnames(temp))]
	for (oo in 1:nexp){
	for (pp in 1:ntreat){
	for (qq in 1:ntf){
	for (rr in 1:nrep){		
		intensity<-data[which(data$experiment==oo & data$treatment==pp & data$repeats==rr),qq]
		interpolationspline<-spline(time,intensity,method=method,xout=xout)
		temp[which(temp$experiment==oo & temp$treatment==pp & temp$repeats==rr),qq]<-interpolationspline$y
		seqtem<-1/(ninter+1)
		temp$time[which(temp$experiment==oo & temp$treatment==pp & temp$repeats==rr)]<-seq(1,ntimeint,by=seqtem)
	}}}}
	temp<-temp[with(temp, order(experiment,treatment,time,repeats)), ]
	if (clstr==0){
		par(ask=TRUE)
	}
	if (clstr==0){
		library(sciplot)
	} else {
		library(sciplot,lib.loc=lib.loc)
	}
	nRow=ceiling(-1+sqrt(2+ntf))
	nCol=2+nRow
	par(mfrow=c(nRow,nCol), mar=c(2,5,0.5,1), oma=c(0,0,0,0))
	pdf(file = paste(folder.output,"/Splineinterpolation",filename,".pdf",sep=""),paper="a4r")
	for(i in 1:ntf){
		lineplot.CI(temp$time, temp[,i], group = temp$treatment,col=c(1:ntreatment),type="b",fixed=TRUE,x.cont=TRUE,
			xaxt="n", ylab=tf_ids[i],fun = function(x) mean(x, na.rm=TRUE),xlab="",legend=FALSE,leg.lab=NULL,
    			ci.fun= function(x) c(mean(x, na.rm=TRUE)-se(x,na.rm=TRUE), mean(x, na.rm=TRUE)+se(x,na.rm=TRUE))) 
			axis(1, at = 1:ntime, labels = timevector)
	}
	plot.new()
	legend("topleft",treatment_ids,col=c(1:ntreatment),lty = c(1,2,3), pch = c(1, 1, 4),
      	merge = TRUE, bg = 'gray90')
	dev.off()
	return (temp)
}

#=====================================================================================


#	wrappingSpline
#=====================================================================================
#=====================================================================================
wrappingSpline<-function(inputParameters=inputParameters,dataAnalysis){
		dataq<-dataAnalysis$dataq
		dataf<-dataAnalysis$dataf
		averageq<-dataAnalysis$averageq
		averagef<-dataAnalysis$averagef
		average_allq<-dataAnalysis$average_allq
		average_allf<-dataAnalysis$average_allf
		stdevq<-dataAnalysis$stdevq
		stdevf<-dataAnalysis$stdevf
		stdev_allq<-dataAnalysis$stdev_allq
		stdev_allf<-dataAnalysis$stdev_allf
		average_timeq<-dataAnalysis$average_timeq
		stdev_timeq<-dataAnalysis$stdev_timeq
		data<-dataAnalysis$data
		clstr<-inputParameters$clstr
		timevector<-inputParameters$timevector
		ntf<-dataAnalysis$ntf
		folder.output<-dataAnalysis$folder.output
		lib.loc<-lib.loc$folder.output
		tf_ids<-dataAnalysis$tf_ids
		ntreatment<-dataAnalysis$ntreatment
		ntime<-dataAnalysis$ntime
		treatment_ids<-dataAnalysis$treatment_ids
		
		dataInt<-list()
		
		# Spline interpolation of the data
		dataq.int<-splineinterpolation(dataq,method="natural",ninter=1,clstr=clstr,filename="dataq",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		dataf.int<-splineinterpolation(dataf,method="natural",ninter=1,clstr=clstr,filename="dataf",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		averageq.int<-splineinterpolation(averageq,method="natural",ninter=1,clstr=clstr,filename="averageq",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		averagef.int<-splineinterpolation(averagef,method="natural",ninter=1,clstr=clstr,filename="averagef",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		average_allq.int<-splineinterpolation(average_allq,method="natural",ninter=1,clstr=clstr,filename="average_allq",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		average_allf.int<-splineinterpolation(average_allf,method="natural",ninter=1,clstr=clstr,filename="average_allf",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		stdevq.int<-splineinterpolation(stdevq,method="natural",ninter=1,clstr=clstr,filename="stdevq",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		stdevf.int<-splineinterpolation(stdevf,method="natural",ninter=1,clstr=clstr,filename="stdevf",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		stdevq_all.int<-splineinterpolation(stdev_allq,method="natural",ninter=1,clstr=clstr,filename="stdevq_all",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		stdevf_all.int<-splineinterpolation(stdev_allf,method="natural",ninter=1,clstr=clstr,filename="stdevf_all",timevector=timevector,ntf=ntf,folder.output=folder.output,
			lib.loc=lib.loc,tf_ids=tf_ids,ntreatment=ntreatment,ntime=ntime,treatment_ids=treatment_ids)
		colnames(dataq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(averageq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdevq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(average_allq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdev_allq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(average_timeq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdev_timeq)[1:ntf]<-colnames(data)[1:ntf]
		colnames(dataq.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(dataf.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(averageq.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(averagef.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(average_allq.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(average_allf.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdevq.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdevf.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdevq_all.int)[1:ntf]<-colnames(data)[1:ntf]
		colnames(stdevf_all.int)[1:ntf]<-colnames(data)[1:ntf]
		dataInt$dataq<-dataq
		dataInt$averageq<-averageq
		dataInt$stdevq<-stdevq
		dataInt$average_allq<-average_allq
		dataInt$stdev_allq<-stdev_allq
		dataInt$average_timeq<-average_timeq
		dataInt$stdev_timeq<-stdev_timeq
		dataInt$dataq.int<-dataq.int
		dataInt$dataf.int<-dataf.int
		dataInt$averageq.int<-averageq.int
		dataInt$averagef.int<-averagef.int
		dataInt$stdevq.int<-stdevq.int
		dataInt$stdevf.int<-stdevf.int
		dataInt$average_allq.int<-average_allq.int
		dataInt$average_allf.int<-average_allf.int
		dataInt$stdevq_all.int<-stdevq_all.int
		dataInt$stdevf_all.int<-stdevf_all.int
		return(dataInt)
	}
	
#=====================================================================================

#	wrappingkmeans   
#=====================================================================================
#=====================================================================================

	wrappingkmeans<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		data.matrix<-as.matrix(t(dataq.int[,1:ntf]))
		data.matrix<-data.matrix[-which((rownames(data.matrix)==name_no_dna)|(rownames(data.matrix)==name_TA)),]
		average.matrix<-as.matrix(t(averageq.int[,1:ntf]))
		average.matrix<-average.matrix[-which((rownames(average.matrix)==name_no_dna)|(rownames(average.matrix)==name_TA)),]
		average_all.matrix<-as.matrix(t(average_allq.int[,1:ntf]))
		average_all.matrix<-average_all.matrix[-which((rownames(average_all.matrix)==name_no_dna)|(rownames(average_all.matrix)==name_TA)),]
		FcvsUntreated<-as.matrix(FcvsUntreated)
		clusters<-function(data=data,scaled=1){
			if(scaled==1){
				data <- t(scale(t(data))) # standardize variables
			}
			# Determine number of clusters
			wss <- (nrow(data)-1)*sum(apply(data,2,var))
			for (i in 2:(nrow(data)-1)) 
				wss[i] <- sum(kmeans(data, centers=i,iter.max = 100)$withinss)
				if (clstr==0){
				plot(1:(nrow(data)-1), wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
				}

			}
			clusters(data.matrix,1)
			clusters(average.matrix,1)
			clusters(average_all.matrix,1)
			clusters(FcvsUntreated,0)
			kmeansTF<-function(data=data,nclusters=nclusters,scaled=scaled){
			if(scaled==1){
				data <- t(scale(t(data))) # standardize variables
			}
			# K-Means Cluster Analysis
			fit <- kmeans(data, nclusters,iter.max = 100) #  cluster solution
			# get cluster means 
			aggregate(data,by=list(fit$cluster),FUN=mean)
			# append cluster assignment
			data <- data.frame(data, cluster=fit$cluster)
			data<-data[order(data$cluster),]
			return(data)
		}
		data.kmeans<-kmeansTF(data.matrix,12,1)
		average.kmeans<-kmeansTF(average.matrix,11,1)
		average_all.kmeans<-kmeansTF(average_all.matrix,12,1)
		FcvsUntreated.kmeans<-kmeansTF(FcvsUntreated,10,0)
		pdf(file = paste(folder.output,"/kmeans",".pdf",sep=""),paper="a4r")
		library(RColorBrewer,lib.loc=lib.loc); 
		mypalette <- c(brewer.pal(9,"Purples")[seq(9,1,by=-1)],brewer.pal(9,"Oranges")) 
		library("gplots")
		heatmap.2(as.matrix(data.kmeans[,(1:ncol(data.kmeans)-1)]), Rowv=F, dendrogram="none", Colv=F, col=mypalette, 
			scale="none", trace="none", key=T) 
		heatmap.2(as.matrix(average.kmeans[,(1:ncol(average.kmeans)-1)]), Rowv=F, dendrogram="none", Colv=F, col=mypalette, 
			scale="row", trace="none", key=T) 
		heatmap.2(as.matrix(average_all.kmeans[,(1:ncol(average_all.kmeans)-1)]), Rowv=F, dendrogram="none", Colv=F, col=mypalette, 
			scale="row", trace="none", key=T) 
		dev.off()
	}
	
#=====================================================================================

#	plotcorrelationmatrix     
#=====================================================================================
#=====================================================================================

plotcorrelationmatrix<-function(corr=corr){
##  Create the colors. The blue colors are designed to be lighter than pure
##  blue by adding equal amounts of red and green.
	red <-c(rep(1,16),0,0,0,0,0)
	green <- c(0,0.2,0.5,0.7,0.9,rep(1,11),0.9,0.7,0.5,0.2,0)
	blue <-c(0,0,0,0,0,rep(1,16))
	colors <- rgb(red = red, green = green, blue = blue )
##  Reverse the columns of the matrix so it will be drawn correctly.
	n = ncol(corr )
	corr2 <- corr[ , n:1 ]
##  Create the image.
	par(mar=c(5,10,10,5))
	image(
   	 z    = corr2,
    	axes = FALSE,
    	col  = colors,
    	zlim = c( -1.0, 1.0 ) )
##  Add labels for the y axis.
	axis(
    	side     = 2,
    	labels   = colnames( corr2 ),
    	at       = seq( 0, 1, length = length( colnames( corr2 ) ) ),
    	cex.axis = 1.2,
    	las      = 2)
##  Add labels for the x axis, but along the top.
	axis(
    	side     = 3,
    	labels   = rownames( corr2 ),
    	at       = seq( 0, 1, length = length( rownames( corr2 ) ) ),
    	cex.axis = 1.2,
    	las      = 2 )
	box()
	par(mar=c(5,4,4,2))
}

#=====================================================================================

#	linearInference   
#=====================================================================================
#=====================================================================================
linearInference<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		# Correlation matrix with all the data
		posTANC<-which((colnames(dataq.int)==name_no_dna)|(colnames(dataq.int)==name_TA))
		dataq.int.sel<-dataq.int[,-posTANC]
		colnames(dataq.int.sel)<-colnames(dataq.int)[-posTANC]
		dataq.int.sel<-data.frame(dataq.int.sel[,1:(which(colnames(dataq.int.sel)=="repeats")-1)],time=dataq.int.sel$time, dataq.int.sel[,(which(colnames(dataq.int.sel)=="experiment")+1):ncol(dataq.int.sel)]*dataq.int.sel$time)
		x<-dataq.int.sel[which(dataq.int.sel$time!=ntime),]	
		colnames(x)<-colnames(dataq.int.sel)
		y<-dataq.int.sel[which(dataq.int.sel$time!=1),1:(which(colnames(dataq.int.sel)=="time")-1)]
		colnames(y)<-colnames(dataq.int.sel)[1:1:(which(colnames(dataq.int.sel)=="time")-1)]
		correlation_matrix<-cor(x,y,method = "pearson",use="na.or.complete")
		colnames(correlation_matrix)<-colnames(dataq.int.sel)[1:1:(which(colnames(dataq.int.sel)=="time")-1)]
		rownames(correlation_matrix)<-colnames(dataq.int.sel)
		par(mfrow=c(1,1))
		if (clstr==0){
			par(ask=TRUE)
		} else {
			par(ask=FALSE)
		}
    	pdf(file = paste(folder.output,"/CorrelationMatrixAllData.pdf",sep=""),paper="a4r")
		plotcorrelationmatrix(correlation_matrix)
		dev.off()
		cor_0.8<-apply(correlation_matrix,c(1,2), function(x) {if(abs(x)<1 && abs(x)>=0.8) {x=1} else x=0})
		if (sum(cor_0.8)>0){
			LinearCorrelation<-NULL
			for (ee in 1:ncol(cor_0.8)){
				for (ff in ee:nrow(cor_0.8)){
					if (cor_0.8[ff,ee]==1 && ff!=ee){
						LinearCorrelation<-rbind(LinearCorrelation,cbind(rownames(cor_0.8)[ff],colnames(cor_0.8)[ee],round(correlation_matrix[ff,ee],2)))
					}
				}
			}
		} else {
			LinearCorrelation<-NULL
			print("No significant correlation found for values using linear correlation")
		}
		if(length(LinearCorrelation)>0){
			LinearCorrelation<-data.frame(LinearCorrelation)
			colnames(LinearCorrelation)<-c("TF1", "TF2","PearsonCorrelation")
			write.table(LinearCorrelation, file = paste(folder.output,"/LinearCorrelation.txt",sep=""), sep = "\t",
     			row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
			panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
				usr <- par("usr"); on.exit(par(usr))
				par(usr = c(0, 1, 0, 1))
				r <- abs(cor(x, y,method = "pearson"))
				txt <- format(c(r, 0.123456789), digits=digits)[1]
				txt <- paste(prefix, txt, sep="")
				if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
					text(0.5, 0.5, txt, cex = cex.cor * r)
			}
			panel.line<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
				cex = 1, col.line = "red") {
				abline(a=0,b=1, col = col.line)
			}
			if (clstr==0){
				par(ask=TRUE)
			}
    		pdf(file = paste(folder.output,"/PairsAllData.pdf",sep=""),paper="a4r")
			pairs(dataq.int.sel,  upper.panel=panel.cor)
			dev.off()
		}
		dataInt$dataq.int.sel<-dataq.int.sel
		# Average data
		# Remove NC and TA 
		posTANC<-which((colnames(averageq.int)=="NC")|(colnames(averageq.int)=="TA"))
		averageq.int.sel<-averageq.int[,-posTANC]
		averageq.int.sel<-data.frame(averageq.int.sel[,1:(which(colnames(averageq.int.sel)=="repeats")-1)],time=averageq.int.sel$time, averageq.int.sel[,(which(colnames(averageq.int.sel)=="experiment")+1):ncol(averageq.int.sel)]*averageq.int.sel$time)
		x<-averageq.int.sel[-which(averageq.int.sel$time==ntime),]	
		y<-averageq.int.sel[-which(averageq.int.sel$time==1),1:(which(colnames(averageq.int.sel)=="time")-1)]
		correlation_matrix_average<-cor(x,y,method = "pearson")
		dataInt$averageq.int.sel<-averageq.int.sel
		par(mfrow=c(1,1))
		if (clstr==0){
			par(ask=TRUE)
		}
    	pdf(file = paste(folder.output,"/CorrelationMatrixAverage.pdf",sep=""),paper="a4r")
		plotcorrelationmatrix(correlation_matrix_average)
		dev.off()
		cor_average_0.8<-apply(correlation_matrix_average,c(1,2), function(x) {if(abs(x)<1 && abs(x)>=0.8) {x=1} else x=0})
		if (sum(cor_average_0.8)>0){
			LinearCorrelation_average<-NULL
			for (ee in 1:ncol(cor_average_0.8)){
				for (ff in ee:nrow(cor_average_0.8)){
					if (cor_average_0.8[ff,ee]==1 && ff!=ee){
					LinearCorrelation_average<-rbind(LinearCorrelation_average,
						cbind(rownames(cor_average_0.8)[ff],colnames(cor_average_0.8)[ee],round(correlation_matrix_average[ff,ee],2)))
					}
				}
			}
		}else {
			LinearCorrelation_average<-NULL
			print("No significant correlation found for average values using linear correlation")
		}
		if (length(LinearCorrelation_average)>0){
			LinearCorrelation_average<-data.frame(LinearCorrelation_average)
			colnames(LinearCorrelation_average)<-c("TF1", "TF2","PearsonCorrelation")
			write.table(LinearCorrelation_average, file =  paste(folder.output,"/LinearCorrelation_average.txt", sep=""),sep = "\t",
     			row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
			panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
    			usr <- par("usr"); on.exit(par(usr))
				par(usr = c(0, 1, 0, 1))
				r <- abs(cor(x, y,method = "pearson"))
				txt <- format(c(r, 0.123456789), digits=digits)[1]
				txt <- paste(prefix, txt, sep="")
				if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
				text(0.5, 0.5, txt, cex = cex.cor * r)
			}
			if (clstr==0){
				par(ask=TRUE)
			}
			pdf(paste(folder.output,"/PairsAverage.pdf",sep=""),paper="a4r")
			pairs(averageq.int.sel, 
				lower.panel=panel.smooth, upper.panel=panel.cor,cex.main=1.5)
			dev.off()
		}
		# Average time
		correlation_matrix_time<-cor(average_timeq[,c(1:ntf,((which(colnames(average_timeq)=="experiment")+1):ncol(average_timeq)))],method = "pearson")
		pdf(paste(folder.output,"/CorrelationMatrixAverageTime.pdf",sep=""),paper="a4r")
		plotcorrelationmatrix(correlation_matrix_time)
		dev.off()
		cor_time_0.8<-apply(correlation_matrix_time,c(1,2), function(x) {if(abs(x)<1 && abs(x)>=0.8) {x=1} else x=0})
		if (sum(cor_time_0.8)>0){
			LinearCorrelation_time<-NULL
			for (ee in 1:ncol(cor_time_0.8)){
				for (ff in ee:nrow(cor_time_0.8)){
					if (cor_time_0.8[ff,ee]==1 && ff!=ee){
					LinearCorrelation_time<-rbind(LinearCorrelation_time,
						cbind(rownames(cor_time_0.8)[ff],rownames(cor_time_0.8)[ee],round(correlation_matrix_time[ff,ee],2)))
				}	
			}
		}
		}else {
			LinearCorrelation_time<-NULL
			print("No significant correlation found for average time values using linear correlation")
		}
		if (length(LinearCorrelation_time)>0){
			LinearCorrelation_time<-data.frame(LinearCorrelation_time)
			colnames(LinearCorrelation_time)<-c("TF1", "TF2","PearsonCorrelation")
			write.table(LinearCorrelation_time, file =  paste(folder.output,"/LinearCorrelation_time.txt", sep=""),sep = "\t",
     			row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
		}
		panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
    		usr <- par("usr"); on.exit(par(usr))
    		par(usr = c(0, 1, 0, 1))
    		r <- abs(cor(x, y,method = "pearson"))
    		txt <- format(c(r, 0.123456789), digits=digits)[1]
    		txt <- paste(prefix, txt, sep="")
    		if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		if (clstr==0){
			par(ask=TRUE)
		}
		pdf(paste(folder.output,"/PairsAverageTime.pdf",sep=""),paper="a4r")
		pairs(average_timeq[,c(1:ntf,((which(colnames(average_timeq)=="experiment")+1):ncol(average_timeq)))], 
		lower.panel=panel.smooth, upper.panel=panel.cor,cex.main=1.5)
		dev.off()
		return(dataInt)
	}
	
#=====================================================================================

#	PLSRInference   
#=====================================================================================
#=====================================================================================

PLSRInference<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		dataq.int.sig<-dataq.int[,c(match(rownames(Significant),colnames(dataq.int)),(which(colnames(dataq.int)=="experiment")+1):ncol(dataq.int))]	
		dataq.int.sig.allinfo<-dataq.int[,c(match(rownames(Significant),colnames(dataq.int)),(which(colnames(dataq.int)=="repeats")):ncol(dataq.int))]	
		dataInt$dataq.int.sig<-dataq.int.sig
		dataInt$dataq.int.sig.allinfo<-dataq.int.sig.allinfo
		# NOTE
			## Just do the analysis with one component and try to get the best model based on that to order your variables and find correlation
			## Implemented TD-PSLR as following for TFs
			##	(ntf(t)-ntf(t-1)) ~ ntf(t-1) 
		## Creation of the matrix to generate the data
		tiniplsr<-ifelse (0 %in% timevector,2,1)
		Y<-as.matrix(dataq.int.sig.allinfo[!is.na(match(dataq.int.sig.allinfo$time,unique(dataq.int.sig.allinfo$time)[tiniplsr:tinitial])),
			(which(colnames(dataq.int.sig.allinfo)=="experiment")+1):ncol(dataq.int.sig.allinfo)])
		if (dim(Y)[2]==1){
			colnames(Y) <- colnames(dataq[length(colnames(dataq))])
		} else {
			colnames(Y)<-colnames(dataq.int.sig.allinfo)[(which(colnames(dataq.int.sig.allinfo)=="experiment")+1):ncol(dataq.int.sig.allinfo)]
		}
		X<-as.matrix(dataq.int.sig.allinfo[!is.na(match(dataq.int.sig.allinfo$time,unique(dataq.int.sig.allinfo$time)[tiniplsr:tinitial])),
			1:nrow(Significant)])
		## Partial Least Square Regression for treatments
		if (clstr==1){
			library(pls, lib.loc=lib.loc)
		} else {
			library(pls)
		}
		par(mar=c(9,9,5,5))
		NET_TDPLSR<-NULL
		for (ii in 1:ncol(Y)){
			data.plsr <- plsr(Y[,ii] ~ X,validation = "LOO",ncomp=1,scale=TRUE)
			if (clstr==0){
				par(ask=TRUE)		
			}
			pdf(paste(folder.output,"/RMSEP",colnames(Y)[ii],".pdf",sep=""),paper="a4r")
			plot(RMSEP(data.plsr), legendpos = "topright",cex.axis=3, cex.lab=3,cex.main=3,
				main=colnames(Y)[ii])#This plots the estimated RMSEPs as functions of the number of components	
			dev.off()
			if (clstr==0){
				par(ask=TRUE)		
			}
			pdf(paste(folder.output,"/PLSR",colnames(Y)[ii],".pdf",sep=""),paper="a4r")
			plot(data.plsr, asp = 1, line = TRUE,cex.axis=3, cex.lab=3,cex.main=3,cex=2,main=colnames(Y)[ii])
			dev.off()
			if (clstr==0){
				par(ask=TRUE)		
			}
			pdf(paste(folder.output,"/Scores",colnames(Y)[ii],".pdf",sep=""),paper="a4r")
			plot(data.plsr, plottype = "scores", comps = 1,cex.axis=3, cex.lab=3,cex.main=3,cex=2,main=colnames(Y)[ii])
			dev.off()
			if (clstr==0){
				par(ask=TRUE)		
			}
			par(mar=c(10,10,5,5))
			if (clstr==0){
				par(ask=TRUE)		
			}
			pdf(paste(folder.output,"/Loadings",colnames(Y)[ii],".pdf",sep=""),paper="a4r")
			plot(data.plsr, "loadings", comps = 1, legendpos = "bottomright",
				xaxt="n",cex.axis=1.5, cex.lab=3,cex.main=3,cex=2,xlab="",ylim=c(-0.6,0.6))
			axis(1,at=1:ncol(X),labels=colnames(X),las=2,cex.axis=1.5, cex.lab=3,cex.main=3,cex=2)
			abline(h = 0)
			abline(h = PLSRCutoffInput,col="red")
			abline(h = -PLSRCutoffInput,col="red")
			dev.off()
			loadings(data.plsr)
			loadings<-loadings(data.plsr)[order(abs(loadings(data.plsr)),decreasing=TRUE),]
			loadings<-loadings[which(abs(loadings)>=PLSRCutoffInput)]
			if (length(loadings)>0){
				interaction<-matrix(0,1,length(loadings))
				for (jj in 1:length(loadings)){
						if (loadings[jj]>0) {
							interaction[jj]<-1 
						}else{
							interaction[jj]<-(-1)
						}
				}
				interaction=t(interaction)
				tmp<-cbind(rep(colnames(Y)[ii],length(loadings)),interaction,names(loadings))
				colnames(tmp)<-c("FromTF","Interaction","ToTF")
				NET_TDPLSR<-data.frame(rbind(NET_TDPLSR,tmp))
			}
		}	
		# For each TF
		posTANC<-which((colnames(dataq)==name_no_dna)|(colnames(dataq)==name_TA))
		dataq.sel<-dataq[,-posTANC]
		colnames(dataq.sel)<-colnames(dataq)[-posTANC]
		dataq.sel_names<-colnames(dataq.sel)[c(1:(which(colnames(dataq.sel)=="repeats")-1),(which(colnames(dataq.sel)=="experiment")+1):ncol(dataq.sel))]
		dataq.sel<-data.frame(dataq.sel[,1:(which(colnames(dataq.sel)=="repeats")-1)],dataq.sel[,(which(colnames(dataq.sel)=="experiment")+1):ncol(dataq.sel)]*dataq.sel$time)
		colnames(dataq.sel)<-dataq.sel_names
		dataInt$dataq.sel<-dataq.sel
		if (PLSRType=="PLSR"){
			x<-as.matrix(dataq.sel)	
			if (tf_ids[ntf]!=name_TA && tf_ids[ntf]!=name_no_dna){
				y<-as.matrix(dataq.sel[,1:(which(colnames(dataq.sel)==tf_ids[ntf]))])
			}else{
				y<-as.matrix(dataq.sel[,1:(which(colnames(dataq.sel)==tf_ids[ntf-1]))])
			}
		} else {
			x<-as.matrix(dataq.sel[which(dataq$time!=ntime),])	
			if (tf_ids[ntf]!=name_TA && tf_ids[ntf]!=name_no_dna){
				y<-as.matrix(dataq.sel[,1:(which(colnames(dataq.sel)==tf_ids[ntf]))])
			}else{
				y<-as.matrix(dataq.sel[,1:(which(colnames(dataq.sel)==tf_ids[ntf-1]))])
			}
			y<-y[which(dataq$time!=ntime),]-y[which(dataq$time!=1),]
		}
		par(mar=c(5,10,10,5))
		for (ii in 1:ncol(y)){
			data.plsr <- plsr(y[,ii] ~ x,validation = "LOO",ncomp=1,scale=TRUE)
			loadings<-loadings(data.plsr)[order(abs(loadings(data.plsr)),decreasing=TRUE),]
			loadings<-loadings[which(abs(loadings)>=PLSRCutoffTF)]
			if (length(loadings)>0){
				interaction<-matrix(0,1,length(loadings))
				for (jj in 1:length(loadings)){
					if (loadings[jj]>0) interaction[jj]<-1 else interaction[jj]<-(-1)
				}
				interaction=t(interaction)
				tmp<-cbind(names(loadings), interaction,rep(colnames(y)[ii],length(loadings)))
				colnames(tmp)<-c("FromTF","Interaction","ToTF")
				NET_TDPLSR<-data.frame(rbind(NET_TDPLSR,tmp))
			}
		}
		if (PLSRType=="PLSR"){
			write.table(NET_TDPLSR, file =  paste(folder.output,"/NET_PLSR.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
				row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		} else {
			write.table(NET_TDPLSR, file =  paste(folder.output,"/NET_TDPLSR.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
				row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		}
		PLSRRes<-list(NET_TDPLSR=NET_TDPLSR,dataInt=dataInt)
		return(PLSRRes)
	}
	
#=====================================================================================	

#     MIInference
#=====================================================================================
#=====================================================================================	
	
	MIInference<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		if (clstr==1){
			library(infotheo,lib.loc=lib.loc)
			library(minet,lib.loc=lib.loc)
		} else {
			library(infotheo)
			library(minet)
		}
		# Build the mutual information matrix for time
		dataq.int.sig.m<-dataq.int.sig
		dataq.int.sig.m[dataq.int$time==1,(which(colnames(dataq.int.sig.m)==rownames(Significant)[nrow(Significant)])+1):ncol(dataq.int.sig.m)]<-0
		tinitmi<-ifelse (0 %in% timevector,2,1)
		dataq.int.sig.m.short<-dataq.int.sig.m[!is.na(match(dataq.int$time,unique(dataq.int$time)[tinitmi:tinitial])),
			c(1:nrow(Significant),(which(colnames(dataq.int.sig.m)==rownames(Significant)[nrow(Significant)])+1):ncol(dataq.int.sig.m))]
		dataq.int.sig.m.all<-dataq.int.sig.m[,c(1:nrow(Significant),
			(which(colnames(dataq.int.sig.m)==rownames(Significant)[nrow(Significant)])+1):ncol(dataq.int.sig.m))]
		dataInt$dataq.int.sig.m.short<-dataq.int.sig.m.short
		dataInt$dataq.int.sig.m.all<-dataq.int.sig.m.all
		# Mutual Information matrix
		ncolMIM<-ncol(dataq.int.sig.m.short)
		nsig<-nrow(Significant)
		MIMdim<-ncolMIM-nsig	
		NET_MIM.all<-NULL
		NET_MIM.short<-NULL
		for (kk in 1:MIMdim){
			indexMIM<-nsig+kk
			datatempshort<-dataq.int.sig.m.short[,c(1:nsig,indexMIM)]
			datatempall<-dataq.int.sig.m.all[,c(1:nsig,indexMIM)]
			if (MIMdim>1){
				mimrowshort<-apply(dataq.int.sig.m.short,1,function(x){
					if (all(x[(nsig+1):ncolMIM]==0)||x[indexMIM]==1){
						r=TRUE
					} else {
						r=FALSE
					}
					return(r)
				})
				mimrowsall<-apply(dataq.int.sig.m.all,1,function(x){
					if (all(x[(nsig+1):ncolMIM]==0)||x[indexMIM]==1){
						r=TRUE
					} else {
						r=FALSE
					}
					return(r)
				})
				MIM<-build.mim(datatempshort[mimrowshort,],estimator = "mi.sg",dis="equalfreq")
				MIM.all<-build.mim(datatempall[mimrowsall,],estimator = "mi.sg",dis="equalfreq")
			} else {
				MIM<-build.mim(datatempshort,estimator = "mi.sg",dis="equalfreq")
				MIM.all<-build.mim(datatempall,estimator = "mi.sg",dis="equalfreq")
			}
			aracne.res<-aracne(MIM, eps=0 )
			clr.res<-clr(MIM)
			mrnet.res<-mrnet(MIM)
			aracne.all.res<-aracne(MIM.all, eps=0 )
			clr.all.res<-clr(MIM.all)
			mrnet.all.res<-mrnet(MIM.all)
			netcon<-function(res,datatemp){
				net.con<-NULL
				for (jj in 1:nrow(Significant)){
					if (res[jj,(nsig+1)]>0){
						interaction=1*sign(cor(datatemp[,(nsig+1)],datatemp[,jj]))
						net.con<-data.frame(rbind(net.con,cbind(FromTF=colnames(res)[(nsig+1)],Interaction=interaction,ToTF=rownames(res)[jj])))
					}
				}
				return(net.con)
			}
			netcon.aracne<-netcon(aracne.res,datatempshort)
			netcon.clr<-netcon(clr.res,datatempshort)
			netcon.mrnet<-netcon(mrnet.res,datatempshort)
			netcon.all.aracne<-netcon(aracne.all.res,datatempall)
			netcon.all.clr<-netcon(clr.all.res,datatempall)
			netcon.all.mrnet<-netcon(mrnet.all.res,datatempall)
			NET_MIM.short<-rbind(NET_MIM.short,cbind(netcon.aracne,METHOD=rep("ARACNE",nrow(netcon.aracne))),cbind(netcon.clr,METHOD=rep("CLR",nrow(netcon.clr))),
				cbind(netcon.mrnet,METHOD=rep("MRNET",nrow(netcon.mrnet))))
			NET_MIM.all<-rbind(NET_MIM.all,cbind(netcon.all.aracne,METHOD=rep("ARACNE",nrow(netcon.all.aracne))),cbind(netcon.all.clr,METHOD=rep("CLR",nrow(netcon.all.clr))),
				cbind(netcon.all.mrnet,METHOD=rep("MRNET",nrow(netcon.all.mrnet))))
		}
		if (MIMdata=="STAT"){
			NET_MIM<-NET_MIM.all
		}
		if (MIMdata=="DYN"){
			NET_MIM<-NET_MIM.short
		}
		if (MIMdata=="BOTH"){
			NET_MIM<-rbind(NET_MIM.short,NET_MIM.all)
		}	
		write.table(NET_MIM, file = paste(folder.output,"/NET_MIM.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		MIRes<-list(NET_MIM=NET_MIM,dataInt=dataInt)
		return(MIRes)
	}
	
#=====================================================================================	


#    BANJOInference
#=====================================================================================
#=====================================================================================


BANJOInference<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		# Banjo Prior Knowledge File
		PKN<-read.table(PKN, header = TRUE, sep = "\t", colClasses = "character")
		PKNList<-list(PKN_Original=PKN)
		ConversionList<-read.table(ConversionList, header = TRUE, sep = "\t", colClasses = "character")
		PKNList[["ConversionList"]]<-ConversionList
		if (!is.null(PKN_INPUTS)){
			PKN_INPUT<-read.table(PKN_INPUTS, header = FALSE, sep = "\t", colClasses = "character")
			PKNList[["PKN_INPUT_Original"]]<-PKN_INPUT
			PKN_INPUT[,1]<-apply(PKN_INPUT,1,function(x){
				if(length(which(ConversionList[,"ORIGINAL"]==x[1]))==0){
					x[1]=x[1]
				} else{
					x[1]<-ConversionList[which(ConversionList[,"ORIGINAL"]==x[1]),"NAME"][1]
				}
			})
			PKN_INPUT[,3]<-apply(PKN_INPUT,1,function(x){
				if(length(which(ConversionList[,"ORIGINAL"]==x[3]))==0){
					x[3]=x[3]
				} else{
					x[3]<-ConversionList[which(ConversionList[,"ORIGINAL"]==x[3]),"NAME"][1]
				}
			})
			PKN_INPUT<-PKN_INPUT[!duplicated(PKN_INPUT),]
			PKN_INPUT["NumberDatabases"]<-rep(1,nrow(PKN_INPUT))
			PKN_INPUT["Databases"]<-rep("GENEGO",nrow(PKN_INPUT))
			colnames(PKN_INPUT)<-c("FromTF","Interaction","ToTF","NumberDatabases","Databases")
			PKN<-rbind(PKN,PKN_INPUT)
		}
		PKNList[["PKN_INPUT"]]<-PKN_INPUT
		PKN_CNO<-PKN
		write.table(PKN_CNO, file = paste(folder.output,"/PKN_ALLPRIOR.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
			col.names = TRUE, qmethod = c("escape", "double"))
		PKNList[["PKN_ALLPRIOR"]]<-PKN_CNO
		dir.create(paste(folder.output,"/Banjo",sep=""))
		test<-NULL
		for (mm in 1: nrow(PKN_CNO)){
			From<-which(ConversionList[,"SYMBOL"]==PKN_CNO[mm,1])
			To<-which(ConversionList[,"SYMBOL"]==PKN_CNO[mm,3])
			if (length(From)==0 && length(To)==0) {
				test<-rbind(test,PKN_CNO[mm,])
			}
			if (length(From)>0 && length(To)==0) {
					test<-rbind(test,
						data.frame(FromTF=ConversionList[From,"NAME"],
								Interaction=rep(PKN_CNO[mm,2],length(From)),
								ToTF=rep(PKN_CNO[mm,3],length(From)),
								NumberDatabases=rep(PKN_CNO[mm,4],length(From)),
								Databases=rep(PKN_CNO[mm,5],length(From))
						)
					)
			}
			if (length(From)==0 && length(To)>0) {
						test<-rbind(test,
							data.frame(FromTF=rep(PKN_CNO[mm,1],length(To)),
								Interaction=rep(PKN_CNO[mm,2],length(To)),
								ToTF=ConversionList[To,"NAME"],
								NumberDatabases=rep(PKN_CNO[mm,4],length(To)),
								Databases=rep(PKN_CNO[mm,5],length(To))
							)
						)
			}
			if (length(From)>0 && length(To)>0){
				for (kk in 1:length(From)){
					test<-rbind(test,
						data.frame(FromTF=rep(ConversionList[From[kk],"NAME"],length(To)),
								Interaction=rep(PKN_CNO[mm,2],length(To)),
								ToTF=ConversionList[To,"NAME"],
								NumberDatabases=rep(PKN_CNO[mm,4],length(To)),
								Databases=rep(PKN_CNO[mm,5],length(To))
						)
					)
				}
			}
		}
		test<-test[!duplicated(test),]
		PKN_CNO<-test
		pos.banjo<-NULL
		for (bb in 1:nrow(PKN_CNO)){
			if (length(which(ConversionList[match(rownames(Significant),ConversionList$ORIGINAL),][,"NAME"]==PKN_CNO[bb,1]))==1){
				if (length(which(ConversionList[match(rownames(Significant),ConversionList$ORIGINAL),][,"NAME"]==PKN_CNO[bb,3]))==1){ 
					pos.banjo<-c(pos.banjo,bb)
				}
			} else {
				if (length(which(colnames(dataq.int.sig)[(nrow(Significant)+1):ncol(dataq.int.sig)]==PKN_CNO[bb,1]))==1){
					if (length(which(ConversionList[match(rownames(Significant),ConversionList$ORIGINAL),][,"NAME"]==PKN_CNO[bb,3]))==1){ 
						pos.banjo<-c(pos.banjo,bb)
					}
				}
			}
		}
		PKN_BANJO<-PKN_CNO[pos.banjo,1:3]
		PKNList[["PKN_BANJO"]]<-PKN_BANJO
		PKNList[["PKN_CNO"]]<-PKN_CNO
		write.table(PKN_BANJO, file = paste(folder.output,"/Banjo/PKN_BANJO.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
			col.names = TRUE, qmethod = c("escape", "double"))
		# Create all the necessary files
		for (j in 1:ntreatment){
			count<-0
			for (k in 1:nexp){
				for (i in 1:nrepeats){
					count<-count+1
					temp<-dataq.int.sig[which(as.logical(dataq.int.sig.allinfo$treatment==j & dataq.int.sig.allinfo$repeats==i & dataq.int.sig.allinfo$experiment==k)),]
					assign(paste(treatment_ids[j],count,sep="_"),temp)
					if (timevector[1]==0){
						temp[1,]=0
					}
					write.table(temp, file = paste(folder.output,"/Banjo/",paste(treatment_ids[j],count,sep="_"),".txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
						col.names = FALSE, qmethod = c("escape", "double"))
				}
			}
		}
		databanjostatic<-dataq.int.sig[-which(dataq.int.sig.allinfo$time==1),]	
		write.table(databanjostatic,file = paste(folder.output,"/Banjo/static.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
			col.names = FALSE, qmethod = c("escape", "double"))
		BANJORes<-PKNList
		return(BANJORes)
	}

#=====================================================================================	

#     PKNProcessing
#=====================================================================================
#=====================================================================================


	PKNProcessing<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		PKNList.names<-names(PKNList)
		for (ii in 1:length(PKNList)){
			assign(PKNList.names[[ii]],PKNList[[ii]])
		}
		InferenceRes.names<-names(InferenceRes)
		for (ii in 1:length(InferenceRes)){
			assign(InferenceRes.names[[ii]],InferenceRes[[ii]])
		}
		sigConversionList<-match(rownames(Significant),ConversionList[,"ORIGINAL"])
		FromTF<-match(PKN_CNO[,1],ConversionList[sigConversionList,"NAME"])
		FromCues<-match(PKN_CNO[,1],colnames(dataq.int.sig)[(nrow(Significant)+1):ncol(dataq.int.sig)])
		ToTF<-match(PKN_CNO[,3],ConversionList[sigConversionList,"NAME"])
		positions<-apply(cbind(FromTF,FromCues,ToTF),1,function(x){
			if ((!is.na(x[1])||!is.na(x[2]))&& !is.na(x[3])){
				r<-TRUE
			} else {
				r<-FALSE
			}
			return(r)
			}
		)
		PKN_CNO<-PKN_CNO[positions,]
		write.table(PKN_CNO, file = paste(folder.output,"/PKN_CNO.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		PKNList[["PKN_CNO_Original"]]<-PKNList[["PKN_CNO"]]
		PKNList[["PKN_CNO_ConversionList"]]<-PKN_CNO
		# Inferred network from ARACNE, CLR and BANJO
		NET_Inferred<-rbind(
			cbind(NET_TDPLSR,METHOD=rep("TDPLSR",nrow(NET_TDPLSR))),
			NET_MIM,
			cbind(NET_BANJO,METHOD=rep("BANJO",nrow(NET_BANJO)))
		)
		ones<-rep(1,nrow(NET_Inferred))
		NET_InferredUnique<-aggregate(ones, by = as.list(NET_Inferred[,1:3]), FUN = sum) 
		NETInferredUniquePaste<-cbind(NET_InferredUnique,apply(NET_InferredUnique,1,function(x){paste(x[1],x[2],x[3],sep="")}))
		NETInferredPaste<-cbind(NET_Inferred,apply(NET_Inferred,1,function(x){paste(x[1],x[2],x[3],sep="")}))
		METHODS<-apply(NETInferredUniquePaste,1,function(x){
			idmethods<-which(NETInferredPaste[,5]==x[5])
			namemethods<-c()
			for (i in 1:length(idmethods)){
				namemethods<-paste(namemethods,NETInferredPaste[idmethods[i],4],sep=" ")
			}
			return(namemethods)
		})
		NET_Inferred<-cbind(NET_InferredUnique,METHODS)
		NET_Inferred[,1]<-apply(NET_Inferred,1,function(x){if(length(which(ConversionList[,"ORIGINAL"]==x[1]))==0){x[1]=x[1]} else{x[1]<-ConversionList[which(ConversionList[,"ORIGINAL"]==x[1]),"NAME"][1]}})
		NET_Inferred[,3]<-apply(NET_Inferred,1,function(x){if(length(which(ConversionList[,"ORIGINAL"]==x[3]))==0){x[3]=x[3]} else{x[3]<-ConversionList[which(ConversionList[,"ORIGINAL"]==x[3]),"NAME"][1]}})
		write.table(NET_Inferred, file = paste(folder.output,"/PKN_Inferred.txt",sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
			row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		NET_Inferred<-NET_Inferred[,1:3]
		InferenceRes[["PKN_Inferred"]]<-NET_Inferred
		sigConversionList<-match(rownames(Significant),ConversionList[,"ORIGINAL"])
		FromTF<-match(NET_Inferred[,1],ConversionList[sigConversionList,"NAME"])
		FromCues<-match(NET_Inferred[,1],colnames(dataq.int.sig)[(nrow(Significant)+1):ncol(dataq.int.sig)])
		ToTF<-match(NET_Inferred[,3],ConversionList[sigConversionList,"NAME"])
		positions<-apply(cbind(FromTF,FromCues,ToTF),1,function(x){
			if ((!is.na(x[1])||!is.na(x[2]))&& !is.na(x[3])){
				r<-TRUE
			} else {
				r<-FALSE
			}
			return(r)
			}
		)
		NET_Inferred<-NET_Inferred[positions,]
		InferenceRes[["NET_Inferred"]]<-NET_Inferred
		#Identify the common edges from prior knowledge and infer methods
		PKN_CNO<-data.frame(rbind(cbind(PKN_CNO[,1:3],source=rep(1,nrow(PKN_CNO))),cbind(NET_Inferred[,1:3], source=rep(2,nrow(NET_Inferred)))))
		# Remove all the duplications
		PKN_CNO<-PKN_CNO[!duplicated(PKN_CNO),]
		testPKN<-data.frame(Name=apply(PKN_CNO,1,function(x){paste(x[1],x[2],x[3],sep="")}),Source=PKN_CNO[,4])
		numbertestPKN<-tapply(testPKN[,2],testPKN[,1],length)
		testPKN<-data.frame(testPKN[match(names(numbertestPKN),testPKN[match(unique(testPKN[,1]),testPKN[,1]),][,1]),],Number=as.vector(numbertestPKN))
		typeedge<-matrix(NA,nrow=nrow(testPKN),ncol=1)
		for (i in 1:nrow(testPKN)){
			if (testPKN[i,2]==1 && testPKN[i,3]==1) {
				typeedge[i]=1				
			}
			if (testPKN[i,2]==2 && testPKN[i,3]==1) {
				typeedge[i]=2				
			}
			if (testPKN[i,3]==2) {
				typeedge[i]=3				
			}
		}
		testPKN<-cbind(testPKN,Type=typeedge)
		PKN_CNO_cytoscape<-cbind(PKN_CNO[match(testPKN[,1],apply(PKN_CNO,1,function(x){paste(x[1],x[2],x[3],sep="")})),1:3],testPKN[,c(-1,-2)])
		write.table(PKN_CNO_cytoscape, file = paste(folder.output,"/PKN_CNO_cytoscape.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
		row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		PKNList[["PKN_CNO_cytoscape"]]<-PKN_CNO_cytoscape
		PKN_CNO_all<-PKN_CNO
		#Remove self-loops
		selfloopind<-NULL
		for (tt in 1:nrow(PKN_CNO)){
			if (PKN_CNO[tt,1]==PKN_CNO[tt,3]){
				selfloopind<-c(selfloopind,tt)
			}
		}	
		PKN_CNO<-PKN_CNO[-selfloopind,1:3]
		PKNList[["PKN_CNO"]]<-PKN_CNO
		PKNSim<-list(PKNList=PKNList,InferenceRes=InferenceRes)
		return(PKNSim)
	}


#=====================================================================================	


#     preCNO
#=====================================================================================
#=====================================================================================

preCNO<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		PKNList.names<-names(PKNList)
		for (ii in 1:length(PKNList)){
			assign(PKNList.names[[ii]],PKNList[[ii]])
		}
		InferenceRes.names<-names(InferenceRes)
		for (ii in 1:length(InferenceRes)){
			assign(InferenceRes.names[[ii]],InferenceRes[[ii]])
		}
		# Create the CNOlist structure
		dir.create(paste(folder.output,"/CNO",sep=""))
		colnames(dataq)[1:ntf]<-colnames(dataq)[1:ntf]
		dataCNO<-dataq[,c(match(rownames(Significant),colnames(dataq)),(ntf+1):ncol(dataq))]
		rownames(dataCNO)<-NULL
		dataCNO<-dataCNO[-which(dataCNO["treatment"]==which(treatment_ids==name_untreated)),]
		# Boolean discretization
		# Without accounting p-value
		dataCNO<-data.frame(
			apply(dataCNO[,1:nrow(Significant)],c(1,2),function(x) {
				if (is.na(x)|| x==0) {
					x=0 
				} else{
				x<-ifelse(x>0,1,(-1))
				}
			 }),
			dataCNO[,(nrow(Significant)+1):ncol(dataCNO)]
		)
		# Accounting p-value
		nocontroltreatment<-length((which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO))
		booltrans<-matrix(0,nrow=nrow(Significant),ncol=(ncol(Significant)-1+nocontroltreatment))
		Significant.nopvalue<-Significant[,-1]
		for (iii in 1:nrow(Significant.nopvalue)){
			for (jjj in 1:nocontroltreatment){
				booltrans[iii,(2+(jjj-1)*ntime):(jjj*ntime)]<-as.numeric(
					Significant.nopvalue[iii,(1+(jjj-1)*(ntime-1)):(jjj*(ntime-1))])
			}
		}
		for (kkk in 1:nrow(Significant)){
			for (lll in 1:1:nocontroltreatment){
				for (mmm in 1:ntime){
					if(booltrans[kkk,mmm+(lll-1)*ntime]==0){
						dataCNO[which((dataCNO[,"time"]==mmm)&
							(dataCNO[,(which(colnames(dataCNO)=="experiment")+lll)]==1)),kkk]<-0
					}
				}
			}
		}		
		colnames(dataCNO)[1:nrow(Significant)]<-rownames(Significant)
		colnames(dataCNO)[1:nrow(Significant)]<-ConversionList[na.omit(match(colnames(dataCNO),ConversionList[,"ORIGINAL"])),"NAME"]
		if (bootstrapping){
            dataCNOboot<-dataCNO
			dataCNOboot[,"repeats"]<-dataCNO[,"repeats"]+nrepeats*(dataCNO[,"experiment"]-1)
			ntrid<-max(dataCNOboot[,"treatment"])
			nSample<-max(dataCNOboot[,"repeats"])
			for (i in 1:nrow(Significant)){
        		for (bb in 1:ntrid){
				orderData<-sample(nSample,replace=T)
					for (aa in 1:nSample){
						dataCNOboot[which(dataCNOboot[,"repeats"]==aa & dataCNOboot[,"treatment"]==bb),i]<-
							dataCNOboot[which(dataCNOboot[,"repeats"]==orderData[aa] & dataCNOboot[,"treatment"]==bb),i]
					}
				}
			}
            dataCNO<-dataCNOboot
		}
		if (permutation){
			if(all(dataCNO[which(dataCNO$time==1),1:nrow(Significant)]==0)){
				dataCNOtmp<-dataCNO[which(dataCNO$time!=1),]
				dataCNOper<-dataCNOtmp
				nSample<-nrow(dataCNOper)
				for (i in 1:nrow(Significant)){
					orderData<-sample(nSample,replace=F)
                	dataCNOper[,i]<-dataCNOtmp[orderData,i]
            	}
				dataCNO[which(dataCNO$time!=1),]<-dataCNOper
			} else {
				dataCNOtmp<-dataCNO
            	dataCNOper<-dataCNO
           		nSample<-nrow(dataCNOper)
            	for (i in 1:nrow(Significant)){
					orderData<-sample(nSample,replace=F)
                	dataCNOper[,i]<-dataCNOtmp[orderData,i]
            	}
            	dataCNO<-dataCNOper
			}
        }
		if (DAPresent){
			namesCues<-c(colnames(dataCNO)[(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],"DA","dummyInh")
			namesStimuli<-c(colnames(dataCNO)[(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],"DA")
		} else {
			namesCues<-c(colnames(dataCNO)[(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],"dummyInh")
			namesStimuli<-colnames(dataCNO)[(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)]
		}
		namesInhibitors<-"dummyInh"
		namesSignals<-colnames(dataCNO)[1:nrow(Significant)]
		treatSignals<-unique(dataCNO[,"treatment"])
		ntreatSignals<-length(treatSignals)
		averageSignals<-matrix(0,nrow=(ntreatSignals*ntime), ncol=nrow(Significant))
		for (aa in 1:ntreatSignals){
			for (bb in 1:nrow(Significant)){
				averageSignals[(1+(aa-1)*ntime):(aa*ntime),bb]<-tapply(
					dataCNO[which(dataCNO[,"treatment"]==treatSignals[aa]),bb],
					dataCNO[which(dataCNO[,"treatment"]==treatSignals[aa]),"time"],
					function(x){round(median(x),0)}
				)
			}
		}
		namesActiveCtrlSignals<-namesSignals[apply(averageSignals,2,function(x){any(x==(-1))})]
		timeSignals<-as.vector(timevector)
		preCNO<-list(dataCNO=dataCNO,averageSignals=averageSignals,namesActiveCtrlSignals=namesActiveCtrlSignals,
			timeSignals=timeSignals)
		if (DAPresent){
			valueCues<-as.matrix(cbind(dataCNO[which(dataCNO$time==1),(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],
				"DA"=matrix(1,nrow=nrow(dataCNO[which(dataCNO$time==1),]),ncol=1),
				"dummyInh"=matrix(0,nrow=nrow(dataCNO[which(dataCNO$time==1),]),ncol=1)))
			valueStimuli<-as.matrix(cbind(dataCNO[which(dataCNO$time==1),(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],
				"DA"=matrix(1,nrow=nrow(dataCNO[which(dataCNO$time==1),]),ncol=1)))
		} else {
			valueCues<-as.matrix(cbind(dataCNO[which(dataCNO$time==1),(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)],
				"dummyInh"=matrix(0,nrow=nrow(dataCNO[which(dataCNO$time==1),]),ncol=1)))
			valueStimuli<-as.matrix(dataCNO[which(dataCNO$time==1),(which(colnames(dataCNO)=="experiment")+1):ncol(dataCNO)])
		}
		valueInhibitors<-matrix(0,nrow = length(which(dataCNO$time==1)),ncol=1)
		colnames(valueInhibitors)<-"dummyInh"
		value.Signals <- list(t0 = matrix(data = 0, nrow = length(which(dataCNO$time==1)),ncol = nrow(Significant)))
		value.Signals$t0<- dataCNO[which(dataCNO$time==1),1:nrow(Significant)]
    	for (i in 2:length(timeSignals)) {
			value.Signals[[i]] <- dataCNO[which(dataCNO$time==i),1:nrow(Significant)]
    	}
		for (aa in 2:length(timeSignals)){
			if(aa==2){
				valueSignals=list(t0=value.Signals$t0,value.Signals[[2]])
			} else {
				valueSignals=list(t0=value.Signals[[aa-1]],value.Signals[[aa]])
			}
			assign(paste("CNOlist",(aa-1),sep="_"), 
				list(namesCues = namesCues, namesStimuli = namesStimuli, 
					namesInhibitors = namesInhibitors, namesSignals = namesSignals,        		
					timeSignals = timeSignals[(aa-1):aa], valueCues = valueCues, valueInhibitors = valueInhibitors, 
					valueStimuli = valueStimuli, valueSignals=valueSignals))		
			tempname<-paste("CNOlist",(aa-1),sep="_")
			preCNO[[tempname]]<-list(namesCues = namesCues, namesStimuli = namesStimuli, 
					namesInhibitors = namesInhibitors, namesSignals = namesSignals,        		
					timeSignals = timeSignals[(aa-1):aa], valueCues = valueCues, valueInhibitors = valueInhibitors, 
					valueStimuli = valueStimuli, valueSignals=valueSignals)
			CNOlist_all<-list(
			namesCues = namesCues, 
			namesStimuli = namesStimuli, 
        		namesInhibitors = namesInhibitors, 
			namesSignals = namesSignals, 
        		timeSignals = timeSignals, 
			valueCues = valueCues, 
			valueInhibitors = valueInhibitors, 
        		valueStimuli = valueStimuli, 
			valueSignals=value.Signals
			)	
			preCNO[["CNOlist_all"]]<-CNOlist_all
		}
		#Add dummy inhibitor
		PKN_CNO<-rbind(PKN_CNO,cbind(FromTF="dummyInh",Interaction=-1,
			ToTF=ConversionList[which(ConversionList[,"ORIGINAL"]==rownames(Significant)[1]),"NAME"]))
		# Add de-activation TF interactions
		if (DAPresent){
			nsig<-nrow(Significant)
			fromTF<-rep("DA",nsig)
			interactionDA<-rep(-1,nsig)
			toTF<-ConversionList[match(rownames(Significant),ConversionList[,"ORIGINAL"]),"NAME"]
			PKN_DA<-cbind(FromTF=fromTF,Interaction=interactionDA,ToTF=toTF)
			PKN_CNO<-rbind(PKN_CNO,PKN_DA)
		}
		write.table(PKN_CNO, file = paste(folder.output,"/PKN_CNO_SimplifiedandInferred.txt", sep=""),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
		row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
		PKNList[["PKN_CNO_SimplifiedandInferred"]]<-PKN_CNO
		PKNList[["PKN_CNO"]]<-PKN_CNO
		preCNO[["PKNList"]]<-PKNList
		return(preCNO)
	}

#=====================================================================================	


#     wrappingCNO
#=====================================================================================
#=====================================================================================
	wrappingCNO<-function(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes,
		dataCNOpre=dataCNOpre){
		inputParameters.names<-names(inputParameters)
		for (ii in 1:length(inputParameters)){
			assign(inputParameters.names[[ii]],inputParameters[[ii]])
		}
		dataAnalysis.names<-names(dataAnalysis)
		for (ii in 1:length(dataAnalysis)){
			assign(dataAnalysis.names[[ii]],dataAnalysis[[ii]])
		}
		dataInt.names<-names(dataInt)
		for (ii in 1:length(dataInt)){
			assign(dataInt.names[[ii]],dataInt[[ii]])
		}
		PKNList.names<-names(PKNList)
		for (ii in 1:length(PKNList)){
			assign(PKNList.names[[ii]],PKNList[[ii]])
		}
		InferenceRes.names<-names(InferenceRes)
		for (ii in 1:length(InferenceRes)){
			assign(InferenceRes.names[[ii]],InferenceRes[[ii]])
		}
		dataCNOpre.names<-names(dataCNOpre)
		for (ii in 1:length(dataCNOpre)){
			assign(dataCNOpre.names[[ii]],dataCNOpre[[ii]])
		}
		# Set up CN0
		if (clstr==1){
			library(CellNOptR,lib.loc=lib.loc)
		} else {
			library(CellNOptR)
		}
		Model<-readSIF(sifFile = paste(folder.output,"/PKN_CNO_SimplifiedandInferred.txt",sep=""))
		dataCNOOpt<-list(Model=Model)
		# Plot CNOlist
		dir.create(paste(folder.output,"/CNO",sep=""))
		setwd(paste(folder.output,"/CNO",sep=""))
		for (bb in 1:(length(timeSignals)-1)){
			plotCNOlistPDF(get(paste("CNOlist_",bb,sep="")),filename=paste("CNOlist_",bb,"plot.pdf",sep=""))
			if (clstr==0){
				plotCNOlist(get(paste("CNOlist_",bb,sep="")))
			}
		}	
		plotCNOlistPDF(CNOlist_all,filename="CNOlist_allplot.pdf")
		if (clstr==0){
			plotCNOlist(CNOlist_all)
		}
		# Compress the model, expand the gates
		valueSignalsIniSim<-as.matrix(CNOlist_all$valueSignals[[1]])
		uniquecond<-unique(apply(CNOlist_all$valueCues,1,function(x){paste(x,collapse="")}))
		for (ii in 1:length(uniquecond)){
			tmppos<-which(uniquecond[ii]==apply(CNOlist_all$valueCues,1,function(x){paste(x,collapse="")}))
			tmpmean<-matrix(apply(valueSignalsIniSim[tmppos,],2,function(x){round(median(x,na.rm=TRUE))}),nrow=1)
			valueSignalsIniSim[tmppos,]<- tmpmean[rep(nrow(tmpmean),length(tmppos)),]
		}
		dataCNOOpt[["valueSignalsIniSim"]]<-valueSignalsIniSim
		finalBString<-c()
		for (bb in 1:(length(timeSignals)-1)){
			CNOlist<-get(paste("CNOlist_",bb,sep=""))
				NPcut<-findNP(CNOlist,Model,verbose=TRUE)
			UpdateValues<-cutNP(CNOlist,valueSignalsIniSim,NPcut)	
			CNOlist<-UpdateValues$CNOlist
			valueSignalsIniSim<-UpdateValues$valueSignalsIniSim
			indices<-indexFinder(CNOlist,Model, verbose = TRUE)
			NCNOindices<-findNONC(Model,indices, verbose = TRUE)
			NCNOcut<-cutNONC2(Model, NCNOindices)
			indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
			NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
			indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
			ignoreList<-which(NCNOcutComp$namesSpecies=="dummyInh")
			if (expandModel){
					allGates<-ifelse(any(rowSums(CNOlist$valueCues,na.rm=TRUE)>1),TRUE,FALSE)
				NCNOcutCompExp<-expandGates3(NCNOcutComp,ignoreList=ignoreList,allGates=allGates,List=CNOlist)
				indicesNCNOcutCompExp<-indexFinder(CNOlist,NCNOcutCompExp)
			} 
			if (DAPresent && !expandModel){
				NCNOcutCompExp<-expandDAGates(Model=NCNOcutComp,List=CNOlist)
				indicesNCNOcutCompExp<-indexFinder(CNOlist,NCNOcutCompExp)
			}
			if(!DAPresent && !expandModel){
				NCNOcutCompExp<-NCNOcutComp
				indicesNCNOcutCompExp<-indicesNCNOcutComp
					NCNOcutCompExp[["SplitANDs"]]<-list(initialReac=c("split1","split2"))
				NCNOcutCompExp[["newANDs"]]<-NULL
			}
			dataCNOOpt[["NCNOcutCompExp"]]<-NCNOcutCompExp
			CtrlSignals<-match(namesActiveCtrlSignals,NCNOcutCompExp$namesSpecies)
			dataCNOOpt[["CtrlSignals"]]<-CtrlSignals
			indicesCtrlSignals<-CtrlSignals[!is.na(CtrlSignals)]
			resECNOlist<-residualError(CNOlist)
			Fields4Sim<-prep4sim(NCNOcutCompExp)
			if (bb==1){
				if(is.null(initStringUser)){
					initBstring<-rep(1,length(NCNOcutCompExp$reacID))
					if (randedgeper=="RAND"){
						randedgeper=runif(1); randedgeper
					}
					initBstring[sample.int(length(NCNOcutCompExp$reacID),replace=FALSE)[1:floor(length(NCNOcutCompExp$reacID)*randedgeper)]]<-0
				} else{
					initStringUser<-initStringUser[,
						match(NCNOcutCompExp$reacID,names(initStringUser))
					]
					initBstring<-as.numeric(initStringUser[1,])
				}	
				nSp <- dim(NCNOcutCompExp$interMat)[1]
    			nReacs <- dim(NCNOcutCompExp$interMat)[2]
                if (is.null(dim(CNOlist$valueStimuli))){
					nCond<-length(CNOlist$valueStimuli)
				} else {
    				nCond <- dim(CNOlist$valueStimuli)[1]
				}
				SimResPrior <- matrix(data = NA, nrow = nCond, ncol = nSp)
				colnames(SimResPrior) <- NCNOcutCompExp$namesSpecies
    			SimResPrior[, indicesNCNOcutCompExp$stimulated] <- CNOlist$valueStimuli
    			SimResPrior[, indicesNCNOcutCompExp$inhibited] <- CNOlist$valueInhibitors
				SimResPrior[,indicesNCNOcutCompExp$signals]<-valueSignalsIniSim
				timeinit<-TRUE
				ChangeSp<-matrix(0,nrow=nCond,ncol=nSp)
				ChangeSp[,indicesNCNOcutCompExp$stimulated]<-1
				ChangeSp[,indicesNCNOcutCompExp$inhibited]<-1 
				SimResultsdata_all<-list(valueSignalsIniSim)               
			} else {
				if(is.null(initStringUser)){
					initBstring<-T1opt$bString
				} else {
					initStringUser<-initStringUser[,
						match(NCNOcutCompExp$reacID,names(initStringUser))
					]
					initBstring<-as.numeric(initStringUser[bb,])
				}			
				timeinit<-FALSE		
				StimuliFac=incPolStim*StimuliFac 
			}
			# Running genetic algorithm
            if (!expandModel){
            	modinterMat<-rbind(NCNOcutCompExp$interMat,
				seq(1,ncol(NCNOcutCompExp$interMat),by=1)
            	)
				OppEdges<-apply(modinterMat,2,function(x){
                  	interMat<-x[-length(x)]
					ipos<-x[length(x)]
					posid<-apply(NCNOcutCompExp$interMat[,-ipos],2,function(y){
						identical(interMat,y)
					})
                  	r<-which(posid==TRUE)
                  	if (length(r)==1){
						t<-r
                  		t<-ifelse(t<ipos,NA,t+1)
					} else {
						t<-NA
					}
					return(t)
				})
				OppEdges<-cbind(seq(1:length(OppEdges)),OppEdges)
					OppEdgesList<-OppEdges[!is.na(OppEdges[,2]),] 
				if (DAPresent){
					DAOppEdgesList<-c()
	     			nstim<-length(indicesNCNOcutCompExp$stimulated)-1
           			for (ii in 1: nstim){
						DAOppEdges<-apply(modinterMat,2,function(x){
                           	notMatx<-NCNOcutCompExp$notMat[
							indicesNCNOcutCompExp$stimulated[ii],x[length(x)]]
							interMatx<-x[indicesNCNOcutCompExp$stimulated[ii]]
							if (interMatx==(-1)&& notMatx==0){
                        		z<-NCNOcutCompExp$interMat[,x[length(x)]]
								z[indicesNCNOcutCompExp$stimulated[nstim+1]]<-(-1)
								posid<-apply(NCNOcutCompExp$interMat,2,function(y){
									identical(y,z)
								})
								r<-which(posid==TRUE)
								r<-ifelse(r==x[length(x)],NA,r)
							} else {
								r<-NA
							}
							return(r)
						})
						DAOppEdges<-cbind(seq(1:length(DAOppEdges)),DAOppEdges)
							DAOppEdges<-DAOppEdges[!is.na(DAOppEdges[,2]),] 
						DAOppEdgesList<-rbind(DAOppEdgesList,DAOppEdges)
					}
					OppEdgesList<-rbind(OppEdgesList,DAOppEdgesList)
				}
			} else {
				OppEdgesList<-NULL
			}
			if (DAPresent){
				DAEdgesList<-c()
	     		nstim<-length(indicesNCNOcutCompExp$stimulated)-1
           		for (ii in 1: nstim){
           			modnotMat<-rbind(NCNOcutCompExp$notMat,
						seq(1,ncol(NCNOcutCompExp$notMat),by=1)
            		)
					DAEdges<-apply(modnotMat,2,function(x){
						if (x[indicesNCNOcutCompExp$stimulated[ii]]==1){
							z<-NCNOcutCompExp$interMat[,x[length(x)]]
							z[indicesNCNOcutCompExp$stimulated[nstim+1]]<-(-1)
							posid<-apply(NCNOcutCompExp$interMat,2,function(y){
								identical(y,z)
							})
							r<-which(posid==TRUE)
							r<-ifelse(r==x[length(x)],NA,r)
						} else {
							r<-NA
						}
						return(r)
					})
					DAEdgesList<-c(DAEdgesList,DAEdges)
				}	
				DAEdgesList<-cbind(rep(seq(1,ncol(NCNOcutCompExp$interMat),by=1),nstim),DAEdgesList)
				DAEdgesList<-DAEdgesList[!is.na(DAEdgesList[,2]),]
			} else {
				DAEdgesList<-NULL
			}
			T1opt<-gaBinaryT3(
				CNOlist = CNOlist,
				Model =NCNOcutCompExp,
				SimList = Fields4Sim,
				indexList =indicesNCNOcutCompExp,
				initBstring = initBstring,
				verbose = TRUE,
				MaxTime = MaxTime,
				maxGens = maxGens,
				RelTol = 0.01,
				PopSize = PopSize, 
				StallGenMax=StallGenMax,
                elitism=elitism,
                Pmutation=Pmutation,
                SelPress=SelPress,
                sizeFac= sizeFac,
				StimuliFac=StimuliFac,
				ResPrior=SimResPrior,
				time0=timeinit,
				outputChangeSp=ChangeSp,
				indCtrlSig=indicesCtrlSignals,
				DAFac=DAFac,
				OppEdgesList=OppEdgesList,
				DAEdgesList=DAEdgesList
			)
			T1opt$results<-T1opt$Results
			if (clstr==1){
				par (ask=FALSE)
			}
			cutAndPlotResultsT3(
				Model=NCNOcutCompExp,
				bString=T1opt$bString,
				SimList=Fields4Sim,
				CNOlist=CNOlist,
				indexList=indicesNCNOcutCompExp,
				plotPDF=TRUE,
				show=FALSE,
				ResPrior=SimResPrior,
				time0=timeinit,
				outputChangeSp=ChangeSp,
				indCtrlSig=indicesCtrlSignals	
			)
			if (clstr==0){
				par(ask=TRUE)
			} else {
				par(ask=FALSE)
			}
			plotFit(optRes=T1opt,filename="evolFit.pdf")
			dev.off()
			write.table(T1opt$Results,file="EvolutionFit.txt",append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
				col.names = TRUE, qmethod = c("escape", "double"))
			writeScaffold(
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				modelOriginal=Model,
				CNOlist=CNOlist
			)
			getScaffold<-getFromNamespace("getSifInfo", "CellNOptR")
			Scaffold<-getScaffold(
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				modelOriginal=Model
			)
			Scaffold.present<-Scaffold$sifFile[which(Scaffold$sifFile[,4]==1),]
			if (!is.null(dim(Scaffold.present)) && any(Scaffold.present[,1]=="dummyInh")){
				Scaffold.present<-Scaffold.present[-which(Scaffold.present[,1]=="dummyInh"),]
			}
			if (!is.null(dim(Scaffold.present)) && dim(Scaffold.present)[1]==0){
				tmpscaffold<-matrix(0,nrow=1,ncol=3)
				colnames(tmpscaffold)<-c("reacInput","inpSign","reacOutput")
				assign(paste("Scaffold.present",bb,sep=""),tmpscaffold)
			}
            if (!is.null(dim(Scaffold.present)) && dim(Scaffold.present)[1]>0){
				assign(paste("Scaffold.present",bb,sep=""),Scaffold.present[,1:3])
			}
            if (is.null(dim(Scaffold.present))){
				assign(paste("Scaffold.present",bb,sep=""),Scaffold.present[1:3])
			}
			writeNetwork(
				modelOriginal=Model,
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=CNOlist
			)
			getNetwork<-getFromNamespace("getNetworkInfo", "CellNOptR")
			Network<-getNetwork(	
				modelOriginal=Model,
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=CNOlist,
				verbose=FALSE
			)
			SimT1Results<-simulateT3(
				CNOlist=CNOlist,
				Model=NCNOcutCompExp,
				bStringT1=T1opt$bString,
				SimList=Fields4Sim,
				indexList=indicesNCNOcutCompExp,
				ResPrior=SimResPrior,
				time0=timeinit,
				outputChangeSp=ChangeSp,
				indCtrlSig=indicesCtrlSignals	
			)
			SimResultdata<-as.matrix(SimT1Results$newInput[,indicesNCNOcutCompExp$signals])
			Expdata<-as.matrix(CNOlist$valueSignals[[2]])
			write.table(SimResultdata,file="SimResultdata.txt",append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
				col.names = TRUE, qmethod = c("escape", "double"))
			write.table(Expdata,file="Expdata.txt",append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
				col.names = TRUE, qmethod = c("escape", "double"))
			SimResultsdata_all[[bb+1]]<-SimResultdata
			namesFiles<-list(
				dataPlot=paste("CNOlist_",bb,"plot.pdf",sep=""),
				evolFit1="evolFit.pdf",
				evolFitT2=NA,
				evolFittable="EvolutionFit.txt",
				simResults1="NCNOcutCompExpSimResultsT1.pdf",
				simResults2=NA,
				SimResultsdata="SimResultdata.txt",
				SimResultsdata="Expdata.txt",
				scaffold="Scaffold.sif",
				scaffoldDot="Scaffold.dot",
				tscaffold="TimesScaffold.EA",
				wscaffold="weightsScaffold.EA",
				PKN="PKN.sif",
				PKNdot="PKN.dot",
				wPKN="TimesPKN.EA",
				nPKN="nodesPKN.NA"
			)
			writeReport(
				modelOriginal=Model,
				modelOpt=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=CNOlist,
				directory=paste("ReportFromT",timevector[bb],"ToT",timevector[bb+1],sep=""),
				namesFiles=namesFiles,
				namesData=list(
					CNOlist=CNOlist,
					model=Model),
				resE=resECNOlist
			)
            ChangeSpPrev<-ChangeSp
			ChangeSp<-SimT1Results$newInput-SimResPrior
			ChangeSp<-apply(ChangeSp,c(1,2),function(x){ifelse(x!=0,1,0)})
			ChangeSp[,indicesNCNOcutCompExp$stimulated]<-1
			ChangeSp[,indicesNCNOcutCompExp$inhibited]<-1 
			ChangeSp[which(ChangeSpPrev==1)]<-1
			SimResPrior<-SimT1Results$newInput
			finalBString<-rbind(finalBString,T1opt$bString)
		}
		colnames(finalBString)<-NCNOcutCompExp$reacID
		dataCNOOpt[["finalBString"]]<-finalBString
		Scaffold.present.allcond<-NULL
		for (cc in 1:(length(timeSignals)-1)){
			scaffoldcc<-get(paste("Scaffold.present",cc,sep=""))
            if(scaffoldcc[1]!=0){
				Scaffold.present.allcond<-rbind(Scaffold.present.allcond,
					scaffoldcc)
			}
		}
		Scaffold.present.allcond<-Scaffold.present.allcond[!duplicated(Scaffold.present.allcond),]
		rownames(Scaffold.present.allcond)<-NULL
		scaffoldfinal<-data.frame(Scaffold.present.allcond)
		scaffoldfinal$Coder<- "final"
		colnames(scaffoldfinal)<-c("FromTF", "Interaction", "ToTF",
			"Coder")
		colnames.scaffoldfinal<-NULL
		for (cc in 1:(length(timeSignals)-1)){
			scaffoldcc<-get(paste("Scaffold.present",cc,sep=""))
            if(scaffoldcc[1]==0){
				condcc<-rep(FALSE,nrow(scaffoldfinal))
			} else {
                if (is.null(dim(scaffoldcc))){
						scaffoldcc <- c(FromTF=scaffoldcc[1],
							Interaction=scaffoldcc[2], ToTF=scaffoldcc[3],
							Coder="cc")
				} else {
					scaffoldcc <- cbind(scaffoldcc,Coder=rep("cc",nrow(scaffoldcc)))
					colnames(scaffoldcc)<-c("FromTF", "Interaction", "ToTF",
						"Coder")
				}
				scaffoldtemp <- rbind(scaffoldfinal[,1:4], scaffoldcc)  
				scaffoldtemp <- scaffoldtemp[,c(4,1,2,3)]          
					dupRows <- dupsBetweenGroups(scaffoldtemp,"Coder")
				scaffoldtempDup<-cbind(scaffoldtemp, dup=dupRows)
				condcc<- subset(scaffoldtempDup, Coder=="final", select=-Coder)[,4]
			}
            scaffoldfinal<-cbind(scaffoldfinal,condcc)
            colnames(scaffoldfinal)[ncol(scaffoldfinal)]<-paste("CondT",
			timeSignals[cc],"ToT",timeSignals[(cc+1)],sep="")
		}
     	colnames(scaffoldfinal)[1:3]<-c("FromTF", "Interaction", "ToTF")
		rownames(scaffoldfinal)<-NULL
		removescaffold<-which(rownames(scaffoldfinal)=="")
		if(length(removescaffold)!=0){
			scaffoldfinal<-scaffoldfinal[-removescaffold,]
		}
		write.table(scaffoldfinal,file="Scaffold.complete.sif",append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
			col.names = TRUE, qmethod = c("escape", "double"))
		dataCNOOpt[["scaffoldfinal"]]<-scaffoldfinal
		write.table(finalBString,file="finalBString.txt",append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
			col.names = TRUE, qmethod = c("escape", "double"))
		if (any(CNOlist_all$namesCues=="dummyInh")){
			namesCues<-CNOlist_all$namesStimuli
			valueCues<-CNOlist_all$valueStimuli
		} else {
			namesCues<-CNOlist_all$namesCues
			valueCues<-CNOlist_all$valueCues
		}
		if (is.null(dim(valueCues))){
			valueCues<-matrix(valueCues,ncol=1)
		}
		UpdateCNO_all<-match(colnames(SimResultsdata_all[[1]]),CNOlist_all$namesSignals)
		indUpdateCNO_all<-UpdateCNO_all[!is.na(UpdateCNO_all)]
		CNOlist_all$namesSignals<-CNOlist_all$namesSignals[indUpdateCNO_all]
		for (dd in 1:length(timeSignals)){
			CNOlist_all$valueSignals[[dd]]<-CNOlist_all$valueSignals[[dd]][,indUpdateCNO_all]
		}
		if (clstr==0){
			plotOptimResults(
				simResults = SimResultsdata_all, 
				expResults = CNOlist_all$valueSignals, 
            	times = CNOlist_all$timeSignals, 
				namesCues = namesCues, 
				namesSignals = CNOlist_all$namesSignals, 
				valueCues =valueCues
			)
		}
		plotOptimResultsPDF(
			simResults = SimResultsdata_all, 
			expResults = CNOlist_all$valueSignals, 
			times = CNOlist_all$timeSignals, 
			filename = "SimulationResults_all.pdf", 
			namesCues = namesCues, 
			namesSignals = CNOlist_all$namesSignals, 
            valueCues = valueCues
		)
		print(randedgeper)
		#txtStop()
		setwd("~")
		#file.copy(paste(tempname,".txt",sep=""),paste(directory,"/",folder.output,sep=""))
		#file.remove(paste(tempname,".txt",sep=""))
		return(dataCNOOpt)
	}

#=====================================================================================	

#     cutNONC2
#=====================================================================================
#=====================================================================================

# The last if loop was not define appropiately in case that there was no reac2remove

cutNONC2<-function (Model, NONCindexes) 
{
    if (length(NONCindexes) == 0) {
        newModel = Model
    }
    else {
        newspecies <- Model$namesSpecies[-NONCindexes]
        newInterMat <- Model$interMat[-NONCindexes, ]
        newNotMat <- Model$notMat[-NONCindexes, ]
        EmptyInOut <- function(x) {
            input <- match(-1, x, nomatch = 0)
            output <- match(1, x, nomatch = 0)
            if ((input == 0) | (output == 0)) {
                return(TRUE)
            }
            else {
                return(FALSE)
            }
        }
        reac2remove <- apply(newInterMat, 2, EmptyInOut)
        if (any(reac2remove)) {
            reac2remove <- which(reac2remove)
            newInterMat <- newInterMat[, -reac2remove]
            newNotMat <- newNotMat[, -reac2remove]
            newreacID <- Model$reacID[-reac2remove]
		newModel <- list(reacID = newreacID, namesSpecies = newspecies, 
            	interMat = newInterMat, notMat = newNotMat)
        }
        else {
            newModel <- list(reacID = Model$reacID, namesSpecies = newspecies, 
                interMat = newInterMat, notMat = newNotMat)
        }
    }
    return(newModel)
}

#     expandGates3
#=====================================================================================
#=====================================================================================

# Change and expand the definition of the ignoreList to account to be a reactant not a
# product to remove

# Definitions
# 	ignoreList ID of the components whose reactions should not be expanded
#     avoid AND in between cues that are not co-present at least one condition

expandGates3<-function (Model, ignoreList = NA, maxInputsPerGate = 2,allGates=TRUE,List=CNOlist) 
{
    if (!is.list(Model)) 
        stop("This function expects as input a Model as output by readSif")
    if (length(Model) == 4) {
        if (all(names(Model) != c("reacID", "namesSpecies", "interMat", 
            "notMat"))) {
            stop("This function expects as input a Model as output by readSif")
        }
    }
    if (length(Model) == 5) {
        if (all(names(Model) != c("reacID", "namesSpecies", "interMat", 
            "notMat", "speciesCompressed"))) {
            stop("This function expects as input a Model as output by readSif")
        }
    }
    SplitANDs <- list(initialReac = c("split1", "split2"))
    splitR <- 1
    andToAdd = c()
    remove.and = c()
    reacs2Ignore = c()
    initialReacN <- length(Model$reacID)
    namesCues<-List$namesCues
    for (ii in 1:length(List$namesCues)){
	namesCues<-c(namesCues,paste("!",List$namesCues[ii],sep=""))
    }
    if (initialReacN == 1) {
        Model$interMat <- as.matrix(Model$interMat)
    }
    if (!is.na(ignoreList[1])) {
        for (s in 1:initialReacN) {
            if (any(Model$interMat[ignoreList, s] == 1)|| any(Model$interMat[ignoreList, s] == -1)) {
                reacs2Ignore = c(reacs2Ignore, s)
            }
        }
    }
    for (r in 1:initialReacN) {
        inNodes <- which(Model$interMat[, r] == -1)
        if (length(inNodes) > 1) {
            if (length(inNodes) > 3) {
                andToAdd = c(andToAdd, r)
            }
            remove.and = c(remove.and, r)
            if (!any(reacs2Ignore == r)) {
                outNode <- which(Model$interMat[, r] == 1)
                newReacs <- matrix(data = 0, nrow = dim(Model$interMat)[1], 
                  ncol = length(inNodes))
                newReacsNots <- matrix(data = 0, nrow = dim(Model$interMat)[1], 
                  ncol = length(inNodes))
                newReacs[outNode, ] <- 1
                newReacsIDs <- rep("a", length(inNodes))
                for (i in 1:length(inNodes)) {
                  newReacs[inNodes[i], i] <- -1
                  newReacsNots[inNodes[i], i] <- Model$notMat[inNodes[i], 
                    r]
                  newReacsIDs[i] <- paste(Model$namesSpecies[inNodes[i]], 
                    "=", Model$namesSpecies[outNode], sep = "")
                  if (Model$notMat[inNodes[i], r] == 1) 
                    newReacsIDs[i] <- paste("!", newReacsIDs[i], 
                      sep = "")
                }
                colnames(newReacs) <- newReacsIDs
                colnames(newReacsNots) <- newReacsIDs
                SplitANDs[[splitR]] <- newReacsIDs
                names(SplitANDs)[splitR] <- Model$reacID[r]
                splitR <- splitR + 1
                Model$notMat <- cbind(Model$notMat, newReacsNots)
                Model$interMat <- cbind(Model$interMat, newReacs)
                Model$reacID <- c(Model$reacID, newReacsIDs)
            } 
        }
    }
    if (length(andToAdd)) {
        toAdd = list()
        toAdd$notMat <- Model$notMat[, andToAdd]
        toAdd$interMat <- Model$interMat[, andToAdd]
        toAdd$reacID <- Model$reacID[andToAdd]
    }
    else {
        toAdd <- NA
    }
    if (length(remove.and)) {
        Model$notMat <- Model$notMat[, -remove.and]
        Model$interMat <- Model$interMat[, -remove.and]
        Model$reacID <- Model$reacID[-remove.and]
    }
    newANDs <- list(finalReac = c("or1", "or2"))
    ANDsadded <- 1
    total.list = 1:length(Model$namesSpecies)
    getlhs <- function(x) {
        spec1 = strsplit(x, "=")[[1]][1]
    }
    getrhs <- function(x) {
        spec2 = strsplit(x, "=")[[1]][2]
    }
    for (sp in total.list) {
        inReacsIndex <- which(Model$interMat[sp, ] == 1)
	  if (length(reacs2Ignore)>0){
        	for (ir2I in 1:length(reacs2Ignore)){
	  		if (reacs2Ignore[ir2I]%in%inReacsIndex){
	  			inReacsIndex<-inReacsIndex[-which(inReacsIndex==reacs2Ignore[ir2I])]
	  		}
	  	}
	  }
        if (length(inReacsIndex) > 1) {
            inReacs <- Model$interMat[, inReacsIndex]
            findInput <- function(x) {
                inNode <- which(x == -1)
            }
            inSp <- apply(inReacs, 2, findInput)
            inSpecies = apply(as.matrix(colnames(inReacs)), 1, 
                getlhs)
            outname = Model$namesSpecies[sp]
            outnames = apply(as.matrix(colnames(inReacs)), 1, 
                getrhs)
            if (length(unique(outnames)) != 1 | outname != outnames[1]) {
                stop("error in expandGates. should not happen here.please report")
            }
            myrownames = rownames(Model$interMat)
            combinations = combn(seq(1, length(inSpecies)), 2)
            for (this in seq(1, dim(combinations)[2])) {
                i = combinations[1, this]
                j = combinations[2, this]
                realname1 = ifelse(substr(inSpecies[i], 1, 1) == 
                  "!", substr(inSpecies[i], 2, 10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1, 1) == 
                  "!", substr(inSpecies[j], 2, 10000), inSpecies[j])
                realnames = c(realname1, realname2)
                if (any(combn(realnames, 2)[1, ] == combn(realnames, 
                  2)[2, ])) {
                  next()
                }
                newcolname = paste(paste(inSpecies[i], inSpecies[j], 
                  sep = "+"), outname, sep = "=")
                if (newcolname %in% colnames(Model$interMat)) {
                  next()
                }
                if (!allGates && all(realnames%in%namesCues)){
			next()
		    }
                Model$reacID <- c(Model$reacID, newcolname)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                values[which(myrownames == realname1)] <- -1
                values[which(myrownames == realname2)] <- -1
                values[which(myrownames == outname)] <- 1
                Model$interMat = cbind(Model$interMat, values)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                if (substr(inSpecies[i], 1, 1) == "!") {
                  values[which(myrownames == realname1)] <- 1
                }
                if (substr(inSpecies[j], 1, 1) == "!") {
                  values[which(myrownames == realname2)] <- 1
                }
                Model$notMat = cbind(Model$notMat, values)
                newreac1 = paste(inSpecies[i], outname, sep = "=")
                newreac2 = paste(inSpecies[j], outname, sep = "=")
                newANDs[[length(newANDs) + 1]] <- c(newreac1, 
                  newreac2)
                names(newANDs)[[length(newANDs)]] <- newcolname
            }
            if (length(inSpecies) >= 3 & maxInputsPerGate >= 
                3) {
                combinations = combn(seq(1, length(inSpecies)), 
                  3)
                indices = seq(1, dim(combinations)[2])
            }
            else {
                indices = seq(length = 0)
            }
            for (this in indices) {
                i = combinations[1, this]
                j = combinations[2, this]
                k = combinations[3, this]
                realname1 = ifelse(substr(inSpecies[i], 1, 1) == 
                  "!", substr(inSpecies[i], 2, 10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1, 1) == 
                  "!", substr(inSpecies[j], 2, 10000), inSpecies[j])
                realname3 = ifelse(substr(inSpecies[k], 1, 1) == 
                  "!", substr(inSpecies[k], 2, 10000), inSpecies[k])
                realnames = c(realname1, realname2, realname3)
                if (any(combn(realnames, 2)[1, ] == combn(realnames, 
                  2)[2, ])) {
                  next()
                }
                newcolname <- paste(paste(inSpecies[i], inSpecies[j], 
                  inSpecies[k], sep = "+"), outname, sep = "=")
                if (newcolname %in% colnames(Model$interMat)) {
                  next()
                }
                if (!allGates && all(realnames%in%namesCues)){
			next()
		    }
                Model$reacID <- c(Model$reacID, newcolname)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                for (name in inSpecies) {
                  realname = ifelse(substr(name, 1, 1) == "!", 
                    substr(name, 2, 10000), name)
                  values[which(myrownames == realname)] <- -1
                }
                values[which(myrownames == outname)] <- 1
                Model$interMat = cbind(Model$interMat, values)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                for (name in inSpecies) {
                  if (substr(name, 1, 1) == "!") {
                    realname = ifelse(substr(name, 1, 1) == "!", 
                      substr(name, 2, 10000), name)
                    values[which(myrownames == realname)] <- 1
                  }
                }
                Model$notMat = cbind(Model$notMat, values)
                newreac1 = paste(inSpecies[i], outname, sep = "=")
                newreac2 = paste(inSpecies[j], outname, sep = "=")
                newreac3 = paste(inSpecies[k], outname, sep = "=")
                newANDs[[length(newANDs) + 1]] <- c(newreac1, 
                  newreac2, newreac3)
                names(newANDs)[[length(newANDs)]] <- newcolname
            }
            if (length(inSpecies) >= 4 & maxInputsPerGate >= 
                4) {
                combinations = combn(seq(1, length(inSpecies)), 
                  4)
                indices = seq(1, dim(combinations)[2])
            }
            else {
                indices = seq(length = 0)
            }
            for (this in indices) {
                i = combinations[1, this]
                j = combinations[2, this]
                k = combinations[3, this]
                l = combinations[4, this]
                realname1 = ifelse(substr(inSpecies[i], 1, 1) == 
                  "!", substr(inSpecies[i], 2, 10000), inSpecies[i])
                realname2 = ifelse(substr(inSpecies[j], 1, 1) == 
                  "!", substr(inSpecies[j], 2, 10000), inSpecies[j])
                realname3 = ifelse(substr(inSpecies[k], 1, 1) == 
                  "!", substr(inSpecies[k], 2, 10000), inSpecies[k])
                realname4 = ifelse(substr(inSpecies[l], 1, 1) == 
                  "!", substr(inSpecies[l], 2, 10000), inSpecies[l])
                realnames = c(realname1, realname2, realname3, 
                  realname4)
                if (any(combn(realnames, 2)[1, ] == combn(realnames, 
                  2)[2, ])) {
                  next()
                }
                newcolname <- paste(paste(inSpecies[i], inSpecies[j], 
                  inSpecies[k], inSpecies[l], sep = "+"), outname, 
                  sep = "=")
                if (newcolname %in% colnames(Model$interMat)) {
                  next()
                }
                if (!allGates && all(realnames%in%namesCues)){
			next()
		    }
                Model$reacID <- c(Model$reacID, newcolname)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                for (name in inSpecies) {
                  realname = ifelse(substr(name, 1, 1) == "!", 
                    substr(name, 2, 10000), name)
                  values[which(myrownames == realname)] <- -1
                }
                values[which(myrownames == outname)] <- 1
                Model$interMat = cbind(Model$interMat, values)
                values = as.matrix(rep(0, length(Model$namesSpecies)))
                colnames(values) <- newcolname
                for (name in inSpecies) {
                  if (substr(name, 1, 1) == "!") {
                    realname = ifelse(substr(name, 1, 1) == "!", 
                      substr(name, 2, 10000), name)
                    values[which(myrownames == realname)] <- 1
                  }
                }
                Model$notMat = cbind(Model$notMat, values)
                newreac1 = paste(inSpecies[i], outname, sep = "=")
                newreac2 = paste(inSpecies[j], outname, sep = "=")
                newreac3 = paste(inSpecies[k], outname, sep = "=")
                newreac4 = paste(inSpecies[l], outname, sep = "=")
                newANDs[[length(newANDs) + 1]] <- c(newreac1, 
                  newreac2, newreac3, newreac4)
                names(newANDs)[[length(newANDs)]] <- newcolname
            }
        }
    }
    if (!is.na(toAdd)) {
        Model$notMat = cbind(Model$notMat, toAdd$notMat)
        Model$interMat = cbind(Model$interMat, toAdd$interMat)
        Model$reacID = c(Model$reacID, toAdd$reacID)
    }
    ModelExp <- Model
    ModelExp$SplitANDs <- SplitANDs
    ModelExp$newANDs <- newANDs
    return(ModelExp)
}


#     expandDAGates
#=====================================================================================
#=====================================================================================

# Create DA and Stimuli specific AND gates

expandDAGates<-function (Model,List=CNOlist) 
{
    if (!is.list(Model)) 
        stop("This function expects as input a Model as output by readSif")
    if (length(Model) == 4) {
        if (all(names(Model) != c("reacID", "namesSpecies", "interMat", 
            "notMat"))) {
            stop("This function expects as input a Model as output by readSif")
        }
    }
    if (length(Model) == 5) {
        if (all(names(Model) != c("reacID", "namesSpecies", "interMat", 
            "notMat", "speciesCompressed"))) {
            stop("This function expects as input a Model as output by readSif")
        }
    }
    SplitANDs <- list(initialReac = c("split1", "split2"))
    newANDs <- list(finalReac = c("or1", "or2"))
    namesStimuli<-List$namesStimuli
    namesStimuli<-namesStimuli[which(namesStimuli!="DA")]
    namesSignals<-List$namesSignals
    for (i in 1:length(namesStimuli)) {
    	for (j in 1:length(namesSignals)){
    		newcolname = paste("!",namesStimuli[i], "+!DA=",
			namesSignals[j], sep = "")
            Model$reacID <- c(Model$reacID, newcolname)
            values = as.matrix(rep(0, length(Model$namesSpecies)))
            colnames(values) <- newcolname
            values[which(Model$namesSpecies == namesStimuli[i])] <- -1
            values[which(Model$namesSpecies == "DA")] <- -1
            values[which(Model$namesSpecies == namesSignals[j])] <- 1
            Model$interMat = cbind(Model$interMat, values)
            values = as.matrix(rep(0, length(Model$namesSpecies)))
            colnames(values) <- newcolname
            values[which(Model$namesSpecies == "DA")] <- 1
            values[which(Model$namesSpecies == namesStimuli[i])] <- 1
            Model$notMat = cbind(Model$notMat, values)
            newreac1 = paste("!",namesStimuli[i], namesSignals[j], sep = "=")
            newreac2 = paste("!DA", namesSignals[j], sep = "=")
            newANDs[[length(newANDs) + 1]] <- c(newreac1, 
            	newreac2)
            names(newANDs)[[length(newANDs)]] <- newcolname
        }
    }
    ModelExp <- Model
    ModelExp$reacID<-Model$reacID[-which(colSums(Model$interMat)==0 & Model$interMat["DA",]==-1)]
    ModelExp$interMat<-Model$interMat[,-which(colSums(Model$interMat)==0 & Model$interMat["DA",]==-1)]
    ModelExp$notMat<-Model$notMat[,-which(colSums(Model$interMat)==0 & Model$interMat["DA",]==-1)]
    ModelExp$SplitANDs <- SplitANDs
    ModelExp$newANDs <- newANDs
    return(ModelExp)
}


#	simulatorT3
#=====================================================================================
#=====================================================================================

# This function is a modification of simulatorT1 from CellNOptR to account for stimuli
# that can be initially activated

# Created by Beatriz Penalver 08/08/2012

simulatorT3<-function (CNOlist, Model, SimList, indexList,ResPrior=ResPrior,time0=time0,
	outputPrevChangeSp=outputChangeSp,indCtrlSig=indCtrlSig) 
{
    if(!is.matrix(SimList$finalCube) || ncol(SimList$finalCube)==1){
    	SimList$finalCube<-matrix(cbind(SimList$finalCube,rep(1,length(SimList$finalCube))),
	   nrow=length(SimList$finalCube), ncol=2)
    	SimList$ixNeg<-matrix(cbind(SimList$ixNeg,rep(FALSE,length(SimList$ixNeg))),
	   nrow=length(SimList$ixNeg), ncol=2)
    	SimList$ignoreCube<-matrix(cbind(SimList$ignoreCube,rep(TRUE,length(SimList$ignoreCube))),
	   nrow=length(SimList$ignoreCube), ncol=2)
    }
    testVal <- 0.001
    initValues <- ResPrior
    valueInhibitors <- 1 - CNOlist$valueInhibitors
    valueInhibitors[which(valueInhibitors == 1)] <- NA
    valueStimulus<-CNOlist$valueStimuli
    valueStimulus[which(valueStimulus == 0)] <- NA
    initValues[, indexList$inhibited] <- valueInhibitors
    initValues[, indexList$stimulated]<-valueStimulus	
    newInput <- initValues
    termCheck1 <- TRUE
    termCheck2 <- TRUE
    count <- 1
    nSp <-length(Model$namesSpecies) 
    if (is.null(dim(CNOlist$valueStimuli))){
	nCond<-length(CNOlist$valueStimuli)
    } else {
    	nCond <- dim(CNOlist$valueStimuli)[1]
    }
    valueSignalsEnd<-as.matrix(CNOlist$valueSignals[[2]])
    if (!is.null(dim(CNOlist$valueCues)) && nrow(CNOlist$valueCues)>1){
    	uniquecond<-unique(apply(CNOlist$valueCues,1,function(x){paste(x,collapse="")}))
    	for (ii in 1:length(uniquecond)){
		tmppos<-which(uniquecond[ii]==apply(CNOlist$valueCues,1,function(x){paste(x,collapse="")}))
		tmpmean<-matrix(apply(valueSignalsEnd[tmppos,],2,function(x){round(median(x,na.rm=TRUE))}),nrow=1)
		valueSignalsEnd[tmppos,]<- tmpmean[rep(nrow(tmpmean),length(tmppos)),]
    	}
    } else {
    	uniquecond<-unique(CNOlist$valueCues)
    	for (ii in 1:length(uniquecond)){
		tmppos<-which(uniquecond[ii]==CNOlist$valueCues)
		tmpmean<-matrix(apply(valueSignalsEnd[tmppos,],2,function(x){round(median(x,na.rm=TRUE))}),nrow=1)
		valueSignalsEnd[tmppos,]<- tmpmean[rep(nrow(tmpmean),length(tmppos)),]
    	}
    }
    predictedEnd<-ResPrior
    predictedEnd[, indexList$stimulated] <- CNOlist$valueStimuli
    predictedEnd[, indexList$inhibited] <- CNOlist$valueInhibitors
    predictedEnd[,indexList$signals]<-valueSignalsEnd
    outputPrevRateChangeSp<-predictedEnd-ResPrior
    outputPrevChangeSp0<-outputPrevChangeSp
    while (termCheck1 && termCheck2) {
 	indicesChangedSp<- which(colSums(outputPrevChangeSp)>0)
	SimListfinaltmp<-SimList$finalCube
	SimListfinaltmp[SimList$ignoreCube]<-NA
    	rxInpChanged<-apply(SimListfinaltmp,1,function(x){
		r<-ifelse(all(x[!is.na(x)] %in% indicesChangedSp),TRUE,FALSE)
	})
      SimListRx<-SimList
	SimListRx$finalCube<-SimList$finalCube[rxInpChanged,]
	SimListRx$ixNeg<-SimList$ixNeg[rxInpChanged,]
 	SimListRx$ignoreCube<-SimList$ignoreCube[rxInpChanged,]
	SimListRx$maxIx<-SimList$maxIx[rxInpChanged]     
	nReacs <- sum(rxInpChanged)
      if (nReacs==1){
		maxIpg<-length(SimListRx$finalCube)
	} 
	if (nReacs>1){
      	maxIpg <- dim(SimListRx$finalCube)[2]  
	}
      if (nReacs==0){
		maxIpg=2
	}
      endIx <- rep(NA, nSp)
      for (i in 1:nSp) {
		endIx[i] <- length(which(SimListRx$maxIx == i))
    	}
   	maxgpo <- max(endIx)
     	if (nReacs>0){
      	outputPrev <- newInput
		outputPrevCor<-newInput
            if (length(indCtrlSig)>1){
    			outputPrevCor[,indCtrlSig]<-apply(outputPrevCor[,indCtrlSig],c(1,2),function(x){
				if(x==0){
					r<-1
				} else {
					r<-x
				}
				return(r)
    			})
		} 
		if (length(indCtrlSig)==1){
			outputPrevCor[which(outputPrevCor[,indCtrlSig]==0),indCtrlSig]<-1
		}
		outputPrevCor<-outputPrevCor*outputPrevChangeSp 
            outputPrevCor[, indexList$stimulated]<-valueStimulus
 		outputPrevCor[, indexList$inhibited] <- valueInhibitors
        	filltempCube <- function(x) {
            	cMatrix <- matrix(data = x, nrow = nReacs, ncol = nCond)
            	cVector <- apply(cMatrix, 1, function(x) {
            	    return(x)
            	})
            	return(cVector)
        	}
        	if (nReacs > 1) {
            	tempStore <- apply(SimListRx$finalCube, 2, function(x) {
                		return(outputPrevCor[, x])
            	})
            	tempStoreChangeSp <- apply(SimListRx$finalCube, 2, function(x) {
                		return(outputPrevChangeSp[, x])
            	})
            	tempIxNeg <- apply(SimListRx$ixNeg, 2, filltempCube)
            	tempIgnore <- apply(SimListRx$ignoreCube, 2, filltempCube)
			tempRate<-outputPrevRateChangeSp[,SimListRx$maxIx]
			tempValue<-ResPrior[,SimListRx$maxIx]
			tempRate<-matrix(tempRate,ncol=1)
			tempValue<-matrix(tempValue,ncol=1)
        	}else {
            	tempStore <- outputPrevCor[, SimListRx$finalCube]
			tempStoreChangeSp <- outputPrevChangeSp[, SimListRx$finalCube]
            	tempIxNeg <- matrix(SimListRx$ixNeg, nrow = nCond, 
                		ncol = length(SimListRx$ixNeg), byrow = TRUE)
            	tempIgnore <- matrix(SimListRx$ignoreCube, nrow = nCond, 
                		ncol = length(SimListRx$ignoreCube), byrow = TRUE)
			tempRate<-outputPrevRateChangeSp[,SimListRx$maxIx]
			tempValue<-ResPrior[,SimListRx$maxIx]
        	}
		tempStore[which(tempStore==(-1))]<-0
        	tempStore[tempIgnore] <- 100
    	  	tempStore[tempIxNeg] <- 1 - tempStore[tempIxNeg]
		tempStore[which(tempStore==0 & !tempIxNeg)]<-NA
            tempStoreChangeSp[tempIgnore]<-NA
            gateType<-function(x){
			r=sum(!is.na(x),na.rm=TRUE)
           		t=sum(any(x,na.rm=TRUE)) 
                  v=sum(any(!x,na.rm=TRUE))
                  gate=4
                  if (r==1 && v==1){
				gate=0
			}
                  if (r==1 && t==1){
				gate=1
			}
			if (r>1 && t==0) {
				gate=2
			}
			if (r>1 && v==0){
				gate=3
			}
                  return(gate)
		}
		outputFinal<-function(x){
			r<-ifelse(!is.na(x[1]) && x[1]==0,-1,x[1])
			u<-r
			if (!is.na(r)){
				Cond2<- x[3]==1 && x[1]==1
				Cond4<- x[3]==3 && !is.na(x[2]) && (x[1]==1 || (x[1]==0 && x[2]>0))
				Cond5<- x[3]==4 && (x[1]==0)
				if (Cond2||Cond4|| Cond5){
					u<-NA
				}
				if (!is.na(u) && !is.na(x[5]) && ((r==(-1) && x[5]==(-1) && x[6]==1)||
			        (r==1 && x[5]==1 && x[6]==(-1)))){
					u=0
				}
			}
			return(u)
		}    
		outputCubeMin <- apply(tempStore, 1,function(x){ifelse(all(is.na(x)),0,min(x))})
            outputCubeSum<-apply(tempStore,1,function(x){r = sum(x)-(100*sum(x==100)) ; return (r)})
            tempCubeGate<-tempIxNeg
            tempCubeGate[tempIgnore]<-NA
            outputCubeGate<-apply(tempCubeGate,1,gateType)
            outputCubeChangeSp<-apply(tempStoreChangeSp,1,function(x){max(x,na.rm=TRUE)})
            outputCubeRate<-tempRate
		outputCubeValue<-tempValue
		tempCubeAll<-cbind(outputCubeMin,outputCubeSum,outputCubeGate,outputCubeChangeSp,outputCubeRate,
			tempValue)
		outputCubeList<-apply(tempCubeAll,1,outputFinal)
		outputCube <- matrix(outputCubeList, nrow = nCond, ncol = nReacs)
            for (s in 1:nSp) {
			if (endIx[s]!=0){
                  	compOR <- function(x) {
					if (all(is.na(x[which(SimListRx$maxIx == s)]))) {
                      			resCube <- NA
                    		} else {
                                    tempres<-x[which(SimListRx$maxIx == s)]
                    			tempsum<- sum(tempres,na.rm = TRUE)
						if (all(tempres[!is.na(tempres)]==0)){
							resCube<-0
						}
						if (all(tempres[!is.na(tempres)]==1)){
							resCube<-1
						}
						if (all(tempres[!is.na(tempres)]==(-1))){
							resCube<-(-1)
						}
						if (!all(tempres[!is.na(tempres)]==0)&& tempsum>0){
							resCube<-1
						}
						if (!all(tempres[!is.na(tempres)]==0) && tempsum<0){
							resCube<-(-1)
						} 
						if (!all(tempres[!is.na(tempres)]==0) && tempsum==0){
							resCube<-NA
						} 
					}
                  		return(resCube)
				} 
		 		newInput[, s] <-apply(outputCube, 1, compOR)
			} else {
				newInput[, s] <- ResPrior[,s]
			}
		}
		for (stim in 1:length(indexList$stimulated)) {
                  dimCNOStim<-is.null(dim(CNOlist$valueStimuli))
			if(dimCNOStim){
				CNOlist$valueStimuli<-matrix(CNOlist$valueStimuli,ncol=1)
			} 
      		stimM <- cbind(CNOlist$valueStimuli[, stim], newInput[, 
                		indexList$stimulated[stim]])
            	maxNA <- function(x) {
                		return(max(x, na.rm = TRUE))
            	}
            	stimV <- apply(stimM, 1, maxNA)
            	newInput[, indexList$stimulated[stim]] <- stimV
      	}
      	newInput[, indexList$inhibited] <- (1 - CNOlist$valueInhibitors)* 
           		newInput[, indexList$inhibited]
      	newInput[is.na(newInput)] <-ResPrior[is.na(newInput)] 
      	outputPrev[is.na(outputPrev)] <- ResPrior[is.na(outputPrev)]
            outputPrevChangeSp<-newInput-ResPrior
		outputPrevChangeSp<-apply(outputPrevChangeSp,c(1,2),function(x){ifelse(x!=0,1,0)})
		outputPrevChangeSp[,indexList$stimulated]<-CNOlist$valueStimuli
		outputPrevChangeSp[,indexList$inhibited]<-CNOlist$valueInhibitors 
		outputPrevChangeSp[which(outputPrevChangeSp0==1)]<-1
     		termCheck1 <- !all(abs(outputPrev - newInput) < testVal)
      	termCheck2 <- (count < (nSp*0.4))
     		count <- count + 1
    	} else {
		newInput<-ResPrior
		outputPrevChangeSp<-newInput-ResPrior
		outputPrevChangeSp<-apply(outputPrevChangeSp,c(1,2),function(x){ifelse(x!=0,1,0)})
		outputPrevChangeSp[,indexList$stimulated]<-CNOlist$valueStimuli
		outputPrevChangeSp[,indexList$inhibited]<-CNOlist$valueInhibitors  
		outputPrevChangeSp[which(outputPrevChangeSp0==1)]<-1
      	termCheck1<-FALSE
            termCheck2<-FALSE
	}
 	indicesChangedSp<- which(colSums(outputPrevChangeSp)>0)
	tempSimList<-SimList$finalCube
	tempSimList[SimList$ignoreCube]<-NA
    	rxInpChanged<-apply(tempSimList,1,function(x){r<-ifelse(all(x[!is.na(x)] %in% indicesChangedSp),TRUE,FALSE)})
   }
   if (nReacs>0){	
	newInput[which(abs(outputPrev - newInput) > testVal)] <- NA
    	#newInput[is.na(newInput)] <-ResPrior[is.na(newInput)]
   }
   simres3<-list(newInput=newInput,rxInpChanged=rxInpChanged)
   return(simres3)
}
#=====================================================================================



#	simulateT3
#=====================================================================================
#=====================================================================================

# This function is a modification of simulatorT1 from CellNOptR to account for stimuli
# that can be initially activated

# Created by Beatriz Penalver 08/08/2012

 simulateT3<-function (CNOlist, Model, bStringT1, SimList, indexList,ResPrior=SimResPrior,time0=time0,
	outputChangeSp=ChangeSp,indCtrlSig=indicesCtrlSignals) 
{
    Modelcut <- Model
    Modelcut$interMat <- Modelcut$interMat[, as.logical(bStringT1)]
    Modelcut$notMat <- Modelcut$notMat[, as.logical(bStringT1)]
    Modelcut$reacID <- Modelcut$reacID[as.logical(bStringT1)]
    SimListcut <- SimList
    SimListcut$finalCube <- SimListcut$finalCube[as.logical(bStringT1), 
        ]
    SimListcut$ixNeg <- SimListcut$ixNeg[as.logical(bStringT1), 
        ]
    SimListcut$ignoreCube <- SimListcut$ignoreCube[as.logical(bStringT1), 
        ]
    SimListcut$maxIx <- SimListcut$maxIx[as.logical(bStringT1)]
    if (is.null(dim(SimListcut$finalCube))) {
        SimListcut$finalCube <- matrix(SimListcut$finalCube, 
            ncol = 1)
        SimListcut$ixNeg <- matrix(SimListcut$ixNeg, ncol = 1)
        SimListcut$ignoreCube <- matrix(SimListcut$ignoreCube, 
            ncol = 1)
    }
    SimRes <- simulatorT3(CNOlist = CNOlist, Model = Modelcut, 
        SimList = SimListcut, indexList = indexList,ResPrior=ResPrior,time0=time0,
	  outputPrevChangeSp=outputChangeSp,indCtrlSig=indCtrlSig)
    return(SimRes)
}
#=====================================================================================



#	cutAndPlotResultsT3
#=====================================================================================
#=====================================================================================

# This function is a modification of simulatorT1 from CellNOptR to account for stimuli
# that can be initially activated

# Created by Beatriz Penalver 08/08/2012

cutAndPlotResultsT3<-function (Model, bString, SimList, CNOlist, indexList, plotPDF = FALSE, 
    tag = NULL, show = TRUE,ResPrior=SimResPrior,time0=timeinit,outputChangeSp=ChangeSp,
    indCtrlSig=indicesCtrlSignals) 
{
    Modelcut <- Model
    Modelcut$interMat <- Modelcut$interMat[, as.logical(bString)]
    Modelcut$notMat <- Modelcut$notMat[, as.logical(bString)]
    Modelcut$reacID <- Modelcut$reacID[as.logical(bString)]
    SimListCut <- cutSimList3(SimList, bString)
    Sim <- simulatorT3(CNOlist = CNOlist, Model = Modelcut, SimList = SimListCut, 
        indexList = indexList,ResPrior=ResPrior,time0=time0,outputPrevChangeSp=outputChangeSp,
	  indCtrlSig=indCtrlSig)
    SimRes <- as.matrix(Sim$newInput[, indexList$signals])
    SimRest0<-as.matrix(ResPrior[, indexList$signals])
    SimResults <- list(t0 = SimRest0, t1 = SimRes)
    expResults <- list(t0 = CNOlist$valueSignals[[1]], t1 = CNOlist$valueSignals[[2]])
    if (any(CNOlist$namesInhibitors=="dummyInh")){
	namesCues<-CNOlist$namesStimuli
	valueCues<-CNOlist$valueStimuli
    } else {
	namesCues<-CNOlist$namesCues
	valueCues<-CNOlist$valueCues
    }
    if (is.null(dim(valueCues))){
	valueCues<-matrix(valueCues,ncol=1)
    }
    if (show == TRUE) {
        plotOptimResults(simResults = SimResults, expResults = expResults, 
            times = CNOlist$timeSignals[1:2], namesCues = namesCues, 
            namesSignals = CNOlist$namesSignals, valueCues = valueCues)
     }
    if (plotPDF == TRUE) {
        if (is.null(tag)) {
            filename <- paste(deparse(substitute(Model)), "SimResultsT1.pdf", 
                sep = "")
        }
        else {
            filename <- paste(tag, "SimResultsT1.pdf", sep = "_")
        }
        plotOptimResultsPDF(simResults = SimResults, expResults = expResults, 
            times = CNOlist$timeSignals[1:2], filename = filename, 
            namesCues = namesCues, namesSignals = CNOlist$namesSignals, 
            valueCues = valueCues)
    }
}
#=====================================================================================



#	gaBinaryT3
#=====================================================================================
#=====================================================================================

# This function is a modification of simulatorT1 from CellNOptR to account for stimuli
# that can be initially activated

# Created by Beatriz Penalver 08/08/2012

gaBinaryT3<-function (CNOlist, Model, SimList, indexList, sizeFac = 1e-04, 
    NAFac = 1, initBstring, PopSize = 50, Pmutation = 0.5, MaxTime = 60, 
    maxGens = 500, StallGenMax = 100, SelPress = 1.2, elitism = 5, 
    RelTol = 0.1, verbose = TRUE,ResPrior=SimResPrior,time0=timeinit,
    outputChangeSp=outputChangeSp,StimuliFac=StimuliFac,indCtrlSig=indicesCtrlSignals,
    DAFac=DAFac,OppEdgesList=OppEdgesList,DAEdgesList=DAEdgesList) 
{
    bLength <- length(initBstring)
    Pop<-matrix(rep(initBstring,PopSize),nrow =PopSize,byrow=1)
    PoppossibleRx <- round(matrix(runif(bLength * (PopSize - 
        1)), nrow = (PopSize - 1), ncol = bLength))
    Pop[-1,]<-PoppossibleRx
    bestbit <- Pop[1, ]
    bestobj <- Inf
    stop <- FALSE
    obj <- rep(0, PopSize)
    g <- 0
    stallGen <- 0
    res <- rbind(c(g, bestobj, toString(bestbit), stallGen, Inf, 
        Inf, toString(bestbit), 0), c(g, bestobj, toString(bestbit), 
        stallGen, Inf, Inf, toString(bestbit), 0))
    colnames(res) <- c("Generation", "Best_score", "Best_bitString", 
        "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
        "Best_bit_Gen", "Iter_time")
    PopTol <- rep(NA, bLength)
    PopTolScores <- NA
    getObj <- function(x) {
        bitString <- x
        ModelCut <- Model
        ModelCut$interMat <- ModelCut$interMat[, as.logical(bitString)]
        ModelCut$notMat <- ModelCut$notMat[, as.logical(bitString)]
        ModelCut$reacID <- ModelCut$reacID[as.logical(bitString)]
        SimListCut <- cutSimList3(SimList, bitString)
        SimResultsAll <- simulatorT3(CNOlist = CNOlist, Model = ModelCut, 
            SimList = SimListCut, indexList = indexList, ResPrior=ResPrior,time0=time0,
		outputPrevChangeSp=outputChangeSp,indCtrlSig=indCtrlSig)
        SimResults<-SimResultsAll$newInput
        bitStringUpdated<-bitString
        bitStringUpdated[as.logical(bitStringUpdated)][!SimResultsAll$rxInpChanged]<-0
        ModelCutUpdated <- Model
        ModelCutUpdated$interMat <- Model$interMat[, as.logical(bitStringUpdated)]
        ModelCutUpdated$notMat <- Model$notMat[, as.logical(bitStringUpdated)]
        ModelCutUpdated$reacID <- Model$reacID[as.logical(bitStringUpdated)]
        SimResultsT0 <- ResPrior
	  nInDATot<-ifelse(DAFac==0,0,length(which(Model$interMat["DA",] == -1)))
	  nInStimTot<-length(which(Model$interMat[indexList$stimulated,] == -1))-nInDATot
	  nInTot<-length(which(Model$interMat == -1))
        Score <- getFit3(SimResults = SimResults, SimResultsT0 = SimResultsT0, 
            CNOlist = CNOlist, Model = ModelCutUpdated, indexList = indexList, 
            timePoint = "t1", sizeFac = sizeFac, NAFac = NAFac, 
            nInTot = nInTot,nInStimTot = nInStimTot,nInDATot = nInDATot,
		StimuliFac=StimuliFac,DAFac=DAFac)
        ScoreAll<-c(Score,bitStringUpdated)
        return(ScoreAll)
    }
    t0 <- Sys.time()
    t <- t0
    while (!stop) {
 	if (!is.null(DAEdgesList) && !is.null(OppEdgesList)){
		if(StimuliFac <= DAFac*sizeFac){
			PopUpdated<-apply(Pop,1,function(x){
				x[DAEdgesList[,2]]<-0
				return(x)
			})
		} else {
			PopUpdated<-apply(Pop,1,function(x){
				x[DAEdgesList[,1]]<-0
				return(x)
			})
		}
    	  }
	  Pop<-t(PopUpdated)
    	  if (!is.null(OppEdgesList)){
		OppEdgesList2<-rbind(OppEdgesList,
			cbind(OppEdgesList[,2],OppEdgesList[,1]))
		PopPrev<-Pop
		PopUpdated<-Pop
            conv<-FALSE
		while(!conv){
			tempOpp<-apply(OppEdgesList2,2,function(x){
				return(PopPrev[,x])
			})
			tempOpp2<-tempOpp[,1]*tempOpp[,2] 
			tempOpp2[which(tempOpp2==1)]<-(-1)
			Poptemp<-matrix(tempOpp2,nrow=PopSize)
			PopUpdated[,OppEdgesList2[,1]]<-PopUpdated[,OppEdgesList2[,1]]+Poptemp
			diffPop<-sum(PopUpdated-PopPrev)
			conv<-ifelse(diffPop==0,TRUE, FALSE)
			PopPrev<-PopUpdated
		}
    	  }
    	  Pop<-PopUpdated
        scoresAll <- t(apply(Pop, 1, getObj))
        scores<-as.vector(scoresAll[,1])
        Pop<-as.matrix(scoresAll[,-1])
        rankP <- order(scores, decreasing = TRUE)
        Pop <- Pop[rankP, ]   
        scores <- scores[rankP]
        fitness <- 2 - SelPress + (2 * (SelPress - 1) * (c(1:PopSize) - 
            1)/(PopSize - 1))
        wheel1 <- cumsum(fitness/sum(fitness))
        breaks <- runif(1) * 1/PopSize
        breaks <- c(breaks, breaks + ((1:(PopSize - 1)))/PopSize)
        sel <- rep(1, PopSize)
        for (i in 1:length(breaks)) {
            sel[i] <- which(wheel1 > breaks[i])[1]
        }
        Pop2 <- Pop[sel, ]
        PSize2 <- dim(Pop2)[1]
        PSize3 <- PopSize - elitism
        mates <- cbind(ceiling(runif(PSize3) * PSize2), ceiling(runif(PSize3) * 
            PSize2))
        InhBit <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
            ncol = bLength)
        InhBit <- InhBit < 0.5
        Pop3par1 <- Pop2[mates[, 1], ]
        Pop3par2 <- Pop2[mates[, 2], ]
        Pop3 <- Pop3par2
        Pop3[InhBit] <- Pop3par1[InhBit]
        MutProba <- matrix(runif((PSize3 * bLength)), nrow = PSize3, 
            ncol = bLength)
        MutProba <- MutProba < Pmutation
        Pop3[MutProba] <- 1 - Pop3[MutProba]
        t <- c(t, Sys.time())
        g <- g + 1
        thisGenBest <- scores[length(scores)]
        thisGenBestBit <- Pop[length(scores), ]
        if (is.na(thisGenBest)) {
            thisGenBest <- min(scores, na.rm = TRUE)
            thisGenBestBit <- Pop[which(scores == thisGenBest)[1], 
                ]
        }
        if (thisGenBest < bestobj) {
            bestobj <- thisGenBest
            bestbit <- thisGenBestBit
            stallGen <- 0
        }
        else {
            stallGen <- stallGen + 1
        }
        resThisGen <- c(g, bestobj, toString(bestbit), stallGen, 
            (mean(scores, na.rm = TRUE)), thisGenBest, toString(thisGenBestBit), 
            as.numeric((t[length(t)] - t[length(t) - 1]), units = "secs"))
        names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", 
            "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", 
            "Best_bit_Gen", "Iter_time")
        if (verbose) 
            print(resThisGen)
        res <- rbind(res, resThisGen)
        Criteria <- c((stallGen > StallGenMax), (as.numeric((t[length(t)] - 
            t[1]), units = "secs") > MaxTime), (g > maxGens))
        if (any(Criteria)) {
            stop <- TRUE
        }
        tolScore <- scores[length(scores)] * RelTol
        TolBs <- which(scores < scores[length(scores)] + tolScore)
        if (length(TolBs) > 0) {
            PopTol <- rbind(PopTol, Pop[TolBs, ])
            PopTolScores <- c(PopTolScores, scores[TolBs])
        }
        if (elitism > 0) {
            Pop <- rbind(Pop3, Pop[(PopSize - elitism + 1):PopSize, 
                ])
        }
        else {
            Pop <- Pop3
        }
    }
    PopTol <- PopTol[-1, ]
    PopTolScores <- PopTolScores[-1]
    TolBs <- which(PopTolScores < scores[length(scores)] + tolScore)
    PopTol <- PopTol[TolBs, ]
    PopTolScores <- PopTolScores[TolBs]
    PopTolT <- cbind(PopTol, PopTolScores)
    PopTolT <- unique(PopTolT, MARGIN = 1)
    if (!is.null(dim(PopTolT))) {
        PopTol <- PopTolT[, 1:(dim(PopTolT)[2] - 1)]
        PopTolScores <- PopTolT[, dim(PopTolT)[2]]
    }
    else {
        PopTol <- PopTolT[1:(length(PopTolT) - 1)]
        PopTolScores <- PopTolT[length(PopTolT)]
    }
    res <- res[3:dim(res)[1], ]
    rownames(res) <- NULL
    return(list(bString = bestbit, Results = res, StringsTol = PopTol, 
        StringsTolScores = PopTolScores))
}
#=====================================================================================

#	cutSimList3
#=====================================================================================
#=====================================================================================

cutSimList3<-function (SimList, bitString) 
{
    bitString<-as.logical(bitString)
    SimListCut <- SimList
    finalCube <- SimListCut$finalCube[bitString, ]
    ixNeg <- SimListCut$ixNeg[bitString, ]
    ignoreCube <- SimListCut$ignoreCube[bitString, ]
    maxIx <- SimListCut$maxIx[bitString]
    if (is.matrix(finalCube) == FALSE) {
        SimListCut$finalCube <- matrix(finalCube, dimnames = list(names(finalCube), 
            1))
        SimListCut$ixNeg <- matrix(ixNeg, dimnames = list(names(ixNeg), 
            1))
        SimListCut$ignoreCube <- matrix(ignoreCube, dimnames = list(names(ignoreCube), 
            1))
        SimListCut$maxIx <- matrix(maxIx, dimnames = list(names(maxIx), 
            1))
    }
    else {
        SimListCut$finalCube <- finalCube
        SimListCut$ixNeg <- ixNeg
        SimListCut$ignoreCube <- ignoreCube
        SimListCut$maxIx <- maxIx
    }
    return(SimListCut)
}

#=====================================================================================

#	dupsBetweenGroups
#=====================================================================================
#=====================================================================================

# From Cookbook R

dupsBetweenGroups <- function (df, idcol) {
    datacols <- setdiff(names(df), idcol)
    sortorder <- do.call(order, df)
    df <- df[sortorder,]
    dupWithin <- duplicated(df)
    dupBetween = rep(NA, nrow(df))
    dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols])
    dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols], fromLast=TRUE) | dupBetween[!dupWithin]
    goodIdx <- !is.na(dupBetween)
    goodVals <- c(NA, dupBetween[goodIdx])
    fillIdx <- cumsum(goodIdx)+1
    dupBetween <- goodVals[fillIdx]
    dupBetween[sortorder] <- dupBetween
    return(dupBetween)
}


#=====================================================================================

#	getFit3
#=====================================================================================
#=====================================================================================
getFit3<-function (SimResults, CNOlist, Model, indexList, timePoint = c("t1", 
    "t2"), sizeFac = 1e-04, NAFac = 1, SimResultsT0 = NA,
    nInTot = nInTot,nInStimTot = nInStimTot,nInDATot = nInDATot,
    StimuliFac=0.001,DAFac=10) 
{
    SimResults <- SimResults[, indexList$signals]
    nSignals<-dim(SimResults)[2]
    if (timePoint == "t1") 
        tPt <- 2
    Diff <- SimResults - CNOlist$valueSignals[[tPt]]
    r <- Diff^2
    deviationPen <- sum(r[!is.na(r)])
    nDataP <- sum(!is.na(CNOlist$valueSignals[[2]]))
    NAFac<-0.1*nDataP*NAFac
    NAPen <- NAFac * sum(is.na(SimResults))
    nDataPts <- dim(CNOlist$valueSignals[[tPt]])[1] * dim(CNOlist$valueSignals[[tPt]])[2]
    nInputs <- length(which(Model$interMat == -1))
    Stimrows<-match(CNOlist$namesStimuli,rownames(Model$interMat))
    DArows<-ifelse(DAFac==0,0,match("DA",rownames(Model$interMat)))
    if (dim(Model$interMat)[2]==1||is.null(dim(Model$interMat))){
	DAPen<-ifelse(is.na(DArows)||DAFac==0,0,1)
	StimuliPen<-ifelse(all(is.na(Stimrows)),0,length(Stimrows))-DAPen
      StimuliPen<-StimuliPen*StimuliFac
	DAPen<-DAFac*sizeFac*DAPen
    } else {
    	nDAInputs<-ifelse(DAFac==0,0,length(which(Model$interMat[DArows,]==-1)))
    	nStiInputs<-length(which(Model$interMat[Stimrows,]==-1))-nDAInputs
	DAPen<-DAFac*sizeFac*nDAInputs
    	StimuliPen<-StimuliFac*nStiInputs
    }
    deviationPen<-deviationPen/nDataPts
    nInDATot<-ifelse(DAFac==0,1,nInDATot)
    sizePen <- (sizeFac * nInputs)/nInTot
    StimuliPen <- StimuliPen/nInTot
    DAPen <- DAPen/nInTot
    NAPen<-NAPen/nDataPts
    score <- deviationPen + NAPen + sizePen+StimuliPen+DAPen
    return(score)
}

#=====================================================================================

# findNP
#=====================================================================================
#=====================================================================================

findNP<-function(CNOlist,Model,verbose=FALSE){

	NPCues<-CNOlist$namesCues[which(is.na(match(CNOlist$namesCues,Model$namesSpecies)))]	
	NPSignals<-CNOlist$namesSignals[which(is.na(match(CNOlist$namesSignals,Model$namesSpecies)))]
	if (verbose==TRUE){
		print(paste("The following Stimulus and/or Inhibitors are not present in the model:",
			NPCues,sep=" "))
		print(paste("The following Signals are not present in the model:",
			NPSignals,sep=" "))
	}
	NPSpecies<-list(NPCues=NPCues,NPSignals=NPSignals)
	return(NPSpecies)
}


# cutNP
#=====================================================================================
#=====================================================================================

cutNP<-function(CNOlist,valueSignalsIniSim,NPcut){
	CNOlist0<-CNOlist
      valueSignalsIniSim0<-valueSignalsIniSim
	NPCues<-NPcut$NPCues
	NPSignals<-NPcut$NPSignals
	if(length(NPCues)!=0 || length(NPSignals)!=0){
		CNOlist$namesCues<-CNOlist0$namesCues[is.na(match(CNOlist0$namesCues,NPCues))]
		CNOlist$namesStimuli<-CNOlist0$namesStimuli[is.na(match(CNOlist0$namesStimuli,NPCues))]
		CNOlist$namesInhibitors<-CNOlist0$namesInhibitors[is.na(match(CNOlist0$namesInhibitors,NPCues))]
		CNOlist$namesSignals<-CNOlist0$namesSignals[is.na(match(CNOlist0$namesSignals,NPSignals))]
		CNOlist$valueCues<-CNOlist0$valueCues[,is.na(match(CNOlist0$namesCues,NPCues))]
		CNOlist$valueStimuli<-CNOlist0$valueStimuli[,is.na(match(CNOlist0$namesStimuli,NPCues))]
		CNOlist$valueInhibitors<-CNOlist0$valueInhibitors[,is.na(match(CNOlist0$namesInhibitors,NPCues))]
		CNOlist$valueSignals$t0<-CNOlist0$valueSignals$t0[,is.na(match(CNOlist0$namesSignals,NPSignals))]
		CNOlist$valueSignals[[2]]<-CNOlist0$valueSignals[[2]][,is.na(match(CNOlist0$namesSignals,NPSignals))]
            valueSignalsIniSim<-valueSignalsIniSim0[,is.na(match(colnames(valueSignalsIniSim0),NPSignals))]
      } else {
		CNOlist<-CNOlist0
		valueSignalsIniSim<-valueSignalsIniSim0
	}
	return(tempcut<-list(CNOlist=CNOlist,valueSignalsIniSim=valueSignalsIniSim))
}

# txtStart.2wd
#=====================================================================================
#=====================================================================================

txtStart.2wd <- function(...)
{
	if(!require(rcom)) {
		print("I couldn't find rcom on your system - trying to install it now")
		install.packages("rcom")
	}
	if(!require(R2wd)) {
		print("I couldn't find R2wd on your system - trying to install it now")
		install.packages("R2wd")
	}
	if(!require(TeachingDemos)) {
		print("I couldn't find TeachingDemos on your system - trying to install it now")
		install.packages("TeachingDemos")
	}	
	my.file <- "temp_ToWord.txt"
	txtStart(my.file,...) # like sink - but catched the commands from the console as well
}

# txtStop.2wd
#=====================================================================================
#=====================================================================================

txtStop.2wd <- function(open.new.file = T, save.word.file = T)
{
	my.file <- "temp_ToWord.txt"
	# The following function is needed to allow the reading of the output of the console into word
	wdBody.anything <- function(output)
	{
		# This function takes the output of an object and prints it line by line into the word document
		# Notice that in many cases you will need to change the text font into courier new roman...
		if(class(output) == "character")
		{
			a <- output
		} else {
			a <- capture.output(output)
		}
		for(i in seq_along(a))
		{
			wdBody(format(a[i]))			
			#wdBody('\n')
		}
	}

      txtStop() # Stop the connection
	console.output <- readLines(con = my.file)
	unlink(my.file)	# erases the temp file that was created to store the console output.
	if(open.new.file) {
		wdGet(T)	# If no word file is open, it will start a new one - can set if to have the file visiable or not
		wdNewDoc("this.doc")	# this creates a new file with "this.doc" name
		}
	wdBody.anything(console.output)	# puts out output into the word file

	if(save.word.file)
	{
		wdSave("This.doc")
	}
}
#=====================================================================================

# getScaffold2
#=====================================================================================
#=====================================================================================

getScaffold2<-function (modelComprExpanded, optimResT1, optimResT2, modelOriginal,boolean) 
{
    bString1 <- optimResT1$bString
	if (!boolean){
		L1<-optimResT1$L
	}
    if (is.na(optimResT2[1])) {
        bString2 <- optimResT1$bString[which(optimResT1$bString == 0)]
		if (!boolean){
			L2 <- optimResT1$L[which(optimResT1$bString == 0)]
		}
    } else {
        bString2 <- optimResT2$bString
		L2<-optimResT2$L
    }
    BStimes <- bString1
    BStimes[which(BStimes == 0)] <- bString2 * 2
	Ltimes <- L1
	Ltimes[which(BStimes == 0)]<-L2
    if (!is.null(dim(optimResT1$stringsTol))) {
        bW1 <- apply(optimResT1$stringsTol, 2, mean)
    }
    else {
        bW1 <- bString1
    }
    if (!is.na(optimResT2[1])) {
        if (!is.null(dim(optimResT2$stringsTol))) {
            bW2 <- apply(optimResT2$stringsTol, 2, mean)
        }
        else {
            bW2 <- bString2
        }
        weightsE <- bW1
        weightsE[which(optimResT1$bString == 0)] <- weightsE[which(optimResT1$bString == 
            0)] + bW2
    } else {
        weightsE <- bW1
    }
    findOutput <- function(x) {
        sp <- which(x == 1)
        sp <- modelComprExpanded$namesSpecies[sp]
    }
    reacOutput <- apply(modelComprExpanded$interMat, 2, findOutput)
    findInput <- function(x) {
        sp <- which(x == -1)
		if (length(sp)==0){ #to account for self-loops
			sp<-which(x==1)
		}
        sp <- modelComprExpanded$namesSpecies[sp]
    }
    reacInput <- apply(modelComprExpanded$interMat, 2, findInput)
    createReac <- function(x) {
        r <- paste(x[1], " (", x[2], ") ", x[3], sep = "")
        return(r)
    }
    if (class(reacInput) != "list") {
        isNeg <- function(x) {
            isNegI <- any(x == 1)
            return(isNegI)
        }
        inpSign <- apply(modelComprExpanded$notMat, 2, isNeg)
        inpSign <- !inpSign
        inpSign[inpSign] <- 1
        inpSign[!inpSign] <- -1
        sifFile <- cbind(reacInput, inpSign, reacOutput)
        EApresent <- apply(sifFile, 1, createReac)
        EApresent <- cbind(EApresent, BStimes)
        EAweights <- cbind(EApresent, weightsE)
        sifFile <- cbind(sifFile, BStimes)
        sifFile <- cbind(sifFile, weightsE)
    } else {
	    if (boolean){
			sifFile <- matrix(0, nrow = 4 * length(reacOutput), ncol = 5)
		} else {
			sifFile <- matrix(0, nrow = 4 * length(reacOutput), ncol = 6)
		}
        nR <- 1
        nANDs <- 1
        for (i in 1:length(reacOutput)) {
            if (length(reacInput[[i]]) == 1) {
                sifFile[nR, 1] <- reacInput[[i]]
                sifFile[nR, 3] <- reacOutput[i]
                sifFile[nR, 2] <- ifelse(any(modelComprExpanded$notMat[, 
                  i] == 1), -1, 1)
                sifFile[nR, 4] <- BStimes[i]
                sifFile[nR, 5] <- weightsE[i]
				if (!boolean){
					sifFile[nR,6]<-Ltimes[i]
				}
                nR <- nR + 1
            } else {
                for (inp in 1:length(reacInput[[i]])) {
                  sifFile[nR, 1] <- reacInput[[i]][inp]
                  sifFile[nR, 3] <- paste("and", nANDs, sep = "")
                  temp_indices = which(reacInput[[i]][inp] == 
                    rownames(modelComprExpanded$notMat))
                  sifFile[nR, 2] <- ifelse(modelComprExpanded$notMat[temp_indices, 
                    i] == 1, -1, 1)
                  sifFile[nR, 4] <- BStimes[i]
                  sifFile[nR, 5] <- weightsE[i]
				  if(!boolean){
					sifFile[nR, 6] <-Ltimes[i]
				  }
                  nR <- nR + 1
                }
                sifFile[nR, 1] <- paste("and", nANDs, sep = "")
                sifFile[nR, 3] <- reacOutput[i]
                sifFile[nR, 2] <- 1
                sifFile[nR, 4] <- BStimes[i]
                sifFile[nR, 5] <- weightsE[i]
				if (!boolean){
					sifFile[nR, 6] <-Ltimes[i]
				}
                nANDs <- nANDs + 1
                nR <- nR + 1
            }
        }
        sifFile <- sifFile[1:(nR - 1), ]
        EApresent <- apply(sifFile[, 1:3], 1, createReac)
        EAweights <- cbind(EApresent, sifFile[, 5])
        EApresent <- cbind(EApresent, sifFile[, 4])
		if (!boolean){
			Lpresent<-cbind(EApresent,sifFile[,6])
		}
	}
    makeEA <- function(x) {
        ea <- paste(x[1], "=", x[2])
        return(ea)
    }
    EApresent <- apply(EApresent, 1, makeEA)
    EAweights <- apply(EAweights, 1, makeEA)
	Lpresent <- apply(Lpresent, 1, makeEA)
	if (boolean){
		return(list(EApresent = EApresent, EAweights = EAweights, 
			sifFile = sifFile))
	} else {
		return(list(EApresent = EApresent, EAweights = EAweights, Lpresent=Lpresent,
			sifFile = sifFile))
	}
}



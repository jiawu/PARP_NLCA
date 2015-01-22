

	rm(list=ls(all=TRUE))

# START COPYING IN TEXT FILE

#	Please indicate the location of the library

	lib.loc="~/R/library"
	#lib.loc="~/R/x86_64-redhat-linux-gnu-library/2.15" #feldspar library location
    #lib.loc="/home/bpe372/usr/lib/R/library" # sapphire location

# Jia: Changed the order for generating temporary file names because this current method doesn't work on quest.

	directory.cluster ="~/PARP_SHORT_MTWT"
	directory = "~/PARP_SHORT_MTWT"



# 	Do not modify 	

	library("TeachingDemos",lib.loc=lib.loc)
	currtime<-paste(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),sep="")
	tempname = tempfile(pattern = paste(currtime,"ver",sep=""), tmpdir = "/projects/p20519/jia_output/PARP_SHORT_MTWT")
	#txtStart(paste(tempname,".txt",sep=""))
	date()
	Sys.Date()
	Sys.time()

#########################################################################################################################################################
####									INPUT DATA														#########
#########################################################################################################################################################

# PLEASE READ IN DETAIL AND FILL UP ALL THE REQUIRED INFORMATION. SEE MANUAL FOR MORE INFORMATION
					
#	Point this to the directory containing all of the data (if you are running in the cluster
#	or using a Mac, please fill the directory.cluster name)


	
# 	Runing in cluster or in a Mac? Please indicate 1 for cluster or if you are running in a Mac.
# 	In that case all the graphs would be saved as pdf.

	clstr<-1
	
#   Indicate the name of the output folder

	folder.output<-c(tempname)
	
#	Indicate the location and name of the source file

	source.location<-"~/Rprofilev3_8_feld.R"

#	Name of NO DNA Row

	name_no_dna <- "NC"

#	Name of Control Row

	name_TA <- "TA"

#	Name of Untreated "Treatment" or BaseLine Treatment

	name_untreated <- "BRCA_WT"

#	Names of your factors

	factornames<-c("Cell_Type")
	
#	Time vector

	timevector<-c(0,3,6,9,24)

#	File that contains the normalized data

	NormalizedData <-"DATA_TABLE_MT2.txt"

#	File that contains the experimental design. The format of the columns should be the following
#		Column1=Number for the condition
#		Column2=TF associated with it
#		Column3=Factor1
#		Column4=Factor2
#		..............
#		Column(N+2)=FactorN

	ExperimentalDesign<-c("BRCA_legend_MT.txt") 

#	File name of the TF that are significant 

	Significant<-"short_mt_sig2.txt"

# 	Are the factors qualitive?

	qfactors<-c(1)

# 	Indicate the levels that can be used as a based line for each of the factors

	baseline<-c("BRCA_WT")

# 	Indicate the name of the prior knowledge network file

	#PKN<-"~/master_PKN1.txt"
	PKN<-"~/PKN_2.txt"

# 	Indicate the name of the CNO file,if available

	PKN_INPUTS<-"~/empty_pkn2.txt"

# 	Indicate the name of the conversion file

	ConversionList<-"~/master_conversion_list.txt"

# 	Do you want your data to be normalized by the control?

	Controlnormalization=1

# PLSR/TD-PLSR paramteres
	
	PLSRType="PLSR" # You can select PLSR or TDPLSR
	PLSRCutoffInput=0.15
    PLSRCutoffTF=0.3

# 	Mutual information parameters
# 	Assumptions for the calculations. Do you want to assume the data to be static for the 
# 	calculations ("STAT"), dynamic ("DYN") or both ("BOTH")?

	MIMdata<-"DYN"

# 	How manay interpolation timepoints would you consider for the initial mutual information method calculation?

	tinitial=3

#	CNO parameters
# 		Account for Deactivation
		DAPresent<-TRUE
# 		expand model
		expandModel=FALSE
# 		Bootstrapping
		bootstrapping=TRUE
#		Permutation
		permutation=FALSE
# 		Population Size
		PopSize=15
# 		Maximum Time
		#MaxTime=360000 (roughly 4 days)
    #maxtime 172800 is 48 hours
		MaxTime= 172800
# 		Maximum number of generations
		#maxGens=8000
    maxGens = 3000
    # 		Maximum number of identical generations
		#StallGenMax=8000
    StallGenMax = 800
# 		Percentaje of non present edges for random start (0 to 1) or "RAND"
		randedgeper=0.1
# 		Elitism
		elitism=1
# 		Mutations
		Pmutation=0.005
# 		Selection Pressure
		SelPress=21.5
# 		sizeFac
		sizeFac=0.001
# 		Stimuli Fac
		StimuliFac=0.005
# 		Deactivation Mechanism Factor Penalty. If DAPresent=FALSE
# 		set it to 0.
		DAFac=50
# 		Increase policy over time for Stimuli Fac
		incPolStim=10
# 		Initial String
		initStringUser<-c()
			#read.table("~/Research/Micheal Wess/Erb2/Activation/FinalData/finalBString.txt", header = TRUE, sep = "\t", 
			#stringsAsFactors = default.stringsAsFactors(),
			#row.names=NULL,check.names=FALSE)


##################################################################################################################
#####									DATA ANALYSIS													##########
##################################################################################################################

	inputParameters<-list(
		directory.cluster = directory.cluster,
		directory = directory,
		clstr=clstr,
		folder.output=folder.output,
		source.location=source.location,
		name_no_dna =name_no_dna,
		name_TA =name_TA,
		name_untreated= name_untreated,
		factornames=factornames,
		timevector=timevector,
		NormalizedData =NormalizedData,
		ExperimentalDesign=ExperimentalDesign,
		Significant=Significant,
		qfactors=qfactors,
		baseline=baseline,
		PKN=PKN,
		PKN_INPUTS=PKN_INPUTS,
		ConversionList=ConversionList,
		Controlnormalization=Controlnormalization,
		PLSRType=PLSRType, 
		PLSRCutoffInput=PLSRCutoffInput,
		PLSRCutoffTF=PLSRCutoffTF,
		MIMdata=MIMdata,
		tinitial=tinitial,
		DAPresent=DAPresent,
		expandModel=expandModel,
		bootstrapping=bootstrapping,
		permutation=permutation,
		PopSize=PopSize,
		MaxTime=MaxTime,
		maxGens=maxGens,
		StallGenMax=StallGenMax,
		randedgeper=randedgeper,
		elitism=elitism,
		Pmutation=Pmutation,
		SelPress=SelPress,
		sizeFac=sizeFac,
		StimuliFac=StimuliFac,
		DAFac=DAFac,
		incPolStim=incPolStim,
		initStringUser=initStringUser)
		
#	DATA PREPROCESSING
##########################################################################################################################################################

	source(source.location)
	if (clstr==0){
			setwd(directory); getwd()
		} else {
			setwd(directory.cluster); getwd()
			directory<-directory.cluster
		}
	old.par <- par(no.readonly = TRUE)
	
	dataAnalysis<-preProcessing(inputParameters=inputParameters)
	dataInt<-wrappingSpline(inputParameters=inputParameters,dataAnalysis=dataAnalysis)
#	wrappingkmeans(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
#		dataInt=dataInt)

#	INFERENCE METHODS AND PRIOR KNOWLEDGE
##########################################################################################################################################################

#	Time-delayed Linear Regression
	TDLRRes<-linearInference(inputParameters=inputParameters,dataAnalysis=dataAnalysis, dataInt=dataInt)
	dataInt<-c(dataInt,TDLRRes$dataInt)

# 	PLSR
	PLSRRes<-PLSRInference(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt)
	dataInt<-c(dataInt,PLSRRes$dataInt)
	InferenceRes<-list(NET_TDPLSR=PLSRRes$NET_TDPLSR)

#	Mutual Information	
	MIRes<-MIInference(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt)
	dataInt<-c(dataInt,MIRes$dataInt)
	InferenceRes[["NET_MIM"]]<-MIRes$NET_MIM

# 	BANJO
	BANJORes<-BANJOInference(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt)
	PKNList<-BANJORes
	# IF YOUR ARE RUNNING FOR THE FIRST TIME, JUST RUN UP TO HERE. ONCE YOU HAVE THE NET_BANJO.txt FILE
    # YOU CAN RUN THE ENTIRE CODE	
	NET_BANJO<-read.table("NET_BANJO_mar2014.txt", header = TRUE, sep = "\t", colClasses = "character")
	InferenceRes[["NET_BANJO"]]<-NET_BANJO
	
#	Prior Knowledge
	PKNSimp<-PKNProcessing(inputParameters=inputParameters,dataAnalysis=dataAnalysis,
		dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes)
	PKNList<-PKNSimp$PKNList
	InferenceRes<-PKNSimp$InferenceRes

#   OPTIMIZATION
##########################################################################################################################################################
	
	#CNO data formating
	dataCNOpre<-preCNO(inputParameters=inputParameters,dataAnalysis=dataAnalysis,dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes)
	PKNList<-dataCNOpre$PKNList
	
	#Running GA and network generation
	dataCNOOpt<-wrappingCNO(inputParameters=inputParameters,dataAnalysis=dataAnalysis,dataInt=dataInt,PKNList=PKNList,InferenceRes=InferenceRes,dataCNOpre=dataCNOpre)
	
	results_location <- paste(dataAnalysis$folder.output, "/","NLCA_cellnetopt_results.rdata",sep="")
	save(dataCNOOpt, file =results_location)

# - - - - - - - - - -- - - - - - - - - - - - - Settings - - - - - - - - - -- - - - - - - - - - - - - 

# File handling
fileName = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/hargrove-OD-22.10.2010.csv"
fileWriteLocation = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/harg-22-10/"
SASFilePref = "LMM_var_S"
SASOutFileCSVnamePref = "LMM_SAS_OUT"

# algorithm variables
ODcut = 0.8
names <- array(c(3,4,5,6, 7))
timeCutLower = 0 # exclude samples for which first BED test is lower than timeCut
timeCutUpper = 140 # exclude samples for which first BED test is higher than timeCut
ExcludeAlreadySeroconverted = TRUE
ExcludeNeverSeroconverted = TRUE
CovMatrix = 1

#SAS specs
SasOpeningString = "data zvit_mixed; infile datalines delimiter=' '; input ID lnTime sqOD @@; datalines;"
SasClosingString = ";proc mixed data=zvit_mixed;model sqOD=lnTime / solution ddfm=kr outpred=pred;repeated / subject=ID r rcorr;random intercept lnTime /type=un subject=ID solution;ods listing exclude solutionr;ods output solutionr=out;run;PROC EXPORT DATA= WORK.Out OUTFILE= \"C:\\Documents and Settings\\heastwood\\My Documents\\SACEMA\\OD cut-off\\harg-22-10\\"
SasClosingString2 = ";proc mixed data=zvit_mixed;model sqOD=lnTime / solution ddfm=kr outpred=pred;repeated / type=sp(pow)(lnTime)  subject=ID r rcorr;random intercept lnTime /type=un subject=ID solution;ods listing exclude solutionr;ods output solutionr=out;run;PROC EXPORT DATA= WORK.Out OUTFILE= \"C:\\Documents and Settings\\heastwood\\My Documents\\SACEMA\\OD cut-off\\harg-22-10\\"

# - - - - - - - - - -- - - - - - - - - - - - - End settings - - - - - - - - - -- - - - - - - - - - - - - 



# Functions
make_scenario <- function(ID.time.OD, cut, nSamples, timeCutLower, timeCutUpper, ExcludeAlreadySeroconverted, ExcludeNeverSeroconverted) {
	# settings                               
	# cut = OD cut off
	# dataNNu.time.OD = data input matrix (person Id, time of test in days) 
	# setting: 										1 = remove already seroconverted at t0, 
	# 												2 = remove slow progressors
	#												3 = remove normal progressors
	#												4 = remove already seroconverted and slow progressors
	
	BEDH <- ID.time.OD
	d<-as.data.frame(BEDH)
	n <- dim(BEDH)[1] - 1
	nnu = d$ID[n]
	omitArray = array(dim=c(nnu,1))
	dataOut <- array(dim=c(nnu,4))
	
	d2 = d
	count <- 0

	# tag already seroconverted at first and subsequent measurements
	if (ExcludeAlreadySeroconverted == TRUE){
		for (j in c(1:nnu)) {
			if (dim(d[d$ID == j & d$OD > cut,])[1] == dim(d[d$ID == j,])[1]){
					d2[d2$ID == j,][2] <- 11111
					d2[d2$ID == j,][3] <- 11111
			}
		}
	}
	
	#tag those who never seroconvert in any measurements
	if (ExcludeNeverSeroconverted == TRUE) {
		for (j in c(1:nnu)) {
			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0)) {
					d2[d2$ID == j,][2] <- 22222
					d2[d2$ID == j,][3] <- 22222
			}
		}
	}

	#tag only n samples
	for (j in c(1:nnu)) {
		if (dim(d[d$ID == j,])[1] < nSamples) {
				d2[d2$ID == j,][2] <- 33333
				d2[d2$ID == j,][3] <- 33333
		}
	}
	
	#tag samples for which the first time is < timeCut
	for (j in c(1:nnu)) {
		if (d[d$ID == j,2][1] <= timeCutLower | d[d$ID == j,2][1] >= timeCutUpper) {
				d2[d2$ID == j,][2] <- 44444
				d2[d2$ID == j,][3] <- 44444
		}
	}
	d3 <- d2[d2$OD != 11111 & d2$OD != 22222 & d2$OD != 33333 & d2$OD != 44444,]
	return(d3)

	}
get_numSamples <- function(dataArrayIn) {
	return(dim(unique(dataArrayIn[1]))[1])
}

if (ExcludeAlreadySeroconverted == TRUE & ExcludeNeverSeroconverted == TRUE) sero = "+"
if (ExcludeAlreadySeroconverted == TRUE & ExcludeNeverSeroconverted == FALSE) sero = "EAST.ENSF"
if (ExcludeAlreadySeroconverted == FALSE & ExcludeNeverSeroconverted == TRUE) sero = "EASF.ENST"
if (ExcludeAlreadySeroconverted == FALSE & ExcludeNeverSeroconverted == FALSE) sero = "EASF.ENSF"

description = paste("_const_OD.", OD, "_TL.", timeCutLower, "_TU.", timeCutUpper, "_",sero, "_COV.", CovMatrix, sep="")
SASFile = paste(SASFilePref, description, sep="")
SASOutFileCSVname = paste(SASOutFileCSVnamePref, description, sep="")

nSampleEnd = dim(names)
fileWriteLocation2 = paste(fileWriteLocation, SASFile, sep="")
fileWriteLocation3 = paste(fileWriteLocation, SASOutFileCSVname, sep="")
BEDH <- read.csv(file=fileName, head=TRUE, sep=" ")
Sdata = list(); Strans = list(); Sunique = list(); Sunique = list(); N <- list()

# write sas start
if (CovMatrix == 1) SasClose = SasClosingString
if (CovMatrix == 2) SasClose = SasClosingString2

#write names and S values to array
write.table(names, file = paste(fileWriteLocation3, ".names", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)

for (i in c(1:nSampleEnd)) {
	Sdata[[i]] <- make_scenario(BEDH, ODcut, i, timeCutLower, timeCutUpper, ExcludeAlreadySeroconverted, ExcludeNeverSeroconverted)
	Strans[[i]] <- Sdata[[i]]
	Strans[[i]][2] = log(Sdata[[i]][2])
	Strans[[i]][3] = sqrt(Sdata[[i]][3])
	Sunique[[i]] = unique(Sdata[[i]][1])
	N[[i]] <- get_numSamples(Sdata[[i]])
	print(paste ("number of samples for ", names[i], "is " , N[[i]]))

	
	# write sas file
		# write sas file
	write.table(SasOpeningString, file = paste(fileWriteLocation2, ".S", names[i], ".SAS", sep=""), sep = " ", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
	write.table(Strans[[i]], file = paste(fileWriteLocation2, ".S", names[i], ".SAS", sep=""), sep = " ", col.names = FALSE, append="TRUE", row.names=FALSE)
	SasClose2 = paste(SasClose, SASOutFileCSVname, ".S", names[i],".csv\" DBMS=CSV REPLACE;RUN;", sep="")
	write.table(SasClose2, file = paste(fileWriteLocation2, ".S", names[i], ".SAS", sep=""), sep = " ", col.names = FALSE, append="TRUE", row.names=FALSE, quote=FALSE)
}

 # - - - - - - - - - -- - - - - - - - - - - - - End code - - - - - - - - - -- - - - - - - - - - - - - 
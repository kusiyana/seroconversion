# - - - - - - - - - -- - - - - - - - - - - - - Settings - - - - - - - - - -- - - - - - - - - - - - - 

# File handling
fileName = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/hargrove-OD-22.10.2010.csv"
fileWriteLocation = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/harg-22-10/"
fileWriteName = "test"

# algorithm variables
cut = 0.8
nSamples = 1 # number of samples to include
timeCutLower = 0 # exclude samples for which first BED test is lower than timeCut
#timeCutUpper = 100 # exclude samples for which first BED test is higher than timeCut
ExcludeAlreadySeroconverted = FALSE	
ExcludeNeverSeroconverted = FALSE

loopLimit = 16 #Set loop if you want to loop through different values and generate data for each setting

#SAS specs
SasOpeningString = "data zvit_mixed; infile datalines delimiter=' '; input ID lnTime sqOD @@; datalines;"
SasClosingString = ";proc mixed data=zvit_mixed;model sqOD=lnTime / solution ddfm=kr outpred=pred;repeated / subject=ID r rcorr;random intercept lnTime /type=un subject=ID solution;ods listing exclude solutionr;ods output solutionr=out;run;PROC EXPORT DATA= WORK.Out OUTFILE= \"C:\\Documents and Settings\\heastwood\\My Documents\\SACEMA\\OD cut-off\\harg-22-10\\"

names <- c("t60", "t70", "t80", "t90", "t100", "t110", "t120", "t130", "t140", "t150", "t160", "t170", "t180", "t190", "t200", "tall")
tvalues <- c(60, 70, 80, 90, 100, 110, 120, 130, 140,150,160, 170, 180, 190, 200, 2500)


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


BEDH <- read.csv(file=fileName, head=TRUE, sep=" ")
Sdata = list(); Strans = list(); Sunique = list(); Sunique = list()

# write sas start
for (i in c(1:loopLimit)) {
	Sdata[[i]] <- make_scenario(BEDH, cut, nSamples, timeCutLower, tvalues[i], ExcludeAlreadySeroconverted, ExcludeNeverSeroconverted)
	Strans[[i]] <- Sdata[[i]]
	Strans[[i]][2] = log(Sdata[[i]][2])
	Strans[[i]][3] = sqrt(Sdata[[i]][3])
	Sunique[[i]] = unique(Sdata[[i]][1])
	numSamples[[i]] <- get_numSamples(Sdata[[i]])
	print(paste ("number of samples for ", names[i], "is " , numSamples[[i]]))
	
	#write csv files
	write.table(Sdata[[i]], file = paste(fileWriteLocation, names[i], ".csv"), sep = " ", col.names = TRUE, append="FALSE", row.names=FALSE)
	write.table(Strans[[i]], file = paste(fileWriteLocation, names[i], "-trans.csv"), sep = " ", col.names = TRUE, append="FALSE", row.names=FALSE)
	
	# write sas file
	write.table(SasOpeningString, file = paste(fileWriteLocation, names[i], ".SAS"), sep = " ", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
	write.table(Strans[[i]], file = paste(fileWriteLocation, names[i], ".SAS"), sep = " ", col.names = FALSE, append="TRUE", row.names=FALSE)
	SasClosingString2 = paste(SasClosingString, "SAS-rand-",names[i],".csv\" DBMS=CSV REPLACE;RUN;")
	write.table(SasClosingString2, file = paste(fileWriteLocation, names[i], ".SAS"), sep = " ", col.names = FALSE, append="TRUE", row.names=FALSE, quote=FALSE)
}

 # - - - - - - - - - -- - - - - - - - - - - - - End code - - - - - - - - - -- - - - - - - - - - - - - 
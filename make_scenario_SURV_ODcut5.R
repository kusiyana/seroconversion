# - - - - - - - - - -- - - - - - - - - - - - - Settings - - - - - - - - - -- - - - - - - - - - - - - 

# File handling
fileName = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/hargrove-OD-22.10.2010.csv"
fileWriteLocation = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/harg-22-10/"
SASFilePref = "SURV_04.08_var_OD"

# algorithm variables
Wmax = 365
ODStart = 0.3
ODStop = 3.0
ODStep = 0.1
nSamples = 4 # ie at least n samples from a client
timeCutLower = 0 # exclude samples for which first BED test is lower than timeCut
timeCutUpper = 120 # exclude samples for which first BED test is higher than timeCut
ExcludeAlreadySeroconverted = FALSE
ExcludeNeverSeroconverted = FALSE
Wmax = 365
# - - - - - - - - - -- - - - - - - - - - - - - End settings - - - - - - - - - -- - - - - - - - - - - - - 



# Functions
make_scenario <- function(ID.time.OD, cut, nSamples, timeCutLower, timeCutUpper, ExcludeAlreadySeroconverted, ExcludeNeverSeroconverted, Wmax) {
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
	countWmaxMat = array();
	countWmax = 0
	d2 = d
	count <- 0
	timeDiff <- array();
	
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
			lastIndex <- dim(d[d$ID == j,])[1]
			timeDiff[j] <- d[d$ID == j,][lastIndex,2] - d[d$ID == j,][1,2]
			
			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & timeDiff[j] < Wmax) {
					d2[d2$ID == j,][2] <- 22222
					d2[d2$ID == j,][3] <- 22222
			}
#			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & timeDiff[j] > Wmax) {
#				countWmax = countWmax + 1
#			}
		}
	countWmaxMat[j] = countWmax	
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

get_recents <- function(Sdata, cut) {
	d<-as.data.frame(Sdata)
	n <- dim(d)[1] - 1
	nnu = d$ID[n]
	count = 0
	continue = TRUE
	for (j in c(1:nnu)) {
		if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & continue == TRUE) {
			continue = FALSE
			count = count + 1
		}
		else {
			continue = TRUE
		}
	}
	return(count)
}

get_alreadyAT <- function(Sdata, cut) {
	d<-as.data.frame(Sdata)
	n <- dim(d)[1] - 1
	nnu = d$ID[n]
	count = 0
	continue = TRUE
	ATArray = d # array(0, c(nnu,1))
	for (j in c(1:nnu)) {
			if (dim(d[d$ID == j & d$OD > cut,])[1] == dim(d[d$ID == j,])[1]){
					ATArray[ATArray$ID == j,][2] <- 11111
					ATArray[ATArray$ID == j,][3] <- 11111
			}
	}
	out <- as.data.frame(ATArray)
			numbers <- which(out$OD ==11111)
		return(unique(out$ID[numbers]))
}

get_neverRT <- function(Sdata, cut, timeDiff) {
	d<-as.data.frame(Sdata)
	n <- dim(d)[1] - 1
	nnu = d$ID[n]
	count = 0
	continue = TRUE
	d2 = d # array(0, c(nnu,1))
		for (j in c(1:nnu)) {
			lastIndex <- dim(d[d$ID == j,])[1]
			timeDiff[j] <- d[d$ID == j,][lastIndex,2] - d[d$ID == j,][1,2]
			
			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & timeDiff[j] < Wmax) {
					d2[d2$ID == j,][2] <- 22222
					d2[d2$ID == j,][3] <- 22222
			}
		}
	out <- as.data.frame(d2)
			numbers <- which(out$OD ==22222)
		return(unique(out$ID[numbers]))
}

get_numWmaxCut <- function(Sdata, Wmax, cut) {
	d<-as.data.frame(Sdata)
	n <- dim(d)[1] - 1
	nnu = d$ID[n]
	timeDiff = array()
	maxCutArray = array(0, c(nnu,1))
		for (j in c(1:nnu)) {
			lastIndex <- dim(d[d$ID == j,])[1]
			timeDiff[j] <- d[d$ID == j,][lastIndex,2] - d[d$ID == j,][1,2]
			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & timeDiff[j] > Wmax) {
				maxCutArray[j] = 1
			}
		}
		out <- as.data.frame(maxCutArray)
		return(which(out ==1))
}

get_indexTnt <- function(Sdata, TntL, TntU) {
	d<-as.data.frame(Sdata)
	d4 = d
	n <- dim(d)[1] - 1
	nnu = d$ID[n]
	timeDiff = array()
	maxCutArray = array(0, c(nnu,1))
		for (j in c(1:nnu)) {
		if (d[d$ID == j,2][1] <= TntL | d[d$ID == j,2][1] >= TntU) {
				d4[d4$ID == j,][2] <- 44444
				d4[d4$ID == j,][3] <- 44444
		}
	}
		out <- as.data.frame(d4)
		numbers <- which(out ==44444)
		return(unique(out$ID[numbers]))
}




if (ExcludeAlreadySeroconverted == TRUE & ExcludeNeverSeroconverted == TRUE) sero = "+"
if (ExcludeAlreadySeroconverted == TRUE & ExcludeNeverSeroconverted == FALSE) sero = "EAST.ENSF"
if (ExcludeAlreadySeroconverted == FALSE & ExcludeNeverSeroconverted == TRUE) sero = "EASF.ENST"
if (ExcludeAlreadySeroconverted == FALSE & ExcludeNeverSeroconverted == FALSE) sero = "EASF.ENSF"

description = paste("_const_S.", nSamples, "_TL.", timeCutLower, "_TU.", timeCutUpper, "_",sero, "_", sep="")
SASFile = paste(SASFilePref, description, sep="")
SASOutFileCSVname = paste(SASFilePref, description, sep="")

fileWriteLocation2 = paste(fileWriteLocation, SASFile, sep="")
fileWriteLocation3 = paste(fileWriteLocation, SASOutFileCSVname, sep="")
BEDH <- read.csv(file=fileName, head=TRUE, sep=" ")
Sdata = list(); Strans = list(); Sunique = list(); Sunique = list(); N <- list(); wMaxCut <- list()
nWMaxCut <- array()
nrecents <- array()
names <- array()

#write names and S values to array

#names = array()
nStart = 1
nEnd = ((ODStop - ODStart)/ODStep) + 1
ODcut = ODStart

for (i in c(nStart:nEnd)) {
	wMaxCut[[i]] <- get_numWmaxCut(BEDH, Wmax, ODcut) 
	nWMaxCut[i] <- length(wMaxCut[[i]])
	Sdata[[i]] <- make_scenario(BEDH, ODcut, nSamples, timeCutLower, timeCutUpper, ExcludeAlreadySeroconverted, ExcludeNeverSeroconverted, Wmax)
	Strans[[i]] <- Sdata[[i]]
	Strans[[i]][2] <- log(Sdata[[i]][2])
	Strans[[i]][3] <- sqrt(Sdata[[i]][3])
	Sunique[[i]] <- unique(Sdata[[i]][1])
	N[[i]] <- get_numSamples(Sdata[[i]])
	print(paste ("number of samples for ", ODcut, "is " , N[[i]]))

	#write survival file
	write.table(Sdata[[i]], file = paste(fileWriteLocation2, "OD.", ODcut, ".csv", sep=""), sep = " ", col.names = TRUE, append="FALSE", row.names=FALSE, quote=FALSE)

	#write wMaxCut file

	write.table(wMaxCut[[i]], file = paste(fileWriteLocation2, "OD.", ODcut, ".", Wmax, ".Wmax", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)

	names[i] <- ODcut
	nrecents[i] <- get_recents(Sdata[[i]], ODcut)
	ODcut <- ODcut + ODStep
}

write.table(N, file = paste(fileWriteLocation2, "OD.nsamples", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
write.table(names, file = paste(fileWriteLocation2, "OD.names", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
write.table(nrecents, file = paste(fileWriteLocation2, "OD.nrecent", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
write.table(nWMaxCut, file = paste(fileWriteLocation2, "OD.nWMaxCut", sep=""), sep = ",", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)


# - - - - - - - - - -- - - - - - - - - - - - - End code - - - - - - - - - -- - - - - - - - - - - - - 
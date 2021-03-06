		# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
		
		#               this file processes linear mixed model (LMM) data generated by SAS and calculates the mean window period.
		#               NB this variant of the script keeps OD constant and varies the cut off time in which OD is considered.
		#               means are stored in "mean" list, ie mean[[j]], where j corresponds to index in "names" array.
		
		# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

		# settings
		fileReadLocation = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/harg-22-10/"
		ODL <- 0.0476
		fileFixedArray <- "LMM_28.7_OUT_var_OD__const_S.3_TL.0_TU.140_+_COV.1.out"
		fileRandomArray <- "LMM_28.7_OUT_var_OD__const_S.3_TL.0_TU.140_+_COV.1"
		ODL <- 0.0476 # for the LMM equation

		# <-------------------------------------------------Functions ----------------------------------------------------------------------------------------->
	
		# boot strap
		bootMean <- function (data=window.elk[,3],B=1000,alpha=c(0.025,0.975)) {
		library(bootstrap)
		n <- length(data)
		theta <- mean(data)
		boot.theta <- numeric(B)
		for(i in 1:B){
			boot.data <- data[sample(1:n,replace=TRUE)]
			boot.theta[i] <- mean(boot.data)
		}
				
		#BCa
		BCa <-bcanon(data,nboot=B,theta=function(x)mean(x),alpha=alpha)$confpoint[,2]
		leng.BCa <- max(BCa)-min(BCa)
		vorm.BCa <- (max(BCa)-theta)/(theta-min(BCa))
		
			
		#lower, upper, bootmean, mean
		return (c(BCa[1], BCa[2],mean(boot.theta),theta))
	}
	
	
		# <-------------------------------------------------End functions ----------------------------------------------------------------------------------------->
	
	#create different OD cutof scenarios and store in list
	Sdata <- list(); surv <- list(); avg <- list(); mean <- list(); var <- list()
	CI <- list(); sasOut <- list(); upper <- list(); lower <- list();SD = list(); 
	CoeffVar=list();tempDiff=list(); count365=list(); bothN <- list(); lenArray <- list()
	
	colours = c("red","blue","green","yellow","black","purple", "cyan","magenta", "orange", "pink", "brown", "violet")
	limitArray = array(); CIout = array()
	fixedArray <- read.table(file = paste(fileReadLocation, fileFixedArray, sep=""), head=TRUE, sep=",")
	namesArray <- read.table(file = paste(fileReadLocation, fileRandomArray, ".names", sep=""), head=FALSE, sep=",")
	loopStart <- 1
	loopLimit <- dim(fixedArray)[1]
	for (j in c(loopStart:loopLimit)) {
		sasOut[[j]] <- read.table(file = paste(fileReadLocation, fileRandomArray, "_OD.", namesArray$V1[j], ".csv", sep=""), head=TRUE, sep=",")
		k <- 1
		l <- 1
		uTemp <- array()
		lTemp <- array()
		len <- dim(sasOut[[j]])[1]
		lenArray[[j]] = len
		while (k < len) {
			uTemp[l] = exp((sqrt(namesArray$V1[j])- fixedArray[j,1] - sasOut[[j]][k,3])/(fixedArray[j,2] + sasOut[[j]][k+1,3]))
			lTemp[l] = exp((sqrt(ODL)- fixedArray[j,1] - sasOut[[j]][k,3])/(fixedArray[j,2] + sasOut[[j]][k+1,3]))
			k <- k + 2
			l <- l + 1
		}

		lower[[j]] <- lTemp;
		upper[[j]] <- uTemp
		mean[[j]] <- mean(uTemp - lTemp)
		SD[[j]] <- sd(uTemp - lTemp)
		CoeffVar[[j]] <- SD[[j]]/mean[[j]]
		tempDiff[[j]] = uTemp - lTemp
		tempDiffY365 = tempDiff[[j]]
		bothN[[j]] = paste( (length(tempDiffY365[tempDiffY365 < 365])), length(tempDiffY365))
		count365[[j]] = length(tempDiffY365[tempDiffY365 < 365])
		CI[[j]] <- bootMean(tempDiffY365[tempDiffY365 < 365]) #modified to create condition that w[i] is less than 365 days 28.07.2011)
		CIout[j] <- paste(namesArray$V1[j], CI[[j]][3], (CI[[j]][3] - CI[[j]][1]), ((CI[[j]][2] - CI[[j]][3])), sep=" ")
	}
	
		#write mean and CI to csv file
		write.table(CIout, file = paste(fileReadLocation, fileRandomArray,".mean.csv", sep=""), sep = " ", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
	
	
	# <-------------------------------------------------End code ----------------------------------------------------------------------------------------->

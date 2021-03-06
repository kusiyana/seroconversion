	fileLocation <- "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/harg-22-10/"
	fileBase <- "SURV_04.08_var_OD_const_S.2_TL.0_TU.120_EASF.ENSF_OD"
	
	#declare functions
	
	makeSurvival <- function(dataNNu.Time.OD, cut, ExcludeSetting = 4) {
		# settings                               
		# cut = OD cut off
		# dataNNu.Time.OD = data input matrix (person Id, Time of test in days) 
		# setting: 										1 = remove already seroconverted at t0, 
		# 												2 = remove slow progressors
		#												3 = remove normal progressors
		#												4 = remove already seroconverted and slow progressors
		
		BEDH <- dataNNu.Time.OD
		d<-as.data.frame(BEDH)
		n <- dim(BEDH)[1] - 1
		nnu = max(d$ID)
		Ubound <- array(dim=nnu)
		Lbound <- array(dim=nnu)
		dataOut <- array(dim=c(nnu,4))
		
		d2 <- d[d$OD > cut,]
		count <- 1
		 
		LboundTemp <- NA
		UboundTemp <- NA
			
		 for (j in c(1:nnu)) {
		 	test1 <- FALSE
		 	test2 <- FALSE
			# were already seroconverted during all measurements
			if (dim(d[d$ID == j & d$OD > cut,])[1] == dim(d[d$ID == j,])[1]) {
				dataOut[j,1] <- j
				dataOut[j,2] <- 0
				dataOut[j,3] <- min(d[d$ID == j & d$OD > cut ,]$t)
				dataOut[j,4] <- 1
				test1 <- TRUE
			}
			
			#never seroconvert in any measurements
			if ((dim(d[d$ID == j & d$OD > cut,])[1] == 0) & test1 == FALSE) {
				dataOut[j,1] <- j
				dataOut[j,2] <-  max(d[d$ID == j & d$OD < cut ,]$t) -  min(d[d$ID == j & d$OD < cut ,]$t)
				dataOut[j,3] <- 365
				dataOut[j,4] <- 2
				test2 <- TRUE
			}
			
			#normal seroconversion during one of the measurements
			if ((dim(d[d$ID == j & d$OD > cut,])[1] > 0) & test2 == FALSE & test1 == FALSE) {
				UboundTemp  <- min(d[d$ID == j & d$OD > cut ,][2])
				tcut = min(d[d$ID == j & d$OD > cut ,][2])
				if (dim((d[d$ID == j & d$OD < cut ,][2]))[1] >= 2) {
					tmax <- dim((d[d$ID == j & d$OD < cut & d$t < tcut ,][2]))[1]
					tarray <- (d[d$ID == j & d$OD < cut & d$t < tcut,])$t
					LboundTemp <- tarray[tmax] - min(tarray)
				}
				if (dim((d[d$ID == j & d$OD < cut ,][2]))[1] == 1) {
					LboundTemp <- 0
				}
			dataOut[j,1] <- j
			dataOut[j,2] <- LboundTemp
			dataOut[j,3] <- UboundTemp
			dataOut[j,4] <-3
			}
		}
		d2 <- as.data.frame(dataOut)
		if (ExcludeSetting == 1) d3 <- d2[d2$V4 != 1,][2:3]       # remove already seroconverted at t0,
		if (ExcludeSetting == 2) d3 <- d2[d2$V4 != 2,][2:3]       # remove slow progressors
		if (ExcludeSetting == 3) d3 <- d2[d2$V4 != 3,][2:3]       # remove normal progressors
		if (ExcludeSetting == 4) d3 <- d2[d2$V4 != 1 & d2$V4 != 2,][2:3]    # remove already seroconverted and slow progressors
		if (ExcludeSetting != 1 & ExcludeSetting != 2 & ExcludeSetting != 3 & ExcludeSetting != 4) d3 <- d2[2:3]
		
		colnames(d3)[1] <- "L"
		colnames(d3)[2] <- "U"
		return(d3)
	} #end of function
	
	turnbull <- function (data,iter) {
		dfram <- data[,1:2]
		Times <- sort(unique(c(dfram$L, dfram$U)))
		nobs <- nrow(dfram)
		ntim <- length(Times)
		alpha <- t(outer(Times, dfram$L, function(x,y)  x > y) *
		           outer(Times, dfram$U, function(x,y)  x <= y)) 
		func <- function(a){
		for(i in 1:length(a)){if(a[i]==0)a[i]=0.0000001}
		                    a} 
		
		Hfunc <- function(Svec, alpha) {
		    pvec <- -diff(c(1,Svec))
		    qvec <- 1/c(alpha %*% pvec)  
		    dvec <- pvec * c(qvec %*% alpha)
		    dvec <- func(dvec)
		    Yvec <- rev(cumsum(rev(dvec)))
		    Yvec <- func(Yvec)
		    out <- cumprod(1-dvec/Yvec)
		    list(dvec=dvec,Yvec=Yvec,out=out) 
		} 
		Sinit <- ((ntim-1):0)/ntim
		Sold <- Sinit
		change <- numeric(iter) 
		Dmat <- matrix(0,ncol=iter,nrow=ntim)	
		Ymat <- matrix(0,ncol=iter,nrow=ntim)
	
		for(i in 1:iter) {
			Snew <- Hfunc(Sold,alpha)$out
			Dmat[,i] <-  Hfunc(Sold,alpha)$dvec
			Ymat[,i] <-  Hfunc(Sold,alpha)$Yvec
			change[i] <- round(sum(abs(Snew-Sold)),7)
			if(change[i]>0.01) {
				cat(i,"th iteration change = ",change[i],"\n")
				Sold <- Snew
			}
		}
		cat("...","\n")
		cat(iter,"th iteration change = ",change[iter],"\n\n") 
		Tabel <- cbind(Times=Times, Surv=round(Sold,4),dvec=Dmat[,iter],yvec=Ymat[,iter])
		list(Tabel=Tabel)
	}
	
	#survival mean function
	makeAvg <- function(survivalData, numPeople) { #NB check to see whether Cari has got numPeople right!
		d<-as.data.frame(survivalData)
		n <- dim(survivalData)[1] #- 1
		avg <- as.data.frame(array(dim=c(n,7)))
	 	colnames(avg)[1] <- "Times"
		colnames(avg)[2] <- "Surv"
		colnames(avg)[3] <- "mean"
		colnames(avg)[4] <- "S"
	 	colnames(avg)[5] <- "Dvec"
	 	colnames(avg)[6] <- "Yvec"
	 	colnames(avg)[7] <- "Var"
		
		 for (j in c(1:n)) {
			k <- j + 1
			avg[j,1] <- d$Times[j]
			avg[j,2] <- d$Surv[j]
		 	avg[j,3] <- (d$Times[k] - d$Times[j]) * d$Surv[j]
		 }
		 
		# calculate S
		avg$mean[is.na(avg$mean)] <- 0
		avg[1,4] <- sum(avg$mean[!is.na(avg$mean)])
	 	 for (k in c(2:n)) {
	 	 	j = k - 1
			avg[k,4] = avg[j,4] - avg[j,3]
		 }
		 
		# calculate Yvec, Dvec + variance
		 for (j in c(n:1)) {
			k <- j + 1
			avg[j,6] <- numPeople * d$Surv[j]
			avg[j,5] <- avg[j,6] - avg[k,6]
			avg[j,7] <- (avg[j,5] + (avg[j,4]^2))/(avg[j,6]*(avg[j,6] - avg[j,5]))
		 }
	
		 return(as.data.frame(avg))
	}
	
	#start code
	
	names <- read.table(file = paste(fileReadLocation, fileBase, ".names", sep=""), head=FALSE, sep=",")
	loopLimit <- dim(names)[1]
	
	#create different OD cutof scenarios and store in list
	Sdata <- list(); surv <- list(); avg <- list(); mean <- list(); var <- list(); CI <- list(); SD <- list();CoeffVar <- list()
	
	colours = c("red","blue","green","yellow","black","purple", "cyan","magenta", "orange", "pink", "brown", "violet")
	limitArray = array()
	meanOut = array()

	for (j in c(1:loopLimit)) {
		BEDH <- read.csv(file= paste(fileLocation, fileBase, ".", names$V1[j], ".csv", sep=""), head=TRUE, sep=" ")
	 	Sdata[[j]] <- makeSurvival(BEDH, names$V1[j], 0)
	 	for (k in c(1:dim(Sdata[[1]])[1])) {
	 		if(Sdata[[j]][k,2] > 365) Sdata[[j]][k,2] <- 365
	 	}
	 	surv[[j]] <- turnbull(Sdata[[j]],200)$Tabel[,1:2]
	 	avg[[j]] <- makeAvg(surv[[j]], max(BEDH$ID))
	 	mean[[j]] <- sum(avg[[j]]$mean[!is.na(avg[[j]]$mean)])
	 	var[[j]] <- sum(avg[[j]]$Var[!is.na(avg[[j]]$Var) & avg[[j]]$Var < Inf])
		SD[[j]] <- sqrt(var[[j]])
		CoeffVar[[j]] <- SD[[j]]/mean[[j]]
	 	CI[[j]] <- c(mean[[j]] - 1.96 * sqrt(var[[j]]), mean[[j]] + 1.96 * sqrt(var[[j]])) 
	 	if (j == 1) plot(surv[[j]][,1],surv[[j]][,2],xlim=c(0,500),ylim=c(0,1),xlab="window period (days)", ylab="estimated exceeding probability",main="Estimation of HIV window period since SC using Turnbull's Algorithm",type="S",col="red")
		else lines(surv[[j]][,1],surv[[j]][,2],col=colours[j],type="S")
		limitArray[j] = paste("OD.",names$V1[j], sep="")
		meanOut[j] <- paste(names$V1[j], mean[[j]], (mean[[j]] - CI[[j]][1]), (CI[[j]][2] - mean[[j]]), SD[[j]], CoeffVar[[j]], sep=" ")
	}
	
legend(x=300,y=0.8,legend=limitArray,col=colours,lty=1)
	
#write mean and CI to csv file
headings = paste("OD", "Mean_R", "-", "+", "SD", "coV", sep = " ")
write.table(headings, file = paste(fileLocation, fileBase,".mean.csv", sep=""), sep = " ", col.names = FALSE, append="FALSE", row.names=FALSE, quote=FALSE)
write.table(meanOut, file = paste(fileLocation, fileBase,".mean.csv", sep=""), sep = " ", col.names = FALSE, append="TRUE", row.names=FALSE, quote=FALSE)
	
		
		
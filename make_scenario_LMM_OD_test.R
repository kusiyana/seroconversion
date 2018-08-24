BEDH <- read.csv(file=fileName, head=TRUE, sep=" ")
ODcut = 0.8
Wcut
d<-as.data.frame(BEDH)
n <- dim(BEDH)[1] - 1
nnu = d$ID[n]
omitArray = array(dim=c(nnu,1))
dataOut <- array(dim=c(nnu,4))
countWmaxMat = array();
countWmax = 0
d2 <- array(dim=c(nnu,4))
d3 = 
count <- 0
timeDiff <- array();


for (j in c(1:nnu)) {
	if (d[d$ID == j,2][1] <= timeCutLower | d[d$ID == j,2][1] >= timeCutUpper) {		
		if ((dim(d[d$ID == j & d$OD < ODcut,])[1] > 0) & (dim(d[d$ID == j & d$OD > ODcut,])[1] > 0) ) { # test for normal conversion above thresh
			d2[d2$ID == j,][2] <- d2[d2$ID == j,][2]
			d2[d2$ID == j,][3] <- d2[d2$ID == j,][2]
			
			# now test for people who hit C AFTER Wcut
			
			tLower = d[d$ID == j,2][1]
			tUpper = d[d$ID == j & d$OD > ODcut,2][1]
			if ((tUpper - tLower) > Wcut) {
				d[d$ID == j & d$OD > 0.8,][2]
			}
			}
			

	}
}
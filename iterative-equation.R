#script to explore John's iterative formula

fileName = "C:/Documents and Settings/heastwood/My Documents/SACEMA/OD cut-off/hargrove-OD-22.10.2010.csv"
BEDH <- read.csv(file=fileName, head=TRUE, sep=" ")
d<-as.data.frame(BEDH)
n <- dim(BEDH)[1] - 1
nnu = d$ID[n]

#number of people who test recent at first test


#get number of people who still test recent after one year
NR1YEAR <-  d[d$time > 365 & d$OD < 0.8,] #n = 23

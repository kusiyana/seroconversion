bootMean <- function (data=window.elk[,3],B=1000,alpha=c(0.05,0.95),nbootsd=1000, nboott=1000, v.nbootg=1000,v.nbootsd=1000,v.nboott=1000) {
	library(bootstrap)
	n <- length(data)
	theta <- mean(data)
	boot.theta <- numeric(B)
	for(i in 1:B){
		boot.data <- data[sample(1:n,replace=TRUE)]
		boot.theta[i] <- mean(boot.data)
	}
	#hist(boot.theta,main="Non-parametric bootstrap histogram",xlab=paste("B=",B),ylab="theta")
	#abline(v=theta,lty=2)
	
	#standaard
	se.boot <- sqrt(var(boot.theta))
	onder.std <- theta - qnorm(alpha[2])*se.boot
	bo.std <- theta - qnorm(alpha[1])*se.boot
	leng.std <- bo.std-onder.std
	vorm.std <- (bo.std-theta)/(theta-onder.std)
	
	#persentiel
	sort.boot <- sort(boot.theta)
	onder.pers <- sort.boot[alpha[1]*B]
	bo.pers <- sort.boot[(alpha[2])*B]
	leng.pers <- bo.pers-onder.pers
	vorm.pers <- (bo.pers-theta)/(theta-onder.pers)
	
	#ABC
	funk <- function(p,x){mean <- sum(p*x)/sum(p)}
	ABC <- abcnon(data,funk,alpha=alpha)$limits[,2]
	leng.ABC <- max(ABC)-min(ABC)
	vorm.ABC <- (max(ABC)-theta)/(theta-min(ABC))
	
	#BCa
	BCa <-bcanon(data,nboot=B,theta=function(x)mean(x),alpha=alpha)$confpoint[,2]
	leng.BCa <- max(BCa)-min(BCa)
	vorm.BCa <- (max(BCa)-theta)/(theta-min(BCa))
	
	#boot-t
	uit <-boott(data,theta=function(x)mean(x),perc=alpha,nbootsd=nbootsd, nboott=nboott)$confpoints
	leng.t <- max(uit)-min(uit)
	vorm.t <- (max(uit)-theta)/(theta-min(uit))
	
	#boot-t(VS=T)
	uit2 <- boott(data,VS=TRUE,theta=function(x)mean(x), perc=alpha, nbootsd=nbootsd, nboott=nboott, v.nbootg=v.nbootg, v.nbootsd=v.nbootsd, v.nboott=v.nboott)$confpoints
	leng.t2 <- max(uit2)-min(uit2)
	vorm.t2 <- (max(uit2)-theta)/(theta-min(uit2))
	
	#basic
	onder.basic <- 2*theta - bo.pers
	bo.basic <- 2*theta - onder.pers
	leng.b <- bo.basic-onder.basic
	vorm.b <- (bo.basic-theta)/(theta-onder.basic)
	Tabel <- data.frame(onder=c(onder.std,onder.pers,ABC[1],BCa[1],uit[1],uit2[1], onder.basic),bo=c(bo.std,bo.pers,ABC[2],BCa[2],uit[2],uit2[2],bo.basic),
	lengte=c(leng.std,leng.pers,leng.ABC,leng.BCa,leng.t,leng.t2,leng.b),
	vorm=c(vorm.std,vorm.pers,vorm.ABC,vorm.BCa,vorm.t,vorm.t2,vorm.b))
	rownames(Tabel) <- c("Standaard","Persentiel","ABC","BCa","Boott","Boot-t(VS=T)","Basic Boot")
	list(boot.mean=mean(boot.theta),Tabel14.2=Tabel)
	return (c(onder.basic, bo.basic))
}
library(xtable)
load('AeDF.Rdata')
AeDF$AeX <- apply(matrix(c(AeDF[[2]], AeDF[[4]], AeDF[[6]], AeDF[[8]]), ncol=4), 1, median)
AeDF$AeY <- apply(matrix(c(AeDF[[3]], AeDF[[5]], AeDF[[7]], AeDF[[9]]), ncol=4), 1, median)
AeDF$Ae <- 0.5*(AeDF$AeX + AeDF$AeY)
antList <- unique(AeDF$Ant)
antList <- sort(as.character(antList[ grep('[C|P]', antList) ]))
antNum <- length(antList)
AeX_mean <- AeY_mean <- AeX_sd <- AeY_sd <- AeD_mean <- AeD_sd <- AeN_mean <- AeN_sd <-numeric(0)
for(ant_index in 1:antNum){
	DF <- AeDF[ (AeDF$Ant == antList[ant_index] ),]
	AeX_median <- median(DF$AeX); AeY_median <- median(DF$AeY)
	flag <- which( (abs(DF$AeX - AeX_median) < 15.0) |(abs(DF$AeY - AeY_median) < 15.0 ) )
    DF <- DF[flag,]
	AeX_mean[ant_index] <- quantile(DF$AeX, 0.75)
	AeY_mean[ant_index] <- quantile(DF$AeY, 0.75)
	AeX_sd[ant_index] <- sd(DF$AeX)
	AeY_sd[ant_index] <- sd(DF$AeY)
	#
	DayDF <- DF[((DF$sunset < 2.5) | (DF$sunset > 18)),]
	NgtDF <- DF[((DF$sunset > 2.5) & (DF$sunset < 18)),]
	AeD_mean[ant_index] <- quantile(DayDF$Ae, 0.75)
	AeN_mean[ant_index] <- quantile(NgtDF$Ae, 0.75)
	AeD_sd[ant_index] <- sd(DayDF$Ae)
	AeN_sd[ant_index] <- sd(NgtDF$Ae) 
}
AeStatDF <- data.frame(Ant=antList, AeXm = AeX_mean, AeXs = AeX_sd, AeYm = AeY_mean, AeYs = AeY_sd, AeDm = AeD_mean, AeDs = AeD_sd, AeNm = AeN_mean, AeNs = AeN_sd)
print(xtable(AeStatDF), include.rownames=F)

pdf('AeStat.pdf')
DayNight <- matrix(c(AeD_mean, AeN_mean - AeD_mean), ncol=2); rownames(DayNight) <- antList; colnames(DayNight) <- c('Night', 'Day')
DayNight <- t(DayNight)
barplot(DayNight, las=2, ylim=c(0,70), ylab='Efficiency [%]', main='Band 3 Day/Night Difference')
xaxis <- (0.2 + 1) * 1:16 - 1/2
arrows( xaxis-0.1, AeD_mean + AeD_sd, xaxis-0.1, AeD_mean - AeD_sd, length=0)
arrows( xaxis+0.1, AeN_mean + AeN_sd, xaxis+0.1, AeN_mean - AeN_sd, length=0)
text(5,65,sprintf('Night - Day = %.1f +- %.1f%%', mean(DayNight['Day',], na.rm=T), sd(DayNight['Day',], na.rm=T)))
dev.off()

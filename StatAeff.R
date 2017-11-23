library(xtable)
load('AeDF.Rdata')
AeDF$AeX <- apply(matrix(c(AeDF[[2]], AeDF[[4]], AeDF[[6]], AeDF[[8]]), ncol=4), 1, median)
AeDF$AeY <- apply(matrix(c(AeDF[[3]], AeDF[[5]], AeDF[[7]], AeDF[[9]]), ncol=4), 1, median)
AeDF$Ae <- 0.5*(AeDF$AeX + AeDF$AeY)
antList <- unique(AeDF$Ant)
antList <- sort(as.character(antList[ grep('[C|P]', antList) ]))
antNum <- length(antList)
AeX_mean <- AeY_mean <- AeX_e <- AeY_e <- AeX_sd <- AeY_sd <- AeD_mean <- AeD_e <- AeD_sd <- AeN_mean <- AeN_e <- AeN_sd <-numeric(0)
for(ant_index in 1:antNum){
	DF <- AeDF[ ((AeDF$Ant == antList[ant_index]) & (AeDF$EL >= 40.0)),]
	weightX <- weightY <- rep(1.0, length(DF$Ant))
	AeX_median <- median(DF$AeX); AeY_median <- median(DF$AeY)
	flag_X <- which( abs(DF$AeX - AeX_median) > 15.0 ); flag_Y <- which( abs(DF$AeY - AeY_median) > 15.0 )
	weightX[flag_X] <- 0.00; weightY[flag_Y] <- 0.00
	EL60 <- DF$EL - 60.0
	Xfit <- lm(DF$AeX ~ EL60, weights=weightX); Yfit <- lm(DF$AeY ~ EL60, weights=weightY)
	AeX_mean[ant_index] <- summary(Xfit)[['coefficients']][1,1]
	AeX_e[ant_index] <- summary(Xfit)[['coefficients']][1,2]
	AeX_sd[ant_index] <- sd(DF$AeX[which(weightX > 0.0)])
	AeY_mean[ant_index] <- summary(Yfit)[['coefficients']][1,1]
	AeY_e[ant_index] <- summary(Yfit)[['coefficients']][1,2]
	AeY_sd[ant_index] <- sd(DF$AeY[which(weightY > 0.0)])
	#
	DayDF <- DF[((DF$sunset < 2.5) | (DF$sunset > 18)),]
	NgtDF <- DF[((DF$sunset > 2.5) & (DF$sunset < 18)),]
	weightD <- rep(1.0, length(DayDF$Ant)); weightN <- rep(1.0, length(NgtDF$Ant))
	AeD_median <- median(DayDF$Ae); AeN_median <- median(NgtDF$Ae)
	flag_D <- which( abs(DayDF$Ae - AeD_median) > 15.0 ); flag_N <- which( abs(NgtDF$Ae - AeN_median) > 15.0 )
	weightD[flag_D] <- 0.00; weightN[flag_N] <- 0.00
	EL60 <- DayDF$EL - 60.0; Dfit <- lm(DayDF$Ae ~ EL60, weights=weightD)
	EL60 <- NgtDF$EL - 60.0; Nfit <- lm(NgtDF$Ae ~ EL60, weights=weightN);
	AeD_mean[ant_index] <- summary(Dfit)[['coefficients']][1,1]
	AeD_e[ant_index] <- summary(Dfit)[['coefficients']][1,2]
	AeD_sd[ant_index] <- sd(DayDF$Ae[which(weightD > 0.0)])
	AeN_mean[ant_index] <- summary(Nfit)[['coefficients']][1,1]
	AeN_e[ant_index] <- summary(Nfit)[['coefficients']][1,2]
	AeN_sd[ant_index] <- sd(NgtDF$Ae[which(weightN > 0.0)])
}
AeStatDF <- data.frame(Ant=antList, AeXm = AeX_mean, AeXe = AeX_e, AeXs = AeX_sd, AeYm = AeY_mean, AeYe = AeY_e, AeYs = AeY_sd, AeDm = AeD_mean, AeDe = AeD_e, AeDs = AeD_sd, AeNm = AeN_mean, AeNe = AeN_e, AeNs = AeN_sd)
print(xtable(AeStatDF), include.rownames=F)

pdf('AeStat.pdf')
DayNight <- matrix(c(AeD_mean, AeN_mean - AeD_mean), ncol=2); rownames(DayNight) <- antList; colnames(DayNight) <- c('Night', 'Day')
DayNight <- t(DayNight)
barplot(DayNight, las=2, ylim=c(0,70), ylab='Efficiency [%]', main='Band 6 Day/Night Difference')
xaxis <- (0.2 + 1) * 1:16 - 1/2
arrows( xaxis-0.1, AeD_mean + AeD_e, xaxis-0.1, AeD_mean - AeD_e, length=0)
arrows( xaxis+0.1, AeN_mean + AeN_e, xaxis+0.1, AeN_mean - AeN_e, length=0)
text(5,65,sprintf('Night - Day = %.1f +- %.1f%%', mean(DayNight['Day',], na.rm=T), sd(DayNight['Day',], na.rm=T)))
dev.off()
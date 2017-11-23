library(xtable)
library(maptools)
ALMA_POS <- matrix(c( -67.755, -23.029 ), nrow=1 )
#-------- Parse arguments
parseArg <- function( args ){
	argNum <- length(args)
	fileNum <- argNum
	return( list(filelist = args[1:argNum]))
}

#-------- Find Calibrator name
findCalibrator <- function( Lines ){
	scalerPointer    <- grep("Scaler", Lines)
	equalizerPointer <- grep("Equalizer", Lines)
	datePointer      <- grep("Flux Calibrator is", Lines)
	scalerName    <- strsplit(Lines[scalerPointer], ' ')[[1]][2]
	equalizerName <- strsplit(Lines[equalizerPointer], ' ')[[1]][2]
	scaleEL <- as.numeric( strsplit(Lines[scalerPointer], ' ')[[1]][5] )
	scalerUTC <- strptime(strsplit(Lines[datePointer], ' ')[[1]][6], "%Y/%m/%d/%H:%M:%S", tz="UTC")
	sunsetUTC <- sunriset(ALMA_POS, as.POSIXct(scalerUTC), POSIXct.out=T, direction="sunset")[[2]]
	return(list(scaler=scalerName, EL=scaleEL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24, equalizer=equalizerName))
}
#-------- Read Aeff Section
readAeffSection <- function(Lines){
	pointer <- grep(" Aeff", Lines) + 1
	antName <- character(0)
	Ae <- list(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))
	while( is.na(strsplit(Lines[pointer], ' ')[[1]][1]) == F){
		antName <- append(antName, strsplit(Lines[pointer], ' ')[[1]][1])
		for(spw in 1:8){
			Ae[[spw]] <- append(Ae[[spw]], as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][spw+2]))
		}
		pointer <- pointer + 1
	}
	names(Ae) <- c('AeX1', 'AeY1', 'AeX2', 'AeY2', 'AeX3', 'AeY3', 'AeX4', 'AeY4')
	return(data.frame(Ant=antName, AeX1=Ae[[1]], AeY1=Ae[[2]], AeX2=Ae[[3]], AeY2=Ae[[4]], AeX3=Ae[[5]], AeY3=Ae[[6]], AeX4=Ae[[7]], AeY4=Ae[[8]]))
}

#-------- Start program
Arguments <- commandArgs(trailingOnly = T)
fileList <- Arguments
#fileList <- c("uid___A002_Xabc3f8_X196f-RB_06-Flux.log")
FMT <- c('Ant', 'AeX1', 'AeY1', 'AeX2', 'AeY2', 'AeX3', 'AeY3', 'AeX4', 'AeY4', 'EL', 'Date', 'sunset', 'EQ')
AeDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(AeDF) <- FMT
for(fileName in fileList){
	cat(fileName); cat('\n')
    fileLines <- readLines(fileName)
	CalList <- findCalibrator(fileLines)
	DF <- readAeffSection(fileLines)
	DF$FL <- CalList$scaler
	DF$EL <- CalList$EL
	DF$Date <- CalList$UTC
	DF$sunset <- CalList$sunset
	DF$EQ <- CalList$equalizer
	DF$File <- fileName
	AeDF <- rbind(AeDF, DF)
}
save(AeDF, file='AeDF.Rdata')
logfileName <- 'AeDF.log'
write(fileList, file=logfileName, append=F)
#-------- Statistics
AeDF$AeX <- apply(matrix(c(AeDF[[2]], AeDF[[4]], AeDF[[6]], AeDF[[8]]), ncol=4), 1, median)
AeDF$AeY <- apply(matrix(c(AeDF[[3]], AeDF[[5]], AeDF[[7]], AeDF[[9]]), ncol=4), 1, median)
antList <- unique(AeDF$Ant); antNum <- length(antList)  # Unique antenna list
FLSList <- unique(AeDF$FL);  FLSNum <- length(FLSList)  # Unique flux-calibrator list
for(ant_index in 1:antNum){
    pdf(sprintf("Ae%s.pdf", antList[ant_index]))
    #-------- EL plot
    plot(AeDF$EL, AeDF$AeX, type='n', ylim=c(20.0, 90.0), xlab='Elevation [deg]', ylab='Efficiency [%]', main=antList[ant_index])
    for(FLS_index in 1:FLSNum){
        Hindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL >= 40.0))
        Lindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL <  40.0))
        points(AeDF$EL[Hindex], AeDF$AeX[Hindex], pch=15, col=FLS_index); points(AeDF$EL[Hindex], AeDF$AeY[Hindex], pch= 2, col=FLS_index)
        points(AeDF$EL[Lindex], AeDF$AeX[Lindex], pch=15, col=FLS_index, cex=0.2); points(AeDF$EL[Lindex], AeDF$AeY[Lindex], pch= 2, col=FLS_index, cex=0.2)
    }
    index <- which((AeDF$Ant == antList[ant_index]) & (AeDF$EL >= 40.0))
	weightX <- weightY <- rep(1.0, length(index))
	AeX_median <- median(AeDF$AeX[index]); AeY_median <- median(AeDF$AeY[index])
	if(length(index) > 4 ){
        flag_X <- which( abs(AeDF$AeX[index] - AeX_median) > 15.0 ); flag_Y <- which( abs(AeDF$AeY[index] - AeY_median) > 15.0 )
        weightX[flag_X] <- 0.00; weightY[flag_Y] <- 0.00
        EL60 <- AeDF$EL[index] - 60.0
        Xfit <- lm(AeDF$AeX[index] ~ EL60, weights=weightX); Yfit <- lm(AeDF$AeY[index] ~ EL60, weights=weightY)
        AeX_mean <- summary(Xfit)[['coefficients']][1,1]; AeX_sd <- summary(Xfit)[['coefficients']][1,2]
        AeY_mean <- summary(Yfit)[['coefficients']][1,1]; AeY_sd <- summary(Yfit)[['coefficients']][1,2]
        legend("topleft", legend=FLSList, col=seq(1,FLSNum), pch=15, cex=0.7)
        legend("topright", legend=c(sprintf("PolX %4.1f +- %.1f %%", AeX_mean, AeX_sd), sprintf("PolY %4.1f +- %.1f %%", AeY_mean, AeY_sd)), col=rep(1, 2), pch=c(15,2))
        write(sprintf("%s: PolX %4.1f +- %.1f %.1f", antList[ant_index], AeX_mean, AeX_sd, sd(AeDF$AeX[index])), file=logfileName, append=T)
        write(sprintf("%s: PolY %4.1f +- %.1f %.1f", antList[ant_index], AeY_mean, AeY_sd, sd(AeDF$AeY[index])), file=logfileName, append=T)
		cat(sprintf('%s  %6.2f  %6.2f\n', antList[ant_index], AeX_mean, AeY_mean))
    } else {
		cat(sprintf('%s  %6.2f  %6.2f\n', antList[ant_index], AeX_median, AeY_median))
	}
    #-------- Date plot
    plot(AeDF$Date, AeDF$AeX, xaxt='n', type='n', ylim=c(20.0, 90.0), xlab='Epoch', ylab='Efficiency [%]', main=antList[ant_index])
    axis.POSIXct(1, at=AeDF$Date[index], format="%Y-%m-%d")
    legTex <- FLSList
    for(FLS_index in 1:FLSNum){
        Hindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL >= 40.0))
        Lindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL <  40.0))
        points(AeDF$Date[Hindex], AeDF$AeX[Hindex], pch=15, col=FLS_index); points(AeDF$Date[Hindex], AeDF$AeY[Hindex], pch= 2, col=FLS_index)
        points(AeDF$Date[Lindex], AeDF$AeX[Lindex], pch=15, col=FLS_index, cex=0.2); points(AeDF$Date[Lindex], AeDF$AeY[Lindex], pch= 2, col=FLS_index, cex=0.2)
        AeS <- median(c(AeDF$AeX[Hindex], AeDF$AeS[Hindex] ))
        legTex[FLS_index] <- sprintf("%.1f%% %s", AeS, FLSList[FLS_index])
    }
    legend("topleft", legend=legTex, col=seq(1,FLSNum), pch=15, cex=0.7)
    #-------- Sunset plot
    plot(AeDF$sunset, AeDF$AeX, type='n', ylim=c(20.0, 90.0), xlab='Hours after sunset', ylab='Efficiency [%]', main=antList[ant_index])
    polygon(c(0.0, 0.0, 2.0, 2.0), c(0.0, 100.0, 100.0, 0.0), col='cyan', border=NA) 
	polygon(c(18.0, 18.0, 24.0, 24.0), c(0.0, 100.0, 100.0, 0.0), col='cyan', border=NA)
    for(FLS_index in 1:FLSNum){
        Hindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL >= 40.0))
        Lindex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$FL == FLSList[FLS_index]) & (AeDF$EL <  40.0))
        points(AeDF$sunset[Hindex], AeDF$AeX[Hindex], pch=15, col=FLS_index); points(AeDF$sunset[Hindex], AeDF$AeY[Hindex], pch= 2, col=FLS_index)
        points(AeDF$sunset[Lindex], AeDF$AeX[Lindex], pch=15, col=FLS_index, cex=0.2); points(AeDF$sunset[Lindex], AeDF$AeY[Lindex], pch= 2, col=FLS_index, cex=0.2)
    }
    DayIndex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$EL >= 40.0) & ((AeDF$sunset < 2.0) | (AeDF$sunset > 18.0)))
    NgtIndex <- which((AeDF$Ant == antList[ant_index]) & (AeDF$EL >= 40.0) & ((AeDF$sunset > 2.0) & (AeDF$sunset < 18.0)))
    if( (length(DayIndex) > 4) & length(NgtIndex > 4)){
        AeDay <- 0.5*(AeDF$AeX[DayIndex] + AeDF$AeY[DayIndex]); flagDay <- which( abs(AeDay - median(AeDay)) > 15.0 ) 
        AeNgt <- 0.5*(AeDF$AeX[NgtIndex] + AeDF$AeY[NgtIndex]); flagNgt <- which( abs(AeNgt - median(AeNgt)) > 15.0 ) 
        EL60D <- AeDF$EL[DayIndex] - 60.0; EL60N <- AeDF$EL[NgtIndex] - 60.0
        WeightDay <- rep(1.0, length(DayIndex)); WeightDay[flagDay] <- 0.0
        WeightNgt <- rep(1.0, length(NgtIndex)); WeightNgt[flagNgt] <- 0.0
        DayFit <- lm(AeDay ~ EL60D, weights=WeightDay); NgtFit <- lm(AeNgt ~ EL60N, weights=WeightNgt) 
        Day_mean <- summary(DayFit)[['coefficients']][1,1];; Day_sd <- summary(DayFit)[['coefficients']][1,2]
        Ngt_mean <- summary(NgtFit)[['coefficients']][1,1];; Ngt_sd <- summary(NgtFit)[['coefficients']][1,2]
        legend("top", legend=c(sprintf("Day %4.1f +- %.1f %%", Day_mean, Day_sd), sprintf("Night %4.1f +- %.1f %%", Ngt_mean, Ngt_sd)), col=c('cyan', 'black'), pch=c(15,15))
        write(sprintf("%s: Day %4.1f +- %.1f %.1f", antList[ant_index], Day_mean, Day_sd, sd(AeDay)), file=logfileName, append=T)
        write(sprintf("%s: Ngt %4.1f +- %.1f %.1f", antList[ant_index], Ngt_mean, Ngt_sd, sd(AeNgt)), file=logfileName, append=T)
    }
    dev.off()
}
#-------- Xtable
FMT <- c('Date', 'UID', 'Antennas')
DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
fileList <- unique(AeDF$File); fileNum <- length(fileList)
for( file_index in 1:fileNum){
    index <- which(AeDF$File == fileList[file_index])
    tempDF <- data.frame( Date = as.character(as.Date(AeDF$Date[index[1]]), format="%Y-%m-%d"), UID = sprintf('%s/%s',strsplit(fileList[file_index], '_+|-')[[1]][3], strsplit(fileList[file_index], '_+|-')[[1]][4]), Antennas=paste(AeDF$Ant[index], collapse=', '))
    DF <- rbind(DF, tempDF)
}
print(xtable(DF), include.rownames=F)

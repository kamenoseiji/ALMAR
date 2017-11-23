#library(xtable)
library(maptools)
ALMA_POS <- matrix(c( -67.755, -23.029 ), nrow=1 )

sourceMatch <- function(sourceName){
	sourceDict <- list(
		c('J0237+2848', 'J0237+288'),
		c('J0238+1636', 'J0238+166'),
		c('J0319+4130', '3c84'),
		c('J0334-4008', 'J0334-401'),
		c('J0423-0120', 'J0423-013'),
		c('J0510+1800', 'J0510+180'),
		c('J0519-4546', 'J0519-454'),
		c('J0522-3627', 'J0522-364'),
		c('J1037-2934', 'J1037-295'),
		c('J1058+0133', 'J1058+015'),
		c('J1107-4449', 'J1107-448'),
		c('J1146+3958', 'J1146+399'),
		c('J1229+0203', '3c273'),
		c('J1256-0547', '3c279'),
		c('J1337-1257', 'J1337-129'),
		c('J1427-4206', 'J1427-421'),
		c('J1517-2422', 'J1517-243'),
		c('J1550+0527', 'J1550+054'),
		c('J1617-5848', 'J1613-586'),
		c('J1642+3948', '3c345'),
		c('J1733-1304', 'J1733-130'), 
		c('J1751+0939', 'J1751+096'), 
		c('J1924-2914', 'J1924-292'),
		c('J2025+3343', 'J2025+337'),
		c('J2056-4714', 'J2056-472'),
		c('J2148+0657', 'J2148+069'),
		c('J2232+1143', 'J2232+117'),
		c('J2253+1608', '3c454.3'),
		c('J2258-2758', 'J2258-279'))
	#
	for(index in 1:length(sourceDict)){ if( !is.na(match(sourceName, sourceDict[[index]]))){ return(sourceDict[[index]][1])} }
	return(sourceName)
}
#-------- Parse arguments
parseArg <- function( args ){
	argNum <- length(args)
	fileNum <- argNum
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-D"){ fileDir <- substring(args[index], 3);  fileNum <- fileNum - 1}
	}
	return( list(fileDir = fileDir, filelist = args[(argNum - fileNum + 1):argNum]))
}

#-------- match time
matchTime <- function(smallSet, largeSet){
	timeNum <- length(smallSet)
	pointer <- integer(timeNum)
	for(index in 1:timeNum){
		pointer[index] <- which.min( abs(largeSet - smallSet[index]) )
	}
	return(pointer)
}

#-------- Find Stokes Parameters
readStokesSection <- function(Lines){
	pointer <- grep(" -----------------------------------------------------------------------------------------", Lines)
	numSource <- length(pointer)
	StokesI <- StokesQ <- StokesU <- StokesV <- numeric(4)
	errI <- errQ <- errU <- errV <- numeric(4)
	I <- Q <- U <- V <- numeric(0)
	eI <- eQ <- eU <- eV <- numeric(0)
	srcList <- character(0); EL <- numeric(0); freq <- numeric(4)
	for(srcIndex in 1:numSource){
		if(length(grep('Only', Lines[(pointer[srcIndex]+1):(pointer[srcIndex]+7)])) > 0){ next }
		if(length(grep('nan', Lines[(pointer[srcIndex]+1):(pointer[srcIndex]+7)])) > 0){ next }
		srcList <- append(srcList, as.character(strsplit(Lines[pointer[srcIndex] -2], '[ |=]+')[[1]][3]) )
		EL <- append(EL, strsplit(Lines[pointer[srcIndex] -2], '[ |=]+')[[1]][3])
		for(spw_index in 1:4){
			lineElements <- strsplit(Lines[pointer[srcIndex] + spw_index], '[ |(|)]+')[[1]]
			freq[spw_index] <- as.numeric(lineElements[2])
			StokesI[spw_index] <- as.numeric(lineElements[4]); errI[spw_index] <- as.numeric(lineElements[5])
			StokesQ[spw_index] <- as.numeric(lineElements[6]); errQ[spw_index] <- as.numeric(lineElements[7])
			StokesU[spw_index] <- as.numeric(lineElements[8]); errU[spw_index] <- as.numeric(lineElements[9])
			StokesV[spw_index] <- as.numeric(lineElements[10]);errV[spw_index] <- as.numeric(lineElements[11])
		}
		FREQ <- freq - mean(freq)
		fitI <- lm(formula= StokesI ~ FREQ, weights=(errI + 0.0005)^(-2))
		fitQ <- lm(formula= StokesQ ~ FREQ, weights=(errQ + 0.0005)^(-2))
		fitU <- lm(formula= StokesU ~ FREQ, weights=(errU + 0.0005)^(-2))
		fitV <- lm(formula= StokesV ~ FREQ, weights=(errV + 0.0005)^(-2))
		I <- append(I, coef(fitI)[[1]]); eI <- append(eI, coef(summary(fitI))[,"Std. Error"][[1]])
		Q <- append(Q, coef(fitQ)[[1]]); eQ <- append(eQ, coef(summary(fitQ))[,"Std. Error"][[1]])
		U <- append(U, coef(fitU)[[1]]); eU <- append(eU, coef(summary(fitU))[,"Std. Error"][[1]])
		V <- append(V, coef(fitV)[[1]]); eV <- append(eV, coef(summary(fitV))[,"Std. Error"][[1]])
		#cat(sprintf('Scan%d %s %f %f %f %f\n', srcIndex, srcList[srcIndex], I[srcIndex], Q[srcIndex], U[srcIndex], V[srcIndex]))
	}
	return(data.frame(Src=as.character(srcList), EL=EL, I=I, Q=Q, U=U, V=V, eI=eI, eQ=eQ, eU=eU, ev=eV))
}

#-------- Find Calibrator name
findCalibrator <- function( Lines ){
	scalerPointer    <- grep("Scaler", Lines)
	equalizerPointer <- grep("Equalizer", Lines)
	datePointer      <- grep("Flux Calibrator is", Lines)
	if( length(scalerPointer) == 0 ){ scalerPointer <- equalizerPointer}
	scalerName    <- strsplit(Lines[scalerPointer], ' ')[[1]][2]
	scaleEL <- as.numeric( strsplit(Lines[scalerPointer], ' ')[[1]][5] )
	equalizerName <- strsplit(Lines[equalizerPointer], ' ')[[1]][2]
	if(length(datePointer) > 0){
		scalerUTC <- strptime(strsplit(Lines[datePointer], ' ')[[1]][6], "%Y/%m/%d/%H:%M:%S", tz="UTC")
	} else {
		scalerUTC <- strptime(strsplit(Lines[equalizerPointer], ' ')[[1]][7], "%Y/%m/%d/%H:%M:%S", tz="UTC")
	}
	sunsetUTC <- sunriset(ALMA_POS, as.POSIXct(scalerUTC), POSIXct.out=T, direction="sunset")[[2]]
	return(list(scaler=scalerName, EL=scaleEL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24, equalizer=equalizerName))
}
#-------- Read Aeff Section
#readAeffSection <- function(Lines){
#	pointer <- grep(" Aeff", Lines) + 1
#	antName <- character(0)
#	Ae <- list(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))
#	while( is.na(strsplit(Lines[pointer], ' ')[[1]][1]) == F){
#		antName <- append(antName, strsplit(Lines[pointer], ' ')[[1]][1])
#		for(spw in 1:8){
#			Ae[[spw]] <- append(Ae[[spw]], as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][spw+2]))
#		}
#		pointer <- pointer + 1
#	}
#	names(Ae) <- c('AeX1', 'AeY1', 'AeX2', 'AeY2', 'AeX3', 'AeY3', 'AeX4', 'AeY4')
#	return(data.frame(Ant=antName, AeX1=Ae[[1]], AeY1=Ae[[2]], AeX2=Ae[[3]], AeY2=Ae[[4]], AeX3=Ae[[5]], AeY3=Ae[[6]], AeX4=Ae[[7]], AeY4=Ae[[8]]))
#}

#-------- Read Trec Section
readTrecSection <- function(Lines){
	pointer <- grep(" Trec:", Lines) + 2
	antName <- character(0)
	Tr <- list(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))
	Ts <- list(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0))
	while(1){
		elements <- strsplit(Lines[pointer], '[ |K]+')[[1]]
		if(elements[1] == ''){ break }
		antName <- append(antName, elements[1])
		for(spw in 1:8){
			Tr[[spw]] <- append(Tr[[spw]], as.numeric(elements[spw+2]))
		}
		pointer <- pointer + 1
	}
	names(Tr) <- c('TrX1', 'TrY1', 'TrX2', 'TrY2', 'TrX3', 'TrY3', 'TrX4', 'TrY4')
	TrsDF <- data.frame(Ant=antName, TrX1=Tr[[1]], TrY1=Tr[[2]], TrX2=Tr[[3]], TrY2=Tr[[4]], TrX3=Tr[[5]], TrY3=Tr[[6]], TrX4=Tr[[7]], TrY4=Tr[[8]])
	#
	pointer <- grep(" Tsys:", Lines) + 1
	for(ant_index in 1:length(antName)){
		elements <- strsplit(Lines[pointer + ant_index], '[ |K]+')[[1]]
		for(spw in 1:8){
			Ts[[spw]] <- append(Ts[[spw]], as.numeric(elements[spw+2]))
		}
	}
	names(Ts) <- c('TsX1', 'TsY1', 'TsX2', 'TsY2', 'TsX3', 'TsY3', 'TsX4', 'TsY4')
	TrsDF$TsX1 <- Ts[[1]]; TrsDF$TsY1 <- Ts[[2]]
	TrsDF$TsX2 <- Ts[[3]]; TrsDF$TsY2 <- Ts[[4]]
	TrsDF$TsX3 <- Ts[[5]]; TrsDF$TsY3 <- Ts[[6]]
	TrsDF$TsX4 <- Ts[[7]]; TrsDF$TsY4 <- Ts[[8]]
	return(TrsDF)
}

#-------- CheckTrex
checkTrec <- function(DF){
	#-------- Check Trec
	for(col_index in 2:9){ negTrecList <- which(TrsDF[[col_index]] < 0.0) }
	for(col_index in 2:9){ negTsysList <- which(TrsDF[[col_index+8]] < 0.0) }
	for(col_index in 2:9){ netTskyList <- which(TrsDF[[col_index+8]] < TrsDF[[col_index]]) }
	return( list(negTrecList, negTsysList, netTskyList) )
}

#-------- Start program
argList <- parseArg(commandArgs(trailingOnly = T))
fileList <- argList$filelist
fileDir <- argList$fileDir

#fileList <- c("uid___A002_Xacc4e4_X1876-RB_03-Flux.log")
FMT <- c('Src', 'EL', 'I', 'Q', 'U', 'V', 'eI', 'eQ', 'eU', 'eV', 'EL')
FLDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(FLDF) <- FMT
APDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(APDF) <- FMT

flagNum <- list()
for(fileName in fileList){
	cat(fileName); cat('\n')
    FLfileLines <- readLines(fileName); APfileLines <- readLines(sprintf('%s%s', fileDir, fileName))
	#-------- Read from FL files
	CalList <- findCalibrator(FLfileLines)
	DF <- readStokesSection(FLfileLines)
	DF$Date <- CalList$UTC
	DF$File <- fileName
	TrsDF <- readTrecSection(FLfileLines)
	flagList <- checkTrec(TrsDF)
	DF$flagNum <- length(flagList[[1]]) + length(flagList[[2]]) + length(flagList[[3]])
	FLDF <- rbind(FLDF, DF)
	#-------- Read from AP files
	CalList <- findCalibrator(APfileLines)
	DF <- readStokesSection(APfileLines)
	DF$Date <- CalList$UTC
	DF$File <- fileName
	TrsDF <- readTrecSection(APfileLines)
	flagList <- checkTrec(TrsDF)
	DF$flagNum <- length(flagList[[1]]) + length(flagList[[2]]) + length(flagList[[3]])
	APDF <- rbind(APDF, DF)
}
FLDF$Src <- as.character(lapply(as.character(FLDF$Src), sourceMatch))
FLDF$P <- sqrt(FLDF$Q^2 + FLDF$U^2)
FLDF$eP <- sqrt(FLDF$eQ^2 + FLDF$eU^2)
FLDF$EVPA <- 0.5* atan2(FLDF$U, FLDF$Q)
FLDF$eEVPA <- FLDF$eP/FLDF$P
APDF$Src <- as.character(lapply(as.character(APDF$Src), sourceMatch))
APDF$P <- sqrt(APDF$Q^2 + APDF$U^2)
APDF$eP <- sqrt(APDF$eQ^2 + APDF$eU^2)
APDF$EVPA <- 0.5* atan2(APDF$U, APDF$Q)
APDF$eEVPA <- APDF$eP/APDF$P
save(FLDF, file='Flux.Rdata')
save(APDF, file='APFlux.Rdata')
BandName <- sprintf('Band-%d', as.numeric(strsplit(FLDF$File[1], '[-|.|_]')[[1]][8]))
#-------- Source List
sourceList <- unique(FLDF$Src)
numSrc <- length(sourceList)
pdf('Flux.pdf', width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))

FluxRatio <- numeric(0)

for(src_index in 1:numSrc){
	if(substr(sourceList[src_index], 1, 1) != 'J'){ next }	# Skip SSO
	
	cat(sourceList[src_index]); cat('\n')

	DFfl <- FLDF[FLDF$Src == sourceList[src_index],]
	DFap <- APDF[APDF$Src == sourceList[src_index],]
	
	
	#-------- Plot 
	plot(DFfl$Date, DFfl$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='Using Flux Calibrator', ylim=c(0, 1.2*max(DFfl$I)))
	color_vector <- rep("#000000FF", length(DFfl$Date))
	color_vector[which(DFfl$flagNum > 0)] <- "#0000FF3F"
	arrows(as.numeric(DFfl$Date), DFfl$I - DFfl$eI, as.numeric(DFfl$Date), DFfl$I + DFfl$eI, angle=90, length=0.0, col=color_vector)
	points(DFfl$Date, DFfl$I, pch=20, col=color_vector)
	
	#-------- Plot polarized flux
	plot(DFap$Date, DFap$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='A Priori Calibration', ylim=c(0, 1.2*max(DFfl$I)))
	color_vector <- rep("#000000FF", length(DFap$Date))
	color_vector[which(DFap$flagNum > 0)] <- "#0000FF3F"
	arrows(as.numeric(DFap$Date), DFap$I - DFap$eI, as.numeric(DFap$Date), DFap$I + DFap$eI, angle=90, length=0.0, col=color_vector)
	points(DFap$Date, DFap$I, pch=20, col=color_vector)
	
	#-------- Plot Ratio
	plot(DFfl$Date, DFfl$EVPA, type='n', xlab='Date', ylab='Ratio - 1.0 (%)', main='A priori / Relative to calibrator', ylim=c(-50,50))
	color_vector <- rep("#000000FF", length(DFfl$Date))
	color_vector[which(DFfl$flagNum > 0)] <- "#0000FF3F"
	abline(h=0, lwd=0.1)
	if(length(DFfl$Date) > length(DFap$Date)){
		pointer <- matchTime( as.numeric(DFap$Date), as.numeric(DFfl$Date) )
		DFfl <- DFfl[pointer,]
	}
	if(length(DFfl$Date) < length(DFap$Date)){
		pointer <- matchTime( as.numeric(DFfl$Date), as.numeric(DFap$Date) )
		DFap <- DFap[pointer,]
	}
	Ratio <- 100.0* (DFap$I / DFfl$I - 1.0)
	error <- 100.0* sqrt((DFfl$eI / DFfl$I)^2 + (DFap$eI / DFap$I)^2)
	arrows(as.numeric(DFfl$Date), Ratio - error, as.numeric(DFfl$Date), Ratio + error, angle=90, length=0.0, col=color_vector)
	points(DFfl$Date, Ratio, pch=20, col=color_vector)
	
	FluxRatio <- c(FluxRatio, Ratio)
	
	mtext(side = 3, line=1, outer=T, text = sprintf('%s %s',sourceList[src_index], BandName), cex=2)
}
dev.off()
par(par.old)

pdf('HistRatio.pdf')
hist(FluxRatio, breaks=c(min(FluxRatio)-1, seq(-50,50,5), max(FluxRatio)+1), xlab='Relative Difference (%)', ylab='Fraction', main=sprintf('%s Flux Ratio (A priori / Relative to calibrator)', BandName), xlim=c(-50,50))
flag <- which( (FluxRatio > -50.0) & (FluxRatio < 50.0))
flag5 <- which( (FluxRatio > -5.0) & (FluxRatio < 5.0))
flag10 <- which( (FluxRatio > -10.0) & (FluxRatio < 10.0))
flag20 <- which( (FluxRatio > -20.0) & (FluxRatio < 20.0))
cat(sprintf('mean=%f  median=%f  sd=%f\n', mean(FluxRatio[flag]), median(FluxRatio[flag]), sd(FluxRatio[flag])))
cat(sprintf('accuracy5 : %4.1f, accracy10 : %4.1f, accuracy20 : %4.1f\n', 100*length(flag5)/length(flag), 100*length(flag10)/length(flag), 100*length(flag20)/length(flag) ))




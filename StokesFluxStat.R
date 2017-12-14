#library(xtable)
#library(maptools)
#ALMA_POS <- matrix(c( -67.755, -23.029 ), nrow=1 )

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
	return( list(filelist = args[1:argNum]))
}

#-------- Find Stokes Parameters
readStokesSection <- function(Lines){
	pointer <- grep("mean ", Lines)
	numSource <- length(pointer)
	StokesI <- StokesQ <- StokesU <- StokesV <- numeric(numSource)
	errI <- errQ <- errU <- errV <- numeric(numSource)
	I <- Q <- U <- V <- numeric(0)
	eI <- eQ <- eU <- eV <- numeric(0)
	srcList <- character(0); EL <- numeric(0)
	srcUTC <- as.Date(as.character(NULL))
	for(srcIndex in 1:numSource){
		srcList <- append(srcList, as.character(strsplit(Lines[pointer[srcIndex] -8], '[ |=]+')[[1]][3]) )
		srcUTC <- append(srcUTC, strptime(strsplit(Lines[pointer[srcIndex] -8], ' +')[[1]][6], "%Y/%m/%d/%H:%M:%S", tz="UTC"))
		EL <- append(EL, as.numeric(strsplit(Lines[pointer[srcIndex] -8], '[ |=]+')[[1]][5]))
		lineElements <- strsplit(Lines[pointer[srcIndex]], '[ |(|)|z]+')[[1]]
		FREQ <- as.numeric(lineElements[3])
		I <- append(I, as.numeric(lineElements[5])); eI <- append(eI, as.numeric(lineElements[6]))
		Q <- append(Q, as.numeric(lineElements[7])); eQ <- append(eQ, as.numeric(lineElements[8]))
		U <- append(U, as.numeric(lineElements[9])); eU <- append(eU, as.numeric(lineElements[10]))
		V <- append(V, as.numeric(lineElements[11])); eV <- append(eV, as.numeric(lineElements[12]))
	}
	return(data.frame(Src=as.character(srcList), Freq=FREQ, EL=EL, I=I, Q=Q, U=U, V=V, eI=eI, eQ=eQ, eU=eU, eV=eV, Date=srcUTC))
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
#	sunsetUTC <- sunriset(ALMA_POS, as.POSIXct(scalerUTC), POSIXct.out=T, direction="sunset")[[2]]
	return(list(scaler=scalerName, EL=scaleEL, UTC=scalerUTC, equalizer=equalizerName))
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

#-------- remove blank lines
removeBlank <- function(Lines){
	lineLength <- nchar(Lines)
	index <- which(lineLength > 1)
	return(Lines[index])
}


#-------- Start program
Arguments <- commandArgs(trailingOnly = T)
fileList <- Arguments
FMT <- c('Src', 'EL', 'I', 'Q', 'U', 'V', 'eI', 'eQ', 'eU', 'eV', 'EL')
FLDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(FLDF) <- FMT

flagNum <- list()
for(fileName in fileList){
	cat(fileName); cat('\n')
    fileLines <- removeBlank(readLines(fileName))
	DF <- readStokesSection(fileLines)
	DF$File <- fileName
	FLDF <- rbind(FLDF, DF)
}
FLDF$Src <- as.character(lapply(as.character(FLDF$Src), sourceMatch))
FLDF$P <- sqrt(FLDF$Q^2 + FLDF$U^2)
FLDF$eP <- sqrt(FLDF$eQ^2 + FLDF$eU^2)
FLDF$EVPA <- 0.5* atan2(FLDF$U, FLDF$Q)
FLDF$eEVPA <- FLDF$eP/FLDF$P
save(FLDF, file='Flux.Rdata')
BandName <- sprintf('Band-%d', as.numeric(strsplit(FLDF$File[1], '[-|.|_]')[[1]][8]))
#-------- Source List
sourceList <- unique(FLDF$Src)
numSrc <- length(sourceList)

if(0){
pdf('Flux.pdf', width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
for(src_index in 1:numSrc){
	DF <- FLDF[FLDF$Src == sourceList[src_index],]
	#-------- Plot Stokes I
	plot(DF$Date, DF$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='Total Flux Density', ylim=c(0, 1.2*max(DF$I)))
	color_vector <- rep("#000000FF", length(DF$Date))
	color_vector[which(DF$flagNum > 0)] <- "#0000FF3F"
	arrows(as.numeric(DF$Date), DF$I - DF$eI, as.numeric(DF$Date), DF$I + DF$eI, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$I, pch=20, col=color_vector)
	
	#-------- Plot polarized flux
	plot(DF$Date, DF$P, type='n', xlab='Date', ylab='Polarized Flux [Jy]', main='Polarized Flux Density', ylim=c(0, 1.2* max(DF$P)))
	arrows(as.numeric(DF$Date), DF$P - DF$eP, as.numeric(DF$Date), DF$P + DF$eP, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$P, pch=20, col=color_vector)
	
	#-------- Plot EVPA
	plot(DF$Date, DF$EVPA, type='n', xlab='Date', ylab='EVPA [deg]', main='Polarization Angle', ylim=c(-90,90))
	abline(h=90, lwd=0.1); abline(h=-90, lwd=0.1)
	arrows(as.numeric(DF$Date), 180*(DF$EVPA - DF$eEVPA)/pi, as.numeric(DF$Date), 180*(DF$EVPA + DF$eEVPA)/pi, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$EVPA*180/pi, pch=20, col=color_vector)
	
	mtext(side = 3, line=1, outer=T, text = sprintf('%s %s',sourceList[src_index], BandName), cex=2)
}
dev.off()
par(par.old)
}

medI <- medQ <- medU <- eI <- eQ <- eU <- numeric(0)
for(src_index in 1:numSrc){
	DF <- FLDF[FLDF$Src == sourceList[src_index],]
	medI[src_index] <- median(DF$I); medQ[src_index] <- median(DF$Q); medU[src_index] <- median(DF$U)
	if(length(DF$eI) == 1){
		eI[src_index] <- DF$eI; eQ[src_index] <- DF$eQ; eU[src_index] <- DF$eU
	} else {
		eI[src_index] <- sd(DF$eI); eQ[src_index] <- sd(DF$eQ); eU[src_index] <- sd(DF$eU)
	}
}
polDF <- data.frame( Src=as.character(sourceList), I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi )
polDF <- polDF[order(polDF$P, decreasing=T),]
#polDF <- polDF[order(polDF$Src, decreasing=F),]
rownames(polDF) <- c(1:nrow(polDF))
save(polDF, file='Pol.Rdata')

sourceList <- grep('^J[0-9]', polDF$Src)
#for(src_index in 1:numSrc){
#	cat(sprintf("%10s  %5.1f  %6.3f  %6.3f\n", polDF$Src[src_index], polDF$I[src_index], polDF$Q[src_index]/polDF$I[src_index], polDF$U[src_index]/polDF#$I[src_index]))
#}

# print(xtable(polDF, digits=3), include.rownames=F)

pdf(sprintf('Flux-%s.pdf', BandName), width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
for(source in sourceList){
	if(source[1] != 'J'){ next }
	DF <- FLDF[FLDF$Src == source,]
	pDF <- polDF[polDF$Src == source,]
	cat(sprintf("%10s  %5.1f  %6.3f  %6.3f\n", pDF$Src, pDF$I, pDF$p, pDF$EVPA))
	#-------- Plot Stokes I
	plot(DF$Date, DF$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='Total Flux Density', ylim=c(0, 1.2*max(DF$I)))
	color_vector <- rep("#000000FF", length(DF$Date))
	color_vector[which(DF$flagNum > 0)] <- "#0000FF3F"
	arrows(as.numeric(DF$Date), DF$I - DF$eI, as.numeric(DF$Date), DF$I + DF$eI, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$I, pch=20, col=color_vector)
	
	#-------- Plot polarized flux
	plot(DF$Date, DF$P, type='n', xlab='Date', ylab='Polarized Flux [Jy]', main='Polarized Flux Density', ylim=c(0, 1.2* max(DF$P)))
	arrows(as.numeric(DF$Date), DF$P - DF$eP, as.numeric(DF$Date), DF$P + DF$eP, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$P, pch=20, col=color_vector)
	
	#-------- Plot EVPA
	plot(DF$Date, DF$EVPA, type='n', xlab='Date', ylab='EVPA [deg]', main='Polarization Angle', ylim=c(-90,90))
	abline(h=90, lwd=0.1); abline(h=-90, lwd=0.1)
	arrows(as.numeric(DF$Date), 180*(DF$EVPA - DF$eEVPA)/pi, as.numeric(DF$Date), 180*(DF$EVPA + DF$eEVPA)/pi, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$EVPA*180/pi, pch=20, col=color_vector)
	
	mtext(side = 3, line=1, outer=T, text = sprintf('%s %s', source, BandName), cex=2)
}
dev.off()
par(par.old)

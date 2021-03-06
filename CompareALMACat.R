library(xtable)
#library(maptools)
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
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-C"){ catFile <- substring(args[index], 3) }
		if(substr(args[index], 1,2) == "-F"){ FLDFFile <- substring(args[index], 3)}
	}
	return(list(catFile = catFile, FLDFFile = FLDFFile))
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
	pointer <- grep("mean", Lines)
	numSource <- length(pointer)
	StokesI <- StokesQ <- StokesU <- StokesV <- numeric(4)
	errI <- errQ <- errU <- errV <- numeric(4)
	I <- Q <- U <- V <- numeric(0)
	eI <- eQ <- eU <- eV <- numeric(0)
	srcList <- character(0); EL <- numeric(0); freq <- numeric(4)
	for(srcIndex in 1:numSource){
		# if(length(grep('Only', Lines[(pointer[srcIndex]+1):(pointer[srcIndex]+7)])) > 0){ next }
		# if(length(grep('nan', Lines[(pointer[srcIndex]+1):(pointer[srcIndex]+7)])) > 0){ next }
		srcList <- append(srcList, sourceMatch(as.character(strsplit(Lines[pointer[srcIndex] - 8], '[ |=]+')[[1]][3])) )
		EL <- append(EL, as.numeric(strsplit(Lines[pointer[srcIndex] - 8], '[ |=]+')[[1]][5]))
		elements <- strsplit(Lines[pointer[srcIndex]], '[ |=|(|)|z]+')[[1]]
		FREQ <- as.numeric(elements[3])
		I <- append(I, as.numeric(elements[5])); eI <- append(eI, as.numeric(elements[6]))
		Q <- append(Q, as.numeric(elements[7])); eQ <- append(eQ, as.numeric(elements[8]))
		U <- append(U, as.numeric(elements[9])); eU <- append(eU, as.numeric(elements[10]))
		V <- append(V, as.numeric(elements[11])); eV <- append(eV, as.numeric(elements[12]))
	}
	return(data.frame(Src=as.character(srcList), EL=EL, I=I, Q=Q, U=U, V=V, eI=eI, eQ=eQ, eU=eU, ev=eV))
}

#-------- UTC of equalizer scan
eqUTC <- function( Lines ){
	equalizerPointer <- grep("Equalizer", Lines)
	return(as.Date(strptime(strsplit(Lines[equalizerPointer], ' ')[[1]][7], "%Y/%m/%d/%H:%M:%S", tz="UTC")))
}

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

#-------- Read ALMA source Catalog
readALMAcat <- function(Lines){
	scanNum <- length(Lines)
	scanDates <- as.Date(as.character(NULL))
	srcName <- as.character(NULL)
	I <- Ierr <- Freq <- numeric(0)
	line_index <- 0
	for(scan_index in 1:scanNum){
		lineElements <- strsplit(Lines[scan_index], '[ |(|)]+')[[1]]
		if(length(lineElements) < 10) next
		if(lineElements[8] != "GHz") next
		if(lineElements[9] != "ALMA") next
		line_index <- line_index + 1
		srcName[line_index] <- as.character(lineElements[10])
		I[line_index] <- as.numeric(lineElements[3])
		if( is.numeric(lineElements[5]) ){ Ierr[line_index] <- as.numeric(lineElements[5]) }
		else { Ierr[line_index] <- 0.2* I[line_index]}
		scanDates[line_index] <- strptime(lineElements[6], "%Y-%m-%d", tz="UTC")
		Freq[line_index] <- as.numeric(lineElements[7])
	}
	
	DF <- data.frame(Src=srcName, I=I, eI=Ierr, Date=scanDates, freq=Freq)
	return(DF[order(DF$Date),])
}
	

#-------- CheckTrex
checkTrec <- function(DF){
	#-------- Check Trec
	for(col_index in 2:9){ negTrecList <- which(TrsDF[[col_index]] < 0.0) }
	for(col_index in 2:9){ negTsysList <- which(TrsDF[[col_index+8]] < 0.0) }
	for(col_index in 2:9){ netTskyList <- which(TrsDF[[col_index+8]] < TrsDF[[col_index]]) }
	return( list(negTrecList, negTsysList, netTskyList) )
}

#-------- Filter catalog by nearest date
nearCat <- function( ObsDF, CatDF ){
	nearestIndex <- 1:length(ObsDF$Date)
	for(index in 1:length(ObsDF$Date)){
		nearestIndex[index] <- which.min( abs(CatDF$Date - ObsDF$Date[index]))[1]
	}
	return( CatDF[nearestIndex,] )
}

catBind <- function(DF){
	dateList <- as.Date(unique(DF$Date))
	dateNum <- length(dateList)
	I <- freq <- numeric(0)
	for(date_index in 1:dateNum){
		index <- which(DF$Date == dateList[date_index])
		I[date_index] <- median(DF$I[index])
		freq[date_index] <- median(DF$freq[index])
	}
	return(data.frame(Src = rep(DF$Src[1], dateNum), I=I, Date=dateList, freq = freq))
}


#-------- Start program
#argList <- parseArg(commandArgs(trailingOnly = T))
argList <- list(FLDFFile = 'Flux.Rdata', catFile = 'ALMAcat.log')
FLDFFile <- argList$FLDFFile
catFile <- argList$catFile
load(FLDFFile)
catDF <- readALMAcat(readLines(catFile))
save(catDF, file='catDF.Rdata')
#if(0){
pdf('FluxComp.pdf')
#-------- Compare between measurements and catalog
representativeFreq <- median(FLDF$Freq)
sourceList <- as.character(unique(FLDF$Src))
relDiffAccum <- numeric(0)
for(source in sourceList){
	cat(sprintf('Statistics for %s\n', source))
	#-------- find source in calibrator catalog
	catSrcDF <- catDF[(catDF$Src == sourceMatch(source)) & (abs(catDF$freq - representativeFreq) < 0.1*representativeFreq),]
	if( length(catSrcDF$I) < 1 ){ next }
	catSrcDF <- catBind(catSrcDF)
	FLSrcDF <- FLDF[FLDF$Src == source,]; FLSrcDF$Date <- as.Date(FLSrcDF$Date)
	plot(catSrcDF$Date, catSrcDF$I, type='n', xlim=c(min(FLSrcDF$Date), max(FLSrcDF$Date)), ylim=c(0.0, 1.2* max(FLSrcDF$I)), xlab='Date', ylab='Stokes I [Jy]', main=sprintf('%s (%.1f GHz)', source, representativeFreq))
	lines(catSrcDF$Date, catSrcDF$I, col='blue')
	points(catSrcDF$Date, catSrcDF$I, pch=20, col='blue', cex=0.5)
	points(as.Date(FLSrcDF$Date), FLSrcDF$I, pch=20, col='red', cex=0.1)
	#-------- Comparison
	nearestCatDF <- nearCat(FLSrcDF, catSrcDF)
	DateFilter1 <- which( abs(FLSrcDF$Date - nearestCatDF$Date) < 0.5 )
	points(FLSrcDF$Date[DateFilter1], FLSrcDF$I[DateFilter1], pch=20, col='red', cex=1)
	relativeDiff <- 100.0* abs(FLSrcDF$I[DateFilter1] - nearestCatDF$I[DateFilter1])/nearestCatDF$I[DateFilter1]
	relDiffAccum <- append(relDiffAccum, 100.0*(FLSrcDF$I[DateFilter1] - nearestCatDF$I[DateFilter1])/nearestCatDF$I[DateFilter1])
	text_sd <- sprintf('Diff. in the same day: Median = %4.1f%%', median(relativeDiff))
	text(min(FLSrcDF$Date), 0.1*max(FLSrcDF$I), text_sd, pos=4)
	#hist( 100.0* (FLSrcDF$I[DateFilter10] - nearestCatDF$I[DateFilter10])/nearestCatDF$I[DateFilter10], xlab='Relative Flux Diff. [%]' )
	#plot( abs(FLSrcDF$Date - nearestCatDF$Date), abs(FLSrcDF$I - nearestCatDF$I)/nearestCatDF$I, pch=20 )
}
H <- hist(relDiffAccum, xlab='Relative Flux Difference [%]', ylab='Frequency', main=sprintf('Diff (a priori cal. - ALMA_SC) / ALMA_SC (%.1f GHz)', representativeFreq), xlim=c(-50,50))
text_sd <- sprintf('Diff. (median) = %4.1f%%', median(abs(relDiffAccum))); text(15, max(H$counts), text_sd, cex=0.7, pos=4)
text_sd <- sprintf('[%4.1f%% : %4.1f%%] in 50%% probability', quantile(relDiffAccum, 0.25), quantile(relDiffAccum, 0.75)); text(15, 0.9* max(H$counts), text_sd, cex=0.7, pos=4)
text_sd <- sprintf('[%4.1f%% : %4.1f%%] in 90%% probability', quantile(relDiffAccum, 0.05), quantile(relDiffAccum, 0.95)); text(15, 0.8* max(H$counts), text_sd, cex=0.7, pos=4)
dev.off()
#}

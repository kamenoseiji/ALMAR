library(parallel)   # multicore parallelization
library(dplyr)
library(suncalc)
library(scales)
numCore = detectCores()
Sys.setenv(TZ="UTC")
ALMA_POS <- matrix(c( -67.755, -23.029 ), nrow=1 )
SSOlist = c('Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus')
#-------- Parse arguments
parseArg <- function( args ){
    fileDF <- read.table(args[1])
    return( as.character(fileDF[[1]]) )
}
#-------- Receiver Band
getBand <- function(fileList){
    return(as.integer(mclapply(fileList, function(fileName){
            bandPointer <- as.integer(regexpr("RB_[0-10]", fileName)[1])
            return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
        },mc.cores=numCore)))
}
#-------- Read UID
readUID <- function( Lines ){
    UIDpointer <- grep("uid___", Lines)
    return( strsplit(Lines[UIDpointer], '[ |,]')[[1]][2] )
}
readTsysDigital <- function( Lines ){
    UIDpointer <- grep("TsysDigitalCorrection", Lines)
    return( ifelse(length(UIDpointer) > 0, strsplit(Lines[UIDpointer], 'TsysDigitalCorrection ')[[1]][2], 'OFF') )
}
#-------- Find Calibrator name
findCalibrator <- function( Lines ){
    SSOListPointer <- grep("Flux Calibrator is", Lines)
    if( length(SSOListPointer) > 0){
        fluxCalName <- strsplit(Lines[SSOListPointer], ' ')[[1]][4]
        scalerUTC <- strptime(strsplit(Lines[SSOListPointer], ' ')[[1]][6], "%Y/%m/%d/%H:%M:%S")
        sunsetUTC <- getSunlightTimes(as.Date(scalerUTC), lat=ALMA_POS[2], lon=ALMA_POS[1])['sunset'][[1]]
        datePointer <- grep(fluxCalName, Lines)[2]
        EL <- as.numeric(strsplit(Lines[datePointer], '[ |=]+')[[1]][5])
		return(list(calibrator=fluxCalName, EL=EL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24))
    } else {
        SSOListPointer <- grep("Aeff", Lines)
    }
    if( length(SSOListPointer) == 0){ return(list(0)) }
    UsedSSOList <- vector(mode='character')
    for(SSO in SSOlist){ if(length(grep(SSO, Lines[SSOListPointer])) > 0){ UsedSSOList[length(UsedSSOList)+1]  <- SSO}}
    if( length(UsedSSOList) == 0){ return(list(0)) }
    for(SSO in UsedSSOList){
        datePointer <- grep(SSO, Lines[4:length(Lines)]) + 3
        fluxCalName <- strsplit(Lines[datePointer], '[ |=]+')[[1]][3]
        EL <- as.numeric(strsplit(Lines[datePointer], '[ |=]+')[[1]][5])
	    scalerUTC <- strptime(strsplit(Lines[datePointer], '[ |=]+')[[1]][7], "%Y/%m/%d/%H:%M:%S")
    	sunsetUTC <- getSunlightTimes(as.Date(scalerUTC), lat=ALMA_POS[2], lon=ALMA_POS[1])['sunset'][[1]]
		return(list(calibrator=fluxCalName, EL=EL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24))
    }
}
#-------- Read Aeff Section
readAeffSection <- function(Lines){
	pointer <- grep(" Aeff", Lines) + 1
    spwLabel <- strsplit(Lines[pointer - 1], '[ |:]+')[[1]]; spwLabel <- spwLabel[(which(spwLabel=='Aeff')+1):length(spwLabel)]
    spwLabel <- spwLabel[grep('SPW', spwLabel)] 
    XspwList <- grep("-X", spwLabel) + 2; YspwList <- grep("-Y", spwLabel) + 2
    FMT <- c('Ant', 'AeX1', 'AeY1', 'AeX2', 'AeY2', 'AeX3', 'AeY3', 'AeX4', 'AeY4')
    DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
	while(length(grep('[--|%]',Lines[pointer])) > 0){
	    if( strsplit(Lines[pointer], ' ')[[1]][2] != ':'){ break }
        if( length(grep('--', Lines[pointer])) != 0 ){pointer <- pointer + 1; next}
        tmpDF <- data.frame(
            Ant = strsplit(Lines[pointer], ' ')[[1]][1],
            #AeX = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList])),
            #AeY = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList])))
            AeX1 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList[1]]),
            AeY1 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList[1]]),
            AeX2 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList[2]]),
            AeY2 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList[2]]),
            AeX3 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList[3]]),
            AeY3 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList[3]]),
            AeX4 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList[4]]),
            AeY4 = as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList[4]]))
        DF <- rbind(DF, tmpDF)
		pointer <- pointer + 1
	}
	return(DF)
}
#-------- Read D-term Section
readDtermSection <- function(Lines){
    FMT <- c('Ant', 'Dx1', 'Dy1', 'Dx2', 'Dy2', 'Dx3', 'Dy3', 'Dx4', 'Dy4', 'Date')
    #---- Date
	pointer <- grep("EL=", Lines)[1]
    BP_UTC <- strptime(strsplit(Lines[pointer], '[ |=]+')[[1]][7], "%Y/%m/%d/%H:%M:%S")
    emptyDF <- data.frame(Ant=NA, Dx1=0.0+0.0i, Dy1=0.0+0.0i, Dx2=0.0+0.0i, Dy2=0.0+0.0i, Dx3=0.0+0.0i, Dy3=0.0+0.0i, Dx4=0.0+0.0i, Dy4=0.0+0.0i, Date=BP_UTC)
    #---- D-term
	pointer <- grep("D-term", Lines)
    if(length(pointer) == 0){ return(emptyDF)}
    spwLabel <- strsplit(Lines[pointer], '[ |:]+')[[1]]
    if( length(spwLabel) < 5){ return(emptyDF) }
    spwLabel <- strsplit(Lines[pointer + 1], '[ |:]+')[[1]]
    spwList <- grep('[Dx|Dy]', spwLabel)
    pointer <- pointer + 2
    while( is.na(strsplit(Lines[pointer], ' ')[[1]][1]) == F){
        Dvec  <- as.complex(strsplit(Lines[pointer], ' +')[[1]][spwList])
        tempDF <- data.frame(Ant = strsplit(Lines[pointer], ' +')[[1]][1])
        tempDF[1,2:9] <- Dvec
        tempDF$Date <- BP_UTC
        colnames(tempDF) <- FMT
        if(is.na(emptyDF[[1]])[1]){
            emptyDF <- tempDF
        } else {
            emptyDF <- rbind(emptyDF, tempDF)
        }
        pointer <- pointer + 1
    }
	return(emptyDF)
}    
#-------- Ae correction
AeCorrect <- function(DF, band=3, thresh=10){
    bandDF <- DF[DF$Band == band,]  # band selection
    bandDF <- bandDF[((abs(bandDF$AeX - median(bandDF$AeX)) < thresh) & (abs(bandDF$AeY - median(bandDF$AeY)) < thresh)),]  # Filter irregular Aeff
    bandDF$Ae <- 0.5* (bandDF$AeX + bandDF$AeY) # combine polarization
    bandDF$sunsetArg <- sinpi((bandDF$sunset + 4.5)/12.0)   # phase since sunset+4.5h
    bandDF$EL45 <- bandDF$EL - 45.0                         # EL reference of 45 deg
    fit <- lm(formula = Ae ~ EL45 + sunsetArg, data=bandDF)   # regression for EL and sunset phase
    AeC <- coef(fit)['EL45'][[1]]* bandDF$EL45 + coef(fit)['sunsetArg'][[1]]* bandDF$sunsetArg  # correction factor
    bandDF$AeC <- bandDF$Ae - AeC       # apply correction
    bandDF$AeX <- bandDF$AeX - AeC      # apply correction
    bandDF$AeY <- bandDF$AeY - AeC      # apply correction
    #-------- correction for different calibrator
    AeRef <- median(bandDF$AeC)
    AeC <- rep(1.0, length(SSOlist))
    names(AeC) <- SSOlist
    for(SSO in SSOlist){
        calDF <- bandDF[bandDF$calibrator == SSO,]
        if(nrow(calDF) < 10){ next }
        AeC[SSO] <- AeRef / median(calDF$AeC)   # correction factor for each calibrator
        text_sd <- sprintf('%10s  : %.3f (%.3f)\n', SSO, median(calDF$AeC), sd(calDF$AeC))
        cat(text_sd)
        bandDF[bandDF$calibrator == SSO,]$AeC <- bandDF[bandDF$calibrator == SSO,]$AeC * AeC[SSO]
        bandDF[bandDF$calibrator == SSO,]$AeX <- bandDF[bandDF$calibrator == SSO,]$AeX * AeC[SSO]
        bandDF[bandDF$calibrator == SSO,]$AeY <- bandDF[bandDF$calibrator == SSO,]$AeY * AeC[SSO]
    }
    bandDF$AeR <- bandDF$AeY / bandDF$AeX
    bandDF$Ae  <- 0.5*(bandDF$AeY + bandDF$AeX)
    write.table(data.frame(AeC), quote=F, col.names=F, sep=',', file=sprintf('SSO.B%d.table', band))
    return( bandDF )
}
#-------- Monthly smoothing
SPL_period <- function(DF, refPeriod, weight=c(0,0)){
    # DF <- data.frame(secSinceRefTIme, Value)
    # refPeriod <- secSinceRefTime to output
    Cadence <- 2* max(diff(refPeriod))
    refTime <- max(refPeriod)
    relSec <- c(min(refPeriod)-Cadence, DF[[1]], refTime+Cadence)
    Value  <- c(median(DF[[2]]), DF[[2]], median(DF[[2]]))
    if(mean(weight) == 0){
        weight <- rep(1.0, length(Value))
    } else {
        weight <- c(1.0, weight, 1.0)
    }
    NumKnots <- round(diff(range(relSec)) / max(diff(relSec)))
    #if(length(Value) < 10){
    if(NumKnots < 4){
        return( data.frame(Date=refPeriod, Value=median(Value) ))
    }
    if(length(Value) > 2*NumKnots){
        SPL <- smooth.spline(relSec, Value, all.knots=F, nknots=NumKnots, w=weight)
    } else {
        SPL <- smooth.spline(relSec, Value, all.knots=F, spar=0.5, w=weight)
    }
    return( data.frame(Date=refPeriod, Value=predict(SPL, refPeriod)$y ))
}
#-------- Read Aeff section and store into a data frame
Log2Aeff <- function(fileName){
    #FMT <- c('Ant', 'AeX', 'AeY', 'Band', 'calibrator', 'EL', 'Date', 'sunset', 'UID', 'TsysCorr')
    FMT <- c('Ant', 'AeX1', 'AeY1', 'AeX2', 'AeY2','AeX3', 'AeY3','AeX4', 'AeY4','Band', 'calibrator', 'EL', 'Date', 'sunset', 'UID', 'TsysCorr')
    fileLines <- readLines(fileName)
    emptyDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]
    colnames(emptyDF) <- FMT
    if(length(fileLines) < 10){ return() }
    if(fileLines[1] == ''){ return() }
    TsysCorr <- readTsysDigital(fileLines)
    #if(TsysCorr == 'OFF'){ return(emptyDF) }
    # UID <- readUID(fileLines)
    UID <- strsplit(fileName, '[ |-]')[[1]][1]
	CalList <- findCalibrator(fileLines)
    if(CalList[[1]][1] == 0){ return() }
	emptyDF  <- readAeffSection(fileLines)
    emptyDF$Band <- getBand(fileName)
    emptyDF$UID <- UID
    emptyDF$TsysCorr <- TsysCorr
    if(CalList[[1]] > 0){
	    emptyDF$calibrator <- CalList$calibrator
	    emptyDF$EL <- CalList$EL
	    emptyDF$Date <- CalList$UTC
	    emptyDF$sunset <- CalList$sunset
    } else {
	    emptyDF$calibrator <- ''
	    emptyDF$EL <- NA
	    emptyDF$Date <- CalList$UTC
	    emptyDF$sunset <- CalList$sunset
    }
    return(emptyDF)
}
#-------- Read D-term section and store into a data frame
Log2Dterm <- function(fileName){
    FMT <- c('Ant', 'Dx1', 'Dy1', 'Dx2', 'Dy2', 'Dx3', 'Dy3', 'Dx4', 'Dy4', 'Date', 'File')
    fileLines <- readLines(fileName)
    emptyDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]
    colnames(emptyDF) <- FMT
    if(length(fileLines) < 10){ return() }
    if(fileLines[1] == ''){ return() }
    emptyDF <- readDtermSection(fileLines)
    emptyDF$File <- fileName
    return(emptyDF)
}

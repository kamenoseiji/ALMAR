library(parallel)   # multicore parallelization
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
	if( length(SSOListPointer) == 0){ SSOListPointer <- grep(" Aeff", Lines) }
	if( length(SSOListPointer) == 0){ return(list(0))}
    UsedSSOList <- SSOlist[SSOlist %in% strsplit(Lines[SSOListPointer], '[ |,]')[[1]]]
    if(exists('scalerUTC')){ rm(scalerUTC) }
    for(SSO in UsedSSOList){
        if(exists('scalerUTC')){ break }
        SSOpointer <- grep(paste(SSO, 'EL'), Lines)
        if(length(SSOpointer) == 0){ next}
        datePointer <- SSOpointer
        while(length(grep('mean', Lines[datePointer]) ) == 0){ datePointer <- datePointer + 1}
        lineData <- strsplit(Lines[datePointer], '[ |(|)|z]+')[[1]]
        StokesI <- as.numeric(lineData[5])
        errI    <- as.numeric(lineData[6])
        if(StokesI / errI < 5.0){ next }
        fluxCalName <- strsplit(Lines[SSOpointer], '[ |=]+')[[1]][3]
        EL <- as.numeric(strsplit(Lines[SSOpointer], '[ |=]+')[[1]][5])
	    scalerUTC <- strptime(strsplit(Lines[SSOpointer], '[ |=]+')[[1]][7], "%Y/%m/%d/%H:%M:%S")
    }
	if(exists('scalerUTC')){
    	sunsetUTC <- getSunlightTimes(as.Date(scalerUTC), lat=ALMA_POS[2], lon=ALMA_POS[1])['sunset'][[1]]
		return(list(calibrator=fluxCalName, EL=EL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24))
	} else {
		return(list(0))
	}
}
#-------- Read Aeff Section
readAeffSection <- function(Lines){
	pointer <- grep(" Aeff", Lines) + 1
    spwLabel <- strsplit(Lines[pointer - 1], '[ |:]+')[[1]]; spwLabel <- spwLabel[(which(spwLabel=='Aeff')+1):length(spwLabel)]
    XspwList <- grep("-X", spwLabel) + 2; YspwList <- grep("-Y", spwLabel) + 2
    FMT <- c('Ant', 'AeX', 'AeY')
    DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
	while(length(grep('%', Lines[pointer])) > 0){
        tmpDF <- data.frame(
            Ant = strsplit(Lines[pointer], ' ')[[1]][1],
            AeX = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList])),
            AeY = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList])))
        DF <- rbind(DF, tmpDF)
		pointer <- pointer + 1
	}
	return(DF)
}
#-------- Read Aeff Section from a priori calibration
readAeApriori <- function(Lines){
    pointer <- grep("Use J", Lines)[2]
    FMT <- c('Ant', 'AeX', 'AeY', 'Band', 'calibrator', 'EL', 'Date', 'sunset', 'UID', 'TsysCorr')
    DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
    if( length(grep('^[A-Z]', Lines[pointer+1])) == 0 ){ return(DF) }
    strList = strsplit(Lines[pointer], '[ |=]+')[[1]]
    tmpDate <- strptime(strList[7], "%Y/%m/%d/%H:%M:%S")
    pointer <- pointer + 1
	while( length(grep('^[A-Z]', Lines[pointer])) > 0){
        antStrList = strsplit(Lines[pointer], ' +')[[1]]
        tmpDF <- data.frame(
            Ant = antStrList[1],
            AeX = as.numeric(antStrList[2]),
            AeY = as.numeric(antStrList[3]),
            calibrator = strList[2],
            EL = as.numeric(strList[5]),
            Date = tmpDate,
            sunset = as.numeric(difftime(tmpDate, getSunlightTimes(as.Date(tmpDate), lat=ALMA_POS[2], lon=ALMA_POS[1])['sunset'][[1]]))%%24)
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
    cat(fileName)
    FMT <- c('Ant', 'AeX', 'AeY', 'Band', 'calibrator', 'EL', 'Date', 'sunset', 'UID', 'TsysCorr')
    fileLines <- readLines(fileName)
    emptyDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]
    colnames(emptyDF) <- FMT
    if(length(fileLines) < 10){ return(emptyDF) }
    if(fileLines[1] == ''){ return(emptyDF) }
    TsysCorr <- readTsysDigital(fileLines)
    if(TsysCorr == 'OFF'){ return(emptyDF) }
    UID <- readUID(fileLines)
    if(fileLines[1] < 10){ return(emptyDF) }
	CalList <- findCalibrator(fileLines)
    if(CalList[[1]] == 0){ return(emptyDF) }
    if(CalList[[1]] > 0){
	    emptyDF  <- readAeffSection(fileLines)
        emptyDF$Band <- getBand(fileName)
	    emptyDF$calibrator <- CalList$calibrator
	    emptyDF$EL <- CalList$EL
	    emptyDF$Date <- CalList$UTC
	    emptyDF$sunset <- CalList$sunset
        emptyDF$UID <- UID
        emptyDF$TsysCorr <- TsysCorr
    } else {
        emptyDF <- readAeApriori(fileLines)
        emptyDF$Band <- getBand(fileName)
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
#-------- Start program
Arguments <- commandArgs(trailingOnly = T)
fileList <- parseArg(Arguments)
#fileList <- parseArg('fileList')
#AeDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(AeDF) <- FMT
DFList <- mclapply(fileList, Log2Aeff, mc.cores=numCore)
AeDF <- do.call("rbind", DFList)
save(AeDF, file='AeDF.Rdata')
AeDF <- AeDF[complete.cases(AeDF$AeX),]
AeDF <- AeDF[complete.cases(AeDF$AeY),]
AeDF <- AeDF[((AeDF$AeX > 25.0) & (AeDF$AeX < 100.0) & (AeDF$AeY > 25.0) & (AeDF$AeY < 100.0)),]
AeDF <- AeDF[complete.cases(AeDF$Date),]
AeDF <- na.omit(AeDF)
AeDF <- AeDF[order(AeDF$Date), ]
#save(AeDF, file='AeDF.Rdata')

DFList <- mclapply(fileList, Log2Dterm, mc.cores=numCore)
DtermDF <- do.call("rbind", DFList)
DtermDF <- DtermDF[complete.cases(DtermDF$Date),]
DtermDF <- na.omit(DtermDF)
DtermDF <- DtermDF[order(DtermDF$Date), ]
DtermDF$Band <- getBand(DtermDF$File)
save(DtermDF, file='Dterm.Rdata')

#-------- Ae table
BandList <- c(1, 3, 4, 6, 7)
pcolors <- c(alpha('orange', 0.5), alpha('purple', 0.5))
lcolors <- c('orange', 'purple')
refTime <- max(AeDF$Date)
MonthSec <- 2629744
refPeriod <- seq(as.numeric(difftime(min(AeDF$Date), refTime, units='sec')), MonthSec, by=MonthSec)
for(Band in BandList){
    BandAeDF <- AeCorrect(AeDF, Band, 10)
    antList <- sort(as.character(unique(BandAeDF$Ant)))
    bandAeDF <- data.frame(Date = c('mean', 'sd', as.character(as.Date(refTime + refPeriod))))
    for(ant in antList){
        BandAntAeDF <- BandAeDF[BandAeDF$Ant == ant,]
        Ae  <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntAeDF$Date, refTime, units='sec')), Value=BandAntAeDF$Ae), refPeriod)
        XYR <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntAeDF$Date, refTime, units='sec')), Value=BandAntAeDF$AeR), refPeriod)
        AeX <- Ae$Value / sqrt(XYR$Value)
        AeY <- Ae$Value * sqrt(XYR$Value) 
        bandAeDF[[paste(ant , '-X', sep='')]] <- c(mean(BandAntAeDF$AeC / sqrt(BandAntAeDF$AeR)), sd(BandAntAeDF$AeC / sqrt(BandAntAeDF$AeR)), AeX)
        bandAeDF[[paste(ant , '-Y', sep='')]] <- c(mean(BandAntAeDF$AeC * sqrt(BandAntAeDF$AeR)), sd(BandAntAeDF$AeC * sqrt(BandAntAeDF$AeR)), AeY)
        pdf(sprintf('Ae-%s-B%d.pdf', ant, Band))
        plot(BandAeDF$Date, BandAeDF$Ae, type='n', ylim=c(0, 100), xlab='Date', ylab='Aeff (%)', main=sprintf('%s Band%d', ant, Band))
        lines(Ae$Date + refTime, Ae$Value / sqrt(XYR$Value), col=lcolors[1], lwd=2)
        lines(Ae$Date + refTime, Ae$Value * sqrt(XYR$Value), col=lcolors[2], lwd=2)
        points(BandAntAeDF$Date, BandAntAeDF$Ae / sqrt(BandAntAeDF$AeR), pch=20, cex=0.5, col=pcolors[1])
        points(BandAntAeDF$Date, BandAntAeDF$Ae * sqrt(BandAntAeDF$AeR), pch=20, cex=0.5, col=pcolors[2])
        legend("bottomleft", legend=c('Pol-X', 'Pol-Y'), col=lcolors, pch=rep(20, 2), lty=rep(1,2))
        dev.off()
    }
    write.table(format(bandAeDF, digits=4), file=sprintf('AeB%d.table', Band), quote=F, row.names=F)
}
#-------- Dterm table
refTime <- max(DtermDF$Date)
refPeriod <- seq(as.numeric(difftime(min(DtermDF$Date), refTime, units='sec')), MonthSec, by=MonthSec)
for(Band in BandList){
    BandDdf <- DtermDF[DtermDF$Band == Band,]
    antList <- sort(as.character(unique(BandDdf$Ant)))
    for(ant in antList){
        BandAntdDF <- BandDdf[BandDdf$Ant == ant,]
        bandAntDdf <- data.frame(Date = as.Date(refTime + refPeriod))
        if(nrow(BandAntdDF) < 3){
            for(BB in c(1,2,3,4)){
                for(pol in c('x','y')){
                    bandAntDdf[[sprintf('%s-BB%d-D%s', ant, BB, pol)]] <- (0.0 + 0.0i)
                }
            }
        } else {
            pdf(sprintf('Dterm.B%d.%s.pdf', Band, ant), width=8, height=11)
            par.old <- par(no.readonly=TRUE)
            par(mfrow=c(4,2), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
            column_index <- 1
            for(BB in c(1,2,3,4)){
                for(pol in c('x','y')){
                    colName <- sprintf('D%s%d', pol, BB)
                    ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    bandAntDdf[[sprintf('%s-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
                    #---- plot
                    column_index <- column_index + 1
                    plot(as.Date(BandAntdDF$Date), Re(BandAntdDF[[column_index]]), pch=21, cex=0.2, ylim=c(-0.1, 0.1), xlab='Date', ylab='D-term', main=colnames(BandAntdDF[column_index]), col=pcolors[1] )
                    points(as.Date(BandAntdDF$Date), Im(BandAntdDF[[column_index]]), pch=21, cex=0.2, col=pcolors[2] )
                    lines(as.Date(bandAntDdf[[1]]), Re(bandAntDdf[[column_index]]), col=lcolors[1], lwd=2)
                    lines(as.Date(bandAntDdf[[1]]), Im(bandAntDdf[[column_index]]), col=lcolors[2], lwd=2)
                    legend("bottomleft", legend=c('Real', 'Imag'), col=lcolors, pch=rep(20, 2), lty=rep(1,2))
                }
            }
            mtext(side = 3, line=1, outer=T, text = sprintf('D-term %s Band-%d', ant, Band), cex=2)
            par(par.old)
            dev.off()
        }
        write.table(format(bandAntDdf, digits=6), file=sprintf('DtermB%d.%s.table', Band, ant), quote=F, row.names=F)
    }
    #---- Dummy antenna
    for(ant in c('CM', 'PM', 'DA', 'DV')){
        BandAntdDF <- subset(BandDdf, grepl(ant,  BandDdf$Ant))
        bandAntDdf <- data.frame(Date = as.Date(refTime + refPeriod))
        for(BB in c(1,2,3,4)){
            for(pol in c('x','y')){
                colName <- sprintf('D%s%d', pol, BB)
                ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                bandAntDdf[[sprintf('%s00-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
            }
        }
        write.table(format(bandAntDdf, digits=6), file=sprintf('DtermB%d.%s00.table', Band, ant), quote=F, row.names=F)
    }
}

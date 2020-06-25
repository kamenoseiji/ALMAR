library(suncalc)
library(scales)
ALMA_POS <- matrix(c( -67.755, -23.029 ), nrow=1 )
SSOlist = c('Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus')
#-------- Parse arguments
parseArg <- function( args ){
    fileDF <- read.table(args[1])
    return( as.character(fileDF[[1]]) )
}
#-------- Find Calibrator name
findCalibrator <- function( Lines ){
	if( length(grep(" Aeff", Lines)) == 0){ return(list(-1))}
    bestR <- -1.0
    for(SSO in SSOlist){
        SSOpointers <- grep(paste(SSO, 'EL'), Lines)
        for( SSOpointer in SSOpointers ){
            if(length(grep('Model I', Lines[SSOpointer+1])) == 0){ next }
        }
        if(length(grep('Model I', Lines[SSOpointer+1])) == 0){ next }
        if(length(SSOpointer) > 0){
            datePointer <- SSOpointer + 3
            StokesI <- modelI <- numeric(0)
            while(length(grep('SPW', Lines[datePointer]) ) > 0){
                lineData <- strsplit(Lines[datePointer], '[ |(|)|z]+')[[1]]
                if(length(lineData) < 13){ datePointer <- datePointer + 1; next }
                StokesI <- append(StokesI, as.numeric(lineData[5]))
                modelI  <- append(modelI,  as.numeric(lineData[13]))
                datePointer <- datePointer + 1
            }
            if(length(StokesI) < 3){ next }
            fluxR <- median(StokesI) / median(modelI)
            if( abs(fluxR - 1.0) < abs(bestR - 1.0) ){
                bestR <- fluxR
                fluxCalName <- strsplit(Lines[SSOpointer], '[ |=]+')[[1]][3]
                EL <- as.numeric(strsplit(Lines[SSOpointer], '[ |=]+')[[1]][5])
	            scalerUTC <- strptime(strsplit(Lines[SSOpointer], '[ |=]+')[[1]][7], "%Y/%m/%d/%H:%M:%S", tz="UTC")
            }
        }
    }
    if( abs(bestR - 1.0) > 0.5 ){ return(list(-1))}
    sunsetUTC <- getSunlightTimes(as.Date(scalerUTC), lat=ALMA_POS[2], lon=ALMA_POS[1])['sunset'][[1]]
	return(list(calibrator=fluxCalName, EL=EL, UTC=scalerUTC, sunset=as.numeric(scalerUTC-sunsetUTC)%%24, fluxR=bestR))
}
#-------- Read Aeff Section
readAeffSection <- function(Lines){
	pointer <- grep(" Aeff", Lines) + 1
    spwLabel <- strsplit(Lines[pointer - 1], '[ |:]+')[[1]]; spwLabel <- spwLabel[(which(spwLabel=='Aeff')+1):length(spwLabel)]
    XspwList <- grep("-X", spwLabel) + 2; YspwList <- grep("-Y", spwLabel) + 2
    FMT <- c('Ant', 'AeX', 'AeY')
    DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
	while( is.na(strsplit(Lines[pointer], ' ')[[1]][1]) == F){
        tmpDF <- data.frame(
            Ant = strsplit(Lines[pointer], ' ')[[1]][1],
            AeX = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][XspwList])),
            AeY = median(as.numeric(strsplit(Lines[pointer], '[ |%]+')[[1]][YspwList])))
        DF <- rbind(DF, tmpDF)
		pointer <- pointer + 1
	}
	return(DF)
}
#-------- Read D-term Section
readDtermSection <- function(Lines){
	pointer <- grep("D-term", Lines)
    if(length(pointer) == 0){ return(-1)}
    spwLabel <- strsplit(Lines[pointer], '[ |:]+')[[1]]
    if( length(spwLabel) < 5){ return(-1) }
    spwLabel <- strsplit(Lines[pointer + 1], '[ |:]+')[[1]]
    spwList <- grep('[Dx|Dy]', spwLabel)
    pointer <- pointer + 2
    FMT <- c('Ant', 'Dx1', 'Dy1', 'Dx2', 'Dy2', 'Dx3', 'Dy3', 'Dx4', 'Dy4')
    DF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DF) <- FMT
    while( is.na(strsplit(Lines[pointer], ' ')[[1]][1]) == F){
        Dline <- gsub('i', 'i ', Lines[pointer])
        Dvec  <- as.complex(strsplit(Dline, ' +')[[1]][spwList])
        tmpDF <- data.frame(Ant = strsplit(Dline, ' +')[[1]][1])
        tmpDF[1,2:9] <- Dvec
        DF <- rbind(DF, tmpDF)
        pointer <- pointer + 1
    }
    names(DF) <- FMT
	return(DF)
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
        text_sd <- sprintf('%10s  : %.2f (%.2f)\n', SSO, median(calDF$AeC), sd(calDF$AeC))
        cat(text_sd)
        bandDF[bandDF$calibrator == SSO,]$AeC <- bandDF[bandDF$calibrator == SSO,]$AeC * AeC[SSO]
        bandDF[bandDF$calibrator == SSO,]$AeX <- bandDF[bandDF$calibrator == SSO,]$AeX * AeC[SSO]
        bandDF[bandDF$calibrator == SSO,]$AeY <- bandDF[bandDF$calibrator == SSO,]$AeY * AeC[SSO]
    }
    bandDF$AeR <- bandDF$AeY / bandDF$AeX
    bandDF$Ae  <- 0.5*(bandDF$AeY + bandDF$AeX)
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
#-------- Start program
#Arguments <- commandArgs(trailingOnly = T)
#fileList <- parseArg(Arguments)
fileList <- parseArg('fileList')
FMT <- c('Ant', 'AeX', 'AeY', 'Band', 'calibrator', 'EL', 'Date', 'sunset', 'fluxR')
AeDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(AeDF) <- FMT
FMT <- c('Ant', 'Dx1', 'Dy1', 'Dx2', 'Dy2', 'Dx3', 'Dy3', 'Dx4', 'Dy4')
DtermDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(DtermDF) <- FMT
for(fileName in fileList){
	cat(fileName); cat('\n')
    fileLines <- readLines(fileName)
	CalList <- findCalibrator(fileLines)
    if(CalList[[1]] == -1){ next }
	DF  <- readAeffSection(fileLines)
	DF$Band <- as.numeric(strsplit(fileName, '_+|-')[[1]][6])
	DF$calibrator <- CalList$calibrator
	DF$EL <- CalList$EL
	DF$Date <- CalList$UTC
	DF$sunset <- CalList$sunset
	DF$fluxR <- CalList$fluxR
	DF$File <- fileName
	AeDF <- rbind(AeDF, DF)
	Ddf <- readDtermSection(fileLines)
    if( length(attributes(Ddf)) != 0){
	    Ddf$Band <- DF$Band
	    Ddf$calibrator <- DF$calibrator
	    Ddf$EL <- DF$EL
	    Ddf$Date <- DF$Date
	    Ddf$sunset <- DF$sunset
	    Ddf$fluxR <- DF$fluxR
	    Ddf$File <- DF$File
	    DtermDF <- rbind(DtermDF, Ddf)
    }
}
AeDF <- AeDF[order(AeDF$Date), ]
DtermDF <- DtermDF[order(DtermDF$Date), ]
save(AeDF, file='AeDF.Rdata')
save(DtermDF, file='Dterm.Rdata')
#-------- Ae table
BandList <- c(3, 6, 7)
pcolors <- c(alpha('orange', 0.5), alpha('purple', 0.5))
lcolors <- c('orange', 'purple')
#refTime <- Sys.time()
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
#refTime <- Sys.time()
refTime <- max(DtermDF$Date)
refPeriod <- seq(as.numeric(difftime(min(DtermDF$Date), refTime, units='sec')), MonthSec, by=MonthSec)
for(Band in BandList){
    BandDdf <- DtermDF[DtermDF$Band == Band,]
    antList <- sort(as.character(unique(BandAeDF$Ant)))
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
            for(BB in c(1,2,3,4)){
                for(pol in c('x','y')){
                    colName <- sprintf('D%s%d', pol, BB)
                    ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    bandAntDdf[[sprintf('%s-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
                }
            }
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
                ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod)
                ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod)
                bandAntDdf[[sprintf('%s00-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
            }
        }
        write.table(format(bandAntDdf, digits=6), file=sprintf('DtermB%d.%s00.table', Band, ant), quote=F, row.names=F)
    }
}

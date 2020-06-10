library(suncalc)
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
save(AeDF, file='AeDF.Rdata')
save(DtermDF, file='Dterm.Rdata')

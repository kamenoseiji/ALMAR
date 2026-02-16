library(parallel)   # multicore parallelization
library(dplyr)
library(VGAM)       # for Rice distribution
Sys.setenv(TZ="UTC")
sysIerr <- 0.005       # temporal Stokes I systematic error
sysPerr <- 0.003       # temporal polarization systematic error
minAntNum <- 5		   # Minimum number of antennas
numCore = detectCores()
cat(sprintf('Number of cores = %d\n', numCore))
sourceMatch <- function(sourceName){
	sourceDict <- list(
		c('J0006-0623', 'J0006-063'),
		c('J0006-0623', '0006-063'),
		c('J0237+2848', 'J0237+288'),
		c('J0238+1636', 'J0238+166'),
		c('J0319+4130', '3c84'),
		c('J0334-4008', 'J0334-401'),
		c('J0423-0120', 'J0423-013'),
		c('J0457+2624', 'J0457+2624_ALMA'),
		c('J0510+1800', 'J0510+180'),
		c('J0519-4546', 'J0519-454'),
		c('J0522-3627', 'J0522-364'),
		c('J0538-4405', 'J0538-440'),
		c('J0750+1231', 'J0750+125'),
		c('J0854+2006', 'J0854+201'),
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
		c('J1819-0258', 'J181917-025807'), 
		c('J1924-2914', 'J1924-292'),
		c('J2025+3343', 'J2025+337'),
		c('J2056-4714', 'J2056-472'),
		c('J2148+0657', 'J2148+069'),
		c('J2232+1143', 'J2232+117'),
		c('J2253+1608', '3c454.3'),
		c('J2258-2758', 'J2258-279'),
		c('Uranus', 'Uranus_1'))
	#
	for(index in 1:length(sourceDict)){ if( !is.na(match(sourceName, sourceDict[[index]]))){ return(sourceDict[[index]][1])} }
	return(sourceName)
}
#-------- Parse arguments
parseArg <- function( args ){
    fileDF <- read.table(args[1])
    return( as.character(fileDF[[1]]) )
}
parseScan <- function(scan){return(as.numeric(strsplit(scan, '[ |(|)|z]+')[[1]][c(3,5,6,7,8,9,10,11,12)]))} # Stokes parameter entry line format
#-------- Read Stokes Parameters
readStokesSection <- function(fileName){
    fileLines <- removeBlank(readLines(fileName, skipNul=TRUE))
    #---- Scan Entries
    scanPointer   <- grep('EL=', fileLines)
    scanEntryList <- fileLines[scanPointer]
    srcVec <- as.character(sapply(scanEntryList, function(scan){ return(strsplit(scan, '[ |=]+')[[1]][3]) }))
    ELVec  <- as.numeric(sapply(scanEntryList, function(scan){ return(strsplit(scan, '[ |=]+')[[1]][5]) }))
    UTCVec <- lapply(scanEntryList, function(scan){ return(  strptime(strsplit(scan, '[ |=]+')[[1]][7], "%Y/%m/%d/%H:%M:%S", tz='UTC'))})
    #---- Stokes Parameters
    StokesPointer <- setdiff(grep('GHz', fileLines) , grep('mean', fileLines))
    StokesEntryList <- fileLines[StokesPointer]
    StokesMatrix <- matrix(sapply(StokesEntryList, parseScan), nrow=9)
    #---- Expand scan to Stokes
    entryNum <- length(StokesPointer)
    srcList <- as.character(entryNum)
    ELList  <- as.numeric(entryNum)
    UTCList <- rep(UTCVec[1][[1]], entryNum)
    for(index in 1:entryNum){
        pointer <- which(scanPointer == max( scanPointer[scanPointer < StokesPointer[index]]))
        srcList[index] <- srcVec[pointer]
        ELList[index]  <- ELVec[pointer]
        UTCList[index] <- UTCVec[pointer][[1]]
    }
    return(data.frame(Src=srcList, Freq=StokesMatrix[1,], EL=ELList, I=StokesMatrix[2,], Q=StokesMatrix[4,], U=StokesMatrix[6,], V=StokesMatrix[8,], eI=StokesMatrix[3,], eQ=StokesMatrix[5,], eU=StokesMatrix[7,], eV=StokesMatrix[9,], Date=UTCList, File=rep(fileName, entryNum)))
}
#-------- Number of antennas
getAntNum <- function(Lines){
    D_pointer <- grep('D-term', Lines)
    if(length(D_pointer) == 0){ return(0) }
	return(length(grep('CM', Lines[1:D_pointer])) + length(grep('PM', Lines[1:D_pointer])) + length(grep('DA', Lines[1:D_pointer])) + length(grep('DV', Lines[1:D_pointer])))
}
#-------- Receiver Band
getBand <- function(fileList){
    return(as.integer(mclapply(fileList, function(fileName){
            bandPointer <- as.integer(regexpr("RB_[0-10]", fileName)[1])
            return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
        },mc.cores=numCore)))
}
#-------- Count record number
countRec <- function(fileName){
    fileLines <- removeBlank(readLines(fileName, skipNul=TRUE))
    return(length(grep('GHz', fileLines)) - length(grep('mean', fileLines)))
}
#-------- remove blank lines
removeBlank <- function(Lines){ return( Lines[which(nchar(Lines) > 1)]) }
#-------- Start program
#Arguments <- commandArgs(trailingOnly = T)
Arguments <- 'fileList'
fileList <- parseArg(Arguments)
#-------- Filter by number of used antennas
fileList <- mclapply(fileList, function(fileName){
    fileLines <- removeBlank(readLines(fileName, skipNul=TRUE))
    if(getAntNum(fileLines) < minAntNum){fileName <- 'FlaggedByAntNum' }
    return(fileName)}, mc.cores=numCore)
fileList <- fileList[fileList != 'FlaggedByAntNum']
#-------- Count number of records
bandID <- getBand(fileList)
recordNum <- sum(as.integer(mclapply(fileList, countRec, mc.cores=numCore)))
#-------- Generage FLDF
DFList <- mclapply(fileList, readStokesSection, mc.cores=numCore)
#FLDF <- do.call("rbind", DFList)
FLDF <- bind_rows(DFList)
#-------- Filter FLDF
FLDF <- na.omit(FLDF)
FLDF <- FLDF[FLDF$I > 2.0* FLDF$eI, ]                       # too large error
FLDF <- FLDF[FLDF$I^2  > FLDF$Q^2 + FLDF$U^2 +  FLDF$V^2,]  # polarization degree
FLDF$Src <- as.character(mclapply(as.character(FLDF$Src), sourceMatch,mc.cores=numCore))
FLDF$eI <- sqrt(FLDF$eI^2 + (sysIerr*FLDF$I)^2)
FLDF$eQ <- sqrt(FLDF$eQ^2 + (sysPerr*FLDF$Q)^2)
FLDF$eU <- sqrt(FLDF$eU^2 + (sysPerr*FLDF$U)^2)
FLDF$eV <- sqrt(FLDF$eV^2 + (sysPerr*FLDF$V)^2)
FLDF$Date <- as.POSIXct(FLDF$Date, tz="GMT")
FLDF <- FLDF[order(FLDF$Date),]
save(FLDF, file='Flux.Rdata')
largeErrorEBs <- unique( FLDF[((substr(FLDF$Src,1,1) == 'J') & (FLDF$eI > 10)),]$File)

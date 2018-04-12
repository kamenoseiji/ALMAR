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
    spw_pointer <- grep("^ SPW[0-99]", Lines)
    spwNum <- ceiling(length(spw_pointer) / numSource)
	StokesI <- StokesQ <- StokesU <- StokesV <- numeric(numSource)
	errI <- errQ <- errU <- errV <- numeric(numSource)
	I <- Q <- U <- V <- numeric(0)
	eI <- eQ <- eU <- eV <- numeric(0)
	srcList <- character(0); EL <- numeric(0)
	srcUTC <- as.Date(as.character(NULL))
	for(srcIndex in 1:numSource){
		srcList <- append(srcList, as.character(strsplit(Lines[pointer[srcIndex] - spwNum - 4], '[ |=]+')[[1]][3]) )
		srcUTC <- append(srcUTC, strptime(strsplit(Lines[pointer[srcIndex] - spwNum - 4], ' +')[[1]][6], "%Y/%m/%d/%H:%M:%S", tz="UTC"))
		EL <- append(EL, as.numeric(strsplit(Lines[pointer[srcIndex] - spwNum - 4], '[ |=]+')[[1]][5]))
		lineElements <- strsplit(Lines[pointer[srcIndex]], '[ |(|)|z]+')[[1]]
		FREQ <- as.numeric(lineElements[3])
		I <- append(I, as.numeric(lineElements[5])); eI <- append(eI, as.numeric(lineElements[6]))
		Q <- append(Q, as.numeric(lineElements[7])); eQ <- append(eQ, as.numeric(lineElements[8]))
		U <- append(U, as.numeric(lineElements[9])); eU <- append(eU, as.numeric(lineElements[10]))
		V <- append(V, as.numeric(lineElements[11])); eV <- append(eV, as.numeric(lineElements[12]))
	}
	return(data.frame(Src=as.character(srcList), Freq=FREQ, EL=EL, I=I, Q=Q, U=U, V=V, eI=eI, eQ=eQ, eU=eU, eV=eV, Date=srcUTC))
}

#-------- remove blank lines
removeBlank <- function(Lines){
	lineLength <- nchar(Lines)
	index <- which(lineLength > 1)
	return(Lines[index])
}


#-------- Start program
Arguments <- commandArgs(trailingOnly = T)
fileList <- as.character(parseArg(Arguments)$filelist)
FMT <- c('Src', 'EL', 'I', 'Q', 'U', 'V', 'eI', 'eQ', 'eU', 'eV', 'EL')
FLDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(FLDF) <- FMT

flagNum <- list()
for(fileName in fileList){
	cat(fileName); cat('\n')
    fileLines <- removeBlank(readLines(fileName))
	DF <- readStokesSection(fileLines)
	DF$File <- fileName
	FLDF <- rbind(FLDF, na.omit(DF))
}
save(FLDF, file='Flux.Rdata')

findSSO <- function(DF, fileName){
    df <- DF[DF$File == fileName,]
    srcList <- unique(df$Src)
    if(length(grep('^J', srcList, invert=TRUE)) > 0){
        #text_sd <- strsplit(fileName, '-')[[1]][1]
        #write(text_sd, file='UIDList', append=TRUE)
        write(fileName, file='UIDList', append=TRUE)
    }
}

parseArg <- function( args ){
    argNum <- length(args)
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-S"){ startDate <- substring(args[index], 3) }
        if(substr(args[index], 1,2) == "-E"){ endDate <- substring(args[index], 3)}
    }
    return(list(startDate = startDate, endDate = endDate))
}

#-------- Start Program
#argList <- parseArg(commandArgs(trailingOnly = T))
if(file.exists('UIDList')){ file.remove('UIDList') }
if(file.exists('Flux.Rdata')){
    load('Flux.Rdata')
} else {
    load(url("https://www.alma.cl/~skameno/AMAPOLA/Flux.Rdata"))
}
#DF <- FLDF[((as.Date(FLDF$Date) >= as.Date(argList$startDate)) & (as.Date(FLDF$Date) <= as.Date(argList$endDate))),]
fileList <- unique(FLDF$File)
for(fileName in fileList){ findSSO(FLDF, fileName) }

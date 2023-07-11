Arguments <- commandArgs(trailingOnly = T)
#-------- Parse arguments
parseArg <- function( args ){
	srcNum <- argNum <- length(args)
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-D"){ refDate <- as.Date(substring(args[index], 3));    srcNum <- srcNum - 1 }
	}
	srcList = args[(argNum - srcNum + 1):argNum]
	return(list(refDate = refDate, srcList = srcList[grep('^J[0-9]',srcList )]))
}
argList <- parseArg(Arguments)
refDate <- argList$refDate
srcList <- argList$srcList
RADDEG <- 180.0/pi
#-------- Load Flux.Rdata from web
#load(url("https://www.alma.cl/~skameno/Grid/Stokes/Flux.Rdata")) 
load("Flux.Rdata") 
FLDF$timeDiff <- as.numeric(difftime(FLDF$Date, refDate, units='days'))
FLDF <- FLDF[abs(FLDF$timeDiff) < 1.01,]
#-------- For each source
for(sourceName in srcList){
	srcDF <- FLDF[FLDF$Src == sourceName,]
	for(index in 1:nrow(srcDF)){
		text_sd <- sprintf('%s : %5.1f GHz : %5.2f (%.3f) %.1f (%.1f) %6.1f (%4.1f) %s\n', srcDF[index,]$Src, srcDF[index,]$Freq, srcDF[index,]$I, srcDF[index,]$eI, srcDF[index,]$P, srcDF[index,]$eP_upper, RADDEG* srcDF[index,]$EVPA, RADDEG* srcDF[index,]$eEVPA, srcDF[index,]$File) 
		cat(text_sd)
	}
}

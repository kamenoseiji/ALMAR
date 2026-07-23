library(RCurl)
library(RColorBrewer)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/ALMAR/refs/heads/master/StatStokes.R", ssl.verifypeer = FALSE)))
#source('../StatStokes.R')
substructureDF <- read.table('../Substructure.txt', header=TRUE)   # Read substructure contamination table
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq<-c(43.0,  75.0, 97.5, 132.0,183.0, 233.0, 343.5,400.0, 650.0, 800.0)
Pthresh7<-c(0.058, 0.060, 0.069, 0.058, 0.086, 0.077, 0.094, 0.153, 0.442, 0.765) # for 7m array, 5-sigma thresholds for polarized flux
Pthresh12 <- 0.5* Pthresh7     # for 12m array 
mthresh <- 0.03                                                                   # polarization degree threshold
SECPERDAY <- 86400
hourPerRad <- 12/pi		# radian to hour angle conversion
RADDEG <- 180.0/pi
DateRange <- 90    # 60-day window
Today <- Sys.time()
#-------- Filter calibrators for band
srcFreqCalibrator <- function(DF, band){
    sourceList <- sort(unique(DF$Src))
    numSrc <- length(sourceList)
    srcDF <- data.frame(Src=sourceList, I=numeric(numSrc), Q=numeric(numSrc), U=numeric(numSrc), V=numeric(numSrc), P=numeric(numSrc), EVPA=numeric(numSrc), eI=numeric(numSrc), eQ=numeric(numSrc), eU=numeric(numSrc), eV=numeric(numSrc), eP=numeric(numSrc), eEVPA=numeric(numSrc), subThresh=numeric(numSrc))
    for(src in sourceList){
        SDF <- DF[DF$Src == src,]
        SDF <- SDF[SDF$eI < 0.5*SDF$I,]
        if(nrow(SDF) < 5){ next }
        srcDF[srcDF$Src == src,][1:13] <- estimateIQUV(SDF, BandFreq[band], Today)
        if( src %in% substructureDF$Source ){
            strdf <- substructureDF[substructureDF$Source == src,]
            srcDF[srcDF$Src == src,]$subThresh  <- sqrt(strdf$StokesQ^2 + strdf$StokesU^2) * exp( strdf$SPIX* log(BandFreq[band]/100.0)) / mthresh
        }
    }
    srcDF$RA  <- pi* (60.0* as.numeric(substring(sourceList, 2, 3)) + as.numeric(substring(sourceList, 4, 5))) / 720.0
    srcDF$DEC <- pi* sign(as.numeric(substring(sourceList, 6, 10)))* (as.numeric(substring(sourceList, 7, 8)) + as.numeric(substring(sourceList, 9, 10))/60.0) / 180.0
    return( srcDF[((srcDF$P - srcDF$eP > Pthresh12[band]) & ((srcDF$P - srcDF$eP)/(srcDF$I + srcDF$eI) > 0.03) & (srcDF$P - srcDF$eP > srcDF$subThresh)),] )     # Filter by polarized flux and polarization degree
}
#-------- Load Flux.Rdata from web
#FluxDataURL <- "https://www.alma.cl/~skameno/AMAPOLA/"
#load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
load("Flux.Rdata")     # Data frame of FLDF
FLDF <- FLDF[as.Date(FLDF$Date) > as.Date(Today) - DateRange,]  # Data frame within DateRange
#-------- Filter quasars
FLDF <- FLDF[substr(FLDF$Src, 1, 1) == 'J',]    # only quasars
FLDF$P  <- sqrt(FLDF$Q^2 + FLDF$U^2)
sourceList <- sort(unique(FLDF$Src))
FLDF$medP <- FLDF$freqRange <- numeric(nrow(FLDF))
for(src in sourceList){ FLDF[FLDF$Src == src,]$medP <- median(FLDF[FLDF$Src == src,]$P)}
for(src in sourceList){ FLDF[FLDF$Src == src,]$freqRange <- diff(range(FLDF[FLDF$Src == src,]$Freq))}
FLDF <- FLDF[FLDF$medP > 0.03,]
FLDF <- FLDF[FLDF$freqRange > 100,]
sourceList <- sort(unique(FLDF$Src))
numSrc <- length(sourceList)
#-------- Loop in frequency band
for(band in seq(1, 7)){
    #-------- Today's IQUV
    srcDF <- na.omit(srcFreqCalibrator(FLDF, band))  # source properties (I, Q, U, V, P, EVPA) at the band
    save(file=sprintf('PolCalBand%d.Rdata', band), srcDF)
}

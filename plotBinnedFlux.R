#-------- plotBinnedFlux.R : R script to plot (Q, U) for specified source, band, and observation period --------
# Usage : Rscript plotQU.R [options]
#   Options:
#        -F[Start Date]  : start date to plot (e.g. -F2016-04-01)
#        -L[End Date]    : end date to plot (e.g. -E2016-09-01)
#        -S[Source Name] : source name (e.g. -SJ0854+2006)
#        -B[band]        : ALMA receiver band (e.g. -B3,6,7)
#        -b[bin]         : Binning days (e.g. -b10)
#        -w[width]       : Binning width in days (e.g. -w20)
#        -L              : Force loading Flux.Rdata from the website
# Requirement:
#   install.packages('plotrix') ,in advance, unless installed
# Caution : No space between the option flag and parameter, i.e. '-F2016-04-01' is OK but '-F 2016-04-01' is invalid.
#
library(xtable)
parseArg <- function( args ){
    argList <- list('2012-01-01', '2026-01-01', 'J1331+3030', '3,7', '10', '15', FALSE)
    names(argList) <- c('ST', 'ET', 'Source', 'Band', 'bin', 'width', 'Load')
    if( !file.exists("Flux.Rdata") ){ argList$load <- TRUE }    # Force downloading Flux.Rdata 
    argNum <- length(args)
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-F"){ argList$ST <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-E"){ argList$ET <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-S"){ argList$Source <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-B"){ argList$Band <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-b"){ argList$bin <- as.numeric(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-w"){ argList$width <- as.numeric(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-L"){ argList$Load <- TRUE}
    }
    return(argList)
}
bindSPW <- function(df){
    centerfreq <- median(df$Freq)
    DF <- data.frame(Date = sort(unique(df$Date)), Freq=centerfreq, I=0, eI=-1, Q=0, eQ=-1, U=0, eU=-1, V=0, eV=-1)
    for(date in DF$Date){
        index <- which(DF$Date == date)
        dayDF <- df[df$Date == date,]
        dayDF$Freq <- dayDF$Freq - centerfreq
        fitI <- lm( formula=I ~ Freq, data=dayDF, weight=1/eI^2); DF[index,]$I <- as.numeric(coef(summary(fitI))[[1]]); DF[index,]$eI <- as.numeric(coef(summary(fitI))[[3]])
        fitQ <- lm( formula=Q ~ Freq, data=dayDF, weight=1/eQ^2); DF[index,]$Q <- as.numeric(coef(summary(fitQ))[[1]]); DF[index,]$eQ <- as.numeric(coef(summary(fitQ))[[3]])
        fitU <- lm( formula=U ~ Freq, data=dayDF, weight=1/eU^2); DF[index,]$U <- as.numeric(coef(summary(fitU))[[1]]); DF[index,]$eU <- as.numeric(coef(summary(fitU))[[3]])
        fitV <- lm( formula=V ~ Freq, data=dayDF, weight=1/eV^2); DF[index,]$V <- as.numeric(coef(summary(fitV))[[1]]); DF[index,]$eV <- as.numeric(coef(summary(fitV))[[3]])
    }
    return( DF )
}
bindDate <- function(df, bin, width){
    DF <- data.frame( Date = seq(min(df$Date), max(df$Date), by=width*86400))
    DF$eV <- DF$V <- DF$eU <- DF$U <- DF$eQ <- DF$Q <- DF$eI <- DF$I <- numeric(nrow(DF))
    #for(refDate in as.Date(DF$Date)){
    for(index in 1:nrow(DF)){
        subDF <- df[ abs(as.Date(df$Date) - as.Date(DF[index,]$Date)) < width,]
        if(nrow(subDF) < 1){
            DF[index,]$I <- NA
        } else {
            subDF$offset <- as.numeric( as.Date(subDF$Date) - as.Date(DF[index,]$Date) )
            #repI <- weightedMean( data.frame(offset
            #repI <- weightedMean(subDF$I, subDF$eI, subDF$weight); DF[index, ]$I <- repI[1]; DF[index, ]$eI <- repI[2]
            #repQ <- weightedMean(subDF$Q, subDF$eQ, subDF$weight); DF[index, ]$Q <- repQ[1]; DF[index, ]$eQ <- repQ[2]
            #repU <- weightedMean(subDF$U, subDF$eU, subDF$weight); DF[index, ]$U <- repU[1]; DF[index, ]$eU <- repU[2]
            #repV <- weightedMean(subDF$V, subDF$eV, subDF$weight); DF[index, ]$V <- repV[1]; DF[index, ]$eV <- repV[2]
        }
    }
    return( na.omit(DF) )
}

weightedMean <- function(value, err, offset){
    
    wvalue <- weight / err^2
    sumWeight <- sum(wvalue)
    represent <- value %*% wvalue / sumWeight
    repErr    <- sqrt( err^2 %*% wvalue^2) / sumWeight
    return( c(represent, repErr) )
}
plotIP <- function(df, color='black'){
    numPoints <- nrow(df)
    numGrad <- 16
    alpha_grad <- seq(0.0, 0.1, length.out=numGrad)
    grid(col='gray', lty='dotted')
    abline(v=0)
    abline(h=0)
    arrows(df$Q, df$U-df$eU, df$Q, df$U+df$eU, length=0, col=color, lwd=0.25)
    arrows(df$Q-df$eQ, df$U, df$Q+df$eQ, df$U, length=0, col=color, lwd=0.25)
    points(df$Q, df$U,pch=20, col=color, cex=4.0e-3/sqrt(df$eQ^2 + df$eU^2))
    arrows(df[1:(numPoints-1),]$Q, df[1:(numPoints-1),]$U, df[2:numPoints,]$Q, df[2:numPoints,]$U, lwd=2, length=0.1, col=color)
    text(df$Q, df$U, substr(as.character(df$Date), 6, 10), col=color, pos=2, cex=0.5)
}
BandRepFreq <- c(40, 86.0, 97.5, 154.9, 183.0, 233.0, 343.4)
BandColors <- c('firebrick4', 'deeppink4', 'orangered4', 'darkgreen', 'turquoise4', 'royalblue4', 'purple4')
#-------- Start
#Arguments <- commandArgs(trailingOnly = TRUE)
Arguments <- strsplit('-F2024-01-01 -E2026-03-03 -SJ1256-0547 -B3,6,7 -b10 -w10', ' ')[[1]]    # for debugging
argList <- parseArg(Arguments)
setwd('./')
#-------- Load Flux.Rdata from web
if( argList$Load ){
    cat('--- Loading Flux.Rdata from the web\n')
    FluxDataURL <- "https://www.alma.cl/~skameno/AMAPOLA/"
    load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
    save(FLDF, file='Flux.Rdata')
} else {
    load('Flux.Rdata')
}
#-------- Filter by the source, date, and band
bandList <- as.numeric(strsplit(argList$Band, ',')[[1]])
for(band in bandList){
    df <- FLDF[((FLDF$Src == argList$Source) & (abs(FLDF$Freq - BandRepFreq[band]) < 12.0) & (as.Date(FLDF$Date) >= as.Date(argList$ST)) & (as.Date(FLDF$Date) <= as.Date(argList$ET))),]
    if(nrow(df) < 2){
        bandList <- bandList[!bandList %in% band]
        next
    }
    #df <- bindSPW(df)
    #df$Src <- argList$Source
    #DF <- bindDate(df, argList$bin, argList$width)
    #DF$P <- sqrt(DF$Q^2 + DF$U^2); DF$eP <- sqrt(DF$eQ^2 + DF$eU^2); DF$EVPA <- 90*atan2(DF$U, DF$Q)/pi; DF$eEVPA <- 90* sqrt(DF$Q^2 * DF$eU^2 + DF$U^2 * DF$eQ^2) / (DF$P)^2/pi
    #assign(sprintf('B%d', band), DF)
    #for(index in 1:nrow(DF)){
    #    text_sd <- sprintf('%s : %5.2f (%4.2f)  %5.2f (%4.2f) %5.1f (%3.1f)\n', as.Date(DF[index,]$Date), DF[index,]$I, DF[index,]$eI,  DF[index,]$P, DF[index,]$eP,  DF[index,]$EVPA, DF[index,]$eEVPA)
    #    cat(text_sd)
    #}
}
#if(0){
##-------- Plot (Q, U)
#dfList <- sprintf('B%d', bandList)
#labels <- dfList
#pdf(sprintf('%s_QU_%s_%s.pdf', argList$Source, argList$ST, argList$ET))
#plot(df$Q, df$U, type='n', xlim=c(-pMax,pMax), ylim=c(-pMax,pMax), xlab='Stokes Q [Jy]', ylab='Stokes U [Jy]', main=sprintf('%s %s to %s', argList$Source, argList$ST, argList$ET))
#for(band in bandList){
#    index <- which(bandList == band)
#    df <- get(dfList[index])
#    plotQU(df, BandColors[band])
#    labels[index] <- sprintf('%.1f GHz', df$Freq[1])
#}
#legend('topleft', legend=labels, col=BandColors[bandList], pch=20)
#dev.off()
#}

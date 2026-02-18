#-------- plotQU.R : R script to plot (Q, U) for specified source, band, and observation period --------
# Usage : Rscript plotQU.R [options]
#   Options:
#        -F[Start Date]  : start date to plot (e.g. -F2016-04-01)
#        -L[End Date]    : end date to plot (e.g. -E2016-09-01)
#        -S[Source Name] : source name (e.g. -SJ0854+2006)
#        -B[band]        : ALMA receiver band (e.g. -B3,6,7)
#        -L              : Force loading Flux.Rdata from the website
# Requirement:
#   install.packages('plotrix') ,in advance, unless installed
# Caution : No space between the option flag and parameter, i.e. '-F2016-04-01' is OK but '-F 2016-04-01' is invalid.
#
library(plotrix)
parseArg <- function( args ){
    argList <- list('2012-01-01', '2026-01-01', 'J1331+3030', '3,6,7', FALSE)
    names(argList) <- c('ST', 'ET', 'Source', 'Band', 'Load')
    if( !file.exists("Flux.Rdata") ){ argList$load <- TRUE }    # Force downloading Flux.Rdata 
    argNum <- length(args)
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-F"){ argList$ST <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-E"){ argList$ET <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-S"){ argList$Source <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-B"){ argList$Band <- as.character(substring(args[index], 3)) }
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
plotQU <- function(df, color='black'){
    numPoints <- nrow(df)
    numGrad <- 16
    alpha_grad <- seq(0.0, 0.1, length.out=numGrad)
    grid(col='gray', lty='dotted')
    abline(v=0)
    abline(h=0)
    arrows(df$Q, df$U-df$eU, df$Q, df$U+df$eU, length=0, col=color, lwd=0.25)
    arrows(df$Q-df$eQ, df$U, df$Q+df$eQ, df$U, length=0, col=color, lwd=0.25)
    for(index in 1:numGrad){
        current_col <- adjustcolor(color, alpha.f=alpha_grad[index])
        draw.ellipse(df$Q, df$U, a=(numGrad - index + 1)*df$eQ/numGrad, b=(numGrad - index + 1)*df$eU/numGrad, col=current_col, border=NA)
    }
    points(df$Q, df$U,pch=20, col=color, cex=4.0e-3/sqrt(df$eQ^2 + df$eU^2))
    arrows(df[1:(numPoints-1),]$Q, df[1:(numPoints-1),]$U, df[2:numPoints,]$Q, df[2:numPoints,]$U, lwd=2, length=0.1, col=color)
    text(df$Q, df$U, substr(as.character(df$Date), 6, 10), col=color, pos=2, cex=0.5)
}
BandRepFreq <- c(40, 86.0, 97.5, 154.9, 183.0, 233.0, 343.4)
BandColors <- c('firebrick4', 'deeppink4', 'orangered4', 'darkgreen', 'turquoise4', 'royalblue4', 'purple4')
#-------- Start
Arguments <- commandArgs(trailingOnly = TRUE)
#Arguments <- strsplit('-F2016-04-01 -E2016-09-01 -SJ0854+2006 -B3,7', ' ')[[1]]    # for debugging
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
pMax <- 0.0     # maximum polarized flux to plot
for(band in bandList){
    df <- FLDF[((FLDF$Src == argList$Source) & (abs(FLDF$Freq - BandRepFreq[band]) < 12.0) & (as.Date(FLDF$Date) >= as.Date(argList$ST)) & (as.Date(FLDF$Date) <= as.Date(argList$ET))),]
    if(nrow(df) < 2){
        bandList <- bandList[!bandList %in% band]
    } else {
        df <- bindSPW(df)
        df$Src <- argList$Source
        pMax <- max(c(pMax, max(sqrt(df$Q^2 + df$U^2))))
        assign(sprintf('B%d', band), df)
    }
}
#-------- Plot (Q, U)
dfList <- sprintf('B%d', bandList)
labels <- dfList
pdf(sprintf('%s_QU_%s_%s.pdf', argList$Source, argList$ST, argList$ET))
plot(df$Q, df$U, type='n', xlim=c(-pMax,pMax), ylim=c(-pMax,pMax), xlab='Stokes Q [Jy]', ylab='Stokes U [Jy]', main=sprintf('%s %s to %s', argList$Source, argList$ST, argList$ET))
for(band in bandList){
    index <- which(bandList == band)
    df <- get(dfList[index])
    plotQU(df, BandColors[band])
    labels[index] <- sprintf('%.1f GHz', df$Freq[1])
}
legend('topleft', legend=labels, col=BandColors[bandList], pch=20)
dev.off()

library(RCurl)
FuncList <- c('date')
for(index in 1:length(FuncList)){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }
setwd('.')
MJDepoch <- doy2mjd(1970,1)
#-------- Execution
StokesDF <- read.table('StokesComp.log', header=T)
StokesDF$Date <- as.POSIXct(StokesDF$mjdSec - MJDepoch*86400, tz="GMT", origin="1970-01-01")
StokesDF$relTime <- StokesDF$mjdSec - median(StokesDF$mjdSec)
StokesDF$secZ <- 1.0 / sin(StokesDF$EL)
StokesDF$Src <- as.character(StokesDF$Src)
StokesDF$relErr <- (StokesDF$aI - StokesDF$I)/StokesDF$I
#-------- Band Separation
B3 <- StokesDF[StokesDF$Band == 3,]
B6 <- StokesDF[StokesDF$Band == 6,]
B7 <- StokesDF[StokesDF$Band == 7,]
colores = c('red', 'blue')
labels = c('Standard', 'a priori')

for(band in c(3,6,7)){
    AFA <- AFS <- numeric(0)    # Astonishing factor
    BandDF <- na.omit(StokesDF[StokesDF$Band == band,])
    png(sprintf('histFluxErrB%d.png', band))
    relErr <- (BandDF$aI - BandDF$I)/BandDF$I
    k = seq(-0.3, 0.33, by=0.01)
    hist(BandDF$relErr, col="#0000ff40", border=F, xlab='(a priori - standard)/standard', breaks=k, main=sprintf('Band-%d', band))
    abline(v=-0.05, lty="dotted"); abline(v=0.05, lty="dotted")
    text_sd <- sprintf("Band-%d : mean=%.3f, sd=%.3f, quantile[0.05 - 0.95] = [%.3f, %.3f]", band, mean(BandDF$relErr), sd(BandDF$relErr), quantile(BandDF$relErr, 0.05), quantile(BandDF$relErr, 0.95))
    print(text_sd)
    dev.off()
    sourceList <- sort(unique(BandDF$Src))
    src_index <- 0
    pdf(sprintf('Flux-Band%d.pdf', band))
    for(source in sourceList){
        SrcDF <- BandDF[BandDF$Src == source,]
        if(nrow(SrcDF) < 10){ next}
        src_index <- src_index + 1
        plot(SrcDF$Date, SrcDF$I, type='n', xlab='Date', ylab='Flux Density [Jy]', ylim=c(1.2*min(SrcDF$I)-0.2*max(SrcDF$I), -0.2*min(SrcDF$I)+1.2*max(SrcDF$I)), main=sprintf('%s Band %d', source, band))
        arrows(SrcDF$Date, SrcDF$I+SrcDF$eI, SrcDF$Date, SrcDF$I-SrcDF$eI, length=0, col='red')
        lines(SrcDF$Date, SrcDF$I, col='red'); points(SrcDF$Date, SrcDF$I, col='red', pch=20)
        arrows(SrcDF$Date, SrcDF$aI+SrcDF$eaI, SrcDF$Date, SrcDF$aI-SrcDF$eaI, length=0, col='blue')
        lines(SrcDF$Date, SrcDF$aI, col='blue'); points(SrcDF$Date, SrcDF$aI, col='blue', pch=20)
        AFS[src_index] <- sd(diff(SrcDF$I))/mean(SrcDF$I); AFA[src_index] <- sd(diff(SrcDF$aI))/mean(SrcDF$I)
        labels <- c( sprintf('AF=%.3f : Standard', AFS[src_index]), sprintf('AF=%.3f : a priori', AFA[src_index]))
        legend("topright", legend=labels, col=colores, pch=20)
    }
    dev.off()
    #if(0){
    k = seq(0, max(c(AFS, AFA))+0.02, by=0.01)
    labels <- c(sprintf('AF(mean)=%.3f : Standard', mean(AFS)), sprintf('AF(mean)=%.3f : a priori', mean(AFA)))
    text_sd <- sprintf('Band %d : AFS = %.3f  AFA = %.3f', band, mean(AFS), mean(AFA)); print(text_sd)
    colores <- c("#0000ff40", "#ff000040")
    png(sprintf('histAF-B%d.png', band))
    hist(AFS, col="#0000ff40", border=F, xlab='Astonishing Factor', breaks=k, main=sprintf('Band-%d Astonishing Factor', band))
    hist(AFA, col="#ff000040", border=F, breaks=k, add=T)
    legend("topright", legend=labels, col=colores, pch=15,cex=1)
    dev.off()
    #}
}
library(RColorBrewer)
library(plotly, warn.conflicts=FALSE)
library(parallel)   # multicore parallelization
#library(doParallel) # multicore parallelization
library(VGAM)       # for Rice distribution
Sys.setenv(TZ="UTC")
ALMA_lat <- -23.029 * pi / 180
#sysIerr <- 0.005       # temporal Stokes I systematic error
sysPerr <- 0.003       # temporal polarization systematic error
minAntNum <- 5		   # Minimum number of antennas
standardFreq <- list(40.0, 80.0, c(91.5,97.5,103.5), 154.9, 183.0, 233.0, 343.4, 410.2)
numCore = detectCores()
#-------- Functions
srcDfFilter  <- function(source){ return(FLDF[FLDF$Src == source,])}
scanDfFilter <- function(scan){ return(FLDF[FLDF$Date == scan,])}
#-------- Get band name
getBand <- function(fileName){
    bandPointer <- regexpr("RB_[0-10]", fileName)
    return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
}
#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
predStokes <- function(df){
    #bandID   <- getBand(df$File[1])
    bandID   <- df$Band[1]
    fitI <- lm(formula=I ~ Freq, data=df, weight=1.0/eI^2)
    fitQ <- lm(formula=Q ~ Freq, data=df, weight=1.0/eQ^2)
    fitU <- lm(formula=U ~ Freq, data=df, weight=1.0/eU^2)
    fitV <- lm(formula=V ~ Freq, data=df, weight=1.0/eV^2)
    newDF <- data.frame(Src=df$Src[1], Freq = standardFreq[[bandID]], EL=df$EL[1])
    pred <- as.numeric(predict(fitI, newDF, interval='confidence', level=0.67)); newDF$I <- matrix(pred, ncol=3)[,1]; newDF$eI <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitQ, newDF, interval='confidence', level=0.67)); newDF$Q <- matrix(pred, ncol=3)[,1]; newDF$eQ <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitU, newDF, interval='confidence', level=0.67)); newDF$U <- matrix(pred, ncol=3)[,1]; newDF$eU <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitV, newDF, interval='confidence', level=0.67)); newDF$V <- matrix(pred, ncol=3)[,1]; newDF$eV <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    newDF$Date <- df$Date[1]
    newDF$File <- df$File[1]
    return(newDF)
}
#-------- Text Formatting
textResult <- function(entry){
    text_sd <- sprintf()
    return(text_sd)
}
getBand <- function(fileList){
    return(as.integer(mclapply(fileList, function(fileName){
            bandPointer <- as.integer(regexpr("RB_[0-10]", fileName)[1])
            return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
        },mc.cores=numCore)))
}
#-------- Starging Process
load('Flux.Rdata')
recNum <- nrow(FLDF)
FLDF <- FLDF[order(FLDF$Date, FLDF$Src, FLDF$Band, FLDF$Freq),]
sliceDF <- function(slice){ return(FLDF[slice[1]:slice[2],])}
diffRange1 <- 1:(recNum-1)
diffRange2 <- 2:recNum
recBorder <- c(1, which( (FLDF$Src[diffRange1] != FLDF$Src[diffRange2]) | diff(FLDF$Date) > 0 | diff(FLDF$Band) > 0 )+1)
recLength <- diff(c(recBorder, recNum+1) )
recMat <- matrix(c(recBorder,   c(tail(recBorder, n=length(recBorder)-1)-1, recNum) ), ncol=2) 
tempDF <- apply(recMat, 1, sliceDF)
result <- mclapply(tempDF, predStokes, mc.cores=numCore)
textDF <- do.call("rbind", result)
textDF <- na.omit(textDF)
textDF$P <- sqrt(textDF$Q^2 + textDF$U^2)
sigmaSQ <- sqrt(textDF$eQ * textDF$eU + (textDF$I* sysPerr)^2)
textDF$eP_lower <- qrice(0.15, sigmaSQ, textDF$P)
textDF[textDF$P < sigmaSQ,]$eP_lower <- 0.0
textDF$eP_upper <- qrice(0.85, sigmaSQ, textDF$P)
textDF$EVPA <- 0.5* atan2(textDF$U, textDF$Q)
textDF$eEVPA <- 0.5* sqrt(textDF$Q^2 * textDF$eU^2 + textDF$U^2 * textDF$eQ^2) / (textDF$P)^2
#---- Output to text data
textDF$Src <- sprintf('%10s ', textDF$Src)
write.table(format(textDF, digits=4), 'amapola.txt', sep='\t', quote=F, col.names=T, row.names=F)
#-------- fiter Band-3 upper and lower sidebands
textDF <- textDF[((textDF$Freq != standardFreq[[3]][1]) & (textDF$Freq != standardFreq[[3]][3])),]
textDF$Band <- getBand(textDF$File)
#-------- Plot time-seried flux densities
textDF$Src <- trimws(textDF$Src)
sourceList <- unique(textDF$Src)
sourceList <- sourceList[grep('^J[0-9]',sourceList)]
sourceList <- sourceList[order(sourceList)]
pdf('AMAPOLA.pdf', width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
labels <- c("B1",        "B3",        "B4",        "B6",        "B7",        "> B7")
colors <- c("#FF80003F", "#E020003F", "#8080003F", "#00C0803F", "#0000FF3F", "#0000003F")
#threshU <-c(41.0,        98.0,        200.0,       280.0,       380.0,       700.0)
#threshL <-c(39.0,        97.0,        120.0,       200.0,       280.0,       380.0)
plotXrange <- range(FLDF$Date)
bandList <- unique(textDF$Band)
bandIndex <- c(1,0,2,3,0,4,5,6,0)
for(sourceName in sourceList){
	DF <- textDF[textDF$Src == sourceName,]
    #-------- Coloring by frequency
    DF$color_vector <- rep("#000000FF", nrow(DF))
    for(band in bandList){
        if(nrow(DF[DF$Band == band,]) == 0){ next }
        DF[DF$Band == band,]$color_vector <- colors[bandIndex[band]]
    }
    DF <- DF[DF$color_vector != "#000000FF",]
	#-------- Plot Stokes I
	plot(DF$Date, DF$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='Total Flux Density', xlim=as.POSIXct(plotXrange), ylim=c(0, 1.2*max(DF$I)))
	arrows(as.numeric(DF$Date), DF$I - DF$eI, as.numeric(DF$Date), DF$I + DF$eI, angle=90, length=0.0, col=DF$color_vector)
	points(DF$Date, DF$I, pch=20, col=DF$color_vector)
    legend("bottomleft", legend=labels, pch=20, col=colors)
    axis.Date(1, at=seq(min(DF$Date), max(DF$Date), "months"), format="%Y-%m", cex=0.5)
	
	#-------- Plot polarized flux
	plot(DF$Date, DF$P, type='n', xlab='Date', ylab='Polarized Flux [Jy]', main='Polarized Flux Density', xlim=as.POSIXct(plotXrange), ylim=c(0, 1.2* max(DF$P)))
	arrows(as.numeric(DF$Date), DF$eP_lower, as.numeric(DF$Date), DF$eP_upper, angle=90, length=0.0, col=DF$color_vector)
	points(DF$Date, DF$P, pch=20, col=DF$color_vector)
    axis.Date(1, at=seq(min(DF$Date), max(DF$Date), "months"), format="%Y-%m", cex=0.5)
	
	#-------- Plot EVPA
	plot(DF$Date, DF$EVPA, type='n', xlab='Date', ylab='EVPA [deg]', main='Polarization Angle', xlim=as.POSIXct(plotXrange), ylim=c(-90,90))
	abline(h=90, lwd=0.1); abline(h=-90, lwd=0.1)
	arrows(as.numeric(DF$Date), 180*(DF$EVPA - DF$eEVPA)/pi, as.numeric(DF$Date), 180*(DF$EVPA + DF$eEVPA)/pi, angle=90, length=0.0, col=DF$color_vector)
	points(DF$Date, DF$EVPA*180/pi, pch=20, col=DF$color_vector)
    axis.Date(1, at=seq(min(DF$Date), max(DF$Date), "months"), format="%Y-%m", cex=0.5)
	mtext(side = 3, line=1, outer=T, text = sourceName, cex=2)

    plot_I <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~I, type="scatter", mode="markers", color=paste("I", freqLabel[1]), colors=bandColor, error_y = list(array=~eI, thickness=1, width=0))
    plot_P <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~P, type="scatter", mode="markers", color=paste("P",freqLabel[1]), colors=bandColor, error_y = list(symmetric=FALSE, array=~errU, arrayminus=~errL, thickness=1, width=0))
    plot_A <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~EVPA*180/pi, type="scatter", mode="markers", color=paste("EVPA",freqLabel[1]), colors=bandColor, error_y = list(array=~eEVPA*180/pi, thickness=1, width=0))
    allPlot <- subplot(plot_I, plot_P, plot_A, nrows=3, shareX=T, titleY=T)
    htmlFile <- sprintf("%s.flux.html", sourceList[src_index])
    htmlwidgets::saveWidget(allPlot, htmlFile, selfcontained=FALSE )
    rm(plot_I); rm(plot_P); rm(plot_A); rm(allPlot); rm(htmlFile)
}
par(par.old)
dev.off()
text_sd <- sprintf('cp -r %s.flux_files common.flux_files', sourceList[1])
system(text_sd)
system('rm -rf J*.flux_files')
for(sourceName in sourceList){
    text_sd <- sprintf('ln -s common.flux_files %s.flux_files', sourceName)
    system(text_sd)
}
if(0){
#-------- Time-series HTML
bandColor <- brewer.pal(numFreq, "Dark2")
for(sourceName in sourceList){
	DF <- textDF[textDF$Src == sourceName,]
    plot_I <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~I, type="scatter", mode="markers", color=paste("I", freqLabel[1]), colors=bandColor, error_y = list(array=~eI, thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_I <- add_trace(plot_I, data=DF[DF$Band == bandList[freq_index],], color=paste("I", freqLabel[freq_index]))}
    }
    plot_I <- layout(plot_I, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Stokes I [Jy]',rangemode='tozero', hoverformat='.3f'), title=sourceList[src_index])

    plot_P <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~P, type="scatter", mode="markers", color=paste("P",freqLabel[1]), colors=bandColor, error_y = list(symmetric=FALSE, array=~errU, arrayminus=~errL, thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_P <- add_trace(plot_P, data=DF[DF$Band == bandList[freq_index],], color=paste("P",freqLabel[freq_index]))}
    }
    plot_P <- layout(plot_P, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Polarized Flux [Jy]',rangemode='tozero', hoverformat='.3f'))

    plot_A <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~EVPA*180/pi, type="scatter", mode="markers", color=paste("EVPA",freqLabel[1]), colors=bandColor, error_y = list(array=~eEVPA*180/pi, thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_A <- add_trace(plot_A, data=DF[DF$Band == bandList[freq_index],], color=paste("EVPA",freqLabel[freq_index]))}
    }
    plot_A <- layout(plot_A, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='EVPA [deg]',range=c(-91,91), hoverformat='.1f'))

    allPlot <- subplot(plot_I, plot_P, plot_A, nrows=3, shareX=T, titleY=T)
    htmlFile <- sprintf("%s.flux.html", sourceList[src_index])
    htmlwidgets::saveWidget(allPlot, htmlFile, selfcontained=FALSE )
    rm(plot_I); rm(plot_P); rm(plot_A); rm(allPlot); rm(htmlFile)
}
text_sd <- sprintf('cp -r %s.flux_files common.flux_files', sourceList[1])
system(text_sd)
system('rm -rf J*.flux_files')
for(src_index in 1:numSrc){
    text_sd <- sprintf('ln -s common.flux_files %s.flux_files', sourceList[src_index])
    system(text_sd)
}
}

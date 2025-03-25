library(RColorBrewer)
library(xtable)
library(plotly, warn.conflicts=FALSE)
library(pandoc)
library(htmlwidgets)   # multicore parallelization
library(parallel)   # multicore parallelization
library(VGAM)       # for Rice distribution
FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
Sys.setenv(TZ="UTC")
ALMA_lat <- -23.029 * pi / 180
Today <- Sys.time()
recentTerm <- 60        # Statistics for recent 60 days
#sysIerr <- 0.005       # temporal Stokes I systematic error
sysPerr <- 0.003       # temporal polarization systematic error
minAntNum <- 5		   # Minimum number of antennas
standardFreq <- list(40.0, 80.0, c(91.5,97.5,103.5), 154.9, 183.0, 233.0, 343.4, 410.2)
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(40.0, 80.0, 97.5, 154.9, 183.0, 233.0, 343.4, 410.2, 650.0, 800.0)
numCore = detectCores()
#-------- Functions
srcDfFilter  <- function(source){ return(FLDF[FLDF$Src == source,])}
scanDfFilter <- function(scan){ return(FLDF[FLDF$Date == scan,])}
#-------- Hour angle to Az, El
ha2azel <-function(ha, phi=ALMA_lat, dec=0){
	# ha2azel : Calculate Azimuth, Elevation, and Parallactic Angle
	# ha  : Hour Angle
	# phi : Latitude of the site [rad]
	# dec : Declination of the source [rad]
	sin_el <- sin(phi)* sin(dec) + cos(phi)* cos(dec)* cos(ha)
	el <- asin(sin_el)
	az <- atan2(cos(dec)*sin(ha), sin(phi)*cos(dec)*cos(ha) - cos(phi)*sin(dec)) + pi
	pa <- atan2(sin(ha), tan(phi)*cos(dec) - sin(dec)* cos(ha))
	return(data.frame(ha=ha, az=az, el=el, pa=pa))
}
#-------- Get band name
getBand <- function(fileName){
    bandPointer <- regexpr("RB_[0-10]", fileName)
    return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
}
#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
predStokes <- function(df){
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
#-------- Plot LST
plotLST <- function(DF, band){
    pDF <- data.frame()
    DF <- DF[abs(DF$DEC - ALMA_lat) > 4.0*pi/180.0,]    # zenith avoidance of 4 degree
    sourceList <- as.character(DF$Src)
    sourceNum <- length(sourceList)
    LST <- seq(0.0, 23.95, 0.05)
    lstNum <- length(LST)
    for(src in sourceList){
        source_index <- which(DF$Src == src)
        HA <- LST*pi/12 - DF$RA[source_index]
        AZEL <- ha2azel( HA, ALMA_lat, DF$DEC[source_index] )
        AZEL$pa <- AZEL$pa + BandPA[band]
        CS <- cos(2.0* AZEL$pa); SN <- sin(2.0* AZEL$pa)
        freqIFact <- (BandFreq[band] / 100)^DF$spixI[source_index]
        freqPFact <- (BandFreq[band] / 100)^DF$spixP[source_index]
        predQ <- DF$Q100[source_index]* freqPFact; predU <- DF$U100[source_index] *freqPFact
        predP <- sqrt(predQ^2 + predU^2); predI <- DF$I100[source_index]* freqIFact
        if((predP < 0.05) | (predP < 0.03* predI))  next
        AZEL$XYcorr <- abs(predU*CS - predQ*SN)
        AZEL[AZEL$el < pi/7.5,]$XYcorr <- NA
        if( nrow(pDF) == 0 ){
            pDF <- data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*180/pi, XYcorr=AZEL$XYcorr)
        } else {
            pDF <- rbind(pDF, data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*180/pi, XYcorr=AZEL$XYcorr))
        }
    }
    return(pDF)
}
#-------- Starging Process
load('Flux.Rdata')
FLDF <- FLDF[FLDF$eI > 1e-5,]
FLDF$Band <- getBand(FLDF$File)
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
#-------- HTML table of source flux
recentDF <- textDF[as.Date(Today) - as.Date(textDF$Date) < recentTerm,]   # Recent 60 days
bandList <- sort(unique(recentDF$Band))
for(band in bandList){
    bandDF <- recentDF[recentDF$Band == band,]
    lastObsIndex <- order(bandDF$Date, decreasing=TRUE)[1]
    sourceList <- unique(bandDF$Src)
    sourceList <- sourceList[grep('^J[0-9]',sourceList)]
    srcDF <- mclapply(sourceList, function(source){
        DF <- bandDF[bandDF$Src == source, ]
        return(data.frame(Src=source, numObs = nrow(DF), I = median(DF$I), eI = sd(DF$I), Q = median(DF$Q), eQ = sd(DF$Q), U = median(DF$U), eU=sd(DF$U)))
    }, mc.cores=numCore)
    srcDF <- na.omit(do.call("rbind", srcDF))
    if(nrow(srcDF) == 0){next}
    srcDF$P  <- sqrt(srcDF$Q^2 + srcDF$U^2)
    srcDF$eP <- sqrt(srcDF$eQ^2 + srcDF$eU^2)
    srcDF$p  <- 100*srcDF$P/srcDF$I
    srcDF$EVPA  <- 90*atan2(srcDF$U, srcDF$Q)/pi
    srcDF$eEVPA <- 90*sqrt(srcDF$Q^2 * srcDF$eU^2 + srcDF$U^2 * srcDF$eQ^2) / (srcDF$P)^2 / pi
    srcDF <- srcDF[order(srcDF$P, decreasing=TRUE),]
    names(srcDF) <- c('Source', '#obs', 'I [Jy]', 'sd(I)', 'Q [Jy]', 'sd(Q)', 'U [Jy]', 'sd(U)', 'P [Jy]', 'sd(P)', '%pol', 'EVPA (deg)', 'sd(EVPA)')
    srcDF$Source <- paste('<a href="', srcDF$Source, '.flux.html" target="_new" >', srcDF$Source, ' </a>', sep='')
    #-------- HTML pol-table
    CaptionText <- paste("<p>", sprintf('Frequency %.1f GHz : %d-day median as of %s / ', BandFreq[band], recentTerm, as.character(Today)),sep='')

    CaptionText <- paste(CaptionText, '<a href="', FluxDataURL, bandDF$File[lastObsIndex], '" target="_new">', 'Last Observation on ', as.character(bandDF$Date[lastObsIndex], tz='UTC'), "</p>", sep='')
    htmlFile <- sprintf('Stokes%.0fGHz.html', BandFreq[band])
    html.head <- paste("<head>", '<link rel="stylesheet" type="text/css" href="https://www.alma.cl/~skameno/resources/amapola.css" />', "</head>", sep='\n')
    html.table <- paste(print(xtable(srcDF, digits=c(0,0,0,3,3,3,3,3,3,3,3,1,1,1)), include.rownames=F, type="html", sanitize.text.function=function(x){x}, htmlFile), collapse="\n")
    html.body <- paste("<body>", CaptionText, html.table, "</body>")
    write(paste(html.head, html.body, sep='\n'), htmlFile)
}
pdf('AMAPOLA.pdf', width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
labels <- c("B1",        "B3",        "B4",        "B6",        "B7",        "> B7")
freqLabels <- c('40 GHz', '97.5 GHz', '154.9 GHz', '233 GHz',   '343.4 GHz', '410.2 GHz')
colors <- c("#FF80003F", "#E020003F", "#8080003F", "#00C0803F", "#0000FF3F", "#0000003F")
plotXrange <- range(FLDF$Date)
bandList <- unique(textDF$Band)
bandIndex <- c(1,0,2,3,0,4,5,6,0)
sourceList <- unique(textDF$Src)
sourceList <- sourceList[grep('^J[0-9]',sourceList)]
sourceList <- sourceList[order(sourceList)]
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

}
par(par.old)
dev.off()
#-------- Time-series HTML
bandList <- sort(unique(FLDF$Band))
numFreq <- length(bandList)
bandColor <- brewer.pal(numFreq, "Dark2")
for(sourceName in sourceList){
	DF <- textDF[textDF$Src == sourceName,]
    plot_I <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~I, type="scatter", mode="markers", color=paste("I", freqLabels[1]), colors=bandColor, error_y = list(array=~eI, thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_I <- add_trace(plot_I, data=DF[DF$Band == bandList[freq_index],], color=paste("I", freqLabels[freq_index]))}
    }
    plot_I <- layout(plot_I, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Stokes I [Jy]',rangemode='tozero', hoverformat='.3f'), title=sourceName)

    plot_P <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~P, type="scatter", mode="markers", color=paste("P",freqLabels[1]), colors=bandColor, error_y = list(symmetric=FALSE, array=~(eP_upper - P), arrayminus=~(P-eP_lower), thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_P <- add_trace(plot_P, data=DF[DF$Band == bandList[freq_index],], color=paste("P",freqLabels[freq_index]))}
    }
    plot_P <- layout(plot_P, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Polarized Flux [Jy]',rangemode='tozero', hoverformat='.3f'))

    plot_A <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~EVPA*180/pi, type="scatter", mode="markers", color=paste("EVPA",freqLabels[1]), colors=bandColor, error_y = list(array=~eEVPA*180/pi, thickness=1, width=0))
    for(freq_index in 2:numFreq){
        if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_A <- add_trace(plot_A, data=DF[DF$Band == bandList[freq_index],], color=paste("EVPA",freqLabels[freq_index]))}
    }
    plot_A <- layout(plot_A, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='EVPA [deg]',range=c(-91,91), hoverformat='.1f'))

    allPlot <- subplot(plot_I, plot_P, plot_A, nrows=3, shareX=T, titleY=T)
    htmlFile <- sprintf("%s.flux.html", sourceName)
    htmlwidgets::saveWidget(allPlot, htmlFile, selfcontained=FALSE )
    rm(plot_I); rm(plot_P); rm(plot_A); rm(allPlot); rm(htmlFile)
}

text_sd <- sprintf('cp -r %s.flux_files common.flux_files', sourceList[1])
system(text_sd)
system('rm -rf J*.flux_files')
for(sourceName in sourceList){
    text_sd <- sprintf('ln -s common.flux_files %s.flux_files', sourceName)
    system(text_sd)
}
#-------- Source 60-day statistics
numSrc <- length(sourceList)
RAList <- (60.0* as.numeric(substring(sourceList, 2,3)) + as.numeric(substring(sourceList, 4,5))) / 720 * pi # RA in [rad]
DecList<- as.numeric(substring(sourceList, 6,8))
DecList<- DecList + sign(DecList)* as.numeric(substring(sourceList, 9,10))/60.0
DecList<- DecList / 180 * pi # DEC in [rad]
I100 <- Q100 <- U100 <- numeric(numSrc)
spixI <- spixP <- rep(NA, numSrc)
for(src_index in 1:numSrc){
    sourceName <- sourceList[src_index]
    srcDF <- textDF[((textDF$Src == sourceName) & (difftime(Today, textDF$Date, units="days") < 60)) , ]
    if(nrow(srcDF) < 6){ next }
    srcDF$deltaDay <- as.numeric(difftime(srcDF$Date, Today))
    bands <- unique(srcDF$Band)
    numBand <- length(bands)
    predI <- predQ <- predU <- eI <- eQ <- eU <- numObs <- freq <- numeric(numBand)
    if(numBand > 1){
        for(band_index in 1:numBand){
            fitDF <- srcDF[srcDF$Band == bands[band_index],]
            fitDF <- fitDF[abs(fitDF$I - median(fitDF$I)) < 100.0* fitDF$eI,]
            numObs[band_index] <- nrow(fitDF)
            BandTimeFit <- TRUE
            if( numObs[band_index] <= 2 ){ BandTimeFit <- FALSE }
            if(BandTimeFit){ if(sd(fitDF$deltaDay) < min(abs(fitDF$deltaDay)) ){ BandTimeFit <- FALSE } }
            if( BandTimeFit ){
                fit <- lm(I ~ deltaDay, weights=1/eI^2/abs(deltaDay + 1), data=fitDF)
                predI[band_index] <- summary(fit)$coefficients[1,'Estimate']
                eI[band_index] <- summary(fit)$coefficients[1,'Std. Error']
                fit <- lm(Q ~ deltaDay, weights=1/eQ^2/abs(deltaDay + 1), data=fitDF)
                predQ[band_index] <- summary(fit)$coefficients[1,'Estimate']
                eQ[band_index] <- summary(fit)$coefficients[1,'Std. Error']
                fit <- lm(U ~ deltaDay, weights=1/eU^2/abs(deltaDay + 1), data=fitDF)
                predU[band_index] <- summary(fit)$coefficients[1,'Estimate']
                eU[band_index] <- summary(fit)$coefficients[1,'Std. Error']
            } else {
                predI[band_index] <- mean(fitDF[fitDF$Band == bands[band_index],]$I)
                predQ[band_index] <- mean(fitDF[fitDF$Band == bands[band_index],]$Q)
                predU[band_index] <- mean(fitDF[fitDF$Band == bands[band_index],]$U)
                eI[band_index]    <- median(fitDF[fitDF$Band == bands[band_index],]$eI)
                eQ[band_index]    <- median(fitDF[fitDF$Band == bands[band_index],]$eQ)
                eU[band_index]    <- median(fitDF[fitDF$Band == bands[band_index],]$eU)
            }
            freq[band_index]  <- median(fitDF[fitDF$Band == bands[band_index],]$Freq)
        }
        fit <- lm( log(predI) ~ log(freq/100), weights=1.0/eI^2 ); spixI[src_index] <- coef(fit)[2]; I100[src_index] <- exp(coef(fit)[1])
        fit <- lm(0.5*log(predQ^2 + predU^2) ~ log(freq/100), weights=1.0/eI^2 ); spixP[src_index] <- coef(fit)[2]
        f100 <- (freq/100.0)^spixP[src_index]
        fit <- lm(predQ ~ f100+0, weights=1.0/eQ^2); Q100[src_index] <- coef(fit)[1]
        fit <- lm(predU ~ f100+0, weights=1.0/eU^2); U100[src_index] <- coef(fit)[1]
    } else {
        I100[src_index] <- median(srcDF$I)
        Q100[src_index] <- median(srcDF$Q)
        U100[src_index] <- median(srcDF$U)
        spixI[src_index] <- -0.7
        spixP[src_index] <- -0.7
    }
}
srcDF <- na.omit(data.frame(Src=sourceList, RA=RAList, DEC=DecList, I100=I100, Q100=Q100, U100=U100, spixI=spixI, spixP=spixP))
write.csv(srcDF[(abs(srcDF$DEC - ALMA_lat) > 4.0*pi/180.0),], 'PolQuery.CSV', row.names=FALSE)
#-------- LST plot
for(band in c(1,3,4,5,6,7,8,9)){
    plotDF <- plotLST(srcDF, band)
    pLST <- plot_ly(data=plotDF, x = ~LST, y = ~XYcorr, type = 'scatter', mode = 'lines', color=~Src, hoverinfo='text', text=~paste(Src, 'EL=',floor(EL)))
    pLST <- layout(pLST, xaxis=list(showgrid=T, title='LST', nticks=24), yaxis=list(showgrid=T, title='XY correlation [Jy]',rangemode='tozero'), title=sprintf('Band-%d Pol-Calibrator Coverage as of %s (30-day statistics)', band, as.character(Today)))
    htmlFile <- sprintf("Band%d_LSTplot.html", band)
    htmlwidgets::saveWidget(pLST, htmlFile, selfcontained=FALSE)
    rm(plotDF)
}

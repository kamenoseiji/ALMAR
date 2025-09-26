library(RColorBrewer)
library(plotly, warn.conflicts=FALSE)
library(pandoc)
library(htmlwidgets)   # multicore parallelization
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq<-c(43.0,  75.0, 97.5, 132.0,183.0, 233.0, 343.5,400.0, 650.0, 800.0)
Pthresh <-c(0.058, 0.060, 0.069, 0.058, 0.086, 0.077, 0.094, 0.153, 0.442, 0.765) # 5-sigma thresholds for polarized flux
# sigma <- function(integ,tsys){
#    ae <- 0.7*pi*9
#    nant <- 10
#    kb <- 1380
#    bw <- 16e6
#    return(2*kb*tsys/(ae*sqrt(nant*(nant-1)*bw*integ*2)))}
# 5.0* sigma( c(150, 200, 200, 300, 300, 300, 600, 900, 1200, 1200), c(55,65,75,77,115,103,178,353,1179,2041) )
mthresh <- 0.03                                                                   # polarization degree threshold
SECPERDAY <- 86400
ALMA_LAT <- -23.029* pi/180.0   # [radian]
ALMA_LONG <- -67.755* pi/180.0  # [radian]
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/180*pi)
minSinEL <- sin(30/180*pi)
hourPerRad <- 12/pi		# radian to hour angle conversion
sessionDuration <- 3/hourPerRad	# 3 hours in [rad]
pointingDuration <- 0.1/hourPerRad	# 0.1 hourss in [rad]
DateRange <- 60    # 60-day window
#-------- # cos(hour angle) when it passes the given EL
EL_HA <- function(sinEL, dec){
	cosHA <- (sinEL - sin_phi* sin(dec)) / (cos_phi* cos(dec))
	return(cosHA)
}
#-------- Hour angle to Az, El
ha2azel <-function(ha, phi=ALMA_LAT, dec=0){
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
#-------- Estimate Stokes parameters by freqneyc and date 
estimateIQUV <- function(DF, refFreq){
    DF$relFreq <- DF$Freq / refFreq
    timeWeightSoftening <- 5* 86400 # 5-day softening
    if( (max(DF$relFreq) < 0.65) | (max(DF$relFreq) / min(DF$relFreq) < 2.0 )){ return( data.frame( Src=DF$Src[1], I=0.0, eI=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, eP=0.0, EVPA=0.0))}
    DF$timeFreqDeparture <- (abs(DF$relTime) + timeWeightSoftening) * (1.0 + abs(DF$relFreq - 1))
    if( (diff(range(DF$relTime)) < 0.25* DateRange* SECPERDAY) | (max(DF$relTime) < -0.25* SECPERDAY)){      # small number of data
        fitI <- lm(formula=log(I) ~ log(relFreq), data=DF, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
	    fitP <- lm(formula=log(P) ~ log(relFreq), data=DF, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
   } else {
        fitI <- lm(formula=log(I) ~ log(relFreq) + relTime, data=DF, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
	    fitP <- lm(formula=log(P) ~ log(relFreq) + relTime, data=DF, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
    }
    weight <- 1.0/(abs(DF$eEVPA) * sqrt(DF$eQ^2 + DF$eU^2)* abs(log(DF$relFreq) + 1.0)^2 * (timeWeightSoftening / abs(DF$relTime + timeWeightSoftening)))
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*DF$EVPA) ) / sum(weight); Twiddle <- Twiddle/abs(Twiddle)
    IQUV <- data.frame(Src=DF$Src[1], I=exp(coef(fitI)[[1]]), eI=exp(coef(fitI)[[1]])* (exp(coef(summary(fitI))[1,2])-1), Q=0.0, U=0.0, V=0.0, P=exp(coef(fitP)[[1]]), eP=exp(coef(fitP)[[1]])* (exp(coef(summary(fitP))[1,2])-1), EVPA=0.5*Arg(Twiddle))
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
	return( IQUV )
}
#-------- Plot LST
plotLST <- function(DF, band){
    pDF <- data.frame()
    DF <- DF[abs(DF$DEC - ALMA_LAT) > 4.0*pi/180.0,]    # zenith avoidance of 4 degree
    sourceList <- as.character(DF$Src)
    sourceNum <- length(sourceList)
    LST <- seq(0.0, 23.95, 0.05)
    lstNum <- length(LST)
    for(source_index in 1:sourceNum){
        src <- sourceList[source_index]
        HA <- LST*pi/12 - DF$RA[source_index]
        AZEL <- ha2azel( HA, ALMA_LAT, DF$DEC[source_index] )
        AZEL$pa <- AZEL$pa + BandPA[band]
        CS <- cos(2.0* AZEL$pa); SN <- sin(2.0* AZEL$pa)
        AZEL$XYcorr <- abs(DF$U[source_index]*CS - DF$Q[source_index]*SN)
        AZEL[AZEL$el < pi/7.5,]$XYcorr <- NA
        if( nrow(pDF) == 0 ){
            pDF <- data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*180/pi, XYcorr=AZEL$XYcorr)
        } else {
            pDF <- rbind(pDF, data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*180/pi, XYcorr=AZEL$XYcorr))
        }
    }
    return(pDF)
}
#-------- Load Flux.Rdata from web
#FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
#load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
load("Flux.Rdata")     # Data frame of FLDF
FLDF <- FLDF[as.Date(FLDF$Date) > Sys.Date() - DateRange,]  # Data frame within DateRange
#-------- Filter quasars
FLDF <- FLDF[substr(FLDF$Src, 1, 1) == 'J',]    # only quasars
FLDF$P  <- sqrt(FLDF$Q^2 + FLDF$U^2)
FLDF$eP <- sqrt(FLDF$eQ^2 + FLDF$eU^2)
FLDF$EVPA  <- 0.5* atan2(FLDF$U, FLDF$Q)
FLDF$eEVPA <- 0.5* sqrt(FLDF$Q^2 * FLDF$eU^2 + FLDF$U^2 * FLDF$eQ^2) / (FLDF$P^2)
FLDF$relTime <- as.numeric(FLDF$Date) - as.numeric(as.POSIXct(Sys.time()))  # Relative second since now
sourceList <- sort(unique(FLDF$Src))
numSrc <- length(sourceList)
FLDF$medP <- FLDF$freqRange <- numeric(nrow(FLDF))
for(src in sourceList){ FLDF[FLDF$Src == src,]$medP <- median(FLDF[FLDF$Src == src,]$P)}
for(src in sourceList){ FLDF[FLDF$Src == src,]$freqRange <- diff(range(FLDF[FLDF$Src == src,]$Freq))}
FLDF <- FLDF[FLDF$medP > 0.03,]
FLDF <- FLDF[FLDF$freqRange > 100,]
sourceList <- sort(unique(FLDF$Src))
numSrc <- length(sourceList)
#-------- Today's IQUV
for(band in seq(1, 7)){
    srcDF <- data.frame(Src = sourceList, RA=numeric(numSrc), DEC=numeric(numSrc), I=numeric(numSrc), eI=numeric(numSrc), Q=numeric(numSrc), U=numeric(numSrc), P=numeric(numSrc), eP=numeric(numSrc), EVPA=numeric(numSrc), LSTmin=numeric(numSrc), LSTmax=numeric(numSrc), png=character(numSrc))
    for(src in sourceList){
        #-------- Filter by polarized flux
        SDF <- FLDF[FLDF$Src == src,]
        SDF <- SDF[SDF$I < median(SDF$I) + 3.0* sd(SDF$I),]
        SDF <- SDF[SDF$P < median(SDF$P) + 3.0* sd(SDF$P),]
        if(nrow(SDF) < 3){ next }
        #cat(sprintf('%s %f\n', src, max(SDF$Freq)))
        IQUV <- estimateIQUV(SDF, BandFreq[band])
        srcDF[srcDF$Src == src,]$I <- IQUV$I
        srcDF[srcDF$Src == src,]$eI<- IQUV$eI
        srcDF[srcDF$Src == src,]$Q <- IQUV$Q
        srcDF[srcDF$Src == src,]$U <- IQUV$U
        srcDF[srcDF$Src == src,]$P <- IQUV$P
        srcDF[srcDF$Src == src,]$eP<- IQUV$eP
        srcDF[srcDF$Src == src,]$EVPA <- IQUV$EVPA
        srcDF[srcDF$Src == src,]$RA  <- pi* (60.0* as.numeric(substring(src, 2, 3)) + as.numeric(substring(src, 4, 5))) / 720.0
        srcDF[srcDF$Src == src,]$DEC <- pi* sign(as.numeric(substring(src, 6, 10)))* (as.numeric(substring(src, 7, 8)) + as.numeric(substring(src, 9, 10))/60.0) / 180.0
    }
    #-------- LST plot
    plotDF <- plotLST(srcDF[srcDF$P - srcDF$eP > Pthresh/2, ], band)
    pLST <- plot_ly(data=plotDF, x = ~LST, y = ~XYcorr, type = 'scatter', mode = 'lines', color=~Src, hoverinfo='text', text=~paste(Src, 'EL=',floor(EL)))
    pLST <- layout(pLST, xaxis=list(showgrid=T, title='LST', nticks=24), yaxis=list(showgrid=T, title='XY correlation [Jy]',rangemode='tozero'), title=sprintf('Band-%d Pol-Calibrator Coverage as of %s (60-day statistics)', band, as.character(Sys.Date())))
    htmlFile <- sprintf("Band%d_LSTplot.html", band)
    htmlwidgets::saveWidget(pLST, htmlFile, selfcontained=FALSE)
    rm(plotDF)
    #-------- LST window plot
    srcDF <- srcDF[(srcDF$P - srcDF$eP) / (srcDF$I + srcDF$eI) > mthresh,]  # filter by polarization degree
    srcDF <- srcDF[srcDF$P - srcDF$eP > Pthresh[band],]  # filter by polarized flux
    srcDF <- srcDF[srcDF$DEC < ALMA_LAT + pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[srcDF$DEC > ALMA_LAT - pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[abs(srcDF$DEC - ALMA_LAT) > 3.0* pi / 180,]  # avoid zenith passage
    srcDF$ELHA <- acos(EL_HA(minSinEL, srcDF$DEC))  # Hour angle above EL limit (30 deg)
    for(index in 1:nrow(srcDF)){
        HA <- seq(-srcDF$ELHA[index], srcDF$ELHA[index], by=0.004)  # Hour angle above EL limit
		HA_min <- min(HA); HA_max <- max(HA) - sessionDuration	# tentative HA range
        sinHA <- sin(HA)
        cosHA <- cos(HA)
        cos_dec <- cos(srcDF$DEC[index])
        sin_dec <- sin(srcDF$DEC[index])
        PA <- atan2(sinHA, sin_phi*cos_dec/cos_phi - sin_dec*cosHA) + BandPA[band]
        XYcorr <- srcDF$U[index]* cos(2.0*PA) - srcDF$Q[index]* sin(2.0*PA)
		XY_overthresh <- HA[which(abs(XYcorr) > Pthresh[band])]		# HA to obtain |XY| > Pthresh
        XY_intercepts <- HA[which(diff(sign(XYcorr)) != 0)]     # hour angles of sign transition
		if(length(XY_overthresh) * length(XY_intercepts) == 0){
            srcDF$P[index] <- NA
            next
		}
		HA_min <- max(HA_min, min(XY_intercepts) - sessionDuration, min(XY_overthresh) - sessionDuration)	# Cover zero-crossing within 3 hours 
		HA_max <- min(HA_max, max(XY_intercepts) - pointingDuration, max(XY_overthresh) - pointingDuration) # 
        if( HA_min > HA_max ){
            srcDF$P[index] <- NA
            next
        }
        srcDF$LSTmin[index] <- HA_min + srcDF$RA[index]   		# fist LST window
        srcDF$LSTmax[index] <- HA_max + srcDF$RA[index]			# last LST window
		pngFile <- sprintf('%s-Band%d-PA.png', srcDF$Src[index], band)
		png(pngFile)
		plot(HA*hourPerRad, XYcorr, type='l', col='red', xlab='Hour Angle [h]', ylab='XY correlation [Jy]', main=sprintf('%s Band%d', srcDF$Src[index], band))
		grid(nx=NULL, ny=NULL, lty=2, col='gray', lwd=1)
		abline(h=Pthresh[band], lty=2); abline(h=-Pthresh[band], lty=2); abline(v=hourPerRad* XY_intercepts)
		lines(c(HA_min, HA_max)*hourPerRad, c(0,0), lwd=4, col='blue')
		dev.off()
		srcDF$png[index] <- pngFile
    }
    srcDF <- na.omit(srcDF)
    cat(sprintf('Band %d : %d sources\n', band, nrow(srcDF)))
    write.csv(srcDF[, c('Src', 'I', 'P', 'EVPA', 'LSTmin', 'LSTmax')], sprintf('PolCalBand%d.csv', band), row.names=FALSE)
	pngFile <- sprintf('Band%d-LST.png', band)
    png(pngFile, width=1536, height=1024)
    par(mar=c(4,8,3,3))
    LSTplot <- barplot(height=rep(NA, nrow(srcDF)), names=srcDF$Src, horiz=TRUE, las=1, xlim=c(0,24), xlab='LST to start [h]', main=sprintf('Polarization Calibrators at Band %d as of %s', band, Sys.Date()), xaxp=c(0, 24, 24))
    grid(nx=24, ny=0, lwd=0.5, lty=1, col='gray')
    for(index in 1:nrow(srcDF)){ abline(h=LSTplot[index], col='black', lwd=0.5)}
    arrows(hourPerRad* srcDF$LSTmin, LSTplot, hourPerRad* srcDF$LSTmax, LSTplot, length=0, lwd=12, col='red')
    for(index in 1:nrow(srcDF)){
        if( srcDF$LSTmin[index] < 0.0){
            arrows(hourPerRad* srcDF$LSTmin[index] + 24.0, LSTplot[index], 24.0, LSTplot[index], length=0, lwd=12, col='red')
            text(hourPerRad*srcDF$LSTmin[index] + 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmin[index] + 24.0), pos=2) 
        } else {
            text(hourPerRad*srcDF$LSTmin[index], LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmin[index]), pos=2) 
        }
        if( srcDF$LSTmax[index] > 2.0*pi){
            arrows(0.0, LSTplot[index], hourPerRad* srcDF$LSTmax[index] - 24.0, LSTplot[index], length=0, lwd=12, col='red')
            text(hourPerRad*srcDF$LSTmax[index] - 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmax[index] - 24.0), pos=4) 
        } else {
            text(hourPerRad*srcDF$LSTmax[index], LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmax[index]), pos=4) 
        }
    }
    dev.off()
}

#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(43.0, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
ALMA_LAT <- -23.029* pi/180.0   # [radian]
ALMA_LONG <- -67.755* pi/180.0  # [radian]
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/180*pi)
minSinEL <- sin(20/180*pi)
#-------- # cos(hour angle) when it passes the given EL
EL_HA <- function(sinEL, dec){
	cosHA <- (sinEL - sin_phi* sin(dec)) / (cos_phi* cos(dec))
	return(cosHA)
}
#-------- Estimate Stokes parameters by freqneyc and date 
estimateIQUV <- function(DF, refFreq){
    DF$relFreq <- DF$Freq / refFreq
    timeWeightSoftening <- 5* 86400 # 5-day softening
    if(diff(range(DF$Freq)) < 100 ){ return( data.frame( Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=0.0))}
    if((length(unique(DF$relTime)) < 4) | (max(DF$relTime) < -15* 86400)){
        fitI <- lm(formula=log(I) ~ log(relFreq), data=DF, weight=(I / eI)^2 * (timeWeightSoftening / abs(relTime + timeWeightSoftening)) )
	    fitP <- lm(formula=log(P) ~ log(relFreq), data=DF, weight=(DF$P/DF$eP)^2 * (timeWeightSoftening / abs(relTime + timeWeightSoftening)))
    } else {
        fitI <- lm(formula=log(I) ~ relTime + log(relFreq), data=DF, weight=(I / eI)^2 * (timeWeightSoftening / abs(relTime + timeWeightSoftening)) )
	    fitP <- lm(formula=log(P) ~ relTime + log(relFreq), data=DF, weight=(P/eP)^2 * (timeWeightSoftening / abs(relTime + timeWeightSoftening)))
    }
    weight <- 1.0/(abs(DF$eEVPA)^2 * abs(log(DF$relFreq) + 1.0)^2 * (timeWeightSoftening / abs(DF$relTime + timeWeightSoftening)))
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*DF$EVPA) ) / sum(weight); Twiddle <- Twiddle/abs(Twiddle)
    IQUV <- data.frame(Src=DF$Src[1], I=exp(coef(fitI)[[1]]), Q=0.0, U=0.0, V=0.0, P=exp(coef(fitP)[[1]]), EVPA=0.5*Arg(Twiddle))
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
	return( IQUV )
}
#-------- Load Flux.Rdata from web
Pthresh <- 0.03    # polarized flux > 30 mJy
DateRange <- 60    # 60-day window
FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
#load("Flux.Rdata")     # Data frame of FLDF
FLDF <- FLDF[as.Date(FLDF$Date) > Sys.Date() - DateRange,]  # Data frame within DateRange
#-------- Filter quasars
FLDF <- FLDF[substr(FLDF$Src, 1, 1) == 'J',]
FLDF$P  <- sqrt(FLDF$Q^2 + FLDF$U^2)
FLDF$eP <- sqrt(FLDF$eQ^2 + FLDF$eU^2)
FLDF$EVPA  <- 0.5* atan2(FLDF$U, FLDF$Q)
FLDF$eEVPA <- 0.5* sqrt(FLDF$Q^2 * FLDF$eU^2 + FLDF$U^2 * FLDF$eQ^2) / (FLDF$P^2)
FLDF$relTime <- as.numeric(FLDF$Date) - as.numeric(as.POSIXct(Sys.time()))  # Relative second since now
sourceList <- sort(unique(FLDF$Src))
numSrc <- length(sourceList)
#-------- Today's IQUV
#for(band in seq(1, 9)){
for(band in seq(1, 7)){
    srcDF <- data.frame(Src = sourceList, RA=numeric(numSrc), DEC=numeric(numSrc), I=numeric(numSrc), Q=numeric(numSrc), U=numeric(numSrc), P=numeric(numSrc), EVPA=numeric(numSrc), LSTmin=numeric(numSrc), LSTmax=numeric(numSrc), png=character(numSrc))
    for(src in sourceList){
        #-------- Filter by polarized flux
        SDF <- FLDF[FLDF$Src == src,]
        IQUV <- estimateIQUV(SDF, BandFreq[band])
        srcDF[srcDF$Src == src,]$I <- IQUV$I
        srcDF[srcDF$Src == src,]$Q <- IQUV$Q
        srcDF[srcDF$Src == src,]$U <- IQUV$U
        srcDF[srcDF$Src == src,]$P <- IQUV$P
        srcDF[srcDF$Src == src,]$EVPA <- IQUV$EVPA
        srcDF[srcDF$Src == src,]$RA  <- pi* (60.0* as.numeric(substring(src, 2, 3)) + as.numeric(substring(src, 4, 5))) / 720.0
        srcDF[srcDF$Src == src,]$DEC <- pi* sign(as.numeric(substring(src, 6, 10)))* (as.numeric(substring(src, 7, 8)) + as.numeric(substring(src, 9, 10))/60.0) / 180.0
    }
    srcDF <- srcDF[srcDF$P > Pthresh,]  # filter by polarized flux
    srcDF <- srcDF[srcDF$DEC < ALMA_LAT + pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[srcDF$DEC > ALMA_LAT - pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[abs(srcDF$DEC - ALMA_LAT) > 3.0* pi / 180,]  # avoid zenith passage
    srcDF$ELHA <- acos(EL_HA(minSinEL, srcDF$DEC))  # Hour angle when EL > 20 deg
    cat(sprintf('Band %d : %d sources\n', band, nrow(srcDF)))
    for(index in 1:nrow(srcDF)){
        HA <- seq(-srcDF$ELHA[index], srcDF$ELHA[index], by=0.004)  # Hour angle for EL > 20 deg
		HA_min <- min(HA); HA_max <- max(HA) - 2/24*2*pi	# tentative HA range
        sinHA <- sin(HA)
        cosHA <- cos(HA)
        cos_dec <- cos(srcDF$DEC[index])
        sin_dec <- sin(srcDF$DEC[index])
        PA <- atan2(sinHA, sin_phi*cos_dec/cos_phi - sin_dec*cosHA) + BandPA[band]
        XYcorr <- srcDF$U[index]* cos(2.0*PA) - srcDF$Q[index]* sin(2.0*PA)
		XY_overthresh <- HA[which(abs(XYcorr) > Pthresh)]		# HA to obtain |XY| > Pthresh
		HA_min <- max(min(XY_overthresh) - 2/24*2*pi, min(HA))
		HA_max <- min(HA_max, max(XY_overthresh) - 2/24*2*pi)				# Cover > Pthresh within 2 hours
        XY_intercepts <- HA[which(diff(sign(XYcorr)) != 0)]     # hour angles of sign transition
		HA_min <- max(HA_min, min(XY_intercepts) - 2/24*2*pi)	# Cover zero-crossing within 2 hours 
		HA_max <- min(HA_max, max(XY_intercepts)) - 0.026		# Cover zero-crossing 10 minutes ago
        if( HA_min > HA_max ){
            srcDF$P[index] <- NA
            next
        }
        srcDF$LSTmin[index] <- HA_min + srcDF$RA[index]   		# fist LST window
        srcDF$LSTmax[index] <- HA_max + srcDF$RA[index]			# last LST window
		pngFile <- sprintf('%s-Band%d-PA.png', srcDF$Src[index], band)
		png(pngFile)
		plot(HA, XYcorr, type='l', col='red', xlab='Hour Angle [rad]', ylab='XY correlation [Jy]', main=sprintf('%s Band%d', srcDF$Src[index], band))
		abline(h=Pthresh, lty=2); abline(h=-Pthresh, lty=2); abline(v=XY_intercepts)
		lines(c(HA_min, HA_max), c(0,0), lwd=4, col='blue')
		dev.off()
		srcDF$png[index] <- pngFile
    }
    srcDF <- na.omit(srcDF)
    write.csv(srcDF[, c('Src', 'I', 'P', 'EVPA', 'LSTmin', 'LSTmax', 'png')], sprintf('PolCalBand%d.csv', band), row.names=FALSE)
}

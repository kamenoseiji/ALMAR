#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(43.0, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
Pthresh  <- c(0.06, 0.07, 0.07, 0.08, 0.10, 0.09, 0.13, 0.3, 1.1, 1.7)     # thresholds for polarized flux
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
#-------- Estimate Stokes parameters by freqneyc and date 
estimateIQUV <- function(DF, refFreq){
    DF$relFreq <- DF$Freq / refFreq
    timeWeightSoftening <- 5* 86400 # 5-day softening
    if( (max(DF$relFreq) < 0.65) | (max(DF$relFreq) / min(DF$relFreq) < 2.0 )){ return( data.frame( Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=0.0))}
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
    IQUV <- data.frame(Src=DF$Src[1], I=exp(coef(fitI)[[1]]), Q=0.0, U=0.0, V=0.0, P=exp(coef(fitP)[[1]]), EVPA=0.5*Arg(Twiddle))
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
	return( IQUV )
}
#-------- Load Flux.Rdata from web
#FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
#load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
load("Flux.Rdata")     # Data frame of FLDF
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
        SDF <- SDF[SDF$eI < 0.1* SDF$I,]
        SDF <- SDF[SDF$eP < 0.5* SDF$P,]
        IQUV <- estimateIQUV(SDF, BandFreq[band])
        srcDF[srcDF$Src == src,]$I <- IQUV$I
        srcDF[srcDF$Src == src,]$Q <- IQUV$Q
        srcDF[srcDF$Src == src,]$U <- IQUV$U
        srcDF[srcDF$Src == src,]$P <- IQUV$P
        srcDF[srcDF$Src == src,]$EVPA <- IQUV$EVPA
        srcDF[srcDF$Src == src,]$RA  <- pi* (60.0* as.numeric(substring(src, 2, 3)) + as.numeric(substring(src, 4, 5))) / 720.0
        srcDF[srcDF$Src == src,]$DEC <- pi* sign(as.numeric(substring(src, 6, 10)))* (as.numeric(substring(src, 7, 8)) + as.numeric(substring(src, 9, 10))/60.0) / 180.0
    }
    srcDF <- srcDF[srcDF$P > Pthresh[band],]  # filter by polarized flux
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
}

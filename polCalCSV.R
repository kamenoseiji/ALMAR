library(RColorBrewer)
library(plotly, warn.conflicts=FALSE)
library(pandoc)
library(htmlwidgets)   # multicore parallelization
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq<-c(43.0,  75.0, 97.5, 132.0,183.0, 233.0, 343.5,400.0, 650.0, 800.0)
Pthresh7<-c(0.058, 0.060, 0.069, 0.058, 0.086, 0.077, 0.094, 0.153, 0.442, 0.765) # for 7m array, 5-sigma thresholds for polarized flux
Pthresh12 <- 0.5* Pthresh7     # for 12m array 
# sigma <- function(integ,tsys){
#    ae <- 0.7*pi*9
#    nant <- 10
#    kb <- 1380
#    bw <- 16e6
#    return(2*kb*tsys/(ae*sqrt(nant*(nant-1)*bw*integ*2)))}
# 5.0* sigma( c(150, 200, 200, 300, 300, 300, 600, 900, 1200, 1200), c(55,65,75,77,115,103,178,353,1179,2041) )
mthresh <- 0.03                                                                   # polarization degree threshold
SECPERDAY <- 86400
hourPerRad <- 12/pi		# radian to hour angle conversion
RADDEG <- 180.0/pi
HAresolution <- 0.0025 # [radian] 0.0025 rad = 34.38 s
ALMA_LAT <- -23.029 / RADDEG  # [radian]
ALMA_LONG <- -67.755/ RADDEG  # [radian]
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/RADDEG)
minSinEL <- sin(30/RADDEG)
sessionDuration <- 1.0/hourPerRad	# at least 1 hour for session, in [rad]
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
    DF <- DF[abs(DF$DEC - ALMA_LAT) > 4.0/RADDEG,]    # zenith avoidance of 4 degree
    sourceList <- as.character(DF$Src)
    sourceNum <- length(sourceList)
    LST <- seq(0.0, 23.95, 0.05)
    lstNum <- length(LST)
    for(source_index in 1:sourceNum){
        src <- sourceList[source_index]
        HA <- LST/hourPerRad - DF$RA[source_index]
        AZEL <- ha2azel( HA, ALMA_LAT, DF$DEC[source_index] )
        AZEL$pa <- AZEL$pa + BandPA[band]
        CS <- cos(2.0* AZEL$pa); SN <- sin(2.0* AZEL$pa)
        AZEL$XYcorr <- abs(DF$U[source_index]*CS - DF$Q[source_index]*SN)
        AZEL[AZEL$el < pi/7.5,]$XYcorr <- NA        # flag EL < 24 deg out
        if( nrow(pDF) == 0 ){
            pDF <- data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*RADDEG, XYcorr=AZEL$XYcorr)
        } else {
            pDF <- rbind(pDF, data.frame(LST=LST, Src=rep(src, lstNum), EL=AZEL$el*180/pi, XYcorr=AZEL$XYcorr))
        }
    }
    return(pDF)
}
#-------- Filter calibrators for band
srcFreqCalibrator <- function(DF, band){
    sourceList <- sort(unique(DF$Src))
    numSrc <- length(sourceList)
    srcDF <- data.frame(Src=sourceList, I=numeric(numSrc), eI=numeric(numSrc), Q=numeric(numSrc), U= numeric(numSrc), V=numeric(numSrc), P=numeric(numSrc), eP=numeric(numSrc), EVPA=numeric(numSrc))
    for(src_index in 1:numSrc){
        src <- sourceList[src_index]
        SDF <- DF[((DF$Src == src) & (DF$I < median(DF$I) + 3.0* sd(DF$I)) & (DF$P < median(DF$P) + 3.0* sd(DF$P))),]   # Filter reliable data
        if(nrow(SDF) < 3){ next }
        srcDF[src_index,] <- estimateIQUV(SDF, BandFreq[band])
    }
    srcDF$RA  <- pi* (60.0* as.numeric(substring(sourceList, 2, 3)) + as.numeric(substring(sourceList, 4, 5))) / 720.0
    srcDF$DEC <- pi* sign(as.numeric(substring(sourceList, 6, 10)))* (as.numeric(substring(sourceList, 7, 8)) + as.numeric(substring(sourceList, 9, 10))/60.0) / 180.0
    return( srcDF[((srcDF$P - srcDF$eP > Pthresh12[band]) & ((srcDF$P - srcDF$eP)/(srcDF$I + srcDF$eI) > 0.03)),] )     # Filter by polarized flux and polarization degree
}
#-------- HA range over threshold
HArange <- function(df, thresh12, thresh7, BPA){
    #df$HA7st1  <- df$HA7st2  <- df$HA7et  <- numeric(1)
    #df$HA12st1 <- df$HA12st2 <- df$HA12et <- numeric(1)
    DF <- df    # for output
    cos_dec <- cos(df$DEC); sin_dec <- sin(df$DEC)
    HApointer  <- which(names(srcDF) == 'HA7et')
    LSTpointer <- which(names(srcDF) == 'LST7et')
    HA <- seq(-df$ELHA, df$ELHA + HAresolution, by=HAresolution)  # Hour angle above EL limit
    sinHA <- sin(HA); cosHA <- cos(HA)
    PA <- atan2(sinHA, sin_phi*cos_dec/cos_phi - sin_dec*cosHA) + BPA   # Parallactic angle + BandPA
    calDF <- data.frame(HA=HA, PA=PA, XYcorr=df$U* cos(2.0*PA) - df$Q* sin(2.0*PA), et12=NA)
    HA_intercepts <- HA[which(diff(sign(calDF$XYcorr)) != 0)]     # hour angles of sign transition
    HA_overthresh12 <- calDF[abs(calDF$XYcorr) > thresh12,]$HA
    calDF <- calDF[calDF$HA < max(HA_intercepts) - pointingDuration,]          # start time must be before the last intercept
    for(index in 1:nrow(calDF)){
        if(calDF[index,]$HA %in% HA_overthresh12){
            calDF[index,]$et12 <- min(HA_intercepts[which(HA_intercepts > calDF[index,]$HA)])
        } else if(calDF[index,]$HA < max(HA_intercepts)){
            calDF[index,]$et12 <- max(min(HA_intercepts[which(HA_intercepts > calDF[index,]$HA)]), min(HA_overthresh12[which(HA_overthresh12 > calDF[index,]$HA)]))
        }
    }
    calDF <- na.omit(calDF)
    et12List <- unique(calDF$et12)
    for(index in seq_along(et12List)){
        DF[index,] <- DF[1,]
        DF$HA12st1[index] <- min(calDF$HA[which(calDF$et12 == et12List[index])])
        DF$HA12st2[index] <- max(calDF$HA[which(calDF$et12 == et12List[index])])
        DF$HA12et[index]  <- et12List[index]
    }




    if(length(HA_intercepts) < 1){ df[HApointer:(HApointer+5)] <- df[LSTpointer:(LSTpointer+5)] <- NA; return(df)}   # No zero-intercept
    XY_transit12 <- HA[which(diff(sign(abs(calDF$XYcorr) - thresh12)) != 0.0)]  # XYcorr straddles threshold
    XY_transit7 <- HA[which(diff(sign(abs(calDF$XYcorr) - thresh7)) != 0.0)]  # XYcorr straddles threshold
    if(length(XY_transit12) < 1){ df[HApointer:(HApointer+5)] <- df[LSTpointer:(LSTpointer+5)] <- NA; return(df)}   # No XY > thresh

	plot((calDF$HA + df$RA)*hourPerRad, calDF$XYcorr, type='l', col='darkgreen', xlab='LST [h]', ylab='XY correlation [Jy]', main=sprintf('%s Band%d', df$Src, band))
	grid(nx=NULL, ny=NULL, lty=2, col='gray', lwd=1)
	abline(h=thresh12, lty=2, col='blue'); abline(h=-thresh12, lty=2, col='blue'); abline(h=thresh7, lty=2, col='red'); abline(h=-thresh7, lty=2, col='red'); abline(v=hourPerRad* (HA_intercepts + df$RA))
    text((min(HA)+df$RA)*hourPerRad+0.1, thresh12, '12-m threshold', col='blue', pos=3); text((min(HA)+df$RA)*hourPerRad+0.1, thresh7, '7-m threshold', col='red', pos=3)
    for(intercept in HA_intercepts){ text((intercept + df$RA)*hourPerRad, min(calDF$XYcorr), sprintf('%.1fh', (intercept + df$RA)*hourPerRad), pos=4, srt=90) }



    #-------- for 12m threshold
    HA_XY12 <- calDF[abs(calDF$XYcorr) > thresh12,]$HA - HAresolution* ceiling(pointingDuration/HAresolution)    # HA range for XY > thresh, 0.1 hour ahead
    calDF$et12 <- ifelse(length(HA_XY12) > 0, max(HA), NA)
    for(index in 1:nrow(calDF)){
        calDF$et12[index] <- ifelse( calDF$HA[index] %in% HA_XY12, max(calDF$HA[index]+sessionDuration, min(HA_intercepts)),  max(calDF$HA[index]+sessionDuration, min(HA_XY12), min(HA_intercepts)))
    }
    flatET12 <- which(diff(calDF$et12) < 0.5*HAresolution)
    lineET12 <- which(diff(calDF$et12) >= 0.5*HAresolution)
    numWindow <- sign(length(flatET12)) + sign(length(lineET12))
    for(row_index in 2:numWindow){ DF[row_index,] <- DF[1,] }
    row_index <- 1
    if(length(flatET12) > 0){
        DF[row_index,]$HA12st1 <- calDF$HA[min(flatET12)]
        DF[row_index,]$HA12st2 <- calDF$HA[max(flatET12)]
        DF[row_index,]$HA12et  <- calDF$et12[max(flatET12)]
        row_index <- row_index + 1
    }
    if(length(lineET12) > 0){
        DF[row_index,]$HA12st1 <- calDF$HA[min(lineET12)]
        DF[row_index,]$HA12st2 <- calDF$HA[max(lineET12)]
        DF[row_index,]$HA12et  <- calDF$et12[max(lineET12)]
    }
    DF$LST12st1 <- DF$HA12st1 + DF$RA
    DF$LST12st2 <- DF$HA12st2 + DF$RA
    DF$LST12et  <- DF$HA12et  + DF$RA
    for(row_index in 1:numWindow){
	    lines(c(DF[row_index,]$LST12st1, DF[row_index,]$LST12st2)*hourPerRad, 0.005*c(row_index-1, row_index-1), lwd=4, col='blue')
	    lines(c(DF[row_index,]$LST12st2, DF[row_index,]$LST12et)*hourPerRad,  0.005*c(row_index-1, row_index-1), lwd=2, lty=2, col='blue')
        text(DF[row_index,]$LST12st1*hourPerRad, 0.01*row_index, sprintf('%.1fh', DF[row_index,]$LST12st1*hourPerRad), pos=1, col='blue')
        text(DF[row_index,]$LST12st2*hourPerRad, 0.0, sprintf('%.1fh', DF[row_index,]$LST12st2*hourPerRad), pos=1, col='blue')
        text(DF[row_index,]$LST12et*hourPerRad, 0.0,  sprintf('%.1fh', DF[row_index,]$LST12et*hourPerRad), pos=1, col='blue')
    }



    #-------- for 7m threshold
    if(length(XY_transit7)  > 0){
        HA_XY7 <- calDF[abs(calDF$XYcorr) > thresh7,]$HA - HAresolution* ceiling(pointingDuration/HAresolution)    # HA range for XY > thresh, 0.1 hour ahead
        calDF$et7 <- ifelse( length(HA_XY7) > 0, max(HA), NA)
        for(index in 1:nrow(calDF)){
            calDF$et7[index]  <- ifelse( calDF$HA[index] %in% HA_XY7 , max(calDF$HA[index]+sessionDuration, min(HA_intercepts)),  max(calDF$HA[index]+sessionDuration, min(HA_XY7), min(HA_intercepts)))
        }
        flatET7 <- which(diff(calDF$et7) < 0.5*HAresolution)
        lineET7 <- which(diff(calDF$et7) >= 0.5*HAresolution)
        if(length(flatET7) > 0){
            df$HA7st1 <- calDF$HA[min(flatET7)]
            df$HA7st2 <- calDF$HA[max(flatET7)]
            df$HA7et  <- calDF$et7[min(flatET7)]
        }
        if(length(lineET7) > 0){
            df$HA7st1 <- append(df$HA7st1, calDF$HA[min(lineET7)])
            df$HA7st2 <- append(df$HA7st1, calDF$HA[max(lineET7)])
            df$HA7et  <- append(df$HA7st1, calDF$et7[min(lineET7)])
        }
    }

    df$HA12min <- max(df$HA12min, min(HA_intercepts) - sessionDuration, min(calDF12$HA) - sessionDuration)	 # Cover zero-crossing within 3 hours 
	df$HA12max <- min(df$HA12max, max(HA_intercepts) - pointingDuration, max(calDF12$HA) - pointingDuration) # 
    df$LSTmin <- df$HA12min + df$RA 
    df$LSTmax <- df$HA12max + df$RA 
    if(nrow(calDF7) > 1){
        XY_transit7 <- HA[which(diff(sign(abs(calDF$XYcorr) - thresh7)) != 0.0)]  # XYcorr straddles threshold
        df$HA7min <- max(df$HA7min, min(HA_intercepts) - sessionDuration, min(calDF7$HA) - sessionDuration)	 # Cover zero-crossing within 3 hours 
	    df$HA7max <- min(df$HA7max, max(HA_intercepts) - pointingDuration, max(calDF7$HA) - pointingDuration) # 
        df$LST7min <- df$HA7min + df$RA 
        df$LST7max <- df$HA7max + df$RA 
    }
	plot((calDF$HA + df$RA)*hourPerRad, calDF$XYcorr, type='l', col='darkgreen', xlab='LST [h]', ylab='XY correlation [Jy]', main=sprintf('%s Band%d', df$Src, band))
	grid(nx=NULL, ny=NULL, lty=2, col='gray', lwd=1)
	abline(h=thresh12, lty=2, col='blue'); abline(h=-thresh12, lty=2, col='blue'); abline(h=thresh7, lty=2, col='red'); abline(h=-thresh7, lty=2, col='red'); abline(v=hourPerRad* (HA_intercepts + df$RA))
    if(nrow(calDF7) > 1){
	    lines(c(df$LST7min, df$LST7max)*hourPerRad, c(0,0), lwd=12, col='red'); text(df$LST7min*hourPerRad, 0.0, sprintf('%.1fh', df$LST7min*hourPerRad), pos=3, col='red'); text(df$LST7max*hourPerRad, 0.0, sprintf('%.1fh', df$LST7max*hourPerRad), pos=3, col='red')
    }
	lines(c(df$LSTmin, df$LSTmax)*hourPerRad, c(0,0), lwd=4, col='blue'); text(df$LSTmin*hourPerRad, 0.0, sprintf('%.1fh', df$LSTmin*hourPerRad), pos=1, col='blue'); text(df$LSTmax*hourPerRad, 0.0, sprintf('%.1fh', df$LSTmax*hourPerRad), pos=1, col='blue')
    text((min(HA)+df$RA)*hourPerRad+0.1, thresh12, '12-m threshold', col='blue', pos=3)
    text((min(HA)+df$RA)*hourPerRad+0.1, thresh7, '7-m threshold', col='red', pos=3)
    for(intercept in HA_intercepts){ text((intercept + df$RA)*hourPerRad, min(calDF$XYcorr), sprintf('%.1fh', (intercept + df$RA)*hourPerRad), pos=4, srt=90) }
    for(transit in XY_transit12){ text((transit + df$RA)*hourPerRad, calDF[calDF$HA == transit,]$XYcorr, sprintf('%.1fh', (transit + df$RA)*hourPerRad), adj=0, col='blue') }
    if(nrow(calDF7) > 1){
        for(transit in XY_transit7){ text((transit + df$RA)*hourPerRad, calDF[calDF$HA == transit,]$XYcorr, sprintf('%.1fh', (transit + df$RA)*hourPerRad), adj=0, col='red') }
    }
    return( df )
}
#-------- Flag LST range
LSTfrag <- function(df){
    LSTrange <- pi* seq(-1, 1440, by=1)/720.0  # 1-min unit, in [rad]
    LSTfrag <- rep(1, length(LSTrange))
    df <- df[!is.na(df$LST7min),]
    df$LST7min <- df$LST7min %% (2*pi)
    df$LST7max <- df$LST7max %% (2*pi)
    for(index in 1:nrow(df)){
        if(df$LST7min[index] < df$LST7max[index]){
            LSTindex <- which((LSTrange > df$LST7min[index]) & (LSTrange < df$LST7max[index]))
        } else {
            LSTindex <- which((LSTrange > df$LST7min[index]) | (LSTrange < df$LST7max[index]))
        }
        LSTfrag[LSTindex] <- 0
    }
    LSTfrag[1] <- LSTfrag[length(LSTfrag)] <- 1
    return(data.frame(LSTbegin = LSTrange[which(diff(LSTfrag) == -1)], LSTend=LSTrange[which(diff(LSTfrag) == 1)]))
}
#-------- Load Flux.Rdata from web
#FluxDataURL <- "https://www.alma.cl/~skameno/AMAPOLA/"
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
#-------- Loop in frequency band
for(band in seq(1, 7)){
    #-------- Today's IQUV
    srcDF <- srcFreqCalibrator(FLDF, band)  # source properties (I, Q, U, V, P, EVPA) at the band
    #-------- XY-LST plot
    plotDF <- plotLST(srcDF[srcDF$P - srcDF$eP > Pthresh12[band], ], band)
    pLST <- plot_ly(data=plotDF, x = ~LST, y = ~XYcorr, type = 'scatter', mode = 'lines', color=~Src, hoverinfo='text', text=~paste(Src, 'EL=',floor(EL)))
    pLST <- layout(pLST, xaxis=list(showgrid=T, title='LST', nticks=24), yaxis=list(showgrid=T, title='XY correlation [Jy]',rangemode='tozero'), title=sprintf('Band-%d Pol-Calibrator Coverage as of %s (60-day statistics)', band, as.character(Sys.Date())))
    htmlFile <- sprintf("Band%d_LSTplot.html", band)
    htmlwidgets::saveWidget(pLST, htmlFile, selfcontained=FALSE)
    rm(plotDF)
    #-------- Calculate HA and LST windows
    srcDF <- srcDF[(srcDF$P - srcDF$eP) / (srcDF$I + srcDF$eI) > mthresh,]  # filter by polarization degree
    srcDF <- srcDF[srcDF$P - srcDF$eP > Pthresh12[band],]  # filter by polarized flux
    srcDF <- srcDF[srcDF$DEC < ALMA_LAT + pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[srcDF$DEC > ALMA_LAT - pi/3,]  # max EL > 30 deg
    srcDF <- srcDF[abs(srcDF$DEC - ALMA_LAT) > 3.0/RADDEG,]  # avoid zenith passage
    srcDF$ELHA <- acos(EL_HA(minSinEL, srcDF$DEC))  # Hour angle [rad] above EL limit (30 deg)
    srcDF$HA12st1 <- srcDF$HA12st2 <- srcDF$HA12et <- srcDF$HA7st1 <- srcDF$HA7st2 <- srcDF$HA7et <- 0.0
    srcDF$LST12st1 <- srcDF$LST12st2 <- srcDF$LST12et <- srcDF$LST7st1 <- srcDF$LST7st2 <- srcDF$LST7et <- 0.0
    srcDF$png <- ''
    for(index in 1:nrow(srcDF)){
		pngFile <- sprintf('%s-Band%d-PA.png', srcDF$Src[index], band)
		png(pngFile, width=1024, height=768)
        srcDF[index,] <- HArange(srcDF[index,], Pthresh12[band], Pthresh7[band], BandPA[band])
		dev.off()
		srcDF$png[index] <- pngFile
    }
    srcDF <- srcDF[!is.na(srcDF$HA12min),]
    LSTwindow <- LSTfrag(srcDF)
    cat(sprintf('Band %d : %d sources\n', band, nrow(srcDF)))
    write.csv(srcDF[, c('Src', 'I', 'P', 'EVPA', 'LSTmin', 'LSTmax')], sprintf('PolCalBand%d.csv', band), row.names=FALSE)
	pngFile <- sprintf('Band%d-LST.png', band)
    png(pngFile, width=1536, height=1024)
    par(mar=c(4,8,3,3))
    LSTplot <- barplot(height=rep(NA, nrow(srcDF)), names=srcDF$Src, horiz=TRUE, las=1, xlim=c(0,24), xlab='LST to start [h]', main=sprintf('Polarization Calibrators at Band %d as of %s', band, Sys.Date()), xaxp=c(0, 24, 24))
    for(index in 1:nrow(LSTwindow)){
        rect( hourPerRad*LSTwindow$LSTbegin[index], 0, hourPerRad*LSTwindow$LSTend[index], max(LSTplot)+1, col='lightgreen', border='transparent')
        text_sd <- sprintf('LST=[%.1fh - %.1fh]', hourPerRad*abs(LSTwindow$LSTbegin[index]), hourPerRad*LSTwindow$LSTend[index])
        text(0.5*hourPerRad*(LSTwindow$LSTbegin[index] + LSTwindow$LSTend[index]), -0.5, text_sd, cex=2, col='darkgreen')
    }
    grid(nx=24, ny=0, lwd=0.5, lty=1, col='gray')
    for(index in 1:nrow(srcDF)){
        segColor <- ifelse(srcDF$P[index] - srcDF$eP[index] > Pthresh7[band], 'red', 'blue')
        segWidth <- min(c(100*srcDF$P[index], 20))
        abline(h=LSTplot[index], col='black', lwd=0.5)
        arrows(hourPerRad* srcDF$LSTmin[index], LSTplot[index], hourPerRad* srcDF$LSTmax[index], LSTplot[index], length=0, lwd=segWidth, col='blue')
        if( srcDF$LSTmin[index] < 0.0){
            arrows(hourPerRad* srcDF$LSTmin[index] + 24.0, LSTplot[index], 24.0, LSTplot[index], length=0, lwd=segWidth, col='blue')
            text(hourPerRad*srcDF$LSTmin[index] + 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmin[index] + 24.0), pos=2) 
        } else {
            text(hourPerRad*srcDF$LSTmin[index], LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmin[index]), pos=2) 
        }
        if( srcDF$LSTmax[index] > 2.0*pi){
            arrows(0.0, LSTplot[index], hourPerRad* srcDF$LSTmax[index] - 24.0, LSTplot[index], length=0, lwd=segWidth, col='blue')
            text(hourPerRad*srcDF$LSTmax[index] - 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmax[index] - 24.0), pos=4) 
        } else {
            text(hourPerRad*srcDF$LSTmax[index], LSTplot[index], sprintf('%.1f', hourPerRad*srcDF$LSTmax[index]), pos=4) 
        }
        legend('right',  legend=c('12-m or 7-m Array','Only 12-m Array','7-m Array Window'), col=c('red','blue', 'lightgreen'), lty=1, lwd=c(2,2,8))
        if( !is.na(srcDF$LST7min[index])){             # Usable for 7m
            arrows(hourPerRad* srcDF$LST7min[index], LSTplot[index], hourPerRad* srcDF$LST7max[index], LSTplot[index], length=0, lwd=segWidth, col='red')
            if( srcDF$LST7min[index] < 0.0){
                arrows(hourPerRad* srcDF$LST7min[index] + 24.0, LSTplot[index], 24.0, LSTplot[index], length=0, lwd=segWidth, col='red')
            }
            if( srcDF$LST7max[index] > 2.0* pi){
                arrows(0.0, LSTplot[index], hourPerRad* srcDF$LST7max[index] - 24.0, LSTplot[index], length=0, lwd=segWidth, col='red')
            }
        }
    }
    dev.off()
}

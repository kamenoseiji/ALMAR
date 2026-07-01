library(RCurl)
library(RColorBrewer)
library(plotly, warn.conflicts=FALSE)
library(pandoc)
library(xtable)
library(htmlwidgets)   # multicore parallelization
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/ALMAR/refs/heads/master/StatStokes.R", ssl.verifypeer = FALSE)))
#source('../StatStokes.R')
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
HAresolution <- 0.0025 # [radian] 0.0025 rad = 34.38 s
ALMA_LAT <- -23.029 / RADDEG  # [radian]
ALMA_LONG <- -67.755/ RADDEG  # [radian]
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/RADDEG)
minSinEL <- sin(30/RADDEG)
sessionDuration <- 1.0/hourPerRad	# at least 1 hour for session, in [rad]
pointingDuration <- 0.2/hourPerRad	# 0.2 hour in [rad]
etMargin <- 0.5/hourPerRad          # 0.5 hour margin to end the session, in [rad] 
DateRange <- 60    # 60-day window
Today <- Sys.time()
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
    srcDF <- data.frame(Src=sourceList, I=numeric(numSrc), Q=numeric(numSrc), U=numeric(numSrc), V=numeric(numSrc), P=numeric(numSrc), EVPA=numeric(numSrc), eI=numeric(numSrc), eQ=numeric(numSrc), eU=numeric(numSrc), eV=numeric(numSrc), eP=numeric(numSrc), eEVPA=numeric(numSrc))
    for(src in sourceList){
        SDF <- DF[((DF$Src == src) & (DF$I < median(DF$I) + 3.0* sd(DF$I)) & (DF$P < median(DF$P) + 3.0* sd(DF$P))),]   # Filter reliable data
        if(nrow(SDF) < 3){ next }
        srcDF[srcDF$Src == src,] <- estimateIQUV(SDF, BandFreq[band], Today)
    }
    srcDF$RA  <- pi* (60.0* as.numeric(substring(sourceList, 2, 3)) + as.numeric(substring(sourceList, 4, 5))) / 720.0
    srcDF$DEC <- pi* sign(as.numeric(substring(sourceList, 6, 10)))* (as.numeric(substring(sourceList, 7, 8)) + as.numeric(substring(sourceList, 9, 10))/60.0) / 180.0
    return( srcDF[((srcDF$P - srcDF$eP > Pthresh12[band]) & ((srcDF$P - srcDF$eP)/(srcDF$I + srcDF$eI) > 0.03)),] )     # Filter by polarized flux and polarization degree
}
#-------- HA range over threshold
HArange <- function(df, thresh, BPA){
    cos_dec <- cos(df$DEC); sin_dec <- sin(df$DEC)
    DF <- data.frame(matrix(rep(NA, 10), nrow=1)); colnames(DF) <- c('Src', 'I', 'P', 'EVPA', 'HA_start1', 'HA_start2', 'HA_end', 'LST_start1', 'LST_start2', 'LST_end')
    #DF$Src <- df$Src; DF$I <- df$I; DF$P <- df$P; DF$EVPA <- df$EVPA
    HA <- seq(-df$ELHA, df$ELHA + HAresolution, by=HAresolution)  # Hour angle above EL limit
    sinHA <- sin(HA); cosHA <- cos(HA)
    PA <- atan2(sinHA, sin_phi*cos_dec/cos_phi - sin_dec*cosHA) + BPA   # Parallactic angle + BandPA
    calDF <- data.frame(HA=HA, PA=PA, XYcorr=df$U* cos(2.0*PA) - df$Q* sin(2.0*PA), et=NA)
    HA_intercepts <- HA[which(diff(sign(calDF$XYcorr)) != 0)]     # hour angles of sign transition
    if(length(HA_intercepts) < 1){ return(na.omit(DF))}
    HA_XY <- calDF[abs(calDF$XYcorr) > thresh,]$HA
    if(length(HA_XY) < 1){ return(na.omit(DF)) }
    #-------- Plot XY vs LST
    XYrange <- range(calDF$XYcorr); XYrange <- c(-0.1*XYrange[1]+1.1*XYrange[2], 1.1*XYrange[1]-0.1*XYrange[2])
    LSTrange <- range((calDF$HA + df$RA)*hourPerRad); LSTrange <- c(-0.1*LSTrange[1]+1.1*LSTrange[2], 1.1*LSTrange[1]-0.1*LSTrange[2])
	plot((calDF$HA + df$RA)*hourPerRad, calDF$XYcorr, type='n', xlab='LST [h]', ylab='XY correlation [Jy]', main=sprintf('%s Band%d as of %s', df$Src, band, as.character(Sys.Date())))
	#grid(nx=NULL, ny=NULL, lty=2, col='gray', lwd=1)
	#abline(h=thresh, lty=2, col='blue'); abline(h=-thresh, lty=2, col='blue')
    polygon( c(LSTrange, rev(LSTrange)), c(-thresh, -thresh, XYrange[2], XYrange[2]), col='#FF000020', border=F)
    polygon( c(LSTrange, rev(LSTrange)), c(thresh, thresh, XYrange[1], XYrange[1]), col='#FF000020', border=F)
    abline(h=0.0, col='gray'); abline(v=hourPerRad* (HA_intercepts + df$RA))
	lines((calDF$HA + df$RA)*hourPerRad, calDF$XYcorr, col='darkgreen', lwd=2)
    for(intercept in HA_intercepts){ text((intercept + df$RA)*hourPerRad, min(calDF$XYcorr), sprintf('%.1fh', (intercept + df$RA)*hourPerRad), pos=4, srt=90) }
    calDF <- calDF[calDF$HA < max(HA_intercepts) - pointingDuration,]          # start time must be before the last intercept
    #-------- HA range for |XY| > thresh
    indexRange <- which(abs(calDF$XYcorr) > thresh)
    for(intercept in sort(HA_intercepts, TRUE)){
        index <- which(calDF$HA[indexRange] < intercept - pointingDuration)
        calDF$et[indexRange[index]] <- intercept + etMargin
    }
    #-------- HA range for |XY| <= thresh12
    indexRange <- which(abs(calDF$XYcorr) <= thresh)
    for(index in indexRange){
        threshCondition    <- which(HA_XY > calDF[index,]$HA)
        interceptCondition <- which(HA_intercepts > calDF[index,]$HA)
        if( length(threshCondition)* length(interceptCondition) > 0){
            calDF[index,]$et <- max(HA_XY[min(threshCondition)], HA_intercepts[min(interceptCondition)]) + etMargin
        }
    }
    calDF <- na.omit(calDF)
    #-------- Summarize LST range into DF
    etList <- sort(unique(calDF$et)); numWindow <- length(etList)
    DF <- data.frame(Src=rep(df$Src, numWindow), I=rep(df$I,numWindow), P=rep(df$P,numWindow), EVPA=rep(df$EVPA,numWindow), HA_start1=rep(NA,numWindow), HA_start2=rep(NA,numWindow), HA_end=etList)
    for(index in 1:nrow(DF)){
        DF[index,]$HA_start1 <- min(calDF[calDF$et == DF[index,]$HA_end,]$HA) - pointingDuration
        DF[index,]$HA_start2 <- max(calDF[calDF$et == DF[index,]$HA_end,]$HA) - pointingDuration
        DF[index,]$HA_end <- max(DF[index,]$HA_end, DF[index,]$HA_start2 + sessionDuration)
    }
	if(numWindow > 1){
        DF <- DF[order(DF$HA_start1),]
        DF[-diff(DF$HA_start2) < 0,]
		numWindow <- nrow(DF)
	}
    DF$LST_start1 <- DF$HA_start1 + df$RA
    DF$LST_start2 <- DF$HA_start2 + df$RA
    DF$LST_end    <- DF$HA_end    + df$RA
    for(row_index in 1:numWindow){
        window_vertical_offset <- 0.1*abs(diff(XYrange))*(row_index-1)
	    lines(c(DF[row_index,]$LST_start1, DF[row_index,]$LST_start2)*hourPerRad, rep(window_vertical_offset,2), lwd=4, col='blue')
	    lines(c(DF[row_index,]$LST_start2, DF[row_index,]$LST_end)*hourPerRad,  rep(window_vertical_offset,2), lwd=0.5, col='blue')
	    lines(c(DF[row_index,]$LST_end, DF[row_index,]$LST_end + 0.25)*hourPerRad,  rep(window_vertical_offset,2), lwd=2, lty=2, col='blue')
	    points(DF[row_index,]$LST_end*hourPerRad, window_vertical_offset, pch=18, cex=2, col='#00000080')
        text(DF[row_index,]$LST_start1*hourPerRad, window_vertical_offset, sprintf('%.1fh', DF[row_index,]$LST_start1*hourPerRad), offset=1, pos=3, col='blue', srt=90)
        text(DF[row_index,]$LST_start2*hourPerRad, window_vertical_offset, sprintf('%.1fh', DF[row_index,]$LST_start2*hourPerRad), offset=1, pos=1, col='blue', srt=-90)
        text(DF[row_index,]$LST_end*hourPerRad, window_vertical_offset,  sprintf('%.1fh', DF[row_index,]$LST_end*hourPerRad), offset=1, pos=1, col='blue', srt=-90)
        text(DF[row_index,]$LST_start1*hourPerRad, window_vertical_offset, 'start', adj=c(-0.2,-1), col='black')
        text(DF[row_index,]$LST_end*hourPerRad, window_vertical_offset,  'end', adj=c(-0.5,-1), col='black')
    }
    return( DF[,c('Src', 'I', 'P', 'EVPA', 'LST_start1', 'LST_start2', 'LST_end')] )
}
#-------- Flag LST range
LSTfrag <- function(df){
    LSTrange <- pi* seq(-1, 1440, by=1)/720.0  # 1-min unit, in [rad]
    LSTfrag <- rep(1, length(LSTrange))
    df <- df[!is.na(df$LST_end),]
    df$LST_start1 <- df$LST_start1 %% (2*pi)
    df$LST_start2 <- df$LST_start2 %% (2*pi)
    for(index in 1:nrow(df)){
        if(df$LST_start1[index] < df$LST_start2[index]){
            LSTindex <- which((LSTrange > df$LST_start1[index]) & (LSTrange < df$LST_start2[index]))
        } else {
            LSTindex <- which((LSTrange > df$LST_start1[index]) | (LSTrange < df$LST_start2[index]))
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
    LST12DF <- na.omit(data.frame(matrix(rep(NA, 7), nrow=1))); names(LST12DF) <- c('Src', 'I', 'P', 'EVPA', 'LST_start1', 'LST_start2', 'LST_end') 
    LST7DF  <- na.omit(data.frame(matrix(rep(NA, 7), nrow=1))); names(LST7DF) <- c('Src', 'I', 'P', 'EVPA', 'LST_start1', 'LST_start2', 'LST_end') 
    #-------- XY-LST plot
    plotDF <- plotLST(na.omit(srcDF[srcDF$P - srcDF$eP > Pthresh12[band], ]), band)
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
    for(index in 1:nrow(srcDF)){
		pngFile <- sprintf('%s-Band%d-PA-12m.png', srcDF$Src[index], band)
		png(pngFile, width=1024, height=768)
        LST12DF <- rbind(LST12DF, HArange(srcDF[index,], Pthresh12[band], BandPA[band]))
		dev.off()
		pngFile <- sprintf('%s-Band%d-PA-7m.png', srcDF$Src[index], band)
		png(pngFile, width=1024, height=768)
        LST7DF <- rbind(LST7DF, HArange(srcDF[index,], Pthresh7[band], BandPA[band]))
		dev.off()
    }
    LSTwindow12 <- LSTfrag(LST12DF)
    LSTwindow7  <- LSTfrag(LST7DF)
	#-------- HTML calibrator table for 12m array
    html.head <- paste("<head>", '<link rel="stylesheet" type="text/css" href="https://www.alma.cl/~skameno/resources/amapola.css" />', "</head>", sep='\n')
    htmlFile <- sprintf('PolCal12m-Band%d.html', band)
    HTMLdf <- LST12DF
    HTMLdf$Src <- paste('<a href="PNG/', sprintf('%s-Band%d',LST12DF$Src, band), '-PA-12m.png" target="_new" > ',  LST12DF$Src,  ' </a>', sep='')
    HTMLdf$EVPA <- RADDEG* HTMLdf$EVPA; HTMLdf$LST_start1 <- hourPerRad* HTMLdf$LST_start1; HTMLdf$LST_start2 <- hourPerRad* HTMLdf$LST_start2; HTMLdf$LST_end <- hourPerRad* HTMLdf$LST_end
    names(HTMLdf) <- c('Source', 'I [Jy]', 'P [Jy]', 'EVPA [deg]', 'LST_start1 [h]', 'LST_start2 [h]', 'LST_end [h]')
    html.table <- paste(print(xtable(HTMLdf, digits=c(0,0,3,3,2,2,2,2)), include.rownames=F, type="html", sanitize.text.function=function(x){x}, htmlFile), collapse="\n")
    CaptionText <- paste("<p>", sprintf('12m-Array Band %d : as of %s', band, as.character(as.Date(Today))),sep='')
    html.body <- paste("<body>", CaptionText, html.table, "</body>")
    write(paste(html.head, html.body, sep='\n'), htmlFile)
	#-------- HTML calibrator table for 7m array
    htmlFile <- sprintf('PolCal7m-Band%d.html', band)
    HTMLdf <- LST7DF
    HTMLdf$Src <- paste('<a href="PNG/', sprintf('%s-Band%d',LST7DF$Src, band), '-PA-7m.png" target="_new" > ',  LST7DF$Src,  ' </a>', sep='')
    HTMLdf$EVPA <- RADDEG* HTMLdf$EVPA; HTMLdf$LST_start1 <- hourPerRad* HTMLdf$LST_start1; HTMLdf$LST_start2 <- hourPerRad* HTMLdf$LST_start2; HTMLdf$LST_end <- hourPerRad* HTMLdf$LST_end
    names(HTMLdf) <- c('Source', 'I [Jy]', 'P [Jy]', 'EVPA [deg]', 'LST_start1 [h]', 'LST_start2 [h]', 'LST_end [h]')
    html.table <- paste(print(xtable(HTMLdf, digits=c(0,0,3,3,2,2,2,2)), include.rownames=F, type="html", sanitize.text.function=function(x){x}, htmlFile), collapse="\n")
    CaptionText <- paste("<p>", sprintf('7m-Array Band %d : as of %s', band, as.character(as.Date(Today))),sep='')
    html.body <- paste("<body>", CaptionText, html.table, "</body>")
    write(paste(html.head, html.body, sep='\n'), htmlFile)
	#CSVdf <- data.frame(Src=LST12DF$Src, I=sprintf(' %7.4f', LST12DF$I), P=sprintf(' %7.4f', LST12DF$P), EVPA=sprintf(' %+7.4f', LST12DF$EVPA), LST_start1=sprintf('    %7.4f', hourPerRad*LST12DF$LST_start1), LST_start2=sprintf('    %7.4f', hourPerRad*LST12DF$LST_start2), LST_end=sprintf(' %7.4f', hourPerRad*LST12DF$LST_end))
	#names(CSVdf) <- c('Src       ', '  I     ', '  P     ', ' EVPA   ', ' LST_start1', ' LST_start2', ' LST_end')
    #cat(sprintf('Band %d 12m array: %d sources\n', band, length(unique(CSVdf$Src))))
    #write.csv(CSVdf, sprintf('PolCal12mBand%d.csv', band), row.names=FALSE, quote=FALSE)
	#CSVdf <- data.frame(Src=LST7DF$Src, I=sprintf(' %7.4f', LST7DF$I), P=sprintf(' %7.4f', LST7DF$P), EVPA=sprintf(' %+7.4f', LST7DF$EVPA), LST_start1=sprintf('    %7.4f', hourPerRad*LST7DF$LST_start1), LST_start2=sprintf('    %7.4f', hourPerRad*LST7DF$LST_start2), LST_end=sprintf(' %7.4f', hourPerRad*LST7DF$LST_end))
	#names(CSVdf) <- c('Src       ', '  I     ', '  P     ', ' EVPA   ', ' LST_start1', ' LST_start2', ' LST_end')
    #cat(sprintf('Band %d  7m array: %d sources\n', band, length(unique(CSVdf$Src))))
    #write.csv(CSVdf, sprintf('PolCal7mBand%d.csv', band), row.names=FALSE, quote=FALSE)
    #-------- Plot 
    uniqueCalibrators <- unique(LST12DF$Src)
	pngFile <- sprintf('Band%d-LST.png', band)
    png(pngFile, width=1536, height=1024)
    par(mar=c(4,8,3,3))
    LSTplot <- barplot(height=rep(NA, length(uniqueCalibrators)), names=uniqueCalibrators, horiz=TRUE, las=1, xlim=c(0,24), xlab='LST to start [h]', main=sprintf('LST windows to start a polarization session in Band %d as of %s', band, Sys.Date()), xaxp=c(0, 24, 24))
    for(index in 1:nrow(LSTwindow7)){
        rect( hourPerRad*LSTwindow7$LSTbegin[index], 0, hourPerRad*LSTwindow7$LSTend[index], max(LSTplot)+1, col='lightgreen', border='transparent')
        text_sd <- sprintf('LST=[%.1fh - %.1fh]', hourPerRad*abs(LSTwindow7$LSTbegin[index]), hourPerRad*LSTwindow7$LSTend[index])
        text(0.5*hourPerRad*(LSTwindow7$LSTbegin[index] + LSTwindow7$LSTend[index]), -0.5, text_sd, cex=2, col='darkgreen')
    }
    grid(nx=24, ny=0, lwd=0.5, lty=1, col='gray')
    #for(calibrator in uniqueCalibrators){
    for(index in seq_along(uniqueCalibrators)){
        calibrator <- uniqueCalibrators[index]
        calDF <- LST12DF[LST12DF$Src == calibrator,]
        segWidth <- max(50*calDF$P[1], 15)
        abline(h=LSTplot[index], col='black', lwd=0.5)
        arrows(hourPerRad* min(calDF$LST_start1), LSTplot[index], hourPerRad* max(calDF$LST_start2), LSTplot[index], length=0, lwd=segWidth, col='blue')
        if( min(calDF$LST_start1) < 0.0){
            arrows(hourPerRad* min(calDF$LST_start1) + 24.0, LSTplot[index], 24.0, LSTplot[index], length=0, lwd=segWidth, col='blue')
            text(hourPerRad*min(calDF$LST_start1) + 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*min(calDF$LST_start1) + 24.0), pos=2) 
        } else {
            text(hourPerRad*min(calDF$LST_start1), LSTplot[index], sprintf('%.1f', hourPerRad*min(calDF$LST_start1)), pos=2) 
        }
        if( max(calDF$LST_start2) > 2.0*pi){
            arrows(0.0, LSTplot[index], hourPerRad* max(calDF$LST_start2) - 24.0, LSTplot[index], length=0, lwd=segWidth, col='blue')
            text(hourPerRad*max(calDF$LST_start2) - 24.0, LSTplot[index], sprintf('%.1f', hourPerRad*max(calDF$LST_start2) - 24.0), pos=4) 
        } else {
            text(hourPerRad*max(calDF$LST_start2), LSTplot[index], sprintf('%.1f', hourPerRad*max(calDF$LST_start2)), pos=4) 
        }
        legend('right',  legend=c('12-m or 7-m Array','Only 12-m Array','7-m Array Window'), col=c('red','blue', 'lightgreen'), lty=1, lwd=c(2,2,8))
        if(calibrator %in% LST7DF$Src ){             # Usable for 7m
            segWidth <- max(30*calDF$P[1], 10)
            calDF <- LST7DF[LST7DF$Src == calibrator,]
            arrows(hourPerRad* min(calDF$LST_start1), LSTplot[index], hourPerRad* max(calDF$LST_start2), LSTplot[index], length=0, lwd=segWidth, col='red')
            if( min(calDF$LST_start1) < 0.0){
                arrows(hourPerRad* min(calDF$LST_start1) + 24.0, LSTplot[index], 24.0, LSTplot[index], length=0, lwd=segWidth, col='red')
            }
            if( max(calDF$LST_start2) > 2.0*pi){
                arrows(0.0, LSTplot[index], hourPerRad* max(calDF$LST_start2) - 24.0, LSTplot[index], length=0, lwd=segWidth, col='red')
            }
        }
    }
    dev.off()
}

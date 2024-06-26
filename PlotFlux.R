library(RColorBrewer)
library(xtable)
library(plotly, warn.conflicts=FALSE)
library(pandoc)
library(htmlwidgets)
library(RCurl)
#-------- Parse arguments
parseArg <- function( args ){
    argList <- list(as.character(Sys.Date()), TRUE, TRUE, TRUE)
    names(argList) <- c('Date', 'Table', 'Source', 'LST')
	if( is.logical(args) == FALSE) return(argList)
    argNum <- length(args)
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-D"){ argList$Date  <- as.character(substring(args[index], 3)) }
        if(substr(args[index], 1,2) == "-T"){ argList$Table <- FALSE }
        if(substr(args[index], 1,2) == "-S"){ argList$Source <- FALSE }
        if(substr(args[index], 1,2) == "-L"){ argList$LST <- FALSE }
    }
    return(argList)
}
#-------- Running input parameters
Arguments <- commandArgs(trailingOnly = TRUE)
argList <- parseArg(Arguments)
Today <- argList$Date
Sys.setenv(TZ="UTC")
ALMA_lat <- -23.029 * pi / 180
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(43.0, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
#-------- Load Flux.Rdata from web
FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
RdataURL    <- "https://www.alma.cl/~skameno/AMAPOLA/"
load("Flux.Rdata")     # Locally generated FLDF
#load(url(paste(RdataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
# (temporal webserver) FluxDataURL <- "https://kamenoseiji.sakura.ne.jp/AMAPOLA/"
# (temporal webserver) load('Flux.Rdata')     # Data frame of FLDF
URL <- "https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R"
Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)

pos <- regexpr("RB",FLDF$File)
FLDF$Band <- as.integer(substr(FLDF$File, pos+3, pos+4))
FLDF$BandPA <- BandPA[FLDF$Band]
FLDF$errU <- FLDF$eP_upper - FLDF$P
FLDF$errL <- FLDF$P - FLDF$eP_lower
FLDF <- FLDF[((FLDF$Band != 3) | (abs(FLDF$Freq - 97.45) < 1.0)),]
#index <- which((FLDF$Band == 3) & (abs(FLDF$Freq - 97.45) > 1.0))
#FLDF <- FLDF[-index,]
#-------- Plot LST
plotLST <- function(DF, band){
    pDF <- data.frame()
	DF <- DF[abs(DF$DEC - ALMA_lat) > 4.0*pi/180.0,]	# zenith avoidance of 4 degree
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
#
#-------- Source List
sourceList <- unique(FLDF$Src)
sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
numSrc <- length(sourceList)
#-------- filter by number of observations
for(src_index in 1:numSrc){
    index <- which(FLDF$Src == sourceList[src_index])
    #if(length(index) < 6){ FLDF <- FLDF[-index,]}
    if(length(index) < 4){ FLDF <- FLDF[-index,]}
}
sourceList <- unique(FLDF$Src)
sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
numSrc <- length(sourceList)
RAList <- (60.0* as.numeric(substring(sourceList, 2,3)) + as.numeric(substring(sourceList, 4,5))) / 720 * pi # RA in [rad]
DecList<- as.numeric(substring(sourceList, 6,8))
DecList<- DecList + sign(DecList)* as.numeric(substring(sourceList, 9,10))/60.0
DecList<- DecList / 180 * pi # DEC in [rad]

#-------- Freq List
bandList <- c(1,3,4,6,7)
FLDF <- FLDF[FLDF$Band %in% bandList,]
numFreq <- length(bandList); freqList <- numeric(numFreq)
for(band_index in 1:length(bandList)){
    freqList[band_index] <- median(FLDF[FLDF$Band == bandList[band_index],]$Freq)
    FLDF[FLDF$Band == bandList[band_index],]$Freq <- freqList[band_index]
}
freqLabel <- sprintf('%.1f GHz', freqList)
#-------- HTML table of source flux 
if(argList$Table){
	for(freq_index in 1:numFreq){
		medI <- medQ <- medU <- eI <- eQ <- eU <- numObs <- rep(NA, numSrc)
    	lastDateIndex <- which.max(as.numeric(FLDF[FLDF$Freq == freqList[freq_index],]$Date))
		for(src_index in 1:numSrc){
			DF <- FLDF[((FLDF$Src == sourceList[src_index]) & (FLDF$Freq == freqList[freq_index]) & (difftime(Today, FLDF$Date, units="days") < 60)) , ]
        	if(nrow(DF) == 0){ next }
			medI[src_index] <- median(DF$I); medQ[src_index] <- median(DF$Q); medU[src_index] <- median(DF$U)
			numObs[src_index] <- length(DF$eI)
			if(numObs[src_index] == 1){
				eI[src_index] <- DF$eI; eQ[src_index] <- DF$eQ; eU[src_index] <- DF$eU
			} else {
				eI[src_index] <- sd(DF$I); eQ[src_index] <- sd(DF$Q); eU[src_index] <- sd(DF$U)
			}
		}
		polDF <- na.omit(data.frame( Src=as.character(sourceList), numObs=numObs, I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi, sdEVPA=90.0*sqrt(eQ^2 + eU^2)/sqrt(medQ^2 + medU^2)/pi ))
		polDF <- polDF[order(polDF$P, decreasing=T),]
		rownames(polDF) <- c(1:nrow(polDF))
		#-------- HTML pol-table
		CaptionText <- paste("<p>", sprintf('Frequency %.1f GHz : 60-day median as of %s / ', freqList[freq_index], as.character(Today)),sep='')
    	CaptionText <- paste(CaptionText, '<a href="', FluxDataURL, FLDF[FLDF$Freq == freqList[freq_index],]$File[lastDateIndex], '" target="_new">', 'Last Observation on ', as.character(FLDF[FLDF$Freq == freqList[freq_index],]$Date[lastDateIndex], tz='UTC'), "</p>", sep='') 
		cat(sprintf('Frequency %.1f GHz : 60-day median as of %s\n', freqList[freq_index], as.character(Today)))
		cat('Source       #obs   I [Jy]   Q [Jy]   U [Jy]    %Pol  EVPA [deg]\n')
		for(index in 1:nrow(polDF)){
			pDF <- polDF[index,]
			cat(sprintf("%10s   (%2d)   %6.2f   %6.2f   %6.2f   %6.1f   %6.1f   %6.1f\n", pDF$Src, pDF$numObs, pDF$I, pDF$Q, pDF$U, pDF$p, pDF$EVPA, pDF$sdEVPA))
		}
		names(polDF) <- c('Source', '#obs', 'I [Jy]', 'sd(I)', 'Q [Jy]', 'sd(Q)', 'U [Jy]', 'sd(U)', 'P [Jy]', 'sd(P)', '%pol', 'EVPA (deg)', 'sd(EVPA)')
		polDF$Source <- paste('<a href="', polDF$Source, '.flux.html" target="_new" >', polDF$Source, ' </a>', sep='')
		htmlFile <- sprintf('Stokes%.0fGHz.html', freqList[freq_index])
		html.head <- paste("<head>", '<link rel="stylesheet" type="text/css" href="https://www.alma.cl/~skameno/resources/amapola.css" />', "</head>", sep='\n')
		html.table <- paste(print(xtable(polDF, digits=c(0,0,0,3,3,3,3,3,3,3,3,1,1,1)), include.rownames=F, type="html", sanitize.text.function=function(x){x}, htmlFile), collapse="\n")
		html.body <- paste("<body>", CaptionText, html.table, "</body>")
		write(paste(html.head, html.body, sep='\n'), htmlFile)
	}
}
#-------- Time-series plots
if(0){
bandColor <- brewer.pal(numFreq, "Dark2")
for(src_index in 1:numSrc){
    rm(DF)
	DF <- FLDF[FLDF$Src == sourceList[src_index],]
	DF$Date <- as.POSIXct(DF$Date, tz='GMT')
	plot_I <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~I, type="scatter", mode="markers", color=paste("I", freqLabel[1]), colors=bandColor, error_y = list(array=~eI, thickness=1, width=0))
    #if(nrow(DF[DF$Band == bandList[2],]) > 0){ plot_I <- add_trace(plot_I, data=DF[DF$Band == bandList[2],], color=paste("I", freqLabel[2]))}
    #if(nrow(DF[DF$Band == bandList[3],]) > 0){ plot_I <- add_trace(plot_I, data=DF[DF$Band == bandList[3],], color=paste("I", freqLabel[3]))}
	for(freq_index in 2:numFreq){
    	if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_I <- add_trace(plot_I, data=DF[DF$Band == bandList[freq_index],], color=paste("I", freqLabel[freq_index]))}
	}
	plot_I <- layout(plot_I, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Stokes I [Jy]',rangemode='tozero', hoverformat='.3f'), title=sourceList[src_index])

	plot_P <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~P, type="scatter", mode="markers", color=paste("P",freqLabel[1]), colors=bandColor, error_y = list(symmetric=FALSE, array=~errU, arrayminus=~errL, thickness=1, width=0))
    #if(nrow(DF[DF$Band == bandList[2],]) > 0){ plot_P <- add_trace(plot_P, data=DF[DF$Band == bandList[2],], color=paste("P",freqLabel[2]))}
    #if(nrow(DF[DF$Band == bandList[3],]) > 0){ plot_P <- add_trace(plot_P, data=DF[DF$Band == bandList[3],], color=paste("P",freqLabel[3]))}
	for(freq_index in 2:numFreq){
    	if(nrow(DF[DF$Band == bandList[freq_index],]) > 0){ plot_P <- add_trace(plot_P, data=DF[DF$Band == bandList[freq_index],], color=paste("P",freqLabel[freq_index]))}
	}
	plot_P <- layout(plot_P, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Polarized Flux [Jy]',rangemode='tozero', hoverformat='.3f'))

	plot_A <- plot_ly(data=DF[DF$Band == bandList[1],], x=~Date, y=~EVPA*180/pi, type="scatter", mode="markers", color=paste("EVPA",freqLabel[1]), colors=bandColor, error_y = list(array=~eEVPA*180/pi, thickness=1, width=0))
    #if(nrow(DF[DF$Band == bandList[2],]) > 0){ plot_A <- add_trace(plot_A, data=DF[DF$Band == bandList[2],], color=paste("EVPA",freqLabel[2]))}
    #if(nrow(DF[DF$Band == bandList[3],]) > 0){ plot_A <- add_trace(plot_A, data=DF[DF$Band == bandList[3],], color=paste("EVPA",freqLabel[3]))}
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
#-------- Source 60-day statistics
I100 <- Q100 <- U100 <- numeric(numSrc)
spixI <- spixP <- rep(NA, numSrc)
for(src_index in 1:numSrc){
	rm(DF)
	DF <- FLDF[((FLDF$Src == sourceList[src_index]) & (difftime(Today, FLDF$Date, units="days") < 60)) , ]
    if(nrow(DF) < 3){ next }
	if(min(abs(as.numeric(difftime(DF$Date, Today)))) > 30){ next }
	bands <- unique(DF$Band)
	numBand <- length(bands)
	predI <- predQ <- predU <- eI <- eQ <- eU <- numObs <- freq <- numeric(numBand)
    if(numBand > 1){
	    for(band_index in 1:numBand){
	    	index <- which(DF$Band == bands[band_index])
            numObs[band_index] <- length(index)
            deltaDay <- as.numeric(difftime(DF[index,]$Date, Today))
            BandTimeFit <- TRUE
	    	if( numObs[band_index] <= 2 ){ BandTimeFit <- FALSE }
            if( BandTimeFit ){
	    	    if(sd(deltaDay) < min(abs(deltaDay)) ){ BandTimeFit <- FALSE }
            }
            if( BandTimeFit ){
	    		deltaDay <- as.numeric(difftime(DF[index,]$Date, Today))
	    		fit <- lm(DF$I[index] ~ deltaDay, weights=1/DF$eI[index]^2/abs(deltaDay + 1))
	    		predI[band_index] <- summary(fit)$coefficients[1,'Estimate']
	    		eI[band_index] <- summary(fit)$coefficients[1,'Std. Error']
	    		fit <- lm(DF$Q[index] ~ deltaDay, weights=1/DF$eQ[index]^2/abs(deltaDay + 1))
	    		predQ[band_index] <- summary(fit)$coefficients[1,'Estimate']
	    		eQ[band_index] <- summary(fit)$coefficients[1,'Std. Error']	
	    		fit <- lm(DF$U[index] ~ deltaDay, weights=1/DF$eU[index]^2/abs(deltaDay + 1))
	    		predU[band_index] <- summary(fit)$coefficients[1,'Estimate']
	    		eU[band_index] <- summary(fit)$coefficients[1,'Std. Error']	
            } else {
	    		predI[band_index] <- mean(DF$I[index]); predQ[band_index] <- mean(DF$Q[index]); predU[band_index] <- mean(DF$U[index])
	    		eI[band_index] <- median(DF$eI[index]); eQ[band_index] <- median(DF$eQ[index]); eU[band_index] <- median(DF$eU[index])
            }
	    	freq[band_index] <- median(DF$Freq[index])
	    }
    	fit <- lm( log(predI) ~ log(freq/100), weights=1.0/eI^2 ); spixI[src_index] <- coef(fit)[2]; I100[src_index] <- exp(coef(fit)[1])
    	fit <- lm(0.5*log(predQ^2 + predU^2) ~ log(freq/100), weights=1.0/eI^2 ); spixP[src_index] <- coef(fit)[2]
    	f100 <- (freq/100.0)^spixP[src_index]
    	fit <- lm(predQ ~ f100+0, weights=1.0/eQ^2); Q100[src_index] <- coef(fit)[1]
    	fit <- lm(predU ~ f100+0, weights=1.0/eU^2); U100[src_index] <- coef(fit)[1]
    } else {
        I100[src_index] <- median(DF$I)
        Q100[src_index] <- median(DF$Q)
        U100[src_index] <- median(DF$U)
        spixI[src_index] <- -0.7
        spixP[src_index] <- -0.7
    }
}
srcDF <- na.omit(data.frame(Src=sourceList, RA=RAList, DEC=DecList, I100=I100, Q100=Q100, U100=U100, spixI=spixI, spixP=spixP))
write.csv(srcDF[(abs(srcDF$DEC - ALMA_lat) > 4.0*pi/180.0),], 'PolQuery.CSV', row.names=FALSE)

for(band in c(1,3,4,5,6,7,8,9)){
	plotDF <- plotLST(srcDF, band)
	pLST <- plot_ly(data=plotDF, x = ~LST, y = ~XYcorr, type = 'scatter', mode = 'lines', color=~Src, hoverinfo='text', text=~paste(Src, 'EL=',floor(EL)))
	pLST <- layout(pLST, xaxis=list(showgrid=T, title='LST', nticks=24), yaxis=list(showgrid=T, title='XY correlation [Jy]',rangemode='tozero'), title=sprintf('Band-%d Pol-Calibrator Coverage as of %s (30-day statistics)', band, as.character(Today)))
	htmlFile <- sprintf("Band%d_LSTplot.html", band)
	htmlwidgets::saveWidget(pLST, htmlFile, selfcontained=FALSE)
	rm(plotDF)
}

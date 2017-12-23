library(xtable)
library(plotly)
library(htmlwidgets)
library(RCurl)
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)

#-------- Load Flux.Rdata from web
load(url("http://www.alma.cl/~skameno/Grid/Stokes/Flux.Rdata"))     # Data frame of FLDF
URL <- "https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R"
Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)

pos <- regexpr("RB",FLDF$File)
FLDF$Band <- as.integer(substr(FLDF$File, pos+3, pos+4))
FLDF$BandPA <- BandPA[FLDF$Band]

#-------- Today
Today <- Sys.Date()

plotLST <- function(DF, band){
	ALMA_lat <- -23.029 * pi / 180
	sourceNum <- length(DF$Src)
	LST <- seq(0.0, 2.0*pi, 0.01)
	plot(LST, LST, type='n', xlim=c(0,24), ylim=c(0,0.5), xlab='LST [h]', ylab='XY response [Jy]', main=sprintf('Band %d', band))
	for(source_index in 1:sourceNum){
		HA <- LST - DF$RA[source_index]
		AZEL <- ha2azel( HA, ALMA_lat, DF$Dec[source_index] )
		AZEL$pa <- AZEL$pa + BandPA[band]
		CS <- cos(2.0* AZEL$pa)
		SN <- sin(2.0* AZEL$pa)
		AZEL$XYcorr <- polDF$U[source_index]*CS - polDF$Q[source_index]*SN
		ELrange <- which(AZEL$el > pi/7.5)	# EL > 24 deg.
		if( max(abs(AZEL$XYcorr[ELrange])) > 0.1 ){
			points(LST[ELrange]*12/pi, abs(AZEL$XYcorr[ELrange]), pch=20, cex=0.5, col=source_index)
			LSTpeak <- LST[ELrange[which.max(abs(AZEL$XYcorr[ELrange]))]]
			text(LSTpeak*12/pi, 0.05, DF$Src[source_index], col=source_index)
		}
	}
}
	
#
#-------- Source List
sourceList <- unique(FLDF$Src)
sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
numSrc <- length(sourceList)
RAList <- (60.0* as.numeric(substring(sourceList, 2,3)) + as.numeric(substring(sourceList, 4,5))) / 720 * pi # RA in [rad]
DecList<- as.numeric(substring(sourceList, 6,8))
DecList<- DecList + sign(DecList)* as.numeric(substring(sourceList, 9,10))/60.0
DecList<- DecList / 180 * pi # DEC in [rad]

#-------- Freq List
freqList <- unique(FLDF$Freq)
bandList <- unique(FLDF$Band)
numFreq <- length(freqList)

#-------- HTML table of source flux 
for(freq_index in 1:numFreq){
	medI <- medQ <- medU <- eI <- eQ <- eU <- numObs <- numeric(0)
	for(src_index in 1:numSrc){
		DF <- FLDF[((FLDF$Src == sourceList[src_index]) & (FLDF$Freq == freqList[freq_index]) & (difftime(Today, FLDF$Date, units="days") < 60)) , ]
		medI[src_index] <- median(DF$I); medQ[src_index] <- median(DF$Q); medU[src_index] <- median(DF$U)
		numObs[src_index] <- length(DF$eI)
		if(numObs[src_index] == 1){
			eI[src_index] <- DF$eI; eQ[src_index] <- DF$eQ; eU[src_index] <- DF$eU
		} else {
			eI[src_index] <- sd(DF$eI); eQ[src_index] <- sd(DF$eQ); eU[src_index] <- sd(DF$eU)
		}
	}
	#polDF <- data.frame( Src=as.character(sourceList), RA=RAList, Dec=DecList, numObs=numObs, I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi )
	polDF <- data.frame( Src=as.character(sourceList), numObs=numObs, I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi )
	polDF <- polDF[order(polDF$P, decreasing=T),]
	rownames(polDF) <- c(1:nrow(polDF))
	#-------- Plot Pol-LST
	#plotLST(polDF[!is.na(polDF$I),], bandList[freq_index])
	CaptionText <- paste("<p>", sprintf('Frequency %.1f GHz : 60-day median as of %s\n', freqList[freq_index], as.character(Today)), "</p>", sep='\n')
	cat(sprintf('Frequency %.1f GHz : 60-day median as of %s\n', freqList[freq_index], as.character(Today)))
	cat('Source       #obs   I [Jy]   Q [Jy]   U [Jy]    %Pol  EVPA [deg]\n')
	for(index in 1:numSrc){
		pDF <- polDF[index,]
		cat(sprintf("%10s   (%2d)   %6.2f   %6.2f   %6.2f   %6.1f    %6.1f\n", pDF$Src, pDF$numObs, pDF$I, pDF$Q, pDF$U, pDF$p, pDF$EVPA))
	}
	names(polDF) <- c('Source', '#obs', 'I (Jy)', 'err I', 'Q (Jy)', 'err Q', 'U (Jy)', 'err U', 'P (Jy)', 'err P', '%pol', 'EVPA (deg)')
	polDF$Source <- paste('<a href="', polDF$Source, '.flux.html" target="_new" >', polDF$Source, ' </a>', sep='')
	htmlFile <- sprintf('Stokes%.0fGHz.html', freqList[freq_index])
	html.head <- paste("<head>", '<link rel="stylesheet" type="text/css" href="http://www.alma.cl/~skameno/resources/amapola.css" />', "</head>", sep='\n')
	html.table <- paste(print(xtable(polDF, digits=c(0,0,0,3,3,3,3,3,3,3,3,1,1)), include.rownames=F, type="html", sanitize.text.function=function(x){x}, htmlFile), collapse="\n")
	html.body <- paste("<body>", CaptionText, html.table, "</body>")
	write(paste(html.head, html.body, sep='\n'), htmlFile)
}

#-------- Time-series plots
for(src_index in 1:numSrc){
	DF <- FLDF[FLDF$Src == sourceList[src_index],]
	DF$Freq <- paste(as.character(DF$Freq), "GHz")
	dataRange <- matrix(c(1.1, -0.1, -0.1, 1.1), nrow=2) %*% c(min(DF$Date), max(DF$Date))
	plot_I <- plot_ly(DF, x=~Date, y=~I, type="scatter", mode="markers", error_y = ~list(array=eI, thickness=1, width=0), color=~Freq, showlegend=F)
	plot_I <- layout(plot_I, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Stokes I [Jy]',rangemode='tozero'), title=sourceList[src_index])
	plot_P <- plot_ly(DF, x=~Date, y=~P, type="scatter", mode="markers", error_y = ~list(array=eP, thickness=1, width=0), color=~Freq, showlegend=F)
	plot_P <- layout(plot_P, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='Polarized Flux [Jy]',rangemode='tozero'))
	plot_A <- plot_ly(DF, x=~Date, y=~EVPA*180/pi, type="scatter", mode="markers", error_y = ~list(array=eEVPA*180/pi, thickness=1, width=0), color=~Freq, showlegend=T)
	plot_A <- layout(plot_A, xaxis=list(showgrid=T, title='Date', range=c(min(DF$Date)-86400, max(DF$Date)+86400)), yaxis=list(showgrid=T, title='EVPA [deg]',range=c(-91,91)))
	allPlot <- subplot(plot_I, plot_P, plot_A, nrows=3, shareX=F, titleY=T)
	htmlFile <- sprintf("%s.flux.html", sourceList[src_index])
	htmlwidgets::saveWidget(allPlot, htmlFile)
}

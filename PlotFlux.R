library(xtable)

#-------- Load Flux.Rdata from web
load(url("http://www.alma.cl/~skameno/Grid/Stokes/Flux.Rdata"))     # Data frame of FLDF

#-------- Today
Today <- Sys.Date()

#BandName <- sprintf('Band-%d', as.numeric(strsplit(FLDF$File[1], '[-|.|_]')[[1]][8]))
#-------- Source List
sourceList <- unique(FLDF$Src)
sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
numSrc <- length(sourceList)
#-------- Freq List
freqList <- unique(FLDF$Freq)
numFreq <- length(freqList)

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
	polDF <- data.frame( Src=as.character(sourceList), numObs=numObs, I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi )
	polDF <- polDF[order(polDF$P, decreasing=T),]
	rownames(polDF) <- c(1:nrow(polDF))
	cat(sprintf('Frequency %.1f GHz : 60-day median by %s\n', freqList[freq_index], as.character(Today)))
	cat('Source       #obs   I [Jy]   Q [Jy]   U [Jy]    %Pol  EVPA [deg]\n')
	for(index in 1:numSrc){
		pDF <- polDF[index,]
		cat(sprintf("%10s   (%2d)   %6.2f   %6.2f   %6.2f   %6.1f    %6.1f\n", pDF$Src, pDF$numObs, pDF$I, pDF$Q, pDF$U, pDF$p, pDF$EVPA))
	}
}

if(0){
	
#for(src_index in 1:numSrc){
#	cat(sprintf("%10s  %5.1f  %6.3f  %6.3f\n", polDF$Src[src_index], polDF$I[src_index], polDF$Q[src_index]/polDF$I[src_index], polDF$U[src_index]/polDF#$I[src_index]))
#}

# print(xtable(polDF, digits=3), include.rownames=F, type = "html")

pdf(sprintf('Flux-%s.pdf', BandName), width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
cat('Source     I [Jy]   Q [Jy]   U [Jy]   %Pol   EVPA [deg]\n')
for(source in sourceList){
	DF <- FLDF[FLDF$Src == source,]
	pDF <- polDF[polDF$Src == source,]
	cat(sprintf("%10s %6.2f %6.2f %6.2f %6.2f %6.2f\n", pDF$Src, pDF$I, pDF$Q, pDF$U, pDF$p, pDF$EVPA))
	#-------- Plot Stokes I
	plot(DF$Date, DF$I, type='n', xlab='Date', ylab='Stokes I [Jy]', main='Total Flux Density', ylim=c(0, 1.2*max(DF$I)))
	color_vector <- rep("#000000FF", length(DF$Date))
	color_vector[which(DF$flagNum > 0)] <- "#0000FF3F"
	arrows(as.numeric(DF$Date), DF$I - DF$eI, as.numeric(DF$Date), DF$I + DF$eI, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$I, pch=20, col=color_vector)
	
	#-------- Plot polarized flux
	plot(DF$Date, DF$P, type='n', xlab='Date', ylab='Polarized Flux [Jy]', main='Polarized Flux Density', ylim=c(0, 1.2* max(DF$P)))
	arrows(as.numeric(DF$Date), DF$P - DF$eP, as.numeric(DF$Date), DF$P + DF$eP, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$P, pch=20, col=color_vector)
	
	#-------- Plot EVPA
	plot(DF$Date, DF$EVPA, type='n', xlab='Date', ylab='EVPA [deg]', main='Polarization Angle', ylim=c(-90,90))
	abline(h=90, lwd=0.1); abline(h=-90, lwd=0.1)
	arrows(as.numeric(DF$Date), 180*(DF$EVPA - DF$eEVPA)/pi, as.numeric(DF$Date), 180*(DF$EVPA + DF$eEVPA)/pi, angle=90, length=0.0, col=color_vector)
	points(DF$Date, DF$EVPA*180/pi, pch=20, col=color_vector)
	mtext(side = 3, line=1, outer=T, text = sprintf('%s %s', source, BandName), cex=2)
}
dev.off()
par(par.old)
}

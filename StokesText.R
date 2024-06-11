library(parallel)   # multicore parallelization
library(doParallel) # multicore parallelization
library(VGAM)       # for Rice distribution
Sys.setenv(TZ="UTC")
#sysIerr <- 0.005       # temporal Stokes I systematic error
#sysPerr <- 0.003       # temporal polarization systematic error
minAntNum <- 5		   # Minimum number of antennas
#numSubBand = c(3,3,3,1,1,1,1,1,1,1) # Number of sub-bands for Band[1-10]
standardFreq <- list(40.0, 80.0, c(91.5,97.5,103.5), 154.9, 183.0, 233.0, 343.4, 410.2)
#-------- Functions
#-------- Get band name
getBand <- function(fileName){
    bandPointer <- regexpr("RB_[0-10]", fileName)
    return(as.integer(substr(fileName, bandPointer+3, bandPointer+4)))
}
srcDfFilter  <- function(source){ return(FLDF[FLDF$Src == source,])}
scanDfFilter <- function(scan){ return(FLDF[FLDF$Date == scan,])}

#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
predStokes <- function(df){
    bandID   <- getBand(df$File[1])
    fitI <- lm(formula=I ~ Freq, data=df, weight=1.0/eI^2)
    fitQ <- lm(formula=Q ~ Freq, data=df, weight=1.0/eQ^2)
    fitU <- lm(formula=U ~ Freq, data=df, weight=1.0/eU^2)
    fitV <- lm(formula=V ~ Freq, data=df, weight=1.0/eV^2)
    newDF <- data.frame(Src=df$Src[1], Freq = standardFreq[[bandID]], EL=df$EL[1])
    pred <- as.numeric(predict(fitI, newDF, interval='confidence', level=0.67)); newDF$I <- matrix(pred, ncol=3)[,1]; newDF$eI <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitQ, newDF, interval='confidence', level=0.67)); newDF$Q <- matrix(pred, ncol=3)[,1]; newDF$eQ <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitU, newDF, interval='confidence', level=0.67)); newDF$U <- matrix(pred, ncol=3)[,1]; newDF$eU <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitV, newDF, interval='confidence', level=0.67)); newDF$V <- matrix(pred, ncol=3)[,1]; newDF$eV <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    return(newDF)
}
load('Flux.Rdata')
FLDF$Band <- getBand(FLDF$File)
recNum <- nrow(FLDF)
index <- order(FLDF$Date, FLDF$Src, FLDF$Band, FLDF$Freq)
diffRange1 <- 1:(recNum-1)
diffRange2 <- 2:recNum
recBorder <- which( (FLDF[index,]$Src[diffRange1] != FLDF[index,]$Src[diffRange2]) | diff(FLDF[index,]$Date) > 0 | diff(FLDF[index,]$Band) > 0 )




recBorder <- which((FLDF[index, ]$Src[1:(length(index)-1)] != FLDF[index, ]$Src[2:length(index)]) || (FLDF[index, ]$Date[1:(length(index)-1)] != FLDF[index, ]$Date[2:length(index)]))










#srcDfList <- mclapply(unique(FLDF$Src), function(source){ return(FLDF[FLDF$Src == source,]) })
srcDfList <- mclapply(unique(FLDF$Src), srcDfFilter)
#scanDFList <- mclapply(unique(FLDF$Date), scanDfFilter)

cores <- detectCores()
cl <- makcCluster(cores)
registerDoParallel(cl)

srcNum <- length(srcDfList)
result <- foreach(src_index = 1:srcNum) %dopar% {
    cat(sprintf('%d : %s %d\n', src_index, srcDfList[[src_index]]$Src[1], length(unique(srcDfList[[src_index]]$Date))))
    scanDfList <- mclapply(unique(srcDfList[[src_index]]$Date), scanDfFilter)
}

for(src_index in 1:srcNum){
    srcDF <- srcDfList[[src_index]]
    cat(sprintf('%d : %s %d\n', src_index, srcDF$Src[1], length(unique(srcDF$Date))))
    scanDfList <- mclapply(unique(srcDF$Date), scanDfFilter)
}
    #tempList <- mclapply( scanDfList, predStokes )
    

for(srcDF in srcDfList){
    cat(srcDF$Src[1])
    scanList <- unique(srcDF$Date)
    scanDfList <- mclapply(scanList, function(scan){ return(srcDF[srcDF$Date == scan,]) })
    tempList <- mclapply( scanDfList, predStokes )
}



scanDFList <- mclapply(unique(FLDF$Date), function(scan){ return(FLDF[FLDF$Date == scan,])})



tempList <- mclapply( srcDFList, function(srcDF){
    scanDFList <- mclapply( unique(srcDF$Date), function(scan){ return(srcDF[srcDF$Date == scan,])})
    return(mclapply(scanDfList, predStokes ))
})
    



tempList <- mclapply( scanDfList, predStokes )
    

return(srcD)



srcDfList <- list()



for(source in sourceList){ srcDfList <- c(srcDfList, list(FLDF[FLDF$Src == source,]))}



srcDF <- srcDfList[[1]]
UniqueDate <- unique(srcDF$Date)
scanDfList <- list()
for(scan in UniqueDate){
    scanDfList <- c(scanDfList, list(srcDF[srcDF$Date == scan,]))
}
tempList <- mclapply( scanDfList, predStokes )


        



FLDF$P <- sqrt(FLDF$Q^2 + FLDF$U^2)
sigmaSQ <- sqrt(FLDF$eQ * FLDF$eU + (FLDF$I* sysPerr)^2)
FLDF$eP_lower <- qrice(0.15, sigmaSQ, FLDF$P)
FLDF[FLDF$P < sigmaSQ,]$eP_lower <- 0.0
FLDF$eP_upper <- qrice(0.85, sigmaSQ, FLDF$P)
FLDF$EVPA <- 0.5* atan2(FLDF$U, FLDF$Q)
FLDF$eEVPA <- 0.5* sqrt(FLDF$Q^2 * FLDF$eU^2 + FLDF$U^2 * FLDF$eQ^2) / (FLDF$P)^2
#---- Output to text data


TextDF <- FLDF[order(FLDF$Date),]
index <- which(abs(TextDF$Freq - 97.45) < 1.0)
TextDF <- TextDF[-index,]
TextDF$Src <- sprintf('%10s ', TextDF$Src)
write.table(format(TextDF, digits=4), 'amapola.txt', sep='\t', quote=F, col.names=T, row.names=F)
#
#-------- Source List
sourceList <- unique(FLDF$Src)
numSrc <- length(sourceList)

medI <- medQ <- medU <- eI <- eQ <- eU <- numeric(0)
for(src_index in 1:numSrc){
	DF <- FLDF[FLDF$Src == sourceList[src_index],]
	medI[src_index] <- median(DF$I); medQ[src_index] <- median(DF$Q); medU[src_index] <- median(DF$U)
	if(length(DF$eI) == 1){
		eI[src_index] <- DF$eI; eQ[src_index] <- DF$eQ; eU[src_index] <- DF$eU
	} else {
		eI[src_index] <- sd(DF$eI); eQ[src_index] <- sd(DF$eQ); eU[src_index] <- sd(DF$eU)
	}
}
polDF <- data.frame( Src=as.character(sourceList), I=medI, eI = eI, Q=medQ, eQ = eQ, U=medU, eU = eU, P=sqrt(medQ^2 + medU^2), eP=sqrt(eQ^2 + eU^2)/medI, p=100.0*sqrt(medQ^2 + medU^2)/medI, EVPA=90.0*atan2(medU, medQ)/pi )
polDF <- polDF[order(polDF$P, decreasing=T),]
rownames(polDF) <- c(1:nrow(polDF))
save(polDF, file='Pol.Rdata')

sourceList <- sort(polDF$Src[grep('^J[0-9]', polDF$Src)])
#-------- Plot time-seried flux densities
pdf('AMAPOLA.pdf', width=8, height=11)
par.old <- par(no.readonly=TRUE)
par(mfrow=c(3,1), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
cat('Source     I [Jy]   Q [Jy]   U [Jy]   %Pol   EVPA [deg]\n')
labels <- c("B1",        "B3",        "B4",        "B6",        "B7",        "> B7")
colors <- c("#FF000000", "#E020003F", "#8080003F", "#00C0803F", "#0000FF3F", "#FFFFFF")
threshU <-c(41.0,        98.0,        200.0,       280.0,       380.0,       700.0)
threshL <-c(39.0,        97.0,        120.0,       200.0,       280.0,       380.0)
plotXrange <- range(FLDF$Date)
for(sourceName in sourceList){
	DF <- FLDF[FLDF$Src == sourceName,]
	pDF <- polDF[polDF$Src == sourceName,]
	cat(sprintf("%10s %6.2f %6.2f %6.2f %6.2f %6.2f\n", pDF$Src, pDF$I, pDF$Q, pDF$U, pDF$p, pDF$EVPA))
    #-------- Coloring by frequency
    DF$color_vector <- rep("#000000FF", nrow(DF))
    for(color_index in 1:length(labels)){
        band_index <- which( (DF$Freq > threshL[color_index]) & (DF$Freq < threshU[color_index]) )
        DF$color_vector[band_index] <- colors[color_index]
    }
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

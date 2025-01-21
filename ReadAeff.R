source('~/ALMAR/ReadAeffLib.R')
#-------- Start program
#Arguments <- commandArgs(trailingOnly = T)
#fileList <- parseArg(Arguments)
fileList <- parseArg('fileList')
#AeDF <- data.frame(matrix(rep(NA, length(FMT)), nrow=1))[numeric(0),]; colnames(AeDF) <- FMT
DFList <- mclapply(fileList, Log2Aeff, mc.cores=numCore)
AeDF <- do.call("rbind", DFList)
save(AeDF, file='AeDF.Rdata')
AeDF <- AeDF[complete.cases(AeDF$AeX),]
AeDF <- AeDF[complete.cases(AeDF$AeY),]
AeDF <- AeDF[((AeDF$AeX > 25.0) & (AeDF$AeX < 100.0) & (AeDF$AeY > 25.0) & (AeDF$AeY < 100.0)),]
AeDF <- AeDF[complete.cases(AeDF$Date),]
AeDF <- na.omit(AeDF)
AeDF <- AeDF[order(AeDF$Date), ]
#save(AeDF, file='AeDF.Rdata')

DFList <- mclapply(fileList, Log2Dterm, mc.cores=numCore)
DtermDF <- do.call("rbind", DFList)
DtermDF <- DtermDF[complete.cases(DtermDF$Date),]
DtermDF <- na.omit(DtermDF)
DtermDF <- DtermDF[order(DtermDF$Date), ]
DtermDF$Band <- getBand(DtermDF$File)
save(DtermDF, file='Dterm.Rdata')

#-------- Ae table
BandList <- c(1, 3, 4, 6, 7)
pcolors <- c(alpha('orange', 0.5), alpha('purple', 0.5))
lcolors <- c('orange', 'purple')
refTime <- max(AeDF$Date)
MonthSec <- 2629744
refPeriod <- seq(as.numeric(difftime(min(AeDF$Date), refTime, units='sec')), MonthSec, by=MonthSec)
for(Band in BandList){
    BandAeDF <- AeCorrect(AeDF, Band, 10)
    antList <- sort(as.character(unique(BandAeDF$Ant)))
    bandAeDF <- data.frame(Date = c('mean', 'sd', as.character(as.Date(refTime + refPeriod))))
    for(ant in antList){
        BandAntAeDF <- BandAeDF[BandAeDF$Ant == ant,]
        Ae  <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntAeDF$Date, refTime, units='sec')), Value=BandAntAeDF$Ae), refPeriod)
        XYR <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntAeDF$Date, refTime, units='sec')), Value=BandAntAeDF$AeR), refPeriod)
        AeX <- Ae$Value / sqrt(XYR$Value)
        AeY <- Ae$Value * sqrt(XYR$Value) 
        bandAeDF[[paste(ant , '-X', sep='')]] <- c(mean(BandAntAeDF$AeC / sqrt(BandAntAeDF$AeR)), sd(BandAntAeDF$AeC / sqrt(BandAntAeDF$AeR)), AeX)
        bandAeDF[[paste(ant , '-Y', sep='')]] <- c(mean(BandAntAeDF$AeC * sqrt(BandAntAeDF$AeR)), sd(BandAntAeDF$AeC * sqrt(BandAntAeDF$AeR)), AeY)
        pdf(sprintf('Ae-%s-B%d.pdf', ant, Band))
        plot(BandAeDF$Date, BandAeDF$Ae, type='n', ylim=c(0, 100), xlab='Date', ylab='Aeff (%)', main=sprintf('%s Band%d', ant, Band))
        lines(Ae$Date + refTime, Ae$Value / sqrt(XYR$Value), col=lcolors[1], lwd=2)
        lines(Ae$Date + refTime, Ae$Value * sqrt(XYR$Value), col=lcolors[2], lwd=2)
        points(BandAntAeDF$Date, BandAntAeDF$Ae / sqrt(BandAntAeDF$AeR), pch=20, cex=0.5, col=pcolors[1])
        points(BandAntAeDF$Date, BandAntAeDF$Ae * sqrt(BandAntAeDF$AeR), pch=20, cex=0.5, col=pcolors[2])
        legend("bottomleft", legend=c('Pol-X', 'Pol-Y'), col=lcolors, pch=rep(20, 2), lty=rep(1,2))
        dev.off()
    }
    write.table(format(bandAeDF, digits=4), file=sprintf('AeB%d.table', Band), quote=F, row.names=F)
}
#-------- Dterm table
refTime <- max(DtermDF$Date)
refPeriod <- seq(as.numeric(difftime(min(DtermDF$Date), refTime, units='sec')), MonthSec, by=MonthSec)
for(Band in BandList){
    BandDdf <- DtermDF[DtermDF$Band == Band,]
    antList <- sort(as.character(unique(BandDdf$Ant)))
    for(ant in antList){
        BandAntdDF <- BandDdf[BandDdf$Ant == ant,]
        bandAntDdf <- data.frame(Date = as.Date(refTime + refPeriod))
        if(nrow(BandAntdDF) < 3){
            for(BB in c(1,2,3,4)){
                for(pol in c('x','y')){
                    bandAntDdf[[sprintf('%s-BB%d-D%s', ant, BB, pol)]] <- (0.0 + 0.0i)
                }
            }
        } else {
            pdf(sprintf('Dterm.B%d.%s.pdf', Band, ant), width=8, height=11)
            par.old <- par(no.readonly=TRUE)
            par(mfrow=c(4,2), oma=c(0, 0, 4, 0), mar=c(4,4,4,4))
            column_index <- 1
            for(BB in c(1,2,3,4)){
                for(pol in c('x','y')){
                    colName <- sprintf('D%s%d', pol, BB)
                    ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                    bandAntDdf[[sprintf('%s-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
                    #---- plot
                    column_index <- column_index + 1
                    plot(as.Date(BandAntdDF$Date), Re(BandAntdDF[[column_index]]), pch=21, cex=0.2, ylim=c(-0.1, 0.1), xlab='Date', ylab='D-term', main=colnames(BandAntdDF[column_index]), col=pcolors[1] )
                    points(as.Date(BandAntdDF$Date), Im(BandAntdDF[[column_index]]), pch=21, cex=0.2, col=pcolors[2] )
                    lines(as.Date(bandAntDdf[[1]]), Re(bandAntDdf[[column_index]]), col=lcolors[1], lwd=2)
                    lines(as.Date(bandAntDdf[[1]]), Im(bandAntDdf[[column_index]]), col=lcolors[2], lwd=2)
                    legend("bottomleft", legend=c('Real', 'Imag'), col=lcolors, pch=rep(20, 2), lty=rep(1,2))
                }
            }
            mtext(side = 3, line=1, outer=T, text = sprintf('D-term %s Band-%d', ant, Band), cex=2)
            par(par.old)
            dev.off()
        }
        write.table(format(bandAntDdf, digits=6), file=sprintf('DtermB%d.%s.table', Band, ant), quote=F, row.names=F)
    }
    #---- Dummy antenna
    for(ant in c('CM', 'PM', 'DA', 'DV')){
        BandAntdDF <- subset(BandDdf, grepl(ant,  BandDdf$Ant))
        bandAntDdf <- data.frame(Date = as.Date(refTime + refPeriod))
        for(BB in c(1,2,3,4)){
            for(pol in c('x','y')){
                colName <- sprintf('D%s%d', pol, BB)
                ReD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Re(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                ImD <- SPL_period(data.frame(relSec=as.numeric(difftime(BandAntdDF$Date, refTime, units='sec')), Value=Im(BandAntdDF[[colName]])), refPeriod, 1.0/(abs(BandAntdDF[[colName]])+0.001))
                bandAntDdf[[sprintf('%s00-BB%d-D%s', ant, BB, pol)]] <- ReD$Value + (0 + 1i)* ImD$Value
            }
        }
        write.table(format(bandAntDdf, digits=6), file=sprintf('DtermB%d.%s00.table', Band, ant), quote=F, row.names=F)
    }
}

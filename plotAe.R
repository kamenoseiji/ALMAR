library(RCurl)
FuncList <- c('date')
for(index in 1:length(FuncList)){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }
setwd('.')
MJDepoch <- doy2mjd(1970,1)
#-------- Execution
AeDF <- read.table('Aeff.log', header=T)
AeDF$Date <- as.POSIXct(AeDF$mjdSec - MJDepoch*86400, tz="GMT", origin="1970-01-01")
AeDF$relTime <- AeDF$mjdSec - median(AeDF$mjdSec)
AeDF$DayTime <- 24* (AeDF$mjdSec %% 1)
AeDF$secZ <- 1.0 / sin(AeDF$EL)
#-------- Band Separation
B3 <- AeDF[AeDF$Band == 3,]
B6 <- AeDF[AeDF$Band == 6,]
B7 <- AeDF[AeDF$Band == 7,]
#-------- Antenna Separation
CMB3 <- na.omit(B3[grep('CM', B3$Ant),])
CMB6 <- na.omit(B6[grep('CM', B3$Ant),])
CMB7 <- na.omit(B7[grep('CM', B7$Ant),])
colores = c('red', 'blue')
labels = c('Ae X', 'Ae Y')
k <- c(0.2, seq(0.4, 0.9, by=0.025), 0.95)
png("AeB3hist.png")
hist(CMB3$AeX, col="#0000ff40", border=F, breaks=k, main="Band-3 Aperture Efficiency", xlab="Ae", xlim=c(0.3, 0.95)) 
hist(CMB3$AeY, col="#ff000040", border=F, breaks=k, add=T)
legend("topright", legend=labels, col=c("#0000ff40", "#ff000040"), pch=15,cex=1)
dev.off()
png("AeB6hist.png")
hist(CMB6$AeX, col="#0000ff40", border=F, breaks=k, main="Band-6 Aperture Efficiency", xlab="Ae",  xlim=c(0.3, 0.95)) 
hist(CMB6$AeY, col="#ff000040", border=F, breaks=k, add=T)
legend("topright", legend=labels, col=c("#0000ff40", "#ff000040"), pch=15,cex=1)
dev.off()
png("AeB7hist.png")
hist(CMB7$AeX, col="#0000ff40", border=F, breaks=k, main="Band-7 Aperture Efficiency", xlab="Ae",  xlim=c(0.3, 0.95)) 
hist(CMB7$AeY, col="#ff000040", border=F, breaks=k, add=T)
legend("topright", legend=labels, col=c("#0000ff40", "#ff000040"), pch=15,cex=1)
dev.off()

png("AeB3secZ.png")
plot(1.0/sin(CMB3$EL), CMB3$AeX, xlim=c(1,2), xlab='SecZ', ylab='Aperture Efficiency', main='Band-3', pch=20, col='red')
points(1.0/sin(CMB3$EL), CMB3$AeY, pch=20, col='blue')
legend("topright", legend=labels, col=colores, pch=20)
dev.off()

png("AeB6secZ.png")
plot(1.0/sin(CMB6$EL), CMB6$AeX, xlim=c(1,2), xlab='SecZ', ylab='Aperture Efficiency', main='Band-6', pch=20, col='red')
points(1.0/sin(CMB6$EL), CMB6$AeY, pch=20, col='blue')
legend("topright", legend=labels, col=colores, pch=20)
dev.off()

png("AeB7secZ.png")
plot(1.0/sin(CMB7$EL), CMB7$AeX, xlim=c(1,2), xlab='SecZ', ylab='Aperture Efficiency', main='Band-7', pch=20, col='red')
points(1.0/sin(CMB7$EL), CMB7$AeY, pch=20, col='blue')
legend("topright", legend=labels, col=colores, pch=20)
dev.off()


antList <- sort(as.character(unique(CMB3$Ant)))
for(band in c(3,6,7)){
    BandDF <- AeDF[AeDF$Band == band,]
    pdf(sprintf('B%dAe.pdf', band))
    for(ant in antList){
        antDF <- BandDF[BandDF$Ant == ant,]
        plot(antDF$Date, antDF$AeX, ylim=c(0.0, 1.0), type='n', xlab='Observation Date', ylab='Aperture Efficiency', main=sprintf('Band %d %s', band, ant))
        points(antDF$Date, antDF$AeX, pch=20, col=colores[1])
        points(antDF$Date, antDF$AeY, pch=20, col=colores[2])
        legend("topright", legend=labels, col=colores, pch=20)
        text_sd <- sprintf('| %s | %.2f (%.2f) | %.2f (%.2f) |', ant, median(antDF$AeX)*100.0, sd(antDF$AeX)*100.0, median(antDF$AeY)*100.0, sd(antDF$AeY)*100.0)
        print(text_sd, quote=F)
    }
    dev.off()
    fit <- lm(BandDF, formula = 0.5*(AeX + AeY) ~ Tau0 + secZ + DayTime)
    text_sd <- summary(fit)
    print(text_sd)
}

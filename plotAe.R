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
B7 <- AeDF[AeDF$Band == 7,]
#-------- Antenna Separation
CMB3 <- B3[grep('CM', B3$Ant),]
CMB7 <- B7[grep('CM', B7$Ant),]
colores = c('red', 'blue')
labels = c('Ae X', 'Ae Y')
antList <- as.character(unique(CMB3$Ant))
for(band in c(3,7)){
    BandDF <- AeDF[AeDF$Band == band,]
    pdf(sprintf('B%dAe.pdf', band))
    for(ant in antList){
        antDF <- BandDF[BandDF$Ant == ant,]
        plot(antDF$Date, antDF$AeX, ylim=c(0.0, 1.0), type='n', xlab='Observation Date', ylab='Aperture Efficiency', main=sprintf('Band %d %s', band, ant))
        points(antDF$Date, antDF$AeX, pch=20, col=colores[1])
        points(antDF$Date, antDF$AeY, pch=20, col=colores[2])
        legend("topright", legend=labels, col=colores, pch=20)
        text_sd <- sprintf('%s  %.2f  %.2f', ant, median(antDF$AeX)*100.0, median(antDF$AeY)*100.0)
        print(text_sd, quote=F)
    }
    dev.off()
    fit <- lm(BandDF, formula = 0.5*(AeX + AeY) ~ relTime + Tau0 + secZ + + DayTime)
    text_sd <- summary(fit)
    print(text_sd)
}

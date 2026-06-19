library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/ALMAR/refs/heads/master/StatStokes.R", ssl.verifypeer = FALSE)))
BandFreq<-c(43.0,  75.0,  97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
timeWindow <- 20*86400  # 20 days
load('Flux.Rdata')
UIDList <- c()
for(freq in BandFreq){
    bandDF <- FLDF[abs(FLDF$Freq - freq) < 10.0,]
    if(nrow(bandDF) < 1000){ next }
    sourceList <- sort(unique(bandDF$Src))
    for(src in sourceList){
        srcDF <- bandDF[bandDF$Src == src,]
        if(nrow(srcDF) < 1000){ next }
        cat(sprintf('%s : %.1f GHz\n', src, freq))
        fluxDF <- data.frame(Date = sort(unique(srcDF$Date)), I=0.0)
        for(index in 1:nrow(fluxDF)){
            subDF <- srcDF[abs(as.numeric(srcDF$Date - fluxDF[index,]$Date)) < timeWindow,]
            weight <- exp(-5*(as.numeric(subDF$Date - fluxDF[index,]$Date)/timeWindow)^2)
            fluxDF[index,]$I <- weight %*% subDF$I / sum(weight)
            #fluxDF[index,]$I <- median(subDF$I)
            if( abs(median(srcDF[srcDF$Date == fluxDF[index,]$Date,]$I) -  fluxDF[index,]$I) > 0.5* fluxDF[index,]$I){  
                #cat(sprintf('%s\n', srcDF[srcDF$Date == fluxDF[index,]$Date,]$File))
                UIDList <- append( UIDList, srcDF[srcDF$Date == fluxDF[index,]$Date,]$File)
            }
        }
    }
}
UIDList <- sort(unique(UIDList))

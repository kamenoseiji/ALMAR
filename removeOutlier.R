library(RCurl)
library(KFAS)
library(arrow)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/ALMAR/refs/heads/master/StatStokes.R", ssl.verifypeer = FALSE)))
#source('../StatStokes.R')
BandFreq<-c(40.0,  75.0,  97.5, 154.9, 183.0, 233.0, 343.4, 400.0, 650.0, 800.0)
timeWindow <- 1000*86400  # 1000 days
Today <- as.numeric(Sys.time())
load('Flux.Rdata')
FLDF$QC <- 1        # Add quality code (0:don't use, 1:usable)
for(freq in BandFreq){
    bandDF <- FLDF[abs(FLDF$Freq - freq) < 10.0,]
    if(nrow(bandDF) < 1000){ next }
    sourceList <- sort(unique(bandDF$Src))
    for(src in sourceList){
        if( substr(src, 1,1) != 'J' ){ next }   # Filter SSOs out
        srcDF <- bandDF[bandDF$Src == src,]
        if(nrow(srcDF) < 100){ next }
        #-------- I, Q, U at the representative frequency
        df <- data.frame(relTime = as.numeric(unique(srcDF$Date) - Today), I=0)
        for(date_index in 1:nrow(df)){
            sdf <- srcDF[as.numeric(srcDF$Date) - Today == df[date_index,]$relTime,]
            sdf$relFreq <- sdf$Freq - freq
            fitI <- lm(formula = I ~ relFreq, data = sdf, weights = 1/eI^2); df$I[date_index] <- coef(summary(fitI))[1]
        }
        #-------- define local level model
        smthDF <- KFASsmooth(data.frame(time=df$relTime, value=df$I), 0.99)
        outlier_index <- which( ((df$I < smthDF$lower) | (df$I > smthDF$upper)))
        for(index in outlier_index){
            outlier_file <- srcDF[as.numeric(srcDF$Date) - Today == df$relTime[index],]$File[1]
            text_sd <- sprintf('%s Iobs= %.1f Jy / Ipred= [%.1f %.1f] Jy : %s\n', src, df$I[index], smthDF$lower[index], smthDF$upper[index], outlier_file)
            cat(text_sd)
            FLDF[as.numeric(FLDF$Date) - Today == df$relTime[index],]$QC <- FLDF[as.numeric(FLDF$Date) - Today == df$relTime[index],]$QC* 0.1   # Outlier flag
        }
    }
}
save(FLDF, file='Flux.Rdata')
write_parquet(FLDF, "Flux.parquet")

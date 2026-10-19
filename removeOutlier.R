library(RCurl)
library(KFAS)
library(arrow)
library(data.table)
library(parallel)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/ALMAR/refs/heads/master/StatStokes.R", ssl.verifypeer = FALSE)))
#source('../StatStokes.R')
numCore = detectCores()
BandFreq<-c(40.0,  75.0,  97.5, 154.9, 183.0, 233.0, 343.4, 400.0, 650.0, 800.0)
minObsNum <- 100  # Minumum number of observations
Today <- as.numeric(Sys.time())
load('Flux.Rdata')
FLDF$QC <- 1        # Add quality code (0:don't use, 1:usable)
for(freq in BandFreq){
    bandDF <- FLDF[abs(FLDF$Freq - freq) < 10.0,]
    if(nrow(bandDF) < minObsNum*10){ next }
    sourceList <- sort(unique(bandDF$Src))
    sourceList <- grep("^.[0-9]", sourceList, value = TRUE)
    bandDF <- bandDF[bandDF$Src %in% sourceList,]
    sourceCounts <- table(bandDF$Src)
    bandDF <- bandDF[bandDF$Src %in% names(sourceCounts[sourceCounts >= minObsNum]),]
    sourceList <- unique(bandDF$Src)
    findOutlier <- function(src){
        srcDF <- bandDF[bandDF$Src == src,]
        setDT(srcDF)
        df <- srcDF[, {
            sdf_relFreq <- Freq - freq
            fitI <- lm(I ~ sdf_relFreq, weights = 1/eI^2)
            .(I = coef(summary(fitI))[1, 1])
        }, by = .(relTime = as.numeric(Date) - Today)]
        #-------- define local level model
        smthDF <- KFASsmooth(data.frame(time=df$relTime, value=df$I), 0.999)
        outlier_index <- which( ((df$I < smthDF$lower) | (df$I > smthDF$upper)))
        return(Today + df$relTime[outlier_index])
    }
    target_dates <- sort(unique(unlist(mclapply(sourceList, findOutlier, mc.cores=numCore))))
    FLDF$QC[as.numeric(FLDF$Date) %in% target_dates] <- 0.1
}
FLDF <- FLDF[FLDF$QC > 0.9,]
save(FLDF, file='Flux.Rdata')
write_parquet(FLDF, "Flux.parquet")

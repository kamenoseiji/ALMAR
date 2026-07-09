library(mgcv)
#-------- Estimate Stokes parameters by frequency and date 
estimateIQUV <- function(DF, refFreq, refDate=Sys.time()){
    DF$relFreq <- DF$Freq / refFreq
    DF$relTime <- as.numeric(DF$Date) - as.numeric(refDate)
    df <- subset(DF[,c('relFreq', 'I', 'eI', 'Q', 'eQ', 'U', 'eU', 'V', 'eV', 'relTime')], is.finite(relFreq) & is.finite(I) & is.finite(eI) & is.finite(relTime) & I > 0 & eI > 0 & relFreq > 0)
    df$P  <- sqrt(df$Q^2 + df$U^2)
    df$eP <- sqrt(df$eQ^2 + df$eU^2)
    df$EVPA  <- 0.5* atan2(df$U, df$Q)
    df$eEVPA <- 0.5* sqrt(df$Q^2 * df$eU^2 + df$U^2 * df$eQ^2) / (df$P^2)
    IQUV <- data.frame(Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=NA, eI=0.0, eQ=0.0, eU=0.0, eV=0.0, eP=0.0, eEVPA=0.0)
    timeWeightSoftening <- 5* 86400 # 5-day softening
    if( (max(df$relFreq) < 0.65) | (max(df$relFreq) / min(df$relFreq) < 2.0 )){ return( IQUV )}
    df$timeFreqDeparture <- (abs(df$relTime) + timeWeightSoftening) * (1.0 + abs(df$relFreq - 1))
    weight <- 1.0/(abs(df$eEVPA) * sqrt(df$eQ^2 + df$eU^2)* abs(log(df$relFreq) + 1.0)^2 * (timeWeightSoftening / abs(df$relTime + timeWeightSoftening)))
    if( diff(range(df$relTime)) > timeWeightSoftening){
        fitI <- lm(formula=log(I) ~ log(relFreq) + relTime, data=df, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
        fitP <- lm(formula=log(P) ~ log(relFreq) + relTime, data=df, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
        df$V <- df$V*df$relFreq^coef(summary(fitI))[2]
        fitV <- lm(formula=V ~ relTime, data=df, weight=(I / eV) * (timeWeightSoftening / timeFreqDeparture))
        IQUV$I <- exp(coef(fitI)[[1]])
        IQUV$P <- exp(coef(fitP)[[1]])
        IQUV$eI <- IQUV$I* coef(summary(fitI))[4]
        IQUV$eP <- IQUV$P* coef(summary(fitP))[4]
    } else {    # insufficient time range
        fitI <- lm(formula=log(I) ~ log(relFreq), data=df, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
        fitP <- lm(formula=log(P) ~ log(relFreq), data=df, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
        fitV <- lm(formula=V ~ relFreq, data=df, weight=(I / eV) * (timeWeightSoftening / timeFreqDeparture))
        IQUV$I <- exp(coef(fitI)[[1]])
        IQUV$P <- exp(coef(fitP)[[1]])
        IQUV$eI <- IQUV$I* coef(summary(fitI))[3]
        IQUV$eP <- IQUV$P* coef(summary(fitP))[3]
    }
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*df$EVPA) ) / sum(weight); Twiddle <- Twiddle/abs(Twiddle)
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
    IQUV$EVPA  <- 0.5* Arg(Twiddle)
    IQUV$eQ <- IQUV$Q* IQUV$eP / IQUV$P
    IQUV$eU <- IQUV$U* IQUV$eP / IQUV$P
    IQUV$eEVPA <- 0.5*sqrt(IQUV$Q^2 * IQUV$eU^2 + IQUV$U^2 * IQUV$eQ^2) / (IQUV$P^2)
    IQUV$V <- coef(summary(fitV))[1]
    IQUV$eV <- coef(summary(fitV))[4]
    return(IQUV)
}
#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
predStokes <- function(df){
    bandID   <- df$Band[1]
    fitI <- lm(formula=I ~ Freq, data=df, weight=1.0/eI^2)
    fitQ <- lm(formula=Q ~ Freq, data=df, weight=1.0/eQ^2)
    fitU <- lm(formula=U ~ Freq, data=df, weight=1.0/eU^2)
    fitV <- lm(formula=V ~ Freq, data=df, weight=1.0/eV^2)
    newDF <- data.frame(Src=df$Src[1], Freq = standardFreq[[bandID]], EL=df$EL[1])
    pred <- as.numeric(predict(fitI, newDF, interval='confidence', level=0.67)); newDF$I <- matrix(pred, ncol=3)[,1]; newDF$eI <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitQ, newDF, interval='confidence', level=0.67)); newDF$Q <- matrix(pred, ncol=3)[,1]; newDF$eQ <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitU, newDF, interval='confidence', level=0.67)); newDF$U <- matrix(pred, ncol=3)[,1]; newDF$eU <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    pred <- as.numeric(predict(fitV, newDF, interval='confidence', level=0.67)); newDF$V <- matrix(pred, ncol=3)[,1]; newDF$eV <- 0.5*(matrix(pred, ncol=3)[,3] - matrix(pred, ncol=3)[,2])
    newDF$Date <- df$Date[1]
    newDF$File <- df$File[1]
    return(newDF)
}
#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
SPfit <- function(df){
    df <- rbind(df, data.frame(relTime=abs(df[1,]$relTime)*seq(0.5, 1.4, by=0.1), Value=rep(median(df$Value), 10), error=rep(sum(df$error), 10)))
    gam_model <- gam(Value ~ s(relTime, bs="cr"), data=df, weights=1/error^2, sp=0.3*abs(df[1,]$relTime))
    pred_gam <- predict(gam_model, newdata=data.frame(relTime=0.0), se.fit=TRUE)
    return( pred_gam$fit  + (0.0 + 1.0i)*pred_gam$se.fit )
}

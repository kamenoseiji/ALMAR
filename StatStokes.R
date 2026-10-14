#library(mgcv)
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
    if( diff(range(df$relTime)) > 1.5*min(abs(df$relTime))){
        fitI <- lm(formula=log(I) ~ log(relFreq) + relTime, data=df, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
        SPIX <- coef(fitI)[[2]] # Spectral Index
        IQUV$I <- exp(coef(fitI)[[1]])
        df$P <- df$P * exp(-SPIX* log(df$relFreq))
        df$eP <- df$eP * exp(-SPIX* log(df$relFreq))
        fitP <- lm(P ~ relTime, data=df, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
        df$V <- df$V*df$relFreq^coef(summary(fitI))[2]
        fitV <- lm(formula=V ~ relTime, data=df, weight=(I / eV) * (timeWeightSoftening / timeFreqDeparture))
        IQUV$P <- coef(fitP)[[1]]
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
SPfit <- function(df){
    #-------- Compact Log-Likelihood Function for smoothing
    scaledTime <- df$relTime / 1e6              # to avoid overlow
    DistMatrix <- as.matrix(dist(scaledTime))
    gp_nll <- function(poly) {      # polynimial parameters
        mu    <- poly[1] + poly[2] * scaledTime
        Sigma <- exp(poly[3]) * exp(-(DistMatrix / exp(poly[4]))^2) + diag(df$error^2)
        L     <- chol(Sigma)
        0.5 * (2 * sum(log(diag(L))) + sum(backsolve(L, df$Value - mu, transpose = TRUE)^2))
    }
    #-------- Fit & Extract Optimal Parameters
    fit <- nlminb(c(0.02, 0.01, -1, 0.01), gp_nll)$par
    b0  <- fit[1]; b1 <- fit[2]; sig2 <- exp(fit[3]); range <- exp(fit[4])
    #-------- Predict at relTime = 0.0
    r_target   <- sig2 * exp(-(abs(scaledTime - 0) / range)^2)
    Sigma_opt  <- sig2 * exp(-(DistMatrix / range)^2) + diag(df$error^2)
    L <- chol(Sigma_opt)
    residual <- df$Value - (b0 + b1 * scaledTime)
    alpha    <- backsolve(L, residual, transpose = TRUE)
    prediction <- as.numeric(b0 + t(r_target) %*% backsolve(L, alpha))
    return(prediction + (0.0 + 1.0i)* sqrt(as.numeric(sig2)))
}

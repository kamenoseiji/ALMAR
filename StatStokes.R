#-------- Estimate Stokes parameters by freqneyc and date 
estimateIQUV <- function(DF, refFreq){
    DF$relFreq <- DF$Freq / refFreq
    timeWeightSoftening <- 5* 86400 # 5-day softening
    if( (max(DF$relFreq) < 0.65) | (max(DF$relFreq) / min(DF$relFreq) < 2.0 )){ return( data.frame( Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=0.0))}
    DF$timeFreqDeparture <- (abs(DF$relTime) + timeWeightSoftening) * (1.0 + abs(DF$relFreq - 1))
    if( (diff(range(DF$relTime)) < 0.25* DateRange* SECPERDAY) | (max(DF$relTime) < -0.25* SECPERDAY)){      # small number of data
        fitI <- lm(formula=log(I) ~ log(relFreq), data=DF, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
        fitP <- lm(formula=log(P) ~ log(relFreq), data=DF, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
    } else {
        fitI <- lm(formula=log(I) ~ log(relFreq) + relTime, data=DF, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
        fitP <- lm(formula=log(P) ~ log(relFreq) + relTime, data=DF, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
    }
    weight <- 1.0/(abs(DF$eEVPA) * sqrt(DF$eQ^2 + DF$eU^2)* abs(log(DF$relFreq) + 1.0)^2 * (timeWeightSoftening / abs(DF$relTime + timeWeightSoftening)))
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*DF$EVPA) ) / sum(weight); Twiddle <- Twiddle/abs(Twiddle)
    IQUV <- data.frame(Src=DF$Src[1], I=exp(coef(fitI)[[1]]), Q=0.0, U=0.0, V=0.0, P=exp(coef(fitP)[[1]]), EVPA=0.5*Arg(Twiddle))
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
    return( IQUV )
}

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
    fitI <- lm(formula=log(I) ~ log(relFreq) + relTime, data=df, weight=(I / eI) * (timeWeightSoftening / timeFreqDeparture))
    fitP <- lm(formula=log(P) ~ log(relFreq) + relTime, data=df, weight=(P / eP) * (timeWeightSoftening / timeFreqDeparture))
    df$V <- df$V*df$relFreq^coef(summary(fitI))[2]
    fitV <- lm(formula=V ~ relTime, data=df, weight=(I / eV) * (timeWeightSoftening / timeFreqDeparture))
    weight <- 1.0/(abs(df$eEVPA) * sqrt(df$eQ^2 + df$eU^2)* abs(log(df$relFreq) + 1.0)^2 * (timeWeightSoftening / abs(df$relTime + timeWeightSoftening)))
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*df$EVPA) ) / sum(weight); Twiddle <- Twiddle/abs(Twiddle)
    IQUV$I <- exp(coef(fitI)[[1]])
    IQUV$eI <- IQUV$I* coef(summary(fitI))[4]
    IQUV$P <- exp(coef(fitP)[[1]])
    IQUV$eP <- IQUV$P* coef(summary(fitP))[4]
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
#library(KFAS)
#KalmanIQUV <- function(DF, refFreq){
#    #-------- Kalman filtering for Stokes I and spectral index
#    datI <- subset(DF[,c("relFreq", "I", "eI", "relTime")], is.finite(relFreq) & is.finite(I) & is.finite(eI) & is.finite(relTime) & I > 0 & eI > 0 & relFreq > 0)
#    #-------- Observation equation : logI = log(I)
#    datI$logI <- log(datI$I); datI$var_logI <- (datI$eI / datI$I)^2
#    #-------- Build state-space model : State vector: a_t = [beta0_t, beta1_t], Observation: y_t = [1, log(Freq_t)] a_t + eps_t
#    Z_array <- array(NA_real_, dim = c(1, 2, nrow(datI)))
#    Z_array[1,1,] <- 1
#    Z_array[1,2,] <- log(datI$relFreq)
#    #-------- Random-walk state evolution
#    build_model <- function(log_q0, log_q1) {
#        Qmat <- diag(c(exp(log_q0), exp(log_q1)))
#        SSModel( datI$logI ~ -1 + SSMcustom(Z = Z_array, T = diag(2), R = diag(2), Q = Qmat, a1 = c(mean(datI$logI), 0), P1 = diag(1e6, 2)), H = array(datI$var_logI, c(1,1,nrow(datI))))
#    }
#    #-------- Maximum likelihood estimation of state variances
#    init_par <- log(c( var(datI$logI, na.rm = TRUE) * 1e-3, 1e-3))
#    fit <- fitSSM( inits = init_par, model = build_model(init_par[1], init_par[2]), updatefn = function(pars, model) {
#            model$Q[,,1] <- diag(exp(pars))
#            model }, method = "BFGS")
#    model_fit <- fit$model
#    #-------- Kalman smoothing
#    kfs <- KFS( model_fit, smoothing = c("state", "signal"))
#    #-------- Estimate state at t = 0
#    beta_hat <- kfs$alphahat
#    last_state <- beta_hat[nrow(datI), ]
#    V_last <- kfs$V[,,nrow(datI)] # State covariance at final epoch
#    V0 <- V_last
#    beta0_now <- last_state[1]
#    alpha_now <- last_state[2]
#    se_beta0_now <- sqrt(V0[1,1]) 
#    se_alpha_now <- sqrt(V0[2,2])
#    S_now <- exp(beta0_now)
#    se_S_now <- S_now * se_beta0_now
#    return(c(S_now, se_S_now, alpha_now, se_alpha_now))
#}

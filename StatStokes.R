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
#-------- Input multiple frequency data and output Stokes parameters at the standard frequency
SPfit <- function(df){
    df <- rbind(df, data.frame(relTime=abs(df[1,]$relTime)*seq(0.5, 1.4, by=0.1), Value=rep(median(df$Value), 10), error=rep(sum(df$error), 10)))
    gam_model <- gam(Value ~ s(relTime, bs="cr"), data=df, weights=1/error^2, sp=0.3*abs(df[1,]$relTime))
    pred_gam <- predict(gam_model, newdata=data.frame(relTime=0.0), se.fit=TRUE)
    return( pred_gam$fit  + (0.0 + 1.0i)*pred_gam$se.fit )
}

#-------- New: Kalman-filter-based estimator for Stokes parameters
# Adds kalman_linear_regression helper and estimateIQUV_kalman function.
# The Kalman filter treats regression coefficients as the state with identity evolution (random walk)
# and uses measurement variances derived from the provided errors.

kalman_linear_regression <- function(y, X, R_vec, Q_scale = 1e-6){
    # y: vector of observations (numeric)
    # X: design matrix (n x p)
    # R_vec: measurement variances (length n)
    # Q_scale: scalar process noise scale for state covariance (diagonal)
    n <- nrow(X)
    p <- ncol(X)
    m <- rep(0, p)                # initial state mean
    P <- diag(1e6, p)             # initial large uncertainty
    Qmat <- diag(Q_scale, p)
    for(i in seq_len(n)){
        xi <- matrix(X[i,], nrow=1)
        # predict (identity evolution)
        m_pred <- m
        P_pred <- P + Qmat
        # measurement update
        Ri <- as.numeric(R_vec[i])
        S <- as.numeric(xi %*% P_pred %*% t(xi) + Ri)
        if(S <= 0) S <- S + 1e-12
        K <- as.numeric((P_pred %*% t(xi)) / S)
        innov <- as.numeric(y[i] - as.numeric(xi %*% m_pred))
        m <- as.numeric(m_pred + K %*% innov)
        P <- (diag(p) - K %*% xi) %*% P_pred
    }
    list(m = m, P = P)
}

estimateIQUV_kalman <- function(DF, refFreq, refDate=Sys.time()){
    # Prepare output with same structure as estimateIQUV
    IQUV <- data.frame(Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=NA,
                       eI=0.0, eQ=0.0, eU=0.0, eV=0.0, eP=0.0, eEVPA=0.0)
    if(nrow(DF) < 2) return(IQUV)

    DF$relFreq <- DF$Freq / refFreq
    DF$relTime <- as.numeric(DF$Date) - as.numeric(refDate)
    df <- subset(DF[,c('relFreq', 'I', 'eI', 'Q', 'eQ', 'U', 'eU', 'V', 'eV', 'relTime')],
                 is.finite(relFreq) & is.finite(I) & is.finite(eI) & is.finite(relTime) & I > 0 & eI > 0 & relFreq > 0)
    if(nrow(df) < 2) return(IQUV)

    # quick checks from original function
    if( (max(df$relFreq) < 0.65) | (max(df$relFreq) / min(df$relFreq) < 2.0 )){ return( IQUV )}

    # design matrix: intercept, log(relFreq), relTime
    X <- cbind(1.0, log(df$relFreq), df$relTime)

    # ----- Intensity (log-space) -----
    yI <- log(df$I)
    RI <- (df$eI / df$I)^2
    resI <- tryCatch(kalman_linear_regression(yI, X, RI, Q_scale = 1e-6), error=function(e) NULL)
    if(!is.null(resI)){
        xi_ref <- c(1.0, 0.0, 0.0)
        logI_ref <- as.numeric(xi_ref %*% resI$m)
        var_logI <- as.numeric(xi_ref %*% resI$P %*% xi_ref)
        IQUV$I <- exp(logI_ref)
        IQUV$eI <- IQUV$I * sqrt(max(var_logI, 0))
    }

    # ----- Polarized intensity P (log-space) -----
    df$P <- sqrt(df$Q^2 + df$U^2)
    df$eP <- sqrt(df$eQ^2 + df$eU^2)
    okP <- which(is.finite(df$P) & is.finite(df$eP) & (df$P > 0) & (df$eP > 0))
    if(length(okP) >= 2){
        yP <- log(df$P[okP])
        RP <- (df$eP[okP] / df$P[okP])^2
        XP <- X[okP, , drop=FALSE]
        resP <- tryCatch(kalman_linear_regression(yP, XP, RP, Q_scale = 1e-6), error=function(e) NULL)
        if(!is.null(resP)){
            xi_ref <- c(1.0, 0.0, 0.0)
            logP_ref <- as.numeric(xi_ref %*% resP$m)
            var_logP <- as.numeric(xi_ref %*% resP$P %*% xi_ref)
            IQUV$P <- exp(logP_ref)
            IQUV$eP <- IQUV$P * sqrt(max(var_logP, 0))
        }
    }

    # compute per-observation EVPA for optional weighting/fallback
    df$EVPA  <- 0.5* atan2(df$U, df$Q)
    df$eEVPA <- 0.5* sqrt(df$Q^2 * df$eU^2 + df$U^2 * df$eQ^2) / (pmax(df$P^2, 1e-12))

    # ----- Q, U, V (linear space) -----
    if(sum(is.finite(df$Q) & is.finite(df$eQ)) >= 2){
        yQ <- df$Q
        RQ <- (df$eQ)^2
        resQ <- tryCatch(kalman_linear_regression(yQ, X, RQ, Q_scale = 1e-6), error=function(e) NULL)
        if(!is.null(resQ)){
            xi_ref <- c(1.0, 0.0, 0.0)
            Q_ref <- as.numeric(xi_ref %*% resQ$m)
            varQ <- as.numeric(xi_ref %*% resQ$P %*% xi_ref)
            IQUV$Q <- Q_ref
            IQUV$eQ <- sqrt(max(varQ, 0))
        }
    }
    if(sum(is.finite(df$U) & is.finite(df$eU)) >= 2){
        yU <- df$U
        RU <- (df$eU)^2
        resU <- tryCatch(kalman_linear_regression(yU, X, RU, Q_scale = 1e-6), error=function(e) NULL)
        if(!is.null(resU)){
            xi_ref <- c(1.0, 0.0, 0.0)
            U_ref <- as.numeric(xi_ref %*% resU$m)
            varU <- as.numeric(xi_ref %*% resU$P %*% xi_ref)
            IQUV$U <- U_ref
            IQUV$eU <- sqrt(max(varU, 0))
        }
    }
    if(sum(is.finite(df$V) & is.finite(df$eV)) >= 2){
        yV <- df$V
        RV <- (df$eV)^2
        resV <- tryCatch(kalman_linear_regression(yV, X, RV, Q_scale = 1e-6), error=function(e) NULL)
        if(!is.null(resV)){
            xi_ref <- c(1.0, 0.0, 0.0)
            V_ref <- as.numeric(xi_ref %*% resV$m)
            varV <- as.numeric(xi_ref %*% resV$P %*% xi_ref)
            IQUV$V <- V_ref
            IQUV$eV <- sqrt(max(varV, 0))
        }
    }

    # If Q and U estimates are available compute EVPA and uncertainty
    if(!is.na(IQUV$Q) & !is.na(IQUV$U) & (IQUV$Q != 0 | IQUV$U != 0)){
        IQUV$EVPA <- 0.5 * atan2(IQUV$U, IQUV$Q)
        denom <- (IQUV$Q^2 + IQUV$U^2)
        if(denom > 0){
            dQ <- -0.5 * IQUV$U / denom
            dU <-  0.5 * IQUV$Q / denom
            varQ <- IQUV$eQ^2
            varU <- IQUV$eU^2
            varEVPA <- dQ^2 * varQ + dU^2 * varU
            IQUV$eEVPA <- sqrt(max(varEVPA, 0))
        }
    } else {
        # fallback: if P estimated but Q/U not individually available, attempt to recover EVPA via weighted complex mean
        if(IQUV$P > 0){
            # prepare weights - avoid zeros/NA
            w_eEVPA <- ifelse(is.finite(df$eEVPA) & df$eEVPA > 0, df$eEVPA, 1.0)
            denom_freq <- (1 + abs(log(df$relFreq)))^2
            denom_err <- pmax(sqrt(df$eQ^2 + df$eU^2), 1e-12)
            weight <- 1.0/( w_eEVPA * denom_err * denom_freq )
            phasors <- (df$Q + 1i * df$U) / pmax(df$P, 1e-12)
            Twiddle <- sum(weight * phasors) / sum(weight)
            Twiddle <- Twiddle / abs(Twiddle)
            IQUV$Q <- IQUV$P * Re(Twiddle)
            IQUV$U <- IQUV$P * Im(Twiddle)
            IQUV$EVPA <- 0.5 * Arg(Twiddle)
            IQUV$eQ <- abs(IQUV$Q) * (IQUV$eP / IQUV$P)
            IQUV$eU <- abs(IQUV$U) * (IQUV$eP / IQUV$P)
            IQUV$eEVPA <- 0.5*sqrt(IQUV$Q^2 * IQUV$eU^2 + IQUV$U^2 * IQUV$eQ^2) / (IQUV$P^2)
        }
    }

    return(IQUV)
}

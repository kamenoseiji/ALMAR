load('Flux.Rdata')
sourceList <- sort(as.character(unique(FLDF$Src)))
pdf('StokesV.pdf')
SrcI <- SrcVfrac <- SrcVe <- numeric(0)
for(source in sourceList){
	SrcDF <- FLDF[FLDF$Src == source,]
	nSmp <- length(SrcDF$I)
	if( nSmp == 1){
		Imean <- SrcDF$I
		Ierr  <- SrcDF$eI
		Pmean <- SrcDF$P/SrcDF$I
		Perr  <- SrcDF$eP/SrcDF$I
		Vmean <- SrcDF$V/SrcDF$I
		Verr  <- SrcDF$eV/SrcDF$I
	} else {
		SrcDF$dEL <- SrcDF$EL - median(SrcDF$EL)
		weight <- 1.0/SrcDF$eI^2; Imean <- weighted.mean(SrcDF$I, weight); Ierr <- sqrt(cov.wt(cbind(SrcDF$I,SrcDF$I), weight)[[1]][1,1] / nSmp)
		weight <- 1.0/SrcDF$eQ^2; Qmean <- weighted.mean(SrcDF$Q, weight); Qerr <- sqrt(cov.wt(cbind(SrcDF$Q,SrcDF$Q), weight)[[1]][1,1] / nSmp)
		weight <- 1.0/SrcDF$eU^2; Umean <- weighted.mean(SrcDF$U, weight); Uerr <- sqrt(cov.wt(cbind(SrcDF$U,SrcDF$U), weight)[[1]][1,1] / nSmp)
		weight <- 1.0/SrcDF$eV^2; Vmean <- weighted.mean(SrcDF$V, weight) / Imean; Verr <- sqrt((cov.wt(cbind(SrcDF$V,SrcDF$V), weight)[[1]][1,1] + min(SrcDF$eV^2))/ nSmp) / Imean
		Pmean <- sqrt(Qmean^2 + Umean^2) / Imean; Perr <- sqrt(Qerr^2 + Uerr^2)/Imean
	}
	SrcI     <- c(SrcI, Imean)
	SrcVfrac <- c(SrcVfrac, Vmean)
	SrcVe    <- c(SrcVe, Verr)
	text_sd <- sprintf('%10s %3d %8.4f (%6.4f)  %4.1f%% (%5.2f%%)  %7.4f%% (%7.4f%%) \n', source, length(SrcDF$I),
				Imean, Ierr, 100.0*Pmean, 100.0*Perr, 100.0*Vmean, 100.0*Verr)
	cat(text_sd)
}
plot(SrcI, 100.0* SrcVfrac, xlab='Stokes I [Jy]', ylab='Fractional Circular Polarization [%]', pch=20, log='x', ylim=c(-3.0, 3.0), main=sprintf('%.1f GHz', FLDF$Freq[1]))
abline(h=-0.5, lty=2); abline(h=0.5, lty=2)
arrows(SrcI, 100.0* (SrcVfrac - SrcVe), SrcI, 100.0* (SrcVfrac + SrcVe), angle=90, length=0.0)
hist(100.0*abs(SrcVfrac), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 100), col = "#ff00ff40", border = "#ff00ff", freq=FALSE, xlim=c(0, 1.5), xlab='| Stokes  V | (%)', main=sprintf('%% Circular Polarization (%.1f GHz)', FLDF$Freq[1]))
hist(abs(SrcVfrac) / SrcVe, breaks=c(0,1,2,3,4,5,100), col = "#0000ff40", border = "#0000ff", freq=FALSE, xlim=c(0, 5.0), xlab='| Stokes  V |/ sigma', main=sprintf('Significance of Circular Polarization (%.1f GHz)', FLDF$Freq[1]))

dev.off()
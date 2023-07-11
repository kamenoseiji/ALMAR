load('~/ALMAR/WORK/Flux.Rdata')

sourceName <- 'J1924-2914'
upperFlux <- 20.0	# [Jy]

DF <- FLDF[FLDF$Src == sourceName,]
outDF <- DF[DF$I > upperFlux,]
fileList <- unique(outDF$File)
cat(fileList)
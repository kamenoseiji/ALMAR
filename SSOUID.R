#-------- Start Program
if(file.exists('UIDList')){ file.remove('UIDList') }
if(file.exists('Flux.Rdata')){
    load('Flux.Rdata')
} else {
    load(url("https://www.alma.cl/~skameno/AMAPOLA/Flux.Rdata"))
}
sourceList <- FLDF$Src
QSO_index <- which((substr(sourceList, 1,1) == 'J') & grepl('[0-9]', substr(sourceList, 2,2)))
fileList <- unique( FLDF[-QSO_index,]$File )
write(file='UIDList', fileList)

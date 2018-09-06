import numpy as np
import glob
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
#-------- indexing
def indexList( refArray, motherArray ):     # Compare two arrays and return matched index
    IL = []
    for currentItem in refArray: IL = IL + np.where( motherArray == currentItem )[0].tolist()
    return IL
#
#-------- Weighted mean of 4 SPWs
def meanFlux(Flux, Ferr):
    if np.min(Ferr) < 1.0e-9 : return np.median(Flux, axis=1), np.median(Ferr, axis=1)
    weight = 1.0 / Ferr**2
    meanFlux = np.sum(Flux* weight, axis=1) / np.sum(weight, axis=1)
    meanErr  = 1.0 / np.sqrt(np.sum(weight, axis=1))
    return meanFlux, meanErr
#
#-------- Each UID
logfile = open('StokesComp.log', 'w')
text_sd = 'Src Band mjdSec EL I Q U V eI eQ eU eV aI aQ aU aV eaI eaQ eaU eaV UID'; logfile.write(text_sd + '\n')
SSODIR = 'SSOCal/'
APRDIR = 'AprioriCal/'
fileList = glob.glob(SSODIR + 'uid*.Aeff.npy')
for index in range(len(fileList)): fileList[index] = fileList[index][7:]
for AeFile in fileList:
    UID = AeFile[0:-15]
    uidBand = AeFile[0:-9]
    bandName = AeFile[-14:-9]
    band = int(AeFile[-11:-9])
    try:
        EL      = np.load(SSODIR + uidBand + '.EL.npy')
    except:
        print '%s is missing from %s' % (UID, SSODIR)
        continue
    antList = np.load(SSODIR + uidBand + '.AntList.npy')
    Flux  = np.load(SSODIR + uidBand + '.Flux.npy')   # Flux[scan, spw, Stokes]
    Ferr  = np.load(SSODIR + uidBand + '.Ferr.npy')   # Ferr[scan, spw, Stokes]
    SourceList = np.load(SSODIR + uidBand + '.Source.npy')
    atmTime    = np.load(SSODIR + uidBand + '.atmTime.npy')
    #
    if 'D' in antList: continue
    try:
        aFlux  = np.load(APRDIR + uidBand + '.Flux.npy')   # aFlux[scan, spw, Stokes]
    except:
        print '%s is missing from %s' % (UID, APRDIR)
        continue
    aFerr  = np.load(APRDIR + uidBand + '.Ferr.npy')   # aFerr[scan, spw, Stokes]
    aSourceList = np.load(APRDIR + uidBand + '.Source.npy')
    scanNum = Flux.shape[0]
    #---- Statistics
    meanF, errF = meanFlux(Flux, Ferr)
    meanA, errA = meanFlux(aFlux, aFerr)
    #---- Output
    for scan_index in range(scanNum):
        if SourceList[scan_index][0] != 'J' : continue
        text_sd = '%s %d %.7f %.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s' % (SourceList[scan_index], band, np.median(atmTime[scan_index]), EL[scan_index], meanF[scan_index,0], meanF[scan_index,1], meanF[scan_index,2], meanF[scan_index,3], errF[scan_index,0], errF[scan_index,1], errF[scan_index,2], errF[scan_index,3], meanA[scan_index,0], meanA[scan_index,1], meanA[scan_index,2], meanA[scan_index,3], errA[scan_index,0], errA[scan_index,1], errA[scan_index,2], errA[scan_index,3], UID)
        logfile.write(text_sd + '\n')
    #
    """
    antMap = indexList(antList, TrxAnt)             # TrxAnt[antMap] to align in AntList
    antNum = len(antMap)
    AePol = np.median(Aeff, axis=2)
    chNum = Tau0.shape[1]; chRange = range(int(0.05*chNum), int(0.95*chNum))
    Tau0  = np.median(Tau0[:,chRange])
    TauE  = np.median(TauE, axis=0)
    Trx = np.median(Trx[:,:,:,chRange], axis=(0,3))
    SSOindex = indexList(np.array(SSOCatalog), SourceList)
    SSOindex = SSOindex[np.argmax(EL[SSOindex])]
    SSO_EL = EL[SSOindex]
    for ant_index in range(antNum):
        antID = antMap[ant_index]
        text_sd = '%s %d %.7f %.3f %.4f %.1f %.1f %.3f %.3f %.3f %s %s' % (antList[ant_index], band, np.median(atmTime[SSOindex]), Tau0, np.std(TauE), Trx[antID,0], Trx[antID,1], AePol[ant_index,0], AePol[ant_index,1], SSO_EL, SourceList[SSOindex], UID)
        logfile.write(text_sd + '\n')
    #
    """
#
logfile.close()    

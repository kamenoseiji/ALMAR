import numpy as np
import glob
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
#-------- indexing
def indexList( refArray, motherArray ):     # Compare two arrays and return matched index
    IL = []
    for currentItem in refArray: IL = IL + np.where( motherArray == currentItem )[0].tolist()
    return IL
#
#-------- Each UID
logfile = open('Aeff.log', 'w')
text_sd = 'Ant Band mjdSec Tau0 TauRMS TrxX TrxY AeX AeY EL SSO UID'; logfile.write(text_sd + '\n')
fileList = glob.glob("uid*.Aeff.npy")
for AeFile in fileList:
    UID = AeFile[0:(len(AeFile)-15)]
    uidBand = AeFile[0:(len(AeFile)-9)]
    bandName = AeFile[(len(AeFile)-14):(len(AeFile)-9)]
    band = int(AeFile[(len(AeFile)-11):(len(AeFile)-9)])
    try:
        EL      = np.load(uidBand + '.EL.npy')
    except:
        continue
    Aeff = np.load(AeFile)
    antList = np.load(uidBand + '.AntList.npy')
    TrxAnt     = np.load(uidBand + '.TrxAnt.npy')
    SourceList = np.load(uidBand + '.Source.npy')
    Tau0       = np.load(uidBand + '.Tau0.npy')
    TauE       = np.load(uidBand + '.TauE.npy')
    Trx        = np.load(uidBand + '.Trx.npy')
    atmTime    = np.load(uidBand + '.atmTime.npy')
    #
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
#
logfile.close()    

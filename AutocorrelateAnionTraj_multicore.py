import numpy as np
import matplotlib.pyplot as plt
import sys
import multiprocessing
from scipy import signal

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    resultSize = int(result.size/2)
    result = result[resultSize:]
    normRange = range(len(result),0,-1)
    result = result/normRange
    return result
    
def runAutoCorrX(ionNo):
    print("Autocorrelating anion in X: %i of %i" %(ionNo+1, noIons))
    autoCorrTempX = autocorr(anionXVel[:,ionNo])
    return autoCorrTempX
def runAutoCorrY(ionNo):
    print("Autocorrelating anion in Y: %i of %i" %(ionNo+1, noIons))
    autoCorrTempY = autocorr(anionYVel[:,ionNo])
    return autoCorrTempY
def runAutoCorrZ(ionNo):
    print("Autocorrelating anion in Z: %i of %i" %(ionNo+1, noIons))
    autoCorrTempZ = autocorr(anionZVel[:,ionNo])
    return autoCorrTempZ
    
def runAutoCorrXDrift(ionNo):
    print("Autocorrelating drift corrected anion in X: %i of %i" %(ionNo+1, noIons))
    autoCorrTempX = autocorr(anionXVelDriftCorr[:,ionNo])
    return autoCorrTempX
def runAutoCorrYDrift(ionNo):
    print("Autocorrelating drift corrected anion in Y: %i of %i" %(ionNo+1, noIons))
    autoCorrTempY = autocorr(anionYVelDriftCorr[:,ionNo])
    return autoCorrTempY
def runAutoCorrZDrift(ionNo):
    print("Autocorrelating drift corrected anion in Z: %i of %i" %(ionNo+1, noIons))
    autoCorrTempZ = autocorr(anionZVelDriftCorr[:,ionNo])
    return autoCorrTempZ

atomMassDictionary = {"C":11.611,"H":1.008,"F":18.998,"N":13.607,"O":15.599,"S":31.666,"X":0.4}
timestep = 1.0 #timestep of trajectory in fs
cationMass =  139.222
anionMass = 280.145
noIons = 256
tSteps = 1000001

fieldStrength = "00"

anionXVel = np.genfromtxt("AnionXTrajectory.txt",delimiter=";"); anionXAvgVel = np.mean(anionXVel)
anionYVel = np.genfromtxt("AnionYTrajectory.txt",delimiter=";"); anionYAvgVel = np.mean(anionYVel)
anionZVel = np.genfromtxt("AnionZTrajectory.txt",delimiter=";"); anionZAvgVel = np.mean(anionZVel)

# Calculate drift corrected velocity arrays
anionXVelDriftCorr = anionXVel - anionXAvgVel
anionYVelDriftCorr = anionYVel - anionYAvgVel
anionZVelDriftCorr = anionZVel - anionZAvgVel

#Pre-allocate autocorrelation vectors
autoCorrAnX = np.zeros([tSteps,noIons]); autoCorrAnY = np.zeros([tSteps,noIons]); autoCorrAnZ = np.zeros([tSteps,noIons])
autoCorrAnXDriftCorr = np.zeros([tSteps,noIons]); autoCorrAnYDriftCorr = np.zeros([tSteps,noIons]); autoCorrAnZDriftCorr = np.zeros([tSteps,noIons])

if __name__ == '__main__':     

    diffpool = multiprocessing.Pool(processes=60)
    
    autoCorrAnX = diffpool.map(runAutoCorrX, range(noIons))
    autoCorrAnY = diffpool.map(runAutoCorrY, range(noIons))
    autoCorrAnZ = diffpool.map(runAutoCorrZ, range(noIons))
    
    autoCorrAnXDriftCorr = diffpool.map(runAutoCorrXDrift, range(noIons))
    autoCorrAnYDriftCorr = diffpool.map(runAutoCorrYDrift, range(noIons))
    autoCorrAnZDriftCorr = diffpool.map(runAutoCorrZDrift, range(noIons))
   
autoCorrAnX = np.sum(autoCorrAnX, axis=0);autoCorrAnX = autoCorrAnX/noIons
autoCorrAnY = np.sum(autoCorrAnY, axis=0);autoCorrAnY = autoCorrAnY/noIons
autoCorrAnZ = np.sum(autoCorrAnZ, axis=0);autoCorrAnZ = autoCorrAnZ/noIons

autoCorrAnXDriftCorr = np.sum(autoCorrAnXDriftCorr, axis=0);autoCorrAnXDriftCorr = autoCorrAnXDriftCorr/noIons
autoCorrAnYDriftCorr = np.sum(autoCorrAnYDriftCorr, axis=0);autoCorrAnYDriftCorr = autoCorrAnYDriftCorr/noIons
autoCorrAnZDriftCorr = np.sum(autoCorrAnZDriftCorr, axis=0);autoCorrAnZDriftCorr = autoCorrAnZDriftCorr/noIons

###Calculate integrals###
anXInt = np.cumsum(autoCorrAnX*timestep/1000)  #factor of 1000 converts timestep in fs to ps
anYInt = np.cumsum(autoCorrAnY*timestep/1000)
anZInt = np.cumsum(autoCorrAnZ*timestep/1000)
anXDriftInt = np.cumsum(autoCorrAnXDriftCorr*timestep/1000)
anYDriftInt = np.cumsum(autoCorrAnYDriftCorr*timestep/1000)
anZDriftInt = np.cumsum(autoCorrAnZDriftCorr*timestep/1000)

#Save autocorrelations
anXOutFile = "Anion_X.txt"
anxOutWrite = open(anXOutFile,"w")
anxOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n") 
for i in range(len(autoCorrAnX)):
    anxOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnX[i],anXInt[i]))
anxOutWrite.close()

anYOutFile = "Anion_Y.txt"
anyOutWrite = open(anYOutFile,"w")
anyOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrAnY)):
    anyOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnY[i],anYInt[i]))
anyOutWrite.close()

anZOutFile = "Anion_Z.txt"
anzOutWrite = open(anZOutFile,"w")
anzOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n") 
for i in range(len(autoCorrAnZ)):
    anzOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnZ[i],anZInt[i]))
anzOutWrite.close()

anXDriftOutFile = "Anion_X_Drift.txt"
anxDriftOutWrite = open(anXDriftOutFile,"w")
anxDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrAnX)):
    anxDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnXDriftCorr[i],anXDriftInt[i]))
anxDriftOutWrite.close()

anYDriftOutFile = "Anion_Y_Drift.txt"
anyDriftOutWrite = open(anYDriftOutFile,"w")
anyDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrAnY)):
    anyDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnYDriftCorr[i],anYDriftInt[i]))
anyDriftOutWrite.close()

anZDriftOutFile = "Anion_Z_drift.txt"
anzDriftOutWrite = open(anZDriftOutFile,"w")
anzDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrAnZ)):
    anzDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrAnZDriftCorr[i],anZDriftInt[i]))
anzDriftOutWrite.close()









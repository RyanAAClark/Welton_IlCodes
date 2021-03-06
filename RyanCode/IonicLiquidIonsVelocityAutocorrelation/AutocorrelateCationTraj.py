'''
Takes the cation centre of mass trajectories generated by "SplitIntoComVelFiles.py" and autocorrelates 
along the x, y, and z vectors using the different trajectories.

Complimentary script to "AutocorrelateAnionTraj.py"

Also calculates the drift corrected autocorrelation along these vectors.
Each autocorrelation is integrated and saved to a file.

THIS IS A VERY TIME CONSUMING SCRIPT!!!
Often it is better to only autocorrelate one direction at a time, even with multithreading.
To do this, comment out the analysis lines for all vectors except the one desired.
'''

import numpy as np
import multiprocessing

# Autocorrelate vector, normalise and cut to half the length 
# (one side of autocorr is all that's needed for this purpose)
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    resultSize = int(result.size/2)
    result = result[resultSize:]
    normRange = range(len(result),0,-1)
    result = result/normRange
    return result

def runAutoCorrX(ionNo):
    print("Autocorrelating cation in X: %i of %i" %(ionNo+1, noIons))
    autoCorrTempX = autocorr(cationXVel[:,ionNo])
    return autoCorrTempX
def runAutoCorrY(ionNo):
    print("Autocorrelating cation in Y: %i of %i" %(ionNo+1, noIons))
    autoCorrTempY = autocorr(cationYVel[:,ionNo])
    return autoCorrTempY
def runAutoCorrZ(ionNo):
    print("Autocorrelating cation in Z: %i of %i" %(ionNo+1, noIons))
    autoCorrTempZ = autocorr(cationZVel[:,ionNo])
    return autoCorrTempZ
    
def runAutoCorrXDrift(ionNo):
    print("Autocorrelating drift corrected cation in X: %i of %i" %(ionNo+1, noIons))
    autoCorrTempX = autocorr(cationXVelDriftCorr[:,ionNo])
    return autoCorrTempX
def runAutoCorrYDrift(ionNo):
    print("Autocorrelating drift corrected cation in Y: %i of %i" %(ionNo+1, noIons))
    autoCorrTempY = autocorr(cationYVelDriftCorr[:,ionNo])
    return autoCorrTempY
def runAutoCorrZDrift(ionNo):
    print("Autocorrelating drift corrected cation in Z: %i of %i" %(ionNo+1, noIons))
    autoCorrTempZ = autocorr(cationZVelDriftCorr[:,ionNo])
    return autoCorrTempZ

timestep = 1.0 #timestep of trajectory in fs, greater than 1 fs not recommended
noIons = 256
tSteps = 1000001

# Import velocities and calculate average velocity along them for drift correction
cationXVel = np.genfromtxt("CationXTrajectory.txt",delimiter=";"); cationXAvgVel = np.mean(cationXVel)
cationYVel = np.genfromtxt("CationYTrajectory.txt",delimiter=";"); cationYAvgVel = np.mean(cationYVel)
cationZVel = np.genfromtxt("CationZTrajectory.txt",delimiter=";"); cationZAvgVel = np.mean(cationZVel)

# Calculate drift corrected velocity arrays
cationXVelDriftCorr = cationXVel - cationXAvgVel
cationYVelDriftCorr = cationYVel - cationYAvgVel
cationZVelDriftCorr = cationZVel - cationZAvgVel

# Pre-allocate autocorrelation vectors
autoCorrCatX = np.zeros([tSteps,noIons])
autoCorrCatY = np.zeros([tSteps,noIons])
autoCorrCatZ = np.zeros([tSteps,noIons])
autoCorrCatXDriftCorr = np.zeros([tSteps,noIons])
autoCorrCatYDriftCorr = np.zeros([tSteps,noIons])
autoCorrCatZDriftCorr = np.zeros([tSteps,noIons])

if __name__ == '__main__': # This if statement is needed for multicore analysis in python, don't ask.....

    diffpool = multiprocessing.Pool(processes=60) # Change processes to the number of cores to run analysis on
    
    autoCorrCatX = diffpool.map(runAutoCorrX, range(noIons))
    autoCorrCatY = diffpool.map(runAutoCorrY, range(noIons))
    autoCorrCatZ = diffpool.map(runAutoCorrZ, range(noIons))
    
    autoCorrCatXDriftCorr = diffpool.map(runAutoCorrXDrift, range(noIons))
    autoCorrCatYDriftCorr = diffpool.map(runAutoCorrYDrift, range(noIons))
    autoCorrCatZDriftCorr = diffpool.map(runAutoCorrZDrift, range(noIons))

# Autocorrelation result has one column per ion, this averages the autocorrelation across all ions
autoCorrCatX = np.sum(autoCorrCatX, axis=0);autoCorrCatX = autoCorrCatX/noIons
autoCorrCatY = np.sum(autoCorrCatY, axis=0);autoCorrCatY = autoCorrCatY/noIons
autoCorrCatZ = np.sum(autoCorrCatZ, axis=0);autoCorrCatZ = autoCorrCatZ/noIons

autoCorrCatXDriftCorr = np.sum(autoCorrCatXDriftCorr, axis=0);autoCorrCatXDriftCorr = autoCorrCatXDriftCorr/noIons
autoCorrCatYDriftCorr = np.sum(autoCorrCatYDriftCorr, axis=0);autoCorrCatYDriftCorr = autoCorrCatYDriftCorr/noIons
autoCorrCatZDriftCorr = np.sum(autoCorrCatZDriftCorr, axis=0);autoCorrCatZDriftCorr = autoCorrCatZDriftCorr/noIons

# Calculate integral of autocorrelation
catXInt = np.cumsum(autoCorrCatX*timestep/1000) # Factor of 1000 converts timestep from fs to ps. 
catYInt = np.cumsum(autoCorrCatY*timestep/1000) # SplitIntoComVelFiles.py already converted distance into pm
catZInt = np.cumsum(autoCorrCatZ*timestep/1000)
catXDriftInt = np.cumsum(autoCorrCatXDriftCorr*timestep/1000)
catYDriftInt = np.cumsum(autoCorrCatYDriftCorr*timestep/1000)
catZDriftInt = np.cumsum(autoCorrCatZDriftCorr*timestep/1000)

# Save autocorrelations to file
catXOutFile = "Cation_X.txt"
catxOutWrite = open(catXOutFile,"w")
catxOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatX)):
    catxOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatX[i],catXInt[i]))
catxOutWrite.close()

catYOutFile = "Cation_Y.txt"
catyOutWrite = open(catYOutFile,"w")
catyOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatY)):
    catyOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatY[i],catYInt[i]))
catyOutWrite.close()

catZOutFile = "Cation_Z.txt"
catzOutWrite = open(catZOutFile,"w")
catzOutWrite.write("Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatZ)):
    catzOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatZ[i],catZInt[i]))
catzOutWrite.close()

catXDriftOutFile = "Cation_X_Drift.txt"
catxDriftOutWrite = open(catXDriftOutFile,"w")
catxDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatX)):
    catxDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatXDriftCorr[i],catXDriftInt[i]))
catxDriftOutWrite.close()

catYDriftOutFile = "Cation_Y_Drift.txt"
catyDriftOutWrite = open(catYDriftOutFile,"w")
catyDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatY)):
    catyDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatYDriftCorr[i],catYDriftInt[i]))
catyDriftOutWrite.close()

catZDriftOutFile = "Cation_Z_drift.txt"
catzDriftOutWrite = open(catZDriftOutFile,"w")
catzDriftOutWrite.write("#Time [fs]; Velocity Autocorrelation [pm^2/ps^2]; Integral [pm^2/ps]\n")
for i in range(len(autoCorrCatZ)):
    catzDriftOutWrite.write("%f;\t%f;\t%f\n" %(i*timestep,autoCorrCatZDriftCorr[i],catZDriftInt[i]))
catzDriftOutWrite.close()


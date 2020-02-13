'''
Calculates the Mean Squared Displacement (MSD), Mean Displacement (MD), and Drift Corrected Mean Squared Displacement (Drift)
of an atom (listed in "atomsToLookFor") in a trajectory.

MD = <x-x0>
MSD = <(x-x0)^2>
Drift = <(x-x0)^2> - <x-x0>^2

<> denotes average across ensemble

Will start up "noProcessers" processes and run on as many processers as it can find.
CARE THAT THIS IS NOT HIGHER THAN THE NUMBER OF PROCESSORS AVAILABLE. THIS WILL SEVERELY AFFECT THE PERFORMANCE NEGATIVELY

Input here is a LAMMPS trajectory file, that has been processed using the PROC function in TRAVIS to save only the atom that
is desired to be studied. 

Currently calculates the MSD, MD and Drift along the x, y, z vectors, as well as the total MSD, MD, Drift.

Uses a manual sliding window technique to calculate the mean displacement and mean square displacement at every lag time
    - Lag time meaning every time seperated by the lag time i.e. time=t -> time=t+lag
    - Range of lag times used are 0-"maxLagTime"
'''

# Import packages
import numpy as np
import sys
import multiprocessing
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Variables to change
File = "trajectory_out.xyz" #Trajectory file containing positions of atoms seperated by the same time
atomsToLookFor = ['#3'] #Atom label to pull out of file for analysis 
                        #WILL ASSUME ALL LABELS CORRESPOND TO THE SAME ATOM IF MULTIPLE DEFINED THIS CODE WILL NOT WORK
nTimeSteps = 20001 #Total number of timesteps in the trajectory file. If unknown, set to a large value and the script will cut to size
timeStep = 500 #Timestep in the file in fs
maxLagTime = 18000 #Maximum size of the sliding window to use to analyse. Recommended is >0.75*nTimeSteps. 
                   #If timesteps are changed, this may need redefining, which will be asked automatically.
noAtoms = 256 #Number of atoms in trajectory to process
noProcessers = 48 #Number of processers available to run analysis on

# Begin timer
sTime = time.time()
print("Start time: " +  time.asctime(time.localtime(time.time())))

## Define fitting functions for data
def Einstein1D(x, a, b):
    return 2*a*x+b
def Einstein3D(x, a, b):
    return 6*a*x+b
def linear(x, a, b):
    return a*x+b

# Define R^2 fitting error function
def rsquared(rawyData,calcyData):
    yBar = np.mean(rawyData)
    sumOfSquareError = 0
    sumOfSquareMean = 0
    for dataPoint in range(len(rawyData)):
        yPoint = rawyData[dataPoint]
        yEvaluated = calcyData[dataPoint]
        sumOfSquareError+=(yPoint-yEvaluated)**2
        sumOfSquareMean+=(yPoint-yBar)**2
    rSquared = 1-(sumOfSquareError/sumOfSquareMean)
    return rSquared

# Define function to calculate MSD and MD
def atomMdMsd(n):
    ssTime = time.time()
    atNoForAnal = n%noAtoms #Determine atom number to analyse
    dirNo = int(n/noAtoms) #Give a number for direction to analyse
    # Pull out the appropriate vectors for analysis depending on the direction number put in.
    if dirNo==0:
        atom = xPos[atNoForAnal,:]
        print("Atom No: %i | X" %(atNoForAnal))  
    elif dirNo==1:
        atom = yPos[atNoForAnal,:]
        print("Atom No: %i | Y" %(atNoForAnal)) 
    elif dirNo==2:
        atom = zPos[atNoForAnal,:]
        print("Atom No: %i | Z" %(atNoForAnal)) 
    elif dirNo==3:
        atomX = xPos[atNoForAnal,:]
        atomY = yPos[atNoForAnal,:]
        atomZ = zPos[atNoForAnal,:]
        print("Atom No: %i | Tot" %(atNoForAnal))
    
    # Pre-allocate cumulative sum arrays for each delay time
    MDCumSum=np.zeros(maxLagTime+1)
    MSDCumSum=np.zeros(maxLagTime+1)
    Count=np.zeros(maxLagTime+1)
    
    # For every delay (deltaT), calculate the total MSD and MD over all sliding windows
    for deltaT in range(maxLagTime+1):
        if deltaT==0:
            MDCumSum[deltaT]+= 0.
            MSDCumSum[deltaT]+= 0.
            Count[deltaT]=1
        elif deltaT<(maxLagTime+1):
            arrayPos=1
            mdArrayCumSum=0
            msdArrayCumSum=0
            while arrayPos<nTimeSteps-deltaT+1:
                
                if dirNo==3:
                    MDtempx = atomX[-arrayPos]-atomX[-(arrayPos+deltaT)] 
                    MDtempy = atomY[-arrayPos]-atomY[-(arrayPos+deltaT)] 
                    MDtempz = atomZ[-arrayPos]-atomZ[-(arrayPos+deltaT)]
                    
                    DISTANCE = (MDtempx**2+MDtempy**2+MDtempz**2)**0.5
                    MSDFROMDISTANCE = (MDtempx**2+MDtempy**2+MDtempz**2)
                    
                else:
                    MDtemp = atom[-arrayPos]-atom[-(arrayPos+deltaT)]
                
                if dirNo==3:
                    mdArrayCumSum+=DISTANCE
                    msdArrayCumSum+=MSDFROMDISTANCE
                else:
                    mdArrayCumSum = mdArrayCumSum+MDtemp
                    msdArrayCumSum = msdArrayCumSum+MDtemp**2
                
                Count[deltaT]+=1 
                arrayPos+=1
            
            MDCumSum[deltaT] = mdArrayCumSum           
            MSDCumSum[deltaT] = msdArrayCumSum

    # Cut any trailing zeros from array
    MDCumSum=MDCumSum[0:deltaT+1]
    MSDCumSum=MSDCumSum[0:deltaT+1]
    Count=Count[0:deltaT+1]
    # Normalise by the number of points
    MD = MDCumSum/Count
    MSD = MSDCumSum/Count 
    
    if dirNo==0:
        print("Atom No: %i | X | Time = %.2f s" %(atNoForAnal,time.time()-ssTime))  
    elif dirNo==1:
        print("Atom No: %i | Y | Time = %.2f s" %(atNoForAnal,time.time()-ssTime)) 
    elif dirNo==2:
        print("Atom No: %i | Z | Time = %.2f s" %(atNoForAnal,time.time()-ssTime)) 
    elif dirNo==3:
        print("Atom No: %i | Tot | Time = %.2f s" %(atNoForAnal,time.time()-ssTime))
    
    return MD,MSD

tStep=0
linecounter=0
xPos=np.zeros([noAtoms,nTimeSteps])
yPos=np.zeros([noAtoms,nTimeSteps])
zPos=np.zeros([noAtoms,nTimeSteps])

# Import trajecotry into 
for line in open(File):
    split=line.split()
    if split[0] in atomsToLookFor:
        AtomNo = linecounter%noAtoms
        xPos[AtomNo,tStep]=float(split[1])*100. # Factor of 100 turns from angstroms to pm
        yPos[AtomNo,tStep]=float(split[2])*100.
        zPos[AtomNo,tStep]=float(split[3])*100.
        linecounter+=1
        
        if AtomNo==noAtoms-1:
            tStep+=1
            sys.stdout.write("\r%i of %i" %(linecounter/noAtoms,nTimeSteps))
            sys.stdout.flush()

sys.stdout.write("\n")
sys.stdout.flush()

# Checks if number of timesteps is correct, if not cuts to appropriate size. Also allows re-entry of lagtime if too high.
if tStep<nTimeSteps:
    print("Actual number of timesteps less than stated value, cutting to size")
    print("Actual number of timesteps = %i" %(tStep))
    xPos = xPos[:,0:tStep]
    yPos = yPos[:,0:tStep]
    zPos = zPos[:,0:tStep]
    nTimeSteps = tStep
    while maxLagTime>nTimeSteps:
        print("Changing the number of timesteps means the lag time is invalid")
        print("Old lag time = %i"%(maxLagTime))
        maxLagTime = input("Please enter a new lag time less than %i: " %(nTimeSteps))

# Normalise x, y, z vectors into distance travelled (not needed, but easier to handle when calculating 3D diffusion)
for ATNO in range(noAtoms):
    xPos[ATNO,:]=xPos[ATNO,:]-xPos[ATNO,0]
    yPos[ATNO,:]=yPos[ATNO,:]-yPos[ATNO,0]
    zPos[ATNO,:]=zPos[ATNO,:]-zPos[ATNO,0]

if __name__ == '__main__':

    diffpool = multiprocessing.Pool(processes=noProcessers)
    
    # Calculate Md and Msd of the atoms. Columns 0-255 are x data, 256-511 y data, 512-767 z data, 768-1023 Total data
    # This is sorted out later into the correct Md and Msd data, labelled by the end part of the array name
    diffAll = diffpool.map(atomMdMsd, range(noAtoms*4))

# Create time array
timeStep = timeStep/1000. # Convert timeStep to ps
timeArray = np.arange(0.,(maxLagTime+1)*timeStep,timeStep) 

# Pre-allocate arrays for sorting
MdX = np.zeros([noAtoms,maxLagTime+1])
MdY = np.zeros([noAtoms,maxLagTime+1])
MdZ = np.zeros([noAtoms,maxLagTime+1])
MdTot = np.zeros([noAtoms,maxLagTime+1])
MsdX = np.zeros([noAtoms,maxLagTime+1])
MsdY = np.zeros([noAtoms,maxLagTime+1])
MsdZ = np.zeros([noAtoms,maxLagTime+1])
MsdTot = np.zeros([noAtoms,maxLagTime+1])

for i in range(noAtoms):
    # Squashes array into Md and Msd per atom and per direction
    tempX = diffAll[i][:][:]
    tempY = diffAll[i+256][:][:]
    tempZ = diffAll[i+512][:][:]
    tempTot = diffAll[i+768][:][:]
    # Extract and save Md per atom
    MdTempX = tempX[0][:]
    MdTempY = tempY[0][:]
    MdTempZ = tempZ[0][:]
    MdTempTot = tempTot[0][:]
    MdX[i,:] = MdTempX
    MdY[i,:] = MdTempY
    MdZ[i,:] = MdTempZ
    MdTot[i,:] = MdTempTot
    # Extract and save Msd per atom
    MsdTempX = tempX[1][:]
    MsdTempY = tempY[1][:]
    MsdTempZ = tempZ[1][:]
    MsdTempTot = tempTot[1][:]
    MsdX[i,:] = MsdTempX
    MsdY[i,:] = MsdTempY
    MsdZ[i,:] = MsdTempZ
    MsdTot[i,:] = MsdTempTot

# Average Md and Msd over all atoms
MdX = np.average(MdX, axis=0)
MdY = np.average(MdY, axis=0)
MdZ = np.average(MdZ, axis=0)
MdTot = np.average(MdTot, axis=0)
MsdX = np.average(MsdX, axis=0)
MsdY = np.average(MsdY, axis=0)
MsdZ = np.average(MsdZ, axis=0)
MsdTot = np.average(MsdTot, axis=0)

# Calculate drift from MSD and MD
DriftX = MsdX-MdX**2
DriftY = MsdY-MdY**2
DriftZ = MsdZ-MdZ**2
DriftTot = MsdTot-MdTot**2

## Fit data ##

# Determine indices to use for fitting
FitSt = int(np.ceil(0.25*maxLagTime))
FitEd = int(np.ceil(0.75*maxLagTime))
    
# Fit Msd Data
msdXFit,msdXCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],MsdX[FitSt:FitEd])
msdYFit,msdYCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],MsdY[FitSt:FitEd])
msdZFit,msdZCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],MsdZ[FitSt:FitEd])
msdTotFit,msdTotCov = curve_fit(Einstein3D,timeArray[FitSt:FitEd],MsdTot[FitSt:FitEd])
# Msd error
msdXError = np.diag(msdXCov)
msdYError = np.diag(msdYCov)
msdZError = np.diag(msdZCov)
msdTotError = np.diag(msdTotCov)
# R Squared
msdXFitRSq = rsquared(MsdX[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],msdXFit[0],msdXFit[1]))
msdXTotRSq = rsquared(MsdX,Einstein1D(timeArray,msdXFit[0],msdXFit[1]))
msdYFitRSq = rsquared(MsdY[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],msdYFit[0],msdYFit[1]))
msdYTotRSq = rsquared(MsdY,Einstein1D(timeArray,msdYFit[0],msdYFit[1]))
msdZFitRSq = rsquared(MsdZ[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],msdZFit[0],msdZFit[1]))
msdZTotRSq = rsquared(MsdZ,Einstein1D(timeArray,msdZFit[0],msdZFit[1]))
msdTotFitRSq = rsquared(MsdTot[FitSt:FitEd],Einstein3D(timeArray[FitSt:FitEd],msdTotFit[0],msdTotFit[1]))
msdTotTotRSq = rsquared(MsdTot,Einstein3D(timeArray,msdTotFit[0],msdTotFit[1]))

# Print MSD output
print("\n")
print("Msd Analysis----------------------------------------")
print("X direction------")
print("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s" %(msdXFit[0],msdXError[0],msdXFit[0]*(10**-12)))
print("R Squared for fit: %f" %(msdXFitRSq))
print("R Squared for total graph: %f" %(msdXTotRSq))

print("Y direction------")
print("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s" %(msdYFit[0],msdYError[0],msdYFit[0]*(10**-12)))
print("R Squared for fit: %f" %(msdYFitRSq))
print("R Squared for total graph: %f" %(msdYTotRSq))

print("Z direction------")
print("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s" %(msdZFit[0],msdZError[0],msdZFit[0]*(10**-12)))
print("R Squared for fit: %f" %(msdZFitRSq))
print("R Squared for total graph: %f" %(msdZTotRSq))

print("All directions---")
print("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s" %(msdTotFit[0],msdTotError[0],msdTotFit[0]*(10**-12)))
print("R Squared for fit: %f" %(msdTotFitRSq))
print("R Squared for total graph: %f" %(msdTotTotRSq))
print("\n")

# Fit Md Data
mdXFit,mdXCov = curve_fit(linear,timeArray[FitSt:FitEd],MdX[FitSt:FitEd])
mdYFit,mdYCov = curve_fit(linear,timeArray[FitSt:FitEd],MdY[FitSt:FitEd])
mdZFit,mdZCov = curve_fit(linear,timeArray[FitSt:FitEd],MdZ[FitSt:FitEd])
mdTotFit,mdTotCov = curve_fit(linear,timeArray[FitSt:FitEd],MdTot[FitSt:FitEd])
# Md error
mdXError = np.diag(mdXCov)
mdYError = np.diag(mdYCov)
mdZError = np.diag(mdZCov)
mdTotError = np.diag(mdTotCov)
# R Squared
mdXFitRSq = rsquared(MdX[FitSt:FitEd],linear(timeArray[FitSt:FitEd],mdXFit[0],mdXFit[1]))
mdXTotRSq = rsquared(MdX,linear(timeArray,mdXFit[0],mdXFit[1]))
mdYFitRSq = rsquared(MdY[FitSt:FitEd],linear(timeArray[FitSt:FitEd],mdYFit[0],mdYFit[1]))
mdYTotRSq = rsquared(MdY,linear(timeArray,mdYFit[0],mdYFit[1]))
mdZFitRSq = rsquared(MdZ[FitSt:FitEd],linear(timeArray[FitSt:FitEd],mdZFit[0],mdZFit[1]))
mdZTotRSq = rsquared(MdZ,linear(timeArray,mdZFit[0],mdZFit[1]))
mdTotFitRSq = rsquared(MdTot[FitSt:FitEd],linear(timeArray[FitSt:FitEd],mdTotFit[0],mdTotFit[1]))
mdTotTotRSq = rsquared(MdTot,linear(timeArray,mdTotFit[0],mdTotFit[1]))

# Print MD output
print("Md Analysis-----------------------------------------")
print("X direction------")
print("Drift Velocity: %f +- %f pm/ps  =  %e m/s" %(mdXFit[0],mdXError[0],mdXFit[0]*(10**-12)))
print("R Squared for fit: %f" %(mdXFitRSq))
print("R Squared for total graph: %f" %(mdXTotRSq))

print("Y direction------")
print("Drift Velocity: %f +- %f pm/ps  =  %e m/s" %(mdYFit[0],mdYError[0],mdYFit[0]*(10**-12)))
print("R Squared for fit: %f" %(mdYFitRSq))
print("R Squared for total graph: %f" %(mdYTotRSq))

print("Z direction------")
print("Drift Velocity: %f +- %f pm/ps  =  %e m/s" %(mdZFit[0],mdZError[0],mdZFit[0]*(10**-12)))
print("R Squared for fit: %f" %(mdZFitRSq))
print("R Squared for total graph: %f" %(mdZTotRSq))

print("All directions---")
print("Drift Velocity: %f +- %f pm/ps  =  %e m/s" %(mdTotFit[0],mdTotError[0],mdTotFit[0]*(10**-12)))
print("R Squared for fit: %f" %(mdTotFitRSq))
print("R Squared for total graph: %f" %(mdTotTotRSq))
print("\n")

# Fit drift corrected Msd data
driftCorrXFit,driftCorrXCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],DriftX[FitSt:FitEd])
driftCorrYFit,driftCorrYCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],DriftY[FitSt:FitEd])
driftCorrZFit,driftCorrZCov = curve_fit(Einstein1D,timeArray[FitSt:FitEd],DriftZ[FitSt:FitEd])
driftCorrTotFit,driftCorrTotCov = curve_fit(Einstein3D,timeArray[FitSt:FitEd],DriftTot[FitSt:FitEd])
# Drift error
driftCorrXError = np.diag(driftCorrXCov)
driftCorrYError = np.diag(driftCorrYCov)
driftCorrZError = np.diag(driftCorrZCov)
driftCorrTotError = np.diag(driftCorrTotCov)
# R Squared
driftCorrXFitRSq = rsquared(DriftX[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],driftCorrXFit[0],driftCorrXFit[1]))
driftCorrXTotRSq = rsquared(DriftX,Einstein1D(timeArray,driftCorrXFit[0],driftCorrXFit[1]))
driftCorrYFitRSq = rsquared(DriftY[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],driftCorrYFit[0],driftCorrYFit[1]))
driftCorrYTotRSq = rsquared(DriftY,Einstein1D(timeArray,driftCorrYFit[0],driftCorrYFit[1]))
driftCorrZFitRSq = rsquared(DriftZ[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],driftCorrZFit[0],driftCorrZFit[1]))
driftCorrZTotRSq = rsquared(DriftZ,Einstein1D(timeArray,driftCorrZFit[0],driftCorrZFit[1]))
driftCorrTotFitRSq = rsquared(DriftTot[FitSt:FitEd],Einstein1D(timeArray[FitSt:FitEd],driftCorrTotFit[0],driftCorrTotFit[1]))
driftCorrTotTotRSq = rsquared(DriftTot,Einstein1D(timeArray,driftCorrTotFit[0],driftCorrTotFit[1]))

# Print drift corrected MSD output
print("Effective Diffusion Analysis------------------------")
print("X direction------")
print("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s" %(driftCorrXFit[0],driftCorrXError[0],driftCorrXFit[0]*(10**-12)))
print("R Squared for fit: %f" %(driftCorrXFitRSq))
print("R Squared for total graph: %f" %(driftCorrXTotRSq))

print("Y direction------")
print("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s" %(driftCorrYFit[0],driftCorrYError[0],driftCorrYFit[0]*(10**-12)))
print("R Squared for fit: %f" %(driftCorrYFitRSq))
print("R Squared for total graph: %f" %(driftCorrYTotRSq))

print("Z direction------")
print("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s" %(driftCorrZFit[0],driftCorrZError[0],driftCorrZFit[0]*(10**-12)))
print("R Squared for fit: %f" %(driftCorrZFitRSq))
print("R Squared for total graph: %f" %(driftCorrZTotRSq))

print("All directions---")
print("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s" %(driftCorrTotFit[0],driftCorrTotError[0],driftCorrTotFit[0]*(10**-12)))
print("R Squared for fit: %f" %(driftCorrTotFitRSq))
print("R Squared for total graph: %f" %(driftCorrTotTotRSq))
print("\n")

# Save MSD raw data and fit
MsdXOutName = "Msd_X.txt"; MsdXOut=open(MsdXOutName,"w"); MsdXOut.write("# Time[ps]; MSD [pm^2]; Fit\n")
MsdYOutName = "Msd_Y.txt"; MsdYOut=open(MsdYOutName,"w"); MsdYOut.write("# Time[ps]; MSD [pm^2]; Fit\n")
MsdZOutName = "Msd_Z.txt"; MsdZOut=open(MsdZOutName,"w"); MsdZOut.write("# Time[ps]; MSD [pm^2]; Fit\n")
MsdTotOutName = "Msd_Tot.txt"; MsdTotOut=open(MsdTotOutName,"w"); MsdTotOut.write("# Time[ps]; MSD [pm^2]; Fit\n")

for i in range(len(timeArray)):
    MsdXOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MsdX[i], Einstein1D(timeArray[i], msdXFit[0], msdXFit[1])))
    MsdYOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MsdY[i], Einstein1D(timeArray[i], msdYFit[0], msdYFit[1])))
    MsdZOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MsdZ[i], Einstein1D(timeArray[i], msdZFit[0], msdZFit[1])))
    MsdTotOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MsdTot[i], Einstein3D(timeArray[i], msdTotFit[0], msdTotFit[1])))

MsdXOut.close(); MsdYOut.close(); MsdZOut.close(); MsdTotOut.close()

# Save MD raw data and fit
MdXOutName = "Md_X.txt"; MdXOut=open(MdXOutName,"w"); MdXOut.write("# Time[ps]; MD [pm]; Fit\n")
MdYOutName = "Md_X.txt"; MdYOut=open(MdYOutName,"w"); MdYOut.write("# Time[ps]; MD [pm]; Fit\n")
MdZOutName = "Md_X.txt"; MdZOut=open(MdZOutName,"w"); MdZOut.write("# Time[ps]; MD [pm]; Fit\n")
MdTotOutName = "Md_X.txt"; MdTotOut=open(MdTotOutName,"w"); MdTotOut.write("# Time[ps]; MD [pm]; Fit\n")

for i in range(len(timeArray)):
    MdXOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MdX[i], linear(timeArray[i], mdXFit[0], mdXFit[1])))
    MdYOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MdY[i], linear(timeArray[i], mdYFit[0], mdYFit[1])))
    MdZOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MdZ[i], linear(timeArray[i], mdZFit[0], mdZFit[1])))
    MdTotOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], MdTot[i], linear(timeArray[i], mdTotFit[0], mdTotFit[1])))

MdXOut.close(); MdYOut.close(); MdZOut.close(); MdTotOut.close()

# Save drift-corrected diffusion raw data and fit
DriftXOutName = "Drift_X.txt"; DriftXOut=open(DriftXOutName,"w"); DriftXOut.write("# Time[ps]; Drift-corrected MSD [pm^2]; Fit\n")
DriftYOutName = "Drift_Y.txt"; DriftYOut=open(DriftYOutName,"w"); DriftYOut.write("# Time[ps]; Drift-corrected MSD [pm^2]; Fit\n")
DriftZOutName = "Drift_Z.txt"; DriftZOut=open(DriftZOutName,"w"); DriftZOut.write("# Time[ps]; Drift-corrected MSD [pm^2]; Fit\n")
DriftTotOutName = "Drift_Tot.txt"; DriftTotOut=open(DriftTotOutName,"w"); DriftTotOut.write("# Time[ps]; Drift-corrected MSD [pm^2]; Fit\n")

for i in range(len(timeArray)):
    DriftXOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], DriftX[i], Einstein1D(timeArray[i], driftCorrXFit[0], driftCorrXFit[1])))
    DriftYOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], DriftY[i], Einstein1D(timeArray[i], driftCorrYFit[0], driftCorrYFit[1])))
    DriftZOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], DriftZ[i], Einstein1D(timeArray[i], driftCorrZFit[0], driftCorrZFit[1])))
    DriftTotOut.write("%10.3f;\t%12.6f;\t%12.6f\n" %(timeArray[i], DriftTot[i], Einstein3D(timeArray[i], driftCorrTotFit[0], driftCorrTotFit[1])))

DriftXOut.close(); DriftYOut.close(); DriftZOut.close(); DriftTotOut.close()

# Print a documant with the summary, which was output to the terminal
OutName = "Summary.txt"
Out=open(OutName,"w")
Out.write("SUMMARY OF RESULTS\n\n")

Out.write("Msd Analysis----------------------------------------\n")
Out.write("X direction------\n")
Out.write("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s\n" %(msdXFit[0],msdXError[0],msdXFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(msdXFitRSq))
Out.write("R Squared for total graph: %f\n" %(msdXTotRSq))
Out.write("\n")
Out.write("Y direction------\n")
Out.write("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s\n" %(msdYFit[0],msdYError[0],msdYFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(msdYFitRSq))
Out.write("R Squared for total graph: %f\n" %(msdYTotRSq))
Out.write("\n")
Out.write("Z direction------\n")
Out.write("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s\n" %(msdZFit[0],msdZError[0],msdZFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(msdZFitRSq))
Out.write("R Squared for total graph: %f\n" %(msdZTotRSq))
Out.write("\n")
Out.write("All directions---\n")
Out.write("Diffusion Coefficient: %f +- %f pm^2/ps  =  %e m^2/s\n" %(msdTotFit[0],msdTotError[0],msdTotFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(msdTotFitRSq))
Out.write("R Squared for total graph: %f\n" %(msdTotTotRSq))
Out.write("\n\n")

Out.write("Md Analysis-----------------------------------------\n")
Out.write("X direction------\n")
Out.write("Drift Velocity: %f +- %f pm/ps  =  %e m/s\n" %(mdXFit[0],mdXError[0],mdXFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(mdXFitRSq))
Out.write("R Squared for total graph: %f\n" %(mdXTotRSq))
Out.write("\n")
Out.write("Y direction------\n")
Out.write("Drift Velocity: %f +- %f pm/ps  =  %e m/s\n" %(mdYFit[0],mdYError[0],mdYFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(mdYFitRSq))
Out.write("R Squared for total graph: %f\n" %(mdYTotRSq))
Out.write("\n")
Out.write("Z direction------\n")
Out.write("Drift Velocity: %f +- %f pm/ps  =  %e m/s\n" %(mdZFit[0],mdZError[0],mdZFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(mdZFitRSq))
Out.write("R Squared for total graph: %f\n" %(mdZTotRSq))
Out.write("\n")
Out.write("All directions---\n")
Out.write("Drift Velocity: %f +- %f pm/ps  =  %e m/s\n" %(mdTotFit[0],mdTotError[0],mdTotFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(mdTotFitRSq))
Out.write("R Squared for total graph: %f\n" %(mdTotTotRSq))
Out.write("\n\n")

Out.write("Effective Diffusion Analysis------------------------\n")
Out.write("X direction------\n")
Out.write("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s\n" %(driftCorrXFit[0],driftCorrXError[0],driftCorrXFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(driftCorrXFitRSq))
Out.write("R Squared for total graph: %f\n" %(driftCorrXTotRSq))
Out.write("\n")
Out.write("Y direction------\n")
Out.write("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s\n" %(driftCorrYFit[0],driftCorrYError[0],driftCorrYFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(driftCorrYFitRSq))
Out.write("R Squared for total graph: %f\n" %(driftCorrYTotRSq))
Out.write("\n")
Out.write("Z direction------\n")
Out.write("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s\n" %(driftCorrZFit[0],driftCorrZError[0],driftCorrZFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(driftCorrZFitRSq))
Out.write("R Squared for total graph: %f\n" %(driftCorrZTotRSq))
Out.write("\n")
Out.write("All directions---\n")
Out.write("Effective Diffusion Coefficient: %f +- %f p^2m/ps  =  %e m^2/s\n" %(driftCorrTotFit[0],driftCorrTotError[0],driftCorrTotFit[0]*(10**-12)))
Out.write("R Squared for fit: %f\n" %(driftCorrTotFitRSq))
Out.write("R Squared for total graph: %f\n" %(driftCorrTotTotRSq))
Out.write("\n\n")

Out.write("Total time for analysis: %10.2f seconds" %(time.time()-sTime))

# Print total time taken to terminal
eTime = time.time()
print("Finished. That took %10.2f seconds" %(eTime-sTime))

# Uncomment to show plots
#plt.plot(timeArray,MdX);plt.plot(timeArray,MdY);plt.plot(timeArray,MdZ);plt.plot(timeArray,MdTot);plt.show()

#plt.plot(timeArray,MsdX);plt.plot(timeArray,MsdY);plt.plot(timeArray,MsdZ);plt.plot(timeArray,MsdTot);plt.show()

#plt.plot(timeArray,DriftX);plt.plot(timeArray,DriftY);plt.plot(timeArray,DriftZ);plt.plot(timeArray,DriftTot);plt.show()

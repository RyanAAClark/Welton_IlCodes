'''
Takes a trajectory file and calculates the distribution of the length and angle of a vector defined by two specified 
markers in the file. The results are then saved to the text files XLength.txt, YLength.txt, ZLength.txt, XAngle.txt, 
YAngle.txt, ZAngle.txt.

The markers that can be used are up to two character strings that identify the origin of the vector and the end of the 
vector.
	NOTE: one character strings can be used, but the markers are required to always be 2 characters long, so they can 
	      be "S "., but not "S"
The switching variable "firstMarker" allows the two markers within the file to be the same.

Takes the file in the variable "File" and searches through it for the "originMarker" and "endMarker" labels. 
This is ideally a .xyz file, structured as such:
	noParticles
	Comment
	MARKER xPosition yPosition zPosition -|
	MARKER xPosition yPosition zPosition  |- repeated noParticles times	
	MARKER xPosition yPosition zPosition -|

These must be sequential, i.e. the file must alternate between originMarker for atom 1, endMarker for atom 1, 
originMarker for atom 2, endMarker for atom 2 etc. It does not matter if other markers are between them, as long as 
they occur in this order. 

The angle data is corrected using a cone correction, allowing a relative occurance similar to g(r), where 1 = uniform 
distribution.

The input "File" can be easily generated using the proc function in TRAVIS, outputting only the atoms (virtual or 
otherwise) desired.

Example:
	-This can be used to calculate the charge arm orientation of an ionic liquid.
	-The charge arm is the vector originating at the centre of mass, and ending at the centre of charge.
	-Therefore the "File" to be analysed contains the centre of mass (a virtual atom) as the originMarker flag, and the 
	 centre of charge (a virtual atom) as the endMarker flag.
	-A "File" should be generated using the proc funtion in TRAVIS for each ion, then this analysis run on each of 
	 those output trajectories.
'''
import numpy as np
import matplotlib.pyplot as plt

# Parameters to change
File = "trajectory_out.xyz"
originMarker = "#2" # Origin point of the vector
endMarker = "#3" # End point of the vector
noParticles = 256 # Number of particles per timestep in the file
nTimeSteps = 20001 # Number of timesteps present in the file
noBinsLength = 100 # Number of bins for the length of the vector array

# Angle calculation between two vectors
def angle(a, b):
    magA = np.linalg.norm(a)
    magB = np.linalg.norm(b)
    dotprod = np.dot(a,b)
    cosAng = dotprod/(magA*magB)
    Ang = np.arccos(cosAng)
    return Ang

## Normalisation for angles
maxAng=180
BigArea = np.zeros(maxAng)
SmallArea = np.zeros(maxAng)
AreaSlice = np.zeros(maxAng)
for angles in range(maxAng):
     AreaSlice[angles] = 2*np.pi*(1-np.pi*np.cos((angles+1)*np.pi/180))-2*np.pi*(1-np.pi*np.cos((angles)*np.pi/180))
 
# Pre-allocate vectors for analysis
xdist = np.zeros(nTimeSteps*noParticles)
ydist = np.zeros(nTimeSteps*noParticles)
zdist = np.zeros(nTimeSteps*noParticles)

xangles = np.zeros(nTimeSteps*noParticles)
yangles = np.zeros(nTimeSteps*noParticles)
zangles = np.zeros(nTimeSteps*noParticles)

atomCounter=0
firstMarker=True

for line in open(File):
    
    if line[:2]==originMarker and firstMarker:
		# Pull out position of the origin of the vector (should be first in list)
        xyz = line.split()
        ogX = float(xyz[1])*100. #Convert from angstrom to pm
        ogY = float(xyz[2])*100.
        ogZ = float(xyz[3])*100.
		firstMarker=False
        
    if line[:2]==endMarker and not firstMarker:
		# Pull out position of the end of the vector (should be second in list)
        xyz = line.split()
        endX = float(xyz[1])*100.
        endY = float(xyz[2])*100.
		endZ = float(xyz[3])*100.
		firstMarker=True
		
		# Calcuate the length of the charge arm in each direction
		xdist[atomCounter] = endX-ogX
		ydist[atomCounter] = endY-ogY
		zdist[atomCounter] = endZ-ogZ
		
		# Describe the charge arm with a single vector (easier for angle calculations)
		vector = np.array([endX-ogX],[endY-ogY],[endZ-ogZ])
		
		# Calculate the angle of the charge arm with respect to each vector of the box
        xangles[atomCounter] = angle(vector,[1,0,0])
		yangles[atomCounter] = angle(vector,[0,1,0])
		zangles[atomCounter] = angle(vector,[0,0,1])
		
        atomCounter+=1
        
        if int(atomCounter/noParticles)*noParticles == atomCounter and int(atomCounter/noParticles)==1:
            sys.stdout.write("%i of %i" %(int(atomCounter/noParticles),nTimeSteps))
        elif int(atomCounter/noParticles)*noParticles == atomCounter and int(atomCounter/noParticles)==nTimeSteps:
            sys.stdout.flush()
            sys.stdout.write("\r%i of %i\n" %(int(atomCounter/noParticles),nTimeSteps))
        elif int(atomCounter/noParticles)*noParticles == atomCounter:
            sys.stdout.flush()
            sys.stdout.write("\r%i of %i" %(int(atomCounter/noParticles),nTimeSteps))

# Calculate length histograms        
histX = np.histogram(xdist, bins=noBinsLength, range=None, normed=None, weights=None, density=False)
histY = np.histogram(ydist, bins=noBinsLength, range=None, normed=None, weights=None, density=False)
histZ = np.histogram(zdist, bins=noBinsLength, range=None, normed=None, weights=None, density=False)

# Collect bin edges (range) and no values in each bin (noPoints)
rangesX = histX[1];noPointsX = histX[0]
rangesY = histY[1];noPointsY = histY[0]
rangesZ = histZ[1];noPointsZ = histZ[0]

# Save length histograms to file
f = open("XLength.txt","w")
f.write("Length [pm];\tHistogram [counts]\n")
for i in range(len(noPointsX)):
    f.write("%f\t%f\n" %(np.mean([rangesX[i],rangesX[i+1]]),noPointsX[i]))
f.close()

f = open("YLength.txt","w")
f.write("Length [pm];\tHistogram [counts]\n")
for i in range(len(noPointsY)):
    f.write("%f\t%f\n" %(np.mean([rangesY[i],rangesY[i+1]]),noPointsY[i]))
f.close()

f = open("ZLength.txt","w")
f.write("Length [pm];\tHistogram [counts]\n")
for i in range(len(noPointsZ)):
    f.write("%f\t%f\n" %(np.mean([rangesZ[i],rangesZ[i+1]]),noPointsZ[i]))
f.close()

# Calculate angle histograms
histx = np.histogram(xangles, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)
histy = np.histogram(yangles, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)
histz = np.histogram(zangles, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)

# Convert bin edges from radians to degrees
xAngleX = histX[1]/np.pi*180
yAngleX = histY[1]/np.pi*180
zAngleX = histZ[1]/np.pi*180

# Apply cone correction to histogram values
xAreaSlice = AreaSlice/np.mean(AreaSlice)*np.mean(histX[0]);xAngleY = histX[0]/xAreaSlice
yAreaSlice = AreaSlice/np.mean(AreaSlice)*np.mean(histY[0]);yAngleY = histY[0]/yAreaSlice
zAreaSlice = AreaSlice/np.mean(AreaSlice)*np.mean(histZ[0]);zAngleY = histZ[0]/ZAreaSlice

# Save angle histograms to file
f = open("XAngle.txt","w")
f.write("Angle [degrees];\tRelative occurrence\n")
for i in range(len(xAngleY)):
    f.write("%f\t%f\n" %(np.mean([xAngleX[i],xAngleX[i+1]]),xAngleY[i]))
f.close()

f = open("YAngle.txt","w")
f.write("Angle [degrees];\tRelative occurrence\n")
for i in range(len(yAngleY)):
    f.write("%f\t%f\n" %(np.mean([yAngleX[i],yAngleX[i+1]]),yAngleY[i]))
f.close()

f = open("ZAngle.txt","w")
f.write("Angle [degrees];\tRelative occurrence\n")
for i in range(len(zAngleY)):
    f.write("%f\t%f\n" %(np.mean([zAngleX[i],zAngleX[i+1]]),zAngleY[i]))
f.close()

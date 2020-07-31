''' 
Reads a .CHARMTRJ file and can extract the centre of mass position, centre of charge position,
charge arm vector or length.

The first line of this file begins with "$$" and contains the data needed to assign arrays.
The second line contains information on the format.
Subsequent lines containing a "#" denote timestep.
The remaining lines are the centre of mass, centre of charge, charge arm vector and charge 
arm length.

Output is arrays with centre of mass, centre of charge, charge arm vectors and length. 
These can be turned off by changing the input flags to negative.
'''

import numpy as np
from tqdm import tqdm

def readCharm(x, COM = True, COC = True, ChArm = True, Len = True):
	# Pulls out no of timesteps and number of ions, and also checks trajectory file is valid
	try:
		trajectoryFile = str(x)
		for line in open(trajectoryFile): 
			if line[0]=="$":
				splitLine = line.split(" ")
				tSteps = int(splitLine[2].strip())
				noIons = int(splitLine[4].strip())
				boxSize = float(splitLine[6].strip())
				extractedParams = [tSteps,noIons,boxSize]
				break
	except:
		print("Error: trajectory file is invalid")
	print("\nReading trajectory")
	# Assign arrays for storing data
	com=np.zeros([tSteps,noIons,3])
	coc=np.zeros([tSteps,noIons,3])
	chArm=np.zeros([tSteps,noIons,3])
	lenChArm=np.zeros([tSteps,noIons])
	# Load progress bar
	pbar = tqdm(total = tSteps)
	# Read though file
	for line in open(trajectoryFile):
		# Ignore first comment on file structure and pull out information on timestep
		if line[0]=="#":
			try:
				ionNo = 0
				lineSplit = line.split(" ")
				tStep = int(lineSplit[2].strip())
				pbar.update(1)
			except:
				pass
		elif line[0]!="$":
			lineSplit = line.split(';') # Split line into components
			# Assign data if flagged to do so
			if COM:
				com[tStep,ionNo,0]=float(lineSplit[0].strip()) #x
				com[tStep,ionNo,1]=float(lineSplit[1].strip()) #y
				com[tStep,ionNo,2]=float(lineSplit[2].strip()) #z
			if COC:
				coc[tStep,ionNo,0]=float(lineSplit[3].strip()) #x
				coc[tStep,ionNo,1]=float(lineSplit[4].strip()) #y
				coc[tStep,ionNo,2]=float(lineSplit[5].strip()) #z
			if ChArm:
				chArm[tStep,ionNo,0]=float(lineSplit[6].strip()) #x
				chArm[tStep,ionNo,1]=float(lineSplit[7].strip()) #y
				chArm[tStep,ionNo,2]=float(lineSplit[8].strip()) #z
			if Len:
				lenChArm=float(lineSplit[9].strip()) #length
			ionNo+=1
	pbar.close()
	return com, coc, chArm, lenChArm, extractedParams
''' 
Reads a LAMMPS trajectory and outputs 2 custom files which contain the centre of mass, centre 
of charge, charge arm vector and charge arm length for analysis.

This requires the generation of lammps trajectory files seperately for each atom. 
This needs to include the element, molecule number, charge and position of each atom.

Example:
dump TRAJECTORY all custom 500 trajectory.lmp element mol q xu yu zu
dump_modify TRAJECTORY sort id element N C C C H H H C C H C C F S N O X X X X X X X X X X X

This will then output two files, one for the cation, one for the anion, with the centre of mass
position, centre of charge position, charge arm vector and charge arm length. Output structure:
"COM(x);COM(y);COM(z);\t COC(x);COC(y);COC(z);\t ChArm(x);ChArm(y);ChArm(z);\t ChArm(length)\n"

This can then be read back into COM, COC, ChArm vectors and length using the script "ReadChArm.py" 
'''
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def strToEmp(x):
	x=[char for char in x] # Turn vector into list
	Elements = []
	# Find elements in string 
	for item in x:
		if item not in Elements:
			Elements += [item]
	Emp = ""
	if "C" in Elements:
		Emp+="C"
		count = sum([char=="C" for char in x])
		if count>1:
			Emp+=str(count)
		Elements.remove("C")
	if "H" in Elements:
		Emp+="H"
		count = sum([char=="H" for char in x])
		if count>1:
			Emp+=str(count)
		Elements.remove("H")
	Elements = sorted(Elements, key=str.upper)
	for element in Elements:
		Emp+=element
		count = sum([char==element for char in x])
		if count>1:
			Emp+=str(count)
	return Emp
	
atomMassDictionary = {"C":11.611,"H":1.008,"F":18.998,"N":13.607,"O":15.599,"S":31.666,"B":10.411,"X":0.4} # Dictionary of masses of atoms

noParticlesNext = False;boxSizeNext = False
totalNoIons = 0
formulaList = ["",""] ; tempFormula = ""
ionChargeList = [0,0] 
atomMassList = [0,0]
tSteps = 0 # Creates initial value which is increased when counting through file

trajectoryFile = "trajectory.lmp"

print("Performing trajectory recognition...")

for line in open(trajectoryFile,'r'): 
	
	splitLine = line.split()
	
	# Pulls box size information out of trajectory
	if boxSizeNext == True:
		boxSize = float(splitLine[1])-float(splitLine[0])
		boxSizeNext = False
		
	# Activates flag for box size measuring
	if "BOUNDS" in splitLine:
		if tSteps == 1:
			boxSizeNext = True
	
	# Pulls out number of particles in simulation in first timestep only
	if noParticlesNext == True: 
		noParticles = int(splitLine[0])
		noParticlesNext = False
	
	# Counts number of timestep and activates flas for number of particle recognition
	if "NUMBER" in splitLine:
		if tSteps == 0:
			noParticlesNext = True
		tSteps+=1
	
	if tSteps==1: # Recognise one ion in first timestep
		
		# Try statement needed as some line splits do not have a second value, leads to errors if removed
		try:
			# Value will increase to maximum, which is total number of ions
			molNo = int(splitLine[1])
			if totalNoIons<molNo:
				totalNoIons = molNo
			# Saves information about molecule 1 (guaranteed to be unique) 
			if molNo==1:
				atomName = splitLine[0]
				atomCharge = float(splitLine[2])
				atomMass = atomMassDictionary[atomName]   
				
				formulaList[0] += atomName
				ionChargeList[0] += atomCharge
				atomMassList[0] += atomMass
				
		except:
			pass
	
	if tSteps==2: # Recognise second ion in second timestep
		
		try:
			molNo = int(splitLine[1])
			# Saves information about last molecule (assumed to be the other ion) 
			if molNo==totalNoIons:
				atomName = splitLine[0]
				atomCharge = float(splitLine[2])
				atomMass = atomMassDictionary[atomName]   
				
				formulaList[1] += atomName
				ionChargeList[1] += atomCharge
				atomMassList[1] += atomMass 
				
		except:
			pass

# Accounts for floating point errors
ionChargeList = [round(num, 5) for num in ionChargeList]
atomMassList  = [round(num, 5) for num in atomMassList]
# Assuming 2 ions, calculates number of ions
noIons = int(totalNoIons/2)

# Determines which ion is cation and ion by higher charge
if ionChargeList[0]>ionChargeList[1]:
	cationIdx = 0
	anionIdx = 1
else:
	cationIdx = 1
	anionIdx = 0

# Assigns ion masses, charges and formula to cation and anion
cationMass = atomMassList[cationIdx]
cationCharge = ionChargeList[cationIdx]
cationFormula = strToEmp(formulaList[cationIdx])

anionMass = atomMassList[anionIdx]
anionCharge = ionChargeList[anionIdx]
anionFormula = strToEmp(formulaList[anionIdx])

# Creates list of atoms to pull out of trajectory
atomList = [char for char in formulaList[0]] 
atomList += [char for char in formulaList[1]] 
atomList = list(dict.fromkeys(atomList))

# Prints output of recognition for error checking
print("\nNumber of timesteps = %i"%(tSteps))

print("\nIon recognition output:\n------------------")
print("Assuming 2 ions...")
print("\nCATION:")
print("Formula = %s"%(cationFormula))
print("Mass = %f"%(cationMass))
print("Charge = %f"%(cationCharge))
print("\nANION:")
print("Formula = %s"%(anionFormula))
print("Mass = %f"%(anionMass))
print("Charge = %f"%(anionCharge))

# Initiate values and flags for file scanning
noLinesPerTStep = 10 # Initates value, this is changed when scanned in the file
# Set flags for recognition in file
firstLine = True;recordNextLine = False
lineCounter = 0
tStep = 0

cationXCOM = np.zeros([tSteps,noIons])
cationYCOM = np.zeros([tSteps,noIons])
cationZCOM = np.zeros([tSteps,noIons])
cationXCOC = np.zeros([tSteps,noIons])
cationYCOC = np.zeros([tSteps,noIons])
cationZCOC = np.zeros([tSteps,noIons])

anionXCOM = np.zeros([tSteps,noIons])
anionYCOM = np.zeros([tSteps,noIons])
anionZCOM = np.zeros([tSteps,noIons])
anionXCOC = np.zeros([tSteps,noIons])
anionYCOC = np.zeros([tSteps,noIons])
anionZCOC = np.zeros([tSteps,noIons])

print("\nReading atom positions")

pbar = tqdm(total = tSteps)

for line in open(trajectoryFile): 
	#    print "Line No: %i" %(lineCounter)

	splitLine = line.split()

	if firstLine and recordNextLine:
		noParticles = int(splitLine[0])
		noLinesPerTStep = int(splitLine[0])+9
		firstLine = False

	if "NUMBER" in splitLine:
		recordNextLine = True
    
	if splitLine[0] in atomList: #Only perform these calculations on atoms in ions
		# Get atom type and mass for centre of mass calculation
		atomType = splitLine[0]
		atomMass = atomMassDictionary[atomType]        
		# molNo is the number molecule i.e. molecule 1, molecule 2 etc..
		molNo = int(splitLine[1])-1
		# Determine if cation or anion
		if molNo < noIons:
			Cation = True
		else:
			Cation = False
			molNo = molNo-noIons # Anion molNo start high so need to be reduced so first anion has molNo = 0
		# Pull charge of molecule fom trajectory
		atomCharge = float(splitLine[2])
		# Pull position. Multiply by 100 to convert position from angstroms to pm
		x = float(splitLine[3].strip())
		y = float(splitLine[4].strip())
		z = float(splitLine[5].strip())
		
		# Use cation flag to add mass to cation or anion cumulative centres of mass and charge
		if Cation == True:
			# centre of mass position = sum (x_i*mass_i)/mass_ion
			cationXCOM[tStep,molNo] = cationXCOM[tStep,molNo] + x*atomMass
			cationYCOM[tStep,molNo] = cationYCOM[tStep,molNo] + y*atomMass
			cationZCOM[tStep,molNo] = cationZCOM[tStep,molNo] + z*atomMass
			
			cationXCOC[tStep,molNo] = cationXCOC[tStep,molNo] + x*atomCharge
			cationYCOC[tStep,molNo] = cationYCOC[tStep,molNo] + y*atomCharge
			cationZCOC[tStep,molNo] = cationZCOC[tStep,molNo] + z*atomCharge
		elif Cation == False:
			anionXCOM[tStep,molNo] = anionXCOM[tStep,molNo] + x*atomMass
			anionYCOM[tStep,molNo] = anionYCOM[tStep,molNo] + y*atomMass
			anionZCOM[tStep,molNo] = anionZCOM[tStep,molNo] + z*atomMass
			
			anionXCOC[tStep,molNo] = anionXCOC[tStep,molNo] + x*atomCharge
			anionYCOC[tStep,molNo] = anionYCOC[tStep,molNo] + y*atomCharge
			anionZCOC[tStep,molNo] = anionZCOC[tStep,molNo] + z*atomCharge
		
#		if molNo==0:
#			print("Atom label = %s | Atom mass = %f | Molecule number = %f | Cation = %s | Atom Charge = %f | x = %f | y = %f | z = %f " %(atomType,atomMass,molNo,Cation,atomCharge,x,y,z))
		
	if lineCounter%noLinesPerTStep == noLinesPerTStep-1:
		tStep+=1

		pbar.update(1)
	lineCounter+=1
    
pbar.close()

print("Normalising vectors")
# Centres of masses need to be normalised by ion mass
cationXCOM = cationXCOM/cationMass
cationYCOM = cationYCOM/cationMass
cationZCOM = cationZCOM/cationMass

anionXCOM = anionXCOM/anionMass
anionYCOM = anionYCOM/anionMass
anionZCOM = anionZCOM/anionMass

# Centres of charge need to be normalised by ion charge
cationXCOC = cationXCOC/cationCharge
cationYCOC = cationYCOC/cationCharge
cationZCOC = cationZCOC/cationCharge

anionXCOC = anionXCOC/anionCharge
anionYCOC = anionYCOC/anionCharge
anionZCOC = anionZCOC/anionCharge

print("Calculating charge arms")
# Calculate charge arm vector for each ion
cationXChArm = cationXCOC - cationXCOM
cationYChArm = cationYCOC - cationYCOM
cationZChArm = cationZCOC - cationZCOM

anionXChArm = anionXCOC - anionXCOM
anionYChArm = anionYCOC - anionYCOM
anionZChArm = anionZCOC - anionZCOM

cationChArmLength = np.zeros([tSteps,noIons])
anionChArmLength = np.zeros([tSteps,noIons])

for steps in range(tSteps):
	for molecules in range(noIons):
		cationChArmLength[steps,molecules] = np.sqrt(cationXChArm[steps,molecules]**2+cationYChArm[steps,molecules]**2+cationZChArm[steps,molecules]**2)
		anionChArmLength[steps,molecules] = np.sqrt(anionXChArm[steps,molecules]**2+anionYChArm[steps,molecules]**2+anionZChArm[steps,molecules]**2)

# Save charge arm data to one .ChArmTrj file for the cation and one .ChArmTrj file for the anion
print("Writing centre of mass/centre of charge trajectories")
catOutFile = "[" + cationFormula + "]in[" + cationFormula + "][" + anionFormula + "].CHARMTRJ"
anOutFile =  "[" + anionFormula + "]in[" + cationFormula + "][" + anionFormula + "].CHARMTRJ"

catOutWrite = open(catOutFile,"w")
anOutWrite = open(anOutFile,"w")

catOutWrite.write("$$ NoTimeSteps= %i NoIons= %i BoxSize= %f\n" %(tSteps,noIons,boxSize))
anOutWrite.write("$$ NoTimeSteps= %i NoIons= %i BoxSize= %f\n" %(tSteps,noIons,boxSize))

catOutWrite.write("#COM coordinates (x,y,z); COC coordinates (x,y,z); Charge arm vector (x,y,z); Charge Arm Length\n")
anOutWrite.write("#COM coordinates (x,y,z); COC coordinates (x,y,z); Charge arm vector (x,y,z); Charge Arm Length\n")

for step in range(tSteps):
	catOutWrite.write("#Timestep = %i\n"%(step))
	anOutWrite.write("#Timestep = %i\n"%(step))
	for molecule in range(noIons):
        
		catOutWrite.write("%f;%f;%f;\t"%(cationXCOM[step,molecule],cationYCOM[step,molecule],cationZCOM[step,molecule]))
		catOutWrite.write("%f;%f;%f;\t"%(cationXCOC[step,molecule],cationYCOC[step,molecule],cationZCOC[step,molecule]))
		catOutWrite.write("%f;%f;%f;\t"%(cationXChArm[step,molecule],cationYChArm[step,molecule],cationZChArm[step,molecule]))
		catOutWrite.write("%f\n"%(cationChArmLength[step,molecule]))
		
		anOutWrite.write("%f;%f;%f;\t"%(anionXCOM[step,molecule],anionYCOM[step,molecule],anionZCOM[step,molecule]))
		anOutWrite.write("%f;%f;%f;\t"%(anionXCOC[step,molecule],anionYCOC[step,molecule],anionZCOC[step,molecule]))
		anOutWrite.write("%f;%f;%f;\t"%(anionXChArm[step,molecule],anionYChArm[step,molecule],anionZChArm[step,molecule]))
		anOutWrite.write("%f\n"%(anionChArmLength[step,molecule]))
		
catOutWrite.close()

anOutWrite.close()


print("Average length of cation charge arm = %f angstrom"%(np.mean(cationChArmLength)))
print("Average length of anion charge arm = %f angstrom"%(np.mean(anionChArmLength)))

# cationChArmLength = np.reshape(cationChArmLength,np.count_nonzero(cationChArmLength))
# anionChArmLength = np.reshape(anionChArmLength,np.count_nonzero(anionChArmLength))
# plt.hist(cationChArmLength,histtype='step')
# plt.hist(anionChArmLength,histtype='step')
# plt.show()



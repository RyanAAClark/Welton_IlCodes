''' 
Part of the processing pipeline to get the autocorrelation of different ions in ionic liquid along 
different dimensions

This requires the generation of lammps trajectory files seperately for each atom. 
This needs to be sorted by id, and include the element, mol, and velocities of each atom.

Example:
#dump TRAJECTORY CATION custom 1 CationTrajectory.lmp element mol vx vy vz
#dump_modify TRAJECTORY sort id element N C C C H H H C C H C C F S N O X X X X X X X X X X X
NB: This example given is for an imidazolium cation defined by the group CATION

This part of the pipeline takes the seperate files and turns them into one centre of mass trajectory 
for each ion resolved along each vector of the box.

The next step of the pipeline uses "AutocorrelateCationTraj.py" and "AutocorrelateAnionTraj.py"
to analyse the centre of mass trajectories generated here
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

atomMassDictionary = {"C":11.611,"H":1.008,"F":18.998,"N":13.607,"O":15.599,"S":31.666,"X":0.4}
cationMass =  139.222
anionMass = 280.145
noIons = 256
tSteps = 1000001 # Change to the total number of timesteps in the file

catFile = "CationTrajectory.lmp"
anFile = "AnionTrajectory.lmp"

noLinesPerTStep = 10 #Creates initial value, this is changed when scanned in the file
firstLine = True
recordNextLine = False
lineCounter = 0
tStep = 0

cationXVel = np.zeros([tSteps,noIons])
cationYVel = np.zeros([tSteps,noIons])
cationZVel = np.zeros([tSteps,noIons])

anionXVel = np.zeros([tSteps,noIons])
anionYVel = np.zeros([tSteps,noIons])
anionZVel = np.zeros([tSteps,noIons])

print("Reading cation trajectory")

for line in open(catFile): 
#    print "Line No: %i" %(lineCounter)
    
    splitLine = line.split()

    if firstLine and recordNextLine:
        noParticles = int(splitLine[0])
        noLinesPerTStep = int(splitLine[0])+9
        firstLine = False
        
    if "NUMBER" in splitLine:
        recordNextLine = True
    
    if splitLine[0] in ["C","N","H","X"]:
        # Get atom type and mass for centre of mass calculation
        atomType = splitLine[0]
        atomMass = atomMassDictionary[atomType]        
        # molNo is the number molecule i.e. cation 1, cation 2 etc..
        molNo = int(splitLine[1])-1
        # Convert velocity from angstroms per femtosecond to metres/second (or pm/ps)
        vx = float(splitLine[2].strip())*100000.
        vy = float(splitLine[3].strip())*100000.
        vz = float(splitLine[4].strip())*100000.
        
        # centre of mass velocity = sum (vel_i*mass_i)/mass_ion
        cationXVel[tStep,molNo] = cationXVel[tStep,molNo] + vx*atomMass
        cationYVel[tStep,molNo] = cationYVel[tStep,molNo] + vy*atomMass
        cationZVel[tStep,molNo] = cationZVel[tStep,molNo] + vz*atomMass
#        print "TimeStep %i | Mol No %i" %(tStep,molNo)
    
    if lineCounter%noLinesPerTStep == noLinesPerTStep-1:
        tStep+=1
        
        sys.stdout.write("\r%i of %i" %(tStep,tSteps))
        sys.stdout.flush()

    lineCounter+=1
    
sys.stdout.write("\n")
sys.stdout.flush()

print("Weighing cation velocities by mass")
cationXVel = cationXVel/cationMass
cationYVel = cationYVel/cationMass
cationZVel = cationZVel/cationMass

# Calculate average cation velocity
cationXAvgVel = np.mean(cationXVel)
cationYAvgVel = np.mean(cationYVel)
cationZAvgVel = np.mean(cationZVel)

firstLine = True    
recordNextLine = False
tStep = 0
noLinesPerTStep = 10
lineCounter = 0

print("Reading anion trajectory")
    
for line in open(anFile):
#    print "Line No: %i" %(lineCounter)
    
    splitLine = line.split()

    if firstLine and recordNextLine:
        noParticles = int(splitLine[0])
        noLinesPerTStep = int(splitLine[0])+9
        firstLine = False
        
    if "NUMBER" in splitLine:
        recordNextLine = True
    
    if splitLine[0] in ["C","F","N","O","S","X"]:
        # Get atom type and mass for centre of mass calculation
        atomType = splitLine[0]
        atomMass = atomMassDictionary[atomType]        
        # molNo is the number molecule i.e. anon 1, anion 2 etc..
        molNo = int(splitLine[1])-257
        # Convert velocity from angstroms per femtosecond to metres/second (or pm/ps)
        vx = float(splitLine[2].strip())*100000.
        vy = float(splitLine[3].strip())*100000.
        vz = float(splitLine[4].strip())*100000.
        
        # centre of mass velocity = sum (vel_i*mass_i)/mass_ion
        anionXVel[tStep,molNo] = anionXVel[tStep,molNo] + vx*atomMass
        anionYVel[tStep,molNo] = anionYVel[tStep,molNo] + vy*atomMass
        anionZVel[tStep,molNo] = anionZVel[tStep,molNo] + vz*atomMass
#        print "TimeStep %i | Mol No %i" %(tStep,molNo)
    
    if lineCounter%noLinesPerTStep == noLinesPerTStep-1:
        tStep+=1
        
        sys.stdout.write("\r%i of %i" %(tStep,tSteps))
        sys.stdout.flush()
    
    lineCounter+=1

sys.stdout.write("\n")
sys.stdout.flush()

print("Weighing anion velocities by mass")
anionXVel = anionXVel/anionMass
anionYVel = anionYVel/anionMass
anionZVel = anionZVel/anionMass

# Calculate average anion velocity
anionXAvgVel = np.mean(anionXVel)
anionYAvgVel = np.mean(anionYVel)
anionZAvgVel = np.mean(anionZVel)

#plt.plot(cationZVel[:,0]);plt.plot(cationZVel[:,1]);plt.plot(cationZVel[:,2]);plt.plot(cationZVel[:,3]);plt.plot(cationZVel[:,4]);plt.plot(cationZVel[:,5]);plt.show()
#plt.plot(anionZVel[:,0]);plt.plot(anionZVel[:,1]);plt.plot(anionZVel[:,2]);plt.plot(anionZVel[:,3]);plt.plot(anionZVel[:,4]);plt.plot(anionZVel[:,5]);plt.show()

avgVelOutput = open("AverageVelocitySummary.txt","w")
avgVelOutput.write("Average cation velocity in x = %f m/s\n" %(cationXAvgVel))
avgVelOutput.write("Average cation velocity in y = %f m/s\n" %(cationYAvgVel))
avgVelOutput.write("Average cation velocity in z = %f m/s\n\n" %(cationZAvgVel))
avgVelOutput.write("Average anion velocity in x = %f m/s\n" %(anionXAvgVel))
avgVelOutput.write("Average anion velocity in y = %f m/s\n" %(anionYAvgVel))
avgVelOutput.write("Average anion velocity in z = %f m/s\n" %(anionZAvgVel))
avgVelOutput.close()

# Save files
print("Writing centre of mass trajectories")
anXOut = "AnionXTrajectory.txt";anXWrite = open(anXOut,"w")
anYOut = "AnionYTrajectory.txt";anYWrite = open(anYOut,"w")
anZOut = "AnionZTrajectory.txt";anZWrite = open(anZOut,"w")

catXOut = "CationXTrajectory.txt";catXWrite = open(catXOut,"w")
catYOut = "CationYTrajectory.txt";catYWrite = open(catYOut,"w")
catZOut = "CationZTrajectory.txt";catZWrite = open(catZOut,"w")
for step in range(tSteps):
    for molecule in range(noIons):
        if molecule != noIons-1:
            anXWrite.write("%f;"%(anionXVel[step,molecule]))
            anYWrite.write("%f;"%(anionYVel[step,molecule]))
            anZWrite.write("%f;"%(anionZVel[step,molecule]))
            
            catXWrite.write("%f;"%(cationXVel[step,molecule]))
            catYWrite.write("%f;"%(cationYVel[step,molecule]))
            catZWrite.write("%f;"%(cationZVel[step,molecule]))
        else:
            anXWrite.write("%f"%(anionXVel[step,molecule]))
            anYWrite.write("%f"%(anionYVel[step,molecule]))
            anZWrite.write("%f"%(anionZVel[step,molecule]))
            
            catXWrite.write("%f"%(cationXVel[step,molecule]))
            catYWrite.write("%f"%(cationYVel[step,molecule]))
            catZWrite.write("%f"%(cationZVel[step,molecule]))

    anXWrite.write("\n")
    anYWrite.write("\n")
    anZWrite.write("\n")
    
    catXWrite.write("\n")
    catYWrite.write("\n")
    catZWrite.write("\n")
    
anXWrite.close()
anYWrite.close()
anZWrite.close()

catXWrite.close()
catYWrite.close()
catZWrite.close()




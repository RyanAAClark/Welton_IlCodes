import numpy as np
import collections

File = "trajectory.lmp"

ReadNextLine = False
firstRepeated = True
writeToFile = True
noTimeSteps = 0
tstepcounter = 0

print("finding no of timesteps")

for line in open(File):
    if "BOX" in line:
        noTimeSteps+=1
        
print(noTimeSteps)

timestepVals = np.zeros(noTimeSteps)

print("finding repeated timesteps")

for line in open(File):
    if ReadNextLine:
        temporary = float(line)
        timestepVals[tstepcounter] = temporary
        tstepcounter += 1
        ReadNextLine = False
        
        
    if "TIMESTEP" in line:
        ReadNextLine = True

print("repeated timesteps are" )
repeatedTSteps =  [item for item, count in collections.Counter(timestepVals).items() if count > 1]

print(repeatedTSteps)

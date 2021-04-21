'''
Python file to count the number of repeated timesteps in a LAMMPS trajectory file,
specified by "File" then create a new trajectory named "newFileName" with the 
repeated steps removed.



    E.g. will turn a trajectory with timesteps of:
        1a, 2a, 3a, 4a, 5a, 3b, 4b, 5b, 6b, 7b, 5c, 6c, 7c, 8c....
		       |___________|  &    |__________| <-Removed
    into:
        1a, 2a, 3b, 4b, 5c, 6c, 7c, 8c....
'''
import numpy as np
from collections import Counter

File = "new_dump.lammpstrj" #File with repeated timesteps
newFileName = "new_new_dump.lammpstrj" #Name for the trajectory to write without any repeated timesteps

timestepVals = []

print("finding repeated timesteps")
notimesteps=0
ReadNextLine=False
# Creates a list of all timesteps in the file
for line in open(File):
    if ReadNextLine:
        temporary = str(int(line))
        notimesteps+=1
        timestepVals.append(temporary)
        ReadNextLine = False
        
    if "TIMESTEP" in line:
        ReadNextLine = True
        
## Finds the repeated values in the list of timesteps
# Creates a counter class (dictionary subclass) with each timestep and the number of times it appears
repeatedTSteps = Counter(timestepVals) 
# Creates another counter class with each timestep appearing only once
tSteps = Counter([i for i in repeatedTSteps])  
# Take one from the count of each timestep
repeatedTSteps.subtract(tSteps)
# Removes empty values (these are timesteps which are not repeated) 
repeatedTSteps += Counter() 
# Converts to list of timesteps, removing all non-repeated timesteps
repeatedTSteps = list(repeatedTSteps.elements())

# Checks the file will no be overwritten
if File==newFileName:
    raise ValueError('Old and new trajectory cannot have the same name!')

# Initiate Variables
ReadNextLine = False
writeToFile = True

newFile = open(newFileName,"w")
print("writing new file")
# Write the trajectory to the new file name, without any of the repeated timesteps 
for line in open(File, "r"):
	
	if ReadNextLine:
		
		currenttStep = str(int((line)))
		
		if currenttStep in repeatedTSteps:
			writeToFile=False
			repeatedTSteps.remove(currenttStep)
		
		else:
			writeToFile=True
		
		ReadNextLine = False
    
	if writeToFile:
		newFile.write(line)
                
	if "TIMESTEP" in line:
		ReadNextLine = True

newFile.close()

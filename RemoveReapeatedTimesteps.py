'''
Python file to count the number of repeated timesteps in a LAMMPS trajectory file, specified by "File" then create a 
new trajectory named "newFileName" with the repeated steps removed.

Can handle up to three repeated timesteps in a row.
    E.g. will turn a trajectory with timesteps of:
        1, 2, 3, 4, 5, 3, 4, 5, 6, 7, 8....
             |________|<-Removed
    into:
        1, 2, 3, 4, 5, 6, 7, 8....

Keeps the second set of positions, assumes that timesteps 3 are the same and that the second set of timesteps are correct.

To debug/see what the script is doing, uncomment the print function at line 117.
This will give:
    TStep - Current timestep
    WTF - Write to file, TRUE to include current line, FALSE to exclude current line
    FR - First Repeated, TRUE before the first repeated timestep, FALSE after the first repeated step
    MR - Middle Repeated, TRUE before the middle of the repeated timestep, FALSE after the middle of the repeated steps
    LR - Last Repeated, TRUE before the last repeated timestep, FALSE after the last repeated step
        FR, MR, LR are reset to zero once the file leaves the range of repeated timesteps
'''
import numpy as np
import collections

File = "trajectory.lmp" #File with repeated timesteps
newFileName = "new_trajectory.lmp" #Name for the trajectory to write without any repeated timesteps
timeDiff = 1000.0 #Difference between the timesteps in the units in the trajectory

# Initialise flags for determining whether to write timesteps to the new file
ReadNextLine = False
firstRepeated = True
lastRepeated = True
midRepeated = True
writeToFile = True
noTimeSteps = 0
tstepcounter = 0

print("finding no of timesteps")

# Finds the number of timesteps in the original trajectory and print to terminal
for line in open(File):
    if "BOX" in line:
        noTimeSteps+=1
        
print(noTimeSteps)

timestepVals = np.zeros(noTimeSteps)

print("finding repeated timesteps")

# Creates a list of all timesteps in the file
for line in open(File):
    if ReadNextLine:
        temporary = float(line)
        timestepVals[tstepcounter] = temporary
        tstepcounter += 1
        ReadNextLine = False
        
    if "TIMESTEP" in line:
        ReadNextLine = True

# Finds the repeated values in the list of timesteps
print("repeated timesteps are" )
repeatedTSteps =  [item for item, count in collections.Counter(timestepVals).items() if count > 1]
print(repeatedTSteps)

if File==newFileName:
    raise ValueError('Old and new trajectory cannot have the same name!')

newFile = open(newFileName,"w")
print("writing new file")
# Write the trajectory to the new file name, without any of the repeated timesteps 
for line in open(File):
            
    if ReadNextLine:
        currenttStep = float(line)

        if currenttStep in repeatedTSteps and currenttStep+timeDiff in repeatedTSteps and currenttStep-timeDiff in repeatedTSteps:
            if midRepeated:
                writeToFile = False
                midRepeated = False
            
            else:
                print("Skipped timestep:%i " %(currenttStep))
                writeToFile = True
                midRepeated = True

        elif currenttStep in repeatedTSteps and currenttStep+timeDiff in repeatedTSteps:
            if firstRepeated:
                firstRepeated=False
                writeToFile = False
                
            else:
                print("Skipped timestep:%i " %(currenttStep))
                writeToFile = True
     
        elif currenttStep in repeatedTSteps and currenttStep-timeDiff in repeatedTSteps:
            
            if lastRepeated:
                lastRepeated = False

            else:
                lastRepeated = True
                firstRepeated = True
                print("Skipped timestep:%i " %(currenttStep))
        
        elif currenttStep in repeatedTSteps:
            if firstRepeated:
                firstRepeated = False
                writeToFile = False
            else:
                print("Skipped timestep:%i " %(currenttStep))
                firstRepeated = True
                writeToFile = True
        ReadNextLine = False
#        print("TStep: %i | WTF = %s | FR = %s | MR = %s | LR = %s " %(currenttStep,writeToFile,firstRepeated,midRepeated,lastRepeated))
    
    if writeToFile:
        newFile.write(line)
                
    if "TIMESTEP" in line:
        ReadNextLine = True

newFile.close()

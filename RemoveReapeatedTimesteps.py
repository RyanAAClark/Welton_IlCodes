import numpy as np
import collections

File = "trajectory.lmp"
newFileName = "new_trajectory.lmp"

ReadNextLine = False
firstRepeated = True
lastRepeated = True
midRepeated = True
writeToFile = True
noTimeSteps = 0
tstepcounter = 0
#lineCounter = 0
timeDiff=1000.0

print("finding no of timesteps")

for line in open(File):
    if "BOX" in line:
        noTimeSteps+=1
        
print(noTimeSteps)
#noTimeSteps = 20016 

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

newFile = open(newFileName,"w")
print("writing new file")
for line in open(File):
            
    if ReadNextLine:
        currenttStep = float(line)
        #print("TStep: %i | WTF = %s | FR = %s" %(currenttStep,writeToFile,firstRepeated))

        if currenttStep in repeatedTSteps and currenttStep+timeDiff in repeatedTSteps and currenttStep-timeDiff in repeatedTSteps:
            #print("4")
            if midRepeated:
                print("4.1")
                writeToFile = False
                midRepeated = False
            
            else:
                print("4.2")
                print("Skipped timestep:%i " %(currenttStep))
                writeToFile = True
                midRepeated = True

        elif currenttStep in repeatedTSteps and currenttStep+timeDiff in repeatedTSteps:
           # print("1")
            if firstRepeated:
                print("1.1")
                firstRepeated=False
                writeToFile = False
                
            else:
                print("1.2")
                print("Skipped timestep:%i " %(currenttStep))
                writeToFile = True
     
        elif currenttStep in repeatedTSteps and currenttStep-timeDiff in repeatedTSteps:
            #print("2")
            
            if lastRepeated:
                lastRepeated = False
                print("2.1")
#                print("---TStep: %i | WTF = %s | FR = %s" %(currenttStep,writeToFile,firstRepeated))

            else:
                print("2.2")
                lastRepeated = True
                firstRepeated = True
                print("Skipped timestep:%i " %(currenttStep))
#            print("---TStep: %i | WTF = %s | FR = %s" %(currenttStep,writeToFile,firstRepeated))
        
        elif currenttStep in repeatedTSteps:
            #print("3")
            if firstRepeated:
                print("3.1")
                firstRepeated = False
                writeToFile = False
            else:
                print("3.2")
                print("Skipped timestep:%i " %(currenttStep))
                firstRepeated = True
                writeToFile = True
        ReadNextLine = False
        print("TStep: %i | WTF = %s | FR = %s | MR = %s | LR = %s " %(currenttStep,writeToFile,firstRepeated,midRepeated,lastRepeated))
    
    if writeToFile:
        newFile.write(line)
        
        
        
        
    if "TIMESTEP" in line:
        ReadNextLine = True


newFile.close()

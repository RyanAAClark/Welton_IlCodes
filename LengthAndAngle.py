import numpy as np
import matplotlib.pyplot as plt
import sys

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
     BigArea[angles] = 2*np.pi*(1-np.pi*np.cos((angles+1)*np.pi/180))
     SmallArea[angles] = 2*np.pi*(1-np.pi*np.cos((angles)*np.pi/180))
 
AreaSlice = BigArea-SmallArea

#FieldStrengths = ["0.000","0.001","0.005"]
FieldStrengths = ["0.000","0.001","0.002","0.005","0.010","0.020","0.050","0.100","0.200","0.300","0.400","0.500"]

for FS in FieldStrengths:
    
    FirstN=True

    File = "Trajectories/" + FS + "_trajectory_out.xyz"
    
    ## Pre-allocate all the vectors to store results
    ## Length
    xLength = np.zeros(20001*256)
    yLength = np.zeros(20001*256)
    zLength = np.zeros(20001*256)
    Length = np.zeros(20001*256)
    ## Angle
    xAngle = np.zeros(20001*256)
    yAngle = np.zeros(20001*256)
    zAngle = np.zeros(20001*256)
    
    ## Counter for inserting into vecotrs
    i=0

    for line in open(File):
        
        if line [:1]=="N" and line [:2]!="No" and FirstN:
            xyz = line.split()
            n1Pos = np.array([float(xyz[1]),float(xyz[2]),float(xyz[3])])
            FirstN=False

        elif line [:1]=="N" and line [:2]!="No":
            FirstN=True
            xyz = line.split()
            n2Pos = np.array([float(xyz[1]),float(xyz[2]),float(xyz[3])])
            
            ## Vector from N1 to N2
            vector = n2Pos-n1Pos
            
            ## Put lengths into vectors
            xLength[i] = vector[0]
            yLength[i] = vector[1]
            zLength[i] = vector[2]
            Length[i] = np.linalg.norm(vector)
            
            ## Calculate angles and put into vectors
            xAngle[i] = angle(vector,[1,0,0])
            yAngle[i] = angle(vector,[0,1,0])
            zAngle[i] = angle(vector,[0,0,1])
            
            ## Increase counter position by 1
            i+=1
            
            # Print timestep when it increases
            if int(i/256)*256 == i and int(i/256)==1:
                sys.stdout.write("%i of 20001" %(int(i/256)))
            elif int(i/256)*256 == i and int(i/256)==20001:
                sys.stdout.flush()
                sys.stdout.write("\r%i of 20001\n" %(int(i/256)))
            elif int(i/256)*256 == i:
                sys.stdout.flush()
                sys.stdout.write("\r%i of 20001" %(int(i/256)))
    
    ## Calculate length histograms        
    xLenHist = np.histogram(xLength, bins=200, range=(-7.5,7.5), normed=None, weights=None, density=False)   
    yLenHist = np.histogram(yLength, bins=200, range=(-7.5,7.5), normed=None, weights=None, density=False)   
    zLenHist = np.histogram(zLength, bins=200, range=(-7.5,7.5), normed=None, weights=None, density=False)   
    LenHist = np.histogram(Length, bins=200, range=(0,7.5), normed=None, weights=None, density=False)   

    ## Calculate angle histograms
    xAngHist = np.histogram(xAngle, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)    
    yAngHist = np.histogram(yAngle, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)
    zAngHist = np.histogram(zAngle, bins=180, range=(0,np.pi), normed=None, weights=None, density=False)
    
    # X values are in hist[1], y values in hist[0]

    ## Save files with length histograms
    f1 = open(FS + "_XLength.txt","w")
    xLengthX = xLenHist[1]
    xLengthY = xLenHist[0]
    for i in range(len(xLengthY)):
        f1.write("%f\t%f\n" %(xLengthX[i],xLengthY[i]))
    f1.close()
    
    f2 = open(FS + "_YLength.txt","w")
    yLengthX = yLenHist[1]
    yLengthY = yLenHist[0]
    for i in range(len(yLengthY)):
        f2.write("%f\t%f\n" %(yLengthX[i],yLengthY[i]))
    f2.close()
    
    f3 = open(FS + "_ZLength.txt","w")
    zLengthX = zLenHist[1]
    zLengthY = zLenHist[0]
    for i in range(len(zLengthY)):
        f3.write("%f\t%f\n" %(zLengthX[i],zLengthY[i]))
    f3.close()
    
    f4 = open(FS + "_Length.txt","w")
    LengthX = LenHist[1]
    LengthY = LenHist[0]
    for i in range(len(LengthY)):
        f4.write("%f\t%f\n" %(LengthX[i],LengthY[i]))
    f4.close()
    
    ## Save files with angle histograms
    f5 = open(FS + "_XAngle.txt","w")
    xAngleX = xAngHist[1]
    xAngleY = xAngHist[0]
    for i in range(len(xAngleY)):
        f5.write("%f\t%f\n" %(180*np.mean([xAngleX[i],xAngleX[i+1]])/np.pi,xAngleY[i]))
    f5.close()
    
    f6 = open(FS + "_YAngle.txt","w")
    yAngleX = yAngHist[1]
    yAngleY = yAngHist[0]
    for i in range(len(yAngleY)):
        f6.write("%f\t%f\n" %(180*np.mean([yAngleX[i],yAngleX[i+1]])/np.pi,yAngleY[i]))
    f6.close()
    
    f7 = open(FS + "_ZAngle.txt","w")
    zAngleX = zAngHist[1]
    zAngleY = zAngHist[0]
    for i in range(len(zAngleY)):
        f7.write("%f\t%f\n" %(180*np.mean([zAngleX[i],zAngleX[i+1]])/np.pi,zAngleY[i]))
    f7.close()
    
    ## Save files with normalised angle histograms
    f5 = open(FS + "_XAngleNorm.txt","w")
    xAngleX = xAngHist[1]
    xAreaSlice=AreaSlice/np.mean(AreaSlice)*np.mean(xAngHist[0])
    xAngleY = xAngHist[0]/AreaSlice
    for i in range(len(xAngleY)):
        f5.write("%f\t%f\n" %(180*np.mean([xAngleX[i],xAngleX[i+1]])/np.pi,xAngleY[i]))
    f5.close()
    
    f6 = open(FS + "_YAngleNorm.txt","w")
    yAngleX = yAngHist[1]
    yAreaSlice=AreaSlice/np.mean(AreaSlice)*np.mean(yAngHist[0])
    yAngleY = yAngHist[0]/AreaSlice
    for i in range(len(yAngleY)):
        f6.write("%f\t%f\n" %(180*np.mean([yAngleX[i],yAngleX[i+1]])/np.pi,yAngleY[i]))
    f6.close()
    
    f7 = open(FS + "_ZAngleNorm.txt","w")
    zAngleX = zAngHist[1]
    zAreaSlice = AreaSlice/np.mean(AreaSlice)*np.mean(zAngHist[0])
    zAngleY = zAngHist[0]/AreaSlice
    for i in range(len(zAngleY)):
        f7.write("%f\t%f\n" %(180*np.mean([zAngleX[i],zAngleX[i+1]])/np.pi,zAngleY[i]))
    f7.close()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



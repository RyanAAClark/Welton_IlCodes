from ReadChArm import readCharm
from WrapArray import wrapTrajectory

com,coc,chArm,lenChArm,extractedParams = readCharm('[C6H11N2X8]in[C6H11N2X8][BF4X5].CHARMTRJ')
com,coc,chArm,lenChArm,_ = readCharm('[BF4X5]in[C6H11N2X8][BF4X5].CHARMTRJ')
[tSteps,noIons,boxSize] = extractedParams

#print(com)

comtStep = com[0,:,:]
coctStep = coc[0,:,:]

print(comtStep)

comtStep = wrapTrajectory((comtStep,coctStep), boxSize)

print(comtStep)

print('Done')
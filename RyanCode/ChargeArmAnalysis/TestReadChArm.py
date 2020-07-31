from ReadChArm import readCharm

com,coc,chArm,lenChArm,extractedParams = readCharm('[C6H11N2X8]in[C6H11N2X8][C3F6NO4S2X10].CHARMTRJ')
com,coc,chArm,lenChArm,_ = readCharm('[C3F6NO4S2X10]in[C6H11N2X8][C3F6NO4S2X10].CHARMTRJ')

[tSteps,noIons,boxSize] = extractedParams

print('Done')
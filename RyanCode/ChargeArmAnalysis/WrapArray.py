''' 
Takes an input trajectory from a .CHARMTRJ file and wraps the centre of mass and centre of charge 
into the box if they are both outside. 

Input should be a list of [centre of mass, centre of charge] positions, along with the box size 
to wrap the molecules into.
'''

import numpy as np

def wrapTrajectory(traj, boxSize):
	
	# Split trajectory into centre of mass and centre of charge vectors
	[COM,COC] = traj
	
	# Split Com and Coc into seperate x,y,z vectors
	COMx = COM[:,0];COMy = COM[:,1];COMz = COM[:,2]
	COCx = COC[:,0];COCy = COC[:,1];COCz = COC[:,2]	
	
	# If either centre of mass or centre of charge is outside box, move both into box
	for i in range(len(COMx)):
		#print("COMx = %f | COCx = %f"%(COMx[i],COCx[i]))
		while COMx[i]>boxSize or COCx[i]>boxSize:
			COMx[i] = COMx[i]-boxSize
			COCx[i] = COCx[i]-boxSize
			#print("New COMx = %f | New COCx = %f"%(COMx[i],COCx[i]))
	
		#print("COMy = %f | COCy = %f"%(COMy[i],COCy[i]))
		while COMy[i]>boxSize or COCy[i]>boxSize:
			COMy[i] = COMy[i]-boxSize
			COCy[i] = COCy[i]-boxSize
			#print("New COMy = %f | New COCy = %f"%(COMy[i],COCy[i]))
	
	
		#print("COMz = %f | COCz = %f"%(COMz[i],COCz[i]))
		while COMz[i]>boxSize or COCz[i]>boxSize:
			COMz[i] = COMz[i]-boxSize
			COCz[i] = COCz[i]-boxSize
			#print("New COMz = %f | New COCz = %f"%(COMz[i],COCz[i]))

		#print("COMx = %f | COCx = %f"%(COMx[i],COCx[i]))
		while COMx[i]<0 or COCx[i]<0:
			COMx[i] = COMx[i]+boxSize
			COCx[i] = COCx[i]+boxSize
			#print("New COMx = %f | New COCx = %f"%(COMx[i],COCx[i]))
	
		#print("COMy = %f | COCy = %f"%(COMy[i],COCy[i]))
		while COMy[i]<0 or COCy[i]<0:
			COMy[i] = COMy[i]+boxSize
			COCy[i] = COCy[i]+boxSize
			#print("New COMy = %f | New COCy = %f"%(COMy[i],COCy[i]))
	
	
		#print("COMz = %f | COCz = %f"%(COMz[i],COCz[i]))
		while COMz[i]<0 or COCz[i]<0:
			COMz[i] = COMz[i]+boxSize
			COCz[i] = COCz[i]+boxSize
			#print("New COMz = %f | New COCz = %f"%(COMz[i],COCz[i]))
	
	# Repack into output
	newCOM = [COMx,COMy,COMz]
	newCOC = [COCx,COCy,COCz]
	
	newtraj = (COM,COC)
	
	return newtraj

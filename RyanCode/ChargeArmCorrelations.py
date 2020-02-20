'''
Calculates the length and angle distributions between charge arms of ions whos centres of mass are 
within a distance criteria.

The angle distribution is a standard angle distribution, taking a histogram of the angle between the two 
charge arm vectors, then applying a cone correction. This means a flat line at a relative abundance of 
one corresponds to a statistically even distribution.

The length distribution looks at the correlation between the length of one charge arm, to the length of
another in a 2D histogram. If there is a tending of this line to y=x then there is a correlation between
the two vectors.
'''
import numpy as np
import sys
import matplotlib.pyplot as plt

timesteps = 20001
ions = 256
distCriteria_CatAn = 950#Distance criteria to put into comparison vector 
distCriteria_CatCat = 1400#If less than this value ion-ion pair will be added to correlation
distCriteria_AnAn = 1300

cationCOMflag = "#4"
cationCOCflag = "#5"
anionCOMflag = "#3"
anionCOCflag = "#4"

cationTraj = "cation_trajectory_out.xyz" 
anionTraj = "anion_trajectory_out.xyz"

def dist1vec(a):
    x = a[0]**2
    y = a[1]**2
    z = a[2]**2
    dist = np.sqrt(x + y + z)
    return dist

def dist2vec(a,b):
    x = (a[0]-b[0])**2
    y = (a[1]-b[1])**2
    z = (a[2]-b[2])**2
    dist = np.sqrt(x + y + z)
    return dist
    
def angle(a, b):
    magA = np.linalg.norm(a)
    magB = np.linalg.norm(b)
    dotprod = np.dot(a,b)
    cosAng = dotprod/(magA*magB)
    Ang = np.arccos(cosAng)
    return Ang

# Cone correction vector for angles
noAngles=180
SmallArea = np.zeros(noAngles)
BigArea = np.zeros(noAngles)
for angles in range(noAngles):
     BigArea[angles] = 2*np.pi*(1-np.pi*np.cos((angles+1)*np.pi/180))
     SmallArea[angles] = 2*np.pi*(1-np.pi*np.cos((angles)*np.pi/180))
ConeCorrectionVector = BigArea-SmallArea

# Pre-allocate cation atom vectors
cationCOM = np.zeros([timesteps*ions,3])
cationChargeArm = np.zeros([timesteps*ions,3])
cationCounter=0

print("Importing cation trajectory...")

# Pull out cation centres of mass, centres of charge and charge arms
for line in open(cationTraj):
    
    if line [:2]==cationCOMflag:
        xyz = line.split()
        comX = float(xyz[1])*100. #Convert from angstrom to pm
        comY = float(xyz[2])*100.
        comZ = float(xyz[3])*100.
        
        cationCOM[cationCounter] = np.array([comX,comY,comZ])

    if line [:2]==cationCOCflag:
        xyz = line.split()
        cocX = float(xyz[1])*100.
        cocY = float(xyz[2])*100.
        cocZ = float(xyz[3])*100.
        
        cationChargeArm[cationCounter] = np.array([cocX-comX,cocY-comY,cocZ-comZ])
        
        cationCounter+=1

# Pre-allocate anion atom vectors        
anionCOM = np.zeros([timesteps*ions,3])
anionChargeArm = np.zeros([timesteps*ions,3])
anionCounter=0

print("Importing anion trajectory")

# Pull out anion centres of mass, centres of charge and charge arms
for line in open(anionTraj):
    
    if line [:2]==anionCOMflag:
        xyz = line.split()
        comX = float(xyz[1])*100. #Convert from angstrom to pm
        comY = float(xyz[2])*100.
        comZ = float(xyz[3])*100.
        
        anionCOM[anionCounter] = np.array([comX,comY,comZ])

    if line [:2]==anionCOCflag:
        xyz = line.split()
        cocX = float(xyz[1])*100.
        cocY = float(xyz[2])*100.
        cocZ = float(xyz[3])*100.
        
        anionChargeArm[anionCounter] = np.array([cocX-comX,cocY-comY,cocZ-comZ])
        
        anionCounter+=1

# Pre-allocate cation-anion, cation-cation and anion-anion comparisons
catan_len_cat = np.zeros(timesteps*(ions**2))
catan_len_an = np.zeros(timesteps*(ions**2))
catan_angle = np.zeros(timesteps*(ions**2))
catan_counter = 0

catcat_len_cat1 = np.zeros(timesteps*(ions**2))
catcat_len_cat2 = np.zeros(timesteps*(ions**2))
catcat_angle = np.zeros(timesteps*(ions**2))
catcat_counter = 0

anan_len_an1 = np.zeros(timesteps*(ions**2))
anan_len_an2 = np.zeros(timesteps*(ions**2))
anan_angle = np.zeros(timesteps*(ions**2))
anan_counter = 0

# Go through data timestep by timestep and perform analyses
for timestep in range(timesteps):
    sys.stdout.write("\r%i of %i" %(timestep+1,timesteps))
    if timestep == timesteps-1:
        sys.stdout.write("\n")
    sys.stdout.flush()
     
    catCOM = cationCOM[timestep*ions:(timestep+1)*ions]
    catChArm = cationChargeArm[timestep*ions:(timestep+1)*ions]
    
    anCOM = anionCOM[timestep*ions:(timestep+1)*ions]
    anChArm = anionChargeArm[timestep*ions:(timestep+1)*ions]
    
    # For every cation...
    for cat in range(ions):
        #...compare with anion
        for an in range(ions):
            cation_com = catCOM[cat]
            cation_charm = catChArm[cat]
            anion_com = anCOM[an]
            anion_charm = anChArm[an]
            
            dist = dist2vec(cation_com,anion_com)
            
            if dist<distCriteria_CatAn:
                catan_len_cat[catan_counter] = dist1vec(cation_charm)
                catan_len_an[catan_counter] = dist1vec(anion_charm)
                
                catan_angle[catan_counter] = angle(cation_charm,anion_charm)
                
                catan_counter+=1
                
        #...compare with cations
        for cat2 in range(ions):
            #Makes sure cation is not the same cation
            if cat!=cat2:
                cation1_com = catCOM[cat]
                cation1_charm = catChArm[cat]
                cation2_com = catCOM[cat2]
                cation2_charm = catChArm[cat2]
                
                dist = dist2vec(cation1_com,cation2_com)

                if dist<distCriteria_CatCat:
                    catcat_len_cat1[catcat_counter] = dist1vec(cation1_charm)
                    catcat_len_cat2[catcat_counter] = dist1vec(cation2_charm)
                    
                    catcat_angle[catcat_counter] = angle(cation1_charm,cation2_charm)
                    
                    catcat_counter+=1
                
    # For every anion...
    for an1 in range(ions):
        #...compare with anion
        for an2 in range(ions):
            if an1!=an2:
                anion1_com = anCOM[an1]
                anion1_charm = anChArm[an1]
                anion2_com = anCOM[an2]
                anion2_charm = anChArm[an2]
                
                dist = dist2vec(anion1_com,anion2_com)
                
                if dist<distCriteria_AnAn:
                    anan_len_an1[anan_counter] = dist1vec(anion1_charm)
                    anan_len_an2[anan_counter] = dist1vec(anion2_charm)
                    
                    anan_angle[anan_counter] = angle(anion1_charm,anion2_charm)
                    
                    anan_counter+=1
                    
# Remove trailing zeros
catan_len_cat = catan_len_cat[0:catan_counter]
catan_len_an = catan_len_an[0:catan_counter]
catan_angle = catan_angle[0:catan_counter]

catcat_len_cat1 = catcat_len_cat1[0:catcat_counter]
catcat_len_cat2 = catcat_len_cat2[0:catcat_counter]
catcat_angle = catcat_angle[0:catcat_counter]

anan_len_an1 = anan_len_an1[0:anan_counter]
anan_len_an2 = anan_len_an2[0:anan_counter]
anan_angle = anan_angle[0:anan_counter]

# Convert angle data to degrees
catan_angle = catan_angle/np.pi*180
catcat_angle = catcat_angle/np.pi*180
anan_angle = anan_angle/np.pi*180

# Histogram angle data
catan_angle_hist, catan_angle_binedges = np.histogram(catan_angle, bins=180, range=(0,180))
catcat_angle_hist, catcat_angle_binedges = np.histogram(catcat_angle, bins=180, range=(0,180))
anan_angle_hist, anan_angle_binedges = np.histogram(anan_angle, bins=180, range=(0,180))
#Manually make bin centres
angle_bin_centres = np.linspace(0.5,179.5,num=180)

# Apply cone correction to angle histograms
catan_angle_norm = ConeCorrectionVector/np.mean(ConeCorrectionVector)*np.mean(catan_angle_hist)
catan_angle_hist = catan_angle_hist/catan_angle_norm

catcat_angle_norm = ConeCorrectionVector/np.mean(ConeCorrectionVector)*np.mean(catcat_angle_hist)
catcat_angle_hist = catcat_angle_hist/catcat_angle_norm

anan_angle_norm = ConeCorrectionVector/np.mean(ConeCorrectionVector)*np.mean(anan_angle_hist)
anan_angle_hist = anan_angle_hist/anan_angle_norm

# Histogram length data
catan_len_zhist, catan_len_xbinedges, catan_len_ybinedges = np.histogram2d(catan_len_cat,catan_len_an,bins=50)
catcat_len_zhist, catcat_len_xbinedges, catcat_len_ybinedges = np.histogram2d(catcat_len_cat1,catcat_len_cat2,bins=50)
anan_len_zhist, anan_len_xbinedges, anan_len_ybinedges = np.histogram2d(anan_len_an1,anan_len_an2,bins=50)

# Turn bin edges into bin centres
catan_len_xbincentre = np.zeros(len(catan_len_xbinedges)-1)
catan_len_ybincentre = np.zeros(len(catan_len_ybinedges)-1)
catcat_len_xbincentre = np.zeros(len(catcat_len_xbinedges)-1)
catcat_len_ybincentre = np.zeros(len(catcat_len_ybinedges)-1)
anan_len_xbincentre = np.zeros(len(anan_len_xbinedges)-1)
anan_len_ybincentre = np.zeros(len(anan_len_ybinedges)-1)

for i in range(len(catan_len_xbincentre)):
    catan_len_xbincentre[i]=(catan_len_xbinedges[i]+catan_len_xbinedges[i+1])/2
for i in range(len(catan_len_ybincentre)):
    catan_len_ybincentre[i]=(catan_len_ybinedges[i]+catan_len_ybinedges[i+1])/2
for i in range(len(catcat_len_xbincentre)):
    catcat_len_xbincentre[i]=(catcat_len_xbinedges[i]+catcat_len_xbinedges[i+1])/2
for i in range(len(catcat_len_ybincentre)):
    catcat_len_ybincentre[i]=(catcat_len_ybinedges[i]+catcat_len_ybinedges[i+1])/2
for i in range(len(anan_len_xbincentre)):
    anan_len_xbincentre[i]=(anan_len_xbinedges[i]+anan_len_xbinedges[i+1])/2
for i in range(len(anan_len_ybincentre)):
    anan_len_ybincentre[i]=(anan_len_ybinedges[i]+anan_len_ybinedges[i+1])/2

# Save data to text files
catan_len_filename = "Cation-Anion_Length_"+str(distCriteria_CatAn)+"pm.txt"
fileopen = open(catan_len_filename,"w")
fileopen.write("\t\t;\t")
for i in range(len(catan_len_xbincentre)):
    fileopen.write("%f;\t"%(catan_len_xbincentre[i]))
fileopen.write("\n")
for i in range(len(catan_len_ybincentre)):
    fileopen.write("%f;\t"%(catan_len_ybincentre[i]))
    for j in range(len(catan_len_xbincentre)):
        fileopen.write("%f;\t"%(catan_len_zhist[j,i]))
    fileopen.write("\n")
    
catcat_len_filename = "Cation-Cation_Length_"+str(distCriteria_CatCat)+"pm.txt"
fileopen = open(catcat_len_filename,"w")
fileopen.write("\t\t;\t")
for i in range(len(catcat_len_xbincentre)):
    fileopen.write("%f;\t"%(catcat_len_xbincentre[i]))
fileopen.write("\n")
for i in range(len(catcat_len_ybincentre)):
    fileopen.write("%f;\t"%(catcat_len_ybincentre[i]))
    for j in range(len(catcat_len_xbincentre)):
        fileopen.write("%f;\t"%(catcat_len_zhist[j,i]))
    fileopen.write("\n")
    
anan_len_filename = "Anion-Anion_Length_"+str(distCriteria_AnAn)+"pm.txt"
fileopen = open(anan_len_filename,"w")
fileopen.write("\t\t;\t")
for i in range(len(anan_len_xbincentre)):
    fileopen.write("%f;\t"%(anan_len_xbincentre[i]))
fileopen.write("\n")
for i in range(len(anan_len_ybincentre)):
    fileopen.write("%f;\t"%(anan_len_ybincentre[i]))
    for j in range(len(anan_len_xbincentre)):
        fileopen.write("%f;\t"%(anan_len_zhist[j,i]))
    fileopen.write("\n")
    
catan_ang_filename = "Cation-Anion_Angle_"+str(distCriteria_CatAn)+"pm.txt"
fileopen = open(catan_ang_filename,"w")
for i in range(len(angle_bin_centres)):
    fileopen.write("%f;\t%f\n"%(angle_bin_centres[i],catan_angle_hist[i]))
fileopen.close()

catcat_ang_filename = "Cation-Cation_Angle_"+str(distCriteria_CatCat)+"pm.txt"
fileopen = open(catcat_ang_filename,"w")
for i in range(len(angle_bin_centres)):
    fileopen.write("%f;\t%f\n"%(angle_bin_centres[i],catcat_angle_hist[i]))
fileopen.close()

anan_ang_filename = "Anion-Anion_Angle_"+str(distCriteria_AnAn)+"pm.txt"
fileopen = open(anan_ang_filename,"w")
for i in range(len(angle_bin_centres)):
    fileopen.write("%f;\t%f\n"%(angle_bin_centres[i],anan_angle_hist[i]))
fileopen.close()
    
# Create and save plots
plt.hist2d(catan_len_cat,catan_len_an,bins=50)
plt.title("Cation length vs anion length")
plt.colorbar()
plt.xlabel("Length [pm]")
plt.ylabel("Length [pm]")
plt.tight_layout()
catan_len_savename = "Cation-Anion_Length_"+str(distCriteria_CatAn)+"pm.png"
plt.savefig(catan_len_savename)
plt.close()

plt.hist2d(catcat_len_cat1,catcat_len_cat2,bins=50)
plt.title("Cation length vs cation length")
plt.colorbar()
plt.xlabel("Length [pm]")
plt.ylabel("Length [pm]")
plt.tight_layout()
catcat_len_savename = "Cation-Cation_Length_"+str(distCriteria_CatCat)+"pm.png"
plt.savefig(catcat_len_savename)
plt.close()

plt.hist2d(anan_len_an1,anan_len_an2,bins=50)
plt.title("Anion length vs anion length")
plt.colorbar()
plt.xlabel("Length [pm]")
plt.ylabel("Length [pm]")
plt.tight_layout()
anan_len_savename = "Anion-Anion_Length_"+str(distCriteria_AnAn)+"pm.png"
plt.savefig(anan_len_savename)
plt.close()

plt.plot(angle_bin_centres,catan_angle_hist)
plt.title("Cation vs anion angle")
plt.xlabel("Angle [degree]")
plt.ylabel("Relative Occurance")
plt.tight_layout()
catan_ang_savename = "Cation-Anion_Angle_"+str(distCriteria_CatAn)+"pm.png"
plt.savefig(catan_ang_savename)
plt.close()

plt.plot(angle_bin_centres,catcat_angle_hist)
plt.title("Cation vs cation angle")
plt.xlabel("Angle [degree]")
plt.ylabel("Relative Occurance")
plt.tight_layout()
catcat_ang_savename = "Cation-Cation_Angle_"+str(distCriteria_CatCat)+"pm.png"
plt.savefig(catcat_ang_savename)
plt.close()


plt.plot(angle_bin_centres,anan_angle_hist)
plt.title("Anion vs anion angle")
plt.xlabel("Angle [degree]")
plt.ylabel("Relative Occurance")
plt.tight_layout()
anan_ang_savename = "Anion-Anion_Angle_"+str(distCriteria_AnAn)+"pm.png"
plt.savefig(anan_ang_savename)
plt.close()








































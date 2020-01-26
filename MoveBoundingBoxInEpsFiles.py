'''
Edit .eps files to allow compatability with LaTeX

Known bug where having the bounding box at the end of an eps file makes it 
unrecognisable to eps2pdf.

This script fins every .eps file within the diretory and all subdirectories 
and moves the bounding box to the top of the file if it is at the bottom.
'''

import os

path = '.'
files = []

# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if '.eps' in file:
            files.append(os.path.join(r, file))
            

for eachFile in files:

    ## Resets edit trigger
    needsEditing=False

    test=open(eachFile,"r").readlines()
    ## test is list with each line as an item
    
    for item in test:
        ## Collect right string (only initiates after needsEditing is triggered)
        if needsEditing==True: 
            if item[0:13]=="%%BoundingBox":
                    newString=item
        ## Checks to see if file needs editing
        if item == "%%BoundingBox: (atend)\n":
            needsEditing=True
    
    ## Edit file if it needs editing
    if needsEditing==True:
        print "Editing file: %s" %(eachFile)
        s = open(eachFile).read()
        s = s.replace("%%BoundingBox: (atend)\n", newString)
        f = open(eachFile, 'w')
        f.write(s)
        f.close()

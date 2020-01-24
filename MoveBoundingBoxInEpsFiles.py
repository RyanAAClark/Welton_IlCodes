import os

path = '.'
files = []

# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if '.eps' in file:
            files.append(os.path.join(r, file))

#for f in files:
#    print(f)

for eachFile in files:

    ## Resets edit trigger
    needsEditing=False

    test=open(eachFile,"r").readlines()
    ## test is list with each line as an item
    
    for item in test:
        ## Collect right string (only initiates after needsEditing is triggered
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

#Dict of dict, each key is a class, with minimum and maximum retention values as second set of keys
#classes = ['CE', 'TG', 'DG', 'FC', 'MG', 'CER', 'HexCer', 'HexCer(OH)', 'PG', 'PE', 'PC', 'SM', 'LPC']

abundance = {}
compoundCount = {}
sampleCompoundCount = {}
sampleCount = 32 #Number of samples
minR = {} #Min. retention time per class
maxR = {} #Max. retention time per class
minM = {} #Min. mass per class
maxM = {} #Max. mass per class
resFac = {} #Response factor per class

f = open('data.txt')
splitLine = f.readline().rstrip().split('\t')
sampleCount = len(splitLine)-3


f = open("ClassSpecs.txt") #Reads file with retention times and mass for each class and loads them into memory
f.readline()
for line in f:
    splitLine = line.rstrip().split('\t')
    abundance[splitLine[1]] = []
    compoundCount[splitLine[1]] = []
    sampleCompoundCount[splitLine[1]] = []
    for i in range(sampleCount):
        abundance[splitLine[1]].append(0)
        compoundCount[splitLine[1]].append(0)
    minR[splitLine[1]] = float(splitLine[2])
    maxR[splitLine[1]] = float(splitLine[3])    
    minM[splitLine[1]] = float(splitLine[4])
    maxM[splitLine[1]] = float(splitLine[5])
    resFac[splitLine[1]] = float(splitLine[6])
    
f = open('data.txt') #Reads in the actual measurement data and counts abundances


splitLine = f.readline().rstrip().split('\t')

sampleGroup = []
groupAbundance = {}
groupCompCount = {}

for i in range(sampleCount):
    sampleGroup.append(splitLine[i+3])

for lipClass in abundance:
    groupAbundance[lipClass] = {}
    groupCompCount[lipClass] = {}
    for i in sampleGroup:
        groupAbundance[lipClass][i] = 0
        groupCompCount[lipClass][i] = 0
        
of = open('compound-abundances.txt', 'w')
splitLine.insert(3, 'Lipid class')
print >> of, '\t'.join(splitLine)

for line in f:
    splitLine = line.rstrip().split('\t')
    retTime = float(splitLine[2])
    mass  = float(splitLine[1])
    for lipClass in abundance:
        if ((retTime > minR[lipClass]) & (retTime < maxR[lipClass])):
            if ((mass > minM[lipClass]) & (mass < maxM[lipClass])):
                for i in range(len(abundance[lipClass])):
                    
                    abundance[lipClass][i] += float(splitLine[i+3])
                    groupAbundance[lipClass][sampleGroup[i]] += float(splitLine[i+3])
                    
                    if ((float(splitLine[i+3]) > 1)):
                        compoundCount[lipClass][i] += 1
                        groupCompCount[lipClass][sampleGroup[i]] += 1
                        
                splitLine.insert(3, lipClass)
                print >> of, '\t'.join(splitLine)
        
of.close()

f = open('LC-abundances.txt', 'w') #Writes output file
line = "Lipid class\t"
for i in range(len(abundance[lipClass])):
    line += "Sample_"+str(i+1)+"_abundance"+"\t"
print >> f, line
for lipClass in abundance:
    line = lipClass+'\t'
    for i in range(len(abundance[lipClass])):
        line += str(abundance[lipClass][i]) + '\t'
    print >> f, line.rstrip()

f.close()


f = open("ClassSpecs.txt") #Reads file with retention times and mass for each class and loads them into memory
f.readline()
for line in f:
    splitLine = line.rstrip().split('\t')
    abundance[splitLine[1]] = []
    compoundCount[splitLine[1]] = []
    sampleCompoundCount[splitLine[1]] = []
    for i in range(sampleCount):
        abundance[splitLine[1]].append(0)
        compoundCount[splitLine[1]].append(0)
    minR[splitLine[1]] = float(splitLine[2])
    maxR[splitLine[1]] = float(splitLine[3])    
    minM[splitLine[1]] = float(splitLine[4])
    maxM[splitLine[1]] = float(splitLine[5])
    resFac[splitLine[1]] = float(splitLine[6])


#Reset tables and do RF-corrected counts
f = open('data.txt') #Reads in the actual measurement data and counts abundances


splitLine = f.readline().rstrip().split('\t')

sampleGroup = []
groupAbundance = {}
groupCompCount = {}

for i in range(sampleCount):
    sampleGroup.append(splitLine[i+3])

for lipClass in abundance:
    groupAbundance[lipClass] = {}
    groupCompCount[lipClass] = {}
    for i in sampleGroup:
        groupAbundance[lipClass][i] = 0
        groupCompCount[lipClass][i] = 0
        


of = open('compound-abundances-RF-corrected.txt', 'w')
splitLine.insert(3, 'Lipid class')
print >> of, '\t'.join(splitLine)

for line in f:
    splitLine = line.rstrip().split('\t')
    retTime = float(splitLine[2])
    mass  = float(splitLine[1])
    for lipClass in abundance:
        if ((retTime > minR[lipClass]) & (retTime < maxR[lipClass])):
            if ((mass > minM[lipClass]) & (mass < maxM[lipClass])):
                outLine = splitLine[:3]
                outLine.append(lipClass)
                for i in range(len(abundance[lipClass])):
                    
                    abundance[lipClass][i] += float(splitLine[i+3])*resFac[lipClass]
                    groupAbundance[lipClass][sampleGroup[i]] += float(splitLine[i+3])*resFac[lipClass]
                    outLine.append(str(float(splitLine[i+3])*resFac[lipClass]))
                    
                    if ((float(splitLine[i+3]) > 1)):
                        compoundCount[lipClass][i] += 1
                        groupCompCount[lipClass][sampleGroup[i]] += 1

                    
                    
                print >> of, '\t'.join(outLine)
        
of.close()

f = open('LC-abundances-RF-corrected.txt', 'w') #Writes output file
line = "Lipid class\t"
for i in range(len(abundance[lipClass])):
    line += "Sample_"+str(i+1)+"_abundance"+"\t"
print >> f, line
for lipClass in abundance:
    line = lipClass+'\t'
    for i in range(len(abundance[lipClass])):
        line += str(abundance[lipClass][i]) + '\t'
    print >> f, line.rstrip()

f.close()


f = open('compoundCounts.txt', 'w')
line = "Lipid class\t"
for i in range(len(abundance[lipClass])):
    line += "Sample "+str(i+1)+"\t"
print >> f, line
for lipClass in abundance:
    line = lipClass+'\t'
    for i in range(len(abundance[lipClass])):
        line += str(compoundCount[lipClass][i]) + '\t'
    print >> f, line.rstrip()

f.close()
 


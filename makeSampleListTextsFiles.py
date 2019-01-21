import os

tag = "p2540"
#inDir = "/export/home/snyder/Stop0LCodeSUSY2345/Stop0LRun2/scripts/"
#inDir = "/direct/usatlas+u/max/Stop0L/Stop0LRun2/data/BNL_PFN/"
inDir = "/direct/usatlas+u/isnyder/Stop0LRun2ASUSY2345/Stop0LRun2/data/BNL_PFN/"
outDir = "../util/sampLists/pnfs/"

fList = [samp for samp in os.listdir(inDir) if tag in samp];
derivation = "SUSY1"

bkgDict = {};
signalDict = {};
dataList = [];

for samp in fList:
    if "data" not in samp:
        if "signal" in samp:
            #signalDict[samp.replace("_", "").replace(tag, "").replace(".txt", "").replace(derivation, "")] = open(samp).read().splitlines();
            signalSamples = open(inDir+samp).read().splitlines()
            
            for i in xrange(len(signalSamples)):
                if signalSamples[i].startswith("#"):
                    signalLines = '_'.join(((signalSamples[i].split(':', 1)[-1]).split('_')[:-2]))
                    signalDict[signalLines] = []
                else: signalDict[signalLines].append(signalSamples[i])
                
        else:
            bkgDict[samp] = open(inDir+samp).read().splitlines(); #.replace("_", "").replace(tag, "").replace(".txt", "").replace(derivation, "")] = open(inDir+samp).read().splitlines();
            bkgSamples = open(inDir+samp).read().splitlines()
                
            for i in xrange(len(bkgSamples)):
                if bkgSamples[i].startswith("#"):
                    bkgLines = '_'.join(((bkgSamples[i].split(':', 1)[-1]).split('_')[:-2]))
                    bkgDict[bkgLines] = []
                else: bkgDict[bkgLines].append(bkgSamples[i])

#print signalDict
for key, value in signalDict.iteritems():
    target = open(outDir+key+".txt", "w")
    for name in value:
        #print name
        target.write(name+"\n")
        
for key, value in bkgDict.iteritems():
    target = open(outDir+key+".txt", "w")
    for name in value:
        target.write(name+"\n")
        

            
#    for samp in fList:
#        if "signal" in samp:
#            signalDict[samp.replace("_", "").replace(tag, "").replace(".txt", "").replace(derivation, "")] = open("sampLists/"+samp).read().splitlines();

#        elif "data" in samp:
#            dataList += open("sampLists/"+samp).read().splitlines();

#        else:
#            bkgDict[samp.replace("_", "").replace(tag, "").replace(".txt", "").replace(derivation, "")] = open("sampLists/"+samp).read().splitlines();

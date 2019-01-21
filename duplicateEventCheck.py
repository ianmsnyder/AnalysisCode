from ROOT import *
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--inDir", help="Where your file is", default="./")
parser.add_argument("--inF", help="Which file you're checking", default="")

args = parser.parse_args()

#inDir = "/usatlas/groups/bnl_local/isnyder/samples/Stop0LRun2-00-01-01/noSystematics/"
#inF = "data15Ext.root"

checkFile = TFile(args.inDir+args.inF)

datatree = checkFile.Get("stop_0LepExt")
#datatree.Scan("RunNumber : EventNumber")

runEventDict = {}

for i in xrange(datatree.GetEntries()):
    datatree.GetEntry(i)
    runNumber = datatree.RunNumber
    eventNumber = datatree.EventNumber

    #print runNumber, eventNumber

    if str(runNumber) not in runEventDict.keys():
        runEventDict[str(runNumber)]=[]
        runEventDict[str(runNumber)].append(eventNumber)
    else:
        runEventDict[str(runNumber)].append(eventNumber)

    if i%1000==0:
        print "Event: ", i, " of ", datatree.GetEntries()
#        print runEventDict.values()
        
print "Performing first check: "
for key, value in runEventDict.iteritems():
    if len(value) != len(set(value)):
        print "There are duplicate events in run number ", key, "  ", len(value), len(set(value))
    else:
        print "All good", len(value), len(set(value))
#Another way to check:
print "Performing second check: "
for key, value in runEventDict.iteritems():
    for val in value:
        #print value.count(val)
        if value.count(val)!=1:
            print "Repeat in: ", key

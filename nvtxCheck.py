import os, sys, code, math
from ROOT import *
from array import array
import plotSettings
import plotting
from plotting import oneDplotter, twoDplotter, profileplotter
import pickle, csv
import collections

#from plotting import oneDplotter

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

##############################################################
##### Define command line arguements  ########################
##############################################################
import argparse
parser = argparse.ArgumentParser()
sys.path.append("./plotting")
from varNamesUnits import varLabels
from setStyle import setStyle
from cutStrings import preCutSR, weightStr, lepWeightStr, srCuts, preCut2LepCR

setStyle()

parser.add_argument("--outDir", help="Where output files go", default="./nvtxChecks/")

parser.add_argument("--inDir", help="Where input root files are", default="/usatlas/groups/bnl_local2/whopkins/samples/Stop0LRun2-00-01-19/")

args = parser.parse_args()
datadict={}
datadict2lep={}
mcdict={}
mcdict2lep={}


dataFileD = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root") 
datatreePeriodD = dataFileD.Get("stop_0Lep")
datadict["datatreePeriodD"]=dataFileD.Get("stop_0Lep")
datadict2lep["datatreePeriodD"]=dataFileD.Get("stop_2Lep")

dataFileE = TFile(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
datatreePeriodE = dataFileE.Get("stop_0Lep")
datadict["datatreePeriodE"] = dataFileE.Get("stop_0Lep")
datadict2lep["datatreePeriodE"] = dataFileE.Get("stop_2Lep")
    
dataFileF = TFile(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
datatreePeriodF = dataFileF.Get("stop_0Lep")
datadict["datatreePeriodF"] = dataFileF.Get("stop_0Lep")
datadict2lep["datatreePeriodF"] = dataFileF.Get("stop_2Lep")
    
dataFileG = TFile(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
datatreePeriodG = dataFileG.Get("stop_0Lep")
datadict["datatreePeriodG"] = dataFileG.Get("stop_0Lep")
datadict2lep["datatreePeriodG"] = dataFileG.Get("stop_2Lep")
    
dataFileH = TFile(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
datatreePeriodH = dataFileH.Get("stop_0Lep")
datadict["datatreePeriodH"] = dataFileH.Get("stop_0Lep")
datadict2lep["datatreePeriodH"] = dataFileH.Get("stop_2Lep")
    
dataFileJ = TFile(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
datatreePeriodJ = dataFileJ.Get("stop_0Lep")
datadict["datatreePeriodJ"] = dataFileJ.Get("stop_0Lep")
datadict2lep["datatreePeriodJ"] = dataFileJ.Get("stop_2Lep")

datafile16PeriodA = TFile(args.inDir+"data16PeriodA_SUSY1_p2667.root")
datatree16PeriodA = datafile16PeriodA.Get("stop_0Lep")
datadict["datatree16PeriodA"] = datafile16PeriodA.Get("stop_0Lep")
datadict2lep["datatree16PeriodA"] = datafile16PeriodA.Get("stop_2Lep")

datafile16PeriodB = TFile(args.inDir+"data16PeriodB123_SUSY1_p2667.root")
datatree16PeriodB = datafile16PeriodB.Get("stop_0Lep")
datadict["datatree16PeriodB"] = datafile16PeriodB.Get("stop_0Lep")
datadict2lep["datatree16PeriodB"] = datafile16PeriodB.Get("stop_2Lep")

datafile16PeriodC = TFile(args.inDir+"data16PeriodC_SUSY1_p2667.root")
datatree16PeriodC = datafile16PeriodC.Get("stop_0Lep")
datadict["datatree16PeriodC"] = datafile16PeriodC.Get("stop_0Lep")
datadict2lep["datatree16PeriodC"] = datafile16PeriodC.Get("stop_2Lep")

datafile16PeriodD = TFile(args.inDir+"data16PeriodD_SUSY1_p2689.root")
datatree16PeriodD = datafile16PeriodD.Get("stop_0Lep")
datadict["datatree16PeriodD"] = datafile16PeriodD.Get("stop_0Lep")
datadict2lep["datatree16PeriodD"] = datafile16PeriodD.Get("stop_2Lep")

datafile16PeriodE = TFile(args.inDir+"data16PeriodE_SUSY1_p2689.root")
datatree16PeriodE = datafile16PeriodE.Get("stop_0Lep")
datadict["datatree16PeriodE"] = datafile16PeriodE.Get("stop_0Lep")
datadict2lep["datatree16PeriodE"] = datafile16PeriodE.Get("stop_2Lep")

datafile16PeriodF = TFile(args.inDir+"data16PeriodF_SUSY1_p2689.root")
datatree16PeriodF = datafile16PeriodF.Get("stop_0Lep")
datadict["datatree16PeriodF"] = datafile16PeriodF.Get("stop_0Lep")
datadict2lep["datatree16PeriodF"] = datafile16PeriodF.Get("stop_2Lep")

datafile16PeriodG = TFile(args.inDir+"data16PeriodG_SUSY1_p2769.root")
datatree16PeriodG = datafile16PeriodG.Get("stop_0Lep")
datadict["datatree16PeriodG"] = datafile16PeriodG.Get("stop_0Lep")
datadict2lep["datatree16PeriodG"] = datafile16PeriodG.Get("stop_2Lep")

datafile16PeriodI = TFile(args.inDir+"data16PeriodI_SUSY1_p2769.root")
datatree16PeriodI = datafile16PeriodI.Get("stop_0Lep")
datadict["datatree16PeriodI"] = datafile16PeriodI.Get("stop_0Lep")
datadict2lep["datatree16PeriodI"] = datafile16PeriodI.Get("stop_2Lep")

datafile16PeriodK = TFile(args.inDir+"data16PeriodK_SUSY1_p2769.root")
datatree16PeriodK = datafile16PeriodK.Get("stop_0Lep")
datadict["datatree16PeriodK"] = datafile16PeriodK.Get("stop_0Lep")
datadict2lep["datatree16PeriodK"] = datafile16PeriodK.Get("stop_2Lep")

datafile16PeriodL = TFile(args.inDir+"data16PeriodL_SUSY1_p2840.root")
datatree16PeriodL = datafile16PeriodL.Get("stop_0Lep")
datadict["datatree16PeriodL"] = datafile16PeriodL.Get("stop_0Lep")
datadict2lep["datatree16PeriodL"] = datafile16PeriodL.Get("stop_2Lep")

####  All the background files

datatreettbarFile = TFile(args.inDir+"ttbarExt.root")
datatreettbar=datatreettbarFile.Get("stop_0LepExt")
mcdict["datatreettbar"] = datatreettbarFile.Get("stop_0LepExt")
mcdict2lep["datatreettbar"] = datatreettbarFile.Get("stop_2LepExt")

datatreeZFile = TFile(args.inDir+"ZExt.root")
datatreeZ=datatreeZFile.Get("stop_0LepExt")
mcdict["datatreeZ"] = datatreeZFile.Get("stop_0LepExt")
mcdict2lep["datatreeZ"] = datatreeZFile.Get("stop_2LepExt")

datatreeWFile = TFile(args.inDir+"WExt.root")
datatreeW=datatreeWFile.Get("stop_0LepExt")
mcdict["datatreeW"] = datatreeWFile.Get("stop_0LepExt")
mcdict2lep["datatreeW"] = datatreeWFile.Get("stop_2LepExt")

datatreettVFile = TFile(args.inDir+"ttVExt.root")
datatreettV=datatreettVFile.Get("stop_0LepExt")
mcdict["datatreettV"] = datatreettVFile.Get("stop_0LepExt")
mcdict2lep["datatreettV"] = datatreettVFile.Get("stop_2LepExt")

datatreesingleTopFile = TFile(args.inDir+"singleTopExt.root")
datatreesingleTop=datatreesingleTopFile.Get("stop_0LepExt")
mcdict["datatreesingleTop"] = datatreesingleTopFile.Get("stop_0LepExt")
mcdict2lep["datatreesingleTop"] = datatreesingleTopFile.Get("stop_2LepExt")

datatreedibosonsFile = TFile(args.inDir+"dibosonsExt.root")
datatreedibosons=datatreedibosonsFile.Get("stop_0LepExt")
mcdict["datatreedibosons"] = datatreedibosonsFile.Get("stop_0LepExt")
mcdict2lep["datatreedibosons"] = datatreedibosonsFile.Get("stop_2LepExt")

#####  2-lep trees:
datatreettbar2lep=datatreettbarFile.Get("stop_2LepExt")
datatreeZ2lep=datatreeZFile.Get("stop_2LepExt")
datatreeW2lep=datatreeWFile.Get("stop_2LepExt")
datatreettV2lep=datatreettVFile.Get("stop_2LepExt")
datatreesingleTop2lep=datatreesingleTopFile.Get("stop_2LepExt")
datatreedibosons2lep=datatreedibosonsFile.Get("stop_2LepExt")



datatreePeriodD2lep = dataFileD.Get("stop_2Lep")

datatreePeriodE2lep = dataFileE.Get("stop_2Lep")
    
datatreePeriodF2lep = dataFileF.Get("stop_2Lep")
    
datatreePeriodG2lep = dataFileG.Get("stop_2Lep")
    
datatreePeriodH2lep = dataFileH.Get("stop_2Lep")
    
datatreePeriodJ2lep = dataFileJ.Get("stop_2Lep")

datatree16PeriodA2lep = datafile16PeriodA.Get("stop_2Lep")

datatree16PeriodB2lep = datafile16PeriodB.Get("stop_2Lep")

datatree16PeriodC2lep = datafile16PeriodC.Get("stop_2Lep")

datatree16PeriodD2lep = datafile16PeriodD.Get("stop_2Lep")


datatree16PeriodE2lep = datafile16PeriodE.Get("stop_2Lep")

datatree16PeriodF2lep = datafile16PeriodF.Get("stop_2Lep")

datatree16PeriodG2lep = datafile16PeriodG.Get("stop_2Lep")

datatree16PeriodI2lep = datafile16PeriodI.Get("stop_2Lep")

datatree16PeriodK2lep = datafile16PeriodK.Get("stop_2Lep")

datatree16PeriodL2lep = datafile16PeriodL.Get("stop_2Lep")

######

branchNames = []
datatreettbarBranchObject = datatreettbar.GetListOfBranches()
for m in range(len(datatreettbarBranchObject)):
    branchNames.append(datatreettbarBranchObject[m].GetName())

#This file is friended with the main branch, so lets get its list of branches

friendObject = datatreettbar.GetListOfFriends()[0]
friendTree = friendObject.GetTree()
friendBranchObject = friendTree.GetListOfBranches()
for i in range(len(friendBranchObject)):
    branchNames.append(friendBranchObject[i].GetName())


usedBranches = []

xLabels = {}
fileNames = {}

MinBins = {}
MaxBins = {}
NumBins = {}

profileYaxismin = {}
profileYaxismax = {}

for i in xrange(len(branchNames)):
    if branchNames[i]=='Met':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["Met"]
        fileNames["%s"%branchNames[i]] = "Met"
        MinBins["%s"%branchNames[i]] = -100
        MaxBins["%s"%branchNames[i]] = 1500
        NumBins["%s"%branchNames[i]] = 160
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 750
    if branchNames[i]=='NJets':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["NJets"]
        fileNames["%s"%branchNames[i]] = "NJets"
        MinBins["%s"%branchNames[i]] = 0
        MaxBins["%s"%branchNames[i]] = 20
        NumBins["%s"%branchNames[i]] = 20
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 12
    if branchNames[i]=='NBJets':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["NBJets"]
        fileNames["%s"%branchNames[i]] = "NBJets"
        MinBins["%s"%branchNames[i]] = 0
        MaxBins["%s"%branchNames[i]] = 20
        NumBins["%s"%branchNames[i]] = 20
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 4
    if branchNames[i]=='MtBMin':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["MtBMin"]
        fileNames["%s"%branchNames[i]] = "MtBMin"
        MinBins["%s"%branchNames[i]] = 0
        MaxBins["%s"%branchNames[i]] = 500
        NumBins["%s"%branchNames[i]] = 100
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 800
    if branchNames[i]=='AntiKt12M':
        usedBranches.append(branchNames[i]+"[0]")
        xLabels["%s"%branchNames[i]+"[0]"] = varLabels["AntiKt12M[0]"]
        fileNames["%s"%branchNames[i]+"[0]"] = "RC12MassLeading"
        MinBins["%s"%branchNames[i]+"[0]"] = -100
        MaxBins["%s"%branchNames[i]+"[0]"] = 1500
        NumBins["%s"%branchNames[i]+"[0]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[0]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[0]"] = 400
    if branchNames[i]=='AntiKt8M':
        usedBranches.append(branchNames[i]+"[0]")
        xLabels["%s"%branchNames[i]+"[0]"] =  varLabels["AntiKt8M[0]"]
        fileNames["%s"%branchNames[i]+"[0]"] = "RC8MassLeading"
        MinBins["%s"%branchNames[i]+"[0]"] = -100
        MaxBins["%s"%branchNames[i]+"[0]"] = 1500
        NumBins["%s"%branchNames[i]+"[0]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[0]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[0]"] = 160
    if branchNames[i]=='JetM':
        usedBranches.append(branchNames[i]+"[0]")
        xLabels["%s"%branchNames[i]+"[0]"] = "RC4Mass (Leading)"
        fileNames["%s"%branchNames[i]+"[0]"] = "RC4MassLeading"
        MinBins["%s"%branchNames[i]+"[0]"] = -100
        MaxBins["%s"%branchNames[i]+"[0]"] = 1500
        NumBins["%s"%branchNames[i]+"[0]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[0]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[0]"] = 50
    if branchNames[i]=='AntiKt12M':
        usedBranches.append(branchNames[i]+"[1]")
        xLabels["%s"%branchNames[i]+"[1]"] = varLabels["AntiKt12M[1]"]
        fileNames["%s"%branchNames[i]+"[1]"] = "RC12MassSubleading"
        MinBins["%s"%branchNames[i]+"[1]"] = -100
        MaxBins["%s"%branchNames[i]+"[1]"] = 1500
        NumBins["%s"%branchNames[i]+"[1]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[1]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[1]"] = 100
    if branchNames[i]=='AntiKt8M':
        usedBranches.append(branchNames[i]+"[1]")
        xLabels["%s"%branchNames[i]+"[1]"] =  varLabels["AntiKt8M[1]"]
        fileNames["%s"%branchNames[i]+"[1]"] = "RC8MassSubleading"
        MinBins["%s"%branchNames[i]+"[1]"] = -100
        MaxBins["%s"%branchNames[i]+"[1]"] = 1500
        NumBins["%s"%branchNames[i]+"[1]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[1]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[1]"] = 60
    if branchNames[i]=='JetM':
        usedBranches.append(branchNames[i]+"[1]")
        xLabels["%s"%branchNames[i]+"[1]"] = "RC4Mass (Subleading)"
        fileNames["%s"%branchNames[i]+"[1]"] = "RC4MassSubleading"
        MinBins["%s"%branchNames[i]+"[1]"] = -100
        MaxBins["%s"%branchNames[i]+"[1]"] = 1500
        NumBins["%s"%branchNames[i]+"[1]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[1]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[1]"] = 60
    if branchNames[i]=='Ht':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["Ht"]
        fileNames["%s"%branchNames[i]] = "Ht"
        MinBins["%s"%branchNames[i]] = -100
        MaxBins["%s"%branchNames[i]] = 1500
        NumBins["%s"%branchNames[i]] = 160
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 2000
    if branchNames[i]=='HtSig':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = varLabels["HtSig"]
        fileNames["%s"%branchNames[i]] = "HtSig"
        MinBins["%s"%branchNames[i]] = 0
        MaxBins["%s"%branchNames[i]] = 100
        NumBins["%s"%branchNames[i]] = 100
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 30
    if branchNames[i]=='MetSig':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = "MetSig"
        fileNames["%s"%branchNames[i]] = "MetSig"
        MinBins["%s"%branchNames[i]] = 0
        MaxBins["%s"%branchNames[i]] = 70
        NumBins["%s"%branchNames[i]] = 7        
        NumBins["%s"%branchNames[i]] = 160
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 40
    if branchNames[i]=='MetSumEt':
        usedBranches.append(branchNames[i])
        xLabels["%s"%branchNames[i]] = "MetSumEt"
        fileNames["%s"%branchNames[i]] = "MetSumEt"
        MinBins["%s"%branchNames[i]] = -100
        MaxBins["%s"%branchNames[i]] = 1500
        NumBins["%s"%branchNames[i]] = 160
        profileYaxismin["%s"%branchNames[i]] = 0
        profileYaxismax["%s"%branchNames[i]] = 1600
    if branchNames[i]=='JetMt':
        usedBranches.append(branchNames[i]+"[NJets-1]")
        xLabels["%s"%branchNames[i]+"[NJets-1]"] = "mT(MET, softest jet)"
        fileNames["%s"%branchNames[i]+"[NJets-1]"] = "MtJetSoftJet"
        MinBins["%s"%branchNames[i]+"[NJets-1]"] = -100
        MaxBins["%s"%branchNames[i]+"[NJets-1]"] = 1500
        NumBins["%s"%branchNames[i]+"[NJets-1]"] = 160
        profileYaxismin["%s"%branchNames[i]+"[NJets-1]"] = 0
        profileYaxismax["%s"%branchNames[i]+"[NJets-1]"] = 200
    if branchNames[i]=='JetPt':
        usedBranches.append(branchNames[i]+"[0]")
        xLabels[branchNames[i]+"[0]"] = varLabels["JetPt[0]"]
        fileNames[branchNames[i]+"[0]"] = "JetPtLeading"
        MinBins[branchNames[i]+"[0]"] = -100
        MaxBins[branchNames[i]+"[0]"] = 1500
        NumBins[branchNames[i]+"[0]"] = 160
        profileYaxismin[branchNames[i]+"[0]"] = 0
        profileYaxismax[branchNames[i]+"[0]"] = 500
    if branchNames[i]=='JetPt':
        usedBranches.append(branchNames[i]+"[1]")
        xLabels[branchNames[i]+"[1]"] = varLabels["JetPt[1]"]
        fileNames[branchNames[i]+"[1]"] = "JetPtSubleading"
        MinBins[branchNames[i]+"[1]"] = -100
        MaxBins[branchNames[i]+"[1]"] = 1500
        NumBins[branchNames[i]+"[1]"] = 160
        profileYaxismin[branchNames[i]+"[1]"] = 0
        profileYaxismax[branchNames[i]+"[1]"] = 500

    else:
        pass

checkDict = {}

for key, value in fileNames.items():
    if value not in checkDict:
        checkDict[value] = [key]
    else:
        checkDict[value].append(key)

checkList = [key for key, values in checkDict.items() if len(values)>1]

if len(checkList)>0:
    print "Repeat values: ", checkList
    sys.exit()

#### Now back to making plots:


drawConditionsDict = {}

#drawConditionsDict[""] = " %s"%dsidStatement
#drawConditionsDict["0Lepton"] = "TopDecayType==0 && AntitopDecayType==0 && %s"%dsidStatement

drawConditionsList = []

for keys, values in drawConditionsDict.iteritems():
    drawConditionsList.append(keys)


#treeNames=[datatreettbar]

#for j in range(len(treeNames)):
for k in range(len(usedBranches)):
    for l in range(len(drawConditionsList)):
        print "BRANCH:  ", usedBranches[k], len(usedBranches[k])
        # oneDplotter(stackPlot=[1]*len(usedBranches[k]), legendNames = ["Met"],
        #             treeName=[datatreettbar],
        #             drawConditions=[1]*len(usedBranches[k]), drawConditionsStatement=[drawConditionsDict[drawConditionsList[l]]],
        #             branchName=[usedBranches[k]],
        #             fileName=fileNames["%s"%usedBranches[k]]+drawConditionsList[l],
        #             numBins=NumBins["%s"%usedBranches[k]], binMin=MinBins["%s"%usedBranches[k]], binMax=MaxBins["%s"%usedBranches[k]],
        #             xAxisLabel=xLabels["%s"%usedBranches[k]], xAxisMin=MinBins["%s"%usedBranches[k]], xAxisMax=MaxBins["%s"%usedBranches[k]],
        #             yAxisMin=0.01, yAxisMax=100, setLogY=1,  yaxisUnit="GeV",
        #             xAxisOffset = .8, yAxisOffset = .7, legendx1 = 0.5, legendx2 = 0.3, legendy1 = 0.7, legendy2 = 0.9, legendDraw = 1,
        #             drawFactor=[1], drawFactorCond=["CombineWeight"])
# MetLegendList=[]
# MetDict={}
# for j in range(0, 6):#dsidList:
#     MetDict["%s"%dsidList[j]] = "dsid == %s"%dsidList[j]
#     MetLegendList.append(dsidList[j])
xMin = 0
xMax = 50
xNumBins = 10
xLabel = "nvtx"

print usedBranches

# #  #start uncommenting here
# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatreePeriodD, datatreePeriodE, datatreePeriodF, datatreePeriodG, datatreePeriodH, datatreePeriodJ], 
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileDataCombined2015", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Data", "#intL=3.2 fb^{-1}, #sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=preCutSR,
#                 legendDraw=0, legendLabels=["20 GeV Jet p_{T}"]#, "35 GeV Jet p_{T}"]
#     )

# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatree16PeriodA, datatree16PeriodB, datatree16PeriodC, datatree16PeriodD, datatree16PeriodE, datatree16PeriodF, datatree16PeriodG, datatree16PeriodI, datatree16PeriodK, datatree16PeriodL],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileDataCombined2016", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Data", "#intL=33.3 fb^{-1}, #sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=preCutSR,
#                 legendDraw=0, legendLabels=["20 GeV Jet p_{T}"]
#     )


# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatreePeriodD, datatreePeriodE, datatreePeriodF, datatreePeriodG, datatreePeriodH, datatreePeriodJ, datatree16PeriodA, datatree16PeriodB, datatree16PeriodC, datatree16PeriodD, datatree16PeriodE, datatree16PeriodF, datatree16PeriodG, datatree16PeriodI, datatree16PeriodK, datatree16PeriodL],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileAllData", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Data", "#intL=36.47 fb^{-1}, #sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=preCutSR,
#                 legendDraw=0, legendLabels=["20 GeV Jet p_{T}"]
#     )

# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatreePeriodD2lep, datatreePeriodE2lep, datatreePeriodF2lep, datatreePeriodG2lep, datatreePeriodH2lep, datatreePeriodJ2lep, datatree16PeriodA2lep, datatree16PeriodB2lep, datatree16PeriodC2lep, datatree16PeriodD2lep, datatree16PeriodE2lep, datatree16PeriodF2lep, datatree16PeriodG2lep, datatree16PeriodI2lep, datatree16PeriodK2lep, datatree16PeriodL2lep], #[datatreePeriodDJetPt35, datatreePeriodEJetPt35, datatreePeriodFJetPt35, datatreePeriodGJetPt35, datatreePeriodHJetPt35, datatreePeriodJJetPt35]],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileAllData2lep", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Data", "#intL=36.47 fb^{-1}, #sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=preCut2LepCR,
#                 legendDraw=0, legendLabels=["20 GeV Jet p_{T}"]
#     )

# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatreettbar, datatreettV, datatreedibosons, datatreeW, datatreeZ, datatreesingleTop],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileMC", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Internal Simulation", "#sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=weightStr+"*("+preCutSR+")",
#     )

# for k in range(len(usedBranches)):
#     profileplotter(
#                 treeName=[datatreettbar2lep, datatreettV2lep, datatreedibosons2lep, datatreeW2lep, datatreeZ2lep, datatreesingleTop2lep],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileMC2lep", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Internal Simulation", "#sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=lepWeightStr+"*("+preCut2LepCR+")",
#     )




###########################################
########## nvtx distribution#########


def setOverUnderFlow(hist, numberBins):
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinContent(numberBins, hist.GetBinContent(numberBins)+hist.GetBinContent(numberBins+1))
    return hist
    
# weightStrNoPileup=weightStr.replace("PileupWeight*","")
# lepWeightStrNoPileup=lepWeightStr.replace("PileupWeight*","")

# nvtxNBins = 50 
# nvtxBinMin = 0
# nvtxBinMax = 50
# nvtxyaxismax=2500
# nvtxyaxismax2lep=1000

# lumi="36.47"

# nvtxMCnopileup = TH1F("nvtxMCnopileup", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# nvtxMCnopileup.SetLineColor(kBlue)
# datatreettbar.Draw("NVtx>>nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# datatreettV.Draw("NVtx>>+nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# datatreesingleTop.Draw("NVtx>>+nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# datatreedibosons.Draw("NVtx>>+nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# datatreeW.Draw("NVtx>>+nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# datatreeZ.Draw("NVtx>>+nvtxMCnopileup", lumi+"*"+weightStrNoPileup+"*("+preCutSR+")")
# nvtxMCnopileup = setOverUnderFlow(nvtxMCnopileup, nvtxNBins)

# nvtxMCpileup = TH1F("nvtxMCpileup", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# nvtxMCpileup.SetLineColor(kRed)
# datatreettbar.Draw("NVtx>>nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# datatreettV.Draw("NVtx>>+nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# datatreesingleTop.Draw("NVtx>>+nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# datatreedibosons.Draw("NVtx>>+nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# datatreeW.Draw("NVtx>>+nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# datatreeZ.Draw("NVtx>>+nvtxMCpileup", lumi+"*"+weightStr+"*("+preCutSR+")")
# nvtxMCpileup = setOverUnderFlow(nvtxMCpileup, nvtxNBins)

# nvtxData = TH1F("nvtxData", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# datatreePeriodD.Draw("NVtx>>nvtxData", preCutSR)
# datatreePeriodE.Draw("NVtx>>+nvtxData", preCutSR)
# datatreePeriodF.Draw("NVtx>>+nvtxData", preCutSR)
# datatreePeriodG.Draw("NVtx>>+nvtxData", preCutSR)
# datatreePeriodH.Draw("NVtx>>+nvtxData", preCutSR)
# datatreePeriodJ.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodA.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodB.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodC.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodD.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodE.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodF.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodG.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodI.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodL.Draw("NVtx>>+nvtxData", preCutSR)
# datatree16PeriodK.Draw("NVtx>>+nvtxData", preCutSR)
# nvtxData = setOverUnderFlow(nvtxData, nvtxNBins)

# nvtxCanvas = TCanvas("nvtxCanvas", "", 800 ,800)
# nvtxMCnopileup.Draw("hist")
# nvtxMCpileup.Draw("same hist")
# nvtxData.Draw("same hist")


# nvtxMCnopileup.GetXaxis().SetTitle("nvtx")
# nvtxMCnopileup.GetYaxis().SetTitle("Events / "+str((nvtxBinMax-nvtxBinMin)/nvtxNBins))
# nvtxMCnopileup.SetMaximum(nvtxyaxismax)
# TGaxis.SetMaxDigits(2)


# nvtxLegend = TLegend(.65, .7, .85, .9)
# nvtxLegend.SetFillColor(0)
# nvtxLegend.SetEntrySeparation(0.5)
# nvtxLegend.SetTextSize(0.025)
# #nvtxLegend.SetTextSize(0.035)
# nvtxLegend.SetShadowColor(10)
# nvtxLegend.SetLineColor(10)
# nvtxLegend.SetFillColor(10)
# nvtxLegend.SetFillStyle(0)
# nvtxLegend.SetBorderSize(0)
# nvtxLegend.AddEntry(nvtxMCnopileup, "Background, No pileup reweighting", "L")
# nvtxLegend.AddEntry(nvtxMCpileup, "Background, with pileup reweighting", "L")
# nvtxLegend.AddEntry(nvtxData, "Data", "L")
# nvtxLegend.Draw()
# nvtxCanvas.SaveAs(args.outDir+"/nvtxDistribution.eps")

# #######################################################
    
# nvtxMCnopileup2lep = TH1F("nvtxMCnopileup2lep", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# nvtxMCnopileup2lep.SetLineColor(kBlue)
# datatreettbar2lep.Draw("NVtx>>nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# datatreettV2lep.Draw("NVtx>>+nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# datatreesingleTop2lep.Draw("NVtx>>+nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# datatreedibosons2lep.Draw("NVtx>>+nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# datatreeW2lep.Draw("NVtx>>+nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# datatreeZ2lep.Draw("NVtx>>+nvtxMCnopileup2lep", lumi+"*"+lepWeightStrNoPileup+"*("+preCut2LepCR+")")
# nvtxMCnopileup2lep = setOverUnderFlow(nvtxMCnopileup2lep, nvtxNBins)

# nvtxMCpileup2lep = TH1F("nvtxMCpileup2lep", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# nvtxMCpileup2lep.SetLineColor(kRed)
# datatreettbar2lep.Draw("NVtx>>nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# datatreettV2lep.Draw("NVtx>>+nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# datatreesingleTop2lep.Draw("NVtx>>+nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# datatreedibosons2lep.Draw("NVtx>>+nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# datatreeW2lep.Draw("NVtx>>+nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# datatreeZ2lep.Draw("NVtx>>+nvtxMCpileup2lep", lumi+"*"+lepWeightStr+"*("+preCut2LepCR+")")
# nvtxMCpileup2lep = setOverUnderFlow(nvtxMCpileup2lep, nvtxNBins)

# nvtxData2lep = TH1F("nvtxData2lep", "", nvtxNBins, nvtxBinMin, nvtxBinMax)
# datatreePeriodD2lep.Draw("NVtx>>nvtxData2lep", preCut2LepCR)
# datatreePeriodE2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatreePeriodF2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatreePeriodG2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatreePeriodH2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatreePeriodJ2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodA2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodB2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodC2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodD2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodE2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodF2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodG2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodI2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodL2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# datatree16PeriodK2lep.Draw("NVtx>>+nvtxData2lep", preCut2LepCR)
# nvtxData2lep = setOverUnderFlow(nvtxData2lep, nvtxNBins)

# nvtxCanvas2lep = TCanvas("nvtxCanvas2lep", "", 800 ,800)
# nvtxMCnopileup2lep.Draw("hist")
# nvtxMCpileup2lep.Draw("same hist")
# nvtxData2lep.Draw("same hist")


# nvtxMCnopileup2lep.GetXaxis().SetTitle("nvtx")
# nvtxMCnopileup2lep.GetYaxis().SetTitle("Events / "+str((nvtxBinMax-nvtxBinMin)/nvtxNBins))
# nvtxMCnopileup2lep.SetMaximum(nvtxyaxismax2lep)
# TGaxis.SetMaxDigits(2)


# nvtxLegend2lep = TLegend(.65, .7, .85, .9)
# nvtxLegend2lep.SetFillColor(0)
# nvtxLegend2lep.SetEntrySeparation(0.5)
# nvtxLegend2lep.SetTextSize(0.025)
# #nvtxLegend2lep.SetTextSize(0.035)
# nvtxLegend2lep.SetShadowColor(10)
# nvtxLegend2lep.SetLineColor(10)
# nvtxLegend2lep.SetFillColor(10)
# nvtxLegend2lep.SetFillStyle(0)
# nvtxLegend2lep.SetBorderSize(0)
# nvtxLegend2lep.AddEntry(nvtxMCnopileup2lep, "Background, No pileup reweighting", "L")
# nvtxLegend2lep.AddEntry(nvtxMCpileup2lep, "Background, with pileup reweighting", "L")
# nvtxLegend2lep.AddEntry(nvtxData2lep, "Data", "L")
# nvtxLegend2lep.Draw()
# nvtxCanvas2lep.SaveAs(args.outDir+"/nvtxDistribution2lep.eps")

histXnumbins=xNumBins
histXmin=xMin
histXmax=xMax

histdictdata={}
makeRootFiles=True
if makeRootFiles:
    #rootfiles={}
    #trees={}
    for j in xrange(20, 85, 5):
        rootfile=TFile("dataPtCut"+str(j)+".root", "recreate")
        testtree=TTree("njets", "njets")
        #trees[str(j)]=TTree("njets", "njet")
        #trees[str(j)].cd()
        nJets=array("i", [0])
        NVTX=array("f", [0])
        testtree.Branch("nJets", nJets, "nJets/I")
        testtree.Branch("NVTX", NVTX, "NVTX/F")
        #histdictdata[str(j)]=TProfile(str(j),"",  histXnumbins, histXmin, histXmax)
        for key,value in datadict.iteritems():
            print key
            for i in xrange(value.GetEntries()):
                value.GetEntry(i)
                jetpt = value.JetPt
                nvtx = value.NVtx
                l=k=0
                #print list(jetpt)
                #print len(jetpt)
                while l<=(len(jetpt)-1):
                    if jetpt[l]>j:
                        k=l
                        #print k
                        #print k, jetpt[k]
                    l+=1
                #print "K: ", k
                nJets[0]=k
                NVTX[0]=nvtx
                testtree.Fill()
                #trees[str(j)].Fill()
                #histdictdata[str(j)].Fill(nvtx, k, 1.)
        rootfile.cd()
        rootfile.Write()
        rootfile.Close()
        del rootfile, testtree, nJets, NVTX


    for j in xrange(20, 85, 5):
        rootfile=TFile("mcPtCut"+str(j)+".root", "recreate")
        testtree=TTree("njets", "njets")
        #trees[str(j)]=TTree("njets", "njet")
        #trees[str(j)].cd()
        nJets=array("i", [0])
        NVTX=array("f", [0])
        testtree.Branch("nJets", nJets, "nJets/I")
        testtree.Branch("NVTX", NVTX, "NVTX/F")
        #histdictdata[str(j)]=TProfile(str(j),"",  histXnumbins, histXmin, histXmax)
        for key,value in mcdict.iteritems():
            print key
            for i in xrange(value.GetEntries()):
                value.GetEntry(i)
                jetpt = value.JetPt
                nvtx = value.NVtx
                l=k=0
                #print list(jetpt)
                #print len(jetpt)
                while l<=(len(jetpt)-1):
                    if jetpt[l]>j:
                        k=l
                        #print k
                        #print k, jetpt[k]
                    l+=1
                #print "K: ", k
                nJets[0]=k
                NVTX[0]=nvtx
                testtree.Fill()
                #trees[str(j)].Fill()
                #histdictdata[str(j)].Fill(nvtx, k, 1.)
        rootfile.cd()
        rootfile.Write()
        rootfile.Close()
        del rootfile, testtree, nJets, NVTX


    for j in xrange(20, 85, 5):
        rootfile=TFile("dataPtCut"+str(j)+"2lep.root", "recreate")
        testtree=TTree("njets", "njets")
        #trees[str(j)]=TTree("njets", "njet")
        #trees[str(j)].cd()
        nJets=array("i", [0])
        NVTX=array("f", [0])
        testtree.Branch("nJets", nJets, "nJets/I")
        testtree.Branch("NVTX", NVTX, "NVTX/F")
        #histdictdata[str(j)]=TProfile(str(j),"",  histXnumbins, histXmin, histXmax)
        for key,value in datadict2lep.iteritems():
            print key
            for i in xrange(value.GetEntries()):
                value.GetEntry(i)
                jetpt = value.JetPt
                nvtx = value.NVtx
                l=k=0
                #print list(jetpt)
                #print len(jetpt)
                while l<=(len(jetpt)-1):
                    if jetpt[l]>j:
                        k=l
                        #print k
                        #print k, jetpt[k]
                    l+=1
                #print "K: ", k
                nJets[0]=k
                NVTX[0]=nvtx
                testtree.Fill()
                #trees[str(j)].Fill()
                #histdictdata[str(j)].Fill(nvtx, k, 1.)
        rootfile.cd()
        rootfile.Write()
        rootfile.Close()
        del rootfile, testtree, nJets, NVTX


    for j in xrange(20, 85, 5):
        rootfile=TFile("mcPtCut"+str(j)+"2lep.root", "recreate")
        testtree=TTree("njets", "njets")
        #trees[str(j)]=TTree("njets", "njet")
        #trees[str(j)].cd()
        nJets=array("i", [0])
        NVTX=array("f", [0])
        testtree.Branch("nJets", nJets, "nJets/I")
        testtree.Branch("NVTX", NVTX, "NVTX/F")
        #histdictdata[str(j)]=TProfile(str(j),"",  histXnumbins, histXmin, histXmax)
        for key,value in mcdict2lep.iteritems():
            print key
            for i in xrange(value.GetEntries()):
                value.GetEntry(i)
                jetpt = value.JetPt
                nvtx = value.NVtx
                l=k=0
                #print list(jetpt)
                #print len(jetpt)
                while l<=(len(jetpt)-1):
                    if jetpt[l]>j:
                        k=l
                        #print k
                        #print k, jetpt[k]
                    l+=1
                #print "K: ", k
                nJets[0]=k
                NVTX[0]=nvtx
                testtree.Fill()
                #trees[str(j)].Fill()
                #histdictdata[str(j)].Fill(nvtx, k, 1.)
        rootfile.cd()
        rootfile.Write()
        rootfile.Close()
        del rootfile, testtree, nJets, NVTX

njetsDataFiles={}
njetsMCFiles={}
njetsDataTrees={}
njetsMCTrees={}

njetsDataFiles2lep={}
njetsMCFiles2lep={}
njetsDataTrees2lep={}
njetsMCTrees2lep={}

for j in xrange(20, 85, 5):
    njetsDataFiles[str(j)]=TFile("dataPtCut"+str(j)+".root")
    njetsDataTrees[str(j)]=njetsDataFiles[str(j)].Get("njets")
    njetsMCFiles[str(j)]=TFile("mcPtCut"+str(j)+".root")
    njetsMCTrees[str(j)]=njetsMCFiles[str(j)].Get("njets")

    njetsDataFiles2lep[str(j)]=TFile("dataPtCut"+str(j)+".root")
    njetsDataTrees2lep[str(j)]=njetsDataFiles2lep[str(j)].Get("njets")
    njetsMCFiles2lep[str(j)]=TFile("mcPtCut"+str(j)+".root")
    njetsMCTrees2lep[str(j)]=njetsMCFiles2lep[str(j)].Get("njets")

    
# profileplotter(
#                 treeName=[[value for key,value in njetsDataTrees],[value for key,value in njetsMCTrees]],
#                 YbranchName="nJets", XbranchName="NVTX",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileDataCombined2015", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Data", "#intL=3.2 fb^{-1}, #sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=preCutSR,
#                 legendDraw=0, legendLabels=["20 GeV Jet p_{T}"]#, "35 GeV Jet p_{T}"]
#     )
    
code.interact(local=locals())

    
#    dataCanvases[str(j)].cd()
#    histdictdata[str(j)].Draw()
#    dataCanvases[str(j)].SaveAs("test.pdf")
        
    
                    
#profileplotter(
#                 treeName=[datatreettbar, datatreettV, datatreedibosons, datatreeW, datatreeZ, datatreesingleTop],
#                 YbranchName=usedBranches[k], XbranchName="NVtx",
#                 fileName=fileNames["%s"%usedBranches[k]]+"NvtxProfileMC", outDir=args.outDir,
#                 setYaxis=1, histYmin=profileYaxismin["%s"%usedBranches[k]], histYmax=profileYaxismax["%s"%usedBranches[k]],
#                 histXnumbins=xNumBins, histXmin=xMin, histXmax=xMax,
#                 xAxisLabel=xLabel, yAxisLabel=xLabels["%s"%usedBranches[k]],
#                 textBox=1, tPaveTextWords=("#bf{#it{ATLAS}} Internal Simulation", "#sqrt{s}=13TeV"), textx1=0,texty1=profileYaxismax["%s"%usedBranches[k]]-.3*profileYaxismax["%s"%usedBranches[k]],textx2=.7*xMax,texty2=profileYaxismax["%s"%usedBranches[k]],
#                 profileCuts=weightStr+"*("+preCutSR+")",
#     )

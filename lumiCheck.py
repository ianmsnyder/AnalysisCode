
# Author: Ian Snyder
# Directions: first run ilumi2histo.py (provided by background forum) to create one or more ilumi files from ilumicalc files.  Example:
#> python ilumi2histo.py --ifile ../data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-005.root
#Then hadd these files together into a file named ilumi2histo.root and run this script as usual.

from ROOT import *
import argparse, math, sys
from cutStrings import preCutSR, preCut1LepCR, preCut2LepCR, srCuts, cr1LepCuts, cr2LepCuts, cr1PhoCuts, vrCuts, weightStr, lepWeightStr, photonWeightStr



parser = argparse.ArgumentParser()

parser.add_argument("--inDir", help="Where your file is", default="/usatlas/groups/bnl_local2/whopkins/samples/Stop0LRun2-00-01-19/")
parser.add_argument("--inF", help="Which file you're checking", default="")

args = parser.parse_args()

def ErrorCalculator(total, a, aError, b, bError):
    if a!=0. and b!=0.:
        return total*math.sqrt(pow(aError/a,2) + pow(bError/b,2))
    else: return 0.
######################################################
############  SRs  ############################
#####################################

#Trying to combine all SRs, VRs, and CRs into one dictionary.  No weights are used so should be ok

allRegionsDict = srCuts.copy()
allRegionsDict.update(cr1LepCuts)
allRegionsDict.update(cr2LepCuts)
allRegionsDict.update(cr1PhoCuts)
allRegionsDict.update(vrCuts)

for srkey, srvalue in allRegionsDict.iteritems(): #srCuts.iteritems():
    print srkey
    weightString = "1."#weightStr
    cutStr = srvalue
    #print weightString, cutStr
    
    ilumiRootFile = TFile("ilumi2histo.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")

    # ilumiRootFile2015 = TFile("ilumi2hist2015.root")
    # ilumiFile2015 = ilumiRootFile.Get("lumi_histo")
    #print type(ilumiFile2015), type(ilumiFile)

    # testHist = TH1F("testHist", "", 1000000, 0, 1000000)
    # testHist.Add(ilumiFile)
    # testHist.Add(ilumiFile2015)
    # c = TCanvas("c", "")
    # #runHisto0Lep = ilumiFile.Add(ilumiFile2015)
    # #runHisto0Lep.Draw()
    # c.cd()
    # testHist.Draw()
    # c.SaveAs("testRunDist.pdf")

    nBins = ilumiFile.GetNbinsX()
    i=0
    while i<nBins:
        if ilumiFile.GetBinContent(nBins)==0.:
            nBins=nBins-i
            i+=1
#            print i
        else: break
    #print "NBins: ", nBins
    #print "Label : ", ilumiFile.GetBinContent(nBins)
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    if newBins<0:
        print "The number of bins has been calculated to be less than 0.  If you made one ilumi root file by hadd'ing two, make sure you added them in the proper order."
        print "Bins information: NBins, bin min, bin max: ", newBins, binMin, binMax
        sys.exit()
    #print newBins, binMin, binMax

    #print "Number of bins: ", nBins

    lumiDict={}
    lumiErrorDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
        lumiErrorDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinError(i)
    #print lumiErrorDict
    
    if srkey in srCuts:
        datatree0lep = TChain("stop_0Lep")
    if srkey in cr1LepCuts:
        datatree0lep = TChain("stop_1Lep")
    if srkey in cr2LepCuts:
        datatree0lep = TChain("stop_2Lep")
    if srkey in cr1PhoCuts:
        datatree0lep = TChain("stop_1Pho")
    if srkey == "ISRTopVR-2b" or srkey == "ISRTopVR-1b" or srkey == "TopVR1" or srkey ==  "TopVR2":
        datatree0lep = TChain("stop_0Lep")
    if srkey ==  "WVR":
        datatree0lep = TChain("stop_1Lep")
    if srkey ==  "ZVRA" or srkey ==  "ZVRB" or srkey ==  "ZVRC" or srkey ==  "ZVRE" or srkey ==  "ZVRF":
        datatree0lep = TChain("stop_2Lep")
    
        
    # data16PeriodAfile = TFile(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    # data16PeriodBfile = TFile(args.inDir+"data16PeriodB_SUSY1_p2667Ext.root")
    # data16PeriodCfile = TFile(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    # data16PeriodDfile = TFile(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    # data16PeriodEfile = TFile(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    # data16PeriodA0lep = data16PeriodAfile.Get("stop_0Lep")
    # data16PeriodB0lep = data16PeriodBfile.Get("stop_0Lep")
    # data16PeriodC0lep = data16PeriodCfile.Get("stop_0Lep")
    # data16PeriodD0lep = data16PeriodDfile.Get("stop_0Lep")
    # data16PeriodE0lep = data16PeriodEfile.Get("stop_0Lep")

    datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    
    datatree0lep.Add(args.inDir+"data16PeriodA_SUSY1_p2667.root")
    datatree0lep.Add(args.inDir+"data16PeriodB123_SUSY1_p2667.root")
    datatree0lep.Add(args.inDir+"data16PeriodC_SUSY1_p2667.root")
    datatree0lep.Add(args.inDir+"data16PeriodD_SUSY1_p2689.root")
    datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689.root")
    datatree0lep.Add(args.inDir+"data16PeriodF_SUSY1_p2689.root")
    datatree0lep.Add(args.inDir+"data16PeriodG_SUSY1_p2769.root")
    datatree0lep.Add(args.inDir+"data16PeriodI_SUSY1_p2769.root")
    datatree0lep.Add(args.inDir+"data16PeriodK_SUSY1_p2769.root")
    datatree0lep.Add(args.inDir+"data16PeriodL_SUSY1_p2840.root")
    
    #datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    #a = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    #b = a.Get("stop_0Lep")
    #b.Scan("RunNumber")

    #Make TChain of all data sets
    #data16PeriodA0lep.Scan("RunNumber")


    #print binMin, binMax
    #print nBins, "Max: ", type(binMin), type(binMax)
    #print nBins, eval(binMin), eval(binMax)
    #print int(binMin)
    # c = TCanvas("c", "")

    #print "Bins: ", newBins, binMin, binMax
    runHisto0Lep = TH1F("runHisto0Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto0Lep", weightString+"*("+cutStr+")")

    # c.cd()
    # runHisto0Lep.Draw()
    # c.SaveAs("testDist.pdf")
    #datatree0lep.Scan("RunNumber")
    
    dataruns = {}
    datarunsError = {}
    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        #print "entry: ", i
        #print str(runHisto0Lep.GetXaxis().(i+binMin))#GetBinLabel(i))
        #print runHisto0Lep.GetBinContent(i)
        dataruns[str(i+binMin)]=runHisto0Lep.GetBinContent(i+1)
        datarunsError[str(i+binMin)] = math.sqrt(dataruns[str(i+binMin)])
    #print "Dat run number: ", dataruns["298967"], dataruns["298966"], dataruns["298968"], dataruns[str(binMin)], lumiDict[str(binMin)]
    #print dataruns

    # for i in xrange(datatree0lep.GetEntries()):
    #     datatree0lep.GetEntry(i)
    #     if str(datatree0lep.RunNumber) in dataruns.keys():
    #         dataruns[str(datatree0lep.RunNumber)]+=1
    #     else:
    #         dataruns[str(datatree0lep.RunNumber)]=1

    lumiNormDict={}
    for i in lumiDict.keys():
        if i!='':
            if i in dataruns.keys():
                if lumiDict[i]!=0:
                    lumiNormDict[i]=dataruns[i]/lumiDict[i]
                else:
                    lumiNormDict[i]=0
            else: print "Key is missing in dataruns: ", i
            
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ################################################
    #### Repeat for  2015 ########################
    #Memory leak potential?

#     ilumiRootFile = TFile("ilumi2hist2015.root")
#     ilumiFile = ilumiRootFile.Get("lumi_histo")


#     nBins = ilumiFile.GetNbinsX()
#     binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
#     binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
#     newBins = binMax-binMin
#     #print "2015 bins: ", newBins, binMin, binMax

#     #print "Number of bins: ", nBins

#     lumiDict2015={}
#     lumiErrorDict2015={}
#     for i in xrange(nBins):
#         #print i
#         lumiDict2015[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
#         lumiErrorDict2015[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinError(i)
#     #print lumiErrorDict

#     if srkey in srCuts:
#         datatree0lep = TChain("stop_0Lep")
#     if srkey in cr1LepCuts:
#         datatree0lep = TChain("stop_1Lep")
#     if srkey in cr2LepCuts:
#         datatree0lep = TChain("stop_2Lep")
#     if srkey in cr1PhoCuts:
#         datatree0lep = TChain("stop_1Pho")
#     if srkey == "ISRTopVR-2b" or srkey == "ISRTopVR-1b" or srkey == "TopVR1" or srkey ==  "TopVR2":
#         datatree0lep = TChain("stop_0Lep")
#     if srkey ==  "WVR":
#         datatree0lep = TChain("stop_1Lep")
#     if srkey ==  "ZVRA" or srkey ==  "ZVRB" or srkey ==  "ZVRC" or srkey ==  "ZVRE" or srkey ==  "ZVRF":
#         datatree0lep = TChain("stop_2Lep")
        
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
#     datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")

#     dataruns2015 = {}

#     #c = TCanvas("c", "")
#     runHisto0Lep = TH1F("runHisto0Lep", "", newBins, binMin, binMax)
#     datatree0lep.Draw("RunNumber>>runHisto0Lep", weightString+"*("+cutStr+")")
#     #c.cd()
#     #runHisto0Lep.Draw()
#     #c.SaveAs("testDist.pdf")

#     datarunsError2015={}
#     for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
#         dataruns2015[str(i+binMin)]=runHisto0Lep.GetBinContent(i+1)
#         datarunsError2015[str(i+binMin)] = math.sqrt(dataruns2015[str(i+binMin)])


#     lumiNormDict2015={}
#     for i in lumiDict2015.keys():
#         #print dataruns.keys()
#         #print i
#         if lumiDict2015[i]!=0:
#             lumiNormDict2015[i]=dataruns2015[i]/lumiDict2015[i]
#         else:
#             lumiNormDict2015[i]=0
#     #print "Data: ", dataruns
#     #print "Lumi: ", lumiDict
#     #print lumiNormDict

#     ##########################################
#     ####### Combine ###############

#     lumiNormDict = lumiNormDict2015.copy()
#     lumiNormDict.update(lumiNormDict2016)
#     del lumiNormDict['']
#     #print lumiNormDict.keys(), "Min: ", min(lumiNormDict.keys()), "Max: ", max(lumiNormDict.keys())

    newBinMin = int(min(lumiNormDict.keys()))
    newBinMax = int(max(lumiNormDict.keys()))
    newNBins = newBinMax-newBinMin
    #print datarunsError
#    dataruns = dataruns2016.copy()
#    dataruns.update(dataruns2015)
#    lumiDict=lumiDict2016.copy()
#    lumiDict.update(lumiDict2015)
#    datarunsError=datarunsError2016.copy()
#    datarunsError.update(datarunsError2015)
#    lumiErrorDict=lumiErrorDict2016.copy()
#    lumiErrorDict.update(lumiErrorDict2015)

    sortedLumiKeys= sorted(lumiNormDict.keys())
#    print sortedLumiKeys, type(sortedLumiKeys)

#    newHist = TH1F("newHist", "", newNBins, newBinMin, newBinMax)
    nBinsRuns=len(sortedLumiKeys)
    newHist=TH1F("newHist", "", nBinsRuns, 1, nBinsRuns+1)
    for i in xrange(1,nBinsRuns+1):
        binLabel=sortedLumiKeys[i-1]
        binValue=lumiNormDict[binLabel]
        newHist.SetBinContent(i, binValue)
        if i==1:
            newHist.GetXaxis().SetBinLabel(i, binLabel)
        if i==nBinsRuns:
            newHist.GetXaxis().SetBinLabel(i, binLabel)
        if i==len(lumiNormDict.keys())+1:
            newHist.GetXaxis().SetBinLabel(i, binLabel)
        #TGaxis.SetMaxDigits(2);
        #TGaxis.fgMaxDigits = 2
        newHist.SetBinError(i, ErrorCalculator(total=binValue, a=lumiDict[binLabel], aError=lumiErrorDict[binLabel], b=dataruns[binLabel], bError=datarunsError[binLabel]))
    
#    for key, value in lumiNormDict.iteritems():
#        newHist.SetBinContent(int(key)-newBinMin, value)
#        #print value
#        #print lumiDict[key]
#        #print lumiErrorDict[key]
#        #print dataruns[key]
#        #print datarunsError[key])
#        newHist.SetBinError(int(key)-newBinMin, ErrorCalculator(total=value, a=lumiDict[key], aError=lumiErrorDict[key], b=dataruns[key], bError=datarunsError[key]) )

    #for i in xrange(newNBins):
    #    if newHist.GetBinContent(i+1)==0:
    #        newHist.
    c = TCanvas("c","")
    newHist.SetMarkerStyle(20)
    newHist.SetMarkerSize(.5)
    gStyle.SetOptStat(0)
    newHist.Draw("P")
    newHist.GetXaxis().SetTitle("Run Number, "+srkey+"      ")
    newHist.GetYaxis().SetTitle("Lumi-normalized yield")
    c.SaveAs("lumiNormPlots/"+srkey+"LumiNormalized.eps")
    del c, newHist, runHisto0Lep

"""
######################################################
############  VRs  ############################
#####################################

for vrkey, vrvalue in vrCuts.iteritems():
    weightString = "1."#weightStr
    cutStr = vrvalue
    print weightString, cutStr
    
    ilumiRootFile = TFile("ilumi2histo.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")

    ilumiRootFile2015 = TFile("ilumi2hist2015.root")
    ilumiFile2015 = ilumiRootFile.Get("lumi_histo")
    print type(ilumiFile2015), type(ilumiFile)

    # testHist = TH1F("testHist", "", 1000000, 0, 1000000)
    # testHist.Add(ilumiFile)
    # testHist.Add(ilumiFile2015)
    # c = TCanvas("c", "")
    # #runHisto0Lep = ilumiFile.Add(ilumiFile2015)
    # #runHisto0Lep.Draw()
    # c.cd()
    # testHist.Draw()
    # c.SaveAs("testRunDist.pdf")

    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_0Lep")
    # data16PeriodAfile = TFile(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    # data16PeriodBfile = TFile(args.inDir+"data16PeriodB_SUSY1_p2667Ext.root")
    # data16PeriodCfile = TFile(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    # data16PeriodDfile = TFile(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    # data16PeriodEfile = TFile(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    # data16PeriodA0lep = data16PeriodAfile.Get("stop_0Lep")
    # data16PeriodB0lep = data16PeriodBfile.Get("stop_0Lep")
    # data16PeriodC0lep = data16PeriodCfile.Get("stop_0Lep")
    # data16PeriodD0lep = data16PeriodDfile.Get("stop_0Lep")
    # data16PeriodE0lep = data16PeriodEfile.Get("stop_0Lep")

    datatree0lep.Add(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodB123_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    #datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    #a = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    #b = a.Get("stop_0Lep")
    #b.Scan("RunNumber")

    #Make TChain of all data sets
    #data16PeriodA0lep.Scan("RunNumber")

    dataruns = {}
    #print binMin, binMax
    #print nBins, "Max: ", type(binMin), type(binMax)
    #print nBins, eval(binMin), eval(binMax)
    #print int(binMin)
    # c = TCanvas("c", "")

    runHisto0Lep = TH1F("runHisto0Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto0Lep", weightString+"*("+cutStr+")")

    # c.cd()
    # runHisto0Lep.Draw()
    # c.SaveAs("testDist.pdf")
    #datatree0lep.Scan("RunNumber")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        #print "entry: ", i
        #print str(runHisto0Lep.GetXaxis().(i+binMin))#GetBinLabel(i))
        #print runHisto0Lep.GetBinContent(i)
        dataruns[str(i+binMin)]=runHisto0Lep.GetBinContent(i+1)
    #print "Dat run number: ", dataruns["298967"], dataruns["298966"], dataruns["298968"], dataruns[str(binMin)], lumiDict[str(binMin)]
    #print dataruns

    # for i in xrange(datatree0lep.GetEntries()):
    #     datatree0lep.GetEntry(i)
    #     if str(datatree0lep.RunNumber) in dataruns.keys():
    #         dataruns[str(datatree0lep.RunNumber)]+=1
    #     else:
    #         dataruns[str(datatree0lep.RunNumber)]=1

    lumiNormDict2016={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2016[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2016[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ################################################
    #### Repeat for  2015 ########################
    #Memory leak potential?

    ilumiRootFile = TFile("ilumi2hist2015.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")


    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print "2015 bins: ", newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_0Lep")

    datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")

    dataruns = {}

    #c = TCanvas("c", "")
    runHisto0Lep = TH1F("runHisto0Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto0Lep", weightString+"*("+cutStr+")")
    #c.cd()
    #runHisto0Lep.Draw()
    #c.SaveAs("testDist.pdf")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        dataruns[str(i+binMin)]=runHisto0Lep.GetBinContent(i+1)


    lumiNormDict2015={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2015[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2015[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ##########################################
    ####### Combine ###############

    lumiNormDict = lumiNormDict2015.copy()
    lumiNormDict.update(lumiNormDict2016)
    del lumiNormDict['']
    #print lumiNormDict.keys(), "Min: ", min(lumiNormDict.keys()), "Max: ", max(lumiNormDict.keys())

    newBinMin = int(min(lumiNormDict.keys()))
    newBinMax = int(max(lumiNormDict.keys()))
    newNBins = newBinMax-newBinMin

    newHist = TH1F("newHist", "", newNBins, newBinMin, newBinMax)
    for key, value in lumiNormDict.iteritems():
        newHist.SetBinContent(int(key)-newBinMin, value)

    #for i in xrange(newNBins):
    #    if newHist.GetBinContent(i+1)==0:
    #        newHist.
    c = TCanvas("c","")
    newHist.SetMarkerStyle(20)
    newHist.SetMarkerSize(.5)
    gStyle.SetOptStat(0)
    newHist.Draw("P")
    newHist.GetXaxis().SetTitle("Run Number")
    newHist.GetYaxis().SetTitle("Lumi-normalized yield")
    c.SaveAs("lumiNormPlots/"+vrkey+"LumiNormalized.eps")

    
######################################################
############  1 lep CRs  ############################
#####################################

for cr1key, cr1value in cr1LepCuts.iteritems():
    weightString = "1."#lepWeightStr
    cutStr = cr1value
    print weightString, cutStr
    
    ilumiRootFile = TFile("ilumi2histo.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")

    ilumiRootFile2015 = TFile("ilumi2hist2015.root")
    ilumiFile2015 = ilumiRootFile.Get("lumi_histo")
    print type(ilumiFile2015), type(ilumiFile)

    # testHist = TH1F("testHist", "", 1000000, 0, 1000000)
    # testHist.Add(ilumiFile)
    # testHist.Add(ilumiFile2015)
    # c = TCanvas("c", "")
    # #runHisto0Lep = ilumiFile.Add(ilumiFile2015)
    # #runHisto0Lep.Draw()
    # c.cd()
    # testHist.Draw()
    # c.SaveAs("testRunDist.pdf")

    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_1Lep")
    # data16PeriodAfile = TFile(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    # data16PeriodBfile = TFile(args.inDir+"data16PeriodB_SUSY1_p2667Ext.root")
    # data16PeriodCfile = TFile(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    # data16PeriodDfile = TFile(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    # data16PeriodEfile = TFile(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    # data16PeriodA0lep = data16PeriodAfile.Get("stop_1Lep")
    # data16PeriodB0lep = data16PeriodBfile.Get("stop_1Lep")
    # data16PeriodC0lep = data16PeriodCfile.Get("stop_1Lep")
    # data16PeriodD0lep = data16PeriodDfile.Get("stop_1Lep")
    # data16PeriodE0lep = data16PeriodEfile.Get("stop_1Lep")

    datatree0lep.Add(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodB123_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    #datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    #a = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    #b = a.Get("stop_1Lep")
    #b.Scan("RunNumber")

    #Make TChain of all data sets
    #data16PeriodA0lep.Scan("RunNumber")

    dataruns = {}
    #print binMin, binMax
    #print nBins, "Max: ", type(binMin), type(binMax)
    #print nBins, eval(binMin), eval(binMax)
    #print int(binMin)
    # c = TCanvas("c", "")

    runHisto1Lep = TH1F("runHisto1Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto1Lep", weightString+"*("+cutStr+")")

    # c.cd()
    # runHisto1Lep.Draw()
    # c.SaveAs("testDist.pdf")
    #datatree0lep.Scan("RunNumber")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        #print "entry: ", i
        #print str(runHisto1Lep.GetXaxis().(i+binMin))#GetBinLabel(i))
        #print runHisto1Lep.GetBinContent(i)
        dataruns[str(i+binMin)]=runHisto1Lep.GetBinContent(i+1)
    #print "Dat run number: ", dataruns["298967"], dataruns["298966"], dataruns["298968"], dataruns[str(binMin)], lumiDict[str(binMin)]
    #print dataruns

    # for i in xrange(datatree0lep.GetEntries()):
    #     datatree0lep.GetEntry(i)
    #     if str(datatree0lep.RunNumber) in dataruns.keys():
    #         dataruns[str(datatree0lep.RunNumber)]+=1
    #     else:
    #         dataruns[str(datatree0lep.RunNumber)]=1

    lumiNormDict2016={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2016[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2016[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ################################################
    #### Repeat for  2015 ########################
    #Memory leak potential?

    ilumiRootFile = TFile("ilumi2hist2015.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")


    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print "2015 bins: ", newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_1Lep")

    datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")

    dataruns = {}

    #c = TCanvas("c", "")
    runHisto1Lep = TH1F("runHisto1Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto1Lep", weightString+"*("+cutStr+")")
    #c.cd()
    #runHisto1Lep.Draw()
    #c.SaveAs("testDist.pdf")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        dataruns[str(i+binMin)]=runHisto1Lep.GetBinContent(i+1)


    lumiNormDict2015={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2015[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2015[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ##########################################
    ####### Combine ###############

    lumiNormDict = lumiNormDict2015.copy()
    lumiNormDict.update(lumiNormDict2016)
    del lumiNormDict['']
    #print lumiNormDict.keys(), "Min: ", min(lumiNormDict.keys()), "Max: ", max(lumiNormDict.keys())

    newBinMin = int(min(lumiNormDict.keys()))
    newBinMax = int(max(lumiNormDict.keys()))
    newNBins = newBinMax-newBinMin

    newHist = TH1F("newHist", "", newNBins, newBinMin, newBinMax)
    for key, value in lumiNormDict.iteritems():
        newHist.SetBinContent(int(key)-newBinMin, value)

    #for i in xrange(newNBins):
    #    if newHist.GetBinContent(i+1)==0:
    #        newHist.
    c = TCanvas("c","")
    newHist.SetMarkerStyle(20)
    newHist.SetMarkerSize(.5)
    gStyle.SetOptStat(0)
    newHist.Draw("P")
    newHist.GetXaxis().SetTitle("Run Number")
    newHist.GetYaxis().SetTitle("Lumi-normalized yield")
    c.SaveAs("lumiNormPlots/"+cr1key+"LumiNormalized.eps")


######################################################
############  2 lep CRs  ############################
#####################################

for cr2key, cr2value in cr2LepCuts.iteritems():
    weightString = "1."#lepWeightStr
    cutStr = cr2value
    print weightString, cutStr
    
    ilumiRootFile = TFile("ilumi2histo.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")

    ilumiRootFile2015 = TFile("ilumi2hist2015.root")
    ilumiFile2015 = ilumiRootFile.Get("lumi_histo")
    print type(ilumiFile2015), type(ilumiFile)

    # testHist = TH1F("testHist", "", 1000000, 0, 1000000)
    # testHist.Add(ilumiFile)
    # testHist.Add(ilumiFile2015)
    # c = TCanvas("c", "")
    # #runHisto0Lep = ilumiFile.Add(ilumiFile2015)
    # #runHisto0Lep.Draw()
    # c.cd()
    # testHist.Draw()
    # c.SaveAs("testRunDist.pdf")

    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_2Lep")
    # data16PeriodAfile = TFile(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    # data16PeriodBfile = TFile(args.inDir+"data16PeriodB_SUSY1_p2667Ext.root")
    # data16PeriodCfile = TFile(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    # data16PeriodDfile = TFile(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    # data16PeriodEfile = TFile(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    # data16PeriodA0lep = data16PeriodAfile.Get("stop_2Lep")
    # data16PeriodB0lep = data16PeriodBfile.Get("stop_2Lep")
    # data16PeriodC0lep = data16PeriodCfile.Get("stop_2Lep")
    # data16PeriodD0lep = data16PeriodDfile.Get("stop_2Lep")
    # data16PeriodE0lep = data16PeriodEfile.Get("stop_2Lep")

    datatree0lep.Add(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodB123_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    #datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    #a = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    #b = a.Get("stop_2Lep")
    #b.Scan("RunNumber")

    #Make TChain of all data sets
    #data16PeriodA0lep.Scan("RunNumber")

    dataruns = {}
    #print binMin, binMax
    #print nBins, "Max: ", type(binMin), type(binMax)
    #print nBins, eval(binMin), eval(binMax)
    #print int(binMin)
    # c = TCanvas("c", "")

    runHisto2Lep = TH1F("runHisto2Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto2Lep", weightString+"*("+cutStr+")")

    # c.cd()
    # runHisto2Lep.Draw()
    # c.SaveAs("testDist.pdf")
    #datatree0lep.Scan("RunNumber")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        #print "entry: ", i
        #print str(runHisto2Lep.GetXaxis().(i+binMin))#GetBinLabel(i))
        #print runHisto2Lep.GetBinContent(i)
        dataruns[str(i+binMin)]=runHisto2Lep.GetBinContent(i+1)
    #print "Dat run number: ", dataruns["298967"], dataruns["298966"], dataruns["298968"], dataruns[str(binMin)], lumiDict[str(binMin)]
    #print dataruns

    # for i in xrange(datatree0lep.GetEntries()):
    #     datatree0lep.GetEntry(i)
    #     if str(datatree0lep.RunNumber) in dataruns.keys():
    #         dataruns[str(datatree0lep.RunNumber)]+=1
    #     else:
    #         dataruns[str(datatree0lep.RunNumber)]=1

    lumiNormDict2016={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2016[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2016[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ################################################
    #### Repeat for  2015 ########################
    #Memory leak potential?

    ilumiRootFile = TFile("ilumi2hist2015.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")


    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print "2015 bins: ", newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_2Lep")

    datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")

    dataruns = {}

    #c = TCanvas("c", "")
    runHisto2Lep = TH1F("runHisto2Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto2Lep", weightString+"*("+cutStr+")")
    #c.cd()
    #runHisto2Lep.Draw()
    #c.SaveAs("testDist.pdf")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        dataruns[str(i+binMin)]=runHisto2Lep.GetBinContent(i+1)


    lumiNormDict2015={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2015[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2015[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ##########################################
    ####### Combine ###############

    lumiNormDict = lumiNormDict2015.copy()
    lumiNormDict.update(lumiNormDict2016)
    del lumiNormDict['']
    #print lumiNormDict.keys(), "Min: ", min(lumiNormDict.keys()), "Max: ", max(lumiNormDict.keys())

    newBinMin = int(min(lumiNormDict.keys()))
    newBinMax = int(max(lumiNormDict.keys()))
    newNBins = newBinMax-newBinMin

    newHist = TH1F("newHist", "", newNBins, newBinMin, newBinMax)
    for key, value in lumiNormDict.iteritems():
        newHist.SetBinContent(int(key)-newBinMin, value)

    #for i in xrange(newNBins):
    #    if newHist.GetBinContent(i+1)==0:
    #        newHist.
    c = TCanvas("c","")
    newHist.SetMarkerStyle(20)
    newHist.SetMarkerSize(.5)
    gStyle.SetOptStat(0)
    newHist.Draw("P")
    newHist.GetXaxis().SetTitle("Run Number")
    newHist.GetYaxis().SetTitle("Lumi-normalized yield")
    c.SaveAs("lumiNormPlots/"+cr2key+"LumiNormalized.eps")


######################################################
############  1 Photon CRs  ############################
#####################################

for p1key, p1value in cr1PhoCuts.iteritems():
    weightString = "1."#lepWeightStr
    cutStr = p1value
    print weightString, cutStr
    
    ilumiRootFile = TFile("ilumi2histo.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")

    ilumiRootFile2015 = TFile("ilumi2hist2015.root")
    ilumiFile2015 = ilumiRootFile.Get("lumi_histo")
    print type(ilumiFile2015), type(ilumiFile)

    # testHist = TH1F("testHist", "", 1000000, 0, 1000000)
    # testHist.Add(ilumiFile)
    # testHist.Add(ilumiFile2015)
    # c = TCanvas("c", "")
    # #runHisto0Lep = ilumiFile.Add(ilumiFile2015)
    # #runHisto0Lep.Draw()
    # c.cd()
    # testHist.Draw()
    # c.SaveAs("testRunDist.pdf")

    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_1Lep")
    # data16PeriodAfile = TFile(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    # data16PeriodBfile = TFile(args.inDir+"data16PeriodB_SUSY1_p2667Ext.root")
    # data16PeriodCfile = TFile(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    # data16PeriodDfile = TFile(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    # data16PeriodEfile = TFile(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    # data16PeriodA0lep = data16PeriodAfile.Get("stop_1Lep")
    # data16PeriodB0lep = data16PeriodBfile.Get("stop_1Lep")
    # data16PeriodC0lep = data16PeriodCfile.Get("stop_1Lep")
    # data16PeriodD0lep = data16PeriodDfile.Get("stop_1Lep")
    # data16PeriodE0lep = data16PeriodEfile.Get("stop_1Lep")

    datatree0lep.Add(args.inDir+"data16PeriodA_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodB123_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodC_SUSY1_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data16PeriodD_SUSY1_p2689Ext.root")
    #datatree0lep.Add(args.inDir+"data16PeriodE_SUSY1_p2689Ext.root")

    #a = TFile(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    #b = a.Get("stop_1Lep")
    #b.Scan("RunNumber")

    #Make TChain of all data sets
    #data16PeriodA0lep.Scan("RunNumber")

    dataruns = {}
    #print binMin, binMax
    #print nBins, "Max: ", type(binMin), type(binMax)
    #print nBins, eval(binMin), eval(binMax)
    #print int(binMin)
    # c = TCanvas("c", "")

    runHisto1Lep = TH1F("runHisto1Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto1Lep", weightString+"*("+cutStr+")")

    # c.cd()
    # runHisto1Lep.Draw()
    # c.SaveAs("testDist.pdf")
    #datatree0lep.Scan("RunNumber")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        #print "entry: ", i
        #print str(runHisto1Lep.GetXaxis().(i+binMin))#GetBinLabel(i))
        #print runHisto1Lep.GetBinContent(i)
        dataruns[str(i+binMin)]=runHisto1Lep.GetBinContent(i+1)
    #print "Dat run number: ", dataruns["298967"], dataruns["298966"], dataruns["298968"], dataruns[str(binMin)], lumiDict[str(binMin)]
    #print dataruns

    # for i in xrange(datatree0lep.GetEntries()):
    #     datatree0lep.GetEntry(i)
    #     if str(datatree0lep.RunNumber) in dataruns.keys():
    #         dataruns[str(datatree0lep.RunNumber)]+=1
    #     else:
    #         dataruns[str(datatree0lep.RunNumber)]=1

    lumiNormDict2016={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2016[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2016[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ################################################
    #### Repeat for  2015 ########################
    #Memory leak potential?

    ilumiRootFile = TFile("ilumi2hist2015.root")
    ilumiFile = ilumiRootFile.Get("lumi_histo")


    nBins = ilumiFile.GetNbinsX()
    binMin = int(ilumiFile.GetXaxis().GetBinLabel(1))
    binMax = int(ilumiFile.GetXaxis().GetBinLabel(nBins))
    newBins = binMax-binMin
    print "2015 bins: ", newBins, binMin, binMax

    print "Number of bins: ", nBins

    lumiDict={}
    for i in xrange(nBins):
        #print i
        lumiDict[str(ilumiFile.GetXaxis().GetBinLabel(i))] = ilumiFile.GetBinContent(i)
    #print lumiDict

    datatree0lep = TChain("stop_1Lep")

    datatree0lep.Add(args.inDir+"data15_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")
    datatree0lep.Add(args.inDir+"data15_13TeV.periodJ.physics_Main.PhysCont.DAOD_SUSY1.grp15_v01_p2667Ext.root")

    dataruns = {}

    #c = TCanvas("c", "")
    runHisto1Lep = TH1F("runHisto1Lep", "", newBins, binMin, binMax)
    datatree0lep.Draw("RunNumber>>runHisto1Lep", weightString+"*("+cutStr+")")
    #c.cd()
    #runHisto1Lep.Draw()
    #c.SaveAs("testDist.pdf")

    for i in xrange(0, newBins):#eval(binMin), eval(binMax)):
        dataruns[str(i+binMin)]=runHisto1Lep.GetBinContent(i+1)


    lumiNormDict2015={}
    for i in lumiDict.keys():
        #print dataruns.keys()
        #print i
        if lumiDict[i]!=0:
            lumiNormDict2015[i]=dataruns[i]/lumiDict[i]
        else:
            lumiNormDict2015[i]=0
    #print "Data: ", dataruns
    #print "Lumi: ", lumiDict
    #print lumiNormDict

    ##########################################
    ####### Combine ###############

    lumiNormDict = lumiNormDict2015.copy()
    lumiNormDict.update(lumiNormDict2016)
    del lumiNormDict['']
    #print lumiNormDict.keys(), "Min: ", min(lumiNormDict.keys()), "Max: ", max(lumiNormDict.keys())

    newBinMin = int(min(lumiNormDict.keys()))
    newBinMax = int(max(lumiNormDict.keys()))
    newNBins = newBinMax-newBinMin

    newHist = TH1F("newHist", "", newNBins, newBinMin, newBinMax)
    for key, value in lumiNormDict.iteritems():
        newHist.SetBinContent(int(key)-newBinMin, value)

    #for i in xrange(newNBins):
    #    if newHist.GetBinContent(i+1)==0:
    #        newHist.
    c = TCanvas("c","")
    newHist.SetMarkerStyle(20)
    newHist.SetMarkerSize(.5)
    gStyle.SetOptStat(0)
    newHist.Draw("P")
    newHist.GetXaxis().SetTitle("Run Number")
    newHist.GetYaxis().SetTitle("Lumi-normalized yield")
    c.SaveAs("lumiNormPlots/"+p1key+"LumiNormalized.eps")

"""

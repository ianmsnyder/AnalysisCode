#!/usr/bin/env python

from ROOT import *#TFile, gROOT, TLorentzVector, TTree, gSystem
from optparse import OptionParser
from array import array
from itertools import combinations
import math


def addTopRecoTrees(inFDir, signalFName, inFName=""):
    gROOT.SetBatch(True);

    inFs = {};
    if inFName == "":
        inFs["ttbar"] = TFile.Open(inFDir+"ttbar_SUSY1_p2540.root", "UPDATE");
        inFs["dibosons"] = TFile.Open(inFDir+"dibosons_SUSY1_p2540.root", "UPDATE");
        inFs["singleTop"] = TFile.Open(inFDir+"singleTop_SUSY1_p2540.root", "UPDATE");
        inFs["ttV"] = TFile.Open(inFDir+"ttV_SUSY1_p2540.root", "UPDATE");
        inFs["W"] = TFile.Open(inFDir+"W_SUSY1_p2540.root", "UPDATE");
        inFs["Z"] = TFile.Open(inFDir+"Z_SUSY1_p2540.root", "UPDATE");

        inFs["signal"] = TFile.Open(inFDir+signalFName, "UPDATE");
    else:
        inFs[inFName.replace(".root", "")] = TFile.Open(inFDir+inFName, "UPDATE");


    for samp in inFs.keys():
        wConstTopM = array('f', [-999]);
        dRAntiKt12AntiKt8 = array('f', [999,999]);
        inFs[samp].cd();

        topTree = TTree("topTree", "topTree");
        topTree.Branch("WConstTopM", wConstTopM, "WConstTopM/F");
        topTree.Branch("DRAntiKt12AntiKt8", dRAntiKt12AntiKt8, "DRAntiKt12AntiKt8[2]/F");

        #MT2 stuff:

        top1BestMass = array("f", [0])
        top2BestMass = array("f", [0])
        topBestMassMT2 = array("f", [0])

        mCT = array("f", [0])
        mt2Variance = array("f", [0])

        topTree.Branch('mt2Variance', mt2Variance, 'mt2Variance/F')
        
        bestWMT2 = array("f", [0])

        topTree.Branch("top1BestMass", top1BestMass, "top1BestMass/F")
        topTree.Branch("top2BestMass", top2BestMass, "top2BestMass/F")
        topTree.Branch("topBestMassMT2", topBestMassMT2, "topBestMassMT2/F")
        topTree.Branch("bestWMT2", bestWMT2, "bestWMT2/F")

        topTree.Branch("mCT", mCT, "mCT/F")

        Mt2TopMass0 = array("f", [0])
        Mt2TopMass1= array("f", [0])
        maxMT2 = array("f", [0])

        topTree.Branch('Mt2TopMass0', Mt2TopMass0, 'Mt2TopMass0/F')
        topTree.Branch('Mt2TopMass1', Mt2TopMass1, 'Mt2TopMass1/F')
        topTree.Branch('maxMT2', maxMT2, 'maxMT2/F')
        
        tree = inFs[samp].Get("stop_0Lep");
        tree.SetBranchStatus("*", False);
        tree.SetBranchStatus("NotInDRMinLeadConstJetMassTop", True);
        tree.SetBranchStatus("DRMinLeadConstJetMassTopPt", True);
        tree.SetBranchStatus("Jet*", True);
        tree.SetBranchStatus("NJets*", True);        
        tree.SetBranchStatus("AntiKt*", True);
        tree.SetBranchStatus("Met*", True);
        #tree.SetBranchStatus("MtBMin", True);
        print samp


        
        entries = tree.GetEntries();
        
        for e in xrange(entries):
            if e%1000==0:
                print e
            tree.GetEntry(e);
            #if len(tree.JetPt) < 4 and tree.AntiKt12M[0] > 120 and tree.AntiKt12M[0] < 60 and tree.MtBMin > 175:
            #    continue;
            if len(tree.AntiKt12M) >= 2 and len(tree.AntiKt8M) >=2:
                for fI in range(2):
                    temp12TLV = TLorentzVector();
                    temp8TLV = TLorentzVector();
                    temp12TLV.SetPtEtaPhiM(tree.AntiKt12Pt[fI], tree.AntiKt12Eta[fI], tree.AntiKt12Phi[fI], tree.AntiKt12M[fI])
                    temp8TLV.SetPtEtaPhiM(tree.AntiKt8Pt[fI], tree.AntiKt8Eta[fI], tree.AntiKt8Phi[fI], tree.AntiKt8M[fI]);
                    dRAntiKt12AntiKt8[fI] = temp12TLV.DeltaR(temp8TLV)
                    
            nSoftSideJets = len(tree.NotInDRMinLeadConstJetMassTop);
            nHeavyJets = 0;
            heavyJetIndexes = []
            bestHeavyJetI = -99;
            for iJet in tree.NotInDRMinLeadConstJetMassTop:
                if tree.JetM[iJet] > 60:
                    nHeavyJets += 1;
                    bestHeavyJetI = iJet;
                    heavyJetIndexes.append(iJet)

            bestDM = 999999;
            if nHeavyJets >= 2:
                for iJet in heavyJetIndexes:
                    if abs(80-tree.JetM[iJet]) < bestDM:
                        bestDM = abs(80-tree.JetM[iJet]);
                        bestHeavyJetI = iJet;
            minDR = 999999;
            minDRI = -99;
            topTLV = TLorentzVector();
            if nHeavyJets > 0 and False:
                topTLV.SetPtEtaPhiM(tree.JetPt[bestHeavyJetI], tree.JetEta[bestHeavyJetI], tree.JetPhi[bestHeavyJetI], tree.JetM[bestHeavyJetI])
                for iJet in tree.NotInDRMinLeadConstJetMassTop:
                    if iJet == bestHeavyJetI:
                        continue;
                    jetTLV = TLorentzVector();
                    jetTLV.SetPtEtaPhiM(tree.JetPt[iJet], tree.JetEta[iJet], tree.JetPhi[iJet], tree.JetM[iJet])
                    dR = jetTLV.DeltaR(topTLV);
                    if dR < minDR:
                        minDR = dR;
                        minDRI = iJet;
                if minDRI > 0:
                    topDaughtTLV = TLorentzVector();
                    topDaughtTLV.SetPtEtaPhiM(tree.JetPt[minDRI], tree.JetEta[minDRI], tree.JetPhi[minDRI], tree.JetM[minDRI])

                    topTLV += topDaughtTLV;

            else:
                if nSoftSideJets <=3:
                    for iJet in tree.NotInDRMinLeadConstJetMassTop:
                        jetTLV = TLorentzVector();
                        jetTLV.SetPtEtaPhiM(tree.JetPt[iJet], tree.JetEta[iJet], tree.JetPhi[iJet], tree.JetM[iJet]);
                        topTLV+=jetTLV;
                else:
                    combos = list(combinations(tree.NotInDRMinLeadConstJetMassTop, 2))
                    bestWCombo = [];
                    bestDW = 999999;
                    wMass = 80;
                    for combo in combos:
                        tempTLV = TLorentzVector();
                        for iJet in combo:
                            jetTLV = TLorentzVector();
                            jetTLV.SetPtEtaPhiM(tree.JetPt[iJet], tree.JetEta[iJet], tree.JetPhi[iJet], tree.JetM[iJet]);
                            tempTLV+=jetTLV;
                        #if tempTLV.Pt() < tree.DRMinLeadConstJetMassTopPt:
                        #    combosWithLowerPt.append(combo)
                        if bestDM > abs(tempTLV.M()-wMass):
                            bestDM = abs(tempTLV.M()-wMass)
                            bestWCombo =  combo
                    wTLV = TLorentzVector();
                    wTLV.SetPtEtaPhiM(tree.JetPt[combo[0]],tree.JetEta[combo[0]],tree.JetPhi[combo[0]],tree.JetM[combo[0]])
                    tempTLV = TLorentzVector();
                    tempTLV.SetPtEtaPhiM(tree.JetPt[combo[1]],tree.JetEta[combo[1]],tree.JetPhi[combo[1]],tree.JetM[combo[1]])
                    wTLV += tempTLV;
                    bestTopDM = 999999;
                    topMass = 173.5;
                    bestTagWeight = -999999;
                    for iJet in tree.NotInDRMinLeadConstJetMassTop:
                        tempTopTLV = TLorentzVector()
                        tempTopTLV+=wTLV;
                        jetTLV = TLorentzVector();
                        jetTLV.SetPtEtaPhiM(tree.JetPt[iJet], tree.JetEta[iJet], tree.JetPhi[iJet], tree.JetM[iJet]);
                        tempTopTLV+=jetTLV
                        if tree.JetMV1Weight[iJet] > bestTagWeight:
                            bestTagWeight =  tree.JetMV1Weight[iJet]
                            topTLV = tempTopTLV
                        # if tree.JetIsBTagged[iJet]:
                        #     tempTopTLV = TLorentzVector()
                        #     tempTopTLV+=wTLV;
                        #     jetTLV = TLorentzVector();
                        #     jetTLV.SetPtEtaPhiM(tree.JetPt[iJet], tree.JetEta[iJet], tree.JetPhi[iJet], tree.JetM[iJet]);
                        #     tempTopTLV+=jetTLV
                        #     if bestTopDM > abs(tempTopTLV.M()-topMass):
                        #         bestTopDM = abs(tempTopTLV.M()-topMass);
                        #         topTLV = tempTopTLV;

                        #print bestDM, bestTopDM, bestTopTLV.M()
                        
            #gSystem.Load("../Stop0LRun2/MT2/MT2_Package/src/libBinnedLik.so")
            
            gSystem.Load("/export/home/snyder/Mt2_2015/MT2_Package/src/libBinnedLik.so")
            MetVector = TLorentzVector()
            MetVector.SetPtEtaPhiE(tree.Met, 0, tree.MetPhi, tree.Met)
            
            if topTLV.M()>0:
                bestWtop1 = TLorentzVector()
                bestWtop1.SetPtEtaPhiM(tree.AntiKt12Pt[0], tree.AntiKt12Eta[0], tree.AntiKt12Phi[0], tree.AntiKt12M[0])
                bestWMT2[0] = ComputeMT2(bestWtop1, topTLV, MetVector, 0., 0.).Compute()
            else:
                bestWMT2[0] = -99

            firstTop = TLorentzVector()
            secondTop = TLorentzVector()
            if len(tree.AntiKt12Pt)>1:
                firstTop.SetPtEtaPhiM(tree.AntiKt12Pt[0], tree.AntiKt12Eta[0], tree.AntiKt12Phi[0], tree.AntiKt12M[0])
                secondTop.SetPtEtaPhiM(tree.AntiKt12Pt[1], tree.AntiKt12Eta[1], tree.AntiKt12Phi[1], tree.AntiKt12M[1])
            mctTop = math.sqrt((firstTop.Et()+secondTop.Et())**2 - (firstTop.Et()-secondTop.Et())**2)
            mCT[0] = mctTop
            # Three best jets to make top
            topMass = 173
            FourVecDict = {}
            MassDict = {}
            MassDiffDict = {}

            for i in xrange(tree.NJets):
                FourVecDict[str(i)] = TLorentzVector()
                FourVecDict[str(i)].SetPtEtaPhiM(tree.JetPt[i], tree.JetEta[i], tree.JetPhi[i], tree.JetM[i])

            # for i in xrange(tree.NJets):
            #     for j in xrange(tree.NJets):
            #         for k in xrange(tree.NJets):
            #             if i!=j and i!=k and j!=k:
            #                     NewVector = TLorentzVector()
            #                     NewVector = FourVecDict[str(i)]+FourVecDict[str(j)]+FourVecDict[str(k)]
            #                     MassDict[str(i)+str(j)+str(k)] = NewVector.M()
            #                     MassDiffDict[str(i)+str(j)+str(k)] = abs(topMass-NewVector.M())
                                

            # Varying MT2
            topMt20 = -99.
            topMt21 = -99.
            varyMt2 = -99.
            mt2List=[]
            if tree.NJets>5:
                jetList = ""
                for i in xrange(tree.NJets):
                    jetList+=str(i)

                #onecombinations = list(combinations(jetList, 1))
                #twocombinations = list(combinations(jetList, 2))
                threecombinations = list(combinations(jetList, 3))

                combinationDict={}
                for i in xrange(len(threecombinations)):
                    combinationDict[str(i)] = TLorentzVector()
                    for j in xrange(len(threecombinations[i])):
                         combinationDict[str(i)]+=FourVecDict[threecombinations[i][j]]
                combinationList = ""
                for i in xrange(len(combinationDict)):
                    combinationList+=str(i)
                allcombinations = list(combinations(combinationList, 2))
                for i in xrange(len(allcombinations)):
                    testMt2 = ComputeMT2(combinationDict[allcombinations[i][0]], combinationDict[allcombinations[i][1]], MetVector, 0., 0.).Compute()
                    NewVector = TLorentzVector()
                    NewVector = combinationDict[allcombinations[i][0]]
                    MassDict[allcombinations[i][0]] = NewVector.M()
                    MassDiffDict[allcombinations[i][0]] = abs(topMass-NewVector.M())
  
                    mt2List.append(testMt2)
                    if testMt2>varyMt2:
                        varyMt2=testMt2
                        topMt20=combinationDict[allcombinations[i][0]].M()
                        topMt21=combinationDict[allcombinations[i][1]].M()
            maxMT2[0] = varyMt2
            if topMt20>topMt21:
                Mt2TopMass0[0]=topMt20
                Mt2TopMass1[0]=topMt21
            else:
                Mt2TopMass1[0]=topMt20
                Mt2TopMass0[0]=topMt21
            if len(mt2List)>2:
                mt2Variance[0] = (max(mt2List)-min(mt2List))/math.sqrt(tree.NJets)
            else: mt2Variance[0]=-99
            wConstTopM[0] = topTLV.M();

            if len(MassDiffDict)>2:
                top1 = min(MassDiffDict, key=MassDiffDict.get)
                top1MassDiff = MassDiffDict[top1]
                del MassDiffDict[top1]
                top2 = min(MassDiffDict, key=MassDiffDict.get)
                top2MassDiff = MassDiffDict[top2]
            
                top1bestmass = combinationDict[str(top1)]#FourVecDict[top1[0]]+FourVecDict[top1[1]]+FourVecDict[top1[2]]
                top2bestmass = combinationDict[str(top2)]#FourVecDict[top2[0]]+FourVecDict[top2[1]]+FourVecDict[top2[2]]


            
                topBestMassMT2[0] = ComputeMT2(top1bestmass, top2bestmass, MetVector, 0., 0.).Compute()

                top1BestMass[0] = top1bestmass.M()
                top2BestMass[0] = top2bestmass.M()
            else:
                topBestMassMT2[0] = -99.
                top1BestMass[0] = -99.
                top2BestMass[0] = -99.

            
            topTree.Fill()
            
#Adding MT2 stuff for this 




            
        topTree.AddFriend(tree.GetName())
        topTree.AddFriend(tree.GetName()+"Ext")
        topTree.Write()
        inFs[samp].Close();

if __name__ == "__main__":
    
    parser=OptionParser();
    parser.add_option("--inFDir",type=str,help="path to input directory or file");
    parser.add_option("--signalFName",type=str,default="mc15_13TeV.387198.MadGraphPythia8EvtGen_A14NNPDF23LO_TT_directTT_800_1.merge.DAOD_SUSY1.e3969_a766_a777_r6282_p2540.root",help="Signal file name");
    parser.add_option("--inFName",type=str,default="",help="Input file name");
    (options,args)=parser.parse_args();
    
    addTopRecoTrees(options.inFDir, options.signalFName, inFName=options.inFName)

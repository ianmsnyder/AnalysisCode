#ifndef Stop0LRun2_Stop0LRun2_H
#define Stop0LRun2_Stop0LRun2_H

#define INIT_VAL -99

// Needed for uint32_t
#include <stdint.h>
#include <EventLoop/Algorithm.h>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "Stop0LRun2/Make_Unique.h"
#include "Stop0LRun2/TwoLeptonProcessor.h"
#include "Stop0LRun2/Ntup4Vec.h"
#include "Stop0LRun2/NtupJetExtraVec.h"
#include "xAODtruthUtils/xAODtruthUtils.h"
//Truth classifying for MET/HT slice
#include "MCTruthClassifier/MCTruthClassifierDefs.h"
#include "MCTruthClassifier/MCTruthClassifier.h"
#ifdef DORJIGSAW
#include "Stop0LRun2/RJigsaw.h"
#endif
#include "JetSubStructureUtils/KtSplittingScale.h"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "Stop0LRun2/NtupSubstructure.h"

// Infrastructure include(s):
#ifdef ROOTCORE
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/TEvent.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#endif

// Histogram related
#include <TH1.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

#include <SUSYTools/SUSYCrossSection.h>
//ROOT override the keyword nullptr and screw it up. Fix that
#undef nullptr

// Other includes
#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"

//#include "MT2_ROOT.h"


#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<vector<int> > >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#pragma link C++ class std::vector < std::vector<double> >+;
#pragma link C++ class std::vector < std::vector<string> >+; 
#endif


// Forward declarations
class GoodRunsListSelectionTool;
class JetCleaningTool;
class JERTool;

namespace ST {
	class SUSYObjDef_xAOD;
}


class Stop0LRun2 : public EL::Algorithm
{

	// put your configuration variables here as public variables.
	// that way they can be set directly from CINT and python.
#ifndef __CINT__
	SUSY::CrossSectionDB *m_xSecDB;  //!
	GoodRunsListSelectionTool *m_grl; //!
	ST::SUSYObjDef_xAOD *m_objTool; //!

#endif // not __CINT__

 public:
	// float cutValue;
	std::string m_outputTreeName;
	std::string m_outputLabel;
	bool m_isData;
	bool m_isAtlFast;
	std::string m_projectTag;
	bool m_isDerived;

	std::string m_configFile;

	bool m_debug;
	bool m_debugTree;

	int m_jesNPSet;

	bool m_syst;
	bool m_doSpecialOR;
	bool m_isMC15b;

	unsigned int m_zeroLepSkimLevel;
	unsigned int m_oneLepSkimLevel;
	unsigned int m_twoLepSkimLevel;

	float m_elPtMin;
	float m_elEtaMax;

	float m_muPtMin;
	float m_muEtaMax;

	float m_jetPtMin;
	float m_jetEtaMax;
	float m_jetSignalPtMin;
	float m_jetSignalEtaMax;
	float m_jetJVTMin;
	float m_mv1Cut;
	float m_ip3dsv1Cut;
	float m_maxBJetEta;
	float m_btagCut;
	float m_truthPtMin;

	double m_dRejet; 
	double m_dRmujet; 
	double m_dRemu; 
	
	int   m_srOffset;

	int   m_finalState;

	bool  m_printEvts;
	int   m_evtToPrint;

	uint64_t m_nEventsProcessed;
	double   m_sumOfWeights;
	double   m_sumOfWeightsSquared;
	int      m_regionSelec;

	std::vector<CP::SystematicSet> m_sysList;

	struct XAODObjects
	{
		xAOD::ElectronContainer* electrons;
		xAOD::ShallowAuxContainer* electrons_aux;
		xAOD::MuonContainer* muons;
		xAOD::ShallowAuxContainer* muons_aux;
		xAOD::JetContainer* jets;
		xAOD::ShallowAuxContainer* jets_aux;

		const xAOD::MissingETContainer* metContainer = 0;
		const xAOD::MissingETContainer* mettrkcontainer = 0;
	};

	// variables that don't get filled at submission time should be
	// protected from being send from the submission node to the worker
	// node (done by the //!)

 public:
	xAOD::TEvent *m_event;  //!
	xAOD::TStore *m_store;  //!

	// this is a standard constructor
	Stop0LRun2();

	// these are the functions inherited from Algorithm
	virtual EL::StatusCode setupJob(EL::Job& job);
	virtual EL::StatusCode fileExecute();
	virtual EL::StatusCode histInitialize();
	virtual EL::StatusCode changeInput(bool firstFile);
	virtual EL::StatusCode initialize();
	virtual EL::StatusCode execute();
	virtual EL::StatusCode postExecute();
	virtual EL::StatusCode finalize();
	virtual EL::StatusCode histFinalize();

	// user defined functions
	void         initTruthVars();
	void         initVars();
	XAODObjects  createObj(std::string systName);
	bool         getTruth();
	void         FillTruthHt(const xAOD::TruthParticleContainer* truthParts);
	void         FillTruthNeutrinoMET(const xAOD::TruthParticleContainer* truthParts);
    /* void         calculateMt2Truth(float topPtTruth, float topEtaTruth, float topPhiTruth, float topETruth, float antitopPtTruth, float antitopEtaTruth, float antitopPhiTruth, float antitopETruth, float metTruth, float metXTruth, float metYTruth); */
    /* void         calculateMt2DRmin(float topPtDRmin, float topEtaDRmin, float topPhiDRmin, float topMDRmin, float antitopPtDRmin, float antitopEtaDRmin, float antitopPhiDRmin, float antitopMDRmin, float metDRmin, float metXDRmin, float metYDRmin); */
    /* void         calculateMt2RCjet(float topPtRCjet, float topEtaRCjet, float topPhiRCjet, float topMRCjet, float antitopPtRCjet, float antitopEtaRCjet, float antitopPhiRCjet, float antitopMRCjet, float metRCjet, float metXRCjet, float metYRCjet); */

	void         FillTruthWb(const xAOD::TruthParticle * particle);
	bool         isFromTopW(const xAOD::TruthParticle* particle);
	int isInChan();
	bool         performSkim(int& channel, unsigned int& systI, const XAODObjects& xAODObjects);
	void         fillCutflow(unsigned int& rI, unsigned int& systI);
	void         setVars();
	void         cleanUp();

	static bool sortBySecond(std::pair<int, float> pair1, std::pair<int, float> pair2) { return pair1.second < pair2.second; };

	std::vector<std::pair<int, float> > getOrderedJetIndexDR(TLorentzVector& refVector, std::vector<int> skipIndices);

	void getProjectTag(std::string sampleName, std::string &projectTag, bool &isData);
	void getTags(std::string sampleName, std::vector<int> &aTags, std::vector<int> &sTags, std::vector<int> &rTags, std::vector<int> &eTags, std::vector<int> &pTags);
	
	// this is needed to distribute the algorithm to the workers
	ClassDef(Stop0LRun2, 1);
 protected:
	void GetExtraJetCollections();
#if DOSUBSTRUCTURE
	void DoTrackSubstructures(XAODObjects& xAODObjects);
	void DoSubstructures();
#endif // DOSUBSTRUCTURE
#if DOPFLOW
	void InitPFlowSUSYTools(ST::SettingDataSource dataSource);
#endif // DOPFLOW
	
 private:
	// Arrays of various objects
	TObjArray *m_electrons; //!
	TObjArray *m_muons; //!
	TObjArray *m_jets; //!

	TFile* m_outputFile; //!
	// Output tree and branch variable definitions
	std::vector<std::vector< TTree* > >     m_outTrees;     //!
	std::vector<std::vector< TH1F* > >      m_cutflow;      //!
	std::vector<std::vector< TH1D* > >      m_cutflowWeighted;      //!
	//Don't use smart pointer. ROOT system own the histograms
	TH1I* m_runNumberNEvents; //!

	std::vector<std::vector<ST::SystInfo> > m_systInfoLists; //!
	std::vector<ST::SystInfo> m_weightSystInfo; //!
	std::map<std::string, double> m_weights; //!

	bool m_isNominal; //!

	bool m_affectsKinematics; //!
	bool m_affectsWeights; //!

	// If necessary (kinematics affected), make a shallow copy with the variation applied
	bool m_syst_affectsElectrons; //!
	bool m_syst_affectsMuons; //!
	bool m_syst_affectsTaus; //!
	bool m_syst_affectsPhotons; //!
	bool m_syst_affectsJets; //!
	bool m_syst_affectsBTag; //!
	// bool m_syst_affectsMET; //!


	// Get the nominal object containers from the event
	// Electrons
	xAOD::ElectronContainer* m_electrons_nominal; //!
	xAOD::ShallowAuxContainer* m_electrons_nominal_aux; //!

	// // Photons
	// xAOD::PhotonContainer* m_photons_nominal; //!
	// xAOD::ShallowAuxContainer* m_photons_nominal_aux; //! 

	// Muons
	xAOD::MuonContainer* m_muons_nominal; //!
	xAOD::ShallowAuxContainer* m_muons_nominal_aux; //!

	// Jets
	xAOD::JetContainer* m_jets_nominal; //!
	xAOD::ShallowAuxContainer* m_jets_nominal_aux; //!

	// // Taus
	// xAOD::TauJetContainer* m_taus_nominal; //!
	// xAOD::ShallowAuxContainer* m_taus_nominal_aux; //!

	// MET
	xAOD::MissingETContainer* m_metcst_nominal; //!
	xAOD::MissingETAuxContainer* m_metcst_nominal_aux; //!

	xAOD::MissingETContainer* m_mettst_nominal; //!
	xAOD::MissingETAuxContainer* m_mettst_nominal_aux; //!

	xAOD::MissingETContainer* m_mettrack_nominal; //!
	xAOD::MissingETAuxContainer* m_mettrack_nominal_aux; //!

	xAOD::MissingETContainer* m_metleptst_nominal; //!
	xAOD::MissingETAuxContainer* m_metleptst_nominal_aux; //!

	TH1I* nEventsProcessedHist;	//!
	TH1D* sumOfWeightsHist;	//!
	TH1D* sumOfWeightsSquaredHist;	//!

	int    m_eventNumber;         //!
	int    m_runNumber;           //!
	int    m_lumiBlock;           //!
	int    m_mcChannelNumber;     //!
	int    m_isZll;               //!
	double m_mcEventWeight;       //!
	double m_xSecWeight;          //!
	double m_pileupWeight;        //!
	double m_elecSF;        //!
	double m_muonSF;        //!
	double m_bTagSF;        //!
	double m_elecSF_nominal;        //!
	double m_muonSF_nominal;        //!
	double m_bTagSF_nominal;        //!
	float  m_avgIntPerXing;       //!

	std::vector<int>       m_aTags;  //!
	std::vector<int>       m_sTags;  //!
	std::vector<int>       m_eTags;  //!
	std::vector<int>       m_rTags;  //!
	std::vector<int>       m_pTags;  //!

	int m_nElectrons; //!
	int m_nSigElectrons; //!
	std::vector<float>     m_elPt;        //!
	std::vector<float>     m_elE; //!
	std::vector<float>     m_elEta;       //!
	std::vector<float>     m_elPhi;       //!
	std::vector<int>       m_elIsSignal;  //!

	int m_nMuons; //!
	int m_nSigMuons; //!
	int m_nCosmicMuons; //!
	int m_nBadMuons; //!
	std::vector<float>     m_muPt;        //!
	std::vector<float>     m_muE; //!
	std::vector<float>     m_muEta;       //!
	std::vector<float>     m_muPhi;       //!
	std::vector<int>       m_muIsSignal;  //!

	int                   m_nJets;        //!
	int                   m_nBadJets;     //!
	int                   m_nBJets;       //!
	std::vector<float>    m_jetPt;        //!
	std::vector<float>    m_jetEMFrac;        //!
	std::vector<float>    m_jetHECFrac;        //!
	std::vector<float>    m_jetE; //!
	std::vector<float>    m_jetEta;       //!
	std::vector<float>    m_jetPhi;       //!
	std::vector<float>    m_jetM; //!
	std::vector<float>    m_jetIsSignal;  //!
	std::vector<float>    m_jetMV1Weight; //!
	std::vector<int>      m_jetIsBTagged; //!
	std::vector<float>    m_jetMt;        //!
	std::vector<float>    m_jetDPhiMet;   //!
	float                 m_jetDPhiMetMin3;       //!
	int                   m_jetLeadTagIndex;      //!
	int                   m_jetSubleadTagIndex;   //!

	float                 m_topDRMin0Pt;  //!
	float                 m_topDRMin0Phi; //!
	float                 m_topDRMin0Eta; //!
	float                 m_topDRMin0M;   //!

	float                 m_topDRMin1Pt;  //!
	float                 m_topDRMin1Phi; //!
	float                 m_topDRMin1Eta; //!
	float                 m_topDRMin1M;   //!

	float                 m_topDRMinB0Pt; //!
	float                 m_topDRMinB0Phi;        //!
	float                 m_topDRMinB0Eta;        //!
	float                 m_topDRMinB0M;  //!
	std::vector<int>      m_topDRMinB0Const;      //!

	float                 m_topDRMinB1Pt; //!
	float                 m_topDRMinB1Phi;        //!
	float                 m_topDRMinB1Eta;        //!
	float                 m_topDRMinB1M;  //!
	std::vector<int>      m_topDRMinB1Const;      //!

	int                   m_nAntiKt8;     //!
	float                 m_antiKt8DijetMassAsym; //!
	std::vector<float>    m_antiKt8Pt;    //!
	std::vector<float>    m_antiKt8Phi;   //!
	std::vector<float>    m_antiKt8Eta;   //!
	std::vector<float>    m_antiKt8M;     //!
	std::vector<int>      m_antiKt8NConst;        //!
	std::vector<std::vector<int> > m_antiKt8Const;        //!

	int                   m_nAntiKt12;    //!
	float                 m_antiKt12DijetMassAsym;        //!
	std::vector<float>    m_antiKt12Pt;   //!
	std::vector<float>    m_antiKt12Phi;  //!
	std::vector<float>    m_antiKt12Eta;  //!
	std::vector<float>    m_antiKt12M;    //!
	std::vector<int>      m_antiKt12NConst;       //!
	std::vector<std::vector<int> > m_antiKt12Const;       //!

	float m_metX; //!
	float m_metY; //!
	float m_metPhi; //!
	float m_met; //!
	float m_metSumEt; //!

	float m_metXOrig; //!
	float m_metYOrig; //!
	float m_metPhiOrig; //!
	float m_metOrig; //!
	float m_metSumEtOrig; //!

	float m_metXTruth; //!
	float m_metYTruth; //!
	float m_metPhiTruth; //!
	float m_metTruth; //!
	float m_metSumEtTruth; //!

    float metXparticle; //!
    float metYparticle; //!
    float metparticle; //!
    
	float m_metXTrk; //!
	float m_metYTrk; //!
	float m_metPhiTrk; //!
	float m_metTrk; //!
	float m_metSumEtTrk; //!

	float m_rawMetXTrk; //!
	float m_rawMetYTrk; //!
	float m_rawMetPhiTrk; //!
	float m_rawMetTrk; //!

	float m_metLepX; //!
	float m_metLepY; //!
	float m_metLepPhi; //!
	float m_metLep; //!
	float m_metLepSumEt; //!
	
	float m_ht; //!
	float m_htSig; //!
	float m_metSig; //!
	float m_mEff; //!

	int m_primVtxNTrks; //!
	int m_nvtx; //!
	float m_mtBMin; //!
	float m_mtNonBMin; //!
	float m_minMt; //!
	float m_mtTauCand; //!
	int m_tauJetNTracks; //!

	bool m_hasPrimVtx; //!

	bool m_passMetTrigger; //!
	bool m_passLepTrigger; //!
	bool m_passPhotonTrigger; //!

	float m_metTriggerPrescale; //!
	float m_lepTriggerPrescale; //!
	float m_photonTriggerPrescale; //!

	bool m_passTrigMatch; //!
	bool m_passOSLep; //!

	std::vector<float> m_lepPt; //!
	std::vector<float> m_lepE; //!
	std::vector<float> m_lepPhi; //!
	std::vector<float> m_lepEta; //!
	std::vector<int> m_lepType; //!

	float m_mtMetLep; //!
	float m_dPhiMetLep; //!
	
	float m_llMass; //!
 
	int m_susyPdgId1; //!
	int m_susyPdgId2; //!

	std::vector<float>    m_mcPt; //!
	std::vector<float>    m_mcPhi;        //!
	std::vector<float>    m_mcEta;        //!
	std::vector<float>    m_mcE;  //!
	std::vector<int>      m_mcPdgId;      //!
	std::vector<int>      m_mcBarcode;    //!
	std::vector<std::vector<int> > m_mcChildren;  //!

	int m_nTruthJets; //!
	std::vector<float> m_truthJetPt; //!
	std::vector<float> m_truthJetPhi; //!
	std::vector<float> m_truthJetEta; //!
	std::vector<float> m_truthJetM; //!

	float m_truthMetNonInt; //!

	float m_dRbbJetHt; //!
	float m_dRbbJet; //!

	float m_mBB; //!
	float m_dRbb; //!

	float m_rawDPhiMetTrackMet; //!
	float m_dPhiMetTrackMet; //!

	float m_allJetMass; //!
	float m_allJetPt; //!
	float m_allJetPhi; //!
	float m_allJetEta; //!

	float m_trueTopPt; //!
	float m_trueTopEta; //!
	float m_trueTopPhi; //!
	float m_trueTopE; //!
	std::vector<int> m_trueTopDaughtPdgId; //!
	std::vector<float> m_trueTopDaughtPt; //!
	std::vector<float> m_trueTopDaughtEta; //!
	std::vector<float> m_trueTopDaughtPhi; //!
	std::vector<float> m_trueTopDaughtE; //!

	float m_trueAntitopPt; //!
	float m_trueAntitopEta; //!
	float m_trueAntitopPhi; //!
	float m_trueAntitopE; //!
	std::vector<int> m_trueAntitopDaughtPdgId; //!
	std::vector<float> m_trueAntitopDaughtPt; //!
	std::vector<float> m_trueAntitopDaughtEta; //!
	std::vector<float> m_trueAntitopDaughtPhi; //!
	std::vector<float> m_trueAntitopDaughtE; //!

	float truthHt;	//!
	float truthNeutrinoMet;	//!

	int m_topDecayType; //!
	int m_antitopDecayType; //!

    float trueMet; //!
    float trueXMet; //!
    float trueYMet; //!
   

    //mt2 stuff:
    float mt2_truth; //!
    float mt2_rcjet; //!
    float mt2_drmin; //!

	const std::string nominalJetCollection = "AntiKt4EMTopoJets";
	const std::vector<std::string> extraJetCollections ={"AntiKt2EMTopoJets", "AntiKt3EMTopoJets"};
	//Need to delay construction, so use (smart) pointer
	std::vector<std::vector<std::unique_ptr<TwoLeptonProcessor> > > twoLeptonProcessors; //!
	std::unique_ptr<Ntup4Vec> truthNeutrinos;	//!
	std::unique_ptr<Ntup4Vec> truthWs;	//!
	std::unique_ptr<Ntup4Vec> truthBottoms;	//!
    std::unique_ptr<Ntup4Vec> trueMetParticle;	//!

	MCTruthClassifier* m_mcTruthClassifier;	//!
#if DORJIGSAW
	std::unique_ptr<RJigsaw> rJigsaw;	//!
#endif
	std::unordered_map<std::string, std::unique_ptr<JetCalibrationTool> > rScanCalibTools;	//!
	std::unique_ptr<JetVertexTaggerTool> jvtTool;	//!
	std::unique_ptr<JetCleaningTool> jetCleaningTool;	//!

	//Extra jet collections
	std::unique_ptr<std::unordered_map<std::string, std::unique_ptr<Ntup4Vec> > > extraJetNtup4VecsMap;	//!
	//Extra jet variables
	std::unique_ptr<std::unordered_map<std::string, std::unique_ptr<NtupJetExtraVec> > > ntupJetExtraMap;	//!

	std::unique_ptr<InDet::InDetTrackSelectionTool> trackSelector;	//!
	std::unique_ptr<JetSubStructureUtils::KtSplittingScale> d12Calculator;	//!
	std::unique_ptr<JetSubStructureUtils::KtSplittingScale> d23Calculator;	//!
	std::unique_ptr<fastjet::contrib::Nsubjettiness> tau1Calculator; //!
	std::unique_ptr<fastjet::contrib::Nsubjettiness> tau2Calculator; //!
	std::unique_ptr<fastjet::contrib::Nsubjettiness> tau3Calculator; //!
	std::unique_ptr<fastjet::contrib::EnergyCorrelator> ECF1Calculator;	//!
	std::unique_ptr<fastjet::contrib::EnergyCorrelator> ECF2Calculator;	//!
	std::unique_ptr<fastjet::contrib::EnergyCorrelator> ECF3Calculator;	//!
	std::unique_ptr<NtupSubstructure> akt4EventSubstructure;	//!
	std::unique_ptr<ST::SUSYObjDef_xAOD> pflowSUSYTools;	//!

	//Don't own these. Don't use smart pointer
	const ToolHandle<CP::IPileupReweightingTool>* prwTool;	//!
};

#endif

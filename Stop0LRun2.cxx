
#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/ToolsMeta.h"
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoop/ProofDriver.h"
#include "Stop0LRun2/Stop0LRun2.h"
#include <fastjet/CompositeJetStructure.hh>
#include "Stop0LRun2/JetFlavourInfo.h"

// Stop0LUtils
#include "Stop0LRun2/Jet.h"
#include "Stop0LRun2/Electron.h"
#include "Stop0LRun2/Muon.h"
#include "Stop0LRun2/Tau.h"
#include "Stop0LRun2/MtCalculator.h"
#include "Stop0LRun2/TopReco.h"
#include "Stop0LRun2/Stop1LepFunctionUtil.h"
#include "Stop0LRun2/MT2_ROOT.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"
#include "xAODBTagging/BTagging.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODCore/ShallowCopy.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "PileupReweighting/TPileupReweighting.h"
// Needed if writing an output xAOD
#include "EventLoop/OutputStream.h"
// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//Use __func__ macro instead of APP_NAME
#define APP_NAME __func__
//#include "CPAnalysisExamples/errorcheck.h"
#include "Stop0LRun2/ErrorCheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#include "TruthUtils/PIDUtils.h"

#include <boost/regex.hpp>

// this is needed to distribute the algorithm to the workers
ClassImp(Stop0LRun2)



Stop0LRun2::Stop0LRun2()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
	m_outputLabel      = "output.01";
	m_outputTreeName   = "OutputTree";
	m_debug            = false;
	m_debugTree        = false;
	m_isData           = false;
	m_isAtlFast        = true;
	m_projectTag       = "mc15_13TeV";
	m_jesNPSet         = 0;
	m_syst             = false;
	m_doSpecialOR      = false;
	m_isMC15b          = false;

	m_maxBJetEta       = 2.5;     // B-jet max eta
	m_btagCut          = -0.4434;	// 77% WP, default of SUSYTools
	m_dRejet           = 0.2;     // Delta R for electron-jet overlap
	m_dRmujet          = 0.2;     // Delta R for muon-jet overlap
	m_dRemu            = 0.01;    // Delta R for muon-electron overlap
	
	m_zeroLepSkimLevel = 0;       // What requirements to make on the output, 0 is output all events
	m_oneLepSkimLevel  = 0;       // What requirements to make on the output, 0 is output all events
	m_twoLepSkimLevel  = 0;       // What requirements to make on the output, 0 is output all events
	m_srOffset         = 5;       // Cutflow offset for SRB
	m_finalState       = 0;       // Final state used for xSec, mainly important for signal
	m_truthPtMin       = 10.;     // Min pT requirement for truth particles
	m_printEvts        = false;   // Print event number flag. Useful for cutflow comparison.
	m_evtToPrint       = -1;      // Print all kinematic info for this event.
	m_regionSelec      = -1;
	
	m_objTool          = 0;
	m_grl              = 0;
	m_xSecDB           = 0;
	m_configFile       = "SUSYTools_Default.conf";
}



EL::StatusCode Stop0LRun2::setupJob(EL::Job& job)
{
	// Here you put code that sets up the job on the submission object
	// so that it is ready to work with your algorithm, e.g. you can
	// request the D3PDReader service or add output files.  Any code you
	// put here could instead also go into the submission script.  The
	// sole advantage of putting it here is that it gets automatically
	// activated/deactivated when you add/remove the algorithm from your
	// job, which may or may not be of value to you.
	job.useXAOD();

	// let's initialize the algorithm to use the xAODRootAccess package
	xAOD::Init("Stop0LRun2").ignore(); // call before opening first file

	// tell EventLoop about our output file:
	EL::OutputStream out(m_outputLabel);
	job.outputAdd(out);
	//Ntuple
	EL::NTupleSvc *ntuple = new EL::NTupleSvc(m_outputLabel);
	job.algsAdd(ntuple);

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::histInitialize()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.

	nEventsProcessedHist = new TH1I("NEventsHist", "NEventsHist", 1, 0, 1);
	sumOfWeightsHist = new TH1D("SumOfWeightsHist", "SumOfWeightsHist", 1, 0, 1);
	sumOfWeightsSquaredHist = new TH1D("SumOfWeightsSquaredHist", "SumOfWeightsSquaredHist", 1, 0, 1);

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::fileExecute()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::changeInput(bool firstFile)
{
	// Here you do everything you need to do when we change input files,
	// e.g. resetting branch addresses on trees.  If you are using
	// D3PDReader or a similar service this method is not needed.

	auto tevent = wk()->xaodEvent();

	std::string sampleName = wk()->metaData()->castString("nc_source");
	getProjectTag(sampleName, m_projectTag, m_isData);

	if (!m_isData) {
		//Read the CutBookkeeper container
		const xAOD::CutBookkeeperContainer* completeCBC = 0;
		if (!tevent->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()) {
			Error(APP_NAME, "Failed to retrieve CutBookkeepers from MetaData! Exiting.");
		}

		// First, let's find the smallest cycle number,
		// i.e., the original first processing step/cycle
		int minCycle = 10000;
		for (auto cbk : *completeCBC) {
			if (minCycle > cbk->cycle()) { minCycle = cbk->cycle(); }
		}
		// Now, let's actually find the right one that contains all the needed info...
		const xAOD::CutBookkeeper* allEventsCBK = 0;
		for (auto cbk : *completeCBC) {
			if (minCycle == cbk->cycle() && cbk->name() == "AllExecutedEvents") {
				allEventsCBK = cbk;
				break;
			}
		}
		if (allEventsCBK) {
			m_nEventsProcessed  = allEventsCBK->nAcceptedEvents();
			m_sumOfWeights        = allEventsCBK->sumOfEventWeights();
			m_sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
			Info(APP_NAME, "CutBookkeepers Accepted %lu SumWei %f sumWei2 %f ", m_nEventsProcessed, m_sumOfWeights, m_sumOfWeightsSquared);
		}
		else { Info(APP_NAME, "No relevent CutBookKeepers found"); }

		const xAOD::EventInfo* eventInfo = 0;
		if (!tevent->retrieve(eventInfo, "EventInfo").isSuccess()) {
			Error("execute()", "Failed to retrieve event info collection. Exiting.");
			return EL::StatusCode::FAILURE;
		}
		auto dsidString = std::to_string(eventInfo->mcChannelNumber());

		nEventsProcessedHist->Fill(dsidString.c_str(), m_nEventsProcessed);
		sumOfWeightsHist->Fill(dsidString.c_str(), m_sumOfWeights);
		sumOfWeightsSquaredHist->Fill(dsidString.c_str(), m_sumOfWeightsSquared);
	}
	else {
		nEventsProcessedHist->Fill("-1", 1);
		sumOfWeightsHist->Fill("-1", 1);
		sumOfWeightsSquaredHist->Fill("-1", 1);
	}
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::initialize()
{
	// Here you do everything that you need to do after the first input
	// file has been connected and before the first event is processed,
	// e.g. create additional histograms based on which variables are
	// available in the input files.  You can also create all of your
	// histograms and trees in here, but be aware that this method
	// doesn't get called if no events are processed.  So any objects
	// you create here won't be available in the output if you have no
	// input events.
	m_event = wk()->xaodEvent();

	Info(APP_NAME, "Number of events = %lli", m_event->getEntries()); // print long long int

	m_store	   = wk()->xaodStore();

	std::string sampleName = wk()->metaData()->castString("nc_source");

	getProjectTag(sampleName, m_projectTag, m_isData);
	getTags(sampleName, m_aTags, m_sTags, m_rTags, m_eTags, m_pTags);

	m_isAtlFast  = m_aTags.size()>0;
	m_isDerived  = sampleName.find("DAOD") != std::string::npos;

	ST::SettingDataSource dataSource = m_isData ? ST::Data : (m_isAtlFast ? ST::AtlfastII : ST::FullSim);

	m_outputFile = wk()->getOutputFile(m_outputLabel);

	std::string mainDir = getenv("ROOTCOREBIN");
	// Setup SUSYTools object
	m_objTool= new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
	m_objTool->msg().setLevel(MSG::ERROR);
	//for (auto&& pair: m_objTool->getPropertyMgr()->getProperties())
	//{
	//	auto tprop = reinterpret_cast<const TProperty<ToolHandle<asg::AsgTool>>*>(pair.second);
	//	if (tprop->pointer() != nullptr)
	//	{
	//		std::cout << "Name " << pair.first << " RTTI " << typeid(*reinterpret_cast<const void*>(tprop->pointer())).name() << std::endl;
	//		//(*tool)->msg().setLevel(MSG::ERROR);
	//	}
	//}
	CHECK(m_objTool->setProperty("DataSource", dataSource));
	CHECK(m_objTool->setProperty("ConfigFile", mainDir+"/data/Stop0LRun2/"+m_configFile));
	CHECK(m_objTool->setProperty("METDoTrkSyst", true));
	CHECK(m_objTool->setProperty("METDoCaloSyst", false));
	CHECK(m_objTool->setProperty("IsMC15b", m_isMC15b));

	// Pile Up Reweighting
	std::vector<std::string> prw_conf;
	prw_conf.push_back(mainDir+"/data/Stop0LRun2/merged_prw.root");
	//The fallback mu distribution. Not ok for publication
	//prw_conf.push_back("dev/PileupReweighting/mc15a_defaults.NotRecommended.prw.root");
	std::vector<std::string> prw_lumicalc;
	prw_lumicalc.push_back(mainDir+"/data/Stop0LRun2/ilumicalc_histograms_None_276262-284484.root");

	CHECK(m_objTool->setProperty("PRWConfigFiles", prw_conf));
	CHECK(m_objTool->setProperty("PRWLumiCalcFiles", prw_lumicalc));
	//CHECK( m_objTool->setProperty("PRWDefaultChannel", 410000) );

	if (m_debug) m_objTool->msg().setLevel(MSG::VERBOSE);
	CHECK(m_objTool->setProperty("DebugMode", m_debug));

	if (m_objTool->initialize() != EL::StatusCode::SUCCESS) {
		Error(APP_NAME, "Cannot intialize SUSYObjDef_xAOD...");
		Error(APP_NAME, "Exiting... ");
		return EL::StatusCode::FAILURE;
	}
	else {
		Info(APP_NAME, "SUSYObjDef_xAOD initialized... ");
	}
	//Max: Manually turn off the debug flag in PRW tool
	auto prwTool = m_objTool->getProperty<ToolHandle<CP::IPileupReweightingTool>>("PileupReweightingTool");
	(*prwTool)->expert()->EnableDebugging(m_debug);
	//dynamic_cast<CP::PileupReweightingTool&>(**prwTool).msg().setLevel(MSG::WARNING);

#if DOPFLOW
	InitPFlowSUSYTools(dataSource);
	if (m_objTool->SUSYToolsInit().isFailure()) {
		Error(APP_NAME, "Failed to initialise tools in SUSYToolsInit()...");
		Error(APP_NAME, "Exiting...");
		return EL::StatusCode::FAILURE;
	}
#endif // DOPFLOW


	if(!m_isData)
		m_xSecDB = new SUSY::CrossSectionDB(mainDir+"/data/SUSYTools/"+m_projectTag+"/", false); // Don't use PathResolver
	truthNeutrinos.reset(new Ntup4Vec());
	truthWs.reset(new Ntup4Vec());
	truthBottoms.reset(new Ntup4Vec());
	trueMetParticle.reset(new Ntup4Vec());

	//Extra jet collections
	extraJetNtup4VecsMap.reset(new std::unordered_map<std::string, std::unique_ptr<Ntup4Vec> >());
	ntupJetExtraMap.reset(new std::unordered_map<std::string, std::unique_ptr<NtupJetExtraVec> >());

	std::vector<ST::SystInfo> systInfoList;
	if (m_syst) {
		systInfoList = m_objTool->getSystInfoList();
	}
	else {
		ST::SystInfo infodef;
		infodef.affectsKinematics = false;
		infodef.affectsWeights = false;
		infodef.affectsType = ST::Unknown;
		systInfoList.push_back(infodef);
	}

	// for (const auto& systInfo : systInfoList) {
	// 	const CP::SystematicSet& syst = systInfo.systset;
	// 	if(systInfo.affectsWeights && !systInfo.affectsKinematics){
	// 	}
	// }

	nEventsProcessedHist->SetDirectory((TDirectory*)m_outputFile);
	sumOfWeightsHist->SetDirectory((TDirectory*)m_outputFile);
	sumOfWeightsSquaredHist->SetDirectory((TDirectory*)m_outputFile);

	std::vector< std::unique_ptr<TwoLeptonProcessor> > twoLepProcessorOneReg;

	for (unsigned int rI = 0; rI<3; rI++) {
		std::vector< TTree*> tempVecTree;
		m_outTrees.push_back(tempVecTree);

		std::vector< TH1F*> tempHist;
		std::vector< TH1D*> tempWHist;
		m_cutflow.push_back(tempHist);
		m_cutflowWeighted.push_back(tempWHist);

		std::vector<ST::SystInfo> tempSystInfo;
		m_systInfoLists.push_back(tempSystInfo);
	}

	twoLeptonProcessors.resize(m_outTrees.size()); //Only [2] is used
	twoLeptonProcessors.push_back(std::move(twoLepProcessorOneReg));

	for (const auto& systInfo : systInfoList) {
		const CP::SystematicSet& syst = systInfo.systset;
		std::string suffix = "";
		if (syst.name() != "")
			suffix = "__SYST__"+syst.name();
		// Initialize histograms and trees here so we can give them names based on the input samples

		bool syst_affectsElectrons = ST::SystObjType::Electron==systInfo.affectsType;
		bool syst_affectsMuons = ST::SystObjType::Muon==systInfo.affectsType;
		bool syst_affectsTaus = ST::SystObjType::Tau==systInfo.affectsType;
		bool syst_affectsPhotons = ST::SystObjType::Photon==systInfo.affectsType;
		bool syst_affectsJets = ST::SystObjType::Jet==systInfo.affectsType;
		bool syst_affectsBTag = ST::SystObjType::BTag==systInfo.affectsType;
		bool syst_affectsEventWeight = ST::SystObjType::EventWeight==systInfo.affectsType;
		bool syst_affectsMET = ST::SystObjType::MET_CST==systInfo.affectsType||ST::SystObjType::MET_TST==systInfo.affectsType||ST::SystObjType::MET_Track==systInfo.affectsType;

		if (systInfo.affectsWeights && !systInfo.affectsKinematics) {
			if ((!syst_affectsTaus && !syst_affectsPhotons) && (syst_affectsJets || syst_affectsMuons || syst_affectsElectrons || syst_affectsMET || syst_affectsBTag || syst_affectsEventWeight)) {
				m_weightSystInfo.push_back(systInfo);
				m_weights[syst.name()] = 1.0;
			}
			continue;
		}

		if ((!syst_affectsTaus && !syst_affectsPhotons) && (syst_affectsJets || syst_affectsMuons || syst_affectsElectrons || syst_affectsMET || syst.name()=="")) {

			if ((!syst_affectsMuons && !syst_affectsElectrons) && (syst_affectsJets|| syst_affectsMET || syst.name()=="")) {
				m_outTrees[0].push_back(new TTree(("stop_0Lep"+suffix).c_str(), "0Lep tree"));
				m_cutflow[0].push_back(new TH1F(("cutflow0Lep"+suffix).c_str(), "", 100, 0, 100));
				m_cutflowWeighted[0].push_back(new TH1D(("cutflowWeighted0Lep"+suffix).c_str(), "", 100, 0, 100));
				m_systInfoLists[0].push_back(systInfo);
			}

			m_outTrees[1].push_back(new TTree(("stop_1Lep"+suffix).c_str(), "1Lep tree"));
			m_outTrees[2].push_back(new TTree(("stop_2Lep"+suffix).c_str(), "2Lep tree"));

			m_cutflow[1].push_back(new TH1F(("cutflow1Lep"+suffix).c_str(), "", 100, 0, 100));
			m_cutflowWeighted[1].push_back(new TH1D(("cutflowWeighted1Lep"+suffix).c_str(), "", 100, 0, 100));
			m_cutflow[2].push_back(new TH1F(("cutflow2Lep"+suffix).c_str(), "", 100, 0, 100));
			m_cutflowWeighted[2].push_back(new TH1D(("cutflowWeighted2Lep"+suffix).c_str(), "", 100, 0, 100));
			std::unique_ptr<TwoLeptonProcessor> tempTwoLepProc;

			twoLeptonProcessors[2].push_back(std::move(tempTwoLepProc));

			m_systInfoLists[1].push_back(systInfo);
			m_systInfoLists[2].push_back(systInfo);
		}
	}
	twoLeptonProcessors[0].resize(m_outTrees[0].size()); //Only [2] is used
	twoLeptonProcessors[1].resize(m_outTrees[1].size()); //Only [2] is used

	if (m_isData) {
		//Number of events by runnummber
		m_runNumberNEvents = new TH1I("RunNumberNEvents", "RunNumberNEvents", 100000, 276000, 376000);
		m_runNumberNEvents->SetDirectory((TDirectory*)m_outputFile);
	}

	for (unsigned int systI = 0; systI<m_cutflow[0].size(); systI++) {
		// SRA
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(1, "All events");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(2, "GRL");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(3, "LAr and Tile error");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(4, "Trigger");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(5, "Primary vertex exists");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(6, "Clean jet");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(7, "Cosmic");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(9, "Lepton veto");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(8, "Clean Muon");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(10, "Trigger match");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(11, "Signal jet (|eta|<2.8) pT (80:80)");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(12, "MET>250");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(13, "NJets>=6");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(14, "dPhiMetJetMin3>pi/5");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(15, "MetTrack>30");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(16, "dPhiMetTrackMet<pi/3");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(17, "TwoBtags");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(18, "tau-veto");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(19, "MtBMin>175");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(20, "50<TopDRminB0<250");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(21, "50<TopDRminB1<400");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(22, "minMt>50");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23, "highMet>350");

		// SRB has an offset in bin index
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+0, "MtBMin>175 (SRB check)");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+1, "AntiKt12DijetMassAsym<0.5");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+2, "80<AntiKt12M[0]");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+3, "60<AntiKt12M[1]<200");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+4, "50<AntiKt8M[0]");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+5, "MtNonBMin>175");
		m_cutflow[0][systI]->GetXaxis()->SetBinLabel(23+m_srOffset+6, "HtSig>17");

		for (int binI=0; binI < m_cutflow[0][systI]->GetXaxis()->GetNbins(); binI++) {
			m_cutflowWeighted[0][systI]->GetXaxis()->SetBinLabel(binI+1, m_cutflow[0][systI]->GetXaxis()->GetBinLabel(binI+1));
		}
		m_cutflow[0][systI]->SetDirectory((TDirectory*)m_outputFile);
		m_cutflowWeighted[0][systI]->SetDirectory((TDirectory*)m_outputFile);
	}

	for (unsigned int systI = 0; systI<m_cutflow[1].size(); systI++) {
		// 1-lepton CR
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(1, "All events");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(2, "GRL");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(3, "LAr and Tile error");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(4, "Trigger");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(5, "Primary vertex exists");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(6, "Clean jet");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(7, "Cosmic");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(9, "1-leptons");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(8, "Clean Muon");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(10, "Trigger match");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(11, "Signal jet (|eta|<2.8) pT (80:80)");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(12, "MET>250");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(13, "NJets>=6");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(14, "dPhiMetJetMin3>pi/5");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(15, "MetTrack>30");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(16, "dPhiMetTrackMet<pi/3");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(17, "TwoBtags");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(18, "tau-veto");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(19, "MtBMin>175");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(20, "50<TopDRminB0<250");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(21, "50<TopDRminB1<400");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(22, "minMt>50");
		m_cutflow[1][systI]->GetXaxis()->SetBinLabel(23, "highMet>350");

		for (int binI=0; binI < m_cutflow[1][systI]->GetXaxis()->GetNbins(); binI++) {
			m_cutflowWeighted[1][systI]->GetXaxis()->SetBinLabel(binI+1, m_cutflow[1][systI]->GetXaxis()->GetBinLabel(binI+1));
		}
		m_cutflow[1][systI]->SetDirectory((TDirectory*)m_outputFile);
		m_cutflowWeighted[1][systI]->SetDirectory((TDirectory*)m_outputFile);
	}

	//for (auto cutflowHist: m_cutflow){
	for (unsigned int systI = 0; systI<m_cutflow[2].size(); systI++) {
		// 2-lepton 
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(1, "All events");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(2, "GRL");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(3, "LAr and Tile error");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(4, "Trigger");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(5, "Primary vertex exists");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(6, "Clean jet");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(7, "Cosmic");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(9, "2 lepton req");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(8, "Clean Muon");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(10, "Trigger match");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(11, "Signal jet (|eta|<2.8) pT (80:80)");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(12, "MET<50");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(13, "OS leps");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(14, "86<Mll<96");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(15, "MetLep>70");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(16, "NJets>=4");
		m_cutflow[2][systI]->GetXaxis()->SetBinLabel(17, "TwoBtags");

		for (int binI=0; binI < m_cutflow[2][systI]->GetXaxis()->GetNbins(); binI++) {
			m_cutflowWeighted[2][systI]->GetXaxis()->SetBinLabel(binI+1, m_cutflow[2][systI]->GetXaxis()->GetBinLabel(binI+1));
		}

		m_cutflow[2][systI]->SetDirectory((TDirectory*)m_outputFile);
		m_cutflowWeighted[2][systI]->SetDirectory((TDirectory*)m_outputFile);
	}

	for (unsigned int rI = 0; rI < m_outTrees.size(); rI++) {
		for (auto& weightMap : m_weights)
			m_outTrees[rI][0]->Branch(("SYST__"+weightMap.first).c_str(), &weightMap.second);

		for (unsigned int systI = 0; systI < m_outTrees[rI].size(); systI++) {
			m_outTrees[rI][systI]->SetDirectory((TDirectory*)m_outputFile);

			m_outTrees[rI][systI]->Branch("IsZll", &m_isZll);
			m_outTrees[rI][systI]->Branch("EventNumber", &m_eventNumber);
			m_outTrees[rI][systI]->Branch("RunNumber", &m_runNumber);
			m_outTrees[rI][systI]->Branch("lb", &m_lumiBlock);
			m_outTrees[rI][systI]->Branch("dsid", &m_mcChannelNumber);
			m_outTrees[rI][systI]->Branch("GenWeight", &m_mcEventWeight);
			m_outTrees[rI][systI]->Branch("XSecWeight", &m_xSecWeight);
			m_outTrees[rI][systI]->Branch("PileupWeight", &m_pileupWeight);
			m_outTrees[rI][systI]->Branch("MuonSF", &m_muonSF);
			m_outTrees[rI][systI]->Branch("ElecSF", &m_elecSF);
			m_outTrees[rI][systI]->Branch("BTagSF", &m_bTagSF);

			m_outTrees[rI][systI]->Branch("AvgIntPerXing", &m_avgIntPerXing);

			m_outTrees[rI][systI]->Branch("aTags", &m_aTags);
			m_outTrees[rI][systI]->Branch("sTags", &m_sTags);
			m_outTrees[rI][systI]->Branch("rTags", &m_rTags);
			m_outTrees[rI][systI]->Branch("eTags", &m_eTags);
			m_outTrees[rI][systI]->Branch("pTags", &m_pTags);

			m_outTrees[rI][systI]->Branch("NElectrons", &m_nElectrons);
			m_outTrees[rI][systI]->Branch("NMuons", &m_nMuons);
			m_outTrees[rI][systI]->Branch("NSigElectrons", &m_nSigElectrons);
			m_outTrees[rI][systI]->Branch("NSigMuons", &m_nSigMuons);

			m_outTrees[rI][systI]->Branch("NBadMuons", &m_nBadMuons);
			m_outTrees[rI][systI]->Branch("NCosmicMuons", &m_nCosmicMuons);

			m_outTrees[rI][systI]->Branch("NJets", &m_nJets);
			m_outTrees[rI][systI]->Branch("NBadJets", &m_nBadJets);
			m_outTrees[rI][systI]->Branch("NBJets", &m_nBJets);
			m_outTrees[rI][systI]->Branch("JetPt", &m_jetPt);
			m_outTrees[rI][systI]->Branch("JetPhi", &m_jetPhi);
			m_outTrees[rI][systI]->Branch("JetEta", &m_jetEta);
			m_outTrees[rI][systI]->Branch("JetEMFrac", &m_jetEMFrac);
			m_outTrees[rI][systI]->Branch("JetHECFrac", &m_jetHECFrac);
			m_outTrees[rI][systI]->Branch("JetM", &m_jetM);
			m_outTrees[rI][systI]->Branch("JetIsSignal", &m_jetIsSignal);
			m_outTrees[rI][systI]->Branch("JetMV1Weight", &m_jetMV1Weight);
			m_outTrees[rI][systI]->Branch("JetIsBTagged", &m_jetIsBTagged);
			m_outTrees[rI][systI]->Branch("JetMt", &m_jetMt);
			m_outTrees[rI][systI]->Branch("JetDPhiMet", &m_jetDPhiMet);
			m_outTrees[rI][systI]->Branch("JetDPhiMetMin3", &m_jetDPhiMetMin3);
			m_outTrees[rI][systI]->Branch("JetLeadTagIndex", &m_jetLeadTagIndex);
			m_outTrees[rI][systI]->Branch("JetSubleadTagIndex", &m_jetSubleadTagIndex);

			m_outTrees[rI][systI]->Branch("TopDRMin0Pt", &m_topDRMin0Pt);
			m_outTrees[rI][systI]->Branch("TopDRMin0Phi", &m_topDRMin0Phi);
			m_outTrees[rI][systI]->Branch("TopDRMin0Eta", &m_topDRMin0Eta);
			m_outTrees[rI][systI]->Branch("TopDRMin0M", &m_topDRMin0M);

			m_outTrees[rI][systI]->Branch("TopDRMin1Pt", &m_topDRMin1Pt);
			m_outTrees[rI][systI]->Branch("TopDRMin1Phi", &m_topDRMin1Phi);
			m_outTrees[rI][systI]->Branch("TopDRMin1Eta", &m_topDRMin1Eta);
			m_outTrees[rI][systI]->Branch("TopDRMin1M", &m_topDRMin1M);

			m_outTrees[rI][systI]->Branch("TopDRMinB0Pt", &m_topDRMinB0Pt);
			m_outTrees[rI][systI]->Branch("TopDRMinB0Phi", &m_topDRMinB0Phi);
			m_outTrees[rI][systI]->Branch("TopDRMinB0Eta", &m_topDRMinB0Eta);
			m_outTrees[rI][systI]->Branch("TopDRMinB0M", &m_topDRMinB0M);

			m_outTrees[rI][systI]->Branch("TopDRMinB1Pt", &m_topDRMinB1Pt);
			m_outTrees[rI][systI]->Branch("TopDRMinB1Phi", &m_topDRMinB1Phi);
			m_outTrees[rI][systI]->Branch("TopDRMinB1Eta", &m_topDRMinB1Eta);
			m_outTrees[rI][systI]->Branch("TopDRMinB1M", &m_topDRMinB1M);

			m_outTrees[rI][systI]->Branch("NAntiKt8", &m_nAntiKt8);
			m_outTrees[rI][systI]->Branch("AntiKt8DijetMassAsym", &m_antiKt8DijetMassAsym);
			m_outTrees[rI][systI]->Branch("AntiKt8Pt", &m_antiKt8Pt);
			m_outTrees[rI][systI]->Branch("AntiKt8Phi", &m_antiKt8Phi);
			m_outTrees[rI][systI]->Branch("AntiKt8Eta", &m_antiKt8Eta);
			m_outTrees[rI][systI]->Branch("AntiKt8M", &m_antiKt8M);
			m_outTrees[rI][systI]->Branch("AntiKt8Const", &m_antiKt8Const);

			m_outTrees[rI][systI]->Branch("NAntiKt12", &m_nAntiKt12);
			m_outTrees[rI][systI]->Branch("AntiKt12DijetMassAsym", &m_antiKt12DijetMassAsym);
			m_outTrees[rI][systI]->Branch("AntiKt12Pt", &m_antiKt12Pt);
			m_outTrees[rI][systI]->Branch("AntiKt12Phi", &m_antiKt12Phi);
			m_outTrees[rI][systI]->Branch("AntiKt12Eta", &m_antiKt12Eta);
			m_outTrees[rI][systI]->Branch("AntiKt12M", &m_antiKt12M);
			m_outTrees[rI][systI]->Branch("AntiKt12Const", &m_antiKt12Const);

			m_outTrees[rI][systI]->Branch("MetX", &m_metX);
			m_outTrees[rI][systI]->Branch("MetY", &m_metY);
			m_outTrees[rI][systI]->Branch("MetPhi", &m_metPhi);
			m_outTrees[rI][systI]->Branch("Met", &m_met);
			m_outTrees[rI][systI]->Branch("MetSumEt", &m_metSumEt);

			m_outTrees[rI][systI]->Branch("PassMetTrig", &m_passMetTrigger);
			m_outTrees[rI][systI]->Branch("MetXOrig", &m_metXOrig);
			m_outTrees[rI][systI]->Branch("MetYOrig", &m_metYOrig);
			m_outTrees[rI][systI]->Branch("MetPhiOrig", &m_metPhiOrig);
			m_outTrees[rI][systI]->Branch("MetOrig", &m_metOrig);
			m_outTrees[rI][systI]->Branch("MetSumEtOrig", &m_metSumEtOrig);

			m_outTrees[rI][systI]->Branch("MetXTruth", &m_metXTruth);
			m_outTrees[rI][systI]->Branch("MetYTruth", &m_metYTruth);
			m_outTrees[rI][systI]->Branch("MetPhiTruth", &m_metPhiTruth);
			m_outTrees[rI][systI]->Branch("MetTruth", &m_metTruth);
			m_outTrees[rI][systI]->Branch("MetSumEtTruth", &m_metSumEtTruth);

			m_outTrees[rI][systI]->Branch("RawMetXTrack", &m_rawMetXTrk);
			m_outTrees[rI][systI]->Branch("RawMetYTrack", &m_rawMetYTrk);
			m_outTrees[rI][systI]->Branch("RawMetTrackPhi", &m_rawMetPhiTrk);
			m_outTrees[rI][systI]->Branch("RawMetTrack", &m_rawMetTrk);

			m_outTrees[rI][systI]->Branch("MetXTrack", &m_metXTrk);
			m_outTrees[rI][systI]->Branch("MetYTrack", &m_metYTrk);
			m_outTrees[rI][systI]->Branch("MetTrackPhi", &m_metPhiTrk);
			m_outTrees[rI][systI]->Branch("MetTrack", &m_metTrk);
			
			m_outTrees[rI][systI]->Branch("MetLepX", &m_metLepX);
			m_outTrees[rI][systI]->Branch("MetLepY", &m_metLepY);
			m_outTrees[rI][systI]->Branch("MetLepPhi", &m_metLepPhi);
			m_outTrees[rI][systI]->Branch("MetLep", &m_metLep);
			m_outTrees[rI][systI]->Branch("MetLepSumEt", &m_metLepSumEt);
			
			m_outTrees[rI][systI]->Branch("Ht", &m_ht);
			m_outTrees[rI][systI]->Branch("HtSig", &m_htSig);
			m_outTrees[rI][systI]->Branch("MetSig", &m_metSig);
			m_outTrees[rI][systI]->Branch("MEff", &m_mEff);

			m_outTrees[rI][systI]->Branch("PrimVertexNTracks", &m_primVtxNTrks);
			m_outTrees[rI][systI]->Branch("NVtx", &m_nvtx);
			
			m_outTrees[rI][systI]->Branch("MtNonBMin", &m_mtNonBMin);
			m_outTrees[rI][systI]->Branch("MtBMin", &m_mtBMin);
			m_outTrees[rI][systI]->Branch("MinMt", &m_minMt);
			m_outTrees[rI][systI]->Branch("MtTauCand", &m_mtTauCand);
			m_outTrees[rI][systI]->Branch("TauJetNTracks", &m_tauJetNTracks);

			m_outTrees[rI][systI]->Branch("mcPt", &m_mcPt);
			m_outTrees[rI][systI]->Branch("mcPhi", &m_mcPhi);
			m_outTrees[rI][systI]->Branch("mcEta", &m_mcEta);
			m_outTrees[rI][systI]->Branch("mcE", &m_mcE);
			m_outTrees[rI][systI]->Branch("mcPdgId", &m_mcPdgId);
			m_outTrees[rI][systI]->Branch("mcBarcode", &m_mcBarcode);
			m_outTrees[rI][systI]->Branch("mcChildren", &m_mcChildren);

			m_outTrees[rI][systI]->Branch("NTruthJets", &m_nTruthJets);
			m_outTrees[rI][systI]->Branch("TruthJetPt", &m_truthJetPt);
			m_outTrees[rI][systI]->Branch("TruthJetPhi", &m_truthJetPhi);
			m_outTrees[rI][systI]->Branch("TruthJetEta", &m_truthJetEta);
			m_outTrees[rI][systI]->Branch("TruthJetM", &m_truthJetM);

			m_outTrees[rI][systI]->Branch("DRBBJetHt", &m_dRbbJetHt);
			m_outTrees[rI][systI]->Branch("DRBBJet", &m_dRbbJet);

			m_outTrees[rI][systI]->Branch("MBB", &m_mBB);
			m_outTrees[rI][systI]->Branch("DRBB", &m_dRbb);

			m_outTrees[rI][systI]->Branch("DPhiMetTrackMet", &m_dPhiMetTrackMet);

			m_outTrees[rI][systI]->Branch("AllJetMass", &m_allJetMass);
			m_outTrees[rI][systI]->Branch("AllJetPt", &m_allJetPt);
			m_outTrees[rI][systI]->Branch("AllJetPhi", &m_allJetPhi);
			m_outTrees[rI][systI]->Branch("AllJetEta", &m_allJetEta);

			m_outTrees[rI][systI]->Branch("TrueTopPt", &m_trueTopPt);
			m_outTrees[rI][systI]->Branch("TrueTopEta", &m_trueTopEta);
			m_outTrees[rI][systI]->Branch("TrueTopPhi", &m_trueTopPhi);
			m_outTrees[rI][systI]->Branch("TrueTopE", &m_trueTopE);

			m_outTrees[rI][systI]->Branch("TrueTopDaughtPdgId", &m_trueTopDaughtPdgId);
			m_outTrees[rI][systI]->Branch("TrueTopDaughtPt", &m_trueTopDaughtPt);
			m_outTrees[rI][systI]->Branch("TrueTopDaughtEta", &m_trueTopDaughtEta);
			m_outTrees[rI][systI]->Branch("TrueTopDaughtPhi", &m_trueTopDaughtPhi);
			m_outTrees[rI][systI]->Branch("TrueTopDaughtE", &m_trueTopDaughtE);

			m_outTrees[rI][systI]->Branch("TrueAntitopPt", &m_trueAntitopPt);
			m_outTrees[rI][systI]->Branch("TrueAntitopEta", &m_trueAntitopEta);
			m_outTrees[rI][systI]->Branch("TrueAntitopPhi", &m_trueAntitopPhi);
			m_outTrees[rI][systI]->Branch("TrueAntitopE", &m_trueAntitopE);

			m_outTrees[rI][systI]->Branch("TrueAntitopDaughtPdgId", &m_trueAntitopDaughtPdgId);
			m_outTrees[rI][systI]->Branch("TrueAntitopDaughtPt", &m_trueAntitopDaughtPt);
			m_outTrees[rI][systI]->Branch("TrueAntitopDaughtEta", &m_trueAntitopDaughtEta);
			m_outTrees[rI][systI]->Branch("TrueAntitopDaughtPhi", &m_trueAntitopDaughtPhi);
			m_outTrees[rI][systI]->Branch("TrueAntitopDaughtE", &m_trueAntitopDaughtE);

			m_outTrees[rI][systI]->Branch("TopDecayType", &m_topDecayType);
			m_outTrees[rI][systI]->Branch("AntitopDecayType", &m_antitopDecayType);

			m_outTrees[rI][systI]->Branch("TruthHt", &truthHt);
			m_outTrees[rI][systI]->Branch("TruthNeutrinoMet", &truthNeutrinoMet);

			m_outTrees[rI][systI]->Branch("PassMetTrigger", &m_passMetTrigger);
			m_outTrees[rI][systI]->Branch("MetTriggerPrescale", &m_metTriggerPrescale);

			m_outTrees[rI][systI]->Branch("PassLepTrigger", &m_passLepTrigger);
			m_outTrees[rI][systI]->Branch("LepTriggerPrescale", &m_lepTriggerPrescale);
			
			m_outTrees[rI][systI]->Branch("PassPhotonTrigger", &m_passPhotonTrigger);
			m_outTrees[rI][systI]->Branch("PhotonTriggerPrescale", &m_photonTriggerPrescale);

			m_outTrees[rI][systI]->Branch("LepTrigMatch", &m_passTrigMatch);
			m_outTrees[rI][systI]->Branch("OSLep", &m_passOSLep);

			m_outTrees[rI][systI]->Branch("LepPt", &m_lepPt);
			m_outTrees[rI][systI]->Branch("LepE", &m_lepE);
			m_outTrees[rI][systI]->Branch("LepPhi", &m_lepPhi);
			m_outTrees[rI][systI]->Branch("LepEta", &m_lepEta);
			m_outTrees[rI][systI]->Branch("LepType", &m_lepType);

			m_outTrees[rI][systI]->Branch("MtMetLep", &m_mtMetLep);
			m_outTrees[rI][systI]->Branch("dPhiMetLep", &m_dPhiMetLep);
			
			m_outTrees[rI][systI]->Branch("LLMass", &m_llMass);
            
			m_outTrees[rI][systI]->Branch("Mt2_Truth", &mt2_truth);			
			m_outTrees[rI][systI]->Branch("Mt2_RCJet", &mt2_rcjet);			
			m_outTrees[rI][systI]->Branch("Mt2_DRMin", &mt2_drmin);			

			truthNeutrinos->ConnectTree(*m_outTrees[rI][systI], "TruthTopNeutrino");
            trueMetParticle->ConnectTree(*m_outTrees[rI][systI], "TruthMetParticles");
			truthWs->ConnectTree(*m_outTrees[rI][systI], "TruthW");
			truthBottoms->ConnectTree(*m_outTrees[rI][systI], "TruthBottom");

			//TwoLeptonProcessor for non-2L-tree are just dummies whose sole purpose is to create the tree branches structure.
			twoLeptonProcessors[rI][systI].reset(new TwoLeptonProcessor(*m_outTrees[rI][systI], *m_objTool));
		}
	}

	if (m_objTool->resetSystematics() != CP::SystematicCode::Ok) {
		Error(APP_NAME, "Cannot reset SUSYTools systematics");
		exit(-2);
	}

	for (auto jetCollection : extraJetCollections)
		{
			auto ntup4Vec = std::make_unique<Ntup4Vec>();
			ntup4Vec->ConnectTree(*m_outTrees[0][0], jetCollection);
			extraJetNtup4VecsMap->emplace(jetCollection, std::move(ntup4Vec));

			auto ntupJetExtraVec =std::make_unique<NtupJetExtraVec>();
			ntupJetExtraVec->ConnectTree(*m_outTrees[0][0], jetCollection);
			ntupJetExtraMap->emplace(jetCollection, std::move(ntupJetExtraVec));
		}
	//Extra jet variable for the nominal
	auto ntupJetExtraVec = std::make_unique<NtupJetExtraVec>();
	ntupJetExtraVec->ConnectTree(*m_outTrees[0][0], nominalJetCollection);
	ntupJetExtraMap->emplace(nominalJetCollection, std::move(ntupJetExtraVec));
	// GRL
	m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
	std::vector<std::string> vecStringGRL;
	vecStringGRL.push_back(mainDir+"/data/Stop0LRun2/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml");
	CHECK(m_grl->setProperty("GoodRunsListVec", vecStringGRL));
	CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
	if (!m_grl->initialize().isSuccess()) { // check this isSuccess
		Error(APP_NAME, "Failed to properly initialize the GRL. Exiting.");
		return EL::StatusCode::FAILURE;
	}

	m_mcTruthClassifier = new MCTruthClassifier("MCTruthClassifier");
#ifdef DORJIGSAW
	rJigsaw = std::make_unique<RJigsaw>();
	rJigsaw->ConnectTree(*m_outTrees[0][0]);
#endif
	for (auto& jetCollection : extraJetCollections)
		{
			rScanCalibTools.emplace(jetCollection, std::unique_ptr<JetCalibrationTool>(new JetCalibrationTool(
				                            "RScanCalibration"+jetCollection,
				                            jetCollection,
				                            "JES_MC15Prerecommendation_Rscan_Aug2015.config",
				                            "JetArea_Origin",
				                            m_isData)));
			CHECK(rScanCalibTools.at(jetCollection)->initializeTool(jetCollection));
		}
	jetCleaningTool.reset(new JetCleaningTool("JetCleaningToolRscan"));
	CHECK(jetCleaningTool->setProperty("CutLevel", "LooseBad"));
	CHECK(jetCleaningTool->initialize());
	jvtTool.reset(new JetVertexTaggerTool("JetVertexTaggerToolRscan"));
	CHECK(jvtTool->setProperty("JVTFileName", "JetMomentTools/JVTlikelihood_20140805.root"));
	CHECK(jvtTool->initialize());

#if DOSUBSTRUCTURE
	trackSelector.reset(new InDet::InDetTrackSelectionTool("TrackSelection"));
	CHECK(trackSelector->setProperty("CutLevel", "TightPrimary"));
	CHECK(trackSelector->initialize());

	d12Calculator.reset(new JetSubStructureUtils::KtSplittingScale(1));
	d23Calculator.reset(new JetSubStructureUtils::KtSplittingScale(2));
	//tau1Calculator.reset(new fastjet::contrib::Nsubjettiness(
	//	1,
	//	fastjet::contrib::                    
	//	fastjet::contrib::Njettiness::unnormalized_measure)  
#endif // DOSUBSTRUCTURE


	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::execute()
{
	// Here you do everything that needs to be done on every single
	// events, e.g. read input variables, apply cuts, and fill
	// histograms and trees.  This is where most of your actual analysis
	// code will go.

	m_eventNumber	    = INIT_VAL;
	m_runNumber       = INIT_VAL;
	m_lumiBlock       = INIT_VAL;
	m_mcChannelNumber = INIT_VAL;
	m_avgIntPerXing   = INIT_VAL;

	// Weights are special, they should be 1
	m_mcEventWeight = 1;
	m_xSecWeight	  = 1;
	m_pileupWeight  = 1;
	m_elecSF = 1;
	m_muonSF = 1;
	m_bTagSF = 1;

	m_elecSF_nominal = 1;
	m_muonSF_nominal = 1;
	m_bTagSF_nominal = 1;

	const xAOD::EventInfo* eventInfo = 0;
	if (!m_event->retrieve(eventInfo, "EventInfo").isSuccess()) {
		Error("execute()", "Failed to retrieve event info collection. Exiting.");
		return EL::StatusCode::FAILURE;
	}
	// Keep track of total number of events that are processed.
	// Events without eventInfo are not considered since this is as if you never downloaded them

	m_eventNumber	  = eventInfo->eventNumber();
	m_runNumber	    = eventInfo->runNumber();
	m_lumiBlock	    = eventInfo->lumiBlock();
	m_avgIntPerXing = eventInfo->averageInteractionsPerCrossing();

	// If we are printing event information for a particular event number just move on
	if (m_eventNumber != m_evtToPrint && m_evtToPrint >= 0)
		{
			return EL::StatusCode::SUCCESS;
		}
	// All events
	for (auto& cutflowVec : m_cutflow) {
		for (auto& cutflow :cutflowVec) {
			cutflow->Fill(0.5);
		}
	}

	// Check if we are dealing with data or MC;
	// if data set generator weights to 1 and exit if event doesn't pass GRL
	if (m_isData) {
		m_xSecWeight = 1.0;

		if (!m_grl->passRunLB(*eventInfo)) {
			return EL::StatusCode::SUCCESS; // go to next event
		}
		// Pass GRL
		for (unsigned int rI = 0; rI < m_cutflow.size(); rI++) {
			for (unsigned int systI = 0; systI<m_cutflow[rI].size(); systI++) {
				m_cutflow[rI][systI]->Fill(1.5);
				m_cutflowWeighted[rI][systI]->Fill(1.5);
			}
		}
		if ((eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error) ||
		    (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error) ||
		    (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)))
			{
				return EL::StatusCode::SUCCESS; // go to the next event
			} // end if event flags check

		// Passed LAr and Tile error
		for (unsigned int rI = 0; rI < m_cutflow.size(); rI++) {
			for (unsigned int systI = 0; systI<m_cutflow[rI].size(); systI++) {
				m_cutflow[rI][systI]->Fill(2.5);
				m_cutflowWeighted[rI][systI]->Fill(2.5);
			}
		}
	}
	else {
		m_mcChannelNumber = eventInfo->mcChannelNumber();
		m_mcEventWeight   = eventInfo->mcEventWeight();
		initTruthVars();
		getTruth();
		m_xSecWeight = m_xSecDB->xsectTimesEff(eventInfo->mcChannelNumber(), m_finalState);
		
		CHECK(m_objTool->ApplyPRWTool());
		m_pileupWeight = m_objTool->GetPileupWeight();

		// All events, GRL, and LAR/Tile
		for (unsigned int rI = 0; rI < m_cutflow.size(); rI++) {
			for (unsigned int systI = 0; systI<m_cutflow[rI].size(); systI++) {
				m_cutflowWeighted[rI][systI]->Fill(0.5, m_mcEventWeight*m_pileupWeight);
				m_cutflow[rI][systI]->Fill(1.5);
				m_cutflowWeighted[rI][systI]->Fill(1.5, m_mcEventWeight*m_pileupWeight);
				m_cutflow[rI][systI]->Fill(2.5);
				m_cutflowWeighted[rI][systI]->Fill(2.5, m_mcEventWeight*m_pileupWeight);
			}
		}
	}

	// Get the nominal object containers from the event
	// Electrons
	m_electrons_nominal = 0;
	m_electrons_nominal_aux = 0;
	CHECK(m_objTool->GetElectrons(m_electrons_nominal, m_electrons_nominal_aux, true));

	// // Photons
	// m_photons_nominal = 0;
	// m_photons_nominal_aux = 0;
	// CHECK( m_objTool->GetPhotons(m_photons_nominal,m_photons_nominal_aux, true) );

	// Muons
	m_muons_nominal = 0;
	m_muons_nominal_aux = 0;
	CHECK(m_objTool->GetMuons(m_muons_nominal, m_muons_nominal_aux, true));

	// Jets
	m_jets_nominal = 0;
	m_jets_nominal_aux = 0;
	CHECK(m_objTool->GetJets(m_jets_nominal, m_jets_nominal_aux, true));

	// // Taus
	// m_taus_nominal = 0;
	// m_taus_nominal_aux = 0;
	// CHECK( m_objTool->GetTaus(m_taus_nominal,m_taus_nominal_aux, true) );

	// MET
	m_metcst_nominal = new xAOD::MissingETContainer;
	m_metcst_nominal_aux = new xAOD::MissingETAuxContainer;
	m_metcst_nominal->setStore(m_metcst_nominal_aux);
	m_metcst_nominal->reserve(10);

	m_mettst_nominal = new xAOD::MissingETContainer;
	m_mettst_nominal_aux = new xAOD::MissingETAuxContainer;
	m_mettst_nominal->setStore(m_mettst_nominal_aux);
	m_mettst_nominal->reserve(10);

	m_mettrack_nominal = new xAOD::MissingETContainer;
	m_mettrack_nominal_aux = new xAOD::MissingETAuxContainer;
	m_mettrack_nominal->setStore(m_mettrack_nominal_aux);
	m_mettrack_nominal->reserve(10);

	m_metleptst_nominal = new xAOD::MissingETContainer;
	m_metleptst_nominal_aux = new xAOD::MissingETAuxContainer;
	m_metleptst_nominal->setStore(m_metleptst_nominal_aux);
	m_metleptst_nominal->reserve(10);

	CHECK(m_store->record(m_metcst_nominal, "MetCSTNominal"));
	CHECK(m_store->record(m_mettst_nominal, "MetTSTNominal"));
	CHECK(m_store->record(m_mettrack_nominal, "MetTRACKNominal"));
	CHECK(m_store->record(m_metcst_nominal_aux, "MetCSTNominalAux."));
	CHECK(m_store->record(m_mettst_nominal_aux, "MetTSTNominalAux."));
	CHECK(m_store->record(m_mettrack_nominal_aux, "MetTRACKNominalAux."));
	
	CHECK(m_store->record(m_metleptst_nominal, "MetLEPTSTNominal"));
	CHECK(m_store->record(m_metleptst_nominal_aux, "MetLEPTSTNominalAux."));

	if (m_isData) {
		//m_runNumberNEvents->Fill(TString::Format("%i", m_runNumber), 1);
		m_runNumberNEvents->Fill(m_runNumber);
	}
	m_isNominal = true;

	for (unsigned int rI = 0; rI<m_outTrees.size(); rI++) {
		for (unsigned int systI = 0; systI<m_outTrees[rI].size(); systI++) {
			const CP::SystematicSet& syst = m_systInfoLists[rI][systI].systset;

			if (syst.name() == "") {
				for (auto& weightSyst : m_weightSystInfo) {
					const CP::SystematicSet& syst = weightSyst.systset;
					bool syst_affectsElectrons = ST::SystObjType::Electron==weightSyst.affectsType;
					bool syst_affectsMuons = ST::SystObjType::Muon==weightSyst.affectsType;
					bool syst_affectsBTag = ST::SystObjType::BTag==weightSyst.affectsType;
					bool syst_affectsEventWeight = ST::SystObjType::EventWeight==weightSyst.affectsType;

					if (m_objTool->applySystematicVariation(syst) != CP::SystematicCode::Ok) {
						Error(APP_NAME, "Cannot configure SUSYTools for systematic var. %s", syst.name().c_str());
					}

					if (rI==1 || rI ==2) {
						if (syst_affectsElectrons)
							m_weights[syst.name()] = m_objTool->GetTotalElectronSF(*m_electrons_nominal);

						if (syst_affectsMuons)
							m_weights[syst.name()] = m_objTool->GetTotalMuonSF(*m_muons_nominal);
					}

					if (syst_affectsBTag)
						m_weights[syst.name()] = m_objTool->BtagSF(m_jets_nominal);

					//if(syst_affectsEventWeight){
					//	
					//}
				}
			}
			else {
				for (auto& weightSyst : m_weightSystInfo) {
					const CP::SystematicSet& syst = weightSyst.systset;
					m_weights[syst.name()] = 1.0;
				}
			}
			// Tell the SUSYObjDef_xAOD which variation to apply
			if (m_objTool->applySystematicVariation(syst) != CP::SystematicCode::Ok) {
				Error(APP_NAME, "Cannot configure SUSYTools for systematic var. %s", syst.name().c_str());
			}
			else {
				if (m_debug) Info(APP_NAME, "Variation \"%s\" configured...", syst.name().c_str());
			}

			m_affectsKinematics = m_systInfoLists[rI][systI].affectsKinematics;
			m_affectsWeights = m_systInfoLists[rI][systI].affectsWeights;
			if (m_affectsKinematics || m_affectsWeights) m_isNominal = false;

			// If necessary (kinematics affected), make a shallow copy with the variation applied
			m_syst_affectsElectrons = ST::testAffectsObject(xAOD::Type::Electron, m_systInfoLists[rI][systI].affectsType);
			m_syst_affectsMuons = ST::testAffectsObject(xAOD::Type::Muon, m_systInfoLists[rI][systI].affectsType);
			m_syst_affectsTaus = ST::testAffectsObject(xAOD::Type::Tau, m_systInfoLists[rI][systI].affectsType);
			m_syst_affectsPhotons = ST::testAffectsObject(xAOD::Type::Photon, m_systInfoLists[rI][systI].affectsType);
			m_syst_affectsJets = ST::testAffectsObject(xAOD::Type::Jet, m_systInfoLists[rI][systI].affectsType);
			m_syst_affectsBTag = ST::testAffectsObject(xAOD::Type::BTag, m_systInfoLists[rI][systI].affectsType);
			//      m_syst_affectsMET = ST::testAffectsObject(xAOD::Type::MissingET, m_systInfoLists[rI][systI].affectsType);

			initVars();

			// Create jets, muons, and electron objects
			auto xAODObjects = createObj(syst.name());

			GetExtraJetCollections();

			// This function tells you what channel we are in, 0=0-lepton, 1=1-lepton, 2=2-lepton. Negative values indicate that no 
			int chanNum = isInChan();

			// If we are in the 1-lepton channel add the lepton to the jet collection
			if (chanNum==1 && (int)rI==chanNum){
				Stop0LUtils::Jet* myJet;
				
				if(m_nSigMuons==1){
					m_passTrigMatch = ((Stop0LUtils::Muon*)(*m_muons)[0])->IsTrigMatched();
					myJet = new Stop0LUtils::Jet(((Stop0LUtils::Muon*)(*m_muons)[0])->Pt()*1000.,
					                             ((Stop0LUtils::Muon*)(*m_muons)[0])->Eta(),
					                             ((Stop0LUtils::Muon*)(*m_muons)[0])->Phi(),
					                             ((Stop0LUtils::Muon*)(*m_muons)[0])->E()*1000.,
					                             0);
					myJet->SetLeptonType(0);				
				}
				else{
					m_passTrigMatch = ((Stop0LUtils::Electron*)(*m_electrons)[0])->IsTrigMatched();
					myJet = new Stop0LUtils::Jet(((Stop0LUtils::Electron*)(*m_electrons)[0])->Pt()*1000.,
					                             ((Stop0LUtils::Electron*)(*m_electrons)[0])->Eta(),
					                             ((Stop0LUtils::Electron*)(*m_electrons)[0])->Phi(),
					                             ((Stop0LUtils::Electron*)(*m_electrons)[0])->E()*1000.,
					                             0);
					myJet->SetLeptonType(1);
				}
				Stop0LUtils::MtCalculator::calcTransverseMass(myJet->Pt(), myJet->Phi(), m_met, m_metPhi, m_dPhiMetLep, m_mtMetLep);
				m_jets->Add(myJet);
				m_jets->Sort();
			}
			else if (chanNum==2 && (int)rI==chanNum){
				if(m_nSigMuons==2){
					// Must be an or since we are using single lepton triggers
					m_passTrigMatch = ((Stop0LUtils::Muon*)(*m_muons)[0])->IsTrigMatched() || ((Stop0LUtils::Muon*)(*m_muons)[1])->IsTrigMatched();
					m_passOSLep = ((Stop0LUtils::Muon*)(*m_muons)[0])->Charge() * ((Stop0LUtils::Muon*)(*m_muons)[1])->Charge() < 0;
					TLorentzVector ll;
					ll += (*(Stop0LUtils::Muon*)(*m_muons)[0]);
					ll += (*(Stop0LUtils::Muon*)(*m_muons)[1]);
					m_llMass = ll.M();
				}
				else{
					m_passTrigMatch = ((Stop0LUtils::Electron*)(*m_electrons)[0])->IsTrigMatched() || ((Stop0LUtils::Electron*)(*m_electrons)[1])->IsTrigMatched();
					m_passOSLep = ((Stop0LUtils::Electron*)(*m_electrons)[0])->Charge() * ((Stop0LUtils::Electron*)(*m_electrons)[1])->Charge() < 0;
					TLorentzVector ll;
					ll += (*(Stop0LUtils::Electron*)(*m_electrons)[0]);
					ll += (*(Stop0LUtils::Electron*)(*m_electrons)[1]);
					m_llMass = ll.M();
				}
			}
			
			// Fill most ntuple variables (some are filled prior to calling this function).
			setVars();

			if (performSkim(chanNum, systI, xAODObjects) && chanNum==(int)rI)
				{
#ifdef DORJIGSAW
					rJigsaw->AnalyzeEvent(*m_jets, m_metX, m_metY);
					//rJigsaw->FillTriggerBits(*m_objTool, *m_event);
#endif
					m_outTrees[rI][systI]->Fill();
				}

			fillCutflow(rI, systI);
			if (m_affectsKinematics) {
				if (m_syst_affectsElectrons) {
					delete xAODObjects.electrons;
					delete xAODObjects.electrons_aux;
				}
				if (m_syst_affectsMuons) {
					delete xAODObjects.muons;
					delete xAODObjects.muons_aux;
				}
				// if(m_syst_affectsTaus) {
				//   delete xAODObjects.taus;
				//   delete xAODObjects.taus_aux;
				// }
				// if(m_syst_affectsPhotons) {
				//   delete xAODObjects.photons;
				//   delete xAODObjects.photons_aux;
				// }
				if (m_syst_affectsJets) {
					delete xAODObjects.jets;
					delete xAODObjects.jets_aux;
				}
			}
			cleanUp();
			m_isNominal = false;
		}
	}
	m_store->clear();

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::postExecute()
{
	// Here you do everything that needs to be done after the main event
	// processing.  This is typically very rare, particularly in user
	// code.  It is mainly used in implementing the NTupleSvc.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::finalize()
{
	// This method is the mirror image of initialize(), meaning it gets
	// called after the last event has been processed on the worker node
	// and allows you to finish up any objects you created in
	// initialize() before they are written to disk.  This is actually
	// fairly rare, since this happens separately for each worker node.
	// Most of the time you want to do your post-processing on the
	// submission node after all your histogram outputs have been
	// merged.  This is different from histFinalize() in that it only
	// gets called on worker nodes that processed input events.

	//Saving the cutflow from 2 leptons processor to the output file
	// for (auto& sysProcessorsSet : twoLeptonProcessors)
	// 	{
	// 		//Get the 2lep processor that connect to 2lep tree in the systematic set.
	// 		auto& twoLeptonProcessor = sysProcessorsSet[2];
	// 		for (auto& hist : twoLeptonProcessor->ReleaseCutFlowHistograms())
	// 			{
	// 				hist.release()->SetDirectory((TDirectory*)m_outputFile);
	// 			}
	// 	}
	if (m_objTool) delete m_objTool;
	if (m_grl) delete m_grl;
	if (m_xSecDB) delete m_xSecDB;

	delete m_mcTruthClassifier;

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode Stop0LRun2::histFinalize()
{
	// This method is the mirror image of histInitialize(), meaning it
	// gets called after the last event has been processed on the worker
	// node and allows you to finish up any objects you created in
	// histInitialize() before they are written to disk.  This is
	// actually fairly rare, since this happens separately for each
	// worker node.  Most of the time you want to do your
	// post-processing on the submission node after all your histogram
	// outputs have been merged.  This is different from finalize() in
	// that it gets called on all worker nodes regardless of whether
	// they processed input events.
	return EL::StatusCode::SUCCESS;
}




//// BEGIN OF CUSTOM METHODS

void Stop0LRun2::initTruthVars() {

	m_metXTruth	 = INIT_VAL;
	m_metYTruth	 = INIT_VAL;
	m_metPhiTruth	 = INIT_VAL;
	m_metTruth	 = INIT_VAL;
	m_metSumEtTruth = INIT_VAL;

	m_isZll = 1;

	m_nTruthJets = 0;
	m_truthJetPt.clear();
	m_truthJetPhi.clear();
	m_truthJetEta.clear();
	m_truthJetM.clear();

	m_trueTopPt = INIT_VAL;
	m_trueTopEta = INIT_VAL;
	m_trueTopPhi = INIT_VAL;
	m_trueTopE = INIT_VAL;
	m_trueTopDaughtPdgId.clear();
	m_trueTopDaughtPt.clear();
	m_trueTopDaughtEta.clear();
	m_trueTopDaughtPhi.clear();
	m_trueTopDaughtE.clear();

	m_trueAntitopPt = INIT_VAL;
	m_trueAntitopEta = INIT_VAL;
	m_trueAntitopPhi = INIT_VAL;
	m_trueAntitopE = INIT_VAL;
	m_trueAntitopDaughtPdgId.clear();
	m_trueAntitopDaughtPt.clear();
	m_trueAntitopDaughtEta.clear();
	m_trueAntitopDaughtPhi.clear();
	m_trueAntitopDaughtE.clear();

}

void Stop0LRun2::initVars() {

	for (auto& weight : m_weights)
		weight.second = 1.0;

	m_lepPt.clear();
	m_lepEta.clear();
	m_lepPhi.clear();
	m_lepE.clear();
	m_lepType.clear();

	m_mtMetLep = INIT_VAL;
	m_dPhiMetLep = INIT_VAL;

	m_nElectrons = 0;
	m_nSigElectrons = 0;
	m_elPt.clear();
	m_elE.clear();
	m_elEta.clear();
	m_elPhi.clear();
	m_elIsSignal.clear();

	m_nMuons = 0;
	m_nBadMuons = 0;
	m_nCosmicMuons = 0;
	m_nSigMuons = 0;
	m_muPt.clear();
	m_muE.clear();
	m_muEta.clear();
	m_muPhi.clear();
	m_muIsSignal.clear();

	m_nJets = 0;
	m_nBJets = 0;
	m_nBadJets = 0;
	m_jetPt.clear();
	m_jetEMFrac.clear();
	m_jetHECFrac.clear();
	m_jetM.clear();
	m_jetEta.clear();
	m_jetPhi.clear();
	m_jetIsSignal.clear();
	m_jetIsBTagged.clear();
	m_jetMV1Weight.clear();
	m_jetMt.clear();
	m_jetDPhiMet.clear();

	m_jetLeadTagIndex = INIT_VAL;
	m_jetSubleadTagIndex = INIT_VAL;
	m_jetDPhiMetMin3 = 0;

	m_topDRMin0Pt = INIT_VAL;
	m_topDRMin0Phi = INIT_VAL;
	m_topDRMin0Eta = INIT_VAL;
	m_topDRMin0M = INIT_VAL;

	m_topDRMin1Pt = INIT_VAL;
	m_topDRMin1Phi = INIT_VAL;
	m_topDRMin1Eta = INIT_VAL;
	m_topDRMin1M = INIT_VAL;

	m_topDRMinB0Pt = INIT_VAL;
	m_topDRMinB0Phi = INIT_VAL;
	m_topDRMinB0Eta = INIT_VAL;
	m_topDRMinB0M = INIT_VAL;
	m_topDRMinB0Const.clear();

	m_topDRMinB1Pt = INIT_VAL;
	m_topDRMinB1Phi = INIT_VAL;
	m_topDRMinB1Eta = INIT_VAL;
	m_topDRMinB1M = INIT_VAL;
	m_topDRMinB1Const.clear();

	m_nAntiKt8 = 0;
	m_antiKt8DijetMassAsym = fabs(INIT_VAL);
	m_antiKt8Pt.clear();
	m_antiKt8Phi.clear();
	m_antiKt8Eta.clear();
	m_antiKt8M.clear();
	m_antiKt8Const.clear();
	m_antiKt8NConst.clear();

	m_nAntiKt12 = 0;
	m_antiKt12DijetMassAsym = fabs(INIT_VAL);
	m_antiKt12Pt.clear();
	m_antiKt12Phi.clear();
	m_antiKt12Eta.clear();
	m_antiKt12M.clear();
	m_antiKt12Const.clear();
	m_antiKt12NConst.clear();

	m_passLepTrigger = false;
	m_passMetTrigger = false;
	m_passPhotonTrigger = false;

	m_lepTriggerPrescale = INIT_VAL;
	m_metTriggerPrescale = INIT_VAL;
	m_photonTriggerPrescale = INIT_VAL;

	m_passTrigMatch = false;
	m_passOSLep = false;
	
	m_metX     = INIT_VAL;
	m_metY     = INIT_VAL;
	m_metPhi   = INIT_VAL;
	m_met	     = INIT_VAL;
	m_metSumEt = INIT_VAL;

	m_metXOrig	 = INIT_VAL;
	m_metYOrig	 = INIT_VAL;
	m_metPhiOrig	 = INIT_VAL;
	m_metOrig	 = INIT_VAL;
	m_metSumEtOrig = INIT_VAL;

	m_metXTrk   = INIT_VAL;
	m_metYTrk   = INIT_VAL;
	m_metPhiTrk = INIT_VAL;
	m_metTrk    = INIT_VAL;
	m_metSumEtTrk = INIT_VAL;

	m_rawMetTrk = INIT_VAL;
	m_rawMetXTrk = INIT_VAL;
	m_rawMetYTrk = INIT_VAL;
	m_rawMetPhiTrk = INIT_VAL;
	m_rawDPhiMetTrackMet = INIT_VAL;

	m_metSig = INIT_VAL;
	m_mEff   = INIT_VAL;
	m_ht	   = 0;
	m_htSig  = INIT_VAL;

	m_primVtxNTrks = 0;

	m_mtBMin	    = INIT_VAL;
	m_mtNonBMin	    = INIT_VAL;
	m_minMt	    = INIT_VAL;
	m_mtTauCand	    = -1.*INIT_VAL;
	m_tauJetNTracks   = INIT_VAL;

	m_mcPt.clear();
	m_mcPhi.clear();
	m_mcEta.clear();
	m_mcE.clear();
	m_mcPdgId.clear();
	m_mcBarcode.clear();
	m_mcChildren.clear();

	m_dRbbJetHt = INIT_VAL;
	m_dRbbJet = INIT_VAL;

	m_mBB = INIT_VAL;
	m_dRbb = INIT_VAL;

	m_dPhiMetTrackMet = -1.*INIT_VAL;

	m_allJetMass = INIT_VAL;
	m_allJetPt = INIT_VAL;
	m_allJetPhi = INIT_VAL;
	m_allJetEta = INIT_VAL;

	m_llMass = INIT_VAL;
	
	m_hasPrimVtx = false;
	m_primVtxNTrks = INIT_VAL;
	m_nvtx = INIT_VAL;
	
	m_electrons = new TObjArray();
	m_electrons->SetOwner(kTRUE);

	m_muons = new TObjArray();
	m_muons->SetOwner(kTRUE);

	m_jets = new TObjArray();
	m_jets->SetOwner(kTRUE);

    mt2_truth = INIT_VAL;
    mt2_drmin = INIT_VAL;
    mt2_rcjet = INIT_VAL;
    
}

bool Stop0LRun2::getTruth() {
	//std::cout << std::endl << "Event number: " << m_eventNumber << std::endl;
	// Get MC truth record and do a little bit of skimming by removing particles with pT < 10

	m_topDecayType = 0;
	m_antitopDecayType = 0;

    TLorentzVector trueMETvector(0., 0., 0., 0.);
    trueMet = INIT_VAL;
    trueXMet = INIT_VAL;
    trueYMet = INIT_VAL;
    trueMetParticle->Reset();


	const xAOD::TruthParticleContainer* truthParts = 0;
	if (m_event->retrieve(truthParts, "TruthParticles").isSuccess() && truthParts->size()>0) {
		CHECK(m_objTool->FindSusyHP(truthParts, m_susyPdgId1, m_susyPdgId2));

		bool isSUSY = false;
		for (unsigned int tI=0; tI< truthParts->size(); ++tI) {
			const xAOD::TruthParticle* tPart = (*truthParts)[tI];

			isSUSY = ((abs(tPart->pdgId())>1000000 && abs(tPart->pdgId())<1000007) || // squarkL
			          (abs(tPart->pdgId())>1000010 && abs(tPart->pdgId())<1000017) || // sleptonL
			          (abs(tPart->pdgId())>2000000 && abs(tPart->pdgId())<2000007) || // squarkR
			          (abs(tPart->pdgId())>2000010 && abs(tPart->pdgId())<2000017) || // sleptonR
			          (abs(tPart->pdgId())>1000020 && abs(tPart->pdgId())<1000040)); // gauginos
			if (isSUSY)
				break;
		}

		if (isSUSY)
			m_finalState = SUSY::finalState(m_susyPdgId1, m_susyPdgId2);

		for (unsigned int tI=0; tI< truthParts->size(); ++tI) {
			const xAOD::TruthParticle* tPart = (*truthParts)[tI];

			if (!xAODtruthUtils::isSelfDecay(tPart) && abs(tPart->pdgId())!=6 && tPart->pt() > 1) {
				m_mcE.push_back(tPart->e());
				m_mcPt.push_back(tPart->pt()/1000.);
				m_mcPhi.push_back(tPart->phi());
				m_mcEta.push_back(tPart->eta());
				m_mcBarcode.push_back(tPart->barcode());
				m_mcPdgId.push_back(tPart->pdgId());
			}

			if (abs(tPart->pdgId())==23 && !xAODtruthUtils::isSelfDecay(tPart) && m_mcChannelNumber >= 410073 && m_mcChannelNumber <= 410075) {
				std::vector< const xAOD::TruthParticle* > descendents;
				xAODtruthUtils::getDescendants(tPart, descendents);
				if (descendents.size() != 2)
					m_isZll = 0;
			}
			if (abs(tPart->pdgId())==6 && !xAODtruthUtils::isSelfDecay(tPart)) {
				std::vector< const xAOD::TruthParticle* > descendents;
				if (tPart->pdgId() == 6) {
					m_trueTopPt = tPart->pt()/1000.;
					m_trueTopEta = tPart->eta();
					m_trueTopPhi = tPart->phi();
					m_trueTopE = tPart->e()/1000.;
				}
				else {
					m_trueAntitopPt = tPart->pt()/1000.;
					m_trueAntitopEta = tPart->eta();
					m_trueAntitopPhi = tPart->phi();
					m_trueAntitopE = tPart->e()/1000.;
				}
				FillTruthWb(tPart);
				xAODtruthUtils::getDescendants(tPart, descendents);
                
                //Truth met vector here:
                
				for (unsigned int dI = 0; dI < descendents.size(); dI++) {
                  if (abs(descendents.at(dI)->pdgId()==1000022)){
                    trueMETvector += descendents.at(dI)->p4();
                    trueMetParticle->push_back(descendents.at(dI)->p4());
                  }  
                  
                  if (descendents.at(dI)->nParents() == 1
					    && abs(descendents.at(dI)->parent(0)->pdgId()) == 24
					    && abs(descendents.at(dI)->pdgId())!=24 &&
					    abs(descendents.at(dI)->pdgId())%2==1) {

                    if (tPart->pdgId() == 6)
							m_topDecayType = abs(descendents.at(dI)->pdgId());
                    else
							m_antitopDecayType = abs(descendents.at(dI)->pdgId());
					}
                  if ((xAODtruthUtils::isFinalParton(descendents.at(dI))) //add truth particle met here?
					    ) {
                  
                    if (tPart->pdgId() == 6) {
							m_trueTopDaughtPdgId.push_back(descendents.at(dI)->pdgId());
							m_trueTopDaughtPt.push_back(descendents.at(dI)->pt()/1000.);
							m_trueTopDaughtEta.push_back(descendents.at(dI)->eta());
							m_trueTopDaughtPhi.push_back(descendents.at(dI)->phi());
							m_trueTopDaughtE.push_back(descendents.at(dI)->e()/1000.);

                            if (abs(descendents.at(dI)->pdgId())==12 || abs(descendents.at(dI)->pdgId())==14 || abs(descendents.at(dI)->pdgId())==16){
							trueMETvector += descendents.at(dI)->p4();
							trueMetParticle->push_back(descendents.at(dI)->p4());
                            }
                            
						}
                    else {
							m_trueAntitopDaughtPdgId.push_back(descendents.at(dI)->pdgId());
							m_trueAntitopDaughtPt.push_back(descendents.at(dI)->pt()/1000.);
							m_trueAntitopDaughtEta.push_back(descendents.at(dI)->eta());
							m_trueAntitopDaughtPhi.push_back(descendents.at(dI)->phi());
							m_trueAntitopDaughtE.push_back(descendents.at(dI)->e()/1000.);

                            if (abs(descendents.at(dI)->pdgId())==12 || abs(descendents.at(dI)->pdgId())==14 || abs(descendents.at(dI)->pdgId())==16){
							trueMETvector += descendents.at(dI)->p4();
							trueMetParticle->push_back(descendents.at(dI)->p4());
                            }
                    }
                  }
				}
                //		}
            //		}

        metXparticle = trueMETvector.Px();
        metYparticle = trueMETvector.Py();
        metparticle = trueMETvector.Pt();
        myversion();
        //float mt2_truth;
        mt2_truth = INIT_VAL;
        
        if (m_trueTopPt>0. && m_trueAntitopPt>0.){
          TLorentzVector visaTruth;
          visaTruth.SetPtEtaPhiE(m_trueTopPt, m_trueTopEta, m_trueTopPhi, m_trueTopE);
          TLorentzVector visbTruth;
          visbTruth.SetPtEtaPhiE(m_trueAntitopPt, m_trueAntitopEta, m_trueAntitopPhi, m_trueAntitopE);
          TLorentzVector metTruthVector;
          metTruthVector.SetPxPyPzE( metXparticle, metYparticle, 0., metparticle);
          mt2_truth = ComputeMT2(visaTruth,visbTruth,metTruthVector,0.,0.).Compute();
        }
        else{
          mt2_truth = -99.;
        }
        //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for Truth is: " << mt2_truth << std::endl;
            }}
    }	
	else {
		return false;
	}
    
	//FillTruthHt(truthParts);
	FillTruthNeutrinoMET(truthParts);
    
    // calculateMt2Truth(m_trueTopPt, m_trueTopEta, m_trueTopPhi, m_trueTopE, m_trueAntitopPt, m_trueAntitopEta, m_trueAntitopPhi, m_trueAntitopE, m_metTruth, m_metXTruth, m_metYTruth);

    // calculateMt2DRmin(m_topDRMin0Pt, m_topDRMin0Eta, m_topDRMin0Phi, m_topDRMin0M, m_topDRMin1Pt, m_topDRMin1Eta, m_topDRMin1Phi, m_topDRMin1M, m_met, m_metX, m_metY);

    // if (m_antiKt12Pt.size()>2) {
    // calculateMt2RCjet(m_antiKt12Pt[0], m_antiKt12Eta[0], m_antiKt12Phi[0], m_antiKt12M[0], m_antiKt12Pt[1], m_antiKt12Eta[1], m_antiKt12Phi[1], m_antiKt12M[1], m_met, m_metX, m_metY);
    //   } 

    
	// Filled all the truth info
	const xAOD::JetContainer* truthJets = 0;
	if (m_event->retrieve(truthJets, "AntiKt4TruthJets").isSuccess()) {
		for (unsigned int jetI=0; jetI<truthJets->size(); ++jetI) {
			m_nTruthJets = truthJets->size();
			const xAOD::Jet* jet = (*truthJets)[jetI];
			m_truthJetPt.push_back(jet->pt());
			m_truthJetEta.push_back(jet->eta());
			m_truthJetPhi.push_back(jet->phi());
			m_truthJetM.push_back(jet->m());
		}
	}
	else {
		return false;
	}

	const xAOD::MissingETContainer* truthMetContainer = 0;
	if (m_event->retrieve(truthMetContainer, "MET_Truth").isSuccess()) {
		xAOD::MissingETContainer::const_iterator it = truthMetContainer->find("NonInt");
		if (it == truthMetContainer->end()) {
			::Error(APP_NAME, "No RefFinal inside MET container");
			return false;
		}
		m_metXTruth	   = (*it)->mpx()/1000.;
		m_metYTruth	   = (*it)->mpy()/1000.;
		m_metPhiTruth   = (*it)->phi();
		m_metSumEtTruth = (*it)->sumet()/1000.;
		m_metTruth = TMath::Sqrt(m_metXTruth*m_metXTruth + m_metYTruth*m_metYTruth);
	}
	else {
		return false;
	}

	return true;
}

void Stop0LRun2::FillTruthHt(const xAOD::TruthParticleContainer* truthParts)
{
	truthHt = 0;
	const xAOD::JetContainer* antiKt4TruthWZJetsContainer = nullptr;
	CHECK(m_event->retrieve(antiKt4TruthWZJetsContainer, "AntiKt4TruthWZJets"));
	for (const auto& jet : *antiKt4TruthWZJetsContainer)
		{
			if (jet->pt() > 35000 && abs(jet->eta()) < 2.5)
				{
					truthHt += jet->pt()/1000;
				}
		}

	for (const auto& particle : *truthParts)
		{
			if ((particle->isElectron() || particle->isMuon())
			    && (particle->pt() > 25000) && (particle->abseta() < 2.5)
			    && (isFromTopW(particle))
			    )
				{
					truthHt += particle->pt()/1000;
				}
		}

}

void Stop0LRun2::FillTruthNeutrinoMET(const xAOD::TruthParticleContainer* truthParts)
{
	truthNeutrinoMet = INIT_VAL;
	truthNeutrinos->Reset();
	TLorentzVector neutrinoMET(0., 0., 0., 0.);
	for (const auto& particle: *truthParts)
		{
			if (particle->isNeutrino())
				{
					if (isFromTopW(particle))
						{
							neutrinoMET += particle->p4();
							truthNeutrinos->push_back(particle->p4());
						}
				}
		}
	truthNeutrinoMet = neutrinoMET.Pt()/1000;
}

void Stop0LRun2::FillTruthWb(const xAOD::TruthParticle * particle)
{
	truthWs->Reset();
	truthBottoms->Reset();
	for (size_t childIndex = 0; childIndex < particle->nChildren(); childIndex++)
		{
			auto child = particle->child(childIndex);
			if (child->isW())
				{
					truthWs->push_back_mev_to_gev(child->p4());
				}
			else if (abs(child->pdgId())==5)
				{
					truthBottoms->push_back_mev_to_gev(child->p4());
				}
		}
}

bool Stop0LRun2::isFromTopW(const xAOD::TruthParticle* particle)
{
	//If we need to classify ourselves
	if (particle->isAvailable<unsigned int>("classifierParticleOrigin"))
		{
			auto origin = particle->auxdata< unsigned int >("classifierParticleOrigin");
			return (origin == MCTruthPartClassifier::WBoson || origin == MCTruthPartClassifier::top);
		}
	else
		{
			auto result = m_mcTruthClassifier->particleTruthClassifier(particle);
			//Second in the pair is particle origin
			auto& type = result.first;
			auto& origin = result.second;
			if (origin == MCTruthPartClassifier::WBoson || origin == MCTruthPartClassifier::top)
				{
					return true;
				}
			else if (origin == MCTruthPartClassifier::NonDefined)
				{
					//MCTruthPartClassifier::ParticleDef particleDef;
					//std::cout << "Event " << m_eventNumber << " MCTruthClassifier failed to classify " << particleDef.sParticleType.at(result.first) << " PdgID: " << particle->pdgId() << std::endl;
					//Failed to classify. Most likely the link to the top is broken. Give it a second try by look at the ancestors
					std::vector<const xAOD::TruthParticle*> ancestors;
					xAODtruthUtils::getAncestors(particle, ancestors);
					if (type == MCTruthPartClassifier::Neutrino)
						{
							if (ancestors[0]->isW())
								{
									return true;
								}
							else if (ancestors[0]->isTau())
								{
									auto tauresult = m_mcTruthClassifier->particleTruthClassifier(ancestors[0]);
									if (tauresult.second == MCTruthPartClassifier::WBoson)
										{
											return true;
										}
									else
										{
											return false;
										}
								}
						}
					return false;
				}
			else
				{
					return false;
				}
		}
}

// Fix mem leaks created by returns without cleanups
Stop0LRun2::XAODObjects Stop0LRun2::createObj(std::string systName)
{
	m_passMetTrigger = m_objTool->IsTrigPassed("HLT_xe70_tc_lcw");
	m_metTriggerPrescale = m_objTool->GetTrigPrescale("HLT_xe70_tc_lcw");

	std::string e24TrigName = "";
	if (!m_isData)
		e24TrigName = "HLT_e24_lhmedium_L1EM18VH";
	else
		e24TrigName = "HLT_e24_lhmedium_L1EM20VH";
	
	std::vector<std::string> elTrigNames = {e24TrigName, "HLT_e60_lhmedium", "HLT_e120_lhloose"};
	std::vector<std::string> passElTrigNames;
	
	bool passElTrigger = false;
	for (auto & trigName: elTrigNames){
		passElTrigger = passElTrigger || m_objTool->IsTrigPassed(trigName);
		if(m_objTool->IsTrigPassed(trigName))
			passElTrigNames.push_back(trigName);
	}
	
	std::vector<std::string> muTrigNames = {"HLT_mu20_iloose_L1MU15", "HLT_mu50"};
	std::vector<std::string> passMuTrigNames;
	
	bool passMuTrigger = false;
	for (auto & trigName: muTrigNames){
		passMuTrigger = passMuTrigger || m_objTool->IsTrigPassed(trigName);
		if(m_objTool->IsTrigPassed(trigName))
			passMuTrigNames.push_back(trigName);
	}
	
	m_passLepTrigger = passMuTrigger || passElTrigger;

	if (m_eventNumber == m_evtToPrint) {
		std::cout << "HLT_xe70_tc_lcw is " << m_objTool->IsTrigPassed("HLT_xe70_tc_lcw") << std::endl << std::endl;
		for (auto & trigName: elTrigNames){
			std::cout << trigName << " is " << m_objTool->IsTrigPassed(trigName) << std::endl;
		}
		std::cout << std::endl;
		
		for (auto & trigName: muTrigNames){
			std::cout << trigName << " is " << m_objTool->IsTrigPassed(trigName) << std::endl;
		}
		std::cout << std::endl;
	}
	
	XAODObjects xAODObjects;
	if (m_debug)
		Info(APP_NAME, "Starting with systematic \"%s\" ", systName.c_str());

	// Generic pointers for either nominal or systematics copy
	xAOD::ElectronContainer* electrons(m_electrons_nominal);
	// xAOD::PhotonContainer* photons(m_photons_nominal);
	xAOD::MuonContainer* muons(m_muons_nominal);
	xAOD::JetContainer* jets(m_jets_nominal);
	// xAOD::TauJetContainer* taus(m_taus_nominal);
	xAOD::MissingETContainer* metcst(m_metcst_nominal);
	xAOD::MissingETContainer* mettst(m_mettst_nominal);
	xAOD::MissingETContainer* mettrack(m_mettrack_nominal);
	xAOD::MissingETContainer* metleptst(m_metleptst_nominal);
	// Aux containers too
	xAOD::ShallowAuxContainer* electrons_aux(m_electrons_nominal_aux);
	// xAOD::ShallowAuxContainer* photons_aux(m_photons_nominal_aux);
	xAOD::ShallowAuxContainer* muons_aux(m_muons_nominal_aux);
	xAOD::ShallowAuxContainer* jets_aux(m_jets_nominal_aux);
	// xAOD::ShallowAuxContainer* taus_aux(m_taus_nominal_aux);
	xAOD::MissingETAuxContainer* metcst_aux(m_metcst_nominal_aux);
	xAOD::MissingETAuxContainer* mettst_aux(m_mettst_nominal_aux);
	xAOD::MissingETAuxContainer* mettrack_aux(m_mettrack_nominal_aux);
	xAOD::MissingETAuxContainer* metleptst_aux(m_metleptst_nominal_aux);

	if (m_affectsKinematics) {
		if (m_syst_affectsElectrons) {
			xAOD::ElectronContainer* electrons_syst(0);
			xAOD::ShallowAuxContainer* electrons_syst_aux(0);
			CHECK(m_objTool->GetElectrons(electrons_syst, electrons_syst_aux));
			electrons = electrons_syst;
			electrons_aux = electrons_syst_aux;

			// CHECK(m_store->record(electrons_syst, "Electrons_"+systName));
			// CHECK(m_store->record(electrons_syst_aux, "Electrons"+systName+"Aux."));
		}

		if (m_syst_affectsMuons) {
			xAOD::MuonContainer* muons_syst(0);
			xAOD::ShallowAuxContainer* muons_syst_aux(0);
			CHECK(m_objTool->GetMuons(muons_syst, muons_syst_aux));
			muons = muons_syst;
			muons_aux = muons_syst_aux;
			// CHECK(m_store->record(muons_syst, "Muons_"+systName));
			// CHECK(m_store->record(muons_syst_aux, "Muons"+systName+"Aux."));
		}

		// if(m_syst_affectsTaus) {
		//   xAOD::TauJetContainer* taus_syst(0);
		//   xAOD::ShallowAuxContainer* taus_syst_aux(0);
		//   CHECK( m_objTool->GetTaus(taus_syst,taus_syst_aux) );
		//   taus = taus_syst;
		//   taus_aux = taus_syst_aux;
		//  CHECK(m_store->record(taus_syst, "Taus_"+systName));
		//	CHECK(m_store->record(taus_syst_aux, "Taus+systName+"Aux."));
		// }

		// if(m_syst_affectsPhotons) {
		//   xAOD::PhotonContainer* photons_syst(0);
		//   xAOD::ShallowAuxContainer* photons_syst_aux(0);
		//   CHECK( m_objTool->GetPhotons(photons_syst,photons_syst_aux) );
		//   photons = photons_syst;
		//   photons_aux = photons_syst_aux;
		// CHECK(m_store->record(photons_syst, "Photons_"+systName));
		// 	CHECK(m_store->record(photons_syst_aux, "Photons"+systName+"Aux."));
		// }

		if (m_syst_affectsJets) {
			xAOD::JetContainer* jets_syst(0);
			xAOD::ShallowAuxContainer* jets_syst_aux(0);
			CHECK(m_objTool->GetJetsSyst(*m_jets_nominal, jets_syst, jets_syst_aux));
			jets = jets_syst;
			jets_aux = jets_syst_aux;
			// CHECK(m_store->record(jets_syst, "Jets_"+systName));
			// CHECK(m_store->record(jets_syst_aux, "Jets"+systName+"Aux."));
		}

		xAOD::MissingETContainer* metcst_syst = new xAOD::MissingETContainer;
		xAOD::MissingETAuxContainer* metcst_syst_aux = new xAOD::MissingETAuxContainer;
		xAOD::MissingETContainer* mettst_syst = new xAOD::MissingETContainer;
		xAOD::MissingETAuxContainer* mettst_syst_aux = new xAOD::MissingETAuxContainer;
		xAOD::MissingETContainer* mettrack_syst = new xAOD::MissingETContainer;
		xAOD::MissingETAuxContainer* mettrack_syst_aux = new xAOD::MissingETAuxContainer;
		xAOD::MissingETContainer* metleptst_syst = new xAOD::MissingETContainer;
		xAOD::MissingETAuxContainer* metleptst_syst_aux = new xAOD::MissingETAuxContainer;

		metcst_syst->setStore(metcst_syst_aux);
		mettst_syst->setStore(mettst_syst_aux);
		mettrack_syst->setStore(mettrack_syst_aux);
		metleptst_syst->setStore(metleptst_syst_aux);

		metcst_syst->reserve(10);
		mettst_syst->reserve(10);
		mettrack_syst->reserve(10);
		metleptst_syst->reserve(10);

		metcst = metcst_syst;
		mettst = mettst_syst;
		mettrack = mettrack_syst;
		metleptst = metleptst_syst;

		metcst_aux = metcst_syst_aux;
		mettst_aux = mettst_syst_aux;
		mettrack_aux = mettrack_syst_aux;
		metleptst_aux = metleptst_syst_aux;

		// CHECK(m_store->record(metcst_syst, "metcst_"+systName));
		// CHECK(m_store->record(mettst_syst, "mettst_"+systName));

		// CHECK(m_store->record(metcst_syst_aux, "metcst"+systName+"Aux."));
		// CHECK(m_store->record(mettst_syst_aux, "mettst"+systName+"Aux."));
	}
	xAOD::IParticleContainer* signalLeptons = new xAOD::IParticleContainer(SG::OwnershipPolicy::VIEW_ELEMENTS);

	//------------
	// OVERLAP REMOVAL
	//------------
	if (m_isNominal || (m_affectsKinematics && (m_syst_affectsElectrons || m_syst_affectsMuons || m_syst_affectsJets))) {
		CHECK(m_objTool->OverlapRemoval(electrons, muons, jets));
		if(m_doSpecialOR)
			Stop0LUtils::OverlapRemoval(electrons, muons, jets, m_dRejet, m_dRmujet, m_dRemu, 0);
	}

	// Fill member electrons that will be stored in the ntuple
	m_electrons->Clear();
	m_muons->Clear();
	m_jets->Clear();

	for (const auto& el : *electrons) {
		if (!el->auxdata< char >("baseline")) continue;
		// Do stuff with baseline electrons without overlap removal requirement below

		if (!(el->auxdata< char >("passOR"))) continue;
		// This does stuff with baseline and overlap removal requirement
		Stop0LUtils::Electron* myElectron = new Stop0LUtils::Electron(el->pt(), el->eta(), el->phi(), el->e(), el->charge());
		if (el->auxdata< char >("signal")){
			m_nSigElectrons++;
			myElectron->SetSignal(true);
			signalLeptons->push_back(el);
			for (auto & trigName: passElTrigNames){
				if(m_objTool->IsTrigMatched(el,trigName)){
					myElectron->SetTrigMatched(true);
					break;
				}
			}
		}
		m_electrons->Add(myElectron);
	}
	

	if (m_isNominal || m_syst_affectsElectrons) {
		m_elecSF = m_objTool->GetTotalElectronSF(*electrons);
	}
	if (m_isNominal) { m_elecSF_nominal = m_elecSF; }
	else if (!m_syst_affectsElectrons) { m_elecSF = m_elecSF_nominal; }

	for (const auto& mu : *muons) {
		if (m_eventNumber == m_evtToPrint) {
			std::cout << "Muons info before everything:" << std::endl;
			std::cout << "(pt,eta,phi,e): (" << mu->pt() << "," << mu->eta() << "," << mu->phi() << "," << mu->e() << ")" << std::endl;
			std::cout << "passBaseline, passOR, isCosmic, isBadMuon: "
			          << (mu->auxdata< char >("baseline")==1) << ","
			          << (mu->auxdata< char >("passOR")==1) << ","
			          << (mu->auxdata< char >("cosmic")==1) << ","
			          << (mu->auxdata< char >("bad")==1) << std::endl << std::endl;
		}
		if (!mu->auxdata< char >("baseline")) continue;

		if (mu->auxdata< char >("bad"))
			m_nBadMuons++;
		// Do stuff with baseline muons without overlap removal requirement below

		if (!mu->auxdata< char >("passOR")) continue;
		// This does stuff with baseline and overlap removal requirement

		// Remove events with cosmics only AFTER overlap removal
		if (mu->auxdata< char >("cosmic")) {
			m_nCosmicMuons++;
			continue;
		}

		Stop0LUtils::Muon* myMuon = new Stop0LUtils::Muon(mu->pt(), mu->eta(), mu->phi(), mu->m(),
		                                                  mu->charge(), mu->auxdata< char >("cosmic"),
		                                                  mu->auxdata< char >("bad"),
		                                                  mu->muonType() == xAOD::Muon::Combined);
		if (mu->auxdata< char >("signal")){
			m_nSigMuons++;
			signalLeptons->push_back(mu);
			myMuon->SetSignal(true);
			for (auto & trigName: passMuTrigNames){
				if(m_objTool->IsTrigMatched(mu,trigName)){
					myMuon->SetTrigMatched(true);
					break;
				}
			}
		}
		m_muons->Add(myMuon);
	}

	if (m_isNominal || m_syst_affectsMuons) {
		m_muonSF = m_objTool->GetTotalMuonSF(*muons, true, true);
	}

	if (m_isNominal) { m_muonSF_nominal = m_muonSF; }
	else if (!m_syst_affectsMuons) { m_muonSF = m_muonSF_nominal; }

	double weight_mv2c20 = INIT_VAL;
	auto& ntupJetExtra = ntupJetExtraMap->at(nominalJetCollection);
	ntupJetExtra->Reset();
	static SG::AuxElement::ConstAccessor<float> getJVT("Jvt");

	for (const auto& jet : *jets) {
		
		if (!jet->auxdata< char >("baseline")) continue;

		if (!jet->auxdata< char >("passOR")) continue;

		// Check for bad jets AFTER overlap removal
		if (jet->auxdata< char >("bad")) {
			m_nBadJets++;
			continue;
		}

		if (!jet->auxdata< char >("signal")) continue;

		std::vector<int> nTrkVec;
		Stop0LUtils::Jet* myJet = new Stop0LUtils::Jet(jet->pt(), jet->eta(), jet->phi(), jet->e(), 0);
		float emfrac = INIT_VAL, hecfrac = INIT_VAL;
		jet->getAttribute(xAOD::JetAttribute::EMFrac, emfrac);
		jet->getAttribute(xAOD::JetAttribute::HECFrac, hecfrac);
		myJet->SetQualityVariables(0, 0, emfrac, 0, hecfrac, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

		myJet->SetSignal(true);
		jet->getAttribute(xAOD::JetAttribute::NumTrkPt500, nTrkVec);
		myJet->SetNTrk(nTrkVec[0]);
		m_jets->Add(myJet);

		myJet->SetBTagged(jet->auxdata< char >("bjet")==1);
		jet->btagging()->MVx_discriminant("MV2c20", weight_mv2c20);
		myJet->SetTagMV1(weight_mv2c20);

		ntupJetExtra->JVT.push_back(getJVT(*jet));
		//(Truth) flavour label
		//https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTagCalib2015#Flavour_Labeling
		ntupJetExtra->FlavourLabel.push_back(xAOD::jetFlavourLabel(jet));
	}

	if (!m_isData) {
		if (m_isNominal || m_syst_affectsBTag || (m_affectsKinematics && m_syst_affectsJets)) {
			m_objTool->BtagSF(jets);
	
			for (const auto& jet : *jets) {
				if (jet->auxdata< char >("baseline") && jet->auxdata< char >("passOR")
				    && !jet->auxdata< char >("bad") && jet->auxdata< char >("signal")) {
					m_bTagSF *= jet->auxdata< double >("effscalefact");
				}
			}
		}
		if(m_isNominal)
			m_bTagSF_nominal = m_bTagSF;
		else if (!(m_syst_affectsBTag || (m_affectsKinematics && m_syst_affectsJets)))
			m_bTagSF = m_bTagSF_nominal;
	}
	
	// Now sort by pT
	m_electrons->Sort();
	m_muons->Sort();
	m_jets->Sort();

	m_nMuons     = m_muons->GetEntries();
	m_nElectrons = m_electrons->GetEntries();

	// MET (re)calculation
	if (m_eventNumber == m_evtToPrint) {
		for (const auto& el : *electrons)
			std::cout << "Electron pT: " << el->pt()<< " Electron eta " <<  el->eta()<< " Electron phi " <<  el->phi() << std::endl;

		for (const auto& mu : *muons)
			std::cout << "Muon pT: " << mu->pt()<< " Muon eta " <<  mu->eta()<< " Muon phi " <<  mu->phi() << std::endl;

		for (const auto& jet : *jets)
			std::cout << "Jet pT: " << jet->pt()<< " Jet eta " <<  jet->eta()<< " Jet phi " <<  jet->phi() << std::endl;
	}

	// ---------
	// MET
	// ---------
	if (m_isNominal || m_affectsKinematics) {
		CHECK(m_objTool->GetMET(*metcst,
		                        jets,
		                        electrons,
		                        muons,
		                        0, // photons,
		                        0, // taus
		                        false));
		CHECK(m_objTool->GetMET(*mettst,
		                        jets,
		                        electrons,
		                        muons,
		                        0, // photons,
		                        0, // taus,
		                        true));
		
		CHECK(m_objTool->GetTrackMET(*mettrack,
		                             jets,
		                             electrons,
		                             muons));

		CHECK(m_objTool->GetMET(*metleptst,
		                             jets,
		                             electrons,
		                             muons,
		                             0,
		                             0,
		                             true,
		                             true,
		                             signalLeptons));
		
	}

	if (m_debug) Info(APP_NAME, "Run/evt: %i/%i. Finished calling GetMET ", m_runNumber, m_eventNumber);
	xAOD::MissingETContainer::const_iterator metIt = mettst->find("Final");
	if (metIt == mettst->end()) {
		::Error(APP_NAME, "Could not find MET final in container");
	}
	m_metX     = (*metIt)->mpx()/1000.;
	m_metY     = (*metIt)->mpy()/1000.;
	m_metPhi   = (*metIt)->phi();
	m_met	     = TMath::Sqrt(m_metX*m_metX + m_metY*m_metY);	// in GeV
	m_metSumEt = (*metIt)->sumet()/1000.;

	xAOD::MissingETContainer::const_iterator metTrkIt = mettrack->find("Track");
	if (metTrkIt == mettrack->end()) {
		::Error(APP_NAME, "Could not find MET track in container");
	}
	m_metXTrk     = (*metTrkIt)->mpx()/1000.;
	m_metYTrk     = (*metTrkIt)->mpy()/1000.;
	m_metPhiTrk   = (*metTrkIt)->phi();
	m_metTrk	     = TMath::Sqrt(m_metXTrk*m_metXTrk + m_metYTrk*m_metYTrk);	// in GeV
	m_metSumEtTrk = (*metTrkIt)->sumet()/1000.;
	
	if (m_metPhiTrk != INIT_VAL && m_metPhi != INIT_VAL)
		m_dPhiMetTrackMet = fabs(TVector2::Phi_mpi_pi(m_metPhiTrk-m_metPhi));

	xAOD::MissingETContainer::const_iterator metLepIt = metleptst->find("Final");
	if (metLepIt == mettst->end()) {
		::Error(APP_NAME, "Could not find MET final in container");
	}
	m_metLepX     = (*metLepIt)->mpx()/1000.;
	m_metLepY     = (*metLepIt)->mpy()/1000.;
	m_metLepPhi   = (*metLepIt)->phi();
	m_metLep	     = TMath::Sqrt(m_metLepX*m_metLepX + m_metLepY*m_metLepY);	// in GeV
	m_metLepSumEt = (*metLepIt)->sumet()/1000.;

	if (m_eventNumber == m_evtToPrint) {
		xAOD::MissingETContainer::const_iterator metItMuon = mettst->find("Muons");

		float metX = (*metItMuon)->mpx()/1000.;
		float metY = (*metItMuon)->mpy()/1000.;
		float met =  TMath::Sqrt(metX*metX + metY*metY);
		std::cout << "Muon met: " << met << std::endl;

		xAOD::MissingETContainer::const_iterator metItJet = mettst->find("RefJet");
		metX = (*metItJet)->mpx()/1000.;
		metY = (*metItJet)->mpy()/1000.;
		met =  TMath::Sqrt(metX*metX + metY*metY);
		std::cout << "Jet met: " << met << std::endl;

		xAOD::MissingETContainer::const_iterator metItEl = mettst->find("RefEle");
		metX = (*metItEl)->mpx()/1000.;
		metY = (*metItEl)->mpy()/1000.;
		met =  TMath::Sqrt(metX*metX + metY*metY);
		std::cout << "Electron met: " << met << std::endl;

		// xAOD::MissingETContainer::const_iterator metItGamma = mettst->find("RefGamma");
		// metX = (*metItGamma)->mpx()/1000.;
		// metY = (*metItGamma)->mpy()/1000.;
		// met =  TMath::Sqrt(metX*metX + metY*metY);
		// std::cout << "Gamma met: " << met << std::endl;

		// xAOD::MissingETContainer::const_iterator metItTau = mettst->find("RefTau");
		// metX = (*metItTau)->mpx()/1000.;
		// metY = (*metItTau)->mpy()/1000.;
		// met =  TMath::Sqrt(metX*metX + metY*metY);
		// std::cout << "Tau met: " << met << std::endl;

	}
	// Get the original MET
	std::string metRefFinalStr = "MET_Reference_AntiKt4EMTopo";
	std::string metTrackStr = "MET_Track";

	const xAOD::MissingETContainer* metContainer = 0;
	if (m_event->retrieve(metContainer, metRefFinalStr.c_str()).isSuccess()) {
		xAOD::MissingETContainer::const_iterator it = metContainer->find("FinalClus");
		if (it == metContainer->end()) {
			::Error(APP_NAME, "No RefFinal inside MET container");
			throw std::runtime_error("No RefFinal inside MET container");
		}
		m_metXOrig	   = (*it)->mpx()/1000.;
		m_metYOrig	   = (*it)->mpy()/1000.;
		m_metPhiOrig   = (*it)->phi();
		m_metSumEtOrig = (*it)->sumet()/1000.;
	}
	else {
		::Error(APP_NAME, "Could not retrieve MET_RefFinal");
		throw std::runtime_error("Could not retrieve MET_RefFinal");
	}

	m_metOrig = TMath::Sqrt(m_metXOrig*m_metXOrig + m_metYOrig*m_metYOrig);

	// Track MET
	const xAOD::MissingETContainer* mettrkcontainer = 0;
	if (m_event->retrieve(mettrkcontainer, metTrackStr.c_str()).isSuccess()) {
		xAOD::MissingETContainer::const_iterator it = mettrkcontainer->find("PVTrack_vx0");
		if (it == mettrkcontainer->end())
			{
				//No vx0
				m_rawMetXTrk = INIT_VAL;
				m_rawMetYTrk = INIT_VAL;
				m_rawMetPhiTrk = INIT_VAL;
				m_rawMetTrk = INIT_VAL;
			}
		else
			{
				m_rawMetXTrk = (*it)->mpx()/1000.;
				m_rawMetYTrk = (*it)->mpy()/1000.;
				m_rawMetPhiTrk = (*it)->phi();
				m_rawMetTrk = TMath::Sqrt(m_rawMetXTrk*m_rawMetXTrk + m_rawMetYTrk*m_rawMetYTrk);
			}
	}
	else {
		::Error(APP_NAME, "Could not retrieve MET_Track");
		throw std::runtime_error("Could not retrieve MET_Track");
	}

	if (m_rawMetPhiTrk != INIT_VAL && m_metPhi != INIT_VAL)
		m_rawDPhiMetTrackMet = fabs(TVector2::Phi_mpi_pi(m_rawMetPhiTrk-m_metPhi));

	// Now get primary verteces
	const xAOD::Vertex* primVx = 0;
	primVx = m_objTool->GetPrimVtx();
	if(primVx){
		m_hasPrimVtx = true;
		m_primVtxNTrks = primVx->nTrackParticles();
	}
	
	const xAOD::VertexContainer* vertices(0);
	if ( m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
		m_nvtx = vertices->size();
	}
	
	xAODObjects.electrons = electrons;
	xAODObjects.electrons_aux = electrons_aux;
	xAODObjects.muons = muons;
	xAODObjects.muons_aux = muons_aux;
	xAODObjects.jets = jets;
	xAODObjects.jets_aux = jets_aux;
	if (m_affectsKinematics) {
		delete metcst;
		delete metcst_aux;
		delete mettst;
		delete mettst_aux;
		delete mettrack;
		delete mettrack_aux;
		delete metleptst;
		delete metleptst_aux;
	}
	return xAODObjects;
}


void Stop0LRun2::setVars() {
	
	// Fill signal electrons, lepType = 0;
	for (int elI=0; elI<m_nElectrons; elI++) {
		Stop0LUtils::Electron*el=(Stop0LUtils::Electron*)(*m_electrons)[elI];
		m_lepPt.push_back(el->Pt());
		m_lepEta.push_back(el->Eta());
		m_lepPhi.push_back(el->Phi());
		m_lepE.push_back(el->E());
		m_lepType.push_back(0);
	}
	
	// Fill signal muons, lepType = 1;
	for (int muI=0; muI<m_nMuons; muI++) {
		Stop0LUtils::Muon*mu=(Stop0LUtils::Muon*)(*m_muons)[muI];
		m_lepPt.push_back(mu->Pt());
		m_lepEta.push_back(mu->Eta());
		m_lepPhi.push_back(mu->Phi());
		m_lepE.push_back(mu->E());
		m_lepType.push_back(1);
	}

	m_ht = 0;
	m_nJets = m_jets->GetEntries();

	TLorentzVector allJets;

	for (int jetI=0; jetI<m_nJets; jetI++) {
		Stop0LUtils::Jet*jet=(Stop0LUtils::Jet*)(*m_jets)[jetI];

		m_jetPt.push_back(jet->Pt());
		m_jetEMFrac.push_back(jet->GetEmfrac());
		m_jetHECFrac.push_back(jet->GetHecf());
		m_jetEta.push_back(jet->Eta());
		m_jetPhi.push_back(jet->Phi());
		m_jetE.push_back(jet->E());
		m_jetM.push_back(jet->M());
		m_jetMV1Weight.push_back(jet->TagMV1());

		if (jet->GetBTagged()) {
			m_nBJets++;

			if (m_jetSubleadTagIndex == INIT_VAL && m_jetLeadTagIndex != INIT_VAL)
				m_jetSubleadTagIndex = jetI;
			if (m_jetLeadTagIndex == INIT_VAL)
				m_jetLeadTagIndex = jetI;
		}

		m_jetIsBTagged.push_back(jet->GetBTagged());
		m_jetIsSignal.push_back(jet->GetSignal());
		m_ht += jet->Pt();
		allJets+=(*jet);
	}

	m_mEff   = m_ht+m_met;
	m_htSig  = m_met/TMath::Sqrt(m_ht);
	m_metSig = 2*m_met/TMath::Sqrt(m_metSumEt);

	if (m_nJets > 0) {
		m_allJetMass = allJets.M();
		m_allJetPt	 = allJets.Pt();
		m_allJetEta	 = allJets.Eta();
		m_allJetPhi	 = allJets.Phi();
	}


	// Calculate mT's and dPhis
	if (m_nJets > 0) {
		Stop0LUtils::MtCalculator mtCalc(m_jets, m_metX, m_metY, m_maxBJetEta);
		m_jetMt          = mtCalc.getAllMts();
		m_jetDPhiMetMin3 = mtCalc.getDPhiMetJetMin3();
		m_jetDPhiMet     = mtCalc.getAllDPhis();
		m_mtBMin         = mtCalc.getMtBMin();
		m_mtNonBMin	     = mtCalc.getMtNonBMin();
		m_minMt          = mtCalc.getMinMt();
		m_mtTauCand	     = mtCalc.getMtTauCand();
		m_tauJetNTracks  = mtCalc.getTauJetNTracks();
	}
	
	// Start top reco
	Stop0LUtils::TopReco topReco(m_jets);

	// Get old DRmin. Maybe add constituent indexes.
	TLorentzVector top0DRMin;
	TLorentzVector top1DRMin;
	int nDRMinTops;
	topReco.GetDRMinTops(top0DRMin, top1DRMin, nDRMinTops);
	if (nDRMinTops > 0) {
		m_topDRMin0Pt  = top0DRMin.Pt();
		m_topDRMin0Eta = top0DRMin.Eta();
		m_topDRMin0Phi = top0DRMin.Phi();
		m_topDRMin0M   = top0DRMin.M();
	}
	if (nDRMinTops > 1) {
		m_topDRMin1Pt  = top1DRMin.Pt();
		m_topDRMin1Eta = top1DRMin.Eta();
		m_topDRMin1Phi = top1DRMin.Phi();
		m_topDRMin1M   = top1DRMin.M();
	}

    if (nDRMinTops > 1) {
      myversion();
      //float mt2_drmin;
      mt2_drmin = INIT_VAL;
      TLorentzVector visaDRmin;
      visaDRmin.SetPtEtaPhiM(m_topDRMin0Pt, m_topDRMin0Eta, m_topDRMin0Phi, m_topDRMin0M);
      TLorentzVector visbDRmin;
      visbDRmin.SetPtEtaPhiM(m_topDRMin1Pt, m_topDRMin1Eta, m_topDRMin1Phi, m_topDRMin1M);
      TLorentzVector metDRminVector;
      metDRminVector.SetPxPyPzE( m_metX, m_metY, 0., m_met);
      mt2_drmin = ComputeMT2(visaDRmin,visbDRmin,metDRminVector,0.,0.).Compute();
    }
    else{
      mt2_drmin = -99.;
    }
    //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for DRMin is: " << mt2_drmin << std::endl;

	// B jet seeded DRmin
	Stop0LUtils::Jet top0DRMinB(0, 0, 0, 0, 0);
	Stop0LUtils::Jet top1DRMinB(0, 0, 0, 0, 0);
	Stop0LUtils::Jet w0DRMinB(0, 0, 0, 0, 0);
	Stop0LUtils::Jet w1DRMinB(0, 0, 0, 0, 0);
	int nDRMinBTops;
	int b1Index = -1;
	int b2Index = -1;
	topReco.GetDRMinBTops(top0DRMinB, top1DRMinB, w0DRMinB, w1DRMinB,
	                      nDRMinBTops, b1Index, b2Index, m_btagCut, m_maxBJetEta);
	if (nDRMinBTops > 0) {
		m_topDRMinB0Pt  = top0DRMinB.Pt();
		m_topDRMinB0Eta = top0DRMinB.Eta();
		m_topDRMinB0Phi = top0DRMinB.Phi();
		m_topDRMinB0M   = top0DRMinB.M();
		m_topDRMinB0Const = top0DRMinB.GetConstIndexes();
	}
	if (nDRMinBTops > 1) {
		m_topDRMinB1Pt  = top1DRMinB.Pt();
		m_topDRMinB1Eta = top1DRMinB.Eta();
		m_topDRMinB1Phi = top1DRMinB.Phi();
		m_topDRMinB1M   = top1DRMinB.M();
		m_topDRMinB1Const = top1DRMinB.GetConstIndexes();
	}


	if (b1Index>=0 && b2Index>=0) {
		TLorentzVector bbVec = (*(Stop0LUtils::Jet*)(*m_jets)[b1Index])+(*(Stop0LUtils::Jet*)(*m_jets)[b2Index]);
		std::vector<int>skipIndexes ={b1Index, b2Index};
		std::vector<std::pair<int, float> > dRBBJets = getOrderedJetIndexDR(bbVec, skipIndexes);
		if (dRBBJets.size()>0) {
			m_dRbbJetHt = dRBBJets.at(0).second/m_ht;
			m_dRbbJet = dRBBJets.at(0).second;
		}

		m_mBB  = bbVec.M();
		m_dRbb = (*(Stop0LUtils::Jet*)(*m_jets)[b1Index]).DeltaR((*(Stop0LUtils::Jet*)(*m_jets)[b2Index]));
	}
	// Reclusetering
	TObjArray * rcJets = new TObjArray();
	rcJets->SetOwner(kTRUE);
	topReco.GetRCTops(rcJets, 1.2);
	for (int rI = 0; rI<rcJets->GetEntries(); rI++) {
		Stop0LUtils::Jet*jet=(Stop0LUtils::Jet*)(*rcJets)[rI];
		m_nAntiKt12++;
		m_antiKt12Pt.push_back(jet->Pt());
		m_antiKt12Eta.push_back(jet->Eta());
		m_antiKt12Phi.push_back(jet->Phi());
		m_antiKt12M.push_back(jet->M());

		m_antiKt12NConst.push_back(jet->GetNconst());
		m_antiKt12Const.push_back(jet->GetConstIndexes());
	}


    if (m_antiKt12Pt.size()>1){
      myversion();
      //float mt2_rcjet;
      mt2_rcjet = INIT_VAL;
      
      if (m_antiKt12Pt[0]>0. && m_antiKt12Pt[1]>0.){
        TLorentzVector visaRCjet;
        visaRCjet.SetPtEtaPhiM(m_antiKt12Pt[0], m_antiKt12Eta[0], m_antiKt12Phi[0], m_antiKt12M[0]);
        TLorentzVector visbRCjet;
        visbRCjet.SetPtEtaPhiM(m_antiKt12Pt[1], m_antiKt12Eta[1], m_antiKt12Phi[1], m_antiKt12M[1]);
        TLorentzVector metRCjetVector;
        metRCjetVector.SetPxPyPzE(m_metX, m_metY, 0., m_met);
        mt2_rcjet = ComputeMT2(visaRCjet,visbRCjet,metRCjetVector,0.,0.).Compute();
      }
      else{
        mt2_rcjet = -99.;
      }
    }
    else{
      mt2_rcjet = -99.;
    }
    //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for RCJets is: " << mt2_rcjet << std::endl;

                                                                                        
    
	if (m_nAntiKt12>1)
		m_antiKt12DijetMassAsym = fabs(m_antiKt12M.at(0)-m_antiKt12M.at(1))/(m_antiKt12M.at(0)+m_antiKt12M.at(1));
	rcJets->Clear();

	rcJets->SetOwner(kTRUE);
	topReco.GetRCTops(rcJets, 0.8);
	for (int rI = 0; rI<rcJets->GetEntries(); rI++) {
		Stop0LUtils::Jet*jet=(Stop0LUtils::Jet*)(*rcJets)[rI];
		m_nAntiKt8++;
		m_antiKt8Pt.push_back(jet->Pt());
		m_antiKt8Eta.push_back(jet->Eta());
		m_antiKt8Phi.push_back(jet->Phi());
		m_antiKt8M.push_back(jet->M());
		m_antiKt8NConst.push_back(jet->GetNconst());
		m_antiKt8Const.push_back(jet->GetConstIndexes());
	}
	if (m_nAntiKt8>1)
		m_antiKt8DijetMassAsym = fabs(m_antiKt8M.at(0)-m_antiKt8M.at(1))/(m_antiKt8M.at(0)+m_antiKt8M.at(1));
	delete rcJets;
}

//mt2 function:

// void Stop0LRun2::calculateMt2Truth(float topPtTruth, float topEtaTruth, float topPhiTruth, float topETruth, float antitopPtTruth, float antitopEtaTruth, float antitopPhiTruth, float antitopETruth, float metTruth, float metXTruth, float metYTruth) {

//   myversion();
//   float mt2_truth;
//   mt2_truth = INIT_VAL;

//   if (topPtTruth>0. && antitopPtTruth>0.){
//     TLorentzVector visaTruth;
//     visaTruth.SetPtEtaPhiE(topPtTruth, topEtaTruth, topPhiTruth, topETruth);
//     TLorentzVector visbTruth;
//     visbTruth.SetPtEtaPhiE(antitopPtTruth, antitopEtaTruth, antitopPhiTruth, antitopETruth);
//     TLorentzVector metTruthVector;
//     metTruthVector.SetPxPyPzE( metXTruth, metYTruth, 0., metTruth);
//     mt2_truth = ComputeMT2(visaTruth,visbTruth,metTruthVector,0.,0.).Compute();
//   }
//   else{
//     mt2_truth = -99.;
//   }
//     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for Truth is: " << mt2_truth << std::endl;

// }

// void Stop0LRun2::calculateMt2DRmin(float topPtDRmin, float topEtaDRmin, float topPhiDRmin, float topMDRmin, float antitopPtDRmin, float antitopEtaDRmin, float antitopPhiDRmin, float antitopMDRmin, float metDRmin, float metXDRmin, float metYDRmin) {

//   myversion();
//   float mt2_drmin;
//   mt2_drmin = INIT_VAL;

//   if (topPtDRmin>0. && antitopPtDRmin>0.){
//     TLorentzVector visaDRmin;
//     visaDRmin.SetPtEtaPhiM(topPtDRmin, topEtaDRmin, topPhiDRmin, topMDRmin);
//     TLorentzVector visbDRmin;
//     visbDRmin.SetPtEtaPhiM(antitopPtDRmin, antitopEtaDRmin, antitopPhiDRmin, antitopMDRmin);
//     TLorentzVector metDRminVector;
//     metDRminVector.SetPxPyPzE( metXDRmin, metYDRmin, 0., metDRmin);
//     mt2_drmin = ComputeMT2(visaDRmin,visbDRmin,metDRminVector,0.,0.).Compute();
//   }
//   else{
//     mt2_drmin = -99.;
//   }
//   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for DRMin is: " << mt2_drmin << std::endl;
// }


// void Stop0LRun2::calculateMt2RCjet(float topPtRCjet, float topEtaRCjet, float topPhiRCjet, float topMRCjet, float antitopPtRCjet, float antitopEtaRCjet, float antitopPhiRCjet, float antitopMRCjet, float metRCjet, float metXRCjet, float metYRCjet) {
//   myversion();
//   float mt2_rcjet;
//   mt2_rcjet = INIT_VAL;

//   if (m_antiKt12Pt.size()>1){
  
//   if (topPtRCjet>0. && antitopPtRCjet>0.){
//     TLorentzVector visaRCjet;
//     visaRCjet.SetPtEtaPhiM(topPtRCjet, topEtaRCjet, topPhiRCjet, topMRCjet);
//     TLorentzVector visbRCjet;
//     visbRCjet.SetPtEtaPhiM(antitopPtRCjet, antitopEtaRCjet, antitopPhiRCjet, antitopMRCjet);
//     TLorentzVector metRCjetVector;
//     metRCjetVector.SetPxPyPzE( metXRCjet, metYRCjet, 0., metRCjet);
//     mt2_rcjet = ComputeMT2(visaRCjet,visbRCjet,metRCjetVector,0.,0.).Compute();
//   }
//   else{
//     mt2_rcjet = -99.;
//   }
//   }
//   else{
//     mt2_rcjet = -99.;
//   }
//     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MT2 for RCJets is: " << mt2_rcjet << std::endl;
// }

void Stop0LRun2::cleanUp() {
	if (m_electrons) delete m_electrons;
	if (m_muons) delete m_muons;
	if (m_jets) delete m_jets;
}

int Stop0LRun2::isInChan() {
	if (m_debugTree)
		return -2;

	std::vector< bool> chanCuts ={
		((m_nElectrons+m_nMuons)==0), // 0 leptons
		((m_nSigMuons==1 && m_nElectrons==0 && m_nMuons==1)||(m_nSigElectrons==1 && m_nElectrons==1 && m_nMuons==0)), // 1 lepton
		((m_nSigMuons==2 && m_nElectrons==0 && m_nMuons==2)||(m_nSigElectrons==2 && m_nElectrons==2 && m_nMuons==0)) // 2 leptons
	};

	// Cuts common to all of the SRs and CRs
	bool basicCuts =
		(m_hasPrimVtx) &&  // Require one primary vertex
		(m_nBadJets==0) &&  // We don't want bad jets in our events
		(m_nCosmicMuons==0) && // Cosmic muon veto
		(m_nBadMuons==0) // Bad muon veto
		;

	for (unsigned int cI=0; cI<chanCuts.size(); cI++) {
		if (chanCuts[cI] && basicCuts)
			return cI;
	}
	// -1 is not in any channel
	return -1;
}


//channel can be change if debug mode, thus reference
bool Stop0LRun2::performSkim(int& channel, unsigned int& systI, const XAODObjects& xAODObjects) {
	// This is debug mode tree mode
	if (channel == -2)
		{
			channel = 0;
			return true;
		}
	// Check the select region. Skip if m_regionSelec >= 0 (< 0 means do all) and the selected region doesn't match
	else if (m_regionSelec >= 0 && m_regionSelec != channel)
		{
			return false;
		}

	switch (channel)
		{
		case 0:
			{
				// 0-lepton cuts
				std::vector<bool> zeroLepCuts ={
					m_passMetTrigger, // Pass MET trigger (SUSYTools)
					true, // Pass trigger match, not yet implemented obviously
					(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.), // Leading jet pT requirement
					(m_met > 200.),
					m_nJets >= 4,
					true, //m_metTrk>30,
					m_dPhiMetTrackMet < TMath::Pi()/3.,
					m_jetDPhiMetMin3 > TMath::Pi()/5.,
					m_nBJets>=1,
					m_mtTauCand < 0,
					m_mtBMin>175.,
					m_minMt>50,
					m_met>350,
				};

				bool skimCut = true;
				if(m_zeroLepSkimLevel > zeroLepCuts.size())
					m_zeroLepSkimLevel = zeroLepCuts.size();
				for (unsigned int cutI=0; cutI<m_zeroLepSkimLevel; cutI++) {
					skimCut = skimCut && zeroLepCuts[cutI];
				}

				return skimCut;
			}
		case 1:
			{
				std::vector<bool> oneLepCuts ={
					m_passLepTrigger,
					m_passTrigMatch, // Pass trigger match, not yet implemented obviously
					(m_met > 150.),
					m_nJets >= 4,
					true, //m_metTrk>30,
					(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.), // Leading jet pT requirement
					m_dPhiMetTrackMet < TMath::Pi()/3.,
					m_jetDPhiMetMin3 > TMath::Pi()/5.,
					m_mtMetLep>40. && m_mtMetLep<120.,
					m_nBJets>=1,
				};
				bool skimCut = true;
				if(m_oneLepSkimLevel > oneLepCuts.size())
					m_oneLepSkimLevel = oneLepCuts.size();
				for (unsigned int cutI=0; cutI<m_oneLepSkimLevel; cutI++) {
					skimCut = skimCut && oneLepCuts[cutI];
				}

				return skimCut;
			}
		case 2:
			{
				std::vector<bool> twoLepCuts ={
					m_passLepTrigger,
					m_passTrigMatch, // Pass trigger match, not yet implemented obviously
					(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.),
					m_met < 50.,
					m_passOSLep, //twoLeptonProcessors[2][systI]->osLep, 
					m_llMass > 86. && m_llMass<96,
					m_metLep > 70, //twoLeptonProcessors[2][systI]->metLep > 70,
					m_nJets >=4,
					m_nBJets>=2,
				};
				bool skimCut = true;
				if(m_twoLepSkimLevel > twoLepCuts.size())
					m_twoLepSkimLevel = twoLepCuts.size();
				for (unsigned int cutI=0; cutI<m_twoLepSkimLevel; cutI++) {
					skimCut = skimCut && twoLepCuts[cutI];
				}

				return skimCut;
				/*return twoLeptonProcessors[channel][systI]->execute(
			                                                    xAODObjects.electrons, xAODObjects.muons, xAODObjects.jets,
			                                                    m_nSigElectrons, m_nSigMuons);*/

			}
					case -2:
			//Debug mode. Always true. But this should already be check, so I'll print a warning
			std::cerr << "Unexpect behavior in " << __func__ << ". Debug channel reach in switch scope while it should be checked before that." << std::endl;
			channel = 0;
			return true;
		default:
			return false;;
		}
}

void Stop0LRun2::fillCutflow(unsigned int& rI, unsigned int& systI) {
	int chanNum = 0;
	std::vector<bool> commonCuts_0lepton ={
		m_passMetTrigger,
		m_hasPrimVtx,
		(m_nBadJets==0),
		(m_nCosmicMuons==0),
		(m_nBadMuons==0),
		((m_nElectrons+m_nMuons)==0),
		true, // Trigger match that isn't yet implemented
		(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.),
		m_met > 250.,
		m_nJets >= 6,
		m_jetDPhiMetMin3 > TMath::Pi()/5.,
		m_metTrk>30,
		m_dPhiMetTrackMet<TMath::Pi()/3.,
		m_nBJets>=2,
		m_mtTauCand < 0,
		m_mtBMin>175.
	};

	std::vector<bool> sraCuts ={
		50.<m_topDRMinB0M && m_topDRMinB0M<250.,
		50.<m_topDRMinB1M && m_topDRMinB1M<400.,
		m_minMt>50,
		m_met>350
	};

	std::vector<bool> srbCuts ={
		m_mtBMin>175., // this is a really dumb sanity check
		m_antiKt12DijetMassAsym < 0.5,
		m_nAntiKt12>0 && 80.<m_antiKt12M.at(0),
		m_nAntiKt12>1 && 60.<m_antiKt12M.at(1) && m_antiKt12M.at(1)<200.,
		m_nAntiKt8>0 && 50.<m_antiKt8M.at(0),
		m_mtNonBMin>175,
		m_htSig>17
	};
	
	std::vector<bool> commonCuts_1lepton ={
		m_passLepTrigger,
		m_hasPrimVtx,
		(m_nBadJets==0),
		(m_nCosmicMuons==0),
		(m_nBadMuons==0),
		((m_nSigMuons==1 && m_nElectrons==0 && m_nMuons==1)||(m_nSigElectrons==1 && m_nElectrons==1 && m_nMuons==0)),
		m_passTrigMatch, // Trigger match that isn't yet implemented
		(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.),
		m_met > 150.,
		m_nJets >= 4,
		m_jetDPhiMetMin3 > TMath::Pi()/5.,
		m_metTrk>30,
		m_dPhiMetTrackMet<TMath::Pi()/3.,
		m_nBJets>=2,
		m_mtBMin>40. && m_mtBMin<120,
	};
	
	std::vector<bool> commonCuts_2lepton ={
		m_passLepTrigger,
		//m_passMetTrigger,
		m_hasPrimVtx,
		(m_nBadJets==0),
		(m_nCosmicMuons==0),
		(m_nBadMuons==0),
		((m_nSigMuons==2 && m_nElectrons==0 && m_nMuons==2)||(m_nSigElectrons==2 && m_nElectrons==2 && m_nMuons==0)),
		m_passTrigMatch, // Trigger match that isn't yet implemented
		(m_nJets>=2 && m_jetPt.at(0)>80. && m_jetPt.at(1)>80.),
		m_met < 50.,
		m_passOSLep, //twoLeptonProcessors[2][systI]->osLep, 
		m_llMass > 86. && m_llMass<96,
		m_metLep > 70, //twoLeptonProcessors[2][systI]->metLep > 70,
		m_nJets >=4,
		m_nBJets>=2,
	};

	std::vector<std::vector<bool> > commonCuts ={commonCuts_0lepton, commonCuts_1lepton, commonCuts_2lepton};
	// Number of cuts that happen before this one
	unsigned int offset = 3;

	std::vector<unsigned int> lepReqIndex = {5,5,5};
	std::vector<unsigned int> bTagIndex = {13, 13, 13};
	float lepWeight = 1;
	float bTagWeight = 1;
	
	std::vector<bool> passCut ={true, true, true};
	//for (unsigned int rI = 0; rI<commonCuts.size(); rI++) {
	for (unsigned int cI = 0; cI<commonCuts[rI].size(); cI++) {
		passCut[rI] = passCut[rI] && commonCuts[rI][cI];
		if (passCut[rI]) {
			if (m_printEvts)
				std::cout << "evtNum, " << m_cutflow[rI][systI]->GetXaxis()->GetBinLabel(offset+cI+1) << ", " << rI
				          << ": " << m_eventNumber << std::endl;
			m_cutflow[rI][systI]->Fill(offset+cI+0.5);
			if(cI >= lepReqIndex.at(rI))
				lepWeight = m_elecSF*m_muonSF;
			if(cI >= bTagIndex.at(rI))
				bTagWeight = m_bTagSF;
			m_cutflowWeighted[rI][systI]->Fill(offset+cI+0.5, m_mcEventWeight*m_pileupWeight*lepWeight*bTagWeight);
		}
		//}
	}

	if (!m_isData && rI==0) {
		bool passSRACut = passCut[0];
		for (unsigned int cI = 0; cI<sraCuts.size(); cI++) {
			passSRACut = passSRACut && sraCuts[cI];
			if (passSRACut) {
				if (m_printEvts)
					std::cout << "evtNum, " << m_cutflow[0][systI]->GetXaxis()->GetBinLabel(offset+cI+commonCuts_0lepton.size()+1)
					          << ", " << chanNum << ": " << m_eventNumber << std::endl;
				m_cutflow[0][systI]->Fill(offset+cI+commonCuts_0lepton.size()+0.5);
				m_cutflowWeighted[0][systI]->Fill(offset+cI+commonCuts_0lepton.size()+0.5, m_mcEventWeight*m_pileupWeight*m_elecSF*m_muonSF*m_bTagSF);
			}
		}

		bool passSRBCut = passCut[0];
		for (unsigned int cI = 0; cI<srbCuts.size(); cI++) {
			passSRBCut = passSRBCut && srbCuts[cI];
			if (passSRBCut) {
				if (m_printEvts)
					std::cout << "evtNum, " << m_cutflow[0][systI]->GetXaxis()->GetBinLabel(offset+cI+commonCuts_0lepton.size()+sraCuts.size()+m_srOffset)
					          << ", " << chanNum << ": " << m_eventNumber << std::endl;
				m_cutflow[0][systI]->Fill(offset+cI+commonCuts_0lepton.size()+sraCuts.size()+m_srOffset-0.5);
				m_cutflowWeighted[0][systI]->Fill(offset+cI+commonCuts_0lepton.size()+sraCuts.size()+m_srOffset-0.5, m_mcEventWeight*m_pileupWeight*m_elecSF*m_muonSF*m_bTagSF);
			}
		}
	}
	// For Cutflow testing
	if (m_eventNumber == m_evtToPrint) {
		std::cout << "nJets: " << m_nJets << std::endl;
		std::cout << "Jet info:" << std::endl;
		for (int jetI=0; jetI<m_nJets; jetI++) {
			Stop0LUtils::Jet*jet=(Stop0LUtils::Jet*)(*m_jets)[jetI];
			std::cout << "jet " << jetI << ": pT=" << jet->Pt() << ", eta="
			          << jet->Eta() << ", phi=" << jet->Phi() << ", m=" << jet->M()
			          << ", ntrk=" << jet->GetNTrk() << ", dphiMet="
			          << fabs(TVector2::Phi_mpi_pi(jet->Phi()-m_metPhi))
			          << ", mv12c2=" << jet->TagMV1()
			          << std::endl;
		}
		std::cout << std::endl;
		std::cout << "Muon info:" << std::endl;
		for (int muI=0; muI<m_nMuons; muI++) {
			Stop0LUtils::Muon *muon=(Stop0LUtils::Muon*)(*m_muons)[muI];
			std::cout << "muon " << muI << ": pT=" << muon->Pt() << ", eta="
			          << muon->Eta() << ", phi=" << muon->Phi() << ", m=" << muon->M()
			          << std::endl;
		}std::cout << std::endl;

		std::cout << "Electron info:" << std::endl;
		for (int elI=0; elI<m_nElectrons; elI++) {
			Stop0LUtils::Electron *electron=(Stop0LUtils::Electron*)(*m_electrons)[elI];
			std::cout << "electron " << elI << ": pT=" << electron->Pt() << ", eta="
			          << electron->Eta() << ", phi=" << electron->Phi() << ", m=" << electron->M()
			          << std::endl;
		}
		std::cout << std::endl;
		std::cout << "Met info" << std::endl;
		std::cout << "MetTotal=" << m_met << ", MetPhi=" << m_metPhi << std::endl;
		std::cout << std::endl;
		std::cout << "Top0 " << ": pT=" << m_topDRMinB0Pt << ", eta="
		          << m_topDRMinB0Eta << ", phi=" << m_topDRMinB0Phi << ", m=" << m_topDRMinB0M
		          << std::endl;
		std::cout << "Top0 constituents:" << std::endl;
		for (unsigned int dI = 0; dI < m_topDRMinB0Const.size(); dI++) {
			Stop0LUtils::Jet *dJet =(Stop0LUtils::Jet*)(*m_jets)[m_topDRMinB0Const[dI]];
			std::cout << "Jet " << dI << ", pt=" << dJet->Pt()
			          << ", eta=" << dJet->Eta() << ", phi=" << dJet->Phi()
			          << ", mv1=" << dJet->TagMV1() << std::endl << std::endl;
		}
		std::cout << "Top1 " << ": pT=" << m_topDRMinB1Pt << ", eta="
		          << m_topDRMinB1Eta << ", phi=" << m_topDRMinB1Phi << ", m=" << m_topDRMinB1M
		          << std::endl;
		std::cout << "Top1 constituents:" << std::endl;
		for (unsigned int dI = 0; dI < m_topDRMinB1Const.size(); dI++) {
			Stop0LUtils::Jet *dJet =(Stop0LUtils::Jet*)(*m_jets)[m_topDRMinB1Const[dI]];
			std::cout << "Jet " << dI << ", pt=" << dJet->Pt()
			          << ", eta=" << dJet->Eta() << ", phi=" << dJet->Phi()
			          << ", mv1=" << dJet->TagMV1() << std::endl << std::endl;
		}
		if (m_nAntiKt12>0) {
			std::cout << "RCTop0 " << ": pT=" << m_antiKt12Pt.at(0) << ", eta="
			          << m_antiKt12Eta.at(0) << ", phi=" << m_antiKt12Phi.at(0) << ", m=" << m_antiKt12M.at(0)
			          << ", nConst=" << m_antiKt12NConst.at(0) << std::endl;
			std::cout << "RCTop0 constituents:" << std::endl;
			for (unsigned int dI = 0; dI < m_antiKt12Const.at(0).size(); dI++) {
				Stop0LUtils::Jet *dJet =(Stop0LUtils::Jet*)(*m_jets)[m_antiKt12Const.at(0)[dI]];
				std::cout << "Jet " << dI << ", pt=" << dJet->Pt()
				          << ", eta=" << dJet->Eta() << ", phi=" << dJet->Phi()
				          << ", mv1=" << dJet->TagMV1() << std::endl << std::endl;
			}
		}
		if (m_nAntiKt12>1) {
			std::cout << "RCTop1 " << ": pT=" << m_antiKt12Pt.at(1) << ", eta="
			          << m_antiKt12Eta.at(1) << ", phi=" << m_antiKt12Phi.at(1) << ", m=" << m_antiKt12M.at(1)
			          << ", nConst=" << m_antiKt12NConst.at(1) << std::endl;

			std::cout << "RCTop1 constituents:" << std::endl;
			for (unsigned int dI = 0; dI < m_antiKt12Const.at(1).size(); dI++) {
				Stop0LUtils::Jet *dJet =(Stop0LUtils::Jet*)(*m_jets)[m_antiKt12Const.at(1)[dI]];
				std::cout << "Jet " << dI << ", pt=" << dJet->Pt()
				          << ", eta=" << dJet->Eta() << ", phi=" << dJet->Phi()
				          << ", mv1=" << dJet->TagMV1() << std::endl << std::endl;
			}
		}
	}
}



std::vector<std::pair<int, float> > Stop0LRun2::getOrderedJetIndexDR(TLorentzVector& refVector, std::vector<int> skipIndices) {
	std::vector<std::pair<int, float> > jetIndexDRPairVector;

	for (int jI = 0; jI < m_nJets; jI++) {
		if (find(skipIndices.begin(), skipIndices.end(), jI) != skipIndices.end()) {
			//Found the index in skip list, skipping.
			continue;
		}
		TLorentzVector jet;
		jet.SetPtEtaPhiM(m_jetPt.at(jI),
		                 m_jetEta.at(jI),
		                 m_jetPhi.at(jI),
		                 m_jetM.at(jI)
		                 );

		jetIndexDRPairVector.push_back(std::make_pair(jI, refVector.DeltaR(jet)));
	}

	std::sort(jetIndexDRPairVector.begin(), jetIndexDRPairVector.end(), Stop0LRun2::sortBySecond);

	return jetIndexDRPairVector;
}

//Code for printing out neutrio truth used to investigate which to keep. Kept here for reference
//void Stop0LRun2::PrintNeutrinoTruth()
//{
//	MCTruthPartClassifier::ParticleDef particleDef;
//	std::cout << "Event " << m_eventNumber <<": Found " << truthNeutrinos->Pt.size() << " neutrinos from Top/W with sum MET = " << truthNeutrinoMet << std::endl;
//	std::cout << "Found " << neutrinos.size() << " raw neutrinos. Print pT, origin, and ancestors PdgID:" << std::endl;
//	for (auto& nu : neutrinos)
//	{
//		auto originaux = nu->auxdata< unsigned int >("classifierParticleOrigin");
//		auto result = mcTruthClassifier->particleTruthClassifier(nu);
//		//Second in the pair is particle origin
//		auto& origin = result.second;
//		std::cout << "pT: " << nu->pt()/1000 << " GeV\tType: " << particleDef.sParticleType[result.first] << " Origin: " << particleDef.sParticleOrigin[origin] << " Original aux: " << particleDef.sParticleOrigin[originaux] << ". Ancestors:" << std::endl;
//		std::vector<const xAOD::TruthParticle*> ancestors;
//		xAODtruthUtils::getAncestors(nu, ancestors);
//		for (auto& ancestor : ancestors)
//		{
//			auto ancestorResult = mcTruthClassifier->particleTruthClassifier(ancestor);
//			std::cout << "\tPdgID: " << ancestor->pdgId() << " Type: " << particleDef.sParticleType[ancestorResult.first] << " Origin: " << particleDef.sParticleOrigin[ancestorResult.second] << std::endl;
//			if (ancestor->isW() || ancestor->isTau())
//			{
//				std::vector<const xAOD::TruthParticle*> grandAncestors;
//				xAODtruthUtils::getAncestors(ancestor, grandAncestors);
//				for (auto& grandAncestor : grandAncestors)
//				{
//					auto grandAncestorResult = mcTruthClassifier->particleTruthClassifier(grandAncestor);
//					std::cout << "\t\tPdgID: " << grandAncestor->pdgId() << " Type: " << particleDef.sParticleType[grandAncestorResult.first] << " Origin: " << particleDef.sParticleOrigin[grandAncestorResult.second] << std::endl;
//				}
//			}
//		}
//	}
//	std::cin.get();
//}


void Stop0LRun2::getProjectTag(std::string sampleName, std::string &projectTag, bool &isData) {
	boost::regex projectTagRegEx("(mc|data)\\d+_\\d+TeV");
	boost::smatch result;

	bool hasProjectTag = boost::regex_search(sampleName, result, projectTagRegEx);
	if (hasProjectTag) {
		projectTag = result[0];
		isData = std::string(result[1]).find("data") != std::string::npos;
	}
}

void Stop0LRun2::getTags(std::string sampleName,
                         std::vector<int> &aTags, std::vector<int> &sTags, std::vector<int> &rTags, std::vector<int> &eTags, std::vector<int> &pTags) {
	boost::regex aTagRegEx("(\\.|_)a(\\d{3,5})");
	boost::regex sTagRegEx("(\\.|_)s(\\d{3,5})");
	boost::regex rTagRegEx("(\\.|_)r(\\d{3,5})");
	boost::regex pTagRegEx("(\\.|_)p(\\d{3,5})");
	boost::regex eTagRegEx("(\\.|_)e(\\d{3,5})");

	boost::sregex_token_iterator iter(sampleName.begin(), sampleName.end(), aTagRegEx, 2);
	boost::sregex_token_iterator end;

	for (; iter != end; ++iter) {
		aTags.push_back(atoi(std::string(*iter).c_str()));
	}

	iter = boost::sregex_token_iterator(sampleName.begin(), sampleName.end(), rTagRegEx, 2);
	for (; iter != end; ++iter) {
		rTags.push_back(atoi(std::string(*iter).c_str()));
	}

	iter = boost::sregex_token_iterator(sampleName.begin(), sampleName.end(), sTagRegEx, 2);
	for (; iter != end; ++iter) {
		sTags.push_back(atoi(std::string(*iter).c_str()));
	}

	iter = boost::sregex_token_iterator(sampleName.begin(), sampleName.end(), eTagRegEx, 2);
	for (; iter != end; ++iter) {
		eTags.push_back(atoi(std::string(*iter).c_str()));
	}

	iter = boost::sregex_token_iterator(sampleName.begin(), sampleName.end(), pTagRegEx, 2);
	for (; iter != end; ++iter) {
		pTags.push_back(atoi(std::string(*iter).c_str()));
	}

}


void Stop0LRun2::GetExtraJetCollections()
{
	for (auto& jetCollectonName : extraJetCollections)
		{
			if (!m_event->contains<xAOD::JetContainer>(jetCollectonName))
				{
					continue;
				}
			const xAOD::JetContainer* jets = 0;

			if (m_event->retrieve(jets, jetCollectonName).isSuccess())
				{
					auto jetsCopy = xAOD::shallowCopyContainer(*jets);
					xAOD::setOriginalObjectLink(*jets, *jetsCopy.first);

					CHECK(m_store->record(jetsCopy.first, jetCollectonName+"ShallowCopiedJets"));
					CHECK(m_store->record(jetsCopy.second, jetCollectonName+"ShallowCopiedJetsAux."));
					auto calibTool = rScanCalibTools.at(jetCollectonName).get();

					extraJetNtup4VecsMap->at(jetCollectonName)->Reset();
					for (auto *ojet : *jets)
						{
							xAOD::Jet * _jet = 0;

							calibTool->calibratedCopy(*ojet, _jet).ignore();
							std::unique_ptr<xAOD::Jet> jet(_jet);
							auto jvt = jvtTool->updateJvt(*jet);
							if (jvt < 0.64)
								{
									continue;
								}
							if (!jetCleaningTool->accept(*jet))
								{
									continue;
								}

							extraJetNtup4VecsMap->at(jetCollectonName)->push_back_mev_to_gev(jet->p4());
						}
				}
		}
}

#if DOSUBSTRUCTURE
void Stop0LRun2::DoTrackSubstructures(XAODObjects& xAODObjects)
{
	// Get input vertex collection
	//const xAOD::VertexContainer* vertexContainer = nullptr;
	//if ( m_store->retrieve(vertexContainer,m_vertexContainer).isFailure()
	//     || vertexContainer == nullptr ) {
	//  ATH_MSG_ERROR("Could not retrieve the VertexContainer from evtStore: "
	//                << m_vertexContainer);
	//  return 1;
	//}
	//
	//// Get the track-vertex association
	//const jet::TrackVertexAssociation* tva = nullptr;
	//if ( evtStore()->retrieve(tva,m_tva).isFailure() || tva==nullptr ) {
	//  ATH_MSG_ERROR("Could not retrieve the TrackVertexAssociation from evtStore: "
	//                << m_tva);
	//  return 2;
	//}
	auto jets = xAODObjects.jets;
	for (auto&& jet : *jets)
		{
			// Retrieve the associated tracks.
			std::vector<const xAOD::TrackParticle*> tracks;
			bool havetracks = jet->getAssociatedObjects("GhostTrack", tracks);
			//if (!havetracks) ATH_MSG_WARNING("Associated tracks not found");
			ATH_MSG_DEBUG("Successfully retrieved track particles");

			//For PFlow jets we will also calculate the same moments, using charged PFO                                                                                                                                                                                                        
			xAOD::Type::ObjectType ctype = jet->rawConstituent(0)->type();
			std::vector<const xAOD::TrackParticle*> pflowTracks;
			bool isPFlowJet = false;
			if (ctype  == xAOD::Type::ParticleFlow) {
				isPFlowJet = true;
				size_t numConstit = jet->numConstituents();
				for (size_t i=0; i<numConstit; i++) {
					const xAOD::PFO* constit = dynamic_cast<const xAOD::PFO*>(jet->rawConstituent(i));
					if (0.0 != constit->charge()) {
						const xAOD::TrackParticle *thisTrack = constit->track(0);//by construction xAOD::PFO can only have one track, in eflowRec usage                                                                                                                                             
						pflowTracks.push_back(thisTrack);
					}//we have a charged PFO                                                                                                                                                                                                                                                      
				}//loop on jet constituents         
			}//yes this jet is made form xAOD::PFO, so we do calculate the pflow moments           

			// For each track cut, get the associated moments
			for (auto&& track : tracks)
				{
					if (!trackSelector->accept(track)) continue;

					using namespace fastjet;
					std::vector<PseudoJet> trackAsJets;
					trackAsJets.emplace_back(track->p4().Px(), track->p4().Py(), track->p4().Pz(), track->p4().T());

					auto trackAsAJet = join(trackAsJets);
				}
		}
}

void Stop0LRun2::DoSubstructures()
{
}
#endif // DOSUBSTRUCTURE


#if DOPFLOW
inline void Stop0LRun2::InitPFlowSUSYTools(ST::SettingDataSource dataSource)
{
	pflowSUSYTools = std::make_unique<ST::SUSYObjDef_xAOD>("SUSYToolsPF");
	CHECK(pflowSUSYTools->setProperty("JetInputType", xAOD::JetInput::EMPFlow));
	CHECK(pflowSUSYTools->setProperty("BtagWP", ""));
	CHECK(pflowSUSYTools->setProperty("BtagWP_OR", ""));
	CHECK(pflowSUSYTools->setProperty("DataSource", dataSource));
	CHECK(pflowSUSYTools->setProperty("MuId", xAOD::Muon::Medium));

	if (m_debug) pflowSUSYTools->msg().setLevel(MSG::VERBOSE);
	CHECK(pflowSUSYTools->setProperty("DebugMode", m_debug));

	if (pflowSUSYTools->SUSYToolsInit().isFailure()) {
		Error(APP_NAME, "Failed to initialise tools in SUSYToolsInit()...");
		Error(APP_NAME, "Exiting...");
		throw std::runtime_error("");
	}

	if (pflowSUSYTools->initialize() != EL::StatusCode::SUCCESS) {
		Error(APP_NAME, "Cannot intialize SUSYToolsPF...");
		Error(APP_NAME, "Exiting... ");
		throw std::runtime_error("");
	}
	else {
		Info(APP_NAME, "SUSYToolsPF initialized... ");
	}
}
#endif // DOPFLOW


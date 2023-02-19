/*************************************************************************

 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

//==================================================================
// Simple class for dn/deta analyses.
// by Beomkyu KIM
//==================================================================

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TGeoGlobalMagField.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliCentrality.h"
#include "AliVMultiplicity.h"
#include "AliPWG0Helper.h"
#include "AliFMDEventInspector.h"
#include "AliMultSelection.h"
#include "AliESDtrack.h"
///#include "AliGenPythiaPlus.h"
#include "AliGenHijingEventHeader.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSMultReconstructor.h"
#include "AliITSsegmentationSPD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"
#include "AliITSMultRecBg.h"
#include "AliESDInputHandlerRP.h"
#include "AliMagF.h"
#include "AliGRPObject.h"
#include "AliMultiplicity.h"
#include "AliVVZERO.h"
#include "AliAnalysisPseudoRapidityDensityTemp.h"
using namespace std;
const Double_t pi = TMath::Pi();
const Int_t kNstable = 20;
const Int_t pdgStable[20] = {
	22,	  // Photon
	11,	  // Electron
	12,	  // Electron Neutrino
	13,	  // Muon
	14,	  // Muon Neutrino
	15,	  // Tau
	16,	  // Tau Neutrino
	211,  // Pion
	321,  // Kaon
	311,  // K0
	130,  // K0s
	310,  // K0l
	2212, // Proton
	2112, // Neutron
	3122, // Lambda_0
	3112, // Sigma Minus
	3222, // Sigma Plus
	3312, // Xsi Minus
	3322, // Xsi0
	3334  // Omega
};
//_________________________________________________________________
AliAnalysisPseudoRapidityDensityTemp::AliAnalysisPseudoRapidityDensityTemp()
	: AliAnalysisTaskSE("AliAnalysisPseudoRapidityDensityTemp"), fOption()
{
	DefineInput(0, TChain::Class());
	//DefineOutput (0, TTree::Class());
	DefineOutput(1, TList::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensityTemp::AliAnalysisPseudoRapidityDensityTemp(
	const char *name, const char *option)
	: AliAnalysisTaskSE(name), fOption(option)
{
	DefineInput(0, TChain::Class());
	//DefineOutput(1, TTree::Class());
	DefineOutput(1, TList::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensityTemp::AliAnalysisPseudoRapidityDensityTemp(
	const AliAnalysisPseudoRapidityDensityTemp &ap)
	: fOption(ap.fOption)
{
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensityTemp &AliAnalysisPseudoRapidityDensityTemp::operator=(
	const AliAnalysisPseudoRapidityDensityTemp &ap)
{
	// assignment operator

	this->~AliAnalysisPseudoRapidityDensityTemp();
	new (this) AliAnalysisPseudoRapidityDensityTemp(ap);
	return *this;
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensityTemp::~AliAnalysisPseudoRapidityDensityTemp()
{
	//delete fHistos;
	delete fTrigger;
	delete fPIDResponse;
	delete fRandom;
}

//___________________________________________________________________
void AliAnalysisPseudoRapidityDensityTemp::UserCreateOutputObjects()

{
	fRandom = new TRandom3;
	fRandom->SetSeed();
	// Histograms container

	// Offline triggers -----------------------------------------------------
	fTrigger = new AliTriggerAnalysis;	 // offline trigger
	fTrigger->SetFMDThreshold(0.3, 0.5); // FMD threshold
	//-----------------------------------------------------------------------

	// TrackCuts for strangeness measure------------------------------------i

	AliESDtrackCuts *tempcutset = new AliESDtrackCuts();
	fTrackCuts.resize(kTrackCutend);
	fTrackCuts.at(0) = *tempcutset;
	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false); 
	tempcutset->SetMaxDCAToVertexXY(2.4);
	tempcutset->SetMaxDCAToVertexZ(3.2);
	tempcutset->SetDCAToVertex2D(kTRUE);
	tempcutset->SetMaxChi2TPCConstrainedGlobal(36);
	tempcutset->SetMaxFractionSharedTPCClusters(0.4);
	fTrackCuts.at(kHybrid) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false); 
	tempcutset->SetMaxDCAToVertexXY(2.4);
	tempcutset->SetMaxDCAToVertexZ(3.2);
	tempcutset->SetDCAToVertex2D(kTRUE);
	tempcutset->SetMaxChi2TPCConstrainedGlobal(36);
	tempcutset->SetMaxFractionSharedTPCClusters(0.4);
	tempcutset->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
	tempcutset->SetRequireITSRefit(kTRUE);
	fTrackCutGC = *tempcutset;

	AliESDtrackCuts *esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	fTrackCuts.at(kITSTPC2011) = *esdTrackCutsH2;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxDCAToVertexZ(1);
	fTrackCuts.at(kITSTPC2011dcazdw) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxDCAToVertexZ(5);
	fTrackCuts.at(kITSTPC2011dcazup) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxDCAToVertexXYPtDep("3*(0.0015+0.0050/pt^1.1)");
	fTrackCuts.at(kITSTPC2011dcardw) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxDCAToVertexXYPtDep("10*(0.0015+0.0050/pt^1.1)");
	fTrackCuts.at(kITSTPC2011dcarup) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMinNCrossedRowsTPC(50);
	fTrackCuts.at(kITSTPC2011nclutpcdw) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMinNCrossedRowsTPC(100);
	fTrackCuts.at(kITSTPC2011nclutpcup) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxChi2PerClusterTPC(3);
	fTrackCuts.at(kITSTPC2011chitpcdw) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxChi2PerClusterTPC(5);
	fTrackCuts.at(kITSTPC2011chitpcup) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxChi2TPCConstrainedGlobal(25);
	fTrackCuts.at(kITSTPC2011globalconsdw) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true,1);
	tempcutset->SetMaxChi2TPCConstrainedGlobal(49);
	fTrackCuts.at(kITSTPC2011globalconsup) = *tempcutset;
	/*
	AliESDtrackCuts dcaz(fTrackCuts.at(kITSTPC2010));
	dcaz.SetMaxDCAToVertexZ(1);
	fTrackCuts.at(kITSTPC2010dcazdw) = dcaz;
	dcaz.SetMaxDCAToVertexZ(5);
	fTrackCuts.at(kITSTPC2010dcazup) = dcaz;

	AliESDtrackCuts dcar(fTrackCuts.at(kITSTPC2010));
	dcar.SetMaxDCAToVertexXYPtDep("10*(0.0026+0.0050/pt^1.01)");
	fTrackCuts.at(kITSTPC2010dcarup) = dcar;
	dcar.SetMaxDCAToVertexXYPtDep("3*(0.0026+0.0050/pt^1.01)");
	fTrackCuts.at(kITSTPC2010dcardw) = dcar;

	AliESDtrackCuts nclutpc(fTrackCuts.at(kITSTPC2010));
	nclutpc.SetMinNCrossedRowsTPC(100);
	fTrackCuts.at(kITSTPC2010nclutpcdw) = nclutpc;
	nclutpc.SetMinNCrossedRowsTPC(130);
	fTrackCuts.at(kITSTPC2010nclutpcup) = nclutpc;

	AliESDtrackCuts chitpc(fTrackCuts.at(kITSTPC2010));
	chitpc.SetMaxChi2PerClusterTPC(3);
	fTrackCuts.at(kITSTPC2010chitpcdw) = chitpc;
	chitpc.SetMaxChi2PerClusterTPC(5);
	fTrackCuts.at(kITSTPC2010chitpcup) = chitpc;

	AliESDtrackCuts globalcons(fTrackCuts.at(kITSTPC2010));
	globalcons.SetMaxChi2TPCConstrainedGlobal(25);
	fTrackCuts.at(kITSTPC2010globalconsdw) = globalcons;
	globalcons.SetMaxChi2TPCConstrainedGlobal(49);
	fTrackCuts.at(kITSTPC2010globalconsup) = globalcons;
	*/

	fHistos = new THistManager("dndeta");

	//auto binType = AxisStr("Type",{"PN","PP","NN","Mixing"});
	//Double1D varcentbin = {0, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 100};
	Double1D varcentbin = {0, 0.01, 0.05,0.1,1,5, 10, 20, 30, 40, 50, 60, 100};
	if (fOption.Contains("MB"))
		varcentbin = {0, 100};
	Double1D varcentbinHeavy = {0, 2.5, 5, 7.5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

	//for (auto i=1; i<=100; i++) varcentbin.push_back(i);
	auto binCent = AxisVar("Cent", IsAA ? varcentbinHeavy : varcentbin);

	auto binEta = AxisFix("Eta", 60, -3, 3);
	auto binPhi = AxisVar("Phi", {0, pi * 2 / 3, pi * 4 / 3, 2 * pi});
	auto binPt = AxisLog("Pt", 20,  0.15, 20, 0.15);
	auto binCentClass = AxisStr("CentClass", {"V0M"});
	//	auto binEventClass = AxisStr("EventClass", {"DATA", "INEL", "INELg0" , "INEL015", "INEL040",  "INEL050", "INEL100", "INEL200","INEL300","INEL04025", "INEL05024","INEL05025"});
	auto binEventClass = AxisStr("EventClass", {"DATA", "INEL", "INELg0" , "INEL300", "INEL100", "INEL200"}); //remove unneccessary bins
	//	auto binTriggClass = AxisStr("TriggClass", {"MBAND015","MBAND040","MBAND050", "MBAND100", "MBAND200", "MBAND300", "HMMBAND" "HMMBAND300"});
	//	auto binTriggClass = AxisStr("TriggClass", {"HMMBAND","HMMBAND300","MBAND050", "MBAND015", "MBAND100", "MBAND300", "MBAND200" "MBAND040"});
	auto binTriggClass = AxisStr("TriggClass", {"HMMBAND","HMMBAND300", "MBAND015", "MBAND300","MBAND100", "MBAND200"}); //remove unneccessary bins
	auto binParType = AxisStr("ParticleType", {"Tracks", "MotherStrange", "Bkg", "Pion", "Kaon", "Proton", "Opar"});
	auto binV0Type = AxisStr("V0Type", {"K0s", "Lambda", "AntiLambda"});
	auto binTrkCutVar = AxisStr("TrkCutVar", {
		"Hybrid", 
		"ITSTPC2011",
		"ITSTPC2010dcazdw", 
        "ITSTPC2010dcazup",
        "ITSTPC2010dcardw", 
        "ITSTPC2010dcarup",
        "ITSTPC2010nclutpcdw", 
        "ITSTPC2010nclutpcup",
        "ITSTPC2010chitpcdw", 
        "ITSTPC2010chitpcup",
        "ITSTPC2010globalconsdw", 
        "ITSTPC2010globalconsup"});

	Double1D binmul = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
	auto binMul = AxisVar("multiplicity", binmul);

	//auto binMass = AxisFix("Mass",200,0,2);

	std::vector<TString> ent = {"All", "PS", "PSpileup", "Goodz", "Goodzcut", "INEL015", "INEL040", "INEL050", "INEL100", "INEL200", "INEL300", "Goodzcut7", "INELg010", "MBANDg010", "INELg0", "PSINELg0", "MBANDg0"};
	auto h = fHistos->CreateTH1("hEventNumbers", "", ent.size(), 0, ent.size());
	fHistos->CreateTH1("quickcheck", "", 10, 0, 10);
	for (auto i = 0u; i < ent.size(); i++)
		h->GetXaxis()->SetBinLabel(i + 1, ent.at(i).Data());
	fHistos->CreateTH2("hPhiEta", "", 180, 0, 2 * pi, 40, -2, 2);
	fHistos->CreateTH2("hPhiEtaCut", "", 180, 0, 2 * pi, 40, -2, 2);
	//Double1D zbins = {-20,-15,-10,-7,-3,0,3,7,10,15,20};
	//Double1D zbins = {-20,-15,-10,-7,-3,0,3,7,10,15,20};
	Double1D zbins = {-15, -10, -7, -6, -5, -4,-3, -2,-1, 0, 1, 2 , 3, 4, 5 , 6 ,7, 10, 15};
	TAxis binZ = AxisVar("Z", zbins);

	fHistos->CreateTH1("hdndeta", "", 60, -6, 6, "s");
	fHistos->CreateTH1("hdndetacut", "", 60, -6, 6, "s");
	fHistos->CreateTH1("hdndetamc", "", 60, -6, 6, "s");
	fHistos->CreateTH1("hdndetaINEL015", "", 1200, -6, 6, "s");
	fHistos->CreateTH1("hdndetaINEL040", "", 1200, -6, 6, "s");
	fHistos->CreateTH1("hdndetaINEL050", "", 1200, -6, 6, "s");
	fHistos->CreateTH1("hdndetaINEL100", "", 1200, -6, 6, "s");
	fHistos->CreateTH1("hdndetaINEL200", "", 1200, -6, 6, "s");
	fHistos->CreateTH1("hGenMultV0AGeo", "", 100000, 0, 1000, "s");
	fHistos->CreateTH1("hGenMultV0CGeo", "", 100000, 0, 1000, "s");
	fHistos->CreateTH1("hGenMultV0MGeo", "", 100000, 0, 1000, "s");
	fHistos->CreateTH1("hGenMultSPDGeo", "", 100000, 0, 1000, "s");
	fForward2d = fHistos->CreateTH2("hfmddetadphi", "", 200, -4, 6, 20, 0, 2 * pi, "s");
	fHistos->CreateTH1("hcentTracklets", "hcentTracklets", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentselection", "hcentselection", 100, 0, 100, "s");
	fHistos->CreateTH1("ztruth", "ztuth", 60, -30, 30, "s");
	fHistos->CreateTH1("zdata", "zdata", 60, -30, 30, "s");
	fHistos->CreateTH1("hchi2", "SPD chi2", 600, 0, 6, "s");
	fHistos->CreateTH1("hcentV0M", "hcentV0M", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0MHighMult", "hcentV0MHighMult", 1000, 0, 1, "s");
	fHistos->CreateTH1("hcentSPD", "hcentSPD", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0Mzcut7", "hcentV0Mzcut7", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentSPDzcut7", "hcentSPDzcut7", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0Mzcut10", "hcentV0Mzcut10", 100, 0., 100, "s");
	fHistos->CreateTH1("hcentSPDzcut10", "hcentSPDzcut10", 100, 0., 100, "s");
	fHistos->CreateTH1("hcentRef05zcut10", "hcentRef05zcut10", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentRef08zcut10", "hcentRef08zcut10", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0Mzcut10MBAND", "hcentV0Mzcut10MBAND:", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentSPDzcut10MBAND", "hcentSPDzcut10MBAND", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentRef05zcut10MBAND", "hcentRef05zcut10MBAND", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentRef08zcut10MBAND", "hcentRef08zcut10MBAND", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0Mzcut20", "hcentV0Mzcut20", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentSPDzcut20", "hcentSPDzcut20", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentV0Mzcut15", "hcentV0Mzcut15", 100, 0, 100, "s");
	fHistos->CreateTH1("hcentSPDzcut15", "hcentSPDzcut15", 100, 0, 100, "s");
	fHistos->CreateTH1("hINELg0Nch", "hINELg0Nch", 300, 0, 300, "s");
	fHistos->CreateTH1("hINT7g0Nch", "hINT7g0Nch", 300, 0, 300, "s");
	fHistos->CreateTH1("hEtTruth", "hEtTruth", 100, 0, 100, "s");
	fHistos->CreateTH1("hEtMCrec", "hEtMCrec", 100, 0, 100, "s");
	fHistos->CreateTH2("hMultResponseV0M", "hMultResponseV0M", 200, 0, 200, 800, 0, 800, "s");
	fHistos->CreateTH2("hMultResponseSPD", "hMultResponseSPD", 200, 0, 200, 200, 0, 200, "s");
	fHistos->CreateTH1("hPhiHybrid", "hPhiHybrid", 360, 0, TMath::Pi(), "s");

	const int nbins = 50;
	Double_t logbins[nbins + 1];
	Double_t low = 1;
	Double_t high = 2000;
	Double_t logbw = (log(high) - log(low)) / nbins;
	for (int ij = 0; ij <= nbins; ij++)
		logbins[ij] = low * exp(ij * logbw);
	fHistos->CreateTH1("hsdmass", "SD mass", nbins, logbins);
	low = 0.001;
	high = 50;
	logbw = (log(high) - log(low)) / nbins;
	for (int ij = 0; ij <= nbins; ij++)
		logbins[ij] = low * exp(ij * logbw);
	fHistos->CreateTH1("hkinept", "Kine only pt", nbins, logbins);
	fHistos->CreateTH1("hmcrecpt", "MC rec pt", nbins, logbins);

	CreateTHnSparse("hrecdndeta", "rec dndeta", 9, {binEventClass, binTriggClass, binCent, binZ, binParType, binTrkCutVar, binEta, binCentClass, binPhi}, "s");
	CreateTHnSparse("hreczvtx", "rec zvtx", 5, {binEventClass, binTriggClass, binCent, binZ, binCentClass}, "s");
	CreateTHnSparse("heventcount", "event count", 3, {binTriggClass, binCent, binCentClass}, "s");

	CreateTHnSparse("hkinedndeta", "kine only dndeta", 8, {binEventClass, binCent, binZ, binParType, binTrkCutVar, binEta, binCentClass, binPhi}, "s");
	CreateTHnSparse("hkinezvtx", "kine zvtx", 4, {binEventClass, binCent, binZ, binCentClass}, "s");

	if (!fOption.Contains("MC")){
	CreateTHnSparse("hrecdndpt", "rec dndpt", 2, {binEta, binPt}, "s");
	CreateTHnSparse("htruedndpt", "rec dndpt", 2, {binEta, binPt}, "s");
	CreateTHnSparse("hdndptres", "response dndpt", 3, {binEta, binPt, binPt}, "s");
	CreateTHnSparse("hdndptfake", "fake dndpt", 2, {binEta, binPt}, "s");
	CreateTHnSparse("hdndptmiss", "miss dndpt", 2, {binEta, binPt}, "s");
	CreateTHnSparse("hrecmult", "rec mult", 1, {binMul }, "s");
	CreateTHnSparse("htruemult", "rec dndpt", 1, {binMul}, "s");
	CreateTHnSparse("hmultres", "response dndpt", 2, {binMul, binMul}, "s");
	CreateTHnSparse("hmultfake", "fake dndpt", 1, {binMul}, "s");
	CreateTHnSparse("hmultmiss", "miss dndpt", 1, {binMul}, "s");
	CreateTHnSparse("hv0mass", "V0 inv mass", 5,
					{binTriggClass, binCent, binV0Type, AxisFix("v0mass", 1200, 0.3, 1.5), binCentClass}, "s");
	}
	CreateTHnSparse("hv0eta", "V0 daughter eta", 5,
					{binTriggClass, binCent, binV0Type, binEta, binCentClass}, "s");

	if (!fOption.Contains("MC")){
	CreateTHnSparse("hkineGeoMult", "", 2, {binCentClass, AxisFix("GeoMult", 30000, 0, 300)}, "s");
	CreateTHnSparse("hkinedNdEtaGeoCent", "", 3, {binCent, binEta, binCentClass}, "s");
	CreateTHnSparse("hcentGeo", "", 2, {binCentClass, binCent}, "s");
	CreateTHnSparse("hMult", "", 2, {binCentClass, AxisFix("Mult", 20000, 0, 2000)}, "s");
	CreateTHnSparse("hMultHigh", "", 2, {binCentClass, AxisFix("Mult", 20000, 0, 2000)}, "s");
	CreateTHnSparse("hMultcent", "", 3, {binCentClass, binCent, AxisFix("Multcent", 20000, 0, 2000)}, "s");
	CreateTHnSparse("hMultHighcent", "", 3, {binCentClass, binCent, AxisFix("highMultcent", 20000, 0, 2000)}, "s");
	CreateTHnSparse("hdNdEtaCent", "", 3, {binCent, binEta, binCentClass}, "s");
	CreateTHnSparse("hcent", "", 2, {binCentClass, binCent}, "s");
	CreateTHnSparse("hPhiEtaTracks", "", 3, {binTrkCutVar, AxisFix("phi", 180, 0, 2 * pi), AxisFix("eta", 40, -2, 2)}, "s");
	}

	PostData(1, fHistos->GetListOfHistograms());
	//PID Combined
	fPIDCombined = new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors(); //Need more update..
	fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
	fPIDCombined->SetDetectorMask(
		AliPIDResponse::kDetTPC |
		AliPIDResponse::kDetTOF |
		AliPIDResponse::kDetITS |
		AliPIDResponse::kDetTRD); //Do we need??

	Double1D fcent(kCentClassBinEnd, -1);
	fCent = fcent;

	if (fOption.Contains("BgInjection"))
		isbginjection = true;
}

void AliAnalysisPseudoRapidityDensityTemp::InitMultReco()
{
	// create mult reconstructor
	if (fMultReco)
		delete fMultReco;
	fMultReco = new AliITSMultRecBg();
	fMultReco->SetCreateClustersCopy(kTRUE);
	fMultReco->SetScaleDThetaBySin2T(true);
	fMultReco->SetNStdDev(25.);
	fMultReco->SetPhiWindow(0.06);
	fMultReco->SetThetaWindow(0.025);
	fMultReco->SetPhiShift(0.0045);
	fMultReco->SetRemoveClustersFromOverlaps(true);
	fMultReco->SetPhiOverlapCut(0.005);
	fMultReco->SetZetaOverlapCut(0.05);
	fMultReco->SetHistOn(kFALSE);
	fMultReco->SetRecType(AliITSMultRecBg::kData);
}

//___________________________________________________________________
void AliAnalysisPseudoRapidityDensityTemp::UserExec(Option_t *)
{

	// Pointer to a event----------------------------------------------------
	AliVEvent *event = InputEvent();
	if (!event)
	{
		Printf("ERROR: Could not retrieve event");
		return;
	}

	// connect to ESD tree --------------------------------------------------
	Int_t runnumber;
	event->IsA() == AliESDEvent::Class()
		? fEvt = dynamic_cast<AliESDEvent *>(event)
		: fEvt = dynamic_cast<AliAODEvent *>(event);
	if (!fEvt)
		return;

	if (IsFirstEvent)
	{
		runnumber = fEvt->GetRunNumber();
		IsFirstEvent = false;
	}


	// ----------------------------------------------------------------------
	// centrality
	Double1D fcent(kCentClassBinEnd, -1);
	fCent = fcent;
	sel = (AliMultSelection *)fEvt->FindListObject("MultSelection");
	//if ( sel->GetEvSelCode() <= 0 ) return;
	//if (! sel ->	GetThisEventPassesTrackletVsCluster() && !fOption.Contains("MC")) return;
	//sel->SetThisEventIsNotAsymmetricInVZERO(true);

	if (sel)
	{
		fCent[kV0M] = sel->GetMultiplicityPercentile("V0M");
		if (fOption.Contains("LHC16q"))
			fCent[kV0M] = sel->GetMultiplicityPercentile("V0A");
	}
	//if(! sel->IsEventSelected() && !fOption.Contains("MC") && !fOption.Contains("LHC10d")) {
	//	fCent[kV0M] = sel->GetEvSelCode();
	//	fCent[kV0A] = sel->GetEvSelCode();
	//	fCent[kV0C] = sel->GetEvSelCode();
	//	fCent[kSPDMult] = sel->GetEvSelCode();
	//}

	if (!fOption.Contains("MC"))
	{
		fHistos->FillTH1("quickcheck", 1);
		if (sel->GetThisEventIsNotPileup())
			fHistos->FillTH1("quickcheck", 2);
		if (sel->GetThisEventIsNotPileupInMultBins())
			fHistos->FillTH1("quickcheck", 3);
		if (sel->GetThisEventHasNoInconsistentVertices())
			fHistos->FillTH1("quickcheck", 4);
		if (sel->GetThisEventPassesTrackletVsCluster())
			fHistos->FillTH1("quickcheck", 5);
		if (/*!sel->GetThisEventIsNotPileup() ||*/ !sel->GetThisEventIsNotPileupInMultBins() || !sel->GetThisEventHasNoInconsistentVertices() || !sel->GetThisEventPassesTrackletVsCluster())
		{
			for (auto i = 0; i < kCentClassBinEnd; i++)
				fCent.at(i) = -1;
		}
		else
		{
			fHistos->FillTH1("hcentV0M", fCent[kV0M]);
		}
	}

	if (sel->IsEventSelected())
	{
	}
	// Pointer to a MC event-------------------------------------------------
	AliMCEvent *mcEvent = MCEvent();
	// ----------------------------------------------------------------------

	Bool_t IsMC = kFALSE;
	TArrayF vtxMC3d(3);

	Double_t eta = -10, phi = -10, pt = -1.;
	AliGenPythiaEventHeader *pythiaGenHeader = NULL;
	AliGenDPMjetEventHeader *dpmHeader = NULL;
	AliGenHijingEventHeader *hijing = NULL;

	Bool_1d bevtc(kECend, false);
	sdweightingfactor = 1.;
	Bool_t issd = 0;
	Bool_t isdd = 0;
	Int_t Nch_INELg0_eta1 = 0;
	Double_t nv0mgeo = 0;
	Double_t nspdgeo = 0;
	if (mcEvent)
	{
		ismc = true;
		pythiaGenHeader =
			dynamic_cast<AliGenPythiaEventHeader *>(mcEvent->GenEventHeader());
		dpmHeader =
			dynamic_cast<AliGenDPMjetEventHeader *>(mcEvent->GenEventHeader());
		hijing =
			dynamic_cast<AliGenHijingEventHeader *>(mcEvent->GenEventHeader());

		if (pythiaGenHeader)
		{ // Pythia6
			//92 SD1, 93 SD2, 94 DD, -2 CD , -1 ND, 91 EL
			Int_t proc = pythiaGenHeader->ProcessType();
			if (proc == 92 || proc == 93)
				issd = true;
			else if (proc == 94)
				isdd = true;
			//else if (proc == -1 || proc == -2) bevtc[kND] = true;

			if (fOption.Contains("PYTHIA8"))
			{
				// 103 SD1, 104 SD2, 105 DD
				if (proc == 103 || proc == 104)
					issd = true;
				else if (proc == 105)
					isdd = true;
			}
		}
		else if (dpmHeader)
		{
			Int_t proc = dpmHeader->ProcessType();
			if (proc == 5 || proc == 6)
				issd = true;
			else if (proc == 7)
				isdd = true;
		}


		bevtc[kINEL] = true;

		if (fOption.Contains("EPOS")){
			bevtc[kINEL] = true; //EPOS doesn't have diff info
		}

		stack = mcEvent->Stack();
		Int_t nPrim = stack->GetNprimary();
		mcEvent->GenEventHeader()->PrimaryVertex(vtxMC3d);
		if (bevtc[kINEL])
			fHistos->FillTH1("ztruth", vtxMC3d[2]);
		for (Int_t i = 0; i < nPrim; i++)
		{
			TParticle *part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE)
				continue;
			eta = part->Eta();
			pt = part->Pt();
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 0.15)
				bevtc[kINEL015] = true;
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 0.5)
				bevtc[kINEL050] = true;
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 0.4)
				bevtc[kINEL040] = true;
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 1)
				bevtc[kINEL100] = true;
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 2)
				bevtc[kINEL200] = true;
			if (bevtc[kINEL] && fabs(eta) < 0.8 && pt > 4) // To Chiara: This is OK because it is the event condition you are saying
				bevtc[kINEL300] = true;
			if (bevtc[kINEL] && fabs(eta) < 2.4 && pt > 0.5)
				bevtc[kINEL05024] = true;
			if (bevtc[kINEL] && fabs(eta) < 2.5 && pt > 0.5)
				bevtc[kINEL05025] = true;
			if (bevtc[kINEL] && fabs(eta) < 2.5 && pt > 0.4)
				bevtc[kINEL04025] = true;
		}

		Double1D mulGeo(kCentClassBinEnd, 0);
		for (Int_t i = 0; i < nPrim; i++)
		{
			TParticle *part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE)
				continue;
			eta = part->Eta();
			phi = part->Phi();
			pt = part->Pt();
			if (bevtc[kINEL] && fabs(eta) < 1.0)
				bevtc[kINELg0] = true;
			if (abs(eta) < 1)
				Nch_INELg0_eta1++;
			fHistos->FillTH1("hdndetamc", eta);
			fHistos->FillTH1("hkinept", pt, 1. / pt);
			if (abs(eta) < 2)
				nspdgeo++;
			if ((eta > -3.7 && eta < -1.7) || (eta > 2.8 && eta < 5.1))
				nv0mgeo++;
		}
		Double1D geocent(kCentClassBinEnd, -1);
		if (abs(vtxMC3d[2]) < 10 && bevtc[kINELg0])
			fHistos->FillTH1("hEtTruth", fCent[kV0M]);
	}
	else
		bevtc[kDATA] = true;

	if (fOption.Contains("MC") && bevtc[kINELg0] && fabs(vtxMC3d[2]) < 10)
	{
		fHistos->FillTH1("hEventNumbers", "INELg010", 1);
		fHistos->FillTH1("hINELg0Nch", Nch_INELg0_eta1, 1);
	}
	if (fOption.Contains("MC") && bevtc[kINELg0])
		fHistos->FillTH1("hEventNumbers", "INELg0", 1);

	AliInputEventHandler *inputHandler = (AliInputEventHandler *)
											 AliAnalysisManager::GetAnalysisManager()
												 ->GetInputEventHandler();
	Bool_t IsMinimumBias = kFALSE;
	fHistos->FillTH1("hEventNumbers", "All", 1);
	Bool_1d btrigc(kTrigend, false);

	IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kINT7;

	if (!fOption.Contains("MC"))
		btrigc[kHMMBAND] = inputHandler->IsEventSelected() & AliVEvent::kHighMultV0;
	else
		btrigc[kHMMBAND] = IsMinimumBias;

	if (fOption.Contains("LHC10") || fOption.Contains("MB"))
	{
		IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kMB;
	}
	//if (!fOption.Contains("MC"))
	//	IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kHighMultV0;

	if (IsMinimumBias)
		fHistos->FillTH1("hEventNumbers", "PS", 1);

	const AliVVertex *trackVtx = fEvt->GetPrimaryVertexTracks();
	const AliVVertex *spdVtx = fEvt->GetPrimaryVertexSPD();

	Bool_t IsGoodVertex = kFALSE;
	Bool_t IsGoodVertexCut = kFALSE;

	double trackz = -50;
	fZ = -50;

	float vtxf[3] = {-50, -50, -50};
	if (spdVtx)
	{
		fZ = spdVtx->GetZ();
		vtxf[0] = spdVtx->GetX();
		vtxf[1] = spdVtx->GetY();
		vtxf[2] = spdVtx->GetZ();
		IsGoodVertex = kTRUE;
		if (spdVtx->GetNContributors() < 1)
			IsGoodVertex = kFALSE;
		else
			IsGoodVertex = true;
	}
	else
		IsGoodVertex = kFALSE;

	if (sel && !sel->GetThisEventHasNoInconsistentVertices())
		IsGoodVertex = false;
	//if (sel && !sel -> GetThisEventHasGoodVertex2016 () && !fOption.Contains("LHC10d")) IsGoodVertex = false;
	if (sel && !sel->GetThisEventHasGoodVertex2016())
		IsGoodVertex = false;

	if (IsMinimumBias && IsGoodVertex)
		fHistos->FillTH1("hEventNumbers", "Goodz", 1);
	if (IsGoodVertex)
		fHistos->FillTH1("zdata", vtxf[2]);

	if (IsGoodVertex && fabs(fZ) < 10.)
	{
		IsGoodVertexCut = kTRUE;
		if (IsMinimumBias)
		{
			fHistos->FillTH1("hEventNumbers", "Goodzcut", 1);
			if (fabs(fZ) < 7)
				fHistos->FillTH1("hEventNumbers", "Goodzcut7", 1);
		}
	}

	if (IsGoodVertexCut)
	{
		const int nESDTracks = fEvt->GetNumberOfTracks();
		Int_t tempmul = 0;
		for (auto it = 0; it < fEvt->GetNumberOfTracks(); it++)
		{
			AliESDtrack *track = (AliESDtrack *)fEvt->GetTrack(it);
			if (!track)
				continue;
			

			if ( fTrackCuts[kITSTPC2011].AcceptTrack(track)){
				if (fabs(track->Eta()) < 0.5)
					tempmul++;
			}

			if (!fTrackCuts[kHybrid].AcceptTrack(track) && !fTrackCutGC.AcceptTrack(track))
				continue;
			//if (!fTrackCuts[kITSTPC2011].AcceptTrack(track))
			//	continue;
			eta = track->Eta();
			pt = track->Pt();
			phi = track->Phi();
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 0.15)
			{
				btrigc[kMBAND015] = true;
			}
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 0.5)
				btrigc[kMBAND050] = true;
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 0.4)
				btrigc[kMBAND040] = true;
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 1)
				btrigc[kMBAND100] = true;
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 2)
				btrigc[kMBAND200] = true;
			if (IsMinimumBias && fabs(eta) < 0.8 && pt > 4)  // To Chiara: This is OK because it is the trigger condition you are talking about
				btrigc[kMBAND300] = true;
			if (btrigc[kHMMBAND] && fabs(eta) < 0.8 && pt > 4)
				btrigc[kHMMBAND300] = true;
		}
		if ( IsMinimumBias ){
		  if (fCent[kV0M]>0 && fCent[kV0M]<0.01) {
		    	if (!fOption.Contains("MC")) FillTHnSparse("hrecmult", {Double_t(tempmul)});
		  }
		}

		for (auto it = 0; it < fEvt->GetNumberOfTracks(); it++)
		{
			AliESDtrack *track = (AliESDtrack *)fEvt->GetTrack(it);
			if (!track)
				continue;
			

			if ( fTrackCuts[kITSTPC2011].AcceptTrack(track)){
				if (fabs(track->Eta()) < 0.5)
					tempmul++;
			}

			if (!fTrackCuts[kHybrid].AcceptTrack(track) && !fTrackCutGC.AcceptTrack(track))
				continue;
			//if (!fTrackCuts[kITSTPC2011].AcceptTrack(track))
			//	continue;
			eta = track->Eta();
			pt = track->Pt();
			phi = track->Phi();

			if ( btrigc[kMBAND200] && fTrackCuts[kITSTPC2011].AcceptTrack(track)){
			  if (!fOption.Contains("MC"))	FillTHnSparse("hrecdndpt", {track->Eta(), track->Pt()});
			}
		}


		this->FillTracks(bevtc, btrigc);
		this->StrangenessMeasure(btrigc);
	}

	//if (IsMinimumBias ){
	for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
	{
	  if (btrigc[itrigc])
	    {
			for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
			{
				FillTHnSparse("heventcount", {Double_t(itrigc), fCent[icentc], Double_t(icentc)});
			}
	    }
	}
	//}

	if (mcEvent)
	{
		for (auto ievtc = 1u; ievtc < kECend; ievtc++)
		{
			if (bevtc[ievtc])
				for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
				{
					FillTHnSparse("hkinezvtx", {Double_t(ievtc), fCent[icentc], vtxMC3d[2], Double_t(icentc)}, sdweightingfactor);
				}
		}
		Int_t nPrim = stack->GetNprimary();
		Int_t pid;
		Int_t label;
		std::vector<Int_t> labels;

		Int_t truetempmul = 0;
		Int_t mctempmul = 0;

		for (auto it = 0; it < fEvt->GetNumberOfTracks(); it++)
		{
			AliESDtrack *track = (AliESDtrack *)fEvt->GetTrack(it);
			if (!track)
				continue;
			if (!fTrackCuts[kITSTPC2011].AcceptTrack(track))
				continue;
			eta = track->Eta();
			phi = track->Phi();
			pt = track->Pt();
			label = fabs(track->GetLabel());
			labels.push_back(label);
			if (fabs(eta) < 0.5)
				mctempmul++;

			TParticle *particle = stack->Particle( (label<nPrim) ? label : 0);
			if (particle)
			{
				TParticle *mother = AliPWG0Helper::FindPrimaryMother(stack, TMath::Abs((label<nPrim) ? label : 0));
				Bool_t isbackground = false;
				if (mother && IsMotherStrangeParticle(mother))
					isbackground = true;
				if (!AliPWG0Helper::IsPrimaryCharged(particle, nPrim))
					isbackground = true;

				if (!isbackground)
				{
					if (IsGoodVertexCut && fabs(vtxMC3d[2]) < 10 && btrigc[kMBAND200] && bevtc[kINEL200] && pt > 2 && particle->Pt() > 2)
					  if (!fOption.Contains("MC")) FillTHnSparse("hdndptres", {eta, pt, particle->Pt()});
				}
				else
				{

					if (IsGoodVertexCut && fabs(vtxMC3d[2]) < 10 && !bevtc[kINEL200] && btrigc[kMBAND200] && pt > 2 && particle->Pt() > 2)
					  if (!fOption.Contains("MC")) FillTHnSparse("hdndptfake", {eta, pt});
				}
				if ((!IsGoodVertexCut || !IsMinimumBias) && fabs(vtxMC3d[2]) < 10 && bevtc[kINEL200] && !btrigc[kMBAND200] && pt > 2 && particle->Pt() > 2)
				  if (!fOption.Contains("MC"))	FillTHnSparse("hdndptmiss", {eta, particle->Pt()});

			}
		}

		for (Int_t i = 0; i < nPrim; i++)
		{
			TParticle *part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE)
				continue;
			eta = part->Eta();
			phi = part->Phi();
			pt = part->Pt();
			if (fabs(eta) < 0.5)
				truetempmul++;
			switch (TMath::Abs(part->GetPdgCode()))
			{
			case 211:
				pid = kPion;
				break;
			case 321:
				pid = kKaon;
				break;
			case 2212:
				pid = kProtonBK;
				break;
			default:
				pid = kOPar;
				break;
			}
			if (fabs(vtxMC3d[2])<10 && bevtc[kINEL200] && pt>2){
			  if (!fOption.Contains("MC"))	FillTHnSparse("htruedndpt", {eta, pt});
				if (std::find(labels.begin(), labels.end(), i) == labels.end())
				{
				  if (!fOption.Contains("MC"))	FillTHnSparse("hdndptmiss", {eta, pt});
				}
			} 


			for (auto ievtc = 1u; ievtc <= kECend; ievtc++)
			{
				Bool_t IsPtSelected = false;
				if (ievtc <= kINELg0)
					IsPtSelected = true;
				else if (ievtc == kINEL015 && pt > 0.15)
					IsPtSelected = true;
				else if (ievtc == kINEL040 && pt > 0.40)
					IsPtSelected = true;
				else if (ievtc == kINEL050 && pt > 0.5)
					IsPtSelected = true;
				//else if (ievtc == kINEL100 && pt > 1)
				else if (ievtc == kINEL100)
					IsPtSelected = true;
				//else if (ievtc == kINEL200 && pt > 2)
				else if (ievtc == kINEL200)
					IsPtSelected = true;
				//else if (ievtc == kINEL300 && pt > 3)
				//	IsPtSelected = true;
				else if (ievtc == kINEL300) // To Chiara: Please find the modification from the commented-out previous condition. You want to fill all charged particles with the event class
					IsPtSelected = true;
				else if (ievtc == kINEL05024 && pt > 0.5)
					IsPtSelected = true;
				else if (ievtc == kINEL05025 && pt > 0.5)
					IsPtSelected = true;
				else if (ievtc == kINEL04025 && pt > 0.4)
					IsPtSelected = true;

				if (bevtc[ievtc])
				{
					if (!IsPtSelected)
						continue;
					for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
					{
						FillTHnSparse("hkinedndeta",
									  {Double_t(ievtc), fCent[icentc], vtxMC3d[2], Double_t(pid), double(kNoTrkCutVar), eta, Double_t(icentc), phi}, sdweightingfactor);
					}
				}
			}
		}

		if (IsGoodVertexCut && fabs(vtxMC3d[2]) < 10 && IsMinimumBias && fCent[kV0M]>0 && fCent[kV0M]<0.01)
		  if (!fOption.Contains("MC")) FillTHnSparse("hmultres", {Double_t(mctempmul), Double_t(truetempmul)});
		if ((!IsGoodVertexCut || !IsMinimumBias) && fabs(vtxMC3d[2]) < 10 && fCent[kV0M]>0 && fCent[kV0M]<0.01)
		  if (!fOption.Contains("MC")) FillTHnSparse("hmultmiss", {Double_t(truetempmul)});
		if (fabs(vtxMC3d[2]) < 10 && fCent[kV0M]>0 && fCent[kV0M]<0.01)
		  if (!fOption.Contains("MC")) FillTHnSparse("htruemult", {Double_t(truetempmul)});
	}

	if (IsMinimumBias)
	{
		if (fOption.Contains("MC") && bevtc[kINELg0])
			fHistos->FillTH1("hEventNumbers", "PSINELg0", 1);
	}

	PostData(1, fHistos->GetListOfHistograms());
}


void AliAnalysisPseudoRapidityDensityTemp::FillTracks(Bool_1d bevtc, Bool_1d btrigc)
{

	const int nESDTracks = fEvt->GetNumberOfTracks();
	Double_t eta, phi, pt = -1.;
	Int_t label;
	Int_t pid;

	for (auto icut = UInt_t(kHybrid); icut < kTrackCutend; icut++)
	{
		for (auto it = 0; it < fEvt->GetNumberOfTracks(); it++)
		{
			AliESDtrack *track = (AliESDtrack *)fEvt->GetTrack(it);
			if (!track)
				continue;
			eta = track->Eta();
			phi = track->Phi();
			pt = track->Pt();
			label = fabs(track->GetLabel());
			if (icut == kHybrid)
			{
				if (!fTrackCuts[kHybrid].AcceptTrack(track) && !fTrackCutGC.AcceptTrack(track))
					continue;
			}
			else
			{
				if (!fTrackCuts[icut].AcceptTrack(track))
					continue;
			}
			//if (icut == kHybrid)
			//	fHistos->FillTH2("hPhiEtaHybrid", phi, eta);

			if (ismc)
			{
				if (label == 0)
					continue;
				if (stack->GetNtrack() < fabs(label))
					continue;
				TParticle *particle = stack->Particle(TMath::Abs(label));
				if (!particle)
					continue;

				switch (TMath::Abs(particle->GetPdgCode()))
				{
				case 211:
					pid = kPion;
					break;
				case 321:
					pid = kKaon;
					break;
				case 2212:
					pid = kProtonBK;
					break;
				default:
					pid = kOPar;
					break;
				}
				TParticle *mother = AliPWG0Helper::FindPrimaryMother(stack, TMath::Abs(label));
				if (IsMotherStrangeParticle(mother))
					pid = kMotherStrange;
				//else if (!stack->IsPhysicalPrimary(TMath::Abs(label)))
				else if (!AliPWG0Helper::IsPrimaryCharged(particle, stack->GetNprimary()))
					pid = kBkg;

				eta = particle->Eta();
				phi = particle->Phi();
				pt = particle->Pt();
			}
			else
			{
				eta = track->Eta();
				phi = track->Phi();
				pt = track->Pt();
				pid = kParDATA;
			}
			for (auto ievtc = 1u; ievtc < kECend; ievtc++)
			{
				for (auto itrigc = UInt_t(kHMMBAND); itrigc < kTrigend; itrigc++)
				{
					if (bevtc[ievtc] && btrigc[itrigc])
					{
						Bool_t IsPtSelected = false;
						if (itrigc == kMBAND015 && pt > 0.15) {
							IsPtSelected = true;
							if (icut == kHybrid && fabs(eta)<0.8)
								fHistos->FillTH1("hPhiHybrid", phi);
						}
						else if (itrigc == kMBAND050 && pt > 0.5){
							IsPtSelected = true;
						}
						else if (itrigc == kMBAND040 && pt > 0.4) {
							IsPtSelected = true;
						}
						//else if (itrigc == kMBAND100 && pt > 1)
						else if (itrigc == kMBAND100)
							IsPtSelected = true;
						//else if (itrigc == kMBAND200 && pt > 2)
						else if (itrigc == kMBAND200)
							IsPtSelected = true;
						//else if (itrigc == kMBAND300 && pt > 3)
						//	IsPtSelected = true;
						else if (itrigc == kMBAND300)  // To Chiara: Please find the modification from the commented-out previous lines. You want to fill all charged tracks with the trigger class
							IsPtSelected = true;
						//else if (itrigc == kHMMBAND300 && pt > 3)
						else if (itrigc == kHMMBAND300){
							IsPtSelected = true;
						}
						else if (itrigc == kHMMBAND) {//I added this (Chiara)
						  IsPtSelected = true;
					        }
						for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
						{
							if (!IsPtSelected)
								continue;
							FillTHnSparse("hrecdndeta",
										  {Double_t(ievtc), Double_t(itrigc), fCent[icentc], fZ, Double_t(pid), double(icut), eta, Double_t(icentc), phi}, sdweightingfactor);
						}
					}
				}
			} // ievtc loop
		}	  // it loop
	}		  // icut loop
	for (auto ievtc = 1u; ievtc < kECend; ievtc++)
	{
		for (auto itrigc = UInt_t(kHMMBAND); itrigc < kTrigend; itrigc++)
		{
			if (bevtc[ievtc] && btrigc[itrigc])
			{
				for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
				{
					FillTHnSparse("hreczvtx", {Double_t(ievtc), Double_t(itrigc), fCent[icentc], fZ, Double_t(icentc)}, sdweightingfactor);
				}
			}
		}
	}
}

void AliAnalysisPseudoRapidityDensityTemp::FillFMDAODtracks()
{

	AliAODForwardMult *fForward = (AliAODForwardMult *)fEvt->FindListObject("Forward");
	TH2D *hist2d = &(fForward->GetHistogram());
	fForward2d->Add(hist2d);
}

Bool_t AliAnalysisPseudoRapidityDensityTemp::FillFMDESDtracks()
{
	AliFMDEventInspector *fEventInspector = new AliFMDEventInspector("event");
	TVector3 ip;
	UInt_t triggers = 0;
	Bool_t lowFlux = kFALSE;
	UShort_t ivz = 0;
	Double_t cent = -1;
	UShort_t nClusters = 0;

	UInt_t found = fEventInspector->Process(dynamic_cast<AliESDEvent *>(fEvt), triggers, lowFlux, ivz, ip, cent, nClusters);
	if (found & AliFMDEventInspector::kNoEvent)
	{
		return false;
	}
	if (found & AliFMDEventInspector::kNoTriggers)
	{
		return false;
	}
	return true;
}

//___________________________________________________________________
//void AliAnalysisPseudoRapidityDensityTemp::FinishTaskOutput()
//{
//OpenFile(1);
//TH1D *fForward1d = (TH1D*)fForward2d->ProjectionX("_x",1,-1,"e");
//TH1D* norm   = fForward2d->ProjectionX("norm", 0, 1, "");
//fForward1d -> Divide(norm);
//fForward1d -> Write();
//}
//___________________________________________________________________
void AliAnalysisPseudoRapidityDensityTemp::Terminate(Option_t *)
{
}
//___________________________________________________________________

Int_t AliAnalysisPseudoRapidityDensityTemp::GetPID(AliPIDResponse *pid, const AliVTrack *trk)
{
	if (!pid)
		return -1; // no pid available

	Double_t prob[AliPID::kSPECIES];
	fPIDCombined->ComputeProbabilities(trk, pid, prob);
	Int_t ipid = -1;
	Double_t iprob = 0;
	for (int i = 0; i < AliPID::kSPECIES; i++)
	{
		if (prob[i] > iprob)
		{
			iprob = prob[i];
			ipid = i;
		}
	}

	return ipid;
}

Bool_t AliAnalysisPseudoRapidityDensityTemp::IsMCTRDtriggered()
{
	Bool_t IsTRD = kFALSE;
	Int_t nTrdTracks = fEvt->GetNumberOfTrdTracks();
	if (nTrdTracks > 0)
	{
		for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack)
		{
			AliESDTrdTrack *trdTrack = (AliESDTrdTrack *)fEvt->GetTrdTrack(iTrack);
			if (!trdTrack)
				continue;

			// simulate HNU
			if ((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) ||
				(trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4))
			{
				IsTRD = kTRUE;
			}

			// simulate HQU
			if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
				trdTrack->GetPID() >= 130 &&
				trdTrack->GetNTracklets() >= 5 &&
				(trdTrack->GetLayerMask() & 1))
			{
				Float_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
				if (sag < 0.2 && sag > -0.2)
				{
					IsTRD = kTRUE;
				}
			}
		}
	}
	return IsTRD;
}

//___________________________________________________________________________________
Float_t AliAnalysisPseudoRapidityDensityTemp::GetInvPtDevFromBC(Int_t b, Int_t c)
{
	//returns d(1/Pt) in c/GeV
	//in case of no gtu simulation -> return maximum 0.5
	if (b == 0 && c == 0)
		return 0.5;
	Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
	tmp += (c & 0xfff);
	Float_t invPtDev = tmp * 0.000001;
	return invPtDev;
}
void AliAnalysisPseudoRapidityDensityTemp::MeasureDiffMass()
{

	Int_t np = stack->GetNprimary();

	Int_t iPart1 = -1;
	Int_t iPart2 = -1;
	Double_t cms = 0;

	if (fOption.Contains("LHC15n"))
		cms = 5020;
	else if (fOption.Contains("LHC10d"))
		cms = 7000;
	Double_t y1 = 1e10;
	Double_t y2 = -1e10;
	for (Int_t i = 0; i < np; ++i)
	{
		TParticle *part = stack->Particle(i);

		Int_t statusCode = part->GetStatusCode();

		// Initial state particle
		if (statusCode != 1)
			continue;

		Int_t pdg = TMath::Abs(part->GetPdgCode());
		Bool_t isStable = kFALSE;
		for (Int_t i1 = 0; i1 < kNstable; i1++)
		{
			if (pdg == pdgStable[i1])
			{
				isStable = kTRUE;
				break;
			}
		}
		if (!isStable)
			continue;

		Double_t y = part->Y();

		if (y < y1)
		{
			y1 = y;
			iPart1 = i;
		}
		if (y > y2)
		{
			y2 = y;
			iPart2 = i;
		}
		if (iPart1 >= 0 && iPart2 >= 0)
		{
			y1 = TMath::Abs(y1);
			y2 = TMath::Abs(y2);

			TParticle *part1 = (TParticle *)stack->Particle(iPart1);
			TParticle *part2 = (TParticle *)stack->Particle(iPart2);

			Int_t pdg1 = part1->GetPdgCode();
			Int_t pdg2 = part2->GetPdgCode();

			Int_t iPart = -1;
			if (pdg1 == 2212 && pdg2 == 2212)
			{
				if (y1 > y2)
					iPart = iPart1;
				else if (y1 < y2)
					iPart = iPart2;
				else
				{
					iPart = iPart1;
					if (fRandom->Uniform(0., 1.) > 0.5)
						iPart = iPart2;
				}
			}
			else if (pdg1 == 2212)
				iPart = iPart1;
			else if (pdg2 == 2212)
				iPart = iPart2;

			Double_t M = -1.;
			if (iPart > 0)
			{
				TParticle *part = (TParticle *)stack->Particle(iPart);
				Double_t E = part->Energy();
				Double_t P = part->P();
				//energy 13 TeV = 13000
				Double_t M2 = (cms - E - P) * (cms - E + P);
				if (M2 > 0)
				{
					M = TMath::Sqrt(M2);
					if (fOption.Contains("DiffTune") && hdifftuneratio->GetBinContent(hdifftuneratio->GetXaxis()->FindBin(M)))
						sdweightingfactor = hdifftuneratio->Interpolate(M);
					if (M > 200)
						sdweightingfactor = 0;
					//sdweightingfactor = hdifftuneratio->GetBinContent(hdifftuneratio->GetXaxis()->FindBin(M));
					fHistos->FillTH1("hsdmass", M, sdweightingfactor);
				}
			}
		} // end of iPart1>=0 && iPart2 >=0
	}
}

void AliAnalysisPseudoRapidityDensityTemp::StrangenessMeasure(Bool_1d btrigc)
{
	Double_t tPrimaryVtxPosition[3];
	Int_t nv0s = 0;
	AliESDEvent *evt = dynamic_cast<AliESDEvent *>(fEvt);
	nv0s = evt->GetNumberOfV0s();
	Int_t lOnFlyStatus = 0, nv0sOn = 0, nv0sOff = 0;
	Double_t lChi2V0 = 0;
	Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
	Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
	Double_t lV0CosineOfPointingAngle = 0;
	Double_t lV0Radius = 0;
	Double_t lV0DecayLength = 0;
	Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
	Double_t lPt = 0, lRapK0s = 0, lRapLambda = 0;
	Double_t lAlphaV0 = 0, lPtArmV0 = 0;
	Double_t tV0Position[3];
	Double_t lMagneticField = 999;
	const AliESDVertex *primaryVtx = evt->GetPrimaryVertex();
	tPrimaryVtxPosition[0] = primaryVtx->GetX();
	tPrimaryVtxPosition[1] = primaryVtx->GetY();
	tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	if (TMath::Abs(tPrimaryVtxPosition[2]) > 10)
		return;

	for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
	{
		AliESDv0 *v0 = evt->GetV0(iV0);
		UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());
		lMagneticField = ((AliESDEvent *)evt)->GetMagneticField();

		AliESDtrack *pTrack = evt->GetTrack(lKeyPos);
		AliESDtrack *nTrack = evt->GetTrack(lKeyNeg);
		if (!pTrack || !nTrack)
		{
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->GetSign() == nTrack->GetSign())
		{
			//cout<< "like sign, continue"<< endl;
			continue;
		}
		// Tracks quality cuts
		if (((pTrack->GetTPCNcls()) < 80) || ((nTrack->GetTPCNcls()) < 80))
			continue;

		// TPC refit condition (done during reconstruction for
		// Offline but not for On-the-fly)
		if (!(pTrack->GetStatus() & AliESDtrack::kTPCrefit))
			continue;
		if (!(nTrack->GetStatus() & AliESDtrack::kTPCrefit))
			continue;

		if (pTrack)
			lDcaPosToPrimVertex =
				TMath::Abs(pTrack->GetD(
					tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], lMagneticField));

		if (nTrack)
			lDcaNegToPrimVertex =
				TMath::Abs(nTrack->GetD(
					tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], lMagneticField));

		lOnFlyStatus = v0->GetOnFlyStatus();
		lChi2V0 = v0->GetChi2V0();
		lDcaV0Daughters = v0->GetDcaV0Daughters();
		lDcaV0ToPrimVertex = v0->GetD(
			tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], tPrimaryVtxPosition[2]);

		lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(
			tPrimaryVtxPosition[0], tPrimaryVtxPosition[1], tPrimaryVtxPosition[2]);
		v0->GetXYZ(tV0Position[0], tV0Position[1], tV0Position[2]);
		lV0Radius = TMath::Sqrt(
			tV0Position[0] * tV0Position[0] + tV0Position[1] * tV0Position[1]);
		lV0DecayLength = TMath::Sqrt(
			TMath::Power(tV0Position[0] - tPrimaryVtxPosition[0], 2) +
			TMath::Power(tV0Position[1] - tPrimaryVtxPosition[1], 2) +
			TMath::Power(tV0Position[2] - tPrimaryVtxPosition[2], 2));

		// Pt:
		lPt = v0->Pt();

		// Armenteros variables: !!
		lAlphaV0 = v0->AlphaV0();
		lPtArmV0 = v0->PtArmV0();

		// Selections:
		if (1)
		{
			if ((lDcaPosToPrimVertex < 0.05) ||
				(lDcaNegToPrimVertex < 0.05) ||
				(lDcaV0Daughters > 0.5) ||
				(lV0CosineOfPointingAngle < 0.99))
				continue;
		}

		for (auto icentc = UInt_t(kCentClassBinBegin); icentc < kCentClassBinEnd; icentc++)
		{
			v0->ChangeMassHypothesis(310);
			for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
				if (btrigc[itrigc])
				  if (!fOption.Contains("MC"))	FillTHnSparse("hv0mass",
								  {Double_t(itrigc), fCent[icentc], Double_t(kK0s), v0->GetEffMass(), Double_t(icentc)});
			if (v0->GetEffMass() > 0.482 && v0->GetEffMass() < 0.509)
			{
				for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
					if (btrigc[itrigc])
					{

						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kK0s), pTrack->Eta(), Double_t(icentc)});
						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kK0s), nTrack->Eta(), Double_t(icentc)});
					}
			}
			v0->ChangeMassHypothesis(3122);
			for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
				if (btrigc[itrigc])
				  if (!fOption.Contains("MC"))		FillTHnSparse("hv0mass",
								  {Double_t(itrigc), fCent[icentc], Double_t(kLambda), v0->GetEffMass(), Double_t(icentc)});
			if (v0->GetEffMass() > 1.11 && v0->GetEffMass() < 1.12)
			{
				for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
					if (btrigc[itrigc])
					{
						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kLambda), pTrack->Eta(), Double_t(icentc)});
						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kLambda), nTrack->Eta(), Double_t(icentc)});
					}
			}

			v0->ChangeMassHypothesis(-3122);
			for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
				if (btrigc[itrigc])
				  if (!fOption.Contains("MC")) FillTHnSparse("hv0mass",
								  {Double_t(itrigc), fCent[icentc], Double_t(kAntiLambda), v0->GetEffMass(), Double_t(icentc)});
			if (v0->GetEffMass() > 1.11 && v0->GetEffMass() < 1.12)
			{
				for (auto itrigc = 1u; itrigc < kTrigend; itrigc++)
					if (btrigc[itrigc])
					{

						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kAntiLambda), pTrack->Eta(), Double_t(icentc)});
						FillTHnSparse("hv0eta",
									  {Double_t(itrigc), fCent[icentc], Double_t(kAntiLambda), nTrack->Eta(), Double_t(icentc)});
					}
			}
		}
	}
}

THnSparse *AliAnalysisPseudoRapidityDensityTemp::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt)
{
	const TAxis *axises[bins.size()];
	for (UInt_t i = 0; i < bins.size(); i++)
		axises[i] = &bins[i];
	THnSparse *h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
	return h;
}

THnSparse *AliAnalysisPseudoRapidityDensityTemp::CreateTHnSparse(TString name, TString title, TString templ, Option_t *opt)
{
	auto o = fHistos->FindObject(templ);
	if (!o)
	{
		cout << "ERROR: no " << templ << endl;
		gSystem->Exit(1);
	}
	auto ht = dynamic_cast<THnSparse *>(o);
	const TAxis *axises[ht->GetNdimensions()];
	for (int i = 0; i < ht->GetNdimensions(); i++)
		axises[i] = ht->GetAxis(i);
	auto h = fHistos->CreateTHnSparse(name, title, ht->GetNdimensions(), axises, opt);
	return h;
}

Long64_t AliAnalysisPseudoRapidityDensityTemp::FillTHnSparse(TString name, std::vector<Double_t> x, Double_t w)
{
	auto hsparse = dynamic_cast<THnSparse *>(fHistos->FindObject(name));
	if (!hsparse)
	{
		cout << "ERROR : no " << name << endl;
		exit(1);
	}
	return FillTHnSparse(hsparse, x, w);
}

Long64_t AliAnalysisPseudoRapidityDensityTemp::FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w)
{
	if (int(x.size()) != h->GetNdimensions())
	{
		cout << "ERROR : wrong sized of array while Fill " << h->GetName() << endl;
		exit(1);
	}
	return h->Fill(&x.front(), w);
}

TAxis AliAnalysisPseudoRapidityDensityTemp::AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax)
{
	TAxis axis(nbin, xmin, xmax);
	axis.SetName(name);
	return axis;
}

TAxis AliAnalysisPseudoRapidityDensityTemp::AxisStr(TString name, std::vector<TString> bin)
{
	TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
	UInt_t i = 1;
	for (auto blabel : bin)
		ax.SetBinLabel(i++, blabel);
	return ax;
}

TAxis AliAnalysisPseudoRapidityDensityTemp::AxisVar(TString name, std::vector<Double_t> bin)
{
	TAxis axis(bin.size() - 1, &bin.front());
	axis.SetName(name);
	return axis;
}

TAxis AliAnalysisPseudoRapidityDensityTemp::AxisLog(TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0)
{
	int binoffset = (xmin0 < 0 || (xmin - xmin0) < 1e-9) ? 0 : 1;
	std::vector<Double_t> bin(nbin + 1 + binoffset, 0);
	double logBW3 = (log(xmax) - log(xmin)) / nbin;
	for (int ij = 0; ij <= nbin; ij++)
		bin[ij + binoffset] = xmin * exp(ij * logBW3);
	TAxis axis(nbin, &bin.front());
	axis.SetName(name);
	return axis;
}

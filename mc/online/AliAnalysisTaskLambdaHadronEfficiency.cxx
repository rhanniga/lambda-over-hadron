#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TStopwatch.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliCFParticle.h"

#include "AliEventPoolManager.h"
#include "AliMultSelection.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskLambdaHadronEfficiency.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskLambdaHadronEfficiency)
//________________________________________________________________________
AliAnalysisTaskLambdaHadronEfficiency::AliAnalysisTaskLambdaHadronEfficiency(const char *name, Float_t multLow, Float_t multHigh)
: AliAnalysisTaskSE(name),
fVevent(0),
fPoolMgr(0x0),
fLSPoolMgr(0x0),
fHHPoolMgr(0x0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fNumTracks(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fRealKDist(0),
fRealPiDist(0),
fRealpDist(0),
fRealeDist(0),
fRealMuonDist(0),
fRealLambdaDist(0),
fRecoLambdaDist(0)
{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
    //printf("\n\n!!!!!!!!!!]\n done with the constructor! \n");
    //fflush(stdout);
    MIN_CROSSED_ROWS_TPC = 70;
    MIN_ROW_CLUSTER_RATIO = 0.8;
    MULT_LOW = multLow;
    MULT_HIGH = multHigh;
    CENT_ESTIMATOR = "V0A";
    // DAUGHTER_TRK_BIT = 16; not used
    ASSOC_TRK_BIT = 1024;
    TRIG_TRK_BIT = 768;
    DAUGHTER_ETA_CUT = 0.8;
    DAUGHTER_MIN_PT = 0.15;

}
//________________________________________________________________________
AliAnalysisTaskLambdaHadronEfficiency::AliAnalysisTaskLambdaHadronEfficiency()
: AliAnalysisTaskSE("default_task"),
fVevent(0),
fPoolMgr(0x0),
fLSPoolMgr(0x0),
fHHPoolMgr(0x0),
fESD(0),
fAOD(0),
fpidResponse(0),
fOutputList(0),
fNevents(0),
fNumTracks(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fRealLambdaDist(0),
fRecoLambdaDist(0)
{
    //Default constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
    MIN_CROSSED_ROWS_TPC = 70;
    MIN_ROW_CLUSTER_RATIO = 0.8;
    MULT_LOW = 0;
    MULT_HIGH = 80;
    CENT_ESTIMATOR = "V0A";
    // DAUGHTER_TRK_BIT = 16; not used
    ASSOC_TRK_BIT = 1024;
    TRIG_TRK_BIT = 768;
    DAUGHTER_ETA_CUT = 0.8;
    DAUGHTER_MIN_PT = 0.15;
}
//________________________________________________________________________
AliAnalysisTaskLambdaHadronEfficiency::~AliAnalysisTaskLambdaHadronEfficiency()
{
    //Destructor
    delete fOutputList;
}
//________________________________________________________________________
void AliAnalysisTaskLambdaHadronEfficiency::UserCreateOutputObjects()
{
    //printf("\n!!!!!\n Starting UserCreateOutputObjects \n\n");
    //fflush(stdout);
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
   
    ////////////////////////////
    // Set-up for Mixed Event //
    ////////////////////////////


    Int_t numVtxZBins = 10;
    //Double_t vtxZBins[11] = {-10.0, -6.15, -3.90, -2.13, -0.59, 0.86, 2.29, 3.77, 5.39, 7.30, 10.0};
    Double_t vtxZBins[11] = {-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
    Int_t numMultBins = 3;
    Double_t multBins[4] = {0.0, 20.0, 50.0, 90.0};


    ////////////////
    //Output list//
    ///////////////
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");

    fNumTracks = new TH1F("fNumTracks", "Number of Tracks/evt", 1000, 0, 1000);
    fOutputList->Add(fNumTracks);
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);
    
    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
    fOutputList->Add(fVtxX);
    
    // Single Particle and inclusive charged hadron histos
    Int_t numbinsSingle[5] = {100, 64, 64, 10, 10};
    Double_t minvalSingle[5] = {0, -3.14159, -2.0, -10.0, 0.0};
    Double_t maxvalSingle[5] = {10.0, 3.14159, 2.0, 10.0, 100.0};

    fRealChargedDist = new THnSparseF("fRealChargedDist", "Real Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealChargedDist);

    fRealKDist = new THnSparseF("fRealKDist", "Real Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealKDist);

    fRealPiDist = new THnSparseF("fRealPiDist", "Real Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealPiDist);

    fRealPiFromLambdaDist = new THnSparseF("fRealPiFromLambdaDist", "Real Pion (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealPiFromLambdaDist);

    fRealpDist = new THnSparseF("fRealpDist", "Real proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealpDist);

    fRealpFromLambdaDist = new THnSparseF("fRealpFromLambdaDist", "Real proton (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealpFromLambdaDist);

    fRealeDist = new THnSparseF("fRealeDist", "Real electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealeDist);

    fRealMuonDist = new THnSparseF("fRealMuonDist", "Real Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRealMuonDist);

    fRecoChargedDist = new THnSparseF("fRecoChargedDist", "Reco Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoChargedDist);

    fRecoKDist = new THnSparseF("fRecoKDist", "Reco Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoKDist);

    fRecoPiDist = new THnSparseF("fRecoPiDist", "Reco Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiDist);

    fRecopDist = new THnSparseF("fRecopDist", "Reco proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecopDist);

    fRecoeDist = new THnSparseF("fRecoeDist", "Reco electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoeDist);

    fRecoMuonDist = new THnSparseF("fRecoMuonDist", "Reco Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoMuonDist);

    fRecoChargedTriggerDist = new THnSparseF("fRecoChargedTriggerDist", "Reco Charged Hadron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoChargedTriggerDist);

    fRecoKTriggerDist = new THnSparseF("fRecoKTriggerDist", "Reco Kaon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoKTriggerDist);
    
    fRecoPiTriggerDist = new THnSparseF("fRecoPiTriggerDist", "Reco Pion Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiTriggerDist);

    fRecopTriggerDist = new THnSparseF("fRecopTriggerDist", "Reco proton Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecopTriggerDist);

    fRecoeTriggerDist = new THnSparseF("fRecoeTriggerDist", "Reco electron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoeTriggerDist);

    fRecoMuonTriggerDist = new THnSparseF("fRecoMuonTriggerDist", "Reco Muon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoMuonTriggerDist);

    // Lambda ditributions used for eff calc
    
    Int_t numbins[6] = {100, 64, 64, 10, 80, 10};
    Double_t minval[6] = {0, -3.14159, -2., -10, 1.06, 0.0};
    Double_t maxval[6] = {10.0, 3.14159,  2.,  10, 1.16, 100.0};

    fRealTotalLambdaDist = new THnSparseF("fRealTotalLambdaDist", "Real (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealTotalLambdaDist);

    fRealLambdaDist = new THnSparseF("fRealLambdaDist", "Real #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealLambdaDist);

    fRealAntiLambdaDist = new THnSparseF("fRealAntiLambdaDist", "Real #bar{#Lambda} distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealAntiLambdaDist);

    fRecoTotalV0LambdaDist = new THnSparseF("fRecoTotalV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoTotalV0LambdaDist);

    fRecoEtaV0LambdaDist = new THnSparseF("fRecoEtaV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaV0LambdaDist);

    fRecoEtaPtV0LambdaDist = new THnSparseF("fRecoEtaPtV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtV0LambdaDist);

    fRecoEtaPtRefitV0LambdaDist = new THnSparseF("fRecoEtaPtRefitV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitV0LambdaDist);

    fRecoEtaPtRefitRowsV0LambdaDist = new THnSparseF("fRecoEtaPtRefitRowsV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit + rows cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitRowsV0LambdaDist);

    fRecoEtaPtRefitRowsRatioV0LambdaDist = new THnSparseF("fRecoEtaPtRefitRowsRatioV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitRowsRatioV0LambdaDist);

    fRecoEtaLambdaDist = new THnSparseF("fRecoEtaLambdaDist", "Reco (#Lambda + #bar{#Lambda})  distribution (eta cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaLambdaDist);

    fRecoEtaPtLambdaDist = new THnSparseF("fRecoEtaPtLambdaDist", "Reco (#Lambda + #bar{#Lambda})  distribution (eta + pt cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtLambdaDist);

    fRecoEtaPtRefitLambdaDist = new THnSparseF("fRecoEtaPtRefitLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitLambdaDist);

    fRecoEtaPtRefitRowsLambdaDist = new THnSparseF("fRecoEtaPtRefitRowsLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitRowsLambdaDist);

    fRecoEtaPtRefitRowsRatioLambdaDist = new THnSparseF("fRecoEtaPtRefitRowsRatioLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoEtaPtRefitRowsRatioLambdaDist);

    fRecoTotalLambdaDist = new THnSparseF("fRecoTotalLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoTotalLambdaDist);

    fRecoLambdaDist = new THnSparseF("fRecoLambdaDist", "Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoLambdaDist);

    fRecoAntiLambdaDist = new THnSparseF("fRecoAntiLambdaDist", "Reco #Lambda (anti) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoAntiLambdaDist);

    fReactionPlane = new TH1D("fReactionPlane", "Reaction Plane Angle; Angle", 64, -3.14159, 3.14159);
    fOutputList->Add(fReactionPlane);

    // Same lambda distributions now with info on daughter filter maps (last two are proton filter map, pion filter map)
    Int_t numbinsFilter[8] = {100, 64, 64, 10, 80, 10, 9, 9};
    Double_t minvalFilter[8] = {0, -3.14159, -2., -10, 1.06, 0.0, 0, 0};
    Double_t maxvalFilter[8] = {10.0, 3.14159,  2.,  10, 1.16, 100.0, 9, 9};

    fRecoTotalLambdaFilterDist = new THnSparseF("fRecoTotalLambdaFilterDist", "Reco (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 8, numbinsFilter, minvalFilter, maxvalFilter);
    fOutputList->Add(fRecoTotalLambdaFilterDist);

    fRealLambdasPerEvent = new TH1F("fRealLambdasPerEvent", "Real Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRealLambdasPerEvent);

    fRecoLambdasPerEvent = new TH1F("fRecoLambdasPerEvent", "Reco Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRecoLambdasPerEvent);

    PostData(1,fOutputList);

}

//_______________________________________________________________________
// UInt_t AliAnalysisTaskLambdaHadronEfficiency::PassProtonCuts(AliAODTrack *track, Double_t TPCnSigma, Double_t TOFnSigma){
//     //returns the level of cuts that the track passed
//     //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
//     UInt_t passLevel = 0;
//     Bool_t pass = kTRUE;
//     pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
//     pass = pass && (track->Pt() >= 0.15);
//     pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
//     pass = pass && (track->GetTPCCrossedRows() > 80);
//     if(pass) passLevel |= TRACK_BIT;    
//     if(TOFnSigma != -999){ //check if there is a TOF signal, but don't care what the signal is
//         passLevel |= TOF_HIT_BIT;
//     }
//     if(TMath::Abs(TPCnSigma) <= 2.0){ // check that kaon passed the TPC nsigma cut
//         passLevel |= TPC_PID_BIT;
//     }
//     if(TMath::Abs(TOFnSigma) <= 2.0){ // check that kaon passed TOF nsigma cut
//         passLevel |= TOF_PID_BIT;
//     }

//     if(TMath::Abs(track->Eta()) <= 0.8) passLevel |= ETA_BIT;
//     if(track->Pt() >= 0.15) passLevel |= PT_BIT;
//     if(track->TestFilterMask(ASSOC_TRK_BIT)) passLevel |= MASK_BIT;
//     if(track->GetTPCCrossedRows() > 80) passLevel |= ROWS_BIT;

//     return passLevel;
// }

//_______________________________________________________________________
// UInt_t AliAnalysisTaskLambdaHadronEfficiency::PassPionCuts(AliAODTrack *track, Double_t TPCnSigma, Double_t TOFnSigma){
//     //returns the level of cuts that the track passed
//     //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
//     UInt_t passLevel = 0;
//     Bool_t pass = kTRUE;
//     pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
//     pass = pass && (track->Pt() >= 0.15);
//     pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
//     pass = pass && (track->GetTPCCrossedRows() > 80);
//     if(pass) passLevel |= TRACK_BIT;    
//     if(TOFnSigma != -999){ //check if there is a TOF signal, but don't care what the signal is
//         passLevel |= TOF_HIT_BIT;
//     }
//     if(TMath::Abs(TPCnSigma) <= 3.0){ // check that kaon passed the TPC nsigma cut
//         passLevel |= TPC_PID_BIT;
//     }
//     if(TMath::Abs(TOFnSigma) <= 3.0){ // check that kaon passed TOF nsigma cut
//         passLevel |= TOF_PID_BIT;
//     }

//     if(TMath::Abs(track->Eta()) <= 0.8) passLevel |= ETA_BIT;
//     if(track->Pt() >= 0.15) passLevel |= PT_BIT;
//     if(track->TestFilterMask(ASSOC_TRK_BIT)) passLevel |= MASK_BIT;
//     if(track->GetTPCCrossedRows() > 80) passLevel |= ROWS_BIT;

//     return passLevel;
// }

//_______________________________________________________________________
uint AliAnalysisTaskLambdaHadronEfficiency::PassDaughterCuts(AliAODTrack *track){

    // Reject the negative ID tracks (TPC constrained to PV, etc.)
    if(track->GetID() < 0) return false;

    uint passLevel = 0;

    // Bit set if track passes eta cut
    if(TMath::Abs(track->Eta()) <= DAUGHTER_ETA_CUT) passLevel |= ETA_BIT;

    // Bit set if track passes pt cut
    if(track->Pt() >= DAUGHTER_MIN_PT) passLevel |= PT_BIT;

    // Bit set if track has TPC refit flag enabled
    if(track->IsOn(AliAODTrack::kTPCrefit)) passLevel |= TPC_REFIT_BIT;

    // Bit set if track has more than min number of crossed rows in TPC (set in constructor, usually 70)
    if(track->GetTPCCrossedRows() > MIN_CROSSED_ROWS_TPC) passLevel |= CROSSED_ROWS_BIT;

    // Bit set if the crossed rows/findable clusters is greater than a certain value (set in constructor, usually 0.8)
    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    if(ratio > MIN_ROW_CLUSTER_RATIO) passLevel |= ROW_CLUSTER_RATIO_BIT;

    return passLevel;

}

Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassAssociatedCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestFilterMask(ASSOC_TRK_BIT);

    return pass;
}

Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassTriggerCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestBit(TRIG_TRK_BIT);

    return pass;
}

//_______________________________________________________________________
// Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassHadronCuts(AliAODTrack *track, Bool_t isTrigger){
//     //returns the level of cuts that the track passed
//     //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
//     Bool_t pass = kTRUE;
//     pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
//     pass = pass && (track->Pt() >= 0.15);
//     if(isTrigger){
//         pass = pass && (track->TestBit(TRIG_TRK_BIT));
//     }else{
//         pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
//     }
//     return pass;
// }


//_______________________________________________________________________
// AliAODTrack* AliAnalysisTaskLambdaHadronEfficiency::GetTrackFromID(AliAODEvent* inputEvent, int trackID) {

//     int numTracks = inputEvent->GetNumberOfTracks();
//     bool trackFound = false;
//     AliAODTrack* returnTrack = 0x0;
//     for(int itrack = 0; itrack < numTracks; itrack++) {
//         auto currentTrack = (AliAODTrack*)inputEvent->GetTrack(itrack);
//         int currentTrackID = currentTrack->GetID();
//         if((currentTrackID == trackID) && !trackFound) {
//             trackFound = true;
//             returnTrack = currentTrack;
//         }    
//         else if((currentTrackID == trackID) && trackFound) {
//             std::cout << "OKAY SO THE TRACK IDs ARE NOT UNIQUE??? WTF IS AN ID THEN???????" << std::endl;
//             std::cout << "THE INTERSECTING ID IS: " << trackID << "\n";
//             std::cout << "track1 stats (pt, eta, phi, filtermap):\n";
//             std::cout << "\t" << returnTrack->Pt() << ", " << returnTrack->Eta() << ", " << returnTrack->Phi() << ", " << returnTrack->GetFilterMap() << std::endl;
//             std::cout << "\t" << currentTrack->Pt() << ", " << currentTrack->Eta() << ", " << currentTrack->Phi() << ", " << currentTrack->GetFilterMap() << std::endl;
//         }
//     }
//     return returnTrack;
// }

//________________________________________________________________________
void AliAnalysisTaskLambdaHadronEfficiency::UserExec(Option_t *){

    //masks for the different cut configurations
    uint maskEta = ETA_BIT;
    uint maskEtaPt = ETA_BIT + PT_BIT;
    uint maskEtaPtRefit = ETA_BIT + PT_BIT + TPC_REFIT_BIT;
    uint maskEtaPtRefitRows = ETA_BIT + PT_BIT + TPC_REFIT_BIT + CROSSED_ROWS_BIT;
    uint maskEtaPtRefitRowsRatio = ETA_BIT + PT_BIT + TPC_REFIT_BIT + CROSSED_ROWS_BIT + ROW_CLUSTER_RATIO_BIT;

    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
     
    

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (fAOD) {
        //printf("fAOD available\n");
        //return;
    }else{
        return;
    }

    // fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    // if(fESD) {
    //     std::cout << "ESD FOUND!" << std::endl;
    // }
    // else {
    //     std::cout << "ESD NOT FOUND" << std::endl;
    // }

    ///////////////////
    //PID initialised//
    //////////////////
   fpidResponse = fInputHandler->GetPIDResponse();

    if(!fpidResponse){
        AliFatal("No PID response loaded!");
    }
    
    ////////////////
    //Event vertex//
    ///////////////
    Int_t ntracks = -999;
    ntracks = fVevent->GetNumberOfTracks();
    
    /////////////////
    //trigger check//
    /////////////////
    fVevent->GetFiredTriggerClasses();
    
    Int_t trigger = -1;
    //Multiplicity stuff
    Double_t multPercentile = 10.0;
    if (fAOD){
        AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
        Double_t multiplicity = header->GetRefMultiplicity();

        fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
        if(fMultSelection){
            multPercentile = fMultSelection->GetMultiplicityPercentile(CENT_ESTIMATOR.Data());
        }else{
            return;
        }
        if(multPercentile < MULT_LOW || multPercentile > MULT_HIGH) return;
    }
    
    fNevents->Fill(0); //all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<3)return;
    fNevents->Fill(1); //events with 3 tracks

    Zvertex = pVtx->GetZ();
    Yvertex = pVtx->GetY();
    Xvertex = pVtx->GetX();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);

    fNumTracks->Fill(ntracks);

    ////////////////////
    //event selection//
    ///////////////////
    if(fabs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut

    Double_t distPoint[6] = {0., 0., 0., 0., 0., 0.};
    Double_t singledistPoint[5] = {0., 0., 0., 0., 0.};
    //Loop over all particles in stack to get real phi, looking for phi->KK
    //AliMCEvent *fMCEvent = dynamic_cast<AliMCEvent*>(InputEvent());
    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
        AliError("Array of MC particles not found");
        return;
    }

    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!fMCHeader) {
        AliError("Could not find MC Header in AOD");
        return;
    }

    Double_t RA = fMCHeader->GetReactionPlaneAngle();
    if(RA > TMath::Pi()){
        RA -= 2.0*TMath::Pi();
    }else if(RA < -1.0*TMath::Pi()){
        RA += 2.0*TMath::Pi();
    }
    fReactionPlane->Fill(RA);

    UInt_t negPassCuts = 0;
    UInt_t posPassCuts = 0;

    Int_t negparPDG = 0;
    Int_t posparPDG = 0;
    Float_t recoPx = 0.0;
    Float_t recoPy = 0.0;
    Float_t recoPz = 0.0;
    Float_t recoP = 0.0;
    Float_t recoE = 0.0;
    Float_t recoM = 0.0;
    Float_t recoEta = 0.0;
    Float_t recoY = 0.0;
    Float_t recoPt = 0.0;
    Float_t recoPhi = 0.0;

    Int_t motherIndex = 0;
    Int_t motherPDG = 0;
    Int_t pdgcode = 0;

    int numRealLambdas = 0;
    int numRecoLambdas = 0;

    std::vector<std::vector<int>> trackIDsV0;

    // quick test loop over V0's to determine filter mask range for daughter tracks
    int numV0s = fAOD->GetNumberOfV0s();
    for(int iv0 = 0; iv0 < numV0s; iv0++) {
        AliAODv0 *vZero = fAOD->GetV0(iv0);
        if(!vZero) continue;

        AliAODTrack *ptrack=(AliAODTrack *)vZero->GetDaughter(0);
        AliAODTrack *ntrack=(AliAODTrack *)vZero->GetDaughter(1);
        
        int plabel = ptrack->GetLabel();
        int nlabel = ntrack->GetLabel();

        if(plabel < 0 || nlabel < 0) continue;

        AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(plabel);
        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(nlabel);

        int posPDG = mcpospart->GetPdgCode();
        int negPDG = mcnegpart->GetPdgCode();

        if(!((posPDG == 2212 && negPDG == -211) || (posPDG == 211 && negPDG == -2212))) continue;

        int posmomlabel = mcpospart->GetMother(); 
        int negmomlabel = mcnegpart->GetMother(); 

        if(posmomlabel < 0 || negmomlabel < 0) continue;

        if(posmomlabel != negmomlabel) {
            continue ;
        }

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(posmomlabel);

        if(!(TMath::Abs(mcmother->GetPdgCode()) == 3122)) continue;


        double distPoint[6];
        double recoPx = ntrack->Px() + ptrack->Px();
        double recoPy = ntrack->Py() + ptrack->Py();
        double recoPz = ntrack->Pz() + ptrack->Pz();
        
        double recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
        double recoE = TMath::Sqrt(ntrack->Px()*ntrack->Px() + ntrack->Py()*ntrack->Py() + ntrack->Pz()*ntrack->Pz() + 0.13957*0.13957) + TMath::Sqrt(ptrack->Px()*ptrack->Px() + ptrack->Py()*ptrack->Py() + ptrack->Pz()*ptrack->Pz() + 0.9383*0.9383);
        double recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
        double recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
        double recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
        double recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
        double recoPhi = TMath::ATan2(recoPy, recoPx);

        if(recoPhi < -TMath::Pi()){
            recoPhi += 2.0*TMath::Pi();
        }else if(recoPhi > TMath::Pi()){
            recoPhi -= 2.0*TMath::Pi();
        }

        distPoint[0] = recoPt;
        distPoint[1] = recoPhi;
        distPoint[2] = recoEta;
        distPoint[3] = Zvertex;
        distPoint[4] = recoM;
        distPoint[5] = multPercentile;
        
        int posTrackID = ptrack->GetID();
        int negTrackID = ntrack->GetID();

        bool alreadyFilled = false;
        for(int i = 0; i < trackIDsV0.size(); i++) {
            if(posTrackID == trackIDsV0[i][0] && negTrackID == trackIDsV0[i][1]) alreadyFilled = true;
        }

        if(alreadyFilled) {
            std::cout << "There are duplicate V0s (whose daughters have exactly the same IDs: " << std::endl;
            std::cout << "\tPositive track ID: " << posTrackID << "; Negative track ID: " << negTrackID << std::endl;
            continue;
        }
        else {
            std::vector<int> v0tracks;
            v0tracks.push_back(posTrackID);
            v0tracks.push_back(negTrackID);
            trackIDsV0.push_back(v0tracks);
        }
        

        fRecoTotalV0LambdaDist->Fill(distPoint);

        uint posPassCuts = PassDaughterCuts(ptrack);
        uint negPassCuts = PassDaughterCuts(ntrack);

        if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
            fRecoEtaV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
            fRecoEtaPtV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
            fRecoEtaPtRefitV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
            fRecoEtaPtRefitRowsV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
            fRecoEtaPtRefitRowsRatioV0LambdaDist->Fill(distPoint);
        }

    }
    
    //first loop over tracks to get pi minuses, find lambda daughters and fill Reco Dist
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();

        if(tracklabel < 0) {
            continue;
        }

        //Get all single particle distributions
        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        pdgcode = mcpart->GetPdgCode();
        Bool_t associatedPass = PassAssociatedCuts(aodnegtrack);
        Bool_t triggerPass = PassTriggerCuts(aodnegtrack);
        singledistPoint[0] = aodnegtrack->Pt();
        singledistPoint[1] = aodnegtrack->Phi();
        if(singledistPoint[1] > TMath::Pi()){
            singledistPoint[1] -= 2.0*TMath::Pi();
        }else if(singledistPoint[1] < -1.0*TMath::Pi()){
            singledistPoint[1] += 2.0*TMath::Pi();
        }
        singledistPoint[2] = aodnegtrack->Eta();
        singledistPoint[3] = Zvertex;
        singledistPoint[4] = multPercentile;

        if(mcpart->IsPhysicalPrimary()){
            if(associatedPass){ 
                if(TMath::Abs(pdgcode) == 211){
                    fRecoPiDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 321){
                    fRecoKDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecopDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoeDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoMuonDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }
            }
            if(triggerPass){
                if(TMath::Abs(pdgcode) == 321){
                    fRecoKTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 211){
                    fRecoPiTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecopTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                    fRecoeTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoMuonTriggerDist->Fill(singledistPoint);
                    fRecoChargedTriggerDist->Fill(singledistPoint);
                }

            }

        }

        negPassCuts = PassDaughterCuts(aodnegtrack);

        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        negparPDG = mcnegpart->GetPdgCode();

        if(negparPDG != -211) continue;

        motherIndex = mcnegpart->GetMother();
        if(motherIndex < 0) continue;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
        motherPDG = mcmother->GetPdgCode();
        if(motherPDG != 3122) continue;

        for(int jtrack = 0; jtrack < ntracks; jtrack++){
            AliVParticle *vpospart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(jtrack));
            AliVTrack *postrack = dynamic_cast<AliVTrack*>(vpospart);
            AliAODTrack *aodpostrack = dynamic_cast<AliAODTrack*>(vpospart);

            Double_t posTPCnSigma = -999;
            Double_t posTOFnSigma = -999;
            if(postrack->Pt() > 0.15){
                posTPCnSigma = fpidResponse->NumberOfSigmasTPC(postrack, AliPID::kProton);
                posTOFnSigma = fpidResponse->NumberOfSigmasTOF(postrack, AliPID::kProton);
            }

            posPassCuts = PassDaughterCuts(aodpostrack);

            Int_t postracklabel = aodpostrack->GetLabel();
            if(postracklabel < 0) continue;

            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(postracklabel);
            posparPDG = mcpospart->GetPdgCode();
            if(posparPDG != 2212) continue;

            if(mcpospart->GetMother() == motherIndex){

                recoPx = aodnegtrack->Px() + aodpostrack->Px();
                recoPy = aodnegtrack->Py() + aodpostrack->Py();
                recoPz = aodnegtrack->Pz() + aodpostrack->Pz();
                
                recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
                recoE = TMath::Sqrt(aodnegtrack->Px()*aodnegtrack->Px() + aodnegtrack->Py()*aodnegtrack->Py() + aodnegtrack->Pz()*aodnegtrack->Pz() + 0.13957*0.13957) + TMath::Sqrt(aodpostrack->Px()*aodpostrack->Px() + aodpostrack->Py()*aodpostrack->Py() + aodpostrack->Pz()*aodpostrack->Pz() + 0.9383*0.9383);
                recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
                recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
                recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
                recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
                recoPhi = TMath::ATan2(recoPy, recoPx);
                if(recoPhi < -TMath::Pi()){
                    recoPhi += 2.0*TMath::Pi();
                }else if(recoPhi > TMath::Pi()){
                    recoPhi -= 2.0*TMath::Pi();
                }

                distPoint[0] = recoP;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;

                fRecoLambdaDist->Fill(distPoint);
                fRecoTotalLambdaDist->Fill(distPoint);

                double filter_distPoint[8];
                filter_distPoint[0] = recoPt;
                filter_distPoint[1] = recoPhi;
                filter_distPoint[2] = recoEta;
                filter_distPoint[3] = Zvertex;
                filter_distPoint[4] = recoM;
                filter_distPoint[5] = multPercentile;
                filter_distPoint[6] = filterMap_map[aodpostrack->GetFilterMap()];
                filter_distPoint[7] = filterMap_map[aodnegtrack->GetFilterMap()];
                
                fRecoTotalLambdaFilterDist->Fill(filter_distPoint);

                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fRecoEtaLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fRecoEtaPtLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
                    fRecoEtaPtRefitLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
                    fRecoEtaPtRefitRowsLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
                    fRecoEtaPtRefitRowsRatioLambdaDist->Fill(distPoint);
                    numRecoLambdas += 1;
                }
            }
        }
    }

    //second loop over tracks to get antiprotons minuses, find lambda daughters and fill Reco Dist
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();
        if(tracklabel < 0) continue;

        //Get pions that came from lambda for lambda reco
        negPassCuts = PassDaughterCuts(aodnegtrack);

        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        negparPDG = mcnegpart->GetPdgCode();
        if(negparPDG != -2212) continue;

        motherIndex = mcnegpart->GetMother();
        if(motherIndex < 0) continue;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
        motherPDG = mcmother->GetPdgCode();
        if(motherPDG != -3122) continue;

        for(int jtrack = 0; jtrack < ntracks; jtrack++){
            AliVParticle *vpospart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(jtrack));
            AliVTrack *postrack = dynamic_cast<AliVTrack*>(vpospart);
            AliAODTrack *aodpostrack = dynamic_cast<AliAODTrack*>(vpospart);

            posPassCuts = PassDaughterCuts(aodpostrack);

            Int_t postracklabel = aodpostrack->GetLabel();
            if(postracklabel < 0) continue;

            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(postracklabel);
            posparPDG = mcpospart->GetPdgCode();
            if(posparPDG != 211) continue;

            if(mcpospart->GetMother() == motherIndex){

                recoPx = aodnegtrack->Px() + aodpostrack->Px();
                recoPy = aodnegtrack->Py() + aodpostrack->Py();
                recoPz = aodnegtrack->Pz() + aodpostrack->Pz();
                
                recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
                recoE = TMath::Sqrt(aodnegtrack->Px()*aodnegtrack->Px() + aodnegtrack->Py()*aodnegtrack->Py() + aodnegtrack->Pz()*aodnegtrack->Pz() + 0.9383*0.9383) + TMath::Sqrt(aodpostrack->Px()*aodpostrack->Px() + aodpostrack->Py()*aodpostrack->Py() + aodpostrack->Pz()*aodpostrack->Pz() + 0.13957*0.13957);
                recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
                recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
                recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
                recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
                recoPhi = TMath::ATan2(recoPy, recoPx);
                if(recoPhi < -TMath::Pi()){
                    recoPhi += 2.0*TMath::Pi();
                }else if(recoPhi > TMath::Pi()){
                    recoPhi -= 2.0*TMath::Pi();
                }
                distPoint[0] = recoP;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;

                fRecoAntiLambdaDist->Fill(distPoint);
                fRecoTotalLambdaDist->Fill(distPoint);

                double filter_distPoint[8];
                filter_distPoint[0] = recoPt;
                filter_distPoint[1] = recoPhi;
                filter_distPoint[2] = recoEta;
                filter_distPoint[3] = Zvertex;
                filter_distPoint[4] = recoM;
                filter_distPoint[5] = multPercentile;
                filter_distPoint[6] = filterMap_map[aodnegtrack->GetFilterMap()];
                filter_distPoint[7] = filterMap_map[aodpostrack->GetFilterMap()];

                fRecoTotalLambdaFilterDist->Fill(filter_distPoint);

                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fRecoEtaLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fRecoEtaPtLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
                    fRecoEtaPtRefitLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
                    fRecoEtaPtRefitRowsLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
                    fRecoEtaPtRefitRowsRatioLambdaDist->Fill(distPoint);
                    numRecoLambdas += 1;
                }
            }
        }
    }

    //loop over MC particles to get original lambda and fill Real dist (and real single particle dist)
    int numCharged = 0;
    for(Int_t imcpart=0; imcpart< fMCArray->GetEntries(); imcpart++){

        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imcpart);
        pdgcode = AODMCtrack->GetPdgCode();

        singledistPoint[0] = AODMCtrack->Pt();
        singledistPoint[1] = AODMCtrack->Phi();
        if(singledistPoint[1] > TMath::Pi()){
            singledistPoint[1] -= 2.0*TMath::Pi();
        }else if(singledistPoint[1] < -1.0*TMath::Pi()){
            singledistPoint[1] += 2.0*TMath::Pi();
        }
        singledistPoint[2] = AODMCtrack->Eta();
        singledistPoint[3] = Zvertex;
        singledistPoint[4] = multPercentile;

        if(AODMCtrack->IsPhysicalPrimary()){
            if(TMath::Abs(pdgcode) == 321){
                fRealKDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 211){
                fRealPiDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealpDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 11){
                fRealeDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 13){
                fRealMuonDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }

        }

        // Get protons and pions that came from lambdas (would not be physical primaries) 
        if(TMath::Abs(pdgcode) == 211){
            motherIndex = AODMCtrack->GetMother();
            if(motherIndex >= 0) {
                AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                motherPDG = mcmother->GetPdgCode();
                if(TMath::Abs(motherPDG) == 3122) fRealPiFromLambdaDist->Fill(singledistPoint);
            } 
        } else if(TMath::Abs(pdgcode) == 2212){
            if(motherIndex >= 0) {
                AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                motherPDG = mcmother->GetPdgCode();
                if(TMath::Abs(motherPDG) == 3122) fRealpFromLambdaDist->Fill(singledistPoint);
            } 
        }

        //select phis
        if(TMath::Abs(pdgcode) != 3122) continue;
        Int_t indexFirstDaughter = 0, indexSecondDaughter = 0;
        indexFirstDaughter = AODMCtrack->GetDaughterFirst();
        indexSecondDaughter = AODMCtrack->GetDaughterLast();

        if(indexFirstDaughter < 0 || indexSecondDaughter < 0) continue;
        AliAODMCParticle* firstDaughter = (AliAODMCParticle*)fMCArray->At(indexFirstDaughter);
        AliAODMCParticle* secondDaughter = (AliAODMCParticle*)fMCArray->At(indexSecondDaughter);

        if(((TMath::Abs(firstDaughter->GetPdgCode()) == 211 && TMath::Abs(secondDaughter->GetPdgCode()) == 2212) ||
            (TMath::Abs(firstDaughter->GetPdgCode()) == 2212 && TMath::Abs(secondDaughter->GetPdgCode()) == 211)) && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) < 0){


            numRealLambdas += 1;
            distPoint[0] = AODMCtrack->Pt();
            distPoint[1] = AODMCtrack->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            distPoint[2] = AODMCtrack->Eta();
            distPoint[3] = Zvertex;
            distPoint[4] = AODMCtrack->GetCalcMass();
            distPoint[5] = multPercentile;
            fRealTotalLambdaDist->Fill(distPoint);

            if(pdgcode == 3122) {
                fRealLambdaDist->Fill(distPoint); 
            }
            else if(pdgcode == -3122) {
                fRealAntiLambdaDist->Fill(distPoint);
            }
            else{
                std::cout << "OH FUCK OH SHIT OH NO, THIS ISN'T A LAMBDA, THE PDG CODE IS INCORRECT: " << pdgcode << std::endl;
            }
        } 
    }
    
    fRealLambdasPerEvent->Fill(numRealLambdas);
    fRecoLambdasPerEvent->Fill(numRecoLambdas);

    PostData(1, fOutputList);
}    
//________________________________________________________________________
void AliAnalysisTaskLambdaHadronEfficiency::Terminate(Option_t *) 
{
    // Draw result to the screen
    // Called once at the end of the query
    printf("terminating task... \n");
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
    
}

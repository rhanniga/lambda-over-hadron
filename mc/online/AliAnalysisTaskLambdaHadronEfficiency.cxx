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
    MULT_LOW = multLow;
    MULT_HIGH = multHigh;
    CENT_ESTIMATOR = "V0A";
    KAON_TRK_BIT = 16;
    ASSOC_TRK_BIT = 16;
    TRIG_TRK_BIT = 768;
    KAON_ETA_CUT = 0.8;

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
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    KAON_TRK_BIT = 16;
    ASSOC_TRK_BIT = 16;
    TRIG_TRK_BIT = 768;
    KAON_ETA_CUT = 0.8;
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

    fRecoVsRealLambdaPtDist = new TH2D("fRecoVsRealLambdaPtDist", "Reco vs. Real #Lambda p_{T} dist", 100, 0, 10, 100, 0, 10);
    fOutputList->Add(fRecoVsRealLambdaPtDist);
    
    // Single Particle and inclusive charged hadron histos
    Int_t numbinsSingle[5] = {95, 64, 64, 10, 10};
    Double_t minvalSingle[5] = {0.5, -3.14159, -2.0, -10.0, 0.0};
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
    
    fTOFPiDist = new THnSparseF("fTOFPiDist", "TOF Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fTOFPiDist);
    
    fTOFProtonDist = new THnSparseF("fTOFProtonDist", "TOF Proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fTOFProtonDist);

    fRecoPiDist = new THnSparseF("fRecoPiDist", "Reco Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiDist);

    fRecoPiFromLambdaDist = new THnSparseF("fRecoPiFromLambdaDist", "Reco Pion (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecoPiFromLambdaDist);

    fRecopFromLambdaDist = new THnSparseF("fRecopFromLambdaDist", "Reco proton (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fOutputList->Add(fRecopFromLambdaDist);

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

    //Real & MC Phi track histos
    
    Int_t numbins[6] = {95, 64, 64, 10, 80, 10};
    Double_t minval[6] = {0.5, -3.14159, -2., -10, 1.06, 0.0};
    Double_t maxval[6] = {10.0, 3.14159,  2.,  10, 1.16, 100.0};

    fRealLambdaDist = new THnSparseF("fRealLambdaDist", "Real #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealLambdaDist);

    fRealNoDecayCutLambdaDist = new THnSparseF("fRealNoDecayCutLambdaDist", "Real #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRealNoDecayCutLambdaDist);

    fRecoLambdaDist = new THnSparseF("fRecoLambdaDist", "Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoLambdaDist);

    fRecoNormalLambdaDist = new THnSparseF("fRecoNormalLambdaDist", "Reco #Lambda (normal) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoNormalLambdaDist);

    fRecoAntiLambdaDist = new THnSparseF("fRecoAntiLambdaDist", "Reco #Lambda (anti) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fRecoAntiLambdaDist);

    fTrackRecoLambdaDist = new THnSparseF("fTrackRecoLambdaDist", "Track Eta, Pt, FilterMask and Crossed Rows Cut Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackRecoLambdaDist);

    fTrackEtaRecoLambdaDist = new THnSparseF("fTrackEtaRecoLambdaDist", "Track Eta Cut (< 0.8) Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackEtaRecoLambdaDist);

    fTrackEtaPtRecoLambdaDist = new THnSparseF("fTrackEtaPtRecoLambdaDist", "Track Eta and Pt Cut (< 0.8, > 0.15 GeV) Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackEtaPtRecoLambdaDist);

    fTrackEtaPtFilterRecoLambdaDist = new THnSparseF("fTrackEtaPtFilterRecoLambdaDist", "Track Eta, Pt and Filtermask Cut (< 0.8, > 0.15 GeV, Filter mask = 16) Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackEtaPtFilterRecoLambdaDist);

    // THIS SHOULD BE THE SAME AS fTrackRecoLambdaDist
    fTrackEtaPtFilterRowsRecoLambdaDist = new THnSparseF("fTrackEtaPtFilterRowsRecoLambdaDist", "Track Eta, Pt, Filtermask and rows Cut (< 0.8, > 0.15 GeV, Filter mask = 16, rows > 80) Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTrackEtaPtFilterRowsRecoLambdaDist);

    fTOFRecoLambdaDist = new THnSparseF("fTOFRecoLambdaDist","Track Cut & TOF hit Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.",6, numbins, minval, maxval);
    fOutputList->Add(fTOFRecoLambdaDist);
    
    fTPCPIDTrackRecoLambdaDist = new THnSparseF("fTPCPIDTrackRecoLambdaDist", "Track Cut & TPC PID 3#sigma Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fTPCPIDTrackRecoLambdaDist);
    
    fTPCPIDRecoLambdaDist = new THnSparseF("fTPCPIDRecoLambdaDist", "Track Cut & TOF hit & TPC PID 3#sigma Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fTPCPIDRecoLambdaDist);

    fPIDRecoLambdaDist = new THnSparseF("fPIDRecoLambdaDist", "Track Cut & TOF hit & TOF&TPC 3#sigma Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fOutputList->Add(fPIDRecoLambdaDist);

    fReactionPlane = new TH1D("fReactionPlane", "Reaction Plane Angle; Angle", 64, -3.14159, 3.14159);
    fOutputList->Add(fReactionPlane);


    // pt Lambda, eta Lambda, pt Proton, pt Pion
    int numBinsDaughters[4] = {128, 128, 128, 128}; 
    double binsMinDaughters[4] = {0, -2, 0, 0};
    double binsMaxDaughters[4] = {10, 2, 10, 10};

    fRealLambdaDaughterDist = new THnSparseF("fRealLambdaDaughterDist", "Real Lambda Daughter (p and #pi) distributions", 4, numBinsDaughters, binsMinDaughters, binsMaxDaughters);
    fOutputList->Add(fRealLambdaDaughterDist);

    fRecoLambdaDaughterDist = new THnSparseF("fRecoLambdaDaughterDist", "Track Cut Reco #Lambda Daughter (p and #pi) distribution", 4, numBinsDaughters, binsMinDaughters, binsMaxDaughters);
    fOutputList->Add(fRecoLambdaDaughterDist);

    fRealLambdasPerEvent = new TH1F("fRealLambdasPerEvent", "Real Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRealLambdasPerEvent);

    fRecoLambdasPerEvent = new TH1F("fRecoLambdasPerEvent", "Reco Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRecoLambdasPerEvent);

    PostData(1,fOutputList);

}

//_______________________________________________________________________
UInt_t AliAnalysisTaskLambdaHadronEfficiency::PassProtonCuts(AliAODTrack *track, Double_t TPCnSigma, Double_t TOFnSigma){
    //returns the level of cuts that the track passed
    //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
    UInt_t passLevel = 0;
    Bool_t pass = kTRUE;
    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);
    pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
    pass = pass && (track->GetTPCCrossedRows() > 80);
    if(pass) passLevel |= TRACK_BIT;    
    if(TOFnSigma != -999){ //check if there is a TOF signal, but don't care what the signal is
        passLevel |= TOF_HIT_BIT;
    }
    if(TMath::Abs(TPCnSigma) <= 2.0){ // check that kaon passed the TPC nsigma cut
        passLevel |= TPC_PID_BIT;
    }
    if(TMath::Abs(TOFnSigma) <= 2.0){ // check that kaon passed TOF nsigma cut
        passLevel |= TOF_PID_BIT;
    }

    if(TMath::Abs(track->Eta()) <= 0.8) passLevel |= ETA_BIT;
    if(track->Pt() >= 0.15) passLevel |= PT_BIT;
    if(track->TestFilterMask(ASSOC_TRK_BIT)) passLevel |= MASK_BIT;
    if(track->GetTPCCrossedRows() > 80) passLevel |= ROWS_BIT;

    return passLevel;
}

//_______________________________________________________________________
UInt_t AliAnalysisTaskLambdaHadronEfficiency::PassPionCuts(AliAODTrack *track, Double_t TPCnSigma, Double_t TOFnSigma){
    //returns the level of cuts that the track passed
    //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
    UInt_t passLevel = 0;
    Bool_t pass = kTRUE;
    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);
    pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
    pass = pass && (track->GetTPCCrossedRows() > 80);
    if(pass) passLevel |= TRACK_BIT;    
    if(TOFnSigma != -999){ //check if there is a TOF signal, but don't care what the signal is
        passLevel |= TOF_HIT_BIT;
    }
    if(TMath::Abs(TPCnSigma) <= 3.0){ // check that kaon passed the TPC nsigma cut
        passLevel |= TPC_PID_BIT;
    }
    if(TMath::Abs(TOFnSigma) <= 3.0){ // check that kaon passed TOF nsigma cut
        passLevel |= TOF_PID_BIT;
    }

    if(TMath::Abs(track->Eta()) <= 0.8) passLevel |= ETA_BIT;
    if(track->Pt() >= 0.15) passLevel |= PT_BIT;
    if(track->TestFilterMask(ASSOC_TRK_BIT)) passLevel |= MASK_BIT;
    if(track->GetTPCCrossedRows() > 80) passLevel |= ROWS_BIT;

    return passLevel;
}

//_______________________________________________________________________
Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassHadronCuts(AliAODTrack *track, Bool_t isTrigger){
    //returns the level of cuts that the track passed
    //cutLevel: 1 = track cuts, 2 = TOF Hit, 4 = TPC PID cut, 8 = TOF PID cut
    Bool_t pass = kTRUE;
    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);
    if(isTrigger){
        pass = pass && (track->TestBit(TRIG_TRK_BIT));
    }else{
        pass = pass && (track->TestFilterMask(ASSOC_TRK_BIT));
    }
    return pass;
}

//________________________________________________________________________
void AliAnalysisTaskLambdaHadronEfficiency::UserExec(Option_t *){

    //masks for the different cut configurations: (track cuts only), (track + TOF hit), (track + TPC PID), (track + TOF hit + TPC PID), (track + TOF hit + TPC PID + TOF PID)
    UInt_t maskEta = ETA_BIT;
    UInt_t maskEtaPt = ETA_BIT + PT_BIT;
    UInt_t maskEtaPtMask = ETA_BIT + PT_BIT + MASK_BIT;
    UInt_t maskEtaPtMaskRows = ETA_BIT + PT_BIT + MASK_BIT + ROWS_BIT;
    UInt_t maskTrackOnly = TRACK_BIT;
    UInt_t maskTrackTOF = TRACK_BIT + TOF_HIT_BIT;
    UInt_t maskTrackTPC = TRACK_BIT + TPC_PID_BIT;
    UInt_t maskTrackTOFTPC = TRACK_BIT + TOF_HIT_BIT + TPC_PID_BIT;
    UInt_t maskTrackPID = TRACK_BIT + TOF_HIT_BIT + TPC_PID_BIT + TOF_PID_BIT;

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

    //TRandom3* randomGen = new TRandom3();
    //Double_t randangle = randomGen->Uniform(0.0, 2.0*TMath::Pi());

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
    
    //first loop over tracks to get pi minuses, find lambda daughters and fill Reco Dist
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();
        // std::cout << tracklabel << " is the track label" << std::endl;
        if(tracklabel < 0) {
            AliAODMCParticle* part = (AliAODMCParticle*)fMCArray->At(TMath::Abs(tracklabel));
            // std::cout << part->GetPdgCode() << " is the PDG code" << std::endl;
            // std::cout << part->GetPdgCode() << std::endl;
            // std::cout << "hello! the track label is: " << tracklabel << std::endl;
            continue;
        }
    
        Double_t negTPCnSigma = -999;
        Double_t negTOFnSigma = -999;
        if(negtrack->Pt() > 0.15){
            negTPCnSigma = fpidResponse->NumberOfSigmasTPC(negtrack, AliPID::kPion);
            negTOFnSigma = fpidResponse->NumberOfSigmasTOF(negtrack, AliPID::kPion);
        }

        //Get all single particle distributions
        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        pdgcode = mcpart->GetPdgCode();
        Bool_t pass = PassHadronCuts(aodnegtrack, kFALSE);
        Bool_t trigpass = kFALSE;
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
            if(pass){ 
                if(TMath::Abs(pdgcode) == 211){
                    Int_t pionPass = PassPionCuts(aodnegtrack, negTPCnSigma, negTOFnSigma);
                    fRecoPiDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                    if((pionPass & maskTrackTOF) == maskTrackTOF){
                        fTOFPiDist->Fill(singledistPoint);
                    }
                }else if(TMath::Abs(pdgcode) == 321){
                    fRecoKDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecopDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                    motherIndex = mcpart->GetMother();
                    if(motherIndex >= 0) {
                        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                        motherPDG = mcmother->GetPdgCode();
                        if(TMath::Abs(motherPDG) == 3122) fRecopFromLambdaDist->Fill(singledistPoint);
                    } 
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoeDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoMuonDist->Fill(singledistPoint);
                    fRecoChargedDist->Fill(singledistPoint);
                }
            }
            if(trigpass){
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

        // Get reco pions and protons that came from lambdas (can't be physical primaries)
        if(pass){ 
            if(TMath::Abs(pdgcode) == 211){
                motherIndex = mcpart->GetMother();
                if(motherIndex >= 0) {
                    AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                    motherPDG = mcmother->GetPdgCode();
                    if(TMath::Abs(motherPDG) == 3122) fRecoPiFromLambdaDist->Fill(singledistPoint);
                } 
            } else if(TMath::Abs(pdgcode) == 2212){
                motherIndex = mcpart->GetMother();
                if(motherIndex >= 0) {
                    AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                    motherPDG = mcmother->GetPdgCode();
                    if(TMath::Abs(motherPDG) == 3122) fRecopFromLambdaDist->Fill(singledistPoint);
                } 
            }
        }

        negPassCuts = PassPionCuts(aodnegtrack,  negTPCnSigma, negTOFnSigma);
        // if(negPassCuts == 0) continue;

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

            posPassCuts = PassProtonCuts(aodpostrack, posTPCnSigma, posTOFnSigma);
            // if(posPassCuts == 0) continue;

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
                distPoint[0] = recoPt;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;
                
                fRecoLambdaDist->Fill(distPoint);
                fRecoNormalLambdaDist->Fill(distPoint);

                fRecoVsRealLambdaPtDist->Fill(recoPt, mcmother->Pt());

                //fill with phi's where daughter kaons pass track cuts
                if(((negPassCuts & maskTrackOnly) == maskTrackOnly) && ((posPassCuts & maskTrackOnly)== maskTrackOnly)){
                    fTrackRecoLambdaDist->Fill(distPoint);
                    numRecoLambdas += 1;
                    if(recoM <= 1.12 && recoM >= 1.11) {
                        double daughterDistPoint[4];
                        daughterDistPoint[0] = recoPt;
                        daughterDistPoint[1] = recoEta;
                        // THE PION IS THE NEGATIVE TRACK, PROTON IS POSITIVE
                        daughterDistPoint[2] = aodpostrack->Pt();
                        daughterDistPoint[3] = aodnegtrack->Pt();
                        fRecoLambdaDaughterDist->Fill(daughterDistPoint);
                    }
                }

                //fill with phi's where daughter kaons pass track cuts and have a TOF hit
                if(((negPassCuts & maskTrackTOF) == maskTrackTOF) && ((posPassCuts & maskTrackTOF)== maskTrackTOF)){
                    fTOFRecoLambdaDist->Fill(distPoint);
                }
               
                //fill with phi's where daughter kaons pass track cuts and TPC PID cuts
                if(((negPassCuts & maskTrackTPC) == maskTrackTPC) && ((posPassCuts & maskTrackTPC)== maskTrackTPC)){
                    fTPCPIDTrackRecoLambdaDist->Fill(distPoint);
                } 
                
                //fill with phi's where daughter kaons pass track cuts and TOF hit and TPC PID cuts 
                if(((negPassCuts & maskTrackTOFTPC) == maskTrackTOFTPC) && ((posPassCuts & maskTrackTOFTPC)== maskTrackTOFTPC)){
                    fTPCPIDRecoLambdaDist->Fill(distPoint);
                }

                //fill with phi daughter kaons that pass track cuts and pass TOF&TPC PID cuts
                if(((negPassCuts & maskTrackPID) == maskTrackPID) && ((posPassCuts & maskTrackPID)== maskTrackPID)){
                    fPIDRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta cuts
                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fTrackEtaRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts
                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fTrackEtaPtRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts and filtermask cuts
                if(((negPassCuts & maskEtaPtMask) == maskEtaPtMask) && ((posPassCuts & maskEtaPtMask)== maskEtaPtMask)){
                    fTrackEtaPtFilterRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts and filtermask cuts and crossed rows
                // NOTE: This should completely match the original track cut distribution
                if(((negPassCuts & maskEtaPtMaskRows) == maskEtaPtMaskRows) && ((posPassCuts & maskEtaPtMaskRows)== maskEtaPtMaskRows)){
                    fTrackEtaPtFilterRowsRecoLambdaDist->Fill(distPoint);
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
    
        Double_t negTPCnSigma = -999;
        Double_t negTOFnSigma = -999;
        if(negtrack->Pt() > 0.15){
            negTPCnSigma = fpidResponse->NumberOfSigmasTPC(negtrack, AliPID::kProton);
            negTOFnSigma = fpidResponse->NumberOfSigmasTOF(negtrack, AliPID::kProton);
        }

        //Get all single particle distributions
        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        if(mcpart->IsPhysicalPrimary()){
            pdgcode = mcpart->GetPdgCode();
            Bool_t pass = PassHadronCuts(aodnegtrack, kFALSE);
            Bool_t trigpass = kFALSE;
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
            if(pass){ 
                if(TMath::Abs(pdgcode) == 2212){
                    Int_t protonPass = PassProtonCuts(aodnegtrack, negTPCnSigma, negTOFnSigma);
                    if((protonPass & maskTrackTOF) == maskTrackTOF){
                        fTOFProtonDist->Fill(singledistPoint);
                    }
                }
            }
        }

        //Get pions that came from lambda for lambda reco
        negPassCuts = PassProtonCuts(aodnegtrack,  negTPCnSigma, negTOFnSigma);
        // if(negPassCuts == 0) continue;

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

            Double_t posTPCnSigma = -999;
            Double_t posTOFnSigma = -999;
            if(postrack->Pt() > 0.15){
                posTPCnSigma = fpidResponse->NumberOfSigmasTPC(postrack, AliPID::kPion);
                posTOFnSigma = fpidResponse->NumberOfSigmasTOF(postrack, AliPID::kPion);
            }

            posPassCuts = PassPionCuts(aodpostrack, posTPCnSigma, posTOFnSigma);
            // if(posPassCuts == 0) continue;

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
                distPoint[0] = recoPt;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;
                
                fRecoLambdaDist->Fill(distPoint);
                fRecoAntiLambdaDist->Fill(distPoint);

                fRecoVsRealLambdaPtDist->Fill(recoPt, mcmother->Pt());

                //fill with phi's where daughter kaons pass track cuts
                if(((negPassCuts & maskTrackOnly) == maskTrackOnly) && ((posPassCuts & maskTrackOnly)== maskTrackOnly)){
                    fTrackRecoLambdaDist->Fill(distPoint);
                    numRecoLambdas += 1;
                    if(recoM <= 1.12 && recoM >= 1.11) {
                        double daughterDistPoint[4];
                        daughterDistPoint[0] = recoPt;
                        daughterDistPoint[1] = recoEta;
                        // THE PION IS THE POSITIVE TRACK, ANTIPROTON IS NEGATIVE
                        daughterDistPoint[2] = aodnegtrack->Pt();
                        daughterDistPoint[3] = aodpostrack->Pt();
                        fRecoLambdaDaughterDist->Fill(daughterDistPoint);
                    }
                }
                //fill with phi's where daughter kaons pass track cuts and have a TOF hit
                if(((negPassCuts & maskTrackTOF) == maskTrackTOF) && ((posPassCuts & maskTrackTOF)== maskTrackTOF)){
                    fTOFRecoLambdaDist->Fill(distPoint);
                }
               
                //fill with phi's where daughter kaons pass track cuts and TPC PID cuts
                if(((negPassCuts & maskTrackTPC) == maskTrackTPC) && ((posPassCuts & maskTrackTPC)== maskTrackTPC)){
                    fTPCPIDTrackRecoLambdaDist->Fill(distPoint);
                } 
                
                //fill with phi's where daughter kaons pass track cuts and TOF hit and TPC PID cuts 
                if(((negPassCuts & maskTrackTOFTPC) == maskTrackTOFTPC) && ((posPassCuts & maskTrackTOFTPC)== maskTrackTOFTPC)){
                    fTPCPIDRecoLambdaDist->Fill(distPoint);
                }

                //fill with phi daughter kaons that pass track cuts and pass TOF&TPC PID cuts
                if(((negPassCuts & maskTrackPID) == maskTrackPID) && ((posPassCuts & maskTrackPID)== maskTrackPID)){
                    fPIDRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta cuts
                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fTrackEtaRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts
                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fTrackEtaPtRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts and filtermask cuts
                if(((negPassCuts & maskEtaPtMask) == maskEtaPtMask) && ((posPassCuts & maskEtaPtMask)== maskEtaPtMask)){
                    fTrackEtaPtFilterRecoLambdaDist->Fill(distPoint);
                }

                // Fill with lambdas whose daughters pass eta and pt cuts and filtermask cuts and crossed rows
                // NOTE: This should completely match the original track cut distribution
                if(((negPassCuts & maskEtaPtMaskRows) == maskEtaPtMaskRows) && ((posPassCuts & maskEtaPtMaskRows)== maskEtaPtMaskRows)){
                    fTrackEtaPtFilterRowsRecoLambdaDist->Fill(distPoint);
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

        //select only lambda that decay to p-pi
        if(((TMath::Abs(firstDaughter->GetPdgCode()) == 211 && TMath::Abs(secondDaughter->GetPdgCode()) == 2212) ||
            (TMath::Abs(firstDaughter->GetPdgCode()) == 2212 && TMath::Abs(secondDaughter->GetPdgCode()) == 211)) && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) < 0 && TMath::Abs(firstDaughter->Eta()) <= 0.8 && TMath::Abs(secondDaughter->Eta()) <= 0.8){

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
            fRealLambdaDist->Fill(distPoint);
        } 

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
            fRealNoDecayCutLambdaDist->Fill(distPoint);


            double daughterDistPoint[4];
            daughterDistPoint[0] = AODMCtrack->Pt();
            daughterDistPoint[1] = AODMCtrack->Eta();
            if(TMath::Abs(firstDaughter->GetPdgCode()) == 211) {
                daughterDistPoint[2] = secondDaughter->Pt();
                daughterDistPoint[3] = firstDaughter->Pt();
            }
            else {
                daughterDistPoint[2] = firstDaughter->Pt();
                daughterDistPoint[3] = secondDaughter->Pt();
            }
            fRealLambdaDaughterDist->Fill(daughterDistPoint);
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
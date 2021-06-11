/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronEfficiency class
                This task is for determining the lambda, trigger, and associated 
                particle efficiency for the AliAnalysisTaskLambdaHadronRatio task
                (can be found in same directory)
                Origin: Ryan Hannigan, January 2021, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronEfficiency_cxx
#define AliAnalysisTaskLambdaHadronEfficiency_cxx

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include <map>

#include "AliAnalysisTaskSE.h"
#include "TString.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TH1F;
class THnSparseF;

class AliEventPoolManager;
class AliESDEvent;
class AliVEvent;
class AliAODEvent;
class AliAODTrack;
class AliAODMCHeader;
class AliMultSelection;
class AliPIDResponse;


class AliAnalysisTaskLambdaHadronEfficiency : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskLambdaHadronEfficiency();
    AliAnalysisTaskLambdaHadronEfficiency(const char *name, Float_t multLow, Float_t multHigh);
    virtual ~AliAnalysisTaskLambdaHadronEfficiency();
    
    virtual void   UserCreateOutputObjects();
    uint PassDaughterCuts(AliAODTrack* track);
    Bool_t PassAssociatedCuts(AliAODTrack* track);
    Bool_t PassTriggerCuts(AliAODTrack* track);
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    
    void SetCentEstimator(TString est){ CENT_ESTIMATOR = est; };
    void SetTrigTrkBit(UInt_t trkbit) { TRIG_TRK_BIT = trkbit; };
    void SetAssocTrkBit(UInt_t trkbit) { ASSOC_TRK_BIT = trkbit; };

    AliAODTrack* GetTrackFromID(AliAODEvent* inputEvent, int trackID);

    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };

private:

    int MIN_CROSSED_ROWS_TPC;
    float MIN_ROW_CLUSTER_RATIO;

    Float_t MULT_LOW; // lower bound of mult
    Float_t MULT_HIGH; // upper bound of mult
    Float_t DAUGHTER_ETA_CUT; // eta cut on daughters
    Float_t DAUGHTER_MIN_PT; // min pt of daughters
    Float_t ASSOC_TRK_BIT; // filter bit for associated hadrons
    Float_t TRIG_TRK_BIT; // "filter" bit for trigger hadrons
    TString CENT_ESTIMATOR; // method used for cent estimator (default V0A)

    // Bit maps for differential efficiency investigation
    uint ETA_BIT = 1 << 0;
    uint PT_BIT = 1 << 1;
    uint TPC_REFIT_BIT = 1 << 2;
    uint CROSSED_ROWS_BIT = 1 << 3;
    uint ROW_CLUSTER_RATIO_BIT = 1 << 4;

    //
    std::map<int, int> filterMap_map = {{0, 0},
                                {1, 1},
                                {2, 2},
                                {5, 3},
                                {21, 4},
                                {128, 5},
                                {512, 6},
                                {3077, 7},
                                {3381, 8}};



    enum{
        kAODanalysis = BIT(20),
    };
      
    AliVEvent   *fVevent;  //!event object
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    AliMultSelection *fMultSelection; //!mult selection
   
    TClonesArray* fMCArray; //!
    AliAODMCHeader* fMCHeader; //!
    TList       *fOutputList; //!Output list

    TH1F        *fNevents;//! no of events
    TH1F        *fNumTracks;//! number of Tracks/evt
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fVtxX;//!Vertex x
    TH1F        *fVtxY;//!Vertex y

    THnSparseF  *fRealTotalLambdaDist;//! Dist of Real lambda and anti lambda
    THnSparseF  *fRealLambdaDist;//! Dist of Real lambda
    THnSparseF  *fRealAntiLambdaDist;//! Dist of Real anti lambda
    
    THnSparseF  *fRecoTotalV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0)
    THnSparseF  *fRecoEtaV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta cut on daughters)
    THnSparseF  *fRecoEtaPtV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt cut on daughters)
    THnSparseF  *fRecoEtaPtRefitV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit cut on daughters)
    THnSparseF  *fRecoEtaPtRefitRowsV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit rows  cut on daughters)
    THnSparseF  *fRecoEtaPtRefitRowsRatioV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit rows ratio cut on daughters)

    THnSparseF  *fRecoTotalLambdaDist;//! Dist of Recon lambda and anti lambda
    THnSparseF  *fRecoEtaLambdaDist;//! Dist of Recon lambda and anti lambda (eta cut on daughters)
    THnSparseF  *fRecoEtaPtLambdaDist;//! Dist of Recon lambda and anti lambda (eta pt cut on daughters)
    THnSparseF  *fRecoEtaPtRefitLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit cut on daughters)
    THnSparseF  *fRecoEtaPtRefitRowsLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit rows  cut on daughters)
    THnSparseF  *fRecoEtaPtRefitRowsRatioLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit rows ratio cut on daughters)

    THnSparseF  *fRecoTotalLambdaFilterDist;//! Dist of Recon lambda and anti lambda with filter bit of daughters
    THnSparseF  *fRecoLambdaDist;//! Dist of Recon lambda
    THnSparseF  *fRecoAntiLambdaDist;//! Dist of Recon anti lambda

    THnSparseF  *fRealChargedDist;//! real charged hadron dist
    THnSparseF  *fRealKDist;//! real charged K dist
    THnSparseF  *fRealPiDist;//! real charged pion dist
    THnSparseF  *fRealPiFromLambdaDist;//! real pions (from lambda) dist
    THnSparseF  *fRealeDist;//! real electron dist
    THnSparseF  *fRealpDist;//! real proton dist
    THnSparseF  *fRealpFromLambdaDist;//! real proton (from lambda) dist
    THnSparseF  *fRealMuonDist;//! real muon dist
   
    THnSparseF  *fRecoChargedDist;//!
    THnSparseF  *fRecoKDist;//!
    THnSparseF  *fRecoPiDist;//!
    THnSparseF  *fRecoeDist;//!
    THnSparseF  *fRecopDist;//!
    THnSparseF  *fRecoMuonDist;//!

    THnSparseF  *fRecoChargedTriggerDist;//!
    THnSparseF  *fRecoKTriggerDist;//!
    THnSparseF  *fRecoPiTriggerDist;//!
    THnSparseF  *fRecoeTriggerDist;//!
    THnSparseF  *fRecopTriggerDist;//!
    THnSparseF  *fRecoMuonTriggerDist;//!

    THnSparseF  *fRealLambdaDaughterDist;//!
    THnSparseF  *fRecoLambdaDaughterDist;//!

    TH1F        *fRealLambdasPerEvent;//!
    TH1F        *fRecoLambdasPerEvent;//!

    TH1D        *fReactionPlane;//!


    ClassDef(AliAnalysisTaskLambdaHadronEfficiency, 1); 
};

#endif
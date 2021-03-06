/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronEfficiency class
                This task is for determining the lambda, trigger, and associated 
                particle efficiency for the AliAnalysisTaskLambdaHadronRatio task
                (can be found in same directory)
                Origin: Ryan Hannigan, January 2021, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronEfficiency_H
#define AliAnalysisTaskLambdaHadronEfficiency_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include <map>

#include "TString.h"

#include "AliAnalysisTaskSE.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TH1F;
class THnSparse;
class TList;

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

    THnSparse  *fRealTotalLambdaDist;//! Dist of Real lambda and anti lambda
    THnSparse  *fRealLambdaDist;//! Dist of Real lambda
    THnSparse  *fRealAntiLambdaDist;//! Dist of Real anti lambda
    
    THnSparse  *fRecoTotalV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0)
    THnSparse  *fRecoEtaV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta cut on daughters)
    THnSparse  *fRecoEtaPtV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt cut on daughters)
    THnSparse  *fRecoEtaPtRefitV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit cut on daughters)
    THnSparse  *fRecoEtaPtRefitRowsV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit rows  cut on daughters)
    THnSparse  *fRecoEtaPtRefitRowsRatioV0LambdaDist;//! Dist of Recon lambda and anti lambda (from v0, eta pt refit rows ratio cut on daughters)

    THnSparse  *fRecoTotalLambdaDist;//! Dist of Recon lambda and anti lambda
    THnSparse  *fRecoEtaLambdaDist;//! Dist of Recon lambda and anti lambda (eta cut on daughters)
    THnSparse  *fRecoEtaPtLambdaDist;//! Dist of Recon lambda and anti lambda (eta pt cut on daughters)
    THnSparse  *fRecoEtaPtRefitLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit cut on daughters)
    THnSparse  *fRecoEtaPtRefitRowsLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit rows  cut on daughters)
    THnSparse  *fRecoEtaPtRefitRowsRatioLambdaDist;//! Dist of Recon lambda and anti lambda ( eta pt refit rows ratio cut on daughters)

    THnSparse  *fRecoTotalLambdaFilterDist;//! Dist of Recon lambda and anti lambda with filter bit of daughters
    THnSparse  *fRecoLambdaDist;//! Dist of Recon lambda
    THnSparse  *fRecoAntiLambdaDist;//! Dist of Recon anti lambda

    THnSparse  *fRealChargedDist;//! real charged hadron dist
    THnSparse  *fRealKDist;//! real charged K dist
    THnSparse  *fRealPiDist;//! real charged pion dist
    THnSparse  *fRealPiFromLambdaDist;//! real pions (from lambda) dist
    THnSparse  *fRealeDist;//! real electron dist
    THnSparse  *fRealpDist;//! real proton dist
    THnSparse  *fRealpFromLambdaDist;//! real proton (from lambda) dist
    THnSparse  *fRealMuonDist;//! real muon dist
  
    THnSparse  *fRecoChargedDist;//!
    THnSparse  *fRecoKDist;//!
    THnSparse  *fRecoPiDist;//!
    THnSparse  *fRecoeDist;//!
    THnSparse  *fRecopDist;//!
    THnSparse  *fRecoMuonDist;//!

    THnSparse  *fRecoChargedTriggerDist;//!
    THnSparse  *fRecoKTriggerDist;//!
    THnSparse  *fRecoPiTriggerDist;//!
    THnSparse  *fRecoeTriggerDist;//!
    THnSparse  *fRecopTriggerDist;//!
    THnSparse  *fRecoMuonTriggerDist;//!

    THnSparse  *fRealLambdaDaughterDist;//!
    THnSparse  *fRecoLambdaDaughterDist;//!

    TH1F        *fRealLambdasPerEvent;//!
    TH1F        *fRecoLambdasPerEvent;//!

    TH1D        *fReactionPlane;//!


    ClassDef(AliAnalysisTaskLambdaHadronEfficiency, 1); 
};

#endif
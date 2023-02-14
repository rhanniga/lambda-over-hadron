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
    unsigned int PassDaughterCuts(AliAODTrack* track);
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
    unsigned int ETA_BIT = 1 << 0;
    unsigned int PT_BIT = 1 << 1;
    unsigned int TPC_REFIT_BIT = 1 << 2;
    unsigned int CROSSED_ROWS_BIT = 1 << 3;
    unsigned int ROW_CLUSTER_RATIO_BIT = 1 << 4;

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
    THnSparse  *fRecoEtaPtRefitRowsRatioLambdaFilterDist;//! Dist of Recon lambda and anti lambda with filter bit of daughters (proton then pion)
    THnSparse  *fRecoEtaPtRefitRowsRatioLambdaDCADist;//! Dist of Recon lambda and anti lambda with DCA of daughters (proton then pion)
    
    THnSparse  *fRecoLambdaDist;//! Dist of Recon lambda
    THnSparse  *fRecoAntiLambdaDist;//! Dist of Recon anti lambda

    THnSparse  *fRecoLambdaWithAODPionDist;//! Dist of Recon lambda with cor AOD pion
    THnSparse  *fRecoLambdaWithAODProtonDist;//! Dist of Recon lambda with cor AOD proton

    THnSparse  *fRealChargedDist;//! real charged hadron dist
    THnSparse  *fRealPrimaryChargedDist;//! real primary charged hadron dist
    THnSparse  *fRealSecondaryChargedDist;//! real secondary charged hadron dist 
    THnSparse  *fRealTriggerDist;//! real trigger hadron dist (charged hadron, not always physical primary)
    THnSparse  *fRealKDist;//! real charged K dist
    THnSparse  *fRealPiDist;//! real charged pion dist
    THnSparse  *fRealPiFromLambdaDist;//! real pions (from lambda) dist
    THnSparse  *fRealeDist;//! real electron dist
    THnSparse  *fRealpDist;//! real proton dist
    THnSparse  *fRealpFromLambdaDist;//! real proton (from lambda) dist
    THnSparse  *fRealMuonDist;//! real muon dist
  
    THnSparse  *fRecoChargedDist;//!
    THnSparse  *fRecoPrimaryChargedDist;//!
    THnSparse  *fRecoSecondaryChargedDist;//!
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

    TH1D        *fPxDifference;//! distribution of difference between real px and track px
    TH1D        *fPxDifferenceFB;//! distribution of difference between real px and track px (tracks have ktrkglobalnodca)
    TH1D        *fPyDifference;//! distribution of difference between real py and track py
    TH1D        *fPzDifference;//! distribution of difference between real pz and track pz
    TH1D        *fPtDifference;//! distribution of difference between real pt (calculated from mc p) and track pt

    // Invariant mass comparison histograms

    TH1D       *fInvMassLambdaResonance; //!
    TH1D       *fInvMassAntiLambdaResonance; //!
    TH1D       *fInvMassLambdaV0; //!
    TH1D       *fInvMassAntiLambdaV0; //!
    TH1D       *fInvMassLambdaDifference; //!
    TH1D       *fInvMassAntiLambdaDifference; //!
    TH1D       *fInvMassLambdaReal; //!
    TH1D       *fInvMassAntiLambdaReal; //!

    TH1D       *fRealPrimaryLambdaPtDist; //!
    TH1D       *fRealSecondaryLambdaPtDist; //!

    // Lambda phi comparisons

    TH1D      *fPhiDifferenceResV0; //!
    TH1D      *fPhiDifferenceResReal; //!
    TH1D      *fPhiDifferenceV0Real; //!
    TH1D      *fPhiV0; //!
    TH1D      *fPhiRes; //!
    TH1D      *fPhiReal; //!

    ClassDef(AliAnalysisTaskLambdaHadronEfficiency, 1); 
};

#endif
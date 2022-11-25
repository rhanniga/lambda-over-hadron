/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronRatio class
                This task is for determining the lambda/hadron (pion) ratio
                in different kinematic regions using a two-particle azimuthal
                correlation method
                Origin: Ryan Hannigan, January 2021, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronRatio_H
#define AliAnalysisTaskLambdaHadronRatio_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH1F;
class THnSparse;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;


class AliAnalysisTaskLambdaHadronRatio : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronRatio();
  AliAnalysisTaskLambdaHadronRatio(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetAssociatedBit(float assocBit);
  void SetCentEstimator(TString estimator);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

private:
  float fMultLow; // lower bound for multiplicity
  float fMultHigh; // upper bound for multiplicity
  float fDaughterBit; // filter bit for daughter particle
  float fAssociatedBit; // filter bit for associated particle
  float fTriggerBit; // filter bit for trigger particle

  TString fCentEstimator;

  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  AliEventPoolManager *fCorPoolMgr_highestPt; //!>! correlation pool manager for highest pt trigger
  
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
  TH1D* fLambdaEff; ///> lambda efficiency

  TH2D* fTriggersAndLambdasPerEvent_All; //!>! triggers and all lambdas per event
  TH2D* fTriggersAndLambdasPerEvent_2_4; //!>! triggers and 2-4 GeV lambdas per event

  TH2D* fTofTest; //!>! tof test

  TH2D* fAssociatedPtEventClass; //!>! Associated hadron pt distribution for different event classes
  TH2D* fTriggerPtEventClass; //!>! Associated hadron pt distribution for different event classes
  TH2D* fLambdaPtEventClass; //!>! Associated hadron pt distribution for different event classes

  TH1D* fMultDistMinBias; //!>! Mult dist (min bias)
  TH1D* fMultDistHHEvent; //!>! Mult dist (events that we could do HH correlation)
  TH1D* fMultDistHLambdaEvent; //!>! Mult dist (events that we could do HH correlation)

  THnSparse* fLooseDist;  //!>! single particle all hadron dist (no cuts at all)
  THnSparse* fTriggerDist;  //!>! single particle trigger dist
  THnSparse* fTriggerDistEff;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDistEff_highestPt;  //!>! single particle trigger dist (corrected for efficiency, highest pt between 4 and 8)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparse* fLambdaDist;  //!>! single particle lambda dist
  THnSparse* fNormalLambdaDist;  //!>! single particle lambda dist (normal lambdas, not anti)
  THnSparse* fAntiLambdaDist;  //!>! single particle lambda dist (anti lambdas)
  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event
  THnSparse* fTriggeredLambdaDistFilterbit;  //!>! single particle lambda dist where daughters have filter bit 16 within a triggered event

  THnSparse* fDphiHLambda;  //!>! hadron-lambda correlation hist
  THnSparse* fDphiHLambdaFilterbit;  //!>! hadron-lambda correlation hist where daughter has filter bit 16
  THnSparse* fDphiHLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparse* fDphiHLambdaEff_highestPt;  //!>! hadron-lambda correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaV0;  //!>! hadron-lambda correlation hist (using v0 finder for lambda)
  THnSparse* fDphiHLambdaRotated;  //!>! hadron-lambda correlation hist with rotated pion
  THnSparse* fDphiHLambdaRotatedPi;  //!>! hadron-lambda correlation hist with daughter rotated by pi
  THnSparse* fDphiHLambdaRotatedProton;  //!>! hadron-lambda correlation hist with rotated proton
  THnSparse* fDphiHLambdaFlipped;  //!>! hadron-lambda correlation hist with flipped pion
  THnSparse* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparse* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHHEff_highestPt;   //!>! hadron-hadron correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiTriggerTrigger;   //!>! trigger-trigger correlation hist
  THnSparse* fDphiHLambdaLS; //!>! hadron-proton+pion like sign correlation hist
  THnSparse* fDphiHLambdaMixed_highestPt; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHHMixed_highestPt; //!>! hadron-hadron mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaLSMixed; //!>! hadron-proton+pion like sign mixed correlation hist
  THnSparse* fDphiTriggerTriggerMixed;   //!>! mixed trigger-trigger correlation hist

  THnSparse* fLambdaDaughterDCA; //!>!

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=false);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist);
  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskLambdaHadronRatio, 3);

};
#endif
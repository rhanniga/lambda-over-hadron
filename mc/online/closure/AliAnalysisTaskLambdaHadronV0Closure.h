/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronV0Closure class
                This is a closure task for the lambda-hadron correlation analysis
                (can be found in same directory)
                Origin: Ryan Hannigan, March 2022, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronV0Closure_H
#define AliAnalysisTaskLambdaHadronV0Closure_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declaration order: standard, ROOT, AliRoot (for pointers in this file only)
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
class AliAODv0;


class AliAnalysisTaskLambdaHadronV0Closure : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronV0Closure();
  AliAnalysisTaskLambdaHadronV0Closure(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronV0Closure();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetAssociatedBit(float assocBit);
  void SetCentEstimator(TString estimator);

  struct AliMotherContainer {
    AliAODv0* vzero;
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
  
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
  TH1D* fLambdaEff; ///> lambda efficiency





  THnSparse* fTriggerDistEff;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event

  THnSparse* fDphiHLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparse* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=false);
  void FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda);
  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);
  bool PassMCTriggerCuts(AliAODMCParticle *mc_track);
  bool PassMCAssociatedCuts(AliAODMCParticle *mc_track);
  bool PassMCLambdaCuts(AliAODMCParticle *mc_track);
  bool IsMCChargedHadron(int pdg_code);

  ClassDef(AliAnalysisTaskLambdaHadronV0Closure, 3);

};
#endif
/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronCloseTrackEff class
                This is a closure task for the lambda-hadron correlation analysis
                (can be found in same directory)
                Origin: Ryan Hannigan, March 2022, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronCloseTrackEff_H
#define AliAnalysisTaskLambdaHadronCloseTrackEff_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declaration order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH1F;
class THnSparse;
class TClonesArray;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;
class AliAODMCParticle;
class AliAODv0;


class AliAnalysisTaskLambdaHadronCloseTrackEff : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronCloseTrackEff();
  AliAnalysisTaskLambdaHadronCloseTrackEff(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronCloseTrackEff();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetCentEstimator(TString estimator);

private:
  float fMultLow; // lower bound for multiplicity
  float fMultHigh; // upper bound for multiplicity
  float fTriggerBit; // filter bit for trigger particle

  TString fCentEstimator;

  AliAODEvent* fAOD; //!>! input event
  TClonesArray* fMCArray; //!>! input MC stack
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager for reco mixed event
  AliEventPoolManager *fMCCorPoolMgr; //!>! correlation pool manager for MC truth mixed event


  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  THnSparse *fRecoHProton; //!>! reco h-proton correlation histogram
  THnSparse *fRecoHProtonMixed; //!>! reco h-proton correlation histogram for mixed event

  THnSparse *fRealHProton; //!>! real h-proton correlation histogram
  THnSparse *fRealHProtonMixed; //!>! real h-proton correlation histogram for mixed event

  void MakeSameHProtonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> proton_list, THnSparse* corr_dist, double bz, double zVtx);
  void MakeMixedHProtonCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> proton_list, THnSparse* corr_dist, double zVtx);

  void MakeSameMCHProtonCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> proton_list, THnSparse* corr_dist, double zVtx);
  void MakeMixedMCHProtonCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> proton_list, THnSparse* corr_dist, double zVtx);

  bool PassDaughterCuts(AliAODTrack *track); // check if the AOD track passes the daughter cuts
  bool PassTriggerCuts(AliAODTrack *track, bool checkMC = false); // check if the AOD track passes the trigger cuts 
  bool PassMCTriggerCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the trigger cuts
  bool IsMCChargedHadron(int pdg_code); // check if the MC particle is a charged hadron

  ClassDef(AliAnalysisTaskLambdaHadronCloseTrackEff, 3);

};
#endif
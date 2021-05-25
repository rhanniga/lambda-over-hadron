#ifndef AliAnalysisTaskLambdaHadronRatio_H
#define AliAnalysisTaskLambdaHadronRatio_H

//All relevant AliPhysics includes (this list will continue to grow):
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH3D.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"
#include "TFile.h"

//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskLambdaHadronRatio : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskLambdaHadronRatio();
  AliAnalysisTaskLambdaHadronRatio(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void TestPrint(TString string);
  void LoadEfficiencies(TFile* inputFile);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

 private:

  float MULT_LOW;
  float MULT_HIGH;
  float DAUGHTER_TRK_BIT;
  float ASSOC_TRK_BIT;
  float TRIG_TRK_BIT;

  TString CENT_ESTIMATOR;
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  
  TH1D* fTriggerEff; //!>! trigger efficiency
  TH1D* fAssociatedEff; //!>! associated efficiency
  TH1D* fLambdaEff; //!>! lambda efficiency

  TH2D* fTriggersAndLambdasPerEvent_All; //!>! triggers and all lambdas per event
  TH2D* fTriggersAndLambdasPerEvent_2_4; //!>! triggers and 2-4 GeV lambdas per event

  THnSparseF* fLooseDist;  //!>! single particle all hadron dist (no cuts at all)
  THnSparseF* fTriggerDist;  //!>! single particle trigger dist
  THnSparseF* fAssociatedHDist;  //!>! single particle associated hadron dist
  THnSparseF* fLambdaDist;  //!>! single particle lambda dist

  THnSparseF* fDphiHLambda;  //!>! hadron-lambda correlation hist
  THnSparseF* fDphiHLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparseF* fDphiHLambdaV0;  //!>! hadron-lambda correlation hist (using v0 finder for lambda)
  THnSparseF* fDphiHLambdaRotated;  //!>! hadron-lambda correlation hist with rotated pion
  THnSparseF* fDphiHLambdaRotatedPi;  //!>! hadron-lambda correlation hist with daughter rotated by pi
  THnSparseF* fDphiHLambdaRotatedProton;  //!>! hadron-lambda correlation hist with rotated proton
  THnSparseF* fDphiHLambdaFlipped;  //!>! hadron-lambda correlation hist with flipped pion
  THnSparseF* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparseF* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparseF* fDphiTriggerTrigger;   //!>! trigger-trigger correlation hist
  THnSparseF* fDphiHLambdaLS; //!>! hadron-proton+pion like sign correlation hist
  THnSparseF* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist
  THnSparseF* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparseF* fDphiHLambdaLSMixed; //!>! hadron-proton+pion like sign mixed correlation hist
  THnSparseF* fDphiTriggerTriggerMixed;   //!>! mixed trigger-trigger correlation hist

  THnSparseF* fLambdaDaughterDCA;

  // THnSparseF* fPid; //!>! histogram to visualize pid cuts
  // THnSparseF* fSignalAnalysis; //!>! histogram to analyze signal with nsigma cuts

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  //hand written functions:

  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist);
  void FillSingleParticleDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist);
  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskLambdaHadronRatio, 0);

};
#endif
#ifndef AliAnalysisTaskStrangeCor_H
#define AliAnalysisTaskStrangeCor_H

//All relevant AliPhysics includes (this list will continue to grow):
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"

//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskStrangeCor : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskStrangeCor();
  AliAnalysisTaskStrangeCor(const char *name);
  virtual ~AliAnalysisTaskStrangeCor();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

 private:
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager

  THnSparseF* fTriggerDist;  //!>! single particle trigger dist
  THnSparseF* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparseF* fDphiHStrangePart;  //!>! hadron-kstar correlation hist
  THnSparseF* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparseF* fDphiHStrangePartLS; //!>! hadron-k pi like sign correlation hist
  THnSparseF* fDphiHStrangePartMixed; //!>! hadron-kstar mixed correlation hist
  THnSparseF* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparseF* fDphiHStrangePartLSMixed; //!>! hadron-k pi like sign mixed correlation hist

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  //hand written functions:

  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist);
  void MakeSameHStrangePartCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> k_star_list, THnSparse* fDphi, double zVtx);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx);
  void MakeMixedHStrangePartCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> k_star_list , THnSparse* fDphi, double zVtx);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx);

  ClassDef(AliAnalysisTaskStrangeCor, 0);

};


#endif

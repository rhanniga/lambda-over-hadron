/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronMC class
                This is the header file for the lambda-hadron correlation analysis in MC (model comparisons, using MC stack only)
                (can be found in same directory)
                Origin: Ryan Hannigan, May 2023, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronMC_H
#define AliAnalysisTaskLambdaHadronMC_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"

#include "AliAnalysisTaskSE.h"

// Forward declaration order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH1F;
class THnSparse;
class TClonesArray;

class AliEventPoolManager;
class AliEventPool;
class AliMCEvent;
class AliAODMCParticle;


class AliAnalysisTaskLambdaHadronMC : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronMC();
  AliAnalysisTaskLambdaHadronMC(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronMC();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

private:

  AliMCEvent* fMCEvent; //!>! input MC stack
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fMCCorPoolMgr; //!>! correlation pool manager for MC mixed event

  TH1D* fMultDist; //!>! distribution for counting charged particles in V0A acceptance (for multiplicity selection)

  THnSparse* fTriggerDist_MC;  //!>! single particle trigger dist (MC truth)
  THnSparse* fAssociatedDist_MC;  //!>! single particle associated hadron dist (MC truth)
  THnSparse* fLambdaDist_MC;  //!>! single particle lambda dist (MC truth)
  THnSparse* fPhiDist_MC;  //!>! single particle phi (1020) dist (MC truth)

  THnSparse* fTriggerDist_MC_no_eta_cut;  //!>! single particle trigger dist (MC_no_eta_cut truth)
  THnSparse* fAssociatedDist_MC_no_eta_cut;  //!>! single particle associated hadron dist (MC_no_eta_cut truth)
  THnSparse* fLambdaDist_MC_no_eta_cut;  //!>! single particle lambda dist (MC_no_eta_cut truth)
  THnSparse* fPhiDist_MC_no_eta_cut;  //!>! single particle phi (1020) dist (MC_no_eta_cut truth)

  THnSparse* fTriggeredTriggerDist_MC;  //!>! single trigger hadron dist within a triggered event (MC truth)
  THnSparse* fTriggeredAssociatedDist_MC;  //!>! single associated hadron dist within a triggered event (MC truth)
  THnSparse* fTriggeredLambdaDist_MC;  //!>! single particle lambda dist within a triggered event (MC truth)
  THnSparse* fTriggeredPhiDist_MC;  //!>! single particle phi (1020) dist within a triggered event (MC truth)

  THnSparse* fDphiHH_MC;   //!>! hadron-hadron correlation hist (MC truth)
  THnSparse* fDphiHLambda_MC;  //!>! hadron-lambda correlation hist (MC truth)
  THnSparse* fDphiHPhi_MC;  //!>! hadron-phi correlation hist (MC truth)

  THnSparse* fDphiHH_MC_no_eta_cut;   //!>! hadron-hadron correlation hist (MC truth)
  THnSparse* fDphiHLambda_MC_no_eta_cut;  //!>! hadron-lambda correlation hist (MC truth)
  THnSparse* fDphiHPhi_MC_no_eta_cut;  //!>! hadron-phi correlation hist (MC truth)

  THnSparse* fDphiHHMixed_MC; //!>! hadron-hadron mixed correlation hist (MC truth)
  THnSparse* fDphiHLambdaMixed_MC; //!>! hadron-lambda mixed correlation hist (MC truth)
  THnSparse* fDphiHPhiMixed_MC; //!>! hadron-phi mixed correlation hist (MC truth)


  void FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist, int multBin);
  void MakeSameMCCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_list, THnSparse* fDphi, double zVtx, int multBin);
  void MakeMixedMCCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> associated_list, THnSparse* fDphi, double zVtx, int multBin);

  bool PassMCTriggerCuts(AliAODMCParticle *mc_track, bool etaCut = true); // check if the MC particle passes the trigger cuts
  bool PassMCAssociatedCuts(AliAODMCParticle *mc_track, bool etaCut = true); // check if the MC particle passes the associated cuts
  bool PassMCLambdaCuts(AliAODMCParticle *mc_track, bool checkPhysicalPrimary = false, bool etaCut = true); // check if the MC particle passes the lambda cuts 
  bool PassMCPhiCuts(AliAODMCParticle *mc_track, bool etaCut = true); // check if the MC particle passes the phi cuts
  bool IsMCChargedHadron(int pdg_code); // check if the MC particle is a charged hadron

  int GetMultBin(int numTracksInV0A);

  ClassDef(AliAnalysisTaskLambdaHadronMC, 3);

};
#endif
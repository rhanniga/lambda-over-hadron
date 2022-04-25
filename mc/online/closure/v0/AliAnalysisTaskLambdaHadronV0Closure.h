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
class TClonesArray;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;
class AliAODMCParticle;
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
    int motherLabel;
  };

private:
  float fMultLow; // lower bound for multiplicity
  float fMultHigh; // upper bound for multiplicity
  float fDaughterBit; // filter bit for daughter particle
  float fAssociatedBit; // filter bit for associated particle
  float fTriggerBit; // filter bit for trigger particle

  TString fCentEstimator;

  AliAODEvent* fAOD; //!>! input event
  TClonesArray* fMCArray; //!>! input MC stack
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  AliEventPoolManager *fMCCorPoolMgr; //!>! correlation pool manager for MC mixed event
  
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
  TH1D* fLambdaEff; ///> lambda efficiency

  THnSparse* fTriggerDistEff;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDistEff_checkMC;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDist;  //!>! single particle trigger dist (not corrected for efficiency, is MC hadron, is MC primary)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist
  THnSparse* fAssociatedHDist_checkMC;  //!>! single particle associated hadron dist (is MC hadron, is MC primary)
  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event
  THnSparse* fLambdaDist;  //!>! single particle lambda dist
  THnSparse* fGuaranteedLambdaDist;  //!>! single particle lambda dist (guaranteed to be lambda from MC pdg)

  THnSparse* fTriggerDist_MC;  //!>! single particle trigger dist (MC truth)
  THnSparse* fAssociatedDist_MC;  //!>! single particle associated hadron dist (MC truth)
  THnSparse* fLambdaDist_MC;  //!>! single particle lambda dist within a triggered event (MC truth)
  THnSparse* fTriggeredLambdaDist_MC;  //!>! single particle lambda dist within a triggered event (MC truth)

  THnSparse* fDphiHLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparse* fDphiHGuaranteedLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected, guaranteed lambda from MC pdg)
  THnSparse* fDphiHLambdaEff_MCKin;  //!>! hadron-lambda correlation hist (efficiency corrected, using MC kinematics)
  THnSparse* fDphiHDaughterProton_MCKin;  //!>! hadron-daughter proton correlation hist
  THnSparse* fDphiHDaughterPion_MCKin;  //!>! hadron-daughter pion correlation hist 
  THnSparse* fDphiHProton;  //!>! hadron-proton correlation hist
  THnSparse* fDphiHPion;  //!>! hadron-pion correlation hist 
  THnSparse* fDphiHLambdaEff_MCKin_physicalPrimary;  //!>! hadron-lambda correlation hist (efficiency corrected, using MC kinematics, trigger and lambda are physical primary)
  THnSparse* fDphiRecoHRealLambdaEff_MCKin_physicalPrimary; //!>! hadron-lambda correlation hist (efficiency corrected, using MC kinematics, trigger and lambda are physical primary)
  THnSparse* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHHEff_checkMC;   //!>! hadron-hadron correlation hist (efficiency corrected, trig and assoc are MC primary)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHLambdaMixed_MCKin; //!>! hadron-lambda mixed correlation hist (MC kinematics)
  THnSparse* fDphiHDaughterProtonMixed_MCKin;  //!>! hadron-daughter proton mixed correlation hist
  THnSparse* fDphiHDaughterPionMixed_MCKin;  //!>! hadron-daughter pion mixed correlation hist 
  THnSparse* fDphiHProtonMixed;  //!>! hadron-proton correlation hist (mixed)
  THnSparse* fDphiHPionMixed;  //!>! hadron-pion correlation hist (mixed)
  THnSparse* fDphiHLambdaMixed_MCKin_physicalPrimary; //!>! hadron-lambda mixed correlation hist (MC kinematics, trigger and lambda are physical primary)
  THnSparse* fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary; //!>! hadron-lambda mixed correlation hist (MC kinematics, trigger and lambda are physical primary) 

  THnSparse* fDphiHLambda_MC;  //!>! hadron-lambda correlation hist (MC truth)
  THnSparse* fDphiHDaughterProton_MC;  //!>! hadron-daughter proton correlation hist (MC truth)
  THnSparse* fDphiHDaughterPion_MC;  //!>! hadron-daughter pion correlation hist (MC truth)
  THnSparse* fDphiHProton_MC;  //!>! hadron-proton correlation hist (MC truth)
  THnSparse* fDphiHPion_MC;  //!>! hadron-pion correlation hist (MC truth)
  THnSparse* fDphiHLambda_MC_physicalPrimary;  //!>! hadron-lambda correlation hist (MC truth, lambda is physical primary)
  THnSparse* fDphiHH_MC;   //!>! hadron-hadron correlation hist (MC truth)
  THnSparse* fDphiHLambdaMixed_MC; //!>! hadron-lambda mixed correlation hist (MC truth)
  THnSparse* fDphiHDaughterProtonMixed_MC; //!>! hadron-daughter proton mixed correlation hist (MC truth)
  THnSparse* fDphiHDaughterPionMixed_MC; //!>! hadron-daughter pion mixed correlation hist (MC truth)
  THnSparse* fDphiHProtonMixed_MC;  //!>! hadron-proton mixed correlation hist (MC truth)
  THnSparse* fDphiHPionMixed_MC;  //!>! hadron-pion mixed correlation hist (MC truth)
  THnSparse* fDphiHLambdaMixed_MC_physicalPrimary; //!>! hadron-lambda mixed correlation hist (MC truth, lambda is physical primary)
  THnSparse* fDphiHHMixed_MC; //!>! hadron-hadron mixed correlation hist (MC truth)

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=false);
  void FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda, bool lambdaEff=true);
  void FillMCMotherDist(std::vector<AliAODMCParticle*> particle_list, float multPercentile, THnSparse* fDist);

  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda);
  void MakeSameRecoHRealLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeSameHLambdaCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeSameHDaughterCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda);
  void MakeMixedHLambdaCorrelations_withMCKin(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeMixedHDaughterCorrelations_withMCKin(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);

  void MakeSameMCHLambdaCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx);
  void MakeSameMCHDaughterCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx);
  void MakeSameMCHHCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx);
  void MakeMixedMCHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx);
  void MakeMixedMCHDaughterCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx);
  void MakeMixedMCHHCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> associated_h_list , THnSparse* fDphi, double zVtx);

  bool PassDaughterCuts(AliAODTrack *track); // check if the AOD track passes the daughter cuts
  bool PassTriggerCuts(AliAODTrack *track, bool checkMC = false); // check if the AOD track passes the trigger cuts 
  bool PassAssociatedCuts(AliAODTrack *track, bool checkMC = false); // check if the AOD track passes the associated cuts
  uint8_t PassV0LambdaCuts(AliAODv0 *v0, bool checkMotherPDG = false, bool checkPhysicalPrimary = false); // check if the AOD v0 passes the lambda cuts (0 = no, 1 = yes and lambda, 2 = yes and anti-lambda, checkMotherPDG verifies if the mother is a lambda or anti-lambda)
  bool PassMCTriggerCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the trigger cuts
  bool PassMCAssociatedCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the associated cuts
  bool PassMCLambdaCuts(AliAODMCParticle *mc_track, bool checkPhysicalPrimary = false); // check if the MC particle passes the lambda cuts 
  bool IsMCChargedHadron(int pdg_code); // check if the MC particle is a charged hadron

  ClassDef(AliAnalysisTaskLambdaHadronV0Closure, 3);

};
#endif
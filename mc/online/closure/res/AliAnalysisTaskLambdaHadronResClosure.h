/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskLambdaHadronResClosure class
                This is a closure task for the lambda-hadron correlation analysis
                (can be found in same directory)
                Origin: Ryan Hannigan, March 2022, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskLambdaHadronResClosure_H
#define AliAnalysisTaskLambdaHadronResClosure_H

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


class AliAnalysisTaskLambdaHadronResClosure : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaHadronResClosure();
  AliAnalysisTaskLambdaHadronResClosure(const char *name);
  virtual ~AliAnalysisTaskLambdaHadronResClosure();
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
  THnSparse* fLambdaDist;  //!>! single particle lambda dist
  THnSparse* fLambdaLSDist;  //!>! single particle lambda dist (LS daughters)
  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event
  THnSparse* fTriggeredLambdaLSDist;  //!>! single particle lambda dist within a triggered event (LS daughters)

  THnSparse* fTriggerDist_MC;  //!>! single particle trigger dist (MC truth)
  THnSparse* fAssociatedDist_MC;  //!>! single particle associated hadron dist (MC truth)
  THnSparse* fLambdaDist_MC;  //!>! single particle lambda dist (MC truth)
  THnSparse* fTriggeredLambdaDist_MC;  //!>! single particle lambda dist within a triggered event (MC truth)

  THnSparse* fDphiHLambdaEff;  //!>! hadron-lambda correlation hist (efficiency corrected)
  THnSparse* fDphiHLambdaEff_MCKin;  //!>! hadron-lambda correlation hist (efficiency corrected, using MC kinematics)
  THnSparse* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHHEff_checkMC;   //!>! hadron-hadron correlation hist (efficiency corrected, trig and assoc are MC primary)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist
  THnSparse* fDphiHLambdaMixed_MCKin; //!>! hadron-lambda mixed correlation hist (using MC kinematics)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist

  THnSparse* fDphiHLambda_MC;  //!>! hadron-lambda correlation hist (MC truth)
  THnSparse* fDphiHH_MC;   //!>! hadron-hadron correlation hist (MC truth)
  THnSparse* fDphiHLambdaMixed_MC; //!>! hadron-lambda mixed correlation hist (MC truth)
  THnSparse* fDphiHHMixed_MC; //!>! hadron-hadron mixed correlation hist (MC truth)

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=false);
  void FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool lambdaEff=true);
  void FillMCMotherDist(std::vector<AliAODMCParticle*> particle_list, float multPercentile, THnSparse* fDist);

  void MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeSameHLambdaCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda);
  void MakeMixedHLambdaCorrelations_withMCKin(AliEventPool *fPool, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);

  void MakeSameMCHLambdaCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx);
  void MakeSameMCHHCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx);
  void MakeMixedMCHLambdaCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx);
  void MakeMixedMCHHCorrelations(AliEventPool *fPool, std::vector<AliAODMCParticle*> associated_h_list , THnSparse* fDphi, double zVtx);

  bool PassDaughterCuts(AliAODTrack *track); // check if the AOD track passes the daughter cuts
  bool PassTriggerCuts(AliAODTrack *track, bool checkMC = false); // check if the AOD track passes the trigger cuts 
  bool PassAssociatedCuts(AliAODTrack *track, bool checkMC = false); // check if the AOD track passes the associated cuts
  uint8_t PassV0LambdaCuts(AliAODv0 *v0); // check if the AOD v0 passes the lambda cuts (0 = no, 1 = yes and lambda, 2 = yes and anti-lambda)
  bool PassMCTriggerCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the trigger cuts
  bool PassMCAssociatedCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the associated cuts
  bool PassMCLambdaCuts(AliAODMCParticle *mc_track); // check if the MC particle passes the lambda cuts 
  bool IsMCChargedHadron(int pdg_code); // check if the MC particle is a charged hadron

  ClassDef(AliAnalysisTaskLambdaHadronResClosure, 3);

};
#endif
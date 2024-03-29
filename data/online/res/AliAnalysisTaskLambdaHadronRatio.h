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

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection
  
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
  TH1D* fLambdaEff; ///> lambda efficiency

  THnSparse* fTriggerDist;  //!>! single particle trigger dist (corrected for efficiency)
  THnSparse* fTriggerDist_highestPt;  //!>! single particle trigger dist (corrected for efficiency, highest pt between 4 and 8)
  THnSparse* fAssociatedHDist;  //!>! single particle associated hadron dist

  THnSparse* fTotalTrackDist;  //!>! single particle total track dist (all AOD tracks in list)
  THnSparse* fCutTrackDist;  //!>! single particle total track dist (AOD tracks that pass daughter cuts)

  THnSparse* fLambdaDist;  //!>! single particle lambda dist (corrected for efficiency)
  THnSparse* fLambdaLSDist;  //!>! single particle lambda dist with like-sign daughters (corrected for efficiency)
  THnSparse* fTriggeredLambdaDist;  //!>! single particle lambda dist within a triggered event (corrected for efficiency)
  THnSparse* fTriggeredLambdaLSDist;  //!>! single particle lambda dist with ligk-sign daughters within a triggered event (corrected for efficiency)

  THnSparse* fDphiHLambda;  //!>! hadron-lambda correlation hist (corrected for efficiency)
  THnSparse* fDphiHLambda_highestPt;  //!>! hadron-lambda correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHH;   //!>! hadron-hadron correlation hist (corrected for efficiency)
  THnSparse* fDphiHH_highestPt;   //!>! hadron-hadron correlation hist (efficiency corrected, highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaMixed_highestPt; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHLambdaMixed; //!>! hadron-lambda mixed correlation hist (highest pt trigger between 4 and 8)
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHHMixed_highestPt; //!>! hadron-hadron mixed correlation hist (highest pt trigger between 4 and 8)

  TH2D*      fTPCnSigmaProton; //!>! TPC n sigma for proton
  TH2D*      fTPCnSigmaPion; //!>! TPC n sigma for pion
  TH2D*      fTOFnSigmaProton; //!>! TOF n sigma for proton
  TH2D*      fTOFnSigmaPion; //!>! TOF n sigma for pion
  TH2D*      fTOFvTPCnSigmaProton; //!>! TOF v TPC n sigma for proton (TPC is X axis, TOF is Y axis)
  TH2D*      fTOFvTPCnSigmaPion; //!>! TOF v TPCn sigma for pion (TPC is X axis, TOF is Y axis)



  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=true);
  void FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool lambda_eff=true);
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
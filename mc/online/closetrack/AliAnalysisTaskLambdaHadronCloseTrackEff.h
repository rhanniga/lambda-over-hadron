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
#include "TArrayF.h"
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


class AliMixingParticle : public AliVParticle {
  public: 
    AliMixingParticle() : AliVParticle(), fPt(0), fEta(0), fPhi(0), fCharge(0), fPosition(0) { }
    AliMixingParticle(Float_t pt, Float_t eta, Float_t phi, Short_t charge, Double_t *position)
    : AliVParticle(), fPt(pt), fEta(eta), fPhi(phi), fCharge(charge), fPosition(position) { }
    virtual ~AliMixingParticle() { }
    
    virtual Double_t Pt()    const { return fPt;      }
    virtual Double_t Phi()   const { return fPhi;     }
    virtual Double_t Eta()   const { return fEta;     }
    virtual Short_t Charge() const { return fCharge;  }
    virtual Double_t Px()    const { return fPt*TMath::Cos(fPhi); }
    virtual Double_t Py()    const { return fPt*TMath::Sin(fPhi); }
    virtual Double_t Pz()    const { return fPt*TMath::SinH(fEta); }
    virtual Double_t Xv() const { return fPosition[0]; }
    virtual Double_t Yv() const { return fPosition[1]; }
    virtual Double_t Zv() const { return fPosition[2]; }

    virtual Double_t P()  const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

    
  private:
    Float_t fPt;
    Float_t fEta;
    Float_t fPhi;
    Float_t fCharge;

    Double_t *fPosition;

    ClassDef(AliMixingParticle, 1);
};

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

  TH1D *fMinDistanceNotFound; //!>! minimum distance of MC tracks not found in reco
  TH1D *fMinDistanceAll; //!>! minimum distance all MC tracks

  bool GetXYZatR(AliAODTrack *track, float radius, double *xyz);
  bool GetXYZatR(AliAODMCParticle *track, float radius, double *xyz);
  bool GetXYZatR(AliMixingParticle *track, float radius, double *xyz);

  void MakeSameHProtonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> proton_list, THnSparse* corr_dist, double zVtx);
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
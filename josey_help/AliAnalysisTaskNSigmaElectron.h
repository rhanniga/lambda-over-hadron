#ifndef AliAnalysisTaskNSigmaElectron_H
#define AliAnalysisTaskNSigmaElectron_H
class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliEventPoolManager;
class AliMultSelection;
class AliAnalysisUtils;
class AliGenEventHeader;
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
#include "AliLog.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"


//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskNSigmaElectron : public AliAnalysisTaskSE {

    
 public:
    
    enum EnhanceSigOrNot {kMB,kEnhance};
    enum pi0etaType {kNoMother, kNoFeedDown, kNotIsPrimary, kLightMesons, kBeauty, kCharm};
    
  AliAnalysisTaskNSigmaElectron();
  AliAnalysisTaskNSigmaElectron(const char *name);
  virtual ~AliAnalysisTaskNSigmaElectron();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

 private:

  TString CENT_ESTIMATOR;
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  THnSparseF* fBigDist;   //!>! big dist to do christina's bidding

  TH2D *fNSigmaElectronBethe; //!>! hopeful dist
  TH2D *fNSigmaPionBethe; //!>! hopeful dist
  TH2D *fNSigmaKaonBethe; //!>! hopeful dist
  TH2D *fNSigmaProtonBethe; //!>! hopeful dist

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  bool PassAssociatedCutsElectrons(AliAODTrack *track);

  ClassDef(AliAnalysisTaskNSigmaElectron, 3);

};
#endif

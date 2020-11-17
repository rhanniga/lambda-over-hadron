#ifndef AliAnalysisTaskLambdaHadronEfficiency_cxx
#define AliAnalysisTaskLambdaHadronEfficiency_cxx

//QA task for EMCAL electron analysis

#include "AliAnalysisTaskSE.h"
#include "AliAODMCHeader.h"
#include "THnSparse.h"
#include "TObject.h"
#include "TRandom3.h"
#include "AliAODTrack.h"

class TH1F;
class AliEventPoolManager;
class THnSparse;
class AliESDEvent;
class AliAODEvent;
class AliAODMCParticle;
class AliMultSelection;

class AliAnalysisTaskLambdaHadronEfficiency : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskLambdaHadronEfficiency();
    AliAnalysisTaskLambdaHadronEfficiency(const char *name, Float_t multLow, Float_t multHigh);
    virtual ~AliAnalysisTaskLambdaHadronEfficiency();
    
    virtual void   UserCreateOutputObjects();
    UInt_t PassProtonCuts(AliAODTrack* track, Double_t TPCnSigma, Double_t TOFnSigma);
    UInt_t PassPionCuts(AliAODTrack* track, Double_t TPCnSigma, Double_t TOFnSigma);
    Bool_t PassHadronCuts(AliAODTrack* track, Bool_t isTrigger);
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    
    void SetKaonTrkBit(Int_t kaonbit){ KAON_TRK_BIT = kaonbit; };
    void SetKaonEtaCut(Float_t eta){ KAON_ETA_CUT = eta; };
    void SetCentEstimator(TString est){ CENT_ESTIMATOR = est; };
    void SetTrigTrkBit(UInt_t trkbit) { TRIG_TRK_BIT = trkbit; };
    void SetAssocTrkBit(UInt_t trkbit) { ASSOC_TRK_BIT = trkbit; };

    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    
private:

    Float_t MULT_LOW;
    Float_t MULT_HIGH;
    Float_t KAON_ETA_CUT;
    Float_t KAON_TRK_BIT;
    Float_t ASSOC_TRK_BIT;
    Float_t TRIG_TRK_BIT;
    TString CENT_ESTIMATOR;

    UInt_t TRACK_BIT = 1UL << 0;
    UInt_t TOF_HIT_BIT = 1UL << 1;
    UInt_t TPC_PID_BIT = 1UL << 2;
    UInt_t TOF_PID_BIT = 1UL << 3;

    enum{
        kAODanalysis = BIT(20),
    };
      
    AliVEvent   *fVevent;  //!event object
    AliEventPoolManager *fPoolMgr; //! Event pool manager for mixed event
    AliEventPoolManager *fLSPoolMgr; //! Event pool manager for LS mixed event
    AliEventPoolManager *fHHPoolMgr; //! Event pool manager for HH
    AliESDEvent *fESD;    //!ESD object
    AliAODEvent *fAOD;    //!AOD object
    AliPIDResponse *fpidResponse; //!pid response
    AliMultSelection *fMultSelection; //!mult selection
   
    TClonesArray* fMCArray; //!
    AliAODMCHeader* fMCHeader; //!
    TList       *fOutputList; //!Output list
    TH1F        *fNevents;//! no of events
    TH1F        *fNumTracks;//! number of Tracks/evt
    TH1F        *fVtxZ;//!Vertex z
    TH1F        *fVtxX;//!Vertex x
    TH1F        *fVtxY;//!Vertex y

    THnSparseF  *fRealLambdaDist;//! Dist of Real lambda
    THnSparseF  *fRealNoDecayCutLambdaDist;//! Dist of Real lambda with no check on decay daughter eta
    THnSparseF  *fRecoLambdaDist;//! Dist of Recon lambda
    THnSparseF  *fRecoNormalLambdaDist;//! Dist of Recon normal lambda
    THnSparseF  *fRecoAntiLambdaDist;//! Dist of Recon anti lambda
    THnSparseF  *fTrackRecoLambdaDist;//! Dist of Recon lambda passing track cuts
    THnSparseF  *fTOFRecoLambdaDist;//! Dist of Recon lambda passing track cuts + TOF hit
    THnSparseF  *fTPCPIDTrackRecoLambdaDist;//! Dist of Recon lambda passing track cuts + TPC PID
    THnSparseF  *fTPCPIDRecoLambdaDist;//! Dist of Recon lambda passing track cuts + TOF hit + TPC PID
    THnSparseF  *fPIDRecoLambdaDist;//! Dist of Recon lambda passing track cuts + TOF&TPC PID 3 sigma cut

    THnSparseF  *fRealChargedDist;//!
    THnSparseF  *fRealKDist;//!
    THnSparseF  *fRealPiDist;//!
    THnSparseF  *fRealeDist;//!
    THnSparseF  *fRealpDist;//!
    THnSparseF  *fRealMuonDist;//!
   
    THnSparseF  *fRecoChargedDist;//!
    THnSparseF  *fRecoKDist;//!
    THnSparseF  *fTOFPiDist;//!
    THnSparseF  *fTOFProtonDist;//!
    THnSparseF  *fRecoPiDist;//!
    THnSparseF  *fRecoeDist;//!
    THnSparseF  *fRecopDist;//!
    THnSparseF  *fRecoMuonDist;//!

    THnSparseF  *fRecoChargedTriggerDist;//!
    THnSparseF  *fRecoKTriggerDist;//!
    THnSparseF  *fRecoPiTriggerDist;//!
    THnSparseF  *fRecoeTriggerDist;//!
    THnSparseF  *fRecopTriggerDist;//!
    THnSparseF  *fRecoMuonTriggerDist;//!

    TH1D        *fReactionPlane;//!


    ClassDef(AliAnalysisTaskLambdaHadronEfficiency, 1); 
};

#endif
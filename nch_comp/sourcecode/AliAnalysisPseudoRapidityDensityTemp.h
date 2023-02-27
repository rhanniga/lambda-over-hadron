/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
//===========================================================

#ifndef ALIANALYSISPSEUDORAPIDITYDENSITYTEMP_H
#define ALIANALYSISPSEUDORAPIDITYDENSITYTEMP_H

#include "TParticle.h"
#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "THistManager.h"
#include <deque>
#include "AliAODForwardMult.h"

class AliVMultiplicity;
class TTree;
class AliAnalysisUtils;
class AliOADBMultSelection;
class AliStack;
class TRandom3;
class AliMultSelection;
class AliMultSelectionTask;
class AliITSMultRecBg;
//class TChain;

class AliAnalysisPseudoRapidityDensityTemp : public AliAnalysisTaskSE
{
public:
    typedef std::vector<Bool_t> Bool_1d;
    typedef std::vector<Double_t> Double1D;

    enum
    {
        kECbegin = 1,
        kDATA = 1,
        kINEL,
        kINELg0,
        kINEL300,
        kINEL100,
        kINEL200,
	    kINEL015,
        kINEL040,
        kINEL050,
        kINEL04025,
        kINEL05024,
        kINEL05025,
        kECend
    };
    enum
    {
        kTrigbegin = 1,
        kHMMBAND=1,
        kHMMBAND300,
        kMBAND015,
        kMBAND300,
        kMBAND100,
        kMBAND200,
        kMBAND040,
        kMBAND050,
        kTrigend
    };
    enum
    {
        kParTypebegin = 1,
        kParDATA = 1,
        kMotherStrange,
        kBkg,
        kPion,
        kKaon,
        kProtonBK,
        kOPar,
        kParTypeend
    };
    enum
    {
        kV0Typebegin = 1,
        kK0s = 1,
        kLambda,
        kAntiLambda,
        kV0Typeend
    };
    enum
    {
        kNoTrkCutVar = 1,
        kHybrid = 1 ,
        kITSTPC2011,
        kITSTPC2011dcazdw, 
        kITSTPC2011dcazup,
        kITSTPC2011dcardw, 
        kITSTPC2011dcarup,
        kITSTPC2011nclutpcdw, 
        kITSTPC2011nclutpcup,
        kITSTPC2011chitpcdw, 
        kITSTPC2011chitpcup,
        kITSTPC2011globalconsdw, 
        kITSTPC2011globalconsup,
        kTrackCutend
    };
    enum
    {
        kCentClassBinBegin = 1,
        kV0M = 1,
        kCentClassBinEnd
    };

    AliAnalysisPseudoRapidityDensityTemp();
    AliAnalysisPseudoRapidityDensityTemp(
        const char *name, const char *option);

    AliAnalysisPseudoRapidityDensityTemp(
        const AliAnalysisPseudoRapidityDensityTemp &ap);

    AliAnalysisPseudoRapidityDensityTemp &operator=(
        const AliAnalysisPseudoRapidityDensityTemp &ap);

    ~AliAnalysisPseudoRapidityDensityTemp();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    void SetOption(char *option) { fOption = option; }
    void SetFilterBit(UInt_t filterbit) { fFilterBit = filterbit; }
    Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk);

    void FillTracklets(Bool_1d bevtc, Bool_1d btrigc);
    void FillTracks(Bool_1d bevtc, Bool_1d btrigc);
    void FillFMDAODtracks();
    void MeasureDiffMass();
    Bool_t FillFMDESDtracks();
    void StrangenessMeasure(Bool_1d btrigc);
    void SetIsAA(Bool_t isaa) { IsAA = isaa; }
    void SetGeoCentHist(std::vector<TH1D> &h)
    {
        fMultGeo = h;
    }
    void SetCentHist(std::vector<TH1D> &h)
    {
        fMult = h;
    }
    void InitMultReco();
    TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
    TAxis AxisVar(TString name, std::vector<Double_t> bin);
    TAxis AxisLog(TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0);
    TAxis AxisStr(TString name, std::vector<TString> bin);
    THnSparse *CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt = "");
    THnSparse *CreateTHnSparse(TString name, TString title, TString templ, Option_t *opt = "");
    Long64_t FillTHnSparse(TString name, std::vector<Double_t> x, Double_t w = 1.);
    Long64_t FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w = 1.);
    Bool_t IsMCTRDtriggered();
    Float_t GetInvPtDevFromBC(Int_t b, Int_t c);

    Bool_t IsMotherStrangeParticle(TParticle *mother)
    {
        if (TMath::Abs(mother->GetPdgCode()) == 3122   //Lambda
            || TMath::Abs(mother->GetPdgCode()) == 310 //K0s
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 3334 //Omega
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 3322 //Xi0
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 3312 //Xi
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 3222 //Sigma
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 3112 //Sigma
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 333 //PHi
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 130 //k0L
                                                       //|| TMath::Abs(mother->GetPdgCode()) == 321 //K

        )
            return kTRUE;
        else
            return kFALSE;
    }

private:
    typedef std::vector<AliVTrack *> tracklist;
    typedef std::deque<tracklist> eventpool;
    typedef std::vector<std::vector<eventpool>> mixingpool;
    typedef std::vector<Int_t> Int_1d;
    typedef std::vector<Double_t> Double_1d;

    TString fOption;

    AliTriggerAnalysis *fTrigger = nullptr; //!
    std::vector<AliESDtrackCuts> fTrackCuts;
    AliESDtrackCuts fTrackCutGC;
    AliVEvent *fEvt = nullptr; //!
    UInt_t fFilterBit;
    Bool_t IsFirstEvent = kTRUE;

    Double_1d fCent;
    Double_t fZ = -30;

    AliPIDResponse *fPIDResponse = nullptr;             //!
    AliPIDCombined *fPIDCombined = nullptr;             //!
    AliAnalysisUtils *fUtils = nullptr;                 //!
    AliOADBMultSelection *fOadbMultSelection = nullptr; //!
    AliMultSelection *sel = nullptr;                    //!
    AliMultSelectionTask *seltask = nullptr;            //!

    //Histograms below are main
    std::vector<std::vector<TH2D *>> fMass2D; //! signbins, centbins
    //Histograms for pT_pair amd pT

    AliStack *stack = nullptr; //!
    Bool_t ismc = false;
    Double_t fptcut = 0.2;
    Double_t fetacut = 0.9;
    Bool_t IsAA = kFALSE;
    THistManager *fHistos = nullptr;           //!
    AliVMultiplicity *fMultiplicity = nullptr; //!
    //TTree*                        fTree=nullptr;//!
    TH2 *fForward2d = nullptr;      //!
    TRandom3 *fRandom = nullptr;    //!
    TH1D *hdifftuneratio = nullptr; //!
    Double_t sdweightingfactor = 1.;
    std::vector<TH1D> fMultGeo;                //[kCentClassBinEnd]
    std::vector<TH1D> fMult;                   //[kCentClassBinEnd]
    AliITSMultRecBg *fMultReco = nullptr; //! mult.reco object
    TTree *fRPTree = nullptr;             //! tree of recpoints
    Bool_t isbginjection = false;

    ClassDef(AliAnalysisPseudoRapidityDensityTemp, 1);
};

#endif

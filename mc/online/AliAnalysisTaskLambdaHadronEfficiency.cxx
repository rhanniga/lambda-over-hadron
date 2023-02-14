/************************************************************************* 
 * Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. * 
 *                                                                        * 
 * Author: Ryan Hannigan
 * Contributors are mentioned in the code where appropriate.              * 
 *                                                                        * 
 * Permission to use, copy, modify and distribute this software and its   * 
 * documentation strictly for non-commercial purposes is hereby granted   * 
 * without fee, provided that the above copyright notice appears in all   * 
 * copies and that both the copyright notice and this permission notice   * 
 * appear in the supporting documentation. The authors make no claims     * 
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  * 
 **************************************************************************/

#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCFParticle.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskLambdaHadronEfficiency.h"

ClassImp(AliAnalysisTaskLambdaHadronEfficiency)

AliAnalysisTaskLambdaHadronEfficiency::AliAnalysisTaskLambdaHadronEfficiency(const char *name, Float_t multLow, Float_t multHigh)
    : AliAnalysisTaskSE(name),
    fVevent(0x0),
    fESD(0x0),
    fAOD(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fOutputList(0x0),
    fNevents(0x0),
    fNumTracks(0x0),
    fVtxZ(0x0),
    fVtxX(0x0),
    fVtxY(0x0),
    fRealTotalLambdaDist(0x0),
    fRealLambdaDist(0x0),
    fRealAntiLambdaDist(0x0),
    fRecoTotalV0LambdaDist(0x0),
    fRecoEtaV0LambdaDist(0x0),
    fRecoEtaPtV0LambdaDist(0x0),
    fRecoEtaPtRefitV0LambdaDist(0x0),
    fRecoEtaPtRefitRowsV0LambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioV0LambdaDist(0x0),
    fRecoTotalLambdaDist(0x0),
    fRecoEtaLambdaDist(0x0),
    fRecoEtaPtLambdaDist(0x0),
    fRecoEtaPtRefitLambdaDist(0x0),
    fRecoEtaPtRefitRowsLambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaFilterDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaDCADist(0x0),
    fRecoTotalLambdaFilterDist(0x0),
    fRecoLambdaDist(0x0),
    fRecoLambdaWithAODPionDist(0x0),
    fRecoLambdaWithAODProtonDist(0x0),
    fRecoAntiLambdaDist(0x0),
    fRealChargedDist(0x0),
    fRealTriggerDist(0x0),
    fRealKDist(0x0),
    fRealPiDist(0x0),
    fRealPiFromLambdaDist(0x0),
    fRealeDist(0x0),
    fRealpDist(0x0),
    fRealpFromLambdaDist(0x0),
    fRealMuonDist(0x0),
    fRecoChargedDist(0x0),
    fRecoKDist(0x0),
    fRecoPiDist(0x0),
    fRecoeDist(0x0),
    fRecopDist(0x0),
    fRecoMuonDist(0x0),
    fRecoChargedTriggerDist(0x0),
    fRecoKTriggerDist(0x0),
    fRecoPiTriggerDist(0x0),
    fRecoeTriggerDist(0x0),
    fRecopTriggerDist(0x0),
    fRecoMuonTriggerDist(0x0),
    fRealLambdaDaughterDist(0x0),
    fRecoLambdaDaughterDist(0x0),
    fRealLambdasPerEvent(0x0),
    fRecoLambdasPerEvent(0x0),
    fReactionPlane(0x0),
    fPxDifference(0x0),
    fPxDifferenceFB(0x0),
    fPyDifference(0x0),
    fPzDifference(0x0),
    fPtDifference(0x0),
    fRealPrimaryLambdaPtDist(0x0),
    fRealSecondaryLambdaPtDist(0x0),
    fInvMassAntiLambdaReal(0x0),
    fInvMassLambdaReal(0x0),
    fInvMassLambdaResonance(0x0),
    fInvMassAntiLambdaResonance(0x0),
    fInvMassLambdaV0(0x0),
    fInvMassAntiLambdaV0(0x0),
    fPhiDifferenceResV0(0x0),
    fPhiDifferenceResReal(0x0),
    fPhiDifferenceV0Real(0x0),
    fPhiV0(0x0),
    fPhiRes(0x0),
    fPhiReal(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    MIN_CROSSED_ROWS_TPC = 70;
    MIN_ROW_CLUSTER_RATIO = 0.8;
    MULT_LOW = multLow;
    MULT_HIGH = multHigh;
    CENT_ESTIMATOR = "V0A";
    ASSOC_TRK_BIT = 1024; // global tracks with tight pt dependent dca cut (selecting primaries)
    TRIG_TRK_BIT = AliAODTrack::kIsHybridGCG; // = 2^20
    DAUGHTER_ETA_CUT = 0.8;
    DAUGHTER_MIN_PT = 0.15;
}

AliAnalysisTaskLambdaHadronEfficiency::AliAnalysisTaskLambdaHadronEfficiency()
    : AliAnalysisTaskSE("default_task"),
    fVevent(0x0),
    fESD(0x0),
    fAOD(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fOutputList(0x0),
    fNevents(0x0),
    fNumTracks(0x0),
    fVtxZ(0x0),
    fVtxX(0x0),
    fVtxY(0x0),
    fRealTotalLambdaDist(0x0),
    fRealLambdaDist(0x0),
    fRealAntiLambdaDist(0x0),
    fRecoTotalV0LambdaDist(0x0),
    fRecoEtaV0LambdaDist(0x0),
    fRecoEtaPtV0LambdaDist(0x0),
    fRecoEtaPtRefitV0LambdaDist(0x0),
    fRecoEtaPtRefitRowsV0LambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioV0LambdaDist(0x0),
    fRecoTotalLambdaDist(0x0),
    fRecoEtaLambdaDist(0x0),
    fRecoEtaPtLambdaDist(0x0),
    fRecoEtaPtRefitLambdaDist(0x0),
    fRecoEtaPtRefitRowsLambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaFilterDist(0x0),
    fRecoEtaPtRefitRowsRatioLambdaDCADist(0x0),
    fRecoTotalLambdaFilterDist(0x0),
    fRecoLambdaDist(0x0),
    fRecoLambdaWithAODPionDist(0x0),
    fRecoLambdaWithAODProtonDist(0x0),
    fRecoAntiLambdaDist(0x0),
    fRealChargedDist(0x0),
    fRealTriggerDist(0x0),
    fRealKDist(0x0),
    fRealPiDist(0x0),
    fRealPiFromLambdaDist(0x0),
    fRealeDist(0x0),
    fRealpDist(0x0),
    fRealpFromLambdaDist(0x0),
    fRealMuonDist(0x0),
    fRecoChargedDist(0x0),
    fRecoKDist(0x0),
    fRecoPiDist(0x0),
    fRecoeDist(0x0),
    fRecopDist(0x0),
    fRecoMuonDist(0x0),
    fRecoChargedTriggerDist(0x0),
    fRecoKTriggerDist(0x0),
    fRecoPiTriggerDist(0x0),
    fRecoeTriggerDist(0x0),
    fRecopTriggerDist(0x0),
    fRecoMuonTriggerDist(0x0),
    fRealLambdaDaughterDist(0x0),
    fRecoLambdaDaughterDist(0x0),
    fRealLambdasPerEvent(0x0),
    fRecoLambdasPerEvent(0x0),
    fReactionPlane(0x0),
    fPxDifference(0x0),
    fPxDifferenceFB(0x0),
    fPyDifference(0x0),
    fPzDifference(0x0),
    fPtDifference(0x0),
    fRealPrimaryLambdaPtDist(0x0),
    fRealSecondaryLambdaPtDist(0x0),
    fInvMassAntiLambdaReal(0x0),
    fInvMassLambdaReal(0x0),
    fInvMassLambdaResonance(0x0),
    fInvMassAntiLambdaResonance(0x0),
    fInvMassLambdaV0(0x0),
    fInvMassAntiLambdaV0(0x0),
    fPhiDifferenceResV0(0x0),
    fPhiDifferenceResReal(0x0),
    fPhiDifferenceV0Real(0x0),
    fPhiV0(0x0),
    fPhiRes(0x0),
    fPhiReal(0x0)

{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    MIN_CROSSED_ROWS_TPC = 70;
    MIN_ROW_CLUSTER_RATIO = 0.8;
    MULT_LOW = 0;
    MULT_HIGH = 20;
    CENT_ESTIMATOR = "V0A";
    ASSOC_TRK_BIT = 1024; // global tracks with tight pt dependent dca cut (selecting primaries)
    TRIG_TRK_BIT = AliAODTrack::kIsHybridGCG; // = 2^20
    DAUGHTER_ETA_CUT = 0.8;
    DAUGHTER_MIN_PT = 0.15;
}

AliAnalysisTaskLambdaHadronEfficiency::~AliAnalysisTaskLambdaHadronEfficiency()
{
    //Destructor
    delete fOutputList;
}

void AliAnalysisTaskLambdaHadronEfficiency::UserCreateOutputObjects()
{
    //printf("\n!!!!!\n Starting UserCreateOutputObjects \n\n");
    //fflush(stdout);
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
   
    ////////////////////////////
    // Set-up for Mixed Event //
    ////////////////////////////


    Int_t numVtxZBins = 10;
    Double_t vtxZBins[11] = {-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
    Int_t numMultBins = 3;
    Double_t multBins[4] = {0.0, 20.0, 50.0, 90.0};


    ////////////////
    //Output list//
    ///////////////
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");

    fNumTracks = new TH1F("fNumTracks", "Number of Tracks/evt", 1000, 0, 1000);
    fOutputList->Add(fNumTracks);
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);

    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
    fOutputList->Add(fVtxX);
    
    // Single Particle and inclusive charged hadron histos
    Int_t numbinsSingle[5] = {100, 64, 20, 10, 10};
    Double_t minvalSingle[5] = {0, -3.14159, -1.0, -10.0, 0.0};
    Double_t maxvalSingle[5] = {10.0, 3.14159, 1.0, 10.0, 100.0};

    fRealChargedDist = new THnSparseF("fRealChargedDist", "Real Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealChargedDist->Sumw2();
    fOutputList->Add(fRealChargedDist);

    fRealPrimaryChargedDist = new THnSparseF("fRealPrimaryChargedDist", "Real PrimaryCharged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealPrimaryChargedDist->Sumw2();
    fOutputList->Add(fRealPrimaryChargedDist);

    fRealSecondaryChargedDist = new THnSparseF("fRealSecondaryChargedDist", "Real SecondaryCharged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealSecondaryChargedDist->Sumw2();
    fOutputList->Add(fRealSecondaryChargedDist);

    fRealTriggerDist = new THnSparseF("fRealTriggerDist", "Real Trigger Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealTriggerDist->Sumw2();
    fOutputList->Add(fRealTriggerDist);

    fRealKDist = new THnSparseF("fRealKDist", "Real Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealKDist->Sumw2();
    fOutputList->Add(fRealKDist);

    fRealPiDist = new THnSparseF("fRealPiDist", "Real Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealPiDist->Sumw2();
    fOutputList->Add(fRealPiDist);

    fRealPiFromLambdaDist = new THnSparseF("fRealPiFromLambdaDist", "Real Pion (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealPiFromLambdaDist->Sumw2();
    fOutputList->Add(fRealPiFromLambdaDist);

    fRealpDist = new THnSparseF("fRealpDist", "Real proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealpDist->Sumw2();
    fOutputList->Add(fRealpDist);

    fRealpFromLambdaDist = new THnSparseF("fRealpFromLambdaDist", "Real proton (from #Lambda) distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealpFromLambdaDist->Sumw2();
    fOutputList->Add(fRealpFromLambdaDist);

    fRealeDist = new THnSparseF("fRealeDist", "Real electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealeDist->Sumw2();
    fOutputList->Add(fRealeDist);

    fRealMuonDist = new THnSparseF("fRealMuonDist", "Real Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRealMuonDist->Sumw2();
    fOutputList->Add(fRealMuonDist);

    fRecoChargedDist = new THnSparseF("fRecoChargedDist", "Reco Charged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoChargedDist->Sumw2();
    fOutputList->Add(fRecoChargedDist);

    fRecoPrimaryChargedDist = new THnSparseF("fRecoPrimaryChargedDist", "Reco PrimaryCharged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoPrimaryChargedDist->Sumw2();
    fOutputList->Add(fRecoPrimaryChargedDist);

    fRecoSecondaryChargedDist = new THnSparseF("fRecoSecondaryChargedDist", "Reco SecondaryCharged Hadron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoSecondaryChargedDist->Sumw2();
    fOutputList->Add(fRecoSecondaryChargedDist);

    fRecoKDist = new THnSparseF("fRecoKDist", "Reco Kaon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoKDist->Sumw2();
    fOutputList->Add(fRecoKDist);

    fRecoPiDist = new THnSparseF("fRecoPiDist", "Reco Pion distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoPiDist->Sumw2();
    fOutputList->Add(fRecoPiDist);

    fRecopDist = new THnSparseF("fRecopDist", "Reco proton distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecopDist->Sumw2();
    fOutputList->Add(fRecopDist);

    fRecoeDist = new THnSparseF("fRecoeDist", "Reco electron distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoeDist->Sumw2();
    fOutputList->Add(fRecoeDist);

    fRecoMuonDist = new THnSparseF("fRecoMuonDist", "Reco Muon distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoMuonDist->Sumw2();
    fOutputList->Add(fRecoMuonDist);

    fRecoChargedTriggerDist = new THnSparseF("fRecoChargedTriggerDist", "Reco Charged Hadron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoChargedTriggerDist->Sumw2();
    fOutputList->Add(fRecoChargedTriggerDist);

    fRecoKTriggerDist = new THnSparseF("fRecoKTriggerDist", "Reco Kaon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoKTriggerDist->Sumw2();
    fOutputList->Add(fRecoKTriggerDist);
    
    fRecoPiTriggerDist = new THnSparseF("fRecoPiTriggerDist", "Reco Pion Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoPiTriggerDist->Sumw2();
    fOutputList->Add(fRecoPiTriggerDist);

    fRecopTriggerDist = new THnSparseF("fRecopTriggerDist", "Reco proton Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecopTriggerDist->Sumw2();
    fOutputList->Add(fRecopTriggerDist);

    fRecoeTriggerDist = new THnSparseF("fRecoeTriggerDist", "Reco electron Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoeTriggerDist->Sumw2();
    fOutputList->Add(fRecoeTriggerDist);

    fRecoMuonTriggerDist = new THnSparseF("fRecoMuonTriggerDist", "Reco Muon Trigger distribution;p_{T};#varphi;#eta;y;Z_{vtx};Multiplicity Percentile", 5, numbinsSingle, minvalSingle, maxvalSingle);
    fRecoMuonTriggerDist->Sumw2();
    fOutputList->Add(fRecoMuonTriggerDist);

    // Lambda ditributions used for eff calc
    
    Int_t numbins[6] = {100, 64, 20, 10, 140, 10};
    Double_t minval[6] = {0, -3.14159, -1., -10, 1.06, 0.0};
    Double_t maxval[6] = {10.0, 3.14159,  1.,  10, 1.2, 100.0};

    fRealTotalLambdaDist = new THnSparseF("fRealTotalLambdaDist", "Real (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRealTotalLambdaDist->Sumw2();
    fOutputList->Add(fRealTotalLambdaDist);

    fRealLambdaDist = new THnSparseF("fRealLambdaDist", "Real #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRealLambdaDist->Sumw2();
    fOutputList->Add(fRealLambdaDist);

    fRealAntiLambdaDist = new THnSparseF("fRealAntiLambdaDist", "Real #bar{#Lambda} distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{KK};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRealAntiLambdaDist->Sumw2();
    fOutputList->Add(fRealAntiLambdaDist);

    fRecoTotalV0LambdaDist = new THnSparseF("fRecoTotalV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoTotalV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoTotalV0LambdaDist);

    fRecoEtaV0LambdaDist = new THnSparseF("fRecoEtaV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaV0LambdaDist);

    fRecoEtaPtV0LambdaDist = new THnSparseF("fRecoEtaPtV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtV0LambdaDist);

    fRecoEtaPtRefitV0LambdaDist = new THnSparseF("fRecoEtaPtRefitV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitV0LambdaDist);

    fRecoEtaPtRefitRowsV0LambdaDist = new THnSparseF("fRecoEtaPtRefitRowsV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit + rows cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitRowsV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsV0LambdaDist);

    fRecoEtaPtRefitRowsRatioV0LambdaDist = new THnSparseF("fRecoEtaPtRefitRowsRatioV0LambdaDist", "Reco (#Lambda + #bar{#Lambda}) from V^{0} distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitRowsRatioV0LambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsRatioV0LambdaDist);

    fRecoEtaLambdaDist = new THnSparseF("fRecoEtaLambdaDist", "Reco (#Lambda + #bar{#Lambda})  distribution (eta cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaLambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaLambdaDist);

    fRecoEtaPtLambdaDist = new THnSparseF("fRecoEtaPtLambdaDist", "Reco (#Lambda + #bar{#Lambda})  distribution (eta + pt cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtLambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtLambdaDist);

    fRecoEtaPtRefitLambdaDist = new THnSparseF("fRecoEtaPtRefitLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitLambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitLambdaDist);

    fRecoEtaPtRefitRowsLambdaDist = new THnSparseF("fRecoEtaPtRefitRowsLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitRowsLambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsLambdaDist);

    fRecoEtaPtRefitRowsRatioLambdaDist = new THnSparseF("fRecoEtaPtRefitRowsRatioLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoEtaPtRefitRowsRatioLambdaDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsRatioLambdaDist);

    fRecoTotalLambdaDist = new THnSparseF("fRecoTotalLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoTotalLambdaDist->Sumw2();
    fOutputList->Add(fRecoTotalLambdaDist);

    fRecoLambdaWithAODPionDist = new THnSparseF("fRecoLambdaWithAODPionDist", "Reco (#Lambda + #bar{#Lambda}) distribution (AOD pion, real proton);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoLambdaWithAODPionDist->Sumw2();
    fOutputList->Add(fRecoLambdaWithAODPionDist);

    fRecoLambdaWithAODProtonDist = new THnSparseF("fRecoLambdaWithAODProtonDist", "Reco (#Lambda + #bar{#Lambda}) distribution (AOD proton, real pion);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoLambdaWithAODProtonDist->Sumw2();
    fOutputList->Add(fRecoLambdaWithAODProtonDist);

    fRecoTotalLambdaDist = new THnSparseF("fRecoTotalLambdaDist", "Reco (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoTotalLambdaDist->Sumw2();
    fOutputList->Add(fRecoTotalLambdaDist);

    fRecoLambdaDist = new THnSparseF("fRecoLambdaDist", "Reco #Lambda distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoLambdaDist->Sumw2();
    fOutputList->Add(fRecoLambdaDist);

    fRecoAntiLambdaDist = new THnSparseF("fRecoAntiLambdaDist", "Reco #Lambda (anti) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 6, numbins, minval, maxval);
    fRecoAntiLambdaDist->Sumw2();
    fOutputList->Add(fRecoAntiLambdaDist);

    fReactionPlane = new TH1D("fReactionPlane", "Reaction Plane Angle; Angle", 64, -3.14159, 3.14159);
    fOutputList->Add(fReactionPlane);

    // Same lambda distributions now with info on daughter dca (last two are proton dca, pion dca)
    Int_t numbinsDCA[8] = {100, 64, 20, 10, 140, 10, 100, 100};
    Double_t minvalDCA[8] = {0, -3.14159, -1, -10, 1.06, 0.0, -10, -10};
    Double_t maxvalDCA[8] = {10.0, 3.14159,  1,  10, 1.2, 100.0, 10, 10};

    fRecoEtaPtRefitRowsRatioLambdaDCADist = new THnSparseF("fRecoEtaPtRefitRowsRatioLambdaDCADist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 8, numbinsDCA, minvalDCA, maxvalDCA);
    fRecoEtaPtRefitRowsRatioLambdaDCADist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsRatioLambdaDCADist);

    // Same lambda distributions now with info on daughter filter maps (last two are proton filter map, pion filter map)
    Int_t numbinsFilter[8] = {100, 64, 20, 10, 140, 10, 9, 9};
    Double_t minvalFilter[8] = {0, -3.14159, -1, -10, 1.06, 0.0, 0, 0};
    Double_t maxvalFilter[8] = {10.0, 3.14159,  1,  10, 1.2, 100.0, 9, 9};

    fRecoTotalLambdaFilterDist = new THnSparseF("fRecoTotalLambdaFilterDist", "Reco (#Lambda + #bar{#Lambda}) distribution;p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 8, numbinsFilter, minvalFilter, maxvalFilter);
    fRecoTotalLambdaFilterDist->Sumw2();
    fOutputList->Add(fRecoTotalLambdaFilterDist);

    fRecoEtaPtRefitRowsRatioLambdaFilterDist = new THnSparseF("fRecoEtaPtRefitRowsRatioLambdaFilterDist", "Reco (#Lambda + #bar{#Lambda}) distribution (eta + pt + refit + rows + ratio cuts on daughters);p_{T};#varphi;#eta;y;Z_{vtx};m_{p#pi};Multiplicity Pctl.", 8, numbinsFilter, minvalFilter, maxvalFilter);
    fRecoEtaPtRefitRowsRatioLambdaFilterDist->Sumw2();
    fOutputList->Add(fRecoEtaPtRefitRowsRatioLambdaFilterDist);

    fRealLambdasPerEvent = new TH1F("fRealLambdasPerEvent", "Real Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRealLambdasPerEvent);

    fRecoLambdasPerEvent = new TH1F("fRecoLambdasPerEvent", "Reco Lambdas per Event", 10, 0, 10);
    fOutputList->Add(fRecoLambdasPerEvent);

    fPxDifference = new TH1D("fPxDifference", "Real p_x - reco p_x", 1000, -1, 1);
    fOutputList->Add(fPxDifference);

    fPxDifferenceFB = new TH1D("fPxDifferenceFB", "Real p_x - reco p_x (filter bit)", 1000, -1, 1);
    fOutputList->Add(fPxDifferenceFB);

    fPyDifference = new TH1D("fPyDifference", "Real p_y - reco p_y", 1000, -1, 1);
    fOutputList->Add(fPyDifference);

    fPzDifference = new TH1D("fPzDifference", "Real p_z - reco p_z", 1000, -1, 1);
    fOutputList->Add(fPzDifference);

    fPtDifference = new TH1D("fPtDifference", "Real p_t - reco p_t", 1000, -1, 1);
    fOutputList->Add(fPtDifference);


    // Invariant mass comparison between resonance and V0
    fInvMassLambdaResonance = new TH1D("fInvMassLambdaResonance", "Invariant Mass of #Lambda Resonance", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassLambdaResonance);

    fInvMassAntiLambdaResonance = new TH1D("fInvMassAntiLambdaResonance", "Invariant Mass of #bar{#Lambda} Resonance", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassAntiLambdaResonance);

    fInvMassLambdaV0 = new TH1D("fInvMassLambdaV0", "Invariant Mass of #Lambda V0", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassLambdaV0);

    fInvMassAntiLambdaV0 = new TH1D("fInvMassAntiLambdaV0", "Invariant Mass of #bar{#Lambda} V0", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassAntiLambdaV0);

    fInvMassLambdaDifference = new TH1D("fInvMassLambdaDifference", "Invariant MassLambda Difference", 1000, -0.2, 0.2);
    fOutputList->Add(fInvMassLambdaDifference);

    fInvMassAntiLambdaDifference = new TH1D("fInvMassAntiLambdaDifference", "Invariant MassAntiLambda Difference", 1000, -0.2, 0.2);
    fOutputList->Add(fInvMassAntiLambdaDifference);

    fInvMassLambdaReal = new TH1D("fInvMassLambdaReal", "Invariant Mass of Real #Lambda", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassLambdaReal);

    fInvMassAntiLambdaReal = new TH1D("fInvMassAntiLambdaReal", "Invariant Mass of Real #bar{#Lambda}", 1000, 1.06, 1.16);
    fOutputList->Add(fInvMassAntiLambdaReal);

    fRealPrimaryLambdaPtDist = new TH1D("fRealPrimaryLambdaPtDist", "Real Primary #Lambda p_{T} Distribution", 100, 0, 10);
    fOutputList->Add(fRealPrimaryLambdaPtDist);

    fRealSecondaryLambdaPtDist = new TH1D("fRealSecondaryLambdaPtDist", "Real Secondary #Lambda p_{T} Distribution", 100, 0, 10);
    fOutputList->Add(fRealSecondaryLambdaPtDist);

    fPhiDifferenceResV0 = new TH1D("fPhiDifferenceResV0", "#phi Difference between Resonance and V0", 128, -3.14159, 3.14159);
    fOutputList->Add(fPhiDifferenceResV0);

    fPhiDifferenceResReal = new TH1D("fPhiDifferenceResReal", "#phi Difference between Resonance and Real", 128, -3.14159, 3.14159);
    fOutputList->Add(fPhiDifferenceResReal);

    fPhiDifferenceV0Real = new TH1D("fPhiDifferenceV0Real", "#phi Difference between V0 and Real", 128, -3.14159, 3.14159);
    fOutputList->Add(fPhiDifferenceV0Real);

    fPhiReal = new TH1D("fPhiReal", "#phi of Real", 64, -3.14159, 3.14159);
    fOutputList->Add(fPhiReal);

    fPhiRes = new TH1D("fPhiRes", "#phi of Resonance", 64, -3.14159, 3.14159);
    fOutputList->Add(fPhiRes);

    fPhiV0 = new TH1D("fPhiV0", "#phi of V0", 64, -3.14159, 3.14159);
    fOutputList->Add(fPhiV0);

    PostData(1,fOutputList);

}

unsigned int AliAnalysisTaskLambdaHadronEfficiency::PassDaughterCuts(AliAODTrack *track){

    // Reject the negative ID tracks (TPC constrained to PV, etc.)
    if(track->GetID() < 0) return false;

    unsigned int passLevel = 0;

    // Bit set if track passes eta cut
    if(TMath::Abs(track->Eta()) < DAUGHTER_ETA_CUT) passLevel |= ETA_BIT;

    // Bit set if track passes pt cut

    if(track->Pt() >= DAUGHTER_MIN_PT) passLevel |= PT_BIT;

    // Bit set if track has TPC refit flag enabled
    if(track->IsOn(AliAODTrack::kTPCrefit)) passLevel |= TPC_REFIT_BIT;

    // Bit set if track has more than min number of crossed rows in TPC (set in constructor, usually 70)
    if(track->GetTPCCrossedRows() > MIN_CROSSED_ROWS_TPC) passLevel |= CROSSED_ROWS_BIT;

    // Bit set if the crossed rows/findable clusters is greater than a certain value (set in constructor, usually 0.8)
    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    if(ratio > MIN_ROW_CLUSTER_RATIO) passLevel |= ROW_CLUSTER_RATIO_BIT;

    return passLevel;

}

Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassAssociatedCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) < 0.8);
    pass = pass && (track->Pt() > 0.15);

    pass = pass && track->TestFilterMask(ASSOC_TRK_BIT);

    return pass;
}

Bool_t AliAnalysisTaskLambdaHadronEfficiency::PassTriggerCuts(AliAODTrack *track){

    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) < 0.8);
    pass = pass && (track->Pt() > 0.15);

    pass = pass && track->TestBit(TRIG_TRK_BIT);

    return pass;
}

void AliAnalysisTaskLambdaHadronEfficiency::UserExec(Option_t *){

    //masks for the different cut configurations
    unsigned int maskEta = ETA_BIT;
    unsigned int maskEtaPt = ETA_BIT + PT_BIT;
    unsigned int maskEtaPtRefit = ETA_BIT + PT_BIT + TPC_REFIT_BIT;
    unsigned int maskEtaPtRefitRows = ETA_BIT + PT_BIT + TPC_REFIT_BIT + CROSSED_ROWS_BIT;
    unsigned int maskEtaPtRefitRowsRatio = ETA_BIT + PT_BIT + TPC_REFIT_BIT + CROSSED_ROWS_BIT + ROW_CLUSTER_RATIO_BIT;

    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
     
    

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) return;

    ///////////////////
    //PID initialised//
    //////////////////

   fpidResponse = fInputHandler->GetPIDResponse();

    if(!fpidResponse){
        AliFatal("No PID response loaded!");
    }
    
    ////////////////
    //Event vertex//
    ///////////////
    Int_t ntracks = -999;
    ntracks = fVevent->GetNumberOfTracks();
    
    /////////////////
    //trigger check//
    /////////////////
    fVevent->GetFiredTriggerClasses();
    
    Int_t trigger = -1;
    //Multiplicity stuff
    Double_t multPercentile = 10.0;
    if (fAOD){
        AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
        Double_t multiplicity = header->GetRefMultiplicity();

        fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
        if(fMultSelection){
            multPercentile = fMultSelection->GetMultiplicityPercentile(CENT_ESTIMATOR.Data());
        }else{
            return;
        }
        if(multPercentile < MULT_LOW || multPercentile > MULT_HIGH) return;
    }
    
    fNevents->Fill(0); //all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<3)return;
    fNevents->Fill(1); //events with 3 tracks

    Zvertex = pVtx->GetZ();
    Yvertex = pVtx->GetY();
    Xvertex = pVtx->GetX();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);

    fNumTracks->Fill(ntracks);

    ////////////////////
    //event selection//
    ///////////////////
    if(fabs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut

    Double_t distPoint[6] = {0., 0., 0., 0., 0., 0.};
    Double_t singledistPoint[5] = {0., 0., 0., 0., 0.};
    //Loop over all particles in stack to get real phi, looking for phi->KK
    //AliMCEvent *fMCEvent = dynamic_cast<AliMCEvent*>(InputEvent());
    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
        AliError("Array of MC particles not found");
        return;
    }

    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!fMCHeader) {
        AliError("Could not find MC Header in AOD");
        return;
    }

    Double_t RA = fMCHeader->GetReactionPlaneAngle();
    if(RA > TMath::Pi()){
        RA -= 2.0*TMath::Pi();
    }else if(RA < -1.0*TMath::Pi()){
        RA += 2.0*TMath::Pi();
    }
    fReactionPlane->Fill(RA);

    UInt_t negPassCuts = 0;
    UInt_t posPassCuts = 0;

    Int_t negparPDG = 0;
    Int_t posparPDG = 0;
    Float_t recoPx = 0.0;
    Float_t recoPy = 0.0;
    Float_t recoPz = 0.0;
    Float_t recoP = 0.0;
    Float_t recoE = 0.0;
    Float_t recoM = 0.0;
    Float_t recoEta = 0.0;
    Float_t recoY = 0.0;
    Float_t recoPt = 0.0;
    Float_t recoPhi = 0.0;

    Int_t motherIndex = 0;
    Int_t motherPDG = 0;
    Int_t pdgcode = 0;

    int numRealLambdas = 0;
    int numRecoLambdas = 0;

    // quick test loop over V0's to determine filter mask range for daughter tracks
    int numV0s = fAOD->GetNumberOfV0s();
    for(int iv0 = 0; iv0 < numV0s; iv0++) {

        AliAODv0 *vZero = fAOD->GetV0(iv0);
        if(!vZero) continue;
        if(vZero->GetOnFlyStatus()) continue;

        AliAODTrack *ptrack=(AliAODTrack *)vZero->GetDaughter(0);
        AliAODTrack *ntrack=(AliAODTrack *)vZero->GetDaughter(1);
        
        int plabel = ptrack->GetLabel();
        int nlabel = ntrack->GetLabel();

        if(plabel < 0 || nlabel < 0) continue;

        AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(plabel);
        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(nlabel);

        int posPDG = mcpospart->GetPdgCode();
        int negPDG = mcnegpart->GetPdgCode();

        if(!((posPDG == 2212 && negPDG == -211) || (posPDG == 211 && negPDG == -2212))) continue;

        int posmomlabel = mcpospart->GetMother(); 
        int negmomlabel = mcnegpart->GetMother(); 

        if(posmomlabel < 0 || negmomlabel < 0) continue;

        if(posmomlabel != negmomlabel) {
            continue ;
        }

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(posmomlabel);

        int motherPDG = mcmother->GetPdgCode();

        if(!(TMath::Abs(motherPDG) == 3122)) continue;


        // THIS IS NOT BEING USED, NOW USE V0 METHODS TO CALCULATE KINEMATIC VARIABLES

         double recoPx = ntrack->Px() + ptrack->Px();
         double recoPy = ntrack->Py() + ptrack->Py();
         double recoPz = ntrack->Pz() + ptrack->Pz();

         double recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
         double recoE;
         if(motherPDG < 0) { 
             recoE = TMath::Sqrt(ptrack->Px()*ptrack->Px() + ptrack->Py()*ptrack->Py() + ptrack->Pz()*ptrack->Pz() + 0.13957*0.13957) + TMath::Sqrt(ntrack->Px()*ntrack->Px() + ntrack->Py()*ntrack->Py() + ntrack->Pz()*ntrack->Pz() + 0.9383*0.9383);
         } 
         else { 
             recoE = TMath::Sqrt(ntrack->Px()*ntrack->Px() + ntrack->Py()*ntrack->Py() + ntrack->Pz()*ntrack->Pz() + 0.13957*0.13957) + TMath::Sqrt(ptrack->Px()*ptrack->Px() + ptrack->Py()*ptrack->Py() + ptrack->Pz()*ptrack->Pz() + 0.9383*0.9383);
         } 
         double recoM = TMath::Sqrt(recoE*recoE - recoP*recoP); 
         double recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy); 
         double recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz)); 
         double recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz)); 
         double recoPhi = TMath::ATan2(recoPy, recoPx);
         if(recoPhi < -TMath::Pi()){ 
             recoPhi += 2.0*TMath::Pi(); 
         }else if(recoPhi > TMath::Pi()){ 
             recoPhi -= 2.0*TMath::Pi(); 
         } 
        double distPoint[6];

        distPoint[0] = vZero->Pt();
        distPoint[1] = vZero->Phi();
        distPoint[2] = vZero->Eta();
        distPoint[3] = Zvertex;
        if(motherPDG < 0) {
            distPoint[4] = vZero->MassAntiLambda();
        }
        else {
            distPoint[4] = vZero->MassLambda();
        }
        distPoint[5] = multPercentile;



        fRecoTotalV0LambdaDist->Fill(distPoint);

        unsigned int posPassCuts = PassDaughterCuts(ptrack);
        unsigned int negPassCuts = PassDaughterCuts(ntrack);

        if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
            fRecoEtaV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
            fRecoEtaPtV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
            fRecoEtaPtRefitV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
            fRecoEtaPtRefitRowsV0LambdaDist->Fill(distPoint);
        }

        if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
            fRecoEtaPtRefitRowsRatioV0LambdaDist->Fill(distPoint);

            std::cout << vZero->RadiusV0() << std::endl;

            float newPhiV0;
            if(vZero->Phi() < -TMath::Pi()){
                newPhiV0 = vZero->Phi() + 2.0*TMath::Pi();
            }
            else if(vZero->Phi() > TMath::Pi()){
                newPhiV0 = vZero->Phi() - 2.0*TMath::Pi();
            }
            else{
                newPhiV0 = vZero->Phi();
            }

            float newPhiReal;
            if(mcmother->Phi() < -TMath::Pi()){
                newPhiReal = mcmother->Phi() + 2.0*TMath::Pi();
            }
            else if(mcmother->Phi() > TMath::Pi()){
                newPhiReal = mcmother->Phi() - 2.0*TMath::Pi();
            }
            else{
                newPhiReal = mcmother->Phi();
            }

            fPhiDifferenceResV0->Fill(recoPhi - newPhiV0);
            fPhiDifferenceResReal->Fill(recoPhi - newPhiReal);
            fPhiDifferenceV0Real->Fill(newPhiV0 - newPhiReal);
            fPhiV0->Fill(newPhiV0);
            fPhiRes->Fill(recoPhi);
            fPhiReal->Fill(newPhiReal);

            if (motherPDG < 0) {
                double massAntiLambdaReal = mcmother->M();
                double massAntiLambdaV0 = vZero->MassAntiLambda();
                double massAntiLambdaResonance = recoM;
                fInvMassAntiLambdaDifference->Fill(massAntiLambdaV0 - massAntiLambdaResonance);
                fInvMassAntiLambdaResonance->Fill(massAntiLambdaResonance);
                fInvMassAntiLambdaV0->Fill(massAntiLambdaV0);
                fInvMassAntiLambdaReal->Fill(massAntiLambdaReal);
            }
            else {
                double massLambdaReal = mcmother->M();
                double massLambdaV0 = vZero->MassLambda();
                double massLambdaResonance = recoM;
                fInvMassLambdaDifference->Fill(massLambdaV0 - massLambdaResonance);
                fInvMassLambdaResonance->Fill(massLambdaResonance);
                fInvMassLambdaV0->Fill(massLambdaV0);
                fInvMassLambdaReal->Fill(massLambdaReal);
            }
        }


    }
    
    std::vector<AliAODMCParticle*> lambdas_with_aod_pion;
    std::vector<AliAODMCParticle*> lambdas_with_aod_proton;

    //first loop over tracks to get pi minuses, find lambda daughters and fill Reco Dist
    // fun padding for git activity smile
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();

        if(tracklabel < 0) {
            continue;
        }

        //Get all single particle distributions
        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        pdgcode = mcpart->GetPdgCode();
        Bool_t associatedPass = PassAssociatedCuts(aodnegtrack);
        Bool_t triggerPass = PassTriggerCuts(aodnegtrack);
        singledistPoint[0] = aodnegtrack->Pt();
        singledistPoint[1] = aodnegtrack->Phi();
        if(singledistPoint[1] > TMath::Pi()){
            singledistPoint[1] -= 2.0*TMath::Pi();
        }else if(singledistPoint[1] < -1.0*TMath::Pi()){
            singledistPoint[1] += 2.0*TMath::Pi();
        }
        singledistPoint[2] = aodnegtrack->Eta();
        singledistPoint[3] = Zvertex;
        singledistPoint[4] = multPercentile;


        if(associatedPass){ 
            if(TMath::Abs(pdgcode) == 211){
                fRecoPiDist->Fill(singledistPoint);
                fRecoChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 321){
                fRecoKDist->Fill(singledistPoint);
                fRecoChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRecopDist->Fill(singledistPoint);
                fRecoChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                fRecoeDist->Fill(singledistPoint);
                fRecoChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRecoMuonDist->Fill(singledistPoint);
                fRecoChargedDist->Fill(singledistPoint);
            }
        }

        unsigned int daughterCuts = PassDaughterCuts(aodnegtrack);
        bool daughterPass = ((daughterCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio);

        if(mcpart->IsPhysicalPrimary()){
            if(daughterPass){ 
                if(TMath::Abs(pdgcode) == 211){
                    fRecoPrimaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 321){
                    fRecoPrimaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecoPrimaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoPrimaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoPrimaryChargedDist->Fill(singledistPoint);
                }
            }
        }
        else {
            if(daughterPass){ 
                if(TMath::Abs(pdgcode) == 211){
                    fRecoSecondaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 321){
                    fRecoSecondaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 2212){
                    fRecoSecondaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 11){
                    fRecoSecondaryChargedDist->Fill(singledistPoint);
                }else if(TMath::Abs(pdgcode) == 13){
                    fRecoSecondaryChargedDist->Fill(singledistPoint);
                }
            }
        }

        if(triggerPass){
            if(TMath::Abs(pdgcode) == 321){
                fRecoKTriggerDist->Fill(singledistPoint);
                fRecoChargedTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 211){
                fRecoPiTriggerDist->Fill(singledistPoint);
                fRecoChargedTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRecopTriggerDist->Fill(singledistPoint);
                fRecoChargedTriggerDist->Fill(singledistPoint);
        }else if(TMath::Abs(pdgcode) == 11){
                fRecoeTriggerDist->Fill(singledistPoint);
                fRecoChargedTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRecoMuonTriggerDist->Fill(singledistPoint);
                fRecoChargedTriggerDist->Fill(singledistPoint);
            }

        }

        negPassCuts = PassDaughterCuts(aodnegtrack);
        // check to see if daughter passes all of our AOD cuts
        if((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) {
            if(TMath::Abs(pdgcode) == 211){
                int motherlabel = mcpart->GetMother();
                if(motherlabel >= 0) {
                    auto mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
                    int motherPDG = mother->GetPdgCode();
                    if(TMath::Abs(motherPDG) == 3122) {
                        lambdas_with_aod_pion.push_back(mother);
                    }
                }
            }
            if(TMath::Abs(pdgcode) == 2212){
                int motherlabel = mcpart->GetMother();
                if(motherlabel >= 0) {
                    auto mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
                    int motherPDG = mother->GetPdgCode();
                    if(TMath::Abs(motherPDG) == 3122) {
                        lambdas_with_aod_proton.push_back(mother);
                    }
                }
            }
        }

        if((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) {
            float px_dif = mcpart->Px() - aodnegtrack->Px();
            float py_dif = mcpart->Py() - aodnegtrack->Py();
            float pz_dif = mcpart->Pz() - aodnegtrack->Pz();
            float pt_dif = mcpart->Pt() - aodnegtrack->Pt();

            fPxDifference->Fill(px_dif);
            fPyDifference->Fill(py_dif);
            fPzDifference->Fill(pz_dif);
            fPtDifference->Fill(pt_dif);
        }

        if(aodnegtrack->TestFilterBit(AliAODTrack::kTrkGlobal)) {
            float px_dif_fb = mcpart->Px() - aodnegtrack->Px();
            fPxDifferenceFB->Fill(px_dif_fb);
        }

        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        negparPDG = mcnegpart->GetPdgCode();

        if(negparPDG != -211) continue;

        motherIndex = mcnegpart->GetMother();
        if(motherIndex < 0) continue;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
        motherPDG = mcmother->GetPdgCode();
        if(motherPDG != 3122) continue;

        if((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) {
            lambdas_with_aod_pion.push_back(mcmother);
        }


        for(int jtrack = 0; jtrack < ntracks; jtrack++){
            AliVParticle *vpospart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(jtrack));
            AliVTrack *postrack = dynamic_cast<AliVTrack*>(vpospart);
            AliAODTrack *aodpostrack = dynamic_cast<AliAODTrack*>(vpospart);

            Double_t posTPCnSigma = -999;
            Double_t posTOFnSigma = -999;
            if(postrack->Pt() > 0.15){
                posTPCnSigma = fpidResponse->NumberOfSigmasTPC(postrack, AliPID::kProton);
                posTOFnSigma = fpidResponse->NumberOfSigmasTOF(postrack, AliPID::kProton);
            }

            posPassCuts = PassDaughterCuts(aodpostrack);

            Int_t postracklabel = aodpostrack->GetLabel();
            if(postracklabel < 0) continue;

            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(postracklabel);
            posparPDG = mcpospart->GetPdgCode();
            if(posparPDG != 2212) continue;

            if(mcpospart->GetMother() == motherIndex){

                recoPx = aodnegtrack->Px() + aodpostrack->Px();
                recoPy = aodnegtrack->Py() + aodpostrack->Py();
                recoPz = aodnegtrack->Pz() + aodpostrack->Pz();
                
                recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
                recoE = TMath::Sqrt(aodnegtrack->Px()*aodnegtrack->Px() + aodnegtrack->Py()*aodnegtrack->Py() + aodnegtrack->Pz()*aodnegtrack->Pz() + 0.13957*0.13957) + TMath::Sqrt(aodpostrack->Px()*aodpostrack->Px() + aodpostrack->Py()*aodpostrack->Py() + aodpostrack->Pz()*aodpostrack->Pz() + 0.9383*0.9383);
                recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
                recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
                recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
                recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
                recoPhi = TMath::ATan2(recoPy, recoPx);
                if(recoPhi < -TMath::Pi()){
                    recoPhi += 2.0*TMath::Pi();
                }else if(recoPhi > TMath::Pi()){
                    recoPhi -= 2.0*TMath::Pi();
                }

                distPoint[0] = recoPt;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;

                fRecoLambdaDist->Fill(distPoint);
                fRecoTotalLambdaDist->Fill(distPoint);

                double filter_distPoint[8];
                filter_distPoint[0] = recoPt;
                filter_distPoint[1] = recoPhi;
                filter_distPoint[2] = recoEta;
                filter_distPoint[3] = Zvertex;
                filter_distPoint[4] = recoM;
                filter_distPoint[5] = multPercentile;
                filter_distPoint[6] = filterMap_map[aodpostrack->GetFilterMap()];
                filter_distPoint[7] = filterMap_map[aodnegtrack->GetFilterMap()];

                double proton_dz[2];
                double proton_covar[3];
                double pion_dz[2];
                double pion_covar[3];
                bool is_protonDCA = aodpostrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., proton_dz, proton_covar);
                bool is_pionDCA = aodnegtrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., pion_dz, pion_covar);
                
                fRecoTotalLambdaFilterDist->Fill(filter_distPoint);

                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fRecoEtaLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fRecoEtaPtLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
                    fRecoEtaPtRefitLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
                    fRecoEtaPtRefitRowsLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
                    fRecoEtaPtRefitRowsRatioLambdaDist->Fill(distPoint);
                    fRecoEtaPtRefitRowsRatioLambdaFilterDist->Fill(filter_distPoint);
                    if(is_pionDCA && is_protonDCA) {
                        double dca_distPoint[8];
                        dca_distPoint[0] = recoPt;
                        dca_distPoint[1] = recoPhi;
                        dca_distPoint[2] = recoEta;
                        dca_distPoint[3] = Zvertex;
                        dca_distPoint[4] = recoM;
                        dca_distPoint[5] = multPercentile;
                        dca_distPoint[6] = proton_dz[0];
                        dca_distPoint[7] = pion_dz[0];
                        fRecoEtaPtRefitRowsRatioLambdaDCADist->Fill(dca_distPoint);
                    }
                    numRecoLambdas += 1;
                }
            }
        }
    }

    //second loop over tracks to get antiprotons, find lambda daughters and fill Reco Dist
    for(int itrack = 0; itrack < ntracks; itrack++){
        AliVParticle *vnegpart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(itrack));
        AliVTrack *negtrack = dynamic_cast<AliVTrack*>(vnegpart);
        AliAODTrack *aodnegtrack = dynamic_cast<AliAODTrack*>(vnegpart);

        Int_t tracklabel = aodnegtrack->GetLabel();
        if(tracklabel < 0) continue;

        //Get pions that came from lambda for lambda reco
        negPassCuts = PassDaughterCuts(aodnegtrack);

        AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(tracklabel);
        negparPDG = mcnegpart->GetPdgCode();
        if(negparPDG != -2212) continue;

        motherIndex = mcnegpart->GetMother();
        if(motherIndex < 0) continue;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
        motherPDG = mcmother->GetPdgCode();
        if(motherPDG != -3122) continue;

        for(int jtrack = 0; jtrack < ntracks; jtrack++){
            AliVParticle *vpospart = dynamic_cast<AliVParticle*>(fVevent->GetTrack(jtrack));
            AliVTrack *postrack = dynamic_cast<AliVTrack*>(vpospart);
            AliAODTrack *aodpostrack = dynamic_cast<AliAODTrack*>(vpospart);

            posPassCuts = PassDaughterCuts(aodpostrack);

            Int_t postracklabel = aodpostrack->GetLabel();
            if(postracklabel < 0) continue;

            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(postracklabel);
            posparPDG = mcpospart->GetPdgCode();
            if(posparPDG != 211) continue;

            if(mcpospart->GetMother() == motherIndex){

                recoPx = aodnegtrack->Px() + aodpostrack->Px();
                recoPy = aodnegtrack->Py() + aodpostrack->Py();
                recoPz = aodnegtrack->Pz() + aodpostrack->Pz();
                
                recoP = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy + recoPz*recoPz);
                recoE = TMath::Sqrt(aodnegtrack->Px()*aodnegtrack->Px() + aodnegtrack->Py()*aodnegtrack->Py() + aodnegtrack->Pz()*aodnegtrack->Pz() + 0.9383*0.9383) + TMath::Sqrt(aodpostrack->Px()*aodpostrack->Px() + aodpostrack->Py()*aodpostrack->Py() + aodpostrack->Pz()*aodpostrack->Pz() + 0.13957*0.13957);
                recoM = TMath::Sqrt(recoE*recoE - recoP*recoP);
                recoPt = TMath::Sqrt(recoPx*recoPx + recoPy*recoPy);
                recoEta = 0.5*TMath::Log((recoP + recoPz)/(recoP -  recoPz));
                recoY = 0.5*TMath::Log((recoE + recoPz)/(recoE - recoPz));
                recoPhi = TMath::ATan2(recoPy, recoPx);
                if(recoPhi < -TMath::Pi()){
                    recoPhi += 2.0*TMath::Pi();
                }else if(recoPhi > TMath::Pi()){
                    recoPhi -= 2.0*TMath::Pi();
                }
                distPoint[0] = recoPt;
                distPoint[1] = recoPhi;
                distPoint[2] = recoEta;
                distPoint[3] = Zvertex;
                distPoint[4] = recoM;
                distPoint[5] = multPercentile;

                fRecoAntiLambdaDist->Fill(distPoint);
                fRecoTotalLambdaDist->Fill(distPoint);

                double filter_distPoint[8];
                filter_distPoint[0] = recoPt;
                filter_distPoint[1] = recoPhi;
                filter_distPoint[2] = recoEta;
                filter_distPoint[3] = Zvertex;
                filter_distPoint[4] = recoM;
                filter_distPoint[5] = multPercentile;
                filter_distPoint[6] = filterMap_map[aodnegtrack->GetFilterMap()];
                filter_distPoint[7] = filterMap_map[aodpostrack->GetFilterMap()];

                double proton_dz[2];
                double proton_covar[3];
                double pion_dz[2];
                double pion_covar[3];
                bool is_protonDCA = aodnegtrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., proton_dz, proton_covar);
                bool is_pionDCA = aodpostrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., pion_dz, pion_covar);

                fRecoTotalLambdaFilterDist->Fill(filter_distPoint);

                if(((negPassCuts & maskEta) == maskEta) && ((posPassCuts & maskEta)== maskEta)){
                    fRecoEtaLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPt) == maskEtaPt) && ((posPassCuts & maskEtaPt)== maskEtaPt)){
                    fRecoEtaPtLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefit) == maskEtaPtRefit) && ((posPassCuts & maskEtaPtRefit)== maskEtaPtRefit)){
                    fRecoEtaPtRefitLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRows) == maskEtaPtRefitRows) && ((posPassCuts & maskEtaPtRefitRows)== maskEtaPtRefitRows)){
                    fRecoEtaPtRefitRowsLambdaDist->Fill(distPoint);
                }

                if(((negPassCuts & maskEtaPtRefitRowsRatio) == maskEtaPtRefitRowsRatio) && ((posPassCuts & maskEtaPtRefitRowsRatio)== maskEtaPtRefitRowsRatio)){
                    fRecoEtaPtRefitRowsRatioLambdaDist->Fill(distPoint);
                    numRecoLambdas += 1;
                    if(is_pionDCA && is_protonDCA) {
                        double dca_distPoint[8];
                        dca_distPoint[0] = recoPt;
                        dca_distPoint[1] = recoPhi;
                        dca_distPoint[2] = recoEta;
                        dca_distPoint[3] = Zvertex;
                        dca_distPoint[4] = recoM;
                        dca_distPoint[5] = multPercentile;
                        dca_distPoint[6] = proton_dz[0];
                        dca_distPoint[7] = pion_dz[0];
                        fRecoEtaPtRefitRowsRatioLambdaDCADist->Fill(dca_distPoint);
                    }
                }
            }
        }
    }

    //loop over MC particles to get original lambda and fill Real dist (and real single particle dist)
    int numCharged = 0;

    for(Int_t imcpart=0; imcpart< fMCArray->GetEntries(); imcpart++){

        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imcpart);
        pdgcode = AODMCtrack->GetPdgCode();

        singledistPoint[0] = AODMCtrack->Pt();
        singledistPoint[1] = AODMCtrack->Phi();
        if(singledistPoint[1] > TMath::Pi()){
            singledistPoint[1] -= 2.0*TMath::Pi();
        }else if(singledistPoint[1] < -1.0*TMath::Pi()){
            singledistPoint[1] += 2.0*TMath::Pi();
        }
        float mcEta = AODMCtrack->Eta();
        singledistPoint[2] = mcEta;
        singledistPoint[3] = Zvertex;
        singledistPoint[4] = multPercentile;

        if(AODMCtrack->IsPhysicalPrimary() && TMath::Abs(mcEta) < 0.8){
            if(TMath::Abs(pdgcode) == 321){
                fRealKDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 211){
                fRealPiDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealpDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 11){
                fRealeDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }else if(TMath::Abs(pdgcode) == 13){
                fRealMuonDist->Fill(singledistPoint);
                fRealChargedDist->Fill(singledistPoint);
                numCharged += 1;
            }
        }

        if(AODMCtrack->IsPhysicalPrimary() && TMath::Abs(mcEta) < 0.8){
            if(TMath::Abs(pdgcode) == 321){
                fRealPrimaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 211){
                fRealPrimaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealPrimaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                fRealPrimaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRealPrimaryChargedDist->Fill(singledistPoint);
            }
        }

        if((!AODMCtrack->IsPhysicalPrimary()) && TMath::Abs(mcEta) < 0.8){
            if(TMath::Abs(pdgcode) == 321){
                fRealSecondaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 211){
                fRealSecondaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealSecondaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                fRealSecondaryChargedDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRealSecondaryChargedDist->Fill(singledistPoint);
            }
        }

        if(AODMCtrack->IsPhysicalPrimary() && TMath::Abs(mcEta) < 0.8 ){
            if(TMath::Abs(pdgcode) == 321){
                fRealTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 211){
                fRealTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 2212){
                fRealTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 11){
                fRealTriggerDist->Fill(singledistPoint);
            }else if(TMath::Abs(pdgcode) == 13){
                fRealTriggerDist->Fill(singledistPoint);
            }
        }

        // Get protons and pions that came from lambdas (would not be physical primaries) 
        if(TMath::Abs(pdgcode) == 211){
            motherIndex = AODMCtrack->GetMother();
            if(motherIndex >= 0) {
                AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                motherPDG = mcmother->GetPdgCode();
                if(TMath::Abs(motherPDG) == 3122) fRealPiFromLambdaDist->Fill(singledistPoint);
            } 
        } else if(TMath::Abs(pdgcode) == 2212){
            if(motherIndex >= 0) {
                AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(motherIndex);
                motherPDG = mcmother->GetPdgCode();
                if(TMath::Abs(motherPDG) == 3122) fRealpFromLambdaDist->Fill(singledistPoint);
            } 
        }

        //select lambdas
        if(TMath::Abs(pdgcode) != 3122) continue;
        Int_t indexFirstDaughter = 0, indexSecondDaughter = 0;
        indexFirstDaughter = AODMCtrack->GetDaughterFirst();
        indexSecondDaughter = AODMCtrack->GetDaughterLast();

        if(indexFirstDaughter < 0 || indexSecondDaughter < 0) continue;
        AliAODMCParticle* firstDaughter = (AliAODMCParticle*)fMCArray->At(indexFirstDaughter);
        AliAODMCParticle* secondDaughter = (AliAODMCParticle*)fMCArray->At(indexSecondDaughter);

        if(((TMath::Abs(firstDaughter->GetPdgCode()) == 211 && TMath::Abs(secondDaughter->GetPdgCode()) == 2212) ||
            (TMath::Abs(firstDaughter->GetPdgCode()) == 2212 && TMath::Abs(secondDaughter->GetPdgCode()) == 211)) && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) < 0){
            
            // check if mother is primary
            if(AODMCtrack->IsPhysicalPrimary()) {
                fRealPrimaryLambdaPtDist->Fill(AODMCtrack->Pt());
            }
            else {
                fRealSecondaryLambdaPtDist->Fill(AODMCtrack->Pt());
            }

            numRealLambdas += 1;
            distPoint[0] = AODMCtrack->Pt();
            distPoint[1] = AODMCtrack->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            float mcLambdaEta = AODMCtrack->Eta();
            distPoint[2] = mcLambdaEta;
            distPoint[3] = Zvertex;
            distPoint[4] = AODMCtrack->GetCalcMass();
            distPoint[5] = multPercentile;

            if(TMath::Abs(mcLambdaEta) < 0.8){
                fRealTotalLambdaDist->Fill(distPoint);

                if(pdgcode == 3122) {
                    fRealLambdaDist->Fill(distPoint); 
                }
                else {
                    fRealAntiLambdaDist->Fill(distPoint);
                }
            }
        } 
    }

    for(auto lambda : lambdas_with_aod_pion) {
        Int_t indexFirstDaughter = 0, indexSecondDaughter = 0;
        indexFirstDaughter = lambda->GetDaughterFirst();
        indexSecondDaughter = lambda->GetDaughterLast();

        if(indexFirstDaughter < 0 || indexSecondDaughter < 0) continue;
        AliAODMCParticle* firstDaughter = (AliAODMCParticle*)fMCArray->At(indexFirstDaughter);
        AliAODMCParticle* secondDaughter = (AliAODMCParticle*)fMCArray->At(indexSecondDaughter);

        if(((TMath::Abs(firstDaughter->GetPdgCode()) == 211 && TMath::Abs(secondDaughter->GetPdgCode()) == 2212) ||
            (TMath::Abs(firstDaughter->GetPdgCode()) == 2212 && TMath::Abs(secondDaughter->GetPdgCode()) == 211)) && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) < 0){
            
            distPoint[0] = lambda->Pt();
            distPoint[1] = lambda->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            float mcLambdaEta = lambda->Eta();
            distPoint[2] = mcLambdaEta;
            distPoint[3] = Zvertex;
            distPoint[4] = lambda->GetCalcMass();
            distPoint[5] = multPercentile;

            if(TMath::Abs(mcLambdaEta) < 0.8){
                fRecoLambdaWithAODPionDist->Fill(distPoint);
            }
        } 
    }

    for(auto lambda : lambdas_with_aod_proton) {
        Int_t indexFirstDaughter = 0, indexSecondDaughter = 0;
        indexFirstDaughter = lambda->GetDaughterFirst();
        indexSecondDaughter = lambda->GetDaughterLast();

        if(indexFirstDaughter < 0 || indexSecondDaughter < 0) continue;
        AliAODMCParticle* firstDaughter = (AliAODMCParticle*)fMCArray->At(indexFirstDaughter);
        AliAODMCParticle* secondDaughter = (AliAODMCParticle*)fMCArray->At(indexSecondDaughter);

        if(((TMath::Abs(firstDaughter->GetPdgCode()) == 211 && TMath::Abs(secondDaughter->GetPdgCode()) == 2212) ||
            (TMath::Abs(firstDaughter->GetPdgCode()) == 2212 && TMath::Abs(secondDaughter->GetPdgCode()) == 211)) && (firstDaughter->GetPdgCode())*(secondDaughter->GetPdgCode()) < 0){
            
            distPoint[0] = lambda->Pt();
            distPoint[1] = lambda->Phi();
            if(distPoint[1] > TMath::Pi()){
                distPoint[1] -= 2.0*TMath::Pi();
            }else if(distPoint[1] < -1.0*TMath::Pi()){
                distPoint[1] += 2.0*TMath::Pi();
            }
            float mcLambdaEta = lambda->Eta();
            distPoint[2] = mcLambdaEta;
            distPoint[3] = Zvertex;
            distPoint[4] = lambda->GetCalcMass();
            distPoint[5] = multPercentile;

            if(TMath::Abs(mcLambdaEta) < 0.8){
                fRecoLambdaWithAODProtonDist->Fill(distPoint);
            }
        } 
    }
    
    fRealLambdasPerEvent->Fill(numRealLambdas);
    fRecoLambdasPerEvent->Fill(numRecoLambdas);

    PostData(1, fOutputList);
}    
//________________________________________________________________________
void AliAnalysisTaskLambdaHadronEfficiency::Terminate(Option_t *) 
{
    // Draw result to the screen
    // Called once at the end of the query
    printf("terminating task... \n");
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
    
}
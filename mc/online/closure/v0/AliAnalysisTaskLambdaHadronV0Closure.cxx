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
#include "TH1D.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCFParticle.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliAODv0.h"

#include "AliAnalysisTaskLambdaHadronV0Closure.h"

ClassImp(AliAnalysisTaskLambdaHadronV0Closure);

AliAnalysisTaskLambdaHadronV0Closure::AliAnalysisTaskLambdaHadronV0Closure() :
    AliAnalysisTaskSE(),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_checkMC(0x0),
    fTriggerDist(0x0),
    fAssociatedHDist(0x0),
    fAssociatedHDist_checkMC(0x0),
    fTriggeredLambdaDist(0x0),
    fLambdaDist(0x0),
    fGuaranteedLambdaDist(0x0),
    fTriggeredLambdaDist_MC(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHGuaranteedLambdaEff(0x0),
    fDphiHLambdaEff_MCKin(0x0),
    fDphiHDaughterProton_MCKin(0x0),
    fDphiHDaughterPion_MCKin(0x0),
    fDphiHLambdaEff_MCKin_physicalPrimary(0x0),
    fDphiRecoHRealLambdaEff_MCKin_physicalPrimary(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_checkMC(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_MCKin(0x0),
    fDphiHDaughterProtonMixed_MCKin(0x0),
    fDphiHDaughterPionMixed_MCKin(0x0),
    fDphiHLambdaMixed_MCKin_physicalPrimary(0x0),
    fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary(0x0),
    fDphiHHMixed(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0),
    fMCCorPoolMgr(0x0),
    fLambdaDist_MC(0x0),
    fTriggerDist_MC(0x0),
    fAssociatedDist_MC(0x0),
    fDphiHLambda_MC(0x0),
    fDphiHLambda_MC_physicalPrimary(0x0),
    fDphiHH_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHLambdaMixed_MC_physicalPrimary(0x0),
    fDphiHHMixed_MC(0x0)
{
}

AliAnalysisTaskLambdaHadronV0Closure::AliAnalysisTaskLambdaHadronV0Closure(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_checkMC(0x0),
    fTriggerDist(0x0),
    fAssociatedHDist(0x0),
    fAssociatedHDist_checkMC(0x0),
    fTriggeredLambdaDist(0x0),
    fLambdaDist(0x0),
    fGuaranteedLambdaDist(0x0),
    fTriggeredLambdaDist_MC(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHGuaranteedLambdaEff(0x0),
    fDphiHLambdaEff_MCKin(0x0),
    fDphiHDaughterProton_MCKin(0x0),
    fDphiHDaughterPion_MCKin(0x0),
    fDphiHLambdaEff_MCKin_physicalPrimary(0x0),
    fDphiRecoHRealLambdaEff_MCKin_physicalPrimary(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_checkMC(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_MCKin(0x0),
    fDphiHDaughterProtonMixed_MCKin(0x0),
    fDphiHDaughterPionMixed_MCKin(0x0),
    fDphiHLambdaMixed_MCKin_physicalPrimary(0x0),
    fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary(0x0),
    fDphiHHMixed(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0),
    fMCCorPoolMgr(0x0),
    fLambdaDist_MC(0x0),
    fTriggerDist_MC(0x0),
    fAssociatedDist_MC(0x0),
    fDphiHLambda_MC(0x0),
    fDphiHLambda_MC_physicalPrimary(0x0),
    fDphiHH_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHLambdaMixed_MC_physicalPrimary(0x0),
    fDphiHHMixed_MC(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronV0Closure::~AliAnalysisTaskLambdaHadronV0Closure()
{
    if(fOutputList) delete fOutputList;
    if(fTriggerEff) delete fTriggerEff;
    if(fAssociatedEff) delete fAssociatedEff;
    if(fLambdaEff) delete fLambdaEff;
}

void AliAnalysisTaskLambdaHadronV0Closure::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {fMultLow, fMultHigh};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    fMCCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fMCCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {200, 16, 20, 10};
    double dist_mins[4] = {0.0, 0, -1, -10};
    double dist_maxes[4] = {20.0, 6.28, 1, 10};

    fTriggerDistEff = new THnSparseF("fTriggerDistEff", "Efficiency Corrected Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDistEff->Sumw2();
    fOutputList->Add(fTriggerDistEff);

    fTriggerDistEff_checkMC = new THnSparseF("fTriggerDistEff_checkMC", "Efficiency Corrected Trigger Hadron Distribution (is MC hadron, is MC physical primary)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDistEff_checkMC->Sumw2();
    fOutputList->Add(fTriggerDistEff_checkMC);

    fTriggerDist = new THnSparseF("fTriggerDist", "Non-Efficiency Corrected Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist->Sumw2();
    fOutputList->Add(fTriggerDist);

    fTriggerDist_MC = new THnSparseF("fTriggerDist_MC", "Trigger Hadron Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist_MC->Sumw2();
    fOutputList->Add(fTriggerDist_MC);

    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist->Sumw2();
    fOutputList->Add(fAssociatedHDist);

    fAssociatedHDist_checkMC = new THnSparseF("fAssociatedHDist_checkMC", "Associated Hadron Distribution (is MC hadron, is MC physical primary)", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist_checkMC->Sumw2();
    fOutputList->Add(fAssociatedHDist_checkMC);

    fAssociatedDist_MC = new THnSparseF("fAssociatedDist_MC", "Associated Hadron Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedDist_MC->Sumw2();
    fOutputList->Add(fAssociatedDist_MC);

    //Mother distribution axes are: Pt, Phi, Eta, Mass, Event multiplicity
    int mother_dist_bins[5] = {100, 16, 20, 100, 10};
    double mother_dist_mins[5] = {0, -3.14, -1, 1.06, 0};
    double mother_dist_maxes[5] = {15, 3.14, 1, 1.16, 100};

    fTriggeredLambdaDist = new THnSparseF("fTriggeredLambdaDist", "Lambda Distribution (with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDist->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist);

    fLambdaDist = new THnSparseF("fLambdaDist", "Lambda Distribution (reco with v0 finder)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fLambdaDist->Sumw2();
    fOutputList->Add(fLambdaDist);

    fGuaranteedLambdaDist = new THnSparseF("fGuaranteedLambdaDist", "Lambda Distribution (reco with v0 finder, checked MC to guarantee lambda)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fGuaranteedLambdaDist->Sumw2();
    fOutputList->Add(fGuaranteedLambdaDist);

    fLambdaDist_MC = new THnSparseF("fLambdaDist_MC", "Lambda Distribution (MC truth)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fLambdaDist_MC->Sumw2();
    fOutputList->Add(fLambdaDist_MC);

    fTriggeredLambdaDist_MC = new THnSparseF("fTriggeredLambdaDist_MC", "Lambda Distribution (MC truth, with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDist_MC->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist_MC);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hl_cor_bins[6] = {18, 16, 16, 20, 100, 10};
    double hl_cor_mins[6] = {1, 0, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hl_cor_maxes[6] = {10, 4, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10};

    fDphiHLambdaEff = new THnSparseF("fDphiHLambdaEff", "Efficiency-corrected Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff->Sumw2();
    fOutputList->Add(fDphiHLambdaEff);

    fDphiHGuaranteedLambdaEff = new THnSparseF("fDphiHGuaranteedLambdaEff", "Efficiency-corrected Hadron-Lambda Correlation Histogram (from guaranteed lambda)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHGuaranteedLambdaEff->Sumw2();
    fOutputList->Add(fDphiHGuaranteedLambdaEff);

    fDphiHLambdaEff_MCKin = new THnSparseF("fDphiHLambdaEff_MCKin", "Efficiency-corrected Hadron-Lambda Correlation Histogram (using MC kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff_MCKin->Sumw2();
    fOutputList->Add(fDphiHLambdaEff_MCKin);

    fDphiHDaughterProton_MCKin = new THnSparseF("fDphiHDaughterProton_MCKin", "Hadron-daughter proton Correlation Histogram (using MC kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterProton_MCKin->Sumw2();
    fOutputList->Add(fDphiHDaughterProton_MCKin);

    fDphiHDaughterPion_MCKin = new THnSparseF("fDphiHDaughterPion_MCKin", "Hadron-daughter pion Correlation Histogram (using MC kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterPion_MCKin->Sumw2();
    fOutputList->Add(fDphiHDaughterPion_MCKin);

    fDphiHLambdaEff_MCKin_physicalPrimary = new THnSparseF("fDphiHLambdaEff_MCKin_physicalPrimary", "Efficiency-corrected Hadron-Lambda Correlation Histogram (using MC kinematics on V0, trigger and lambda are physical primary)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff_MCKin_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiHLambdaEff_MCKin_physicalPrimary);

    fDphiRecoHRealLambdaEff_MCKin_physicalPrimary = new THnSparseF("fDphiRecoHRealLambdaEff_MCKin_physicalPrimary", "Efficiency-corrected recoHadron-realLambda Correlation Histogram (using MC kinematics on V0, trigger and lambda are physical primary)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiRecoHRealLambdaEff_MCKin_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiRecoHRealLambdaEff_MCKin_physicalPrimary);

    fDphiHLambda_MC = new THnSparseF("fDphiHLambda_MC", "Hadron-Lambda Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambda_MC->Sumw2();
    fOutputList->Add(fDphiHLambda_MC);

    fDphiHDaughterProton_MC = new THnSparseF("fDphiHDaughterProton_MC", "Hadron-DaughterProton Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterProton_MC->Sumw2();
    fOutputList->Add(fDphiHDaughterProton_MC);

    fDphiHDaughterPion_MC = new THnSparseF("fDphiHDaughterPion_MC", "Hadron-DaughterPion Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterPion_MC->Sumw2();
    fOutputList->Add(fDphiHDaughterPion_MC);

    fDphiHLambda_MC_physicalPrimary = new THnSparseF("fDphiHLambda_MC_physicalPrimary", "Hadron-Lambda Correlation Histogram (MC truth, primary lambdas)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambda_MC_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiHLambda_MC_physicalPrimary);

    fDphiHLambdaMixed = new THnSparseF("fDphiHLambdaMixed", "Mixed Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed);

    fDphiHDaughterProtonMixed_MCKin = new THnSparseF("fDphiHDaughterProtonMixed_MCKin", "Mixed Hadron-daughter proton Correlation Histogram (using MC kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterProtonMixed_MCKin->Sumw2();
    fOutputList->Add(fDphiHDaughterProtonMixed_MCKin);

    fDphiHDaughterPionMixed_MCKin = new THnSparseF("fDphiHDaughterPionMixed_MCKin", "Mixed Hadron-daughter proton Correlation Histogram (using MC kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterPionMixed_MCKin->Sumw2();
    fOutputList->Add(fDphiHDaughterPionMixed_MCKin);

    fDphiHLambdaMixed_MCKin = new THnSparseF("fDphiHLambdaMixed_MCKin", "Mixed Hadron-Lambda Correlation Histogram (using Mc kinematics on V0)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MCKin->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MCKin);

    fDphiHLambdaMixed_MCKin_physicalPrimary = new THnSparseF("fDphiHLambdaMixed_MCKin_physicalPrimary", "Mixed Hadron-Lambda Correlation Histogram (using MC kinematics on V0, trigger and lambda are physical primary)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MCKin_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MCKin_physicalPrimary);

    fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary = new THnSparseF("fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary", "Mixed Hadron-Lambda Correlation Histogram (using MC kinematics on V0, reco trigger and real lambda are physical primary)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary);

    fDphiHLambdaMixed_MC = new THnSparseF("fDphiHLambdaMixed_MC", "Mixed Hadron-Lambda Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MC->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MC);

    fDphiHDaughterProtonMixed_MC = new THnSparseF("fDphiHDaughterProtonMixed_MC", "Hadron-DaughterProtonMixed Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterProtonMixed_MC->Sumw2();
    fOutputList->Add(fDphiHDaughterProtonMixed_MC);

    fDphiHDaughterPionMixed_MC = new THnSparseF("fDphiHDaughterPionMixed_MC", "Hadron-DaughterPionMixed Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHDaughterPionMixed_MC->Sumw2();
    fOutputList->Add(fDphiHDaughterPionMixed_MC);

    fDphiHLambdaMixed_MC_physicalPrimary = new THnSparseF("fDphiHLambdaMixed_MC_physicalPrimary", "Mixed Hadron-Lambda Correlation Histogram (MC truth, primary lambdas)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MC_physicalPrimary->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MC_physicalPrimary);


    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {18, 18, 16, 20, 10};
    double hh_cor_mins[5] = {1, 1, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {10, 10, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff->Sumw2();
    fOutputList->Add(fDphiHHEff);

    fDphiHHEff_checkMC = new THnSparseF("fDphiHHEff_checkMC", "Efficiency corrected Hadron-Hadron Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff_checkMC->Sumw2();
    fOutputList->Add(fDphiHHEff_checkMC);

    fDphiHProton = new THnSparseF("fDphiHProton", "Efficiency corrected Hadron-Proton Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHProton->Sumw2();
    fOutputList->Add(fDphiHProton);

    fDphiHPion = new THnSparseF("fDphiHPion", "Efficiency corrected Hadron-Pion Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHPion->Sumw2();
    fOutputList->Add(fDphiHPion);

    fDphiHH_MC = new THnSparseF("fDphiHH_MC", "Hadron-Hadron Correlation Histogram (MC truth)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHH_MC->Sumw2();
    fOutputList->Add(fDphiHH_MC);

    fDphiHProton_MC = new THnSparseF("fDphiHProton_MC", "Efficiency corrected Hadron-Proton_MC Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHProton_MC->Sumw2();
    fOutputList->Add(fDphiHProton_MC);

    fDphiHPion_MC = new THnSparseF("fDphiHPion_MC", "Efficiency corrected Hadron-Pion_MC Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHPion_MC->Sumw2();
    fOutputList->Add(fDphiHPion_MC);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed->Sumw2();
    fOutputList->Add(fDphiHHMixed);

    fDphiHProtonMixed = new THnSparseF("fDphiHProtonMixed", "Efficiency corrected Hadron-ProtonMixed Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHProtonMixed->Sumw2();
    fOutputList->Add(fDphiHProtonMixed);

    fDphiHPionMixed = new THnSparseF("fDphiHPionMixed", "Efficiency corrected Hadron-PionMixed Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHPionMixed->Sumw2();
    fOutputList->Add(fDphiHPionMixed);

    fDphiHHMixed_MC = new THnSparseF("fDphiHHMixed_MC", "Mixed Hadron-Hadron Correlation Histogram (MC truth)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed_MC->Sumw2();
    fOutputList->Add(fDphiHHMixed_MC);

    fDphiHProtonMixed_MC = new THnSparseF("fDphiHProtonMixed_MC", "Efficiency corrected Hadron-ProtonMixed_MC Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHProtonMixed_MC->Sumw2();
    fOutputList->Add(fDphiHProtonMixed_MC);

    fDphiHPionMixed_MC = new THnSparseF("fDphiHPionMixed_MC", "Efficiency corrected Hadron-PionMixed_MC Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHPionMixed_MC->Sumw2();
    fOutputList->Add(fDphiHPionMixed_MC);

    PostData(1, fOutputList);

}

void AliAnalysisTaskLambdaHadronV0Closure::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff)
{
    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        bool in_pt_range = (particle->Pt() < 10 && particle->Pt() > 0.5);
        if(trig_eff && in_pt_range) {
            int trigBin = fTriggerEff->FindBin(particle->Pt());
            double trigEff = fTriggerEff->GetBinContent(trigBin);
            double triggerScale = 1.0/trigEff;
            fDist->Fill(dist_points, triggerScale);
        }
        else{
            fDist->Fill(dist_points);
        }

    }
}

void AliAnalysisTaskLambdaHadronV0Closure::FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist)
{
    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        fDist->Fill(dist_points);
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda, bool lambdaEff)
{
    double dist_points[5]; //Pt, Phi, Eta, M, event multiplicity
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].vzero;
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        if(isAntiLambda) {
            dist_points[3] = particle->MassAntiLambda();
        }
        else{
            dist_points[3] = particle->MassLambda();
        }
        dist_points[4] = multPercentile;
        bool in_pt_range = (particle->Pt() < 10 && particle->Pt() > 0.5);
        if(lambdaEff && in_pt_range) {
            int lambdaBin = fLambdaEff->FindBin(particle->Pt());
            double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
            double lambdaScale = 1.0/lambdaEff;
            fDist->Fill(dist_points, lambdaScale);
        }
        else{
            fDist->Fill(dist_points);
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::FillMCMotherDist(std::vector<AliAODMCParticle*> particle_list, float multPercentile, THnSparse* fDist)
{
    double dist_points[5]; //Pt, Phi, Eta, M, event multiplicity
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = particle->M();
        dist_points[4] = multPercentile;
        fDist->Fill(dist_points);
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskLambdaHadronV0Closure::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskLambdaHadronV0Closure::SetAssociatedBit(float associatedBit) {
    fAssociatedBit = associatedBit;
}

void AliAnalysisTaskLambdaHadronV0Closure::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}

void AliAnalysisTaskLambdaHadronV0Closure::LoadEfficiencies(TString filePath) {
    TFile* effFile = TFile::Open(filePath);

    if(!effFile) {
        AliFatal("NULL INPUT FILE WHEN LOADING EFFICIENCIES, EXITING");
    }

    fLambdaEff = (TH1D*) effFile->Get("fLambdaV0Eff")->Clone("fLambdaV0EffClone");
    if(!fLambdaEff) {
        AliFatal("UNABLE TO FIND LAMBDA EFF, EXITING");
    }
    
    fAssociatedEff = (TH1D*) effFile->Get("fAssociatedEff")->Clone("fAssociatedEffClone");
    if(!fAssociatedEff) {
        AliFatal("UNABLE TO FIND ASSOCIATED EFF, EXITING");
    }

    fTriggerEff = (TH1D*) effFile->Get("fTriggerEff")->Clone("fTriggerEffClone");
    if(!fTriggerEff) {
        AliFatal("UNABLE TO FIND TRIGGER EFF, EXITING");
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetID() == lambda.daughter1ID) || (trigger->GetID() == lambda.daughter2ID)) continue;

            dphi_point[1] = lambda.vzero->Pt();
            dphi_point[2] = trigger->Phi() - lambda.vzero->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - lambda.vzero->Eta();
            if(isAntiLambda) dphi_point[4] = lambda.vzero->MassAntiLambda();
            else dphi_point[4] = lambda.vzero->MassLambda();
            dphi_point[5] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (lambda.vzero->Pt() < 10 && lambda.vzero->Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int lambdaBin = fLambdaEff->FindBin(lambda.vzero->Pt());
                double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
                double lambdaScale = 1.0/lambdaEff;
                double totalScale = triggerScale*lambdaScale;
                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}
void AliAnalysisTaskLambdaHadronV0Closure::MakeSameRecoHRealLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetLabel() == lambda->GetDaughterFirst()) || (trigger->GetID() == lambda->GetDaughterLast())) continue;

            dphi_point[1] = lambda->Pt();
            dphi_point[2] = trigger->Phi() - lambda->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - lambda->Eta();
            dphi_point[4] = lambda->M();
            dphi_point[5] = zVtx;

            bool in_pt_range = (trigger->Pt() < 10 && trigger->Pt() > 0.5);

            if(eff && in_pt_range) {
                // In this function we only want to use the efficiency for the trigger
                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                fDphi->Fill(dphi_point, triggerScale);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeSameHLambdaCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {

            auto lambda = lambda_list[i];
            AliAODTrack *posTrack=(AliAODTrack *)lambda.vzero->GetDaughter(0);
                    
            int plabel = posTrack->GetLabel();
            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(plabel);
            int mlabel_pos = mcpospart->GetMother();
            AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(mlabel_pos);

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetID() == lambda.daughter1ID) || (trigger->GetID() == lambda.daughter2ID)) continue;

            dphi_point[1] = mcmother->Pt();
            dphi_point[2] = trigger->Phi() - mcmother->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - mcmother->Eta();
            dphi_point[4] = mcmother->M();
            dphi_point[5] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (mcmother->Pt() < 10 && mcmother->Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int lambdaBin = fLambdaEff->FindBin(mcmother->Pt());
                double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
                double lambdaScale = 1.0/lambdaEff;
                double totalScale = triggerScale*lambdaScale;
                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeSameHDaughterCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {

            auto lambda = lambda_list[i];

            AliAODTrack *posTrack=(AliAODTrack *)lambda.vzero->GetDaughter(0);
            int plabel = posTrack->GetLabel();
            AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(plabel);
            int mlabel_pos = mcpospart->GetMother();

            AliAODTrack *negTrack=(AliAODTrack *)lambda.vzero->GetDaughter(1);
            int nlabel = negTrack->GetLabel();
            AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(nlabel);
            int mlabel_neg = mcnegpart->GetMother();

            AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(mlabel_pos);
            if(mcmother->Pt() > 4.0 || mcmother->Pt() < 2.0) continue;

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetID() == lambda.daughter1ID) || (trigger->GetID() == lambda.daughter2ID)) continue;

            if(mcpospart->GetPdgCode() == 211) {
                dphi_point[1] = mcpospart->Pt();
                dphi_point[2] = trigger->Phi() - mcpospart->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - mcpospart->Eta();
                dphi_point[4] = mcpospart->M();
                dphi_point[5] = zVtx;
                fDphi_pion->Fill(dphi_point);
            }

            if(mcnegpart->GetPdgCode() == -211) {
                dphi_point[1] = mcnegpart->Pt();
                dphi_point[2] = trigger->Phi() - mcnegpart->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - mcnegpart->Eta();
                dphi_point[4] = mcnegpart->M();
                dphi_point[5] = zVtx;
                fDphi_pion->Fill(dphi_point);
            }

            if(mcpospart->GetPdgCode() == 2212) {
                dphi_point[1] = mcpospart->Pt();
                dphi_point[2] = trigger->Phi() - mcpospart->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - mcpospart->Eta();
                dphi_point[4] = mcpospart->M();
                dphi_point[5] = zVtx;
                fDphi_proton->Fill(dphi_point);
            }

            if(mcnegpart->GetPdgCode() == -2212) {
                dphi_point[1] = mcnegpart->Pt();
                dphi_point[2] = trigger->Phi() - mcnegpart->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - mcnegpart->Eta();
                dphi_point[4] = mcnegpart->M();
                dphi_point[5] = zVtx;
                fDphi_proton->Fill(dphi_point);
            }

        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeSameMCHLambdaCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];

            int first_daughter_index = lambda->GetDaughterFirst();
            int second_daughter_index = lambda->GetDaughterLast();

            // guaranteed to exist since lambda already passes cuts at this point
            AliAODMCParticle* first_daughter = (AliAODMCParticle*)fMCArray->At(first_daughter_index);
            AliAODMCParticle* second_daughter = (AliAODMCParticle*)fMCArray->At(second_daughter_index);

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetLabel() == first_daughter->GetLabel()) || (trigger->GetLabel() == second_daughter->GetLabel())) continue;

            dphi_point[1] = lambda->Pt();
            dphi_point[2] = trigger->Phi() - lambda->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - lambda->Eta();
            dphi_point[4] = lambda->M();
            dphi_point[5] = zVtx;
            fDphi->Fill(dphi_point);
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeSameMCHDaughterCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];
            if(lambda->Pt() > 4.0 || lambda->Pt() < 2.0) continue;

            int first_daughter_index = lambda->GetDaughterFirst();
            int second_daughter_index = lambda->GetDaughterLast();

            // GetDaughterFirst always returns proton
            // GetDaughterLast always returns pion
            AliAODMCParticle* dproton = (AliAODMCParticle*)fMCArray->At(first_daughter_index);
            AliAODMCParticle* dpion = (AliAODMCParticle*)fMCArray->At(second_daughter_index);

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetLabel() == dproton->GetLabel()) || (trigger->GetLabel() == dpion->GetLabel())) continue;

            // fill proton dist first
            dphi_point[1] = dproton->Pt();
            dphi_point[2] = trigger->Phi() - dproton->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - dproton->Eta();
            dphi_point[4] = dproton->M();
            dphi_point[5] = zVtx;
            fDphi_proton->Fill(dphi_point);

            // now fill pion dist
            dphi_point[1] = dpion->Pt();
            dphi_point[2] = trigger->Phi() - dpion->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }
            dphi_point[3] = trigger->Eta() - dpion->Eta();
            dphi_point[4] = dpion->M();
            dphi_point[5] = zVtx;
            fDphi_pion->Fill(dphi_point);
        }
    }
}


void AliAnalysisTaskLambdaHadronV0Closure::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_h_list.size(); i++) {
            auto associate = associated_h_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;

                int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;

                double totalScale = triggerScale*associatedScale;

                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}
void AliAnalysisTaskLambdaHadronV0Closure::MakeSameMCHHCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_h_list.size(); i++) {
            auto associate = associated_h_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));

            fDphi->Fill(dphi_point);
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)lambda_list.size(); j++) {
                auto lambda = lambda_list[j];

                dphi_point[1] = lambda.vzero->Pt();
                dphi_point[2] = trigger->Phi() - lambda.vzero->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - lambda.vzero->Eta();
                
                if(isAntiLambda) dphi_point[4] = lambda.vzero->MassAntiLambda();
                else dphi_point[4] = lambda.vzero->MassLambda();

                dphi_point[5] = zVtx;
                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (lambda.vzero->Pt() < 10 && lambda.vzero->Pt() > 0.5));
                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int lambdaBin = fLambdaEff->FindBin(lambda.vzero->Pt());
                    double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
                    double lambdaScale = 1.0/lambdaEff;
                    double totalScale = triggerScale*lambdaScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedHLambdaCorrelations_withMCKin(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)lambda_list.size(); j++) {
                auto lambda = (AliAODMCParticle*)fMCArray->At(lambda_list[j].motherLabel);

                dphi_point[1] = lambda->Pt();
                dphi_point[2] = trigger->Phi() - lambda->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - lambda->Eta();
                dphi_point[4] = lambda->M();
                dphi_point[5] = zVtx;
                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (lambda->Pt() < 10 && lambda->Pt() > 0.5));
                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int lambdaBin = fLambdaEff->FindBin(lambda->Pt());
                    double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
                    double lambdaScale = 1.0/lambdaEff;
                    double totalScale = triggerScale*lambdaScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedHDaughterCorrelations_withMCKin(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list , THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)lambda_list.size(); j++) {
                auto lambda = (AliAODMCParticle*)fMCArray->At(lambda_list[j].motherLabel);
                // it just so happens GetDaughterFirst always returns proton and GetDaughterLast always returns pion
                if(lambda->Pt() > 4.0 || lambda->Pt() < 2.0) continue;
                auto dproton = (AliAODMCParticle*)fMCArray->At(lambda->GetDaughterFirst());
                auto dpion = (AliAODMCParticle*)fMCArray->At(lambda->GetDaughterLast());

                // fill proton dist first
                dphi_point[1] = dproton->Pt();
                dphi_point[2] = trigger->Phi() - dproton->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - dproton->Eta();
                dphi_point[4] = dproton->M();
                dphi_point[5] = zVtx;
                fDphi_proton->Fill(dphi_point);

                // now fill pion dist
                dphi_point[1] = dpion->Pt();
                dphi_point[2] = trigger->Phi() - dpion->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - dpion->Eta();
                dphi_point[4] = dpion->M();
                dphi_point[5] = zVtx;
                fDphi_pion->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedMCHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> lambda_list , THnSparse* fDphi, double zVtx)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)lambda_list.size(); j++) {
                auto lambda = lambda_list[j];

                dphi_point[1] = lambda->Pt();
                dphi_point[2] = trigger->Phi() - lambda->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - lambda->Eta();
                dphi_point[4] = lambda->M();
                dphi_point[5] = zVtx;
                fDphi->Fill(dphi_point);
            }
        }
    }
}
void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedMCHDaughterCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> lambda_list , THnSparse* fDphi_proton, THnSparse* fDphi_pion, double zVtx)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)lambda_list.size(); j++) {
                auto lambda = lambda_list[j];
                if(lambda->Pt() > 4.0 || lambda->Pt() < 2.0) continue;
                // it just so happens GetDaughterFirst always returns proton and GetDaughterLast always returns pion
                auto dproton = (AliAODMCParticle*)fMCArray->At(lambda->GetDaughterFirst());
                auto dpion = (AliAODMCParticle*)fMCArray->At(lambda->GetDaughterLast());

                // fill proton dist first
                dphi_point[1] = dproton->Pt();
                dphi_point[2] = trigger->Phi() - dproton->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - dpion->Eta();
                dphi_point[4] = dpion->M();
                dphi_point[5] = zVtx;
                fDphi_pion->Fill(dphi_point);

                // now fill pion dist
                dphi_point[1] = dpion->Pt();
                dphi_point[2] = trigger->Phi() - dpion->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - dpion->Eta();
                dphi_point[4] = dpion->M();
                dphi_point[5] = zVtx;
                fDphi_pion->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_h_list.size(); j++) {
                auto associate = associated_h_list[j];

                dphi_point[1] = associate->Pt();
                dphi_point[2] = trigger->Phi() - associate->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - associate->Eta();
                dphi_point[4] = zVtx;

                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (associate->Pt() < 10 && associate->Pt() > 0.5));

                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                    double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                    double associatedScale = 1.0/associatedEff;
                    double totalScale = triggerScale*associatedScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronV0Closure::MakeMixedMCHHCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx)
{
    double dphi_point[5];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_h_list.size(); j++) {
                auto associate = associated_h_list[j];

                dphi_point[1] = associate->Pt();
                dphi_point[2] = trigger->Phi() - associate->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - associate->Eta();
                dphi_point[4] = zVtx;

                fDphi->Fill(dphi_point);
            }
        }
    }
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassDaughterCuts(AliAODTrack *track){

    if(track->GetID() < 0) return false;

    bool pass = true;

    pass = pass && (TMath::Abs(track->Eta()) < 0.8);
    pass = pass && (track->Pt() > 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    return pass;
}

uint8_t AliAnalysisTaskLambdaHadronV0Closure::PassV0LambdaCuts(AliAODv0 *v0, bool checkMotherPDG, bool checkPhysicalPrimary) {

    if(v0->GetOnFlyStatus()) return 0;
    if(!(TMath::Abs(v0->Eta()) < 0.8)) return 0;

    AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);
    AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);

    if(!PassDaughterCuts(ptrack)) return 0;
    if(!PassDaughterCuts(ntrack)) return 0;
        
    int plabel = ptrack->GetLabel();
    int nlabel = ntrack->GetLabel();

    if(plabel < 0 || nlabel < 0) return 0;

    AliAODMCParticle* mcpospart = (AliAODMCParticle*)fMCArray->At(plabel);
    AliAODMCParticle* mcnegpart = (AliAODMCParticle*)fMCArray->At(nlabel);

    if(checkMotherPDG) {
        int mlabel_pos = mcpospart->GetMother();
        int mlabel_neg = mcnegpart->GetMother();

        if(mlabel_pos < 0 || mlabel_neg < 0) return 0;
        if(mlabel_pos != mlabel_neg) return 0;

        AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(mlabel_pos);

        int momPDG = mcmother->GetPdgCode();
        if(TMath::Abs(momPDG) != 3122) return 0;

        if(checkPhysicalPrimary) {
            if(!mcmother->IsPhysicalPrimary()) return 0;
        }
    }

    int posPDG = mcpospart->GetPdgCode();
    int negPDG = mcnegpart->GetPdgCode();

    if(posPDG == 2212 && negPDG == -211) {
        return 1;
    } else if(posPDG == 211 && negPDG == -2212) {
        return 2;
    }
    else {
        return 0;
    }
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassAssociatedCuts(AliAODTrack *track, bool checkMC){

    if(!(TMath::Abs(track->Eta()) < 0.8)) return false;
    if(!(track->Pt() > 0.15)) return false;
    if(!track->TestFilterMask(fAssociatedBit)) return false;

    if(checkMC) {
        int label = track->GetLabel();
        if(label < 0) return false;

        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(label);

        int pdg = mcpart->GetPdgCode();
        if(!IsMCChargedHadron(pdg)) return false;
        if(!mcpart->IsPhysicalPrimary()) return false;
    }

    return true;
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassTriggerCuts(AliAODTrack *track, bool checkMC){

    if(!(TMath::Abs(track->Eta()) < 0.8)) return false;
    if(!(track->Pt() > 0.15)) return false;
    if(!track->TestBit(fTriggerBit)) return false;

    if(checkMC) {
        int label = track->GetLabel();
        if(label < 0) return false;

        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(label);

        int pdg = mcpart->GetPdgCode();
        if(!IsMCChargedHadron(pdg)) return false;
        if(!mcpart->IsPhysicalPrimary()) return false;
    }

    return true;
}

bool AliAnalysisTaskLambdaHadronV0Closure::IsMCChargedHadron(int pdg_code) {
    // checks to see if pdg code matches pion, proton, kaon, electron or muon
    if((TMath::Abs(pdg_code) == 321)
        || (TMath::Abs(pdg_code) == 211)
        || (TMath::Abs(pdg_code) == 2212)
        || (TMath::Abs(pdg_code) == 11)
        || (TMath::Abs(pdg_code) == 13)) return true;

    return false;
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassMCTriggerCuts(AliAODMCParticle *mc_particle){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false; // for now trigger is physical primary, could change
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassMCAssociatedCuts(AliAODMCParticle *mc_particle){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false;
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronV0Closure::PassMCLambdaCuts(AliAODMCParticle *mc_particle, bool checkPhysicalPrimary){

    if(!(TMath::Abs(mc_particle->GetPdgCode()) == 3122)) return false;
    if(checkPhysicalPrimary) {
        if(!mc_particle->IsPhysicalPrimary()) return false;
    }
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;

    int first_daughter_index = 0;
    int second_daughter_index = 0;

    first_daughter_index = mc_particle->GetDaughterFirst();
    second_daughter_index = mc_particle->GetDaughterLast();

    if(first_daughter_index < 0 || second_daughter_index < 0) return false;

    AliAODMCParticle* first_daughter = (AliAODMCParticle*)fMCArray->At(first_daughter_index);
    AliAODMCParticle* second_daughter = (AliAODMCParticle*)fMCArray->At(second_daughter_index);

    // make sure lambda decays into p-pi
    if(!(((TMath::Abs(first_daughter->GetPdgCode()) == 211 && TMath::Abs(second_daughter->GetPdgCode()) == 2212) 
        ||  (TMath::Abs(first_daughter->GetPdgCode()) == 2212 && TMath::Abs(second_daughter->GetPdgCode()) == 211)) 
        && (first_daughter->GetPdgCode())*(second_daughter->GetPdgCode()) < 0)) return false;
    
    return true;
}
        

void AliAnalysisTaskLambdaHadronV0Closure::UserExec(Option_t*)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
    }

    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
        AliError("Array of MC particles not found");
        return;
    }

    fpidResponse = fInputHandler->GetPIDResponse();

    //Event cuts
    TString cent_estimator = fCentEstimator;
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < fMultLow || multPercentile > fMultHigh) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();

    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;

    std::vector<AliAODTrack*> trigger_list_checkMC;
    std::vector<AliAODTrack*> associated_h_list_checkMC;

    std::vector<AliAODTrack*> associated_proton_list;
    std::vector<AliAODTrack*> associated_pion_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    //MC trigger list used for event mixing
    TObjArray* fMixedMCTrackObjArray = new TObjArray;
    fMixedMCTrackObjArray->SetOwner(kTRUE);

    // Bool to keep track if the event has a high-pt (> 4 GeV) trigger
    bool is_triggered_event = false;

    float maxTrigPt = 0;
    AliAODTrack* maxTrigger = 0x0;

    int NCharged = 0;


    // RECO SAME EVENT SECTION

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);
            if(triggerPart->Pt() > 4) is_triggered_event = true;
            if(triggerPart->Pt() > 4 && triggerPart->Pt() < 8) {
                if(triggerPart->Pt() > maxTrigPt) {
                    maxTrigPt = triggerPart->Pt();
                    maxTrigger = track;
                }
            }
        }

        if(PassAssociatedCuts(track)) {
            associated_h_list.push_back(track);
        }

        if(PassTriggerCuts(track, true)) {
            trigger_list_checkMC.push_back(track);
        }

        if(PassAssociatedCuts(track, true)) {
            associated_h_list_checkMC.push_back(track);
        }

        if(PassDaughterCuts(track)) {
            int mc_label = track->GetLabel();
            if(mc_label >= 0) {
                auto mc_part = (AliAODMCParticle*)fMCArray->At(mc_label);
                if(TMath::Abs(mc_part->GetPdgCode()) == 2212) {
                    associated_proton_list.push_back(track);
                }
                if(TMath::Abs(mc_part->GetPdgCode()) == 211) {
                    associated_pion_list.push_back(track);
                }
            }
        }
    }

    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> antilambda_list;
    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list;
    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> antilambda_list_checkMotherPDG;
    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list_checkMotherPDG;
    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> antilambda_list_checkMotherPDG_isPrimary;
    std::vector<AliAnalysisTaskLambdaHadronV0Closure::AliMotherContainer> lambda_list_checkMotherPDG_isPrimary;


    // V0 SECTION
    int numV0s = fAOD->GetNumberOfV0s();
    for(int i = 0; i < numV0s; i++) {
        AliAODv0 *v0 = fAOD->GetV0(i);

        AliAODTrack *posTrack=(AliAODTrack *)v0->GetDaughter(0);
        AliAODTrack *negTrack=(AliAODTrack *)v0->GetDaughter(1);

        if(PassV0LambdaCuts(v0) == 1) {
            AliMotherContainer lambda;
            lambda.vzero = v0;
            lambda.daughter1ID = posTrack->GetID();
            lambda.daughter2ID = negTrack->GetID();
            // set mother label to 1 when not requiring MC check
            lambda.motherLabel = -1;
            lambda_list.push_back(lambda);
        }

        if(PassV0LambdaCuts(v0) == 2) {
            AliMotherContainer antilambda;
            antilambda.vzero = v0;
            antilambda.daughter1ID = posTrack->GetID();
            antilambda.daughter2ID = negTrack->GetID();
            // set mother label to 1 when not requiring MC check
            antilambda.motherLabel = -1;
            antilambda_list.push_back(antilambda);
        }

        if(PassV0LambdaCuts(v0, true) == 1) {
            AliMotherContainer lambda;
            auto mc_pos = (AliAODMCParticle*)fMCArray->At(posTrack->GetLabel());
            int mother_label = mc_pos->GetMother();
            lambda.vzero = v0;
            lambda.daughter1ID = posTrack->GetID();
            lambda.daughter2ID = negTrack->GetID();
            lambda.motherLabel = mother_label;
            lambda_list_checkMotherPDG.push_back(lambda);
        }

        if(PassV0LambdaCuts(v0, true) == 2) {
            AliMotherContainer antilambda;
            auto mc_pos = (AliAODMCParticle*)fMCArray->At(posTrack->GetLabel());
            int mother_label = mc_pos->GetMother();
            antilambda.vzero = v0;
            antilambda.daughter1ID = posTrack->GetID();
            antilambda.daughter2ID = negTrack->GetID();
            antilambda.motherLabel = mother_label;
            antilambda_list_checkMotherPDG.push_back(antilambda);
        }

        if(PassV0LambdaCuts(v0, true, true) == 1) {
            AliMotherContainer lambda;
            auto mc_pos = (AliAODMCParticle*)fMCArray->At(posTrack->GetLabel());
            int mother_label = mc_pos->GetMother();
            lambda.vzero = v0;
            lambda.daughter1ID = posTrack->GetID();
            lambda.daughter2ID = negTrack->GetID();
            lambda.motherLabel = mother_label;
            lambda_list_checkMotherPDG_isPrimary.push_back(lambda);
        }

        if(PassV0LambdaCuts(v0, true, true) == 2) {
            AliMotherContainer antilambda;
            auto mc_pos = (AliAODMCParticle*)fMCArray->At(posTrack->GetLabel());
            int mother_label = mc_pos->GetMother();
            antilambda.vzero = v0;
            antilambda.daughter1ID = posTrack->GetID();
            antilambda.daughter2ID = negTrack->GetID();
            antilambda.motherLabel = mother_label;
            antilambda_list_checkMotherPDG_isPrimary.push_back(antilambda);
        }

    }


    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDistEff, true);
    FillSingleParticleDist(trigger_list, primZ, fTriggerDist, false);
    FillSingleParticleDist(trigger_list_checkMC, primZ, fTriggerDistEff_checkMC, true);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);
    FillSingleParticleDist(associated_h_list_checkMC, primZ, fAssociatedHDist);
    FillMotherDist(lambda_list, primZ, fLambdaDist, false);
    FillMotherDist(antilambda_list, primZ, fLambdaDist, true);
    FillMotherDist(lambda_list_checkMotherPDG, primZ, fGuaranteedLambdaDist, false);
    FillMotherDist(antilambda_list_checkMotherPDG, primZ, fGuaranteedLambdaDist, true);

    // Filling our single particle lambda distribution histogram:
    if(is_triggered_event) FillMotherDist(lambda_list, multPercentile, fTriggeredLambdaDist, false);
    if(is_triggered_event) FillMotherDist(antilambda_list, multPercentile, fTriggeredLambdaDist, true);

    MakeSameHLambdaCorrelations(trigger_list, antilambda_list, fDphiHLambdaEff, primZ, true, true);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambdaEff, primZ, true, false);
    MakeSameHLambdaCorrelations(trigger_list, antilambda_list_checkMotherPDG, fDphiHGuaranteedLambdaEff, primZ, true, true);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_checkMotherPDG, fDphiHGuaranteedLambdaEff, primZ, true, false);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list, antilambda_list_checkMotherPDG, fDphiHLambdaEff_MCKin, primZ, true);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list, lambda_list_checkMotherPDG, fDphiHLambdaEff_MCKin, primZ, true);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list_checkMC, antilambda_list_checkMotherPDG_isPrimary, fDphiHLambdaEff_MCKin_physicalPrimary, primZ, true);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list_checkMC, lambda_list_checkMotherPDG_isPrimary, fDphiHLambdaEff_MCKin_physicalPrimary, primZ, true);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, true);
    MakeSameHHCorrelations(trigger_list_checkMC, associated_h_list_checkMC, fDphiHHEff_checkMC, primZ, true);

    MakeSameHHCorrelations(trigger_list_checkMC, associated_proton_list, fDphiHProton, primZ, true);
    MakeSameHHCorrelations(trigger_list_checkMC, associated_pion_list, fDphiHPion, primZ, true);

    MakeSameHDaughterCorrelations_withMCKin(trigger_list, antilambda_list_checkMotherPDG, fDphiHDaughterProton_MCKin, fDphiHDaughterPion_MCKin, primZ);
    MakeSameHDaughterCorrelations_withMCKin(trigger_list, lambda_list_checkMotherPDG, fDphiHDaughterProton_MCKin, fDphiHDaughterPion_MCKin, primZ);


    // MC SAME EVENT SECTION

    std::vector<AliAODMCParticle*> real_trigger_list;
    std::vector<AliAODMCParticle*> real_associated_list;
    std::vector<AliAODMCParticle*> real_proton_list;
    std::vector<AliAODMCParticle*> real_pion_list;
    std::vector<AliAODMCParticle*> real_lambda_list;
    std::vector<AliAODMCParticle*> real_lambda_list_physicalPrimary;

    bool is_triggered_event_MC = false;

    for(int mc_index = 0; mc_index < fMCArray->GetEntries(); mc_index++){

        AliAODMCParticle *mc_particle = (AliAODMCParticle*)fMCArray->At(mc_index);

        if(PassMCTriggerCuts(mc_particle)) {
            real_trigger_list.push_back(mc_particle);
            if(mc_particle->Pt() > 4) is_triggered_event_MC = true;
            AliCFParticle *trigger_particle = new AliCFParticle(mc_particle->Pt(), mc_particle->Eta(), mc_particle->Phi(), mc_particle->Charge(), 0);
            fMixedMCTrackObjArray->Add(trigger_particle);
        }
        if(PassMCAssociatedCuts(mc_particle)) real_associated_list.push_back(mc_particle);
        if(PassMCLambdaCuts(mc_particle)) real_lambda_list.push_back(mc_particle);
        if(PassMCLambdaCuts(mc_particle, true)) real_lambda_list_physicalPrimary.push_back(mc_particle);
        if(TMath::Abs(mc_particle->Eta()) < 0.8 && mc_particle->Pt() > 0.15) {
            if(TMath::Abs(mc_particle->GetPdgCode()) == 2212) real_proton_list.push_back(mc_particle);
            if(TMath::Abs(mc_particle->GetPdgCode()) == 211) real_pion_list.push_back(mc_particle);
        }
    }

    FillSingleMCParticleDist(real_trigger_list, primZ, fTriggerDist_MC);
    FillSingleMCParticleDist(real_associated_list, primZ, fAssociatedDist_MC);
    FillMCMotherDist(real_lambda_list, primZ, fLambdaDist_MC);
    if(is_triggered_event_MC) FillMCMotherDist(real_lambda_list, multPercentile, fTriggeredLambdaDist_MC);

    MakeSameMCHLambdaCorrelations(real_trigger_list, real_lambda_list, fDphiHLambda_MC, primZ);
    MakeSameMCHLambdaCorrelations(real_trigger_list, real_lambda_list_physicalPrimary, fDphiHLambda_MC_physicalPrimary, primZ);
    MakeSameMCHHCorrelations(real_trigger_list, real_associated_list, fDphiHH_MC, primZ);

    MakeSameMCHHCorrelations(real_trigger_list, real_proton_list, fDphiHProton_MC, primZ);
    MakeSameMCHHCorrelations(real_trigger_list, real_pion_list, fDphiHPion_MC, primZ);

    MakeSameMCHDaughterCorrelations(real_trigger_list, real_lambda_list, fDphiHDaughterProton_MC, fDphiHDaughterPion_MC, primZ);

    MakeSameRecoHRealLambdaCorrelations(trigger_list, real_lambda_list, fDphiRecoHRealLambdaEff_MCKin_physicalPrimary, primZ, true);

    // MIXED EVENT SECTION (added to very end to correctly do a mixture of reco/real correlations)

    if(lambda_list.size() > 0  || associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fCorPool->IsReady()) {
                MakeMixedHLambdaCorrelations(fCorPool, antilambda_list, fDphiHLambdaMixed, primZ, true, true);
                MakeMixedHLambdaCorrelations(fCorPool, lambda_list, fDphiHLambdaMixed, primZ, true, false);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, antilambda_list_checkMotherPDG, fDphiHLambdaMixed_MCKin, primZ, true);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, lambda_list_checkMotherPDG, fDphiHLambdaMixed_MCKin, primZ, true);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, antilambda_list_checkMotherPDG_isPrimary, fDphiHLambdaMixed_MCKin_physicalPrimary, primZ, true);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, lambda_list_checkMotherPDG_isPrimary, fDphiHLambdaMixed_MCKin_physicalPrimary, primZ, true);
                MakeMixedMCHLambdaCorrelations(fCorPool, real_lambda_list_physicalPrimary, fDphiRecoHRealLambdaMixed_MCKin_physicalPrimary, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_proton_list, fDphiHProtonMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_pion_list, fDphiHPionMixed, primZ);
                MakeMixedHDaughterCorrelations_withMCKin(fCorPool, antilambda_list_checkMotherPDG, fDphiHDaughterProtonMixed_MCKin, fDphiHDaughterPionMixed_MCKin, primZ);
                MakeMixedHDaughterCorrelations_withMCKin(fCorPool, lambda_list_checkMotherPDG, fDphiHDaughterProtonMixed_MCKin, fDphiHDaughterPionMixed_MCKin, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }


    if(real_associated_list.size() > 0 || real_lambda_list.size() > 0) {
        AliEventPool *fMCCorPool = fMCCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fMCCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fMCCorPool->IsReady()) {
                MakeMixedMCHLambdaCorrelations(fMCCorPool, real_lambda_list, fDphiHLambdaMixed_MC, primZ);
                MakeMixedMCHLambdaCorrelations(fMCCorPool, real_lambda_list_physicalPrimary, fDphiHLambdaMixed_MC_physicalPrimary, primZ);
                MakeMixedMCHHCorrelations(fMCCorPool, real_associated_list, fDphiHHMixed_MC, primZ);
                MakeMixedMCHHCorrelations(fMCCorPool, real_proton_list, fDphiHProtonMixed_MC, primZ);
                MakeMixedMCHHCorrelations(fMCCorPool, real_pion_list, fDphiHPionMixed_MC, primZ);
                MakeMixedMCHDaughterCorrelations(fMCCorPool, real_lambda_list, fDphiHDaughterProtonMixed_MC, fDphiHDaughterPionMixed_MC, primZ);
            }
            if(fMixedMCTrackObjArray->GetEntries() > 0) {
                fMCCorPool->UpdatePool(fMixedMCTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronV0Closure::Terminate(Option_t *option)
{
}
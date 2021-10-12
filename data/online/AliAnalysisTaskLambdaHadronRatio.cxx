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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
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

#include "AliAnalysisTaskLambdaHadronRatio.h"

ClassImp(AliAnalysisTaskLambdaHadronRatio);

AliAnalysisTaskLambdaHadronRatio::AliAnalysisTaskLambdaHadronRatio() :
    AliAnalysisTaskSE(),
    fAOD(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggersAndLambdasPerEvent_All(0x0),
    fTriggersAndLambdasPerEvent_2_4(0x0),
    fLooseDist(0x0),
    fTriggerDist(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_highestPt(0x0),
    fAssociatedHDist(0x0),
    fLambdaDist(0x0),
    fTriggeredLambdaDist(0x0),
    fTriggeredLambdaDistFilterbit(0x0),
    fDphiHLambda(0x0),
    fDphiHLambdaFilterbit(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_highestPt(0x0),
    fDphiHLambdaV0(0x0),
    fDphiHLambdaRotated(0x0),
    fDphiHLambdaRotatedPi(0x0),
    fDphiHLambdaRotatedProton(0x0),
    fDphiHLambdaFlipped(0x0),
    fDphiHH(0x0),
    fDphiHHEff(0x0),
    fDphiTriggerTrigger(0x0),
    fDphiHLambdaLS(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_highestPt(0x0),
    fDphiHHMixed(0x0),
    fDphiHHMixed_highestPt(0x0),
    fDphiHLambdaLSMixed(0x0),
    fDphiTriggerTriggerMixed(0x0),
    fLambdaDaughterDCA(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0),
    fTofTest(0x0),
    fAssociatedPtEventClass(0x0),
    fLambdaPtEventClass(0x0),
    fTriggerPtEventClass(0x0)
{
}

AliAnalysisTaskLambdaHadronRatio::AliAnalysisTaskLambdaHadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggersAndLambdasPerEvent_All(0x0),
    fTriggersAndLambdasPerEvent_2_4(0x0),
    fLooseDist(0x0),
    fTriggerDist(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_highestPt(0x0),
    fAssociatedHDist(0x0),
    fLambdaDist(0x0),
    fTriggeredLambdaDist(0x0),
    fTriggeredLambdaDistFilterbit(0x0),
    fDphiHLambda(0x0),
    fDphiHLambdaFilterbit(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_highestPt(0x0),
    fDphiHLambdaV0(0x0),
    fDphiHLambdaRotated(0x0),
    fDphiHLambdaRotatedPi(0x0),
    fDphiHLambdaRotatedProton(0x0),
    fDphiHLambdaFlipped(0x0),
    fDphiHH(0x0),
    fDphiHHEff(0x0),
    fDphiTriggerTrigger(0x0),
    fDphiHLambdaLS(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_highestPt(0x0),
    fDphiHHMixed(0x0),
    fDphiHHMixed_highestPt(0x0),
    fDphiHLambdaLSMixed(0x0),
    fDphiTriggerTriggerMixed(0x0),
    fLambdaDaughterDCA(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0),
    fTofTest(0x0),
    fAssociatedPtEventClass(0x0),
    fLambdaPtEventClass(0x0),
    fTriggerPtEventClass(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronRatio::~AliAnalysisTaskLambdaHadronRatio()
{
    if(fOutputList) delete fOutputList;
    if(fTriggerEff) delete fTriggerEff;
    if(fAssociatedEff) delete fAssociatedEff;
    if(fLambdaEff) delete fLambdaEff;
}

void AliAnalysisTaskLambdaHadronRatio::UserCreateOutputObjects()
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

    fTriggersAndLambdasPerEvent_All = new TH2D("fTriggersAndLambdasPerEvent_All", "Triggers and Lambdas per event (all p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndLambdasPerEvent_All);

    fTriggersAndLambdasPerEvent_2_4 = new TH2D("fTriggersAndLambdasPerEvent_2_4", "Triggers and Lambdas per event (2-4 p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndLambdasPerEvent_2_4);


    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {200, 16, 20, 10};
    double dist_mins[4] = {0.0, 0, -1, -10};
    double dist_maxes[4] = {20.0, 6.28, 1, 10};

    fLooseDist = new THnSparseF("fLooseDist", "All Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fLooseDist->Sumw2();
    fOutputList->Add(fLooseDist);

    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist->Sumw2();
    fOutputList->Add(fTriggerDist);

    fTriggerDistEff = new THnSparseF("fTriggerDistEff", "Efficiency Corrected Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDistEff->Sumw2();
    fOutputList->Add(fTriggerDistEff);

    fTriggerDistEff_highestPt = new THnSparseF("fTriggerDistEff_highestPt", "Efficiency Corrected Highest p_{t} Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDistEff_highestPt->Sumw2();
    fOutputList->Add(fTriggerDistEff_highestPt);

    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist->Sumw2();
    fOutputList->Add(fAssociatedHDist);

    //Mother distribution axes are: Pt, Phi, Eta, Mass, Event multiplicity
    int mother_dist_bins[5] = {100, 16, 20, 100, 10};
    double mother_dist_mins[5] = {0, -3.14, -1, 1.06, 0};
    double mother_dist_maxes[5] = {15, 3.14, 1, 1.16, 100};

    fLambdaDist = new THnSparseF("fLambdaDist", "Lambda Distribution", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fLambdaDist->Sumw2();
    fOutputList->Add(fLambdaDist);

    fTriggeredLambdaDist = new THnSparseF("fTriggeredLambdaDist", "Lambda Distribution (with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDist->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist);

    fTriggeredLambdaDistFilterbit = new THnSparseF("fTriggeredLambdaDistFilter", "Lambda Distribution (with triggered event, filterbit 16 on daughters)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDistFilterbit->Sumw2();
    fOutputList->Add(fTriggeredLambdaDistFilterbit);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hl_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    double hl_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hl_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10};

    fDphiHLambda = new THnSparseF("fDphiHLambda", "Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambda->Sumw2();
    fOutputList->Add(fDphiHLambda);

    fDphiHLambdaFilterbit = new THnSparseF("fDphiHLambdaFilterbit", "Hadron-Lambda Correlation Histogram (daughter has filter bit kTrkGlobalNoDCA) ", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaFilterbit->Sumw2();
    fOutputList->Add(fDphiHLambdaFilterbit);

    fDphiHLambdaEff = new THnSparseF("fDphiHLambdaEff", "Efficiency-corrected Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff->Sumw2();
    fOutputList->Add(fDphiHLambdaEff);

    fDphiHLambdaV0 = new THnSparseF("fDphiHLambdaV0", "Hadron-Lambda (using V0) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaV0->Sumw2();
    fOutputList->Add(fDphiHLambdaV0);

    fDphiHLambdaRotated = new THnSparseF("fDphiHLambdaRotated", "Hadron-Lambda (rotated) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaRotated->Sumw2();
    fOutputList->Add(fDphiHLambdaRotated);

    fDphiHLambdaRotatedPi = new THnSparseF("fDphiHLambdaRotatedPi", "Hadron-Lambda (rotated pi) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaRotatedPi->Sumw2();
    fOutputList->Add(fDphiHLambdaRotatedPi);

    fDphiHLambdaRotatedProton = new THnSparseF("fDphiHLambdaRotatedProton", "Hadron-Lambda (proton rotated) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaRotatedProton->Sumw2();
    fOutputList->Add(fDphiHLambdaRotatedProton);

    fDphiHLambdaFlipped = new THnSparseF("fDphiHLambdaFlipped", "Hadron-Lambda (flipped) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaFlipped->Sumw2();
    fOutputList->Add(fDphiHLambdaFlipped);

    fDphiHLambdaLS = new THnSparseF("fDphiHLambdaLS", "Hadron-Lambda LS Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaLS->Sumw2();
    fOutputList->Add(fDphiHLambdaLS);

    fDphiHLambdaMixed = new THnSparseF("fDphiHLambdaMixed", "Mixed Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed);

    fDphiHLambdaLSMixed = new THnSparseF("fDphiHLambdaLSMixed", "Mixed Hadron-Lambda LS Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaLSMixed->Sumw2();
    fOutputList->Add(fDphiHLambdaLSMixed);


    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {20, 20, 16, 20, 10};
    double hh_cor_mins[5] = {2, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 12, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHH->Sumw2();
    fOutputList->Add(fDphiHH);

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff->Sumw2();
    fOutputList->Add(fDphiHHEff);

    fDphiTriggerTrigger = new THnSparseF("fDphiTriggerTrigger", "Trigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTrigger->Sumw2();
    fOutputList->Add(fDphiTriggerTrigger);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed->Sumw2();
    fOutputList->Add(fDphiHHMixed);

    fDphiTriggerTriggerMixed = new THnSparseF("fDphiTriggerTriggerMixed", "MixedTrigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTriggerMixed->Sumw2();
    fOutputList->Add(fDphiTriggerTriggerMixed);


    //axes are pion 0: dca pion 1: dca proton
    //              2: pt pion  3: pt proton
    //              4: pt L     5: mass L
    int dca_bins[6] = {100, 100, 20, 20, 20, 50};
    double dca_mins[6] = {-2.4, -2.4, 0, 0, 0, 1.08};
    double dca_maxes[6] = {2.4, 2.4, 10, 10, 10, 1.16};

    fLambdaDaughterDCA = new THnSparseF("fLambdaDaughterDCA", "#Lambda^{0} daughter DCA dist", 6, dca_bins, dca_mins, dca_maxes);
    fLambdaDaughterDCA->Sumw2();
    fOutputList->Add(fLambdaDaughterDCA);

    fTofTest = new TH2D("fTofTest", "Beta vs P test hist", 1000, 0, 10, 230, 0, 2.3);
    fOutputList->Add(fTofTest);

    fAssociatedPtEventClass = new TH2D("fAssociatedPtEventClass", "Associated p_{T} dist in different event classes", 3, 0, 3, 500, 0, 10);
    fOutputList->Add(fAssociatedPtEventClass);

    fTriggerPtEventClass = new TH2D("fTriggerPtEventClass", "Trigger p_{T} dist in different event classes", 3, 0, 3, 500, 0, 10);
    fOutputList->Add(fTriggerPtEventClass);

    fLambdaPtEventClass = new TH2D("fLambdaPtEventClass", "Lambda p_{T} dist in different event classes", 3, 0, 3, 500, 0, 10);
    fOutputList->Add(fLambdaPtEventClass);

    fMultDistMinBias = new TH1D("fMultDistMinBias", "Multiplicty distribution (min bias, all events)",  100, 0, 100);
    fOutputList->Add(fMultDistMinBias);

    fMultDistHHEvent = new TH1D("fMultDistHHEvent", "Multiplicty distribution (events that h-h correlation performed", 100, 0, 100);
    fOutputList->Add(fMultDistHHEvent);

    fMultDistHLambdaEvent = new TH1D("fMultDistHLambdaEvent", "Multiplicty distribution (events that h-#Lambda correlation performed", 100, 0, 100);
    fOutputList->Add(fMultDistHLambdaEvent);

    PostData(1, fOutputList);
}


void AliAnalysisTaskLambdaHadronRatio::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff)
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

void AliAnalysisTaskLambdaHadronRatio::FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist)
{
    double dist_points[5]; //Pt, Phi, Eta, M, event multiplicity
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = particle.M();
        dist_points[4] = multPercentile;
        fDist->Fill(dist_points);
    }
}

AliAnalysisTaskLambdaHadronRatio::AliMotherContainer AliAnalysisTaskLambdaHadronRatio::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{
    AliAnalysisTaskLambdaHadronRatio::AliMotherContainer mom;

    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));

    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

AliAnalysisTaskLambdaHadronRatio::AliMotherContainer AliAnalysisTaskLambdaHadronRatio::RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle)
{
    AliAnalysisTaskLambdaHadronRatio::AliMotherContainer mom;
    // Rotating track1
    TVector3 track1Vector(track1->Px(), track1->Py(), track1->Pz());
    track1Vector.RotateZ(angle);
    mom.particle.SetPx(track1Vector(0) + track2->Px());
    mom.particle.SetPy(track1Vector(1) + track2->Py());
    mom.particle.SetPz(track1Vector(2) + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

void AliAnalysisTaskLambdaHadronRatio::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskLambdaHadronRatio::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskLambdaHadronRatio::SetAssociatedBit(float associatedBit) {
    fAssociatedBit = associatedBit;
}

void AliAnalysisTaskLambdaHadronRatio::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}

void AliAnalysisTaskLambdaHadronRatio::LoadEfficiencies(TString filePath) {
    TFile* effFile = TFile::Open(filePath);

    if(!effFile) {
        AliFatal("NULL INPUT FILE WHEN LOADING EFFICIENCIES, EXITING");
    }

    fLambdaEff = (TH1D*) effFile->Get("fLambdaEff")->Clone("fLambdaEffClone");
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

AliAnalysisTaskLambdaHadronRatio::AliMotherContainer AliAnalysisTaskLambdaHadronRatio::FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{
    AliAnalysisTaskLambdaHadronRatio::AliMotherContainer mom;
    // Flipping track1
    mom.particle.SetPx(-track1->Px() + track2->Px());
    mom.particle.SetPy(-track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

void AliAnalysisTaskLambdaHadronRatio::MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetID() == lambda.daughter1ID) || (trigger->GetID() == lambda.daughter2ID)) continue;

            dphi_point[1] = lambda.particle.Pt();
            dphi_point[2] = trigger->Phi() - lambda.particle.Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - lambda.particle.Eta();
            dphi_point[4] = lambda.particle.M();
            dphi_point[5] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (lambda.particle.Pt() < 10 && lambda.particle.Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int lambdaBin = fLambdaEff->FindBin(lambda.particle.Pt());
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

void AliAnalysisTaskLambdaHadronRatio::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
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

void AliAnalysisTaskLambdaHadronRatio::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = j+1; i < (int)trigger_list.size(); i++) {
            auto associate = trigger_list[i];

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
                int associatedBin = fTriggerEff->FindBin(associate->Pt());
                double associatedEff = fTriggerEff->GetBinContent(associatedBin);
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

void AliAnalysisTaskLambdaHadronRatio::MakeMixedHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff)
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

                dphi_point[1] = lambda.particle.Pt();
                dphi_point[2] = trigger->Phi() - lambda.particle.Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - lambda.particle.Eta();
                dphi_point[4] = lambda.particle.M();
                dphi_point[5] = zVtx;
                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (lambda.particle.Pt() < 10 && lambda.particle.Pt() > 0.5));
                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int lambdaBin = fLambdaEff->FindBin(lambda.particle.Pt());
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

void AliAnalysisTaskLambdaHadronRatio::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
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

bool AliAnalysisTaskLambdaHadronRatio::PassDaughterCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    return pass;
}

bool AliAnalysisTaskLambdaHadronRatio::PassAssociatedCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestFilterMask(fAssociatedBit);

    return pass;
}

Bool_t AliAnalysisTaskLambdaHadronRatio::PassTriggerCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestBit(fTriggerBit);

    return pass;
}

void AliAnalysisTaskLambdaHadronRatio::UserExec(Option_t*)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
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

    std::vector<AliAODTrack*> unlikelyProton_list;
    std::vector<AliAODTrack*> unlikelyPion_list;
    std::vector<AliAODTrack*> proton_list;
    std::vector<AliAODTrack*> antiProton_list;
    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> filterbit_proton_list;
    std::vector<AliAODTrack*> filterbit_antiProton_list;
    std::vector<AliAODTrack*> filterbit_piPlus_list;
    std::vector<AliAODTrack*> filterbit_piMinus_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> associated_h_list_2_4;
    std::vector<AliAODTrack*> all_hadron_list;
    std::vector<AliAODTrack*> k_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    // Bool to keep track if the event has a high-pt (> 4 GeV) trigger
    bool is_triggered_event = false;


    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //List for comparison with cuts/filter bits
        all_hadron_list.push_back(track);

        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);
            if(triggerPart->Pt() > 4) is_triggered_event = true;
        }

        if(PassAssociatedCuts(track)) {
            associated_h_list.push_back(track);
            if(track->Pt() < 4 && track->Pt() > 2) {
                associated_h_list_2_4.push_back(track);
            }

        }

        if(track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) {

            double filterbit_TPCNSigmaPion = 1000;
            double filterbit_TOFNSigmaPion = 1000;

            filterbit_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            filterbit_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

            if(TMath::Abs(filterbit_TPCNSigmaPion) <= 3 && (TMath::Abs(filterbit_TOFNSigmaPion) <= 3 || filterbit_TOFNSigmaPion == 1000)) {

                if(track->Charge() == 1){
                    filterbit_piPlus_list.push_back(track);
                }
                else {
                    filterbit_piMinus_list.push_back(track);
                }
            }

            double filterbit_TPCNSigmaProton = 1000;
            double filterbit_TOFNSigmaProton = 1000;


            filterbit_TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
            filterbit_TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(track, AliPID::kProton);

            if(TMath::Abs(filterbit_TPCNSigmaProton) <= 2 && (TMath::Abs(filterbit_TOFNSigmaProton) <= 2 || filterbit_TOFNSigmaProton == 1000)) {

                if(track->Charge() == 1){
                    filterbit_proton_list.push_back(track);
                }
                else {
                    filterbit_antiProton_list.push_back(track);
                }
            }
        } 
            
        if(PassDaughterCuts(track)) {
            double TPCNSigmaPion = 1000;
            double TOFNSigmaPion = 1000;

            TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

            double time = track->GetTOFsignal() - fpidResponse->GetTOFResponse().GetStartTime(track->P());
            double length = track->GetIntegratedLength();
            double v = length/time;
            double c = 0.0288782;
            double beta = v/c;

            fTofTest->Fill(track->P(), beta);

            if(TOFNSigmaPion != 1000 && track->Charge() != 1) unlikelyPion_list.push_back(track);

            if(TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000)) {

                if(track->Charge() == 1){
                    piPlus_list.push_back(track);
                }
                else {
                    piMinus_list.push_back(track);
                }
            }

            double TPCNSigmaProton = 1000;
            double TOFNSigmaProton = 1000;


            TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
            TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(track, AliPID::kProton);

            if(TOFNSigmaProton != 1000 && track->Charge() == 1) unlikelyProton_list.push_back(track);

            if(TMath::Abs(TPCNSigmaProton) <= 2 && (TMath::Abs(TOFNSigmaProton) <= 2 || TOFNSigmaProton == 1000)) {

                if(track->Charge() == 1){
                    proton_list.push_back(track);
                }
                else {
                    antiProton_list.push_back(track);
                }
            }
        }
    }

    //Making list of possible lambdas (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_filterbit_daughters;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_v0;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_signal_region;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_signal_region_2_4;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedPion;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedProton;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedPi;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_Flipped;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_LS;

    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) proton_list.size(); j++) {
            AliMotherContainer lambda = DaughtersToMother(piMinus_list[i], proton_list[j], 0.1396, 0.9383);
            auto pion = piMinus_list[i];
            auto proton = proton_list[j];
            
            double pion_dz[2];
            double pion_covar[3];

            double proton_dz[2];
            double proton_covar[3];

            bool is_pionDCA = piMinus_list[i]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., pion_dz, pion_covar);
            bool is_protonDCA = proton_list[j]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., proton_dz, proton_covar);

            if(is_pionDCA && is_protonDCA) {
                double fillArray[6] = {pion_dz[0], proton_dz[0], pion->Pt(), proton->Pt(), lambda.particle.Pt(), lambda.particle.M()};
                fLambdaDaughterDCA->Fill(fillArray);
            }

            AliMotherContainer lambda_RotatedPi = RotatedDaughtersToMother(piMinus_list[i], proton_list[j], 0.1396, 0.9383, TMath::Pi());
            AliMotherContainer lambda_Flipped = FlippedDaughtersToMother(piMinus_list[i], proton_list[j], 0.1396, 0.9383);
            lambda_list.push_back(lambda);
            lambda_list_RotatedPi.push_back(lambda_RotatedPi);
            lambda_list_Flipped.push_back(lambda_Flipped);

            AliMotherContainer lambda_Rotated;
            AliMotherContainer lambda_RotatedProton;
            for(int k = 1; k < 12; k++) {
                lambda_Rotated = RotatedDaughtersToMother(piMinus_list[i], proton_list[j], 0.1396, 0.9383, (2*TMath::Pi()*k)/12);
                lambda_RotatedProton = RotatedDaughtersToMother(proton_list[j], piMinus_list[i], 0.9383, 0.1396, (2*TMath::Pi()*k)/12);
                lambda_list_RotatedPion.push_back(lambda_Rotated);
                lambda_list_RotatedProton.push_back(lambda_RotatedProton);

            }
        }
    }

    for(int i = 0; i < (int)filterbit_piMinus_list.size(); i++) {
        for(int j = 0; j < (int) filterbit_proton_list.size(); j++) {
            AliMotherContainer filterbit_lambda = DaughtersToMother(filterbit_piMinus_list[i], filterbit_proton_list[j], 0.1396, 0.9383);
            lambda_list_filterbit_daughters.push_back(filterbit_lambda);
        }
    }

    for(int i = 0; i < (int)filterbit_piPlus_list.size(); i++) {
        for(int j = 0; j < (int) filterbit_antiProton_list.size(); j++) {
            AliMotherContainer filterbit_lambda = DaughtersToMother(filterbit_piPlus_list[i], filterbit_antiProton_list[j], 0.1396, 0.9383);
            lambda_list_filterbit_daughters.push_back(filterbit_lambda);
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) antiProton_list.size(); j++) {
            AliMotherContainer lambda = DaughtersToMother(piPlus_list[i], antiProton_list[j], 0.1396, 0.9383);
            AliMotherContainer lambda_RotatedPi = RotatedDaughtersToMother(piPlus_list[i], antiProton_list[j], 0.1396, 0.9383, TMath::Pi());
            AliMotherContainer lambda_Flipped = FlippedDaughtersToMother(piPlus_list[i], antiProton_list[j], 0.1396, 0.9383);
            lambda_list.push_back(lambda);
            lambda_list_RotatedPi.push_back(lambda_RotatedPi);
            lambda_list_Flipped.push_back(lambda_Flipped);

            AliMotherContainer lambda_Rotated;
            AliMotherContainer lambda_RotatedProton;
            for(int k = 1; k < 12; k++) {
                lambda_Rotated = RotatedDaughtersToMother(piPlus_list[i], antiProton_list[j], 0.1396, 0.9383, (2*TMath::Pi()*k)/12);
                lambda_RotatedProton = RotatedDaughtersToMother(antiProton_list[j], piPlus_list[i], 0.9383, 0.1396, (2*TMath::Pi()*k)/12);
                lambda_list_RotatedPion.push_back(lambda_Rotated);
                lambda_list_RotatedProton.push_back(lambda_RotatedProton);

            }

        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) proton_list.size(); j++) {
            if(piPlus_list[i]->GetID() == proton_list[j]->GetID()) continue;
            AliMotherContainer lambda = DaughtersToMother(piPlus_list[i], proton_list[j], 0.1396, 0.9383);
            lambda_list_LS.push_back(lambda);
        }
    }

    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) antiProton_list.size(); j++) {
            if(piMinus_list[i]->GetID() == antiProton_list[j]->GetID()) continue;
            AliMotherContainer lambda = DaughtersToMother(piMinus_list[i], antiProton_list[j], 0.1396, 0.9383);
            lambda_list_LS.push_back(lambda);
        }
    }


    for(int i = 0; i < (int)lambda_list.size(); i++) {
        if(lambda_list[i].particle.M() < 1.125 && lambda_list[i].particle.M() > 1.105) {
            lambda_list_signal_region.push_back(lambda_list[i]);
            if(lambda_list[i].particle.Pt() < 4 && lambda_list[i].particle.Pt() > 2) {
                lambda_list_signal_region_2_4.push_back(lambda_list[i]);
            }
        }
    }


    // V0 SECTION

    int numV0s = fAOD->GetNumberOfV0s();
    for(int i = 0; i < numV0s; i++) {
        AliAODv0 *v0 = fAOD->GetV0(i);
        if(v0->GetOnFlyStatus()) continue;

        AliAODTrack* posTrack = (AliAODTrack*) v0->GetDaughter(0);
        AliAODTrack* negTrack = (AliAODTrack*) v0->GetDaughter(1);

        // Occasionally returns null, not quite sure why...
        if(!posTrack || !negTrack) continue;
        if(!(PassDaughterCuts(posTrack) && PassDaughterCuts(negTrack))) continue;


        double TPCNSigmaProton = 1000;
        double TOFNSigmaProton = 1000;
        double TPCNSigmaPion = 1000;
        double TOFNSigmaPion = 1000;

        TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton);
        TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(posTrack, AliPID::kProton);
        TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion);
        TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(negTrack, AliPID::kPion);

        bool isNegTrackPion = TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000);
        bool isPosTrackProton = TMath::Abs(TPCNSigmaProton) <= 2 && (TMath::Abs(TOFNSigmaProton) <= 2 || TOFNSigmaProton == 1000);

        if(isNegTrackPion && isPosTrackProton) {
            auto lambda = DaughtersToMother(negTrack, posTrack, 0.1396, 0.9383);
            lambda_list_v0.push_back(lambda);
        }
    }


    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDist);
    FillSingleParticleDist(trigger_list, primZ, fTriggerDistEff, true);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);
    FillSingleParticleDist(all_hadron_list, primZ, fLooseDist);

    // Filling our single particle lambda distribution histogram:
    if(is_triggered_event) FillMotherDist(lambda_list, multPercentile, fTriggeredLambdaDist);
    if(is_triggered_event) FillMotherDist(lambda_list_filterbit_daughters, multPercentile, fTriggeredLambdaDistFilterbit);
    FillMotherDist(lambda_list, multPercentile, fLambdaDist);

    // Filling all of our correlation histograms
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambda, primZ, false);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_filterbit_daughters, fDphiHLambdaFilterbit, primZ, false);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambdaEff, primZ, true);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_v0, fDphiHLambdaV0, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedPi, fDphiHLambdaRotatedPi, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_Flipped, fDphiHLambdaFlipped, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedPion, fDphiHLambdaRotated, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedProton, fDphiHLambdaRotatedProton, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ, false);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, true);
    MakeSameTriggerTriggerCorrelations(trigger_list, fDphiTriggerTrigger, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_LS, fDphiHLambdaLS, primZ);

    fTriggersAndLambdasPerEvent_All->Fill(trigger_list.size(), lambda_list_signal_region.size());
    fTriggersAndLambdasPerEvent_2_4->Fill(trigger_list.size(), lambda_list_signal_region_2_4.size());

    if(is_triggered_event) {
        for(auto part : trigger_list) {
            fTriggerPtEventClass->Fill(0.1, part->Pt());
        }
        for(auto part : associated_h_list) {
            fAssociatedPtEventClass->Fill(0.1, part->Pt());
        }
        for(auto part : lambda_list_signal_region) {
            fLambdaPtEventClass->Fill(0.1, part.particle.Pt());
        }
        if(associated_h_list_2_4.size()) {
            for(auto part : trigger_list) {
                fTriggerPtEventClass->Fill(1.1, part->Pt());
            }
            for(auto part : associated_h_list) {
                fAssociatedPtEventClass->Fill(1.1, part->Pt());
            }
            for(auto part : lambda_list_signal_region) {
                fLambdaPtEventClass->Fill(1.1, part.particle.Pt());
            }
            if(lambda_list_signal_region_2_4.size()) {
                for(auto part : trigger_list) {
                    fTriggerPtEventClass->Fill(2.1, part->Pt());
                }
                for(auto part : associated_h_list) {
                    fAssociatedPtEventClass->Fill(2.1, part->Pt());
                }
                for(auto part : lambda_list_signal_region) {
                    fLambdaPtEventClass->Fill(2.1, part.particle.Pt());
                }
            }
        }
    }

    if(lambda_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHLambdaCorrelations(fCorPool, lambda_list, fDphiHLambdaMixed, primZ);
                MakeMixedHLambdaCorrelations(fCorPool, lambda_list_LS, fDphiHLambdaLSMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, trigger_list, fDphiTriggerTriggerMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronRatio::Terminate(Option_t *option)
{
}
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

#include "AliAnalysisTaskLambdaHadronResClosure.h"

ClassImp(AliAnalysisTaskLambdaHadronResClosure);

AliAnalysisTaskLambdaHadronResClosure::AliAnalysisTaskLambdaHadronResClosure() :

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
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_MCKin(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_checkMC(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_MCKin(0x0),
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
    fDphiHH_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHHMixed_MC(0x0)
{
}

AliAnalysisTaskLambdaHadronResClosure::AliAnalysisTaskLambdaHadronResClosure(const char *name) :
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
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_MCKin(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_checkMC(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_MCKin(0x0),
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
    fDphiHH_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHHMixed_MC(0x0)

{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronResClosure::~AliAnalysisTaskLambdaHadronResClosure()
{
    if(fOutputList) delete fOutputList;
    if(fTriggerEff) delete fTriggerEff;
    if(fAssociatedEff) delete fAssociatedEff;
    if(fLambdaEff) delete fLambdaEff;
}

void AliAnalysisTaskLambdaHadronResClosure::UserCreateOutputObjects()
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
    int mother_dist_bins[5] = {100, 16, 20, 140, 10};
    double mother_dist_mins[5] = {0, -3.14, -1, 1.06, 0};
    double mother_dist_maxes[5] = {15, 3.14, 1, 1.2, 100};

    fTriggeredLambdaDist = new THnSparseF("fTriggeredLambdaDist", "Lambda Distribution (with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDist->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist);

    fLambdaDist = new THnSparseF("fLambdaDist", "Lambda Distribution (reco with v0 finder)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fLambdaDist->Sumw2();
    fOutputList->Add(fLambdaDist);

    fLambdaDist_MC = new THnSparseF("fLambdaDist_MC", "Lambda Distribution (MC truth)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fLambdaDist_MC->Sumw2();
    fOutputList->Add(fLambdaDist_MC);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hl_cor_bins[6] = {9, 10, 16, 20, 140, 10};
    double hl_cor_mins[6] = {1.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hl_cor_maxes[6] = {10.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 1.2, 10};

    fDphiHLambdaEff = new THnSparseF("fDphiHLambdaEff", "Efficiency-corrected Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff->Sumw2();
    fOutputList->Add(fDphiHLambdaEff);

    fDphiHLambdaEff_MCKin = new THnSparseF("fDphiHLambdaEff_MCKin", "Efficiency-corrected Hadron-Lambda Correlation Histogram (using MC kinematics on lambda)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff_MCKin->Sumw2();
    fOutputList->Add(fDphiHLambdaEff_MCKin);

    fDphiHLambda_MC = new THnSparseF("fDphiHLambda_MC", "Hadron-Lambda Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambda_MC->Sumw2();
    fOutputList->Add(fDphiHLambda_MC);

    fDphiHLambdaMixed = new THnSparseF("fDphiHLambdaMixed", "Mixed Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed);

    fDphiHLambdaMixed_MCKin = new THnSparseF("fDphiHLambdaMixed_MCKin", "Mixed_MCKin Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MCKin->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MCKin);

    fDphiHLambdaMixed_MC = new THnSparseF("fDphiHLambdaMixed_MC", "Mixed Hadron-Lambda Correlation Histogram (MC truth)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_MC->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MC);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {9, 10, 16, 20, 10};
    double hh_cor_mins[5] = {1, 1, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {10, 6, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff->Sumw2();
    fOutputList->Add(fDphiHHEff);

    fDphiHHEff_checkMC = new THnSparseF("fDphiHHEff_checkMC", "Efficiency corrected Hadron-Hadron Correlation Histogram (trig, assoc physical MC prim)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff_checkMC->Sumw2();
    fOutputList->Add(fDphiHHEff_checkMC);

    fDphiHH_MC = new THnSparseF("fDphiHH_MC", "Hadron-Hadron Correlation Histogram (MC truth)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHH_MC->Sumw2();
    fOutputList->Add(fDphiHH_MC);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed->Sumw2();
    fOutputList->Add(fDphiHHMixed);

    fDphiHHMixed_MC = new THnSparseF("fDphiHHMixed_MC", "Mixed Hadron-Hadron Correlation Histogram (MC truth)", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed_MC->Sumw2();
    fOutputList->Add(fDphiHHMixed_MC);

    PostData(1, fOutputList);

}

AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer AliAnalysisTaskLambdaHadronResClosure::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{
    AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer mom;

    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));

    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

void AliAnalysisTaskLambdaHadronResClosure::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff)
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

void AliAnalysisTaskLambdaHadronResClosure::FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist)
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

void AliAnalysisTaskLambdaHadronResClosure::FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda, bool lambdaEff)
{
    double dist_points[5]; //Pt, Phi, Eta, M, event multiplicity
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = particle.M();
        dist_points[4] = multPercentile;
        bool in_pt_range = (particle.Pt() < 10 && particle.Pt() > 0.5);
        if(lambdaEff && in_pt_range) {
            int lambdaBin = fLambdaEff->FindBin(particle.Pt());
            double lambdaEff = fLambdaEff->GetBinContent(lambdaBin);
            double lambdaScale = 1.0/lambdaEff;
            fDist->Fill(dist_points, lambdaScale);
        }
        else{
            fDist->Fill(dist_points);
        }
    }
}

void AliAnalysisTaskLambdaHadronResClosure::FillMCMotherDist(std::vector<AliAODMCParticle*> particle_list, float multPercentile, THnSparse* fDist)
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

void AliAnalysisTaskLambdaHadronResClosure::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskLambdaHadronResClosure::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskLambdaHadronResClosure::SetAssociatedBit(float associatedBit) {
    fAssociatedBit = associatedBit;
}

void AliAnalysisTaskLambdaHadronResClosure::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}

void AliAnalysisTaskLambdaHadronResClosure::LoadEfficiencies(TString filePath) {
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

void AliAnalysisTaskLambdaHadronResClosure::MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {
            auto lambda = lambda_list[i];
            auto lambda_mc = (AliAODMCParticle*)fMCArray->At(lambda.motherLabel);

            //Make sure trigger isn't one of the daughters of lambda
            if((trigger->GetID() == lambda.daughter1ID) || (trigger->GetID() == lambda.daughter2ID)) {
                continue;
            }

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
void AliAnalysisTaskLambdaHadronResClosure::MakeSameHLambdaCorrelations_withMCKin(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)lambda_list.size(); i++) {

            auto lambda = lambda_list[i];
            AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(lambda.motherLabel);

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

void AliAnalysisTaskLambdaHadronResClosure::MakeSameMCHLambdaCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> lambda_list, THnSparse* fDphi, double zVtx)
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

            //Make sure trigger isn't one of the daughters of lambda, this is not necessarry if trigger is physical primary (oh well)
            if((trigger->GetLabel() == first_daughter_index) || (trigger->GetLabel() == second_daughter_index)) {
                continue;
            }

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


void AliAnalysisTaskLambdaHadronResClosure::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
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
void AliAnalysisTaskLambdaHadronResClosure::MakeSameMCHHCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx)
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

void AliAnalysisTaskLambdaHadronResClosure::MakeMixedHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
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

void AliAnalysisTaskLambdaHadronResClosure::MakeMixedHLambdaCorrelations_withMCKin(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff)
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
                AliAODMCParticle* mcmother = (AliAODMCParticle*)fMCArray->At(lambda.motherLabel);

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
}

void AliAnalysisTaskLambdaHadronResClosure::MakeMixedMCHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> lambda_list , THnSparse* fDphi, double zVtx)
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

void AliAnalysisTaskLambdaHadronResClosure::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
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

void AliAnalysisTaskLambdaHadronResClosure::MakeMixedMCHHCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> associated_h_list, THnSparse* fDphi, double zVtx)
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

bool AliAnalysisTaskLambdaHadronResClosure::PassDaughterCuts(AliAODTrack *track){

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

bool AliAnalysisTaskLambdaHadronResClosure::PassAssociatedCuts(AliAODTrack *track, bool checkMC){

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

bool AliAnalysisTaskLambdaHadronResClosure::PassTriggerCuts(AliAODTrack *track, bool checkMC){

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

bool AliAnalysisTaskLambdaHadronResClosure::IsMCChargedHadron(int pdg_code) {
    // checks to see if pdg code matches pion, proton, kaon, electron or muon
    if((TMath::Abs(pdg_code) == 321)
        || (TMath::Abs(pdg_code) == 211)
        || (TMath::Abs(pdg_code) == 2212)
        || (TMath::Abs(pdg_code) == 11)
        || (TMath::Abs(pdg_code) == 13)) return true;

    return false;
}

bool AliAnalysisTaskLambdaHadronResClosure::PassMCTriggerCuts(AliAODMCParticle *mc_particle){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false; // for now trigger is physical primary, could change
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronResClosure::PassMCAssociatedCuts(AliAODMCParticle *mc_particle){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false;
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronResClosure::PassMCLambdaCuts(AliAODMCParticle *mc_particle){

    if(!(TMath::Abs(mc_particle->GetPdgCode()) == 3122)) return false;
    // if(!(mc_particle->IsPhysicalPrimary())) return false; // testing, testing, 1 2 3
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
        

void AliAnalysisTaskLambdaHadronResClosure::UserExec(Option_t*)
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
    
    // The following lists contain tracks that pass daughter cuts
    std::vector<AliAODTrack*> piminus_list;
    std::vector<AliAODTrack*> piplus_list;
    std::vector<AliAODTrack*> proton_list;
    std::vector<AliAODTrack*> antiproton_list;

    std::vector<AliAODTrack*> trigger_list_checkMC;
    std::vector<AliAODTrack*> associated_h_list_checkMC;

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
            if(mc_label < 0) continue;
            AliAODMCParticle* mc_part = (AliAODMCParticle*)fMCArray->At(mc_label);
            if(mc_part->GetPdgCode() == 211) {
                piplus_list.push_back(track);
            }
            if(mc_part->GetPdgCode() == -211) {
                piminus_list.push_back(track);
            }
            if(mc_part->GetPdgCode() == 2212) {
                proton_list.push_back(track);
            }
            if(mc_part->GetPdgCode() == -2212) {
                antiproton_list.push_back(track);
            }
        }
    }

    //Making list of possible lambdas (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> antilambda_list;
    std::vector<AliAnalysisTaskLambdaHadronResClosure::AliMotherContainer> lambda_list;

    // For now we guarantee the p-pi pair came from an actual lambda, and that the mothers of each daughter are the same
    for(int i = 0; i < (int)proton_list.size(); i++) {
        for(int j = 0; j < (int) piminus_list.size(); j++) {
            auto proton = proton_list[i];
            auto piminus = piminus_list[j];
            if(proton->GetID() == piminus->GetID()) continue;
            auto proton_mc = (AliAODMCParticle*)fMCArray->At(proton->GetLabel());
            auto piminus_mc = (AliAODMCParticle*)fMCArray->At(piminus->GetLabel());
            int mlabel_pos = proton_mc->GetMother();
            int mlabel_neg = piminus_mc->GetMother();
            if(mlabel_pos < 0 || mlabel_neg < 0) continue;
            if(mlabel_pos != mlabel_neg) continue;
            auto mother = (AliAODMCParticle*)fMCArray->At(mlabel_pos);
            if(mother->GetPdgCode() != 3122) continue;
            AliMotherContainer lambda = DaughtersToMother(proton, piminus, 0.9383, 0.1396);
            lambda.motherLabel = mlabel_pos;
            lambda_list.push_back(lambda);
        }
    }

    for(int i = 0; i < (int)antiproton_list.size(); i++) {
        for(int j = 0; j < (int) piplus_list.size(); j++) {
            auto antiproton = antiproton_list[i];
            auto piplus = piplus_list[j];
            if(antiproton->GetID() == piplus->GetID()) continue;
            auto antiproton_mc = (AliAODMCParticle*)fMCArray->At(antiproton->GetLabel());
            auto piplus_mc = (AliAODMCParticle*)fMCArray->At(piplus->GetLabel());
            int mlabel_pos = antiproton_mc->GetMother();
            int mlabel_neg = piplus_mc->GetMother();
            if(mlabel_pos < 0 || mlabel_neg < 0) continue;
            if(mlabel_pos != mlabel_neg) continue;
            auto mother = (AliAODMCParticle*)fMCArray->At(mlabel_pos);
            if(mother->GetPdgCode() != -3122) continue;
            AliMotherContainer antilambda = DaughtersToMother(antiproton, piplus, 0.9383, 0.1396);
            antilambda.motherLabel = mlabel_pos;
            antilambda_list.push_back(antilambda);
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

    // Filling our single particle lambda distribution histogram:
    if(is_triggered_event) FillMotherDist(lambda_list, multPercentile, fTriggeredLambdaDist, false);
    if(is_triggered_event) FillMotherDist(antilambda_list, multPercentile, fTriggeredLambdaDist, true);

    MakeSameHLambdaCorrelations(trigger_list, antilambda_list, fDphiHLambdaEff, primZ, true, true);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambdaEff, primZ, true, false);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list, antilambda_list, fDphiHLambdaEff_MCKin, primZ, true);
    MakeSameHLambdaCorrelations_withMCKin(trigger_list, lambda_list, fDphiHLambdaEff_MCKin, primZ, true);

    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, true);
    MakeSameHHCorrelations(trigger_list_checkMC, associated_h_list_checkMC, fDphiHHEff_checkMC, primZ, true);

    if(/*lambda_list.size() > 0 && */ associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fCorPool->IsReady()) {
                MakeMixedHLambdaCorrelations(fCorPool, antilambda_list, fDphiHLambdaMixed, primZ, true, true);
                MakeMixedHLambdaCorrelations(fCorPool, lambda_list, fDphiHLambdaMixed, primZ, true, false);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, antilambda_list, fDphiHLambdaMixed_MCKin, primZ, true);
                MakeMixedHLambdaCorrelations_withMCKin(fCorPool, lambda_list, fDphiHLambdaMixed_MCKin, primZ, true);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }
    


    std::vector<AliAODMCParticle*> real_trigger_list;
    std::vector<AliAODMCParticle*> real_associated_list;
    std::vector<AliAODMCParticle*> real_lambda_list;

    for(int mc_index = 0; mc_index < fMCArray->GetEntries(); mc_index++){

        AliAODMCParticle *mc_particle = (AliAODMCParticle*)fMCArray->At(mc_index);

        if(PassMCTriggerCuts(mc_particle)) {
            real_trigger_list.push_back(mc_particle);
            AliCFParticle *trigger_particle = new AliCFParticle(mc_particle->Pt(), mc_particle->Eta(), mc_particle->Phi(), mc_particle->Charge(), 0);
            fMixedMCTrackObjArray->Add(trigger_particle);
        }
        if(PassMCAssociatedCuts(mc_particle)) real_associated_list.push_back(mc_particle);
        if(PassMCLambdaCuts(mc_particle)) real_lambda_list.push_back(mc_particle);

    }

    FillSingleMCParticleDist(real_trigger_list, primZ, fTriggerDist_MC);
    FillSingleMCParticleDist(real_associated_list, primZ, fAssociatedDist_MC);
    FillMCMotherDist(real_lambda_list, primZ, fLambdaDist_MC);

    MakeSameMCHLambdaCorrelations(real_trigger_list, real_lambda_list, fDphiHLambda_MC, primZ);
    MakeSameMCHHCorrelations(real_trigger_list, real_associated_list, fDphiHH_MC, primZ);
    


    if(real_associated_list.size() > 0 /*&& real_lambda_list.size() > 0*/) {
        AliEventPool *fMCCorPool = fMCCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fMCCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fMCCorPool->IsReady()) {
                MakeMixedMCHLambdaCorrelations(fMCCorPool, real_lambda_list, fDphiHLambdaMixed_MC, primZ);
                MakeMixedMCHHCorrelations(fMCCorPool, real_associated_list, fDphiHHMixed_MC, primZ);
            }
            if(fMixedMCTrackObjArray->GetEntries() > 0) {
                fMCCorPool->UpdatePool(fMixedMCTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronResClosure::Terminate(Option_t *option)
{
}
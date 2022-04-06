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
    fCorPoolMgr_highestPt(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_highestPt(0x0),
    fAssociatedHDist(0x0),
    fTriggeredLambdaDist(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_highestPt(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_highestPt(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_highestPt(0x0),
    fDphiHHMixed(0x0),
    fDphiHHMixed_highestPt(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0)
{
}

AliAnalysisTaskLambdaHadronRatio::AliAnalysisTaskLambdaHadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fCorPoolMgr_highestPt(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fLambdaEff(0x0),
    fTriggerDistEff(0x0),
    fTriggerDistEff_highestPt(0x0),
    fAssociatedHDist(0x0),
    fTriggeredLambdaDist(0x0),
    fDphiHLambdaEff(0x0),
    fDphiHLambdaEff_highestPt(0x0),
    fDphiHHEff(0x0),
    fDphiHHEff_highestPt(0x0),
    fDphiHLambdaMixed(0x0),
    fDphiHLambdaMixed_highestPt(0x0),
    fDphiHHMixed(0x0),
    fDphiHHMixed_highestPt(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0)

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

    fCorPoolMgr_highestPt = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr_highestPt->SetTargetValues(trackDepth, 0.1, 5);

    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {200, 16, 20, 10};
    double dist_mins[4] = {0.0, 0, -1, -10};
    double dist_maxes[4] = {20.0, 6.28, 1, 10};

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

    fTriggeredLambdaDist = new THnSparseF("fTriggeredLambdaDist", "Lambda Distribution (with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredLambdaDist->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hl_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    double hl_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hl_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10};

    fDphiHLambdaEff = new THnSparseF("fDphiHLambdaEff", "Efficiency-corrected Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff->Sumw2();
    fOutputList->Add(fDphiHLambdaEff);

    fDphiHLambdaEff_highestPt = new THnSparseF("fDphiHLambdaEff_highestPt", "Efficiency-corrected Hadron-Lambda Correlation Histogram (highest pt trigger)", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaEff_highestPt->Sumw2();
    fOutputList->Add(fDphiHLambdaEff_highestPt);

    fDphiHLambdaMixed = new THnSparseF("fDphiHLambdaMixed", "Mixed Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed);

    fDphiHLambdaMixed_highestPt = new THnSparseF("fDphiHLambdaMixed_highestPt", "Mixed Hadron-Lambda Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHLambdaMixed_highestPt->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_highestPt);


    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {16, 10, 16, 20, 10};
    double hh_cor_mins[5] = {4, 1, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 6, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff->Sumw2();
    fOutputList->Add(fDphiHHEff);

    fDphiHHEff_highestPt = new THnSparseF("fDphiHHEff_highestPt", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff_highestPt->Sumw2();
    fOutputList->Add(fDphiHHEff_highestPt);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed->Sumw2();
    fOutputList->Add(fDphiHHMixed);

    fDphiHHMixed_highestPt = new THnSparseF("fDphiHHMixed_highestPt", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed_highestPt->Sumw2();
    fOutputList->Add(fDphiHHMixed_highestPt);

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

void AliAnalysisTaskLambdaHadronRatio::FillMotherDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist, bool isAntiLambda)
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
        fDist->Fill(dist_points);
    }
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

void AliAnalysisTaskLambdaHadronRatio::MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
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

void AliAnalysisTaskLambdaHadronRatio::MakeMixedHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx, bool eff, bool isAntiLambda)
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

    if(track->GetID() < 0) return false;

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

    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> trigger_list_highestPt; // not actually a list but too lazy to rewrite correlation function
    std::vector<AliAODTrack*> associated_h_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    TObjArray* fMixedTrackObjArray_highestPt = new TObjArray;
    fMixedTrackObjArray_highestPt->SetOwner(kTRUE);

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
    }

    if(maxTrigger){
        trigger_list_highestPt.push_back(maxTrigger); 
        AliCFParticle *triggerPart_highestPt = new AliCFParticle(maxTrigger->Pt(), maxTrigger->Eta(), maxTrigger->Phi(), maxTrigger->Charge(), 0);
        fMixedTrackObjArray_highestPt->Add(triggerPart_highestPt);
    }

    //Making list of possible lambdas (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> antilambda_list;

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


        double pos_TPCNSigmaProton = 1000;
        double pos_TOFNSigmaProton = 1000;
        double neg_TPCNSigmaPion = 1000;
        double neg_TOFNSigmaPion = 1000;

        double neg_TPCNSigmaProton = 1000;
        double neg_TOFNSigmaProton = 1000;
        double pos_TPCNSigmaPion = 1000;
        double pos_TOFNSigmaPion = 1000;

        pos_TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton);
        pos_TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(posTrack, AliPID::kProton);
        neg_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion);
        neg_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(negTrack, AliPID::kPion);

        neg_TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(negTrack, AliPID::kProton);
        neg_TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(negTrack, AliPID::kProton);
        pos_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion);
        pos_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(posTrack, AliPID::kPion);

        bool isNegTrackPion = TMath::Abs(neg_TPCNSigmaPion) <= 3 && (TMath::Abs(neg_TOFNSigmaPion) <= 3 || neg_TOFNSigmaPion == 1000);
        bool isPosTrackProton = TMath::Abs(pos_TPCNSigmaProton) <= 2 && (TMath::Abs(pos_TOFNSigmaProton) <= 2 || pos_TOFNSigmaProton == 1000);

        bool isPosTrackPion = TMath::Abs(pos_TPCNSigmaPion) <= 3 && (TMath::Abs(pos_TOFNSigmaPion) <= 3 || pos_TOFNSigmaPion == 1000);
        bool isNegTrackProton = TMath::Abs(neg_TPCNSigmaProton) <= 2 && (TMath::Abs(neg_TOFNSigmaProton) <= 2 || neg_TOFNSigmaProton == 1000);

        if((isNegTrackPion && isPosTrackProton)) {
            AliMotherContainer lambda;
            lambda.vzero = v0;
            lambda.daughter1ID = posTrack->GetID();
            lambda.daughter2ID = negTrack->GetID();
            lambda_list.push_back(lambda);
        }

        if((isPosTrackPion && isNegTrackProton)) {
            AliMotherContainer antilambda;
            antilambda.vzero = v0;
            antilambda.daughter1ID = posTrack->GetID();
            antilambda.daughter2ID = negTrack->GetID();
            antilambda_list.push_back(antilambda);
        }

    }


    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDistEff, true);
    FillSingleParticleDist(trigger_list_highestPt, primZ, fTriggerDistEff_highestPt, true);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);

    // Filling our single particle lambda distribution histogram:
    if(is_triggered_event) FillMotherDist(lambda_list, multPercentile, fTriggeredLambdaDist, false);
    if(is_triggered_event) FillMotherDist(antilambda_list, multPercentile, fTriggeredLambdaDist, true);

    MakeSameHLambdaCorrelations(trigger_list, antilambda_list, fDphiHLambdaEff, primZ, true, true);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambdaEff, primZ, true, false);

    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, true);

    // Highest pt trigger correlations
    MakeSameHLambdaCorrelations(trigger_list_highestPt, antilambda_list, fDphiHLambdaEff_highestPt, primZ, true, true);
    MakeSameHLambdaCorrelations(trigger_list_highestPt, lambda_list, fDphiHLambdaEff_highestPt, primZ, true, false);
    MakeSameHHCorrelations(trigger_list_highestPt, associated_h_list, fDphiHHEff_highestPt, primZ, true);

    if(lambda_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        AliEventPool *fCorPool_highestPt = fCorPoolMgr_highestPt->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fCorPool->IsReady()) {
                MakeMixedHLambdaCorrelations(fCorPool, antilambda_list, fDphiHLambdaMixed, primZ, true, true);
                MakeMixedHLambdaCorrelations(fCorPool, lambda_list, fDphiHLambdaMixed, primZ, true, false);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
            }
            if(fCorPool_highestPt->IsReady()) {
                MakeMixedHLambdaCorrelations(fCorPool_highestPt, antilambda_list, fDphiHLambdaMixed_highestPt, primZ, true, true);
                MakeMixedHLambdaCorrelations(fCorPool_highestPt, lambda_list, fDphiHLambdaMixed_highestPt, primZ, true, false);
                MakeMixedHHCorrelations(fCorPool_highestPt, associated_h_list, fDphiHHMixed_highestPt, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
            if(fMixedTrackObjArray_highestPt->GetEntries() > 0) {
                fCorPool_highestPt->UpdatePool(fMixedTrackObjArray_highestPt);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronRatio::Terminate(Option_t *option)
{
}

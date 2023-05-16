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
#include "AliEventPoolManager.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCFParticle.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliAODv0.h"

#include "AliAnalysisTaskLambdaHadronMC.h"

ClassImp(AliAnalysisTaskLambdaHadronMC);

AliAnalysisTaskLambdaHadronMC::AliAnalysisTaskLambdaHadronMC() :
    AliAnalysisTaskSE(),
    fMCEvent(0x0),
    fOutputList(0x0),
    fMCCorPoolMgr(0x0),
    fTriggerDist_MC(0x0),
    fAssociatedDist_MC(0x0),
    fLambdaDist_MC(0x0),
    fPhiDist_MC(0x0),
    fTriggeredTriggerDist_MC(0x0),
    fTriggeredAssociatedDist_MC(0x0),
    fTriggeredLambdaDist_MC(0x0),
    fTriggeredPhiDist_MC(0x0),
    fDphiHH_MC(0x0),
    fDphiHLambda_MC(0x0),
    fDphiHPhi_MC(0x0),
    fDphiHHMixed_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHPhiMixed_MC(0x0)
{
}

AliAnalysisTaskLambdaHadronMC::AliAnalysisTaskLambdaHadronMC(const char *name) :
    AliAnalysisTaskSE(name),
    fMCEvent(0x0),
    fOutputList(0x0),
    fMCCorPoolMgr(0x0),
    fTriggerDist_MC(0x0),
    fAssociatedDist_MC(0x0),
    fLambdaDist_MC(0x0),
    fPhiDist_MC(0x0),
    fTriggeredTriggerDist_MC(0x0),
    fTriggeredAssociatedDist_MC(0x0),
    fTriggeredLambdaDist_MC(0x0),
    fTriggeredPhiDist_MC(0x0),
    fDphiHH_MC(0x0),
    fDphiHLambda_MC(0x0),
    fDphiHPhi_MC(0x0),
    fDphiHHMixed_MC(0x0),
    fDphiHLambdaMixed_MC(0x0),
    fDphiHPhiMixed_MC(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronMC::~AliAnalysisTaskLambdaHadronMC()
{
    if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskLambdaHadronMC::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {0, 100.0};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fMCCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fMCCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);


    fMultDist = new TH1D("fMultDist", "Charged part. in V0A acceptance", 1000, 0, 1000);
    fOutputList->Add(fMultDist);

    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {100, 16, 100, 10};
    double dist_mins[4] = {0.0, 0, -10, -10};
    double dist_maxes[4] = {10.0, 6.28, 10, 10};

    fTriggerDist_MC = new THnSparseF("fTriggerDist_MC", "Trigger Hadron Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist_MC->Sumw2();
    fOutputList->Add(fTriggerDist_MC);

    fAssociatedDist_MC = new THnSparseF("fAssociatedDist_MC", "Associated Hadron Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedDist_MC->Sumw2();
    fOutputList->Add(fAssociatedDist_MC);

    fLambdaDist_MC = new THnSparseF("fLambdaDist_MC", "Lambda Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fLambdaDist_MC->Sumw2();
    fOutputList->Add(fLambdaDist_MC);

    fPhiDist_MC = new THnSparseF("fPhiDist_MC", "Phi(1020) Distribution (MC truth)", 4, dist_bins, dist_mins, dist_maxes);
    fPhiDist_MC->Sumw2();
    fOutputList->Add(fPhiDist_MC);

    fTriggerDist_MC_no_eta_cut = new THnSparseF("fTriggerDist_MC_no_eta_cut", "Trigger Hadron Distribution (MC_no_eta_cut truth)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist_MC_no_eta_cut->Sumw2();
    fOutputList->Add(fTriggerDist_MC_no_eta_cut);

    fAssociatedDist_MC_no_eta_cut = new THnSparseF("fAssociatedDist_MC_no_eta_cut", "Associated Hadron Distribution (MC_no_eta_cut truth)", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedDist_MC_no_eta_cut->Sumw2();
    fOutputList->Add(fAssociatedDist_MC_no_eta_cut);

    fLambdaDist_MC_no_eta_cut = new THnSparseF("fLambdaDist_MC_no_eta_cut", "Lambda Distribution (MC_no_eta_cut truth)", 4, dist_bins, dist_mins, dist_maxes);
    fLambdaDist_MC_no_eta_cut->Sumw2();
    fOutputList->Add(fLambdaDist_MC_no_eta_cut);

    fPhiDist_MC_no_eta_cut = new THnSparseF("fPhiDist_MC_no_eta_cut", "Phi(1020) Distribution (MC_no_eta_cut truth)", 4, dist_bins, dist_mins, dist_maxes);
    fPhiDist_MC_no_eta_cut->Sumw2();
    fOutputList->Add(fPhiDist_MC_no_eta_cut);

    fTriggeredTriggerDist_MC = new THnSparseF("fTriggeredTriggerDist_MC", "Trigger Distribution (MC truth, with triggered event)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggeredTriggerDist_MC->Sumw2();
    fOutputList->Add(fTriggeredTriggerDist_MC);

    fTriggeredAssociatedDist_MC = new THnSparseF("fTriggeredAssociatedDist_MC", "Associated Distribution (MC truth, with triggered event)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggeredAssociatedDist_MC->Sumw2();
    fOutputList->Add(fTriggeredAssociatedDist_MC);

    fTriggeredLambdaDist_MC = new THnSparseF("fTriggeredLambdaDist_MC", "Lambda Distribution (MC truth, with triggered event)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggeredLambdaDist_MC->Sumw2();
    fOutputList->Add(fTriggeredLambdaDist_MC);

    fTriggeredPhiDist_MC = new THnSparseF("fTriggeredPhiDist_MC", "Phi(1020) Distribution (MC truth, with triggered event)", 4, dist_bins, dist_mins, dist_maxes);
    fTriggeredPhiDist_MC->Sumw2();
    fOutputList->Add(fTriggeredPhiDist_MC);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int cor_bins[5] = {18, 16, 16, 20, 10};
    double cor_mins[5] = {1, 0, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double cor_maxes[5] = {10, 4, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHH_MC = new THnSparseF("fDphiHH_MC", "Hadron-Hadron Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHH_MC->Sumw2();
    fOutputList->Add(fDphiHH_MC);

    fDphiHLambda_MC = new THnSparseF("fDphiHLambda_MC", "Hadron-Lambda Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHLambda_MC->Sumw2();
    fOutputList->Add(fDphiHLambda_MC);

    fDphiHPhi_MC = new THnSparseF("fDphiHPhi_MC", "Hadron-Phi Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHPhi_MC->Sumw2();
    fOutputList->Add(fDphiHPhi_MC);

    fDphiHHMixed_MC = new THnSparseF("fDphiHHMixed_MC", "Mixed Hadron-H Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHHMixed_MC->Sumw2();
    fOutputList->Add(fDphiHHMixed_MC);

    fDphiHLambdaMixed_MC = new THnSparseF("fDphiHLambdaMixed_MC", "Mixed Hadron-Lambda Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHLambdaMixed_MC->Sumw2();
    fOutputList->Add(fDphiHLambdaMixed_MC);

    fDphiHPhiMixed_MC = new THnSparseF("fDphiHPhiMixed_MC", "Mixed Hadron-Phi(1020) Correlation Histogram (MC truth)", 5, cor_bins, cor_mins, cor_maxes);
    fDphiHPhiMixed_MC->Sumw2();
    fOutputList->Add(fDphiHPhiMixed_MC);


    PostData(1, fOutputList);

}

void AliAnalysisTaskLambdaHadronMC::FillSingleMCParticleDist(std::vector<AliAODMCParticle*> particle_list, double zVtx, THnSparse* fDist)
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

void AliAnalysisTaskLambdaHadronMC::MakeSameMCCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> associated_list, THnSparse* fDphi, double zVtx)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_list.size(); i++) {
            auto associate = associated_list[i];

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

void AliAnalysisTaskLambdaHadronMC::MakeMixedMCCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> associated_list, THnSparse* fDphi, double zVtx)
{
    double dphi_point[5];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_list.size(); j++) {
                auto associate = associated_list[j];

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

bool AliAnalysisTaskLambdaHadronMC::IsMCChargedHadron(int pdg_code) {
    // checks to see if pdg code matches pion, proton, kaon, electron or muon
    if((TMath::Abs(pdg_code) == 321)
        || (TMath::Abs(pdg_code) == 211)
        || (TMath::Abs(pdg_code) == 2212)
        || (TMath::Abs(pdg_code) == 11)
        || (TMath::Abs(pdg_code) == 13)) return true;

    return false;
}

bool AliAnalysisTaskLambdaHadronMC::PassMCTriggerCuts(AliAODMCParticle *mc_particle, bool etaCut){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false; 
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8) && etaCut) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronMC::PassMCAssociatedCuts(AliAODMCParticle *mc_particle, bool etaCut){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false;
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8) && etaCut) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

bool AliAnalysisTaskLambdaHadronMC::PassMCLambdaCuts(AliAODMCParticle *mc_particle, bool checkPhysicalPrimary, bool etaCut){

    if(!(TMath::Abs(mc_particle->GetPdgCode()) == 3122)) return false;
    if(checkPhysicalPrimary) {
        if(!mc_particle->IsPhysicalPrimary()) return false;
    }
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8) && etaCut) return false;

    int first_daughter_index = 0;
    int second_daughter_index = 0;

    first_daughter_index = mc_particle->GetDaughterFirst();
    second_daughter_index = mc_particle->GetDaughterLast();

    if(first_daughter_index < 0 || second_daughter_index < 0) return false;

    AliAODMCParticle* first_daughter = (AliAODMCParticle*)fMCEvent->GetTrack(first_daughter_index);
    AliAODMCParticle* second_daughter = (AliAODMCParticle*)fMCEvent->GetTrack(second_daughter_index);

    // make sure lambda decays into p-pi
    if(!(((TMath::Abs(first_daughter->GetPdgCode()) == 211 && TMath::Abs(second_daughter->GetPdgCode()) == 2212) 
        ||  (TMath::Abs(first_daughter->GetPdgCode()) == 2212 && TMath::Abs(second_daughter->GetPdgCode()) == 211)) 
        && (first_daughter->GetPdgCode())*(second_daughter->GetPdgCode()) < 0)) return false;
    
    return true;
}

bool AliAnalysisTaskLambdaHadronMC::PassMCPhiCuts(AliAODMCParticle *mc_particle, bool etaCut){

    if(!(TMath::Abs(mc_particle->GetPdgCode()) == 333)) return false;
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8) && etaCut) return false;

    int first_daughter_index = 0;
    int second_daughter_index = 0;

    first_daughter_index = mc_particle->GetDaughterFirst();
    second_daughter_index = mc_particle->GetDaughterLast();

    if(first_daughter_index < 0 || second_daughter_index < 0) return false;

    AliAODMCParticle* first_daughter = (AliAODMCParticle*)fMCEvent->GetTrack(first_daughter_index);
    AliAODMCParticle* second_daughter = (AliAODMCParticle*)fMCEvent->GetTrack(second_daughter_index);

    // make sure lambda decays into p-pi
    if(!((TMath::Abs(first_daughter->GetPdgCode()) == 321 && TMath::Abs(second_daughter->GetPdgCode()) == 321) 
        && (first_daughter->GetPdgCode())*(second_daughter->GetPdgCode()) < 0)) return false;
    
    return true;
}
        

void AliAnalysisTaskLambdaHadronMC::UserExec(Option_t*)
{

    fMCEvent  = dynamic_cast<AliMCEvent*> (MCEvent());
    if(!fMCEvent){
        Printf("No MC particle branch found");
        return;
    }

    // using 50 as mult percentile until I have actual multiplicity distribution
    double multPercentile = 50;

    AliVVertex * primVertex = (AliVVertex* ) fMCEvent->GetPrimaryVertex();
    double primZ = primVertex->GetZ();
    if (TMath::Abs(primZ) > 10) return;


    //MC trigger list used for event mixing
    TObjArray* fMixedMCTrackObjArray = new TObjArray;
    fMixedMCTrackObjArray->SetOwner(kTRUE);


    // MC SAME EVENT SECTION
    std::vector<AliAODMCParticle*> real_trigger_list;
    std::vector<AliAODMCParticle*> real_associated_list;
    std::vector<AliAODMCParticle*> real_lambda_list;
    std::vector<AliAODMCParticle*> real_phi_list;

    std::vector<AliAODMCParticle*> real_trigger_list_no_eta_cut;
    std::vector<AliAODMCParticle*> real_associated_list_no_eta_cut;
    std::vector<AliAODMCParticle*> real_lambda_list_no_eta_cut;
    std::vector<AliAODMCParticle*> real_phi_list_no_eta_cut;


    // bool to keep track if event has a trigger particle with pT > 4 GeV/c
    bool is_triggered_event_MC = false;

    int numTracksInV0A = 0;

    int numTracks = fMCEvent->GetNumberOfTracks();

    for(int mc_index = 0; mc_index < numTracks; mc_index++){

        AliAODMCParticle *mc_particle = (AliAODMCParticle*)fMCEvent->GetTrack(mc_index);

        if(PassMCTriggerCuts(mc_particle)) {
            real_trigger_list.push_back(mc_particle);
            AliCFParticle *trigger_particle = new AliCFParticle(mc_particle->Pt(), mc_particle->Eta(), mc_particle->Phi(), mc_particle->Charge(), 0);
            fMixedMCTrackObjArray->Add(trigger_particle);

            if(mc_particle->Pt() > 4) is_triggered_event_MC = true;
        }

        if(PassMCAssociatedCuts(mc_particle)) real_associated_list.push_back(mc_particle);
        if(PassMCLambdaCuts(mc_particle)) real_lambda_list.push_back(mc_particle);
        if(PassMCPhiCuts(mc_particle)) real_phi_list.push_back(mc_particle);

        if(PassMCTriggerCuts(mc_particle, false)) real_trigger_list_no_eta_cut.push_back(mc_particle);
        if(PassMCAssociatedCuts(mc_particle, false)) real_associated_list_no_eta_cut.push_back(mc_particle);
        if(PassMCLambdaCuts(mc_particle, false, false)) real_lambda_list_no_eta_cut.push_back(mc_particle);
        if(PassMCPhiCuts(mc_particle, false)) real_phi_list_no_eta_cut.push_back(mc_particle);

        if(mc_particle->Eta() > 2.8 && mc_particle->Eta() < 5.1 && mc_particle->IsPhysicalPrimary() && mc_particle->Charge()) numTracksInV0A++;
    }

    fMultDist->Fill(numTracksInV0A);

    FillSingleMCParticleDist(real_trigger_list, primZ, fTriggerDist_MC);
    FillSingleMCParticleDist(real_associated_list, primZ, fAssociatedDist_MC);
    FillSingleMCParticleDist(real_lambda_list, primZ, fLambdaDist_MC);
    FillSingleMCParticleDist(real_phi_list, primZ, fPhiDist_MC);

    FillSingleMCParticleDist(real_trigger_list_no_eta_cut, primZ, fTriggerDist_MC_no_eta_cut);
    FillSingleMCParticleDist(real_associated_list_no_eta_cut, primZ, fAssociatedDist_MC_no_eta_cut);
    FillSingleMCParticleDist(real_lambda_list_no_eta_cut, primZ, fLambdaDist_MC_no_eta_cut);
    FillSingleMCParticleDist(real_phi_list_no_eta_cut, primZ, fPhiDist_MC_no_eta_cut);

    if(is_triggered_event_MC)  {
        FillSingleMCParticleDist(real_trigger_list, primZ, fTriggeredTriggerDist_MC);
        FillSingleMCParticleDist(real_associated_list, primZ, fTriggeredAssociatedDist_MC);
        FillSingleMCParticleDist(real_lambda_list, primZ, fTriggeredLambdaDist_MC);
        FillSingleMCParticleDist(real_phi_list, primZ, fTriggeredPhiDist_MC);
    }



    MakeSameMCCorrelations(real_trigger_list, real_associated_list, fDphiHH_MC, primZ);
    MakeSameMCCorrelations(real_trigger_list, real_lambda_list, fDphiHLambda_MC, primZ);
    MakeSameMCCorrelations(real_trigger_list, real_phi_list, fDphiHPhi_MC, primZ);


    if(real_associated_list.size() > 0 || real_lambda_list.size() > 0) {
        AliEventPool *fMCCorPool = fMCCorPoolMgr->GetEventPool(50, primZ);
        if(!fMCCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fMCCorPool->IsReady()) {
                MakeMixedMCCorrelations(fMCCorPool, real_associated_list, fDphiHHMixed_MC, primZ);
                MakeMixedMCCorrelations(fMCCorPool, real_lambda_list, fDphiHLambdaMixed_MC, primZ);
                MakeMixedMCCorrelations(fMCCorPool, real_phi_list, fDphiHPhiMixed_MC, primZ);
            }
            if(fMixedMCTrackObjArray->GetEntries() > 0) {
                fMCCorPool->UpdatePool(fMixedMCTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronMC::Terminate(Option_t *option)
{
}
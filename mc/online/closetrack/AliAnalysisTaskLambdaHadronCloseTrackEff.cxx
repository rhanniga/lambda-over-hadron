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

#include "AliAnalysisTaskLambdaHadronCloseTrackEff.h"

ClassImp(AliAnalysisTaskLambdaHadronCloseTrackEff);

AliAnalysisTaskLambdaHadronCloseTrackEff::AliAnalysisTaskLambdaHadronCloseTrackEff() :
    AliAnalysisTaskSE(),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fTriggerBit(0.0),
    fCorPoolMgr(0x0),
    fMCCorPoolMgr(0x0)
{
}

AliAnalysisTaskLambdaHadronCloseTrackEff::AliAnalysisTaskLambdaHadronCloseTrackEff(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fTriggerBit(0.0),
    fCorPoolMgr(0x0),
    fMCCorPoolMgr(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronCloseTrackEff::~AliAnalysisTaskLambdaHadronCloseTrackEff()
{
    if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::UserCreateOutputObjects()
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


    //Correlation axes are: q (relative momenta), p (avg momenta), dphi, zvtx
    int hp_cor_bins[5] = {100, 100, 100, 32, 10};
    double hp_cor_mins[5] = {0, 0, 0, -1.0*TMath::Pi()/2.0, -10};
    double hp_cor_maxes[5] = {10, 10, 10, 3.0*TMath::Pi()/2.0, 10};

    fRecoHProton = new THnSparseF("fRecoHProton", "fRecoHProton", 4, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRecoHProton->Sumw2();
    fRecoHProton->GetAxis(0)->SetTitle("close (r=85cm)");
    fRecoHProton->GetAxis(1)->SetTitle("mid (r=167.5cm)");
    fRecoHProton->GetAxis(2)->SetTitle("far (r=250cm)");
    fRecoHProton->GetAxis(3)->SetTitle("dphi");
    fRecoHProton->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRecoHProton);

    fRecoHProtonMixed = new THnSparseF("fRecoHProtonMixed", "fRecoHProtonMixed", 4, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRecoHProtonMixed->Sumw2();
    fRecoHProtonMixed->GetAxis(0)->SetTitle("close (r=85cm)");
    fRecoHProtonMixed->GetAxis(1)->SetTitle("mid (r=167.5cm)");
    fRecoHProtonMixed->GetAxis(2)->SetTitle("far (r=250cm)");
    fRecoHProtonMixed->GetAxis(3)->SetTitle("dphi");
    fRecoHProtonMixed->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRecoHProtonMixed);

    fRealHProton = new THnSparseF("fRealHProton", "fRealHProton", 4, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRealHProton->Sumw2();
    fRealHProton->GetAxis(0)->SetTitle("close (r=85cm)");
    fRealHProton->GetAxis(1)->SetTitle("mid (r=167.5cm)");
    fRealHProton->GetAxis(2)->SetTitle("far (r=250cm)");
    fRealHProton->GetAxis(3)->SetTitle("dphi");
    fRealHProton->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRealHProton);

    fRealHProtonMixed = new THnSparseF("fRealHProtonMixed", "fRealHProtonMixed", 4, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRealHProtonMixed->Sumw2();
    fRealHProtonMixed->GetAxis(0)->SetTitle("close (r=85cm)");
    fRealHProtonMixed->GetAxis(1)->SetTitle("mid (r=167.5cm)");
    fRealHProtonMixed->GetAxis(2)->SetTitle("far (r=250cm)");
    fRealHProtonMixed->GetAxis(3)->SetTitle("dphi");
    fRealHProtonMixed->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRealHProtonMixed);

    PostData(1, fOutputList);

}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}


void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeSameHProtonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> proton_list, THnSparse* corr_dist, double bz, double zVtx)
{
    double corr_point[4];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        // for now we put the pt cuts here, maybe at axes later
        if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;

        for(int i = 0; i < (int)proton_list.size(); i++) {
            auto proton = proton_list[i];
            // for now we put the pt cuts here, maybe at axes later
            if(proton->Pt() < 2 || proton->Pt() > 4) continue;

            // make sure current trigger isn't current proton
            if(trigger->GetLabel() == proton->GetLabel()) continue;

            // getting q and p

            double close_point_trigger[3];
            trigger->GetXYZatR(85, bz, close_point_trigger);

            double mid_point_trigger[3];
            trigger->GetXYZatR(167.5, bz, mid_point_trigger);

            double far_point_trigger[3];
            trigger->GetXYZatR(250, bz, far_point_trigger);

            double close_point_proton[3];
            proton->GetXYZatR(85, bz, close_point_proton);

            double mid_point_proton[3];
            proton->GetXYZatR(167.5, bz, mid_point_proton);

            double far_point_proton[3];
            proton->GetXYZatR(250, bz, far_point_proton);

            double close_distance = sqrt(pow(close_point_trigger[0] - close_point_proton[0], 2) + pow(close_point_trigger[1] - close_point_proton[1], 2) + pow(close_point_trigger[2] - close_point_proton[2], 2));
            double mid_distance = sqrt(pow(mid_point_trigger[0] - mid_point_proton[0], 2) + pow(mid_point_trigger[1] - mid_point_proton[1], 2) + pow(mid_point_trigger[2] - mid_point_proton[2], 2));
            double far_distance = sqrt(pow(far_point_trigger[0] - far_point_proton[0], 2) + pow(far_point_trigger[1] - far_point_proton[1], 2) + pow(far_point_trigger[2] - far_point_proton[2], 2));

            corr_point[0] = close_distance;
            corr_point[1] = mid_distance;
            corr_point[2] = far_distance;


            corr_point[3] = trigger->Phi() - proton->Phi();

            if(corr_point[3] < -TMath::Pi()/2.0) {
                corr_point[3] += 2.0*TMath::Pi();
            }
            else if(corr_point[3] > 3.0*TMath::Pi()/2.0) {
                corr_point[3] -= 2.0*TMath::Pi();
            }

            corr_point[4] = zVtx;

            corr_dist->Fill(corr_point);

        }
    }
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeSameMCHProtonCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> proton_list, THnSparse* corr_dist, double zVtx)
{
    double corr_point[4];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        // for now we put the pt cuts here, maybe at axes later
        if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;

        for(int i = 0; i < (int)proton_list.size(); i++) {
            auto proton = proton_list[i];
            // for now we put the pt cuts here, maybe at axes later
            if(proton->Pt() < 2 || proton->Pt() > 4) continue;

            // make sure current trigger isn't current proton
            if(trigger->GetLabel() == proton->GetLabel()) continue;

            // getting q and p
            corr_point[0] = TMath::Sq(trigger->Px() - proton->Px()) + 
                            TMath::Sq(trigger->Py() - proton->Py()) +
                            TMath::Sq(trigger->Pz() - proton->Pz());

            corr_point[0] = TMath::Sqrt(corr_point[0]);

            corr_point[1] = TMath::Sq((trigger->Px() + proton->Px()) / 2.0) +
                            TMath::Sq((trigger->Py() + proton->Py()) / 2.0) +
                            TMath::Sq((trigger->Pz() + proton->Pz()) / 2.0);

            corr_point[1] = TMath::Sqrt(corr_point[1]);

            corr_point[2] = trigger->Phi() - proton->Phi();

            if(corr_point[2] < -TMath::Pi()/2.0) {
                corr_point[2] += 2.0*TMath::Pi();
            }
            else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                corr_point[2] -= 2.0*TMath::Pi();
            }

            corr_point[3] = zVtx;

            corr_dist->Fill(corr_point);

        }
    }
}



void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeMixedHProtonCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> proton_list , THnSparse* corr_dist, double zVtx)
{
    
    double corr_point[4];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            // for now we put the pt cuts here, maybe at axes later
            if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;
            for(int j = 0; j < (int)proton_list.size(); j++) {
                auto proton = proton_list[j];
                // for now we put the pt cuts here, maybe at axes later
                if(proton->Pt() < 2 || proton->Pt() > 4) continue;

                float tpx = trigger->Pt() * TMath::Cos(trigger->Phi());
                float tpy = trigger->Pt() * TMath::Sin(trigger->Phi());
                float tpz = trigger->Pt() * TMath::SinH(trigger->Eta());

                corr_point[0] = TMath::Sq(tpx - proton->Px()) + 
                                TMath::Sq(tpy - proton->Py()) +
                                TMath::Sq(tpz - proton->Pz());

                corr_point[0] = TMath::Sqrt(corr_point[0]);

                corr_point[1] = TMath::Sq((tpx + proton->Px()) / 2.0) +
                                TMath::Sq((tpy + proton->Py()) / 2.0) +
                                TMath::Sq((tpz + proton->Pz()) / 2.0);

                corr_point[1] = TMath::Sqrt(corr_point[1]);

                corr_point[2] = trigger->Phi() - proton->Phi();

                if(corr_point[2] < -TMath::Pi()/2.0) {
                    corr_point[2] += 2.0*TMath::Pi();
                }
                else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                    corr_point[2] -= 2.0*TMath::Pi();
                }

                corr_point[3] = zVtx;
                corr_dist->Fill(corr_point);

            }
        }
    }
}
void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeMixedMCHProtonCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> proton_list , THnSparse* corr_dist, double zVtx)
{
    
    double corr_point[4];
    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            // for now we put the pt cuts here, maybe at axes later
            if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;
            for(int j = 0; j < (int)proton_list.size(); j++) {
                auto proton = proton_list[j];
                // for now we put the pt cuts here, maybe at axes later
                if(proton->Pt() < 2 || proton->Pt() > 4) continue;

                float tpx = trigger->Pt() * TMath::Cos(trigger->Phi());
                float tpy = trigger->Pt() * TMath::Sin(trigger->Phi());
                float tpz = trigger->Pt() * TMath::SinH(trigger->Eta());

                corr_point[0] = TMath::Sq(tpx - proton->Px()) + 
                                TMath::Sq(tpy - proton->Py()) +
                                TMath::Sq(tpz - proton->Pz());

                corr_point[0] = TMath::Sqrt(corr_point[0]);

                corr_point[1] = TMath::Sq((tpx + proton->Px()) / 2.0) +
                                TMath::Sq((tpy + proton->Py()) / 2.0) +
                                TMath::Sq((tpz + proton->Pz()) / 2.0);

                corr_point[1] = TMath::Sqrt(corr_point[1]);

                corr_point[2] = trigger->Phi() - proton->Phi();

                if(corr_point[2] < -TMath::Pi()/2.0) {
                    corr_point[2] += 2.0*TMath::Pi();
                }
                else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                    corr_point[2] -= 2.0*TMath::Pi();
                }

                corr_point[3] = zVtx;
                corr_dist->Fill(corr_point);

            }
        }
    }
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassDaughterCuts(AliAODTrack *track){

    if(track->GetID() < 0) return false;

    bool pass = true;

    pass = pass && (TMath::Abs(track->Eta()) < 0.8);
    pass = pass && (track->Pt() > 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    int label = track->GetLabel();

    if(label < 0) return false;

    AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(label);

    // we are only concerned with secondary protons
    if(mcpart->IsPhysicalPrimary()) return false;

    return pass;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassTriggerCuts(AliAODTrack *track, bool checkMC){

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

bool AliAnalysisTaskLambdaHadronCloseTrackEff::IsMCChargedHadron(int pdg_code) {
    // checks to see if pdg code matches pion, proton, kaon, electron or muon
    if((TMath::Abs(pdg_code) == 321)
        || (TMath::Abs(pdg_code) == 211)
        || (TMath::Abs(pdg_code) == 2212)
        || (TMath::Abs(pdg_code) == 11)
        || (TMath::Abs(pdg_code) == 13)) return true;

    return false;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassMCTriggerCuts(AliAODMCParticle *mc_particle){

    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false; // for now trigger is physical primary, could change
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::UserExec(Option_t*)
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
    std::vector<AliAODTrack*> proton_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    //MC trigger list used for event mixing
    TObjArray* fMixedMCTrackObjArray = new TObjArray;
    fMixedMCTrackObjArray->SetOwner(kTRUE);

    // RECO SAME EVENT SECTION
    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);
            }

        if(PassDaughterCuts(track)) {
            int mc_label = track->GetLabel();
            if(mc_label >= 0) {
                auto mc_part = (AliAODMCParticle*)fMCArray->At(mc_label);
                if(TMath::Abs(mc_part->GetPdgCode()) == 2212) {
                    // guaranteed secondary protons that pass loose daughter cuts
                    proton_list.push_back(track);
                }
            }
        }
    }

    MakeSameHProtonCorrelations(trigger_list, proton_list, fRecoHProton, fAOD->GetMagneticField(), primZ);

    // MC SAME EVENT SECTION

    std::vector<AliAODMCParticle*> real_trigger_list;
    std::vector<AliAODMCParticle*> real_proton_list;

    for(int mc_index = 0; mc_index < fMCArray->GetEntries(); mc_index++){

        AliAODMCParticle *mc_particle = (AliAODMCParticle*)fMCArray->At(mc_index);

        if(PassMCTriggerCuts(mc_particle)) {
            real_trigger_list.push_back(mc_particle);
            AliCFParticle *trigger_particle = new AliCFParticle(mc_particle->Pt(), mc_particle->Eta(), mc_particle->Phi(), mc_particle->Charge(), 0);
            fMixedMCTrackObjArray->Add(trigger_particle);
        }
        if(TMath::Abs(mc_particle->Eta()) < 0.8 && mc_particle->Pt() > 0.15) {
            // again selecting real protons, but not physical primaries
            if(TMath::Abs(mc_particle->GetPdgCode()) == 2212 && !mc_particle->IsPhysicalPrimary()) real_proton_list.push_back(mc_particle);
        }
    }

    MakeSameMCHProtonCorrelations(real_trigger_list, real_proton_list, fRealHProton, primZ);

    // MIXED EVENT SECTION (added to very end to correctly do a mixture of reco/real correlations)

    if(proton_list.size() > 0 ) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fCorPool->IsReady()) {
                MakeMixedHProtonCorrelations(fCorPool, proton_list, fRecoHProtonMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }


    if(real_proton_list.size() > 0) {
        AliEventPool *fMCCorPool = fMCCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fMCCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fMCCorPool->IsReady()) {
                MakeMixedMCHProtonCorrelations(fMCCorPool, real_proton_list, fRealHProtonMixed, primZ);
            }
            if(fMixedMCTrackObjArray->GetEntries() > 0) {
                fMCCorPool->UpdatePool(fMixedMCTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::Terminate(Option_t *option)
{
}

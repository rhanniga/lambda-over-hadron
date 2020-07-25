#include "AliAnalysisTaskLambdaHadronRatio.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskLambdaHadronRatio;
ClassImp(AliAnalysisTaskLambdaHadronRatio);

static const int centLow = 0;
static const int centHigh = 100;

AliAnalysisTaskLambdaHadronRatio::AliAnalysisTaskLambdaHadronRatio() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHLambda{0},
    fDphiHLambdaRotated{0},
    fDphiHLambdaRotatedProton{0},
    fDphiHLambdaRotatedPi{0},
    fDphiHLambdaFlipped{0},
    fDphiHH{0},
    fDphiHLambdaLS{0},
    fDphiHLambdaMixed{0},
    fDphiHHMixed{0},
    fDphiHLambdaLSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{

}

AliAnalysisTaskLambdaHadronRatio::AliAnalysisTaskLambdaHadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHLambda{0},
    fDphiHLambdaRotated{0},
    fDphiHLambdaRotatedProton{0},
    fDphiHLambdaRotatedPi{0},
    fDphiHLambdaFlipped{0},
    fDphiHH{0},
    fDphiHLambdaLS{0},
    fDphiHLambdaMixed{0},
    fDphiHHMixed{0},
    fDphiHLambdaLSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

}

AliAnalysisTaskLambdaHadronRatio::~AliAnalysisTaskLambdaHadronRatio()
{

    if(fOutputList) delete fOutputList;

}

void AliAnalysisTaskLambdaHadronRatio::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:

    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {centLow, centHigh};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    fTriggersAndLambdasPerEvent_All = new TH2D("fTriggersAndLambdasPerEvent_All", "Triggers and Lambdas per event (all p_{T})", 10, 0, 10, 10, 0, 10);
    fTriggersAndLambdasPerEvent_2_4 = new TH2D("fTriggersAndLambdasPerEvent_2_4", "Triggers and Lambdas per event (2-4 p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndLambdasPerEvent_All);
    fOutputList->Add(fTriggersAndLambdasPerEvent_2_4);


    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {100, 16, 20, 10};
    double dist_mins[4] = {0, 0, -1, -10};
    double dist_maxes[4] = {15, 6.28, 1, 10};

    fLooseDist = new THnSparseF("fLooseDist", "All Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fLambdaDist = new THnSparseF("fLambdaDist", "Lambda Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fLooseDist);
    fOutputList->Add(fTriggerDist);
    fOutputList->Add(fAssociatedHDist);
    fOutputList->Add(fLambdaDist);
    

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Associated Inv Mass (Lambda only), Zvtx
    int hk_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    double hk_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hk_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10};

    int hh_cor_bins[5] = {20, 20, 16, 20, 10};
    double hh_cor_mins[5] = {2, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 12, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHLambda = new THnSparseF("fDphiHLambda", "Hadron-Lambda Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaV0 = new THnSparseF("fDphiHLambdaV0", "Hadron-Lambda (using V0) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaRotated = new THnSparseF("fDphiHLambdaRotated", "Hadron-Lambda (rotated) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaRotatedPi = new THnSparseF("fDphiHLambdaRotatedPi", "Hadron-Lambda (rotated pi) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaRotatedProton = new THnSparseF("fDphiHLambdaRotatedProton", "Hadron-Lambda (proton rotated) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaFlipped = new THnSparseF("fDphiHLambdaFlipped", "Hadron-Lambda (flipped) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTrigger = new THnSparseF("fDphiTriggerTrigger", "Trigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHLambdaLS = new THnSparseF("fDphiHLambdaLS", "Hadron-Lambda LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHLambdaMixed = new THnSparseF("fDphiHLambdaMixed", "Mixed Hadron-Lambda Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTriggerMixed = new THnSparseF("fDphiTriggerTriggerMixed", "MixedTrigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHLambdaLSMixed = new THnSparseF("fDphiHLambdaLSMixed", "Mixed Hadron-Lambda LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fOutputList->Add(fDphiHLambda);
    fOutputList->Add(fDphiHLambdaV0);
    fOutputList->Add(fDphiHLambdaRotated);
    fOutputList->Add(fDphiHLambdaRotatedPi);
    fOutputList->Add(fDphiHLambdaRotatedProton);
    fOutputList->Add(fDphiHLambdaFlipped);
    fOutputList->Add(fDphiHH);
    fOutputList->Add(fDphiTriggerTrigger);
    fOutputList->Add(fDphiHLambdaLS);
    fOutputList->Add(fDphiHLambdaMixed);
    fOutputList->Add(fDphiHHMixed);
    fOutputList->Add(fDphiTriggerTriggerMixed);
    fOutputList->Add(fDphiHLambdaLSMixed);

    //axes are pion DCA proton DCA lambda mass
    fLambdaDaughterDCA = new TH3D("fLambdaDaughterDCA", "#Lambda^{0} daughter DCA dist", 100, -5, 5, 100, -5, 5, 100, 1.08, 1.16);
    fOutputList->Add(fLambdaDaughterDCA);

    PostData(1, fOutputList);

}

void AliAnalysisTaskLambdaHadronRatio::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist)
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

void AliAnalysisTaskLambdaHadronRatio::FillSingleParticleDist(std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = zVtx;
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

void AliAnalysisTaskLambdaHadronRatio::MakeSameHLambdaCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list, THnSparse* fDphi, double zVtx)
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
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskLambdaHadronRatio::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskLambdaHadronRatio::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx)
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
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskLambdaHadronRatio::MakeMixedHLambdaCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list , THnSparse* fDphi, double zVtx)
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
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskLambdaHadronRatio::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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

void AliAnalysisTaskLambdaHadronRatio::UserExec(Option_t*)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        std::cout << "THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!" << std::endl;
        std::abort();
    }

    fpidResponse = fInputHandler->GetPIDResponse();


    //Event cuts
    TString cent_estimator = "V0A";
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < centLow || multPercentile > centHigh) return;

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
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> all_hadron_list;
    std::vector<AliAODTrack*> k_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //List for comparison with cuts/filter bits
        all_hadron_list.push_back(track);

        //Filter for trigger particles
        bool trigFilter = track->TestBit(AliAODTrack::kIsHybridGCG);
        if(trigFilter) {

            if(track->Pt() > 4 && track->Pt() < 8 && TMath::Abs(track->Eta()) < 1) {
                trigger_list.push_back(track);
                AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
                fMixedTrackObjArray->Add(triggerPart);
           }
        }

        //Filter for associated particles

        bool assocFilter = track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA);
        if(assocFilter) {

            if(track->Pt() > 0.15 && TMath::Abs(track->Eta()) < 1) {

                associated_h_list.push_back(track);

                double TPCNSigmaPion = 1000;
                double TOFNSigmaPion = 1000;

                TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
                TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

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

                // double pid_array[7] = {track->P(), track->GetTPCsignal(), track->GetTOFsignal(), TPCNSigmaPion, TOFNSigmaPion, TPCNSigmaProton, TOFNSigmaProton};
                // fPid->Fill(pid_array);
            
            }
        }
    }

    //Making list of possible lambdas (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_v0;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_signal_region;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_signal_region_2_4;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedPion;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedProton;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_RotatedPi;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_Flipped;
    std::vector<AliAnalysisTaskLambdaHadronRatio::AliMotherContainer> lambda_list_LS;

    // for(int i = 0; i < (int)unlikelyPion_list.size(); i++) {
    //     for(int j = 0; j < (int) unlikelyProton_list.size(); j++) {

    //         double TPCNSigmaPion = 1000;
    //         double TOFNSigmaPion = 1000;
    //         double TPCNSigmaProton = 1000;
    //         double TOFNSigmaProton = 1000;

    //         TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(unlikelyPion_list[i], AliPID::kPion);
    //         TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(unlikelyPion_list[i], AliPID::kPion);
    //         TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(unlikelyProton_list[j], AliPID::kProton);
    //         TOFNSigmaProton = fpidResponse->NumberOfSigmasTOF(unlikelyProton_list[j], AliPID::kProton);

    //         AliMotherContainer mother = DaughtersToMother(unlikelyPion_list[i], unlikelyProton_list[j], 0.1396, 0.9383);
    //         double mass = mother.particle.M();
    //         double pt = mother.particle.Pt();

    //         double signal_array[6] = {TPCNSigmaPion, TOFNSigmaPion, TPCNSigmaProton, TOFNSigmaProton, mass, pt};
    //         fSignalAnalysis->Fill(signal_array);
    //     }
    // }


    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) proton_list.size(); j++) {
            AliMotherContainer lambda = DaughtersToMother(piMinus_list[i], proton_list[j], 0.1396, 0.9383);
            
            double pion_dz[2];
            double pion_covar[3];

            double proton_dz[2];
            double proton_covar[3];

            bool is_pionDCA = piMinus_list[i]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., pion_dz, pion_covar);
            bool is_protonDCA = proton_list[j]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., proton_dz, proton_covar);

            if(is_pionDCA && is_protonDCA) {
                fLambdaDaughterDCA->Fill(pion_dz[0], proton_dz[0], lambda.particle.M());
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

        // int posTrackID = v0->GetPosID();
        // int negTrackID = v0->GetNegID();
        // AliAODTrack* posTrack = (AliAODTrack*)fAOD->GetTrack(posTrackID);
        // AliAODTrack* negTrack = (AliAODTrack*)fAOD->GetTrack(negTrackID);

        AliAODTrack* posTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
        AliAODTrack* negTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);

        // Occasionally returns null, not quite sure why...
        if(!posTrack || !negTrack) continue;

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
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);
    FillSingleParticleDist(all_hadron_list, primZ, fLooseDist);
    FillSingleParticleDist(lambda_list_signal_region, primZ, fLambdaDist);

    // Filling all of our correlation histograms
    MakeSameHLambdaCorrelations(trigger_list, lambda_list, fDphiHLambda, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_v0, fDphiHLambdaV0, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedPi, fDphiHLambdaRotatedPi, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_Flipped, fDphiHLambdaFlipped, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedPion, fDphiHLambdaRotated, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_RotatedProton, fDphiHLambdaRotatedProton, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ);
    MakeSameTriggerTriggerCorrelations(trigger_list, fDphiTriggerTrigger, primZ);
    MakeSameHLambdaCorrelations(trigger_list, lambda_list_LS, fDphiHLambdaLS, primZ);


    fTriggersAndLambdasPerEvent_All->Fill(trigger_list.size(), lambda_list_signal_region.size());
    fTriggersAndLambdasPerEvent_2_4->Fill(trigger_list.size(), lambda_list_signal_region_2_4.size());


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

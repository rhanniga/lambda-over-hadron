#include "AliAnalysisTaskStrangeCor.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskStrangeCor;
ClassImp(AliAnalysisTaskStrangeCor);

static const int centLow = 20;
static const int centHigh = 50;

AliAnalysisTaskStrangeCor::AliAnalysisTaskStrangeCor() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHStrangePart{0},
    fDphiHH{0},
    fDphiHStrangePartLS{0},
    fDphiHStrangePartMixed{0},
    fDphiHHMixed{0},
    fDphiHStrangePartLSMixed{0},
    fCorPoolMgr{0}
{

}

AliAnalysisTaskStrangeCor::AliAnalysisTaskStrangeCor(const char *name) :
    AliAnalysisTaskSE(name),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHStrangePart{0},
    fDphiHH{0},
    fDphiHStrangePartLS{0},
    fDphiHStrangePartMixed{0},
    fDphiHHMixed{0},
    fDphiHStrangePartLSMixed{0},
    fCorPoolMgr{0}
{

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

}

AliAnalysisTaskStrangeCor::~AliAnalysisTaskStrangeCor()
{

    if(fOutputList) delete fOutputList;

}

void AliAnalysisTaskStrangeCor::UserCreateOutputObjects()
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

    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {15, 16, 20, 10};
    double dist_mins[4] = {0, 0, -1, -10};
    double dist_maxes[4] = {15, 6.28, 1, 10};

    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fTriggerDist);
    fOutputList->Add(fAssociatedHDist);



    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Associated Inv Mass (K* only)
    int hk_cor_bins[6] = {8, 4, 16, 20, 100, 10};
    double hk_cor_mins[6] = {4.0, 2, -1.0*TMath::Pi()/2.0, -2.0, 1.06, -10};
    double hk_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 1.16, 10};

    int hh_cor_bins[5] = {8, 4, 16, 20, 10};
    double hh_cor_mins[5] = {4, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 10};



    fDphiHStrangePart = new THnSparseF("fDphiHStrangePart", "Hadron-K* Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHStrangePartLS = new THnSparseF("fDphiHStrangePartLS", "Hadron-K* LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHStrangePartMixed = new THnSparseF("fDphiHStrangePartMixed", "Mixed Hadron-K* Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHStrangePartLSMixed = new THnSparseF("fDphiHStrangePartLSMixed", "Mixed Hadron-K* LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fOutputList->Add(fDphiHStrangePart);
    fOutputList->Add(fDphiHH);
    fOutputList->Add(fDphiHStrangePartLS);
    fOutputList->Add(fDphiHStrangePartMixed);
    fOutputList->Add(fDphiHHMixed);
    fOutputList->Add(fDphiHStrangePartLSMixed);

    PostData(1, fOutputList);

}

void AliAnalysisTaskStrangeCor::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist)
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

AliAnalysisTaskStrangeCor::AliMotherContainer AliAnalysisTaskStrangeCor::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskStrangeCor::AliMotherContainer mom;
    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

void AliAnalysisTaskStrangeCor::MakeSameHStrangePartCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> strange_part_list, THnSparse* fDphi, double zVtx)
{

    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)strange_part_list.size(); i++) {
            auto strange_part = strange_part_list[i];

            //Make sure trigger isn't one of the daughters of kstar
            if((trigger->GetID() == strange_part.daughter1ID) || (trigger->GetID() == strange_part.daughter2ID)) continue;

            dphi_point[1] = strange_part.particle.Pt();
            dphi_point[2] = trigger->Phi() - strange_part.particle.Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - strange_part.particle.Eta();
            dphi_point[4] = strange_part.particle.M();
            dphi_point[5] = zVtx;
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskStrangeCor::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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

void AliAnalysisTaskStrangeCor::MakeMixedHStrangePartCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> strange_part_list , THnSparse* fDphi, double zVtx)
{

    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)strange_part_list.size(); j++) {
                auto strange_part = strange_part_list[j];

                dphi_point[1] = strange_part.particle.Pt();
                dphi_point[2] = trigger->Phi() - strange_part.particle.Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - strange_part.particle.Eta();
                dphi_point[4] = strange_part.particle.M();
                dphi_point[5] = zVtx;
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskStrangeCor::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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

void AliAnalysisTaskStrangeCor::UserExec(Option_t*)
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

    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> kPlus_list;
    std::vector<AliAODTrack*> kMinus_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> k_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //Filter for trigger particles
        if(track->TestBit(AliAODTrack::kIsHybridGCG)) {

            if(track->Pt() > 4 && TMath::Abs(track->Eta()) < 1) {
                trigger_list.push_back(track);
                AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
                fMixedTrackObjArray->Add(triggerPart);
           }
        }

        //Filter for associated particles
        if(true) {

            if(track->Pt() > 0.5 && TMath::Abs(track->Eta()) < 1) {

                associated_h_list.push_back(track);

                double TPCNSigmaPion = 1000;
                double TOFNSigmaPion = 1000;

                TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
                TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

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

                if(TMath::Abs(TPCNSigmaProton) <= 2 && (TMath::Abs(TOFNSigmaProton) <= 2 || TOFNSigmaProton == 1000)) {
                    if(track->Charge() == 1){
                        kPlus_list.push_back(track);
                        k_list.push_back(track);
                    }
                    else {
                        kMinus_list.push_back(track);
                        k_list.push_back(track);
                    }
                }
            }
        }
    }

    //Making list of possible kStars (have to do +/- for K or pi):

    std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> kStar_list;
    std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> kStar_list_LS;
    std::vector<AliAnalysisTaskStrangeCor::AliMotherContainer> kStar_list_LS_same;

    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) kPlus_list.size(); j++) {
            AliMotherContainer kStar = DaughtersToMother(piMinus_list[i], kPlus_list[j], 0.1396, 0.9383);
            kStar_list.push_back(kStar);
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) kMinus_list.size(); j++) {
            AliMotherContainer kStar = DaughtersToMother(piPlus_list[i], kMinus_list[j], 0.1396, 0.9383);
            kStar_list.push_back(kStar);
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) kPlus_list.size(); j++) {
            if(piPlus_list[i]->GetID() == kPlus_list[j]->GetID()) continue;
            AliMotherContainer kStar = DaughtersToMother(piPlus_list[i], kPlus_list[j], 0.1396, 0.9383);
            kStar_list_LS.push_back(kStar);
        }
    }

    for(int i = 0; i < (int)piMinus_list.size(); i++) {
        for(int j = 0; j < (int) kMinus_list.size(); j++) {
            if(piMinus_list[i]->GetID() == kMinus_list[j]->GetID()) continue;
            AliMotherContainer kStar = DaughtersToMother(piMinus_list[i], kMinus_list[j], 0.1396, 0.9383);
            kStar_list_LS.push_back(kStar);
        }
    }

    //Filling all of our single particle distribution histograms:

    FillSingleParticleDist(trigger_list, primZ, fTriggerDist);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);

    //Filling all of our correlation histograms

    MakeSameHStrangePartCorrelations(trigger_list, kStar_list, fDphiHStrangePart, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ);
    MakeSameHStrangePartCorrelations(trigger_list, kStar_list_LS, fDphiHStrangePartLS, primZ);


    if(kStar_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHStrangePartCorrelations(fCorPool, kStar_list, fDphiHStrangePartMixed, primZ);
                MakeMixedHStrangePartCorrelations(fCorPool, kStar_list_LS, fDphiHStrangePartLSMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);

}

void AliAnalysisTaskStrangeCor::Terminate(Option_t *option)
{
}







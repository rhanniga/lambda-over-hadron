
#include "AliAnalysisTaskNSigmaElectron.h"
#include "AliKFParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "TGeoGlobalMagField.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"


#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskNSigmaElectron;
ClassImp(AliAnalysisTaskNSigmaElectron);


AliAnalysisTaskNSigmaElectron::AliAnalysisTaskNSigmaElectron() :
    AliAnalysisTaskSE()
{
}

AliAnalysisTaskNSigmaElectron::AliAnalysisTaskNSigmaElectron(const char *name) :
    AliAnalysisTaskSE(name)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskNSigmaElectron::~AliAnalysisTaskNSigmaElectron()
{
    if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskNSigmaElectron::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);


    // the only distribution we need

    // bins are momentum, tpc signal, nsigmaelec, nsigmapi, nsigmak, nsigmap
    int nBins[6] = {100, 120, 200, 200, 200, 200};
    double mins[6] = {0, 0, -10, -10, -10, -10};
    double maxes[6] = {10, 120, 10, 10, 10, 10};

    fBigDist = new THnSparseF("fBigDist", "fBigDist", 6, nBins, mins, maxes);
    fBigDist->Sumw2();
    fOutputList->Add(fBigDist);

    fNSigmaElectronBethe = new TH2D("fNSigmaElectronBethe", "fNSigmaElectronBethe", 200, 0, 10, 200, -10, 10);
    fOutputList->Add(fNSigmaElectronBethe);

    fNSigmaPionBethe = new TH2D("fNSigmaPionBethe", "fNSigmaPionBethe", 200, 0, 10, 200, -10, 10);
    fOutputList->Add(fNSigmaPionBethe);

    fNSigmaKaonBethe = new TH2D("fNSigmaKaonBethe", "fNSigmaKaonBethe", 200, 0, 10, 200, -10, 10);
    fOutputList->Add(fNSigmaKaonBethe);

    fNSigmaProtonBethe = new TH2D("fNSigmaProtonBethe", "fNSigmaProtonBethe", 200, 0, 10, 200, -10, 10);
    fOutputList->Add(fNSigmaProtonBethe);
    
    PostData(1, fOutputList);

}

bool AliAnalysisTaskNSigmaElectron::PassAssociatedCutsElectrons(AliAODTrack *track){

  Bool_t pass = kTRUE;
    
  pass = pass && (track->Pt() >= 0.15);
  pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
   
  pass = pass && track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA);//
    
    
  return pass;

}

void AliAnalysisTaskNSigmaElectron::UserExec(Option_t*)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
    }


    fpidResponse = fInputHandler->GetPIDResponse();

    //Event cuts
    TString cent_estimator = "V0A";
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < 0.0 || multPercentile > 100.0) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();


    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;


        //------------------ THIS IS THE SECTION THAT RYAN ADDED ----------------------

        double TOFNSigmaElec = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);

        // we're going to make the nsigma plot look even BETTER by requiring our track pass the
        // associated cuts
        if(PassAssociatedCutsElectrons(track) && TMath::Abs(TOFNSigmaElec) <= 3) {

          double momentum = track->P();
          double TPCSignal = track->GetTPCsignal();
          double TPCNSigmaElec = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
          double TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
          double TPCNSigmaKaon = fpidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
          double TPCNSigmaProton = fpidResponse->NumberOfSigmasTPC(track, AliPID::kProton);

          double fillArray[6] = {momentum, TPCSignal, TPCNSigmaElec, TPCNSigmaPion, TPCNSigmaKaon, TPCNSigmaProton};

          fBigDist->Fill(fillArray);

          AliTPCPIDResponse TPCResponse = fpidResponse->GetTPCResponse();

          double expectedSigmaElectron = TPCResponse.GetExpectedSigma(track, AliPID::kElectron);

          double expectedSignalElectron = TPCResponse.GetExpectedSignal(track, AliPID::kElectron);
          double expectedSignalPion = TPCResponse.GetExpectedSignal(track, AliPID::kPion);
          double expectedSignalKaon = TPCResponse.GetExpectedSignal(track, AliPID::kKaon);
          double expectedSignalProton = TPCResponse.GetExpectedSignal(track, AliPID::kProton);

          double nSigmaBetheElectron = -(expectedSignalElectron - expectedSignalElectron)/expectedSigmaElectron;
          double nSigmaBethePion = -(expectedSignalElectron - expectedSignalPion)/expectedSigmaElectron;
          double nSigmaBetheKaon = -(expectedSignalElectron - expectedSignalKaon)/expectedSigmaElectron;
          double nSigmaBetheProton = -(expectedSignalElectron - expectedSignalProton)/expectedSigmaElectron;

          fNSigmaElectronBethe->Fill(momentum, nSigmaBetheElectron);
          fNSigmaPionBethe->Fill(momentum, nSigmaBethePion);
          fNSigmaKaonBethe->Fill(momentum, nSigmaBetheKaon);
          fNSigmaProtonBethe->Fill(momentum, nSigmaBetheProton);

        }
    }

  PostData(1, fOutputList);//writes to output list
}

void AliAnalysisTaskNSigmaElectron::Terminate(Option_t *option)
{
}
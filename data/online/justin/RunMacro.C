#include "AliAnalysisTaskHadronPhiCorr_current.h"

void LoadLibraries();

void RunMacro(){

   // Firstly, set some variables
   const char* launch = "local"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "full"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kTRUE; //TRUE = merging done on grid, FALSE = merge happens locally   
   Int_t cyclenumber = 1;
   Bool_t debug = kTRUE;
   char* work_dir = "PhiCorrelations_LHC16q";
   //char* output_dir = "output_2020_07_30_hphi_alltrig_0_100";
   char* output_dir = "output_2021_09_14_hh_test_0_20";
   Int_t ttl = 50000;
   Int_t noffiles = 50;
   Int_t runcycle[]={0,16,31};
   //Int_t runcycle[] ={0,31};
   //Int_t runcycle[]={0,1,12,27};
   //Int_t runcycle[]={0,27};
   Bool_t UseParfiles = kFALSE;

   // load libraries
   LoadLibraries();

   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
   gROOT->ProcessLine(".include $PWD/.");

   gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_PHYSICS/PWGCF/Correlations/Base -I$ALICE_PHYSICS/PWGCF/Correlations -g ");

   //printf("\n!!!!!!!!!!!!!!!!!!!!!!\n AliAnalysis Manager \n\n");
   AliAnalysisManager *mgr = new AliAnalysisManager("PhiAnalysis");
   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);

    //switch on aliphysicsselection
    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE, kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C")))); //"kTRUE, kTRUE" for MC, "kFALSE, kTRUE" for Data 
    
    gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")));
    //Only set true for MC
    Bool_t isMC = kTRUE;
    
    gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")));


    TFile* testfile = TFile::Open("./justin_eff.root");
    std::cout << testfile << " alskdfglaisjdgflaskdfgalskcbvas'gkhawetlkjgqsdfkjg" << std::endl;
    AliAnalysisTaskHadronPhiCorr_current *task1 = reinterpret_cast<AliAnalysisTaskHadronPhiCorr_current*>(gInterpreter->ProcessLine(Form(".x %s(kTRUE, 0.0, 20.0)", gSystem->ExpandPathName("AddTaskHadronPhiCorr_current.C"))));

    task1->SetIsMCTrue(kFALSE);
    task1->SetIsMCKaon(kFALSE);
    task1->SetIsMCKTrack(kFALSE);
    task1->SetUseAccpt(kTRUE);
    task1->SetSingleTrigger(kFALSE);
    task1->SetSelectTrigger(kFALSE);
    task1->SetHighestTriggerOnly(kFALSE);
    task1->SetKaonEtaCut(0.8);
    task1->SetKaonTPCCut(3.0);
    task1->SetKaonTOFCut(3.0);
    task1->SetTOFVeto(kFALSE);
    task1->SetKaonTrkBit(1024);
    task1->SetAssocTrkBit(1024);
    task1->SetTrigTrkBit(AliAODTrack::kIsHybridGCG);
    task1->SetZVertexMin(-10.0);
    task1->SetZVertexMax(10.0);
    task1->SetZVertexNbins(10);
    task1->SetCentEstimator("V0A");
    task1->SetEtaPhiRegion(0);
    
    TF1* phiEff = (TF1*)(testfile->Get("phiFit")->Clone("phiEff"));
    TF1* hEff = (TF1*)(testfile->Get("hFit")->Clone("hEff"));
    TF1* trigEff = (TF1*)(testfile->Get("trigFit")->Clone("trigEff"));

    task1->LoadEfficiencies(phiEff, hEff, trigEff);

    if (!mgr->InitAnalysis())
     return;

   printf("mgr initialized analysis\n");
   mgr->PrintStatus();
   mgr->SetUseProgressBar(1, 25);
   fprintf(stdout, "\n!!!!!!!!!!!!!\nAbout to launch analysis... \n");
    TChain *chain = new TChain("aodTree");
    chain->Add("~/Wonderland/native/data/pPb_5_tev_1.root");
    chain->Add("~/Wonderland/native/data/pPb_5_tev_69.root");
    chain->Add("~/Wonderland/native/data/pPb_5_tev_420.root");
    mgr->StartAnalysis("local", chain);
}

//---------------------------------------
void LoadLibraries()
{
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libOADB");
    gSystem->Load("libCDB");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libSTEER");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0rec");
  gSystem->Load("libCORRFW");

  gSystem->Load("libPWGTools");
    gSystem->Load("libPWGHFhfe");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
    gSystem->Load("libPWGCFCorrelationsBase.so");

  gSystem->Load("libEMCALbase.so");
  gSystem->Load("libEMCALUtils.so");
  gSystem->Load("libEMCALrec.so");
  //  gSystem->Load("libPWG4CaloCalib.so");
  gSystem->Load("libPWGCaloTrackCorrBase.so");

  gSystem->Load("libpythia6.so");

  printf("!!!!!!!!!!!!!! loaded all libraries\n\n");
  //    if(use_parFiles)
  //    {
  // //     AliAnalysisAlien::SetupPar("PWGflowBase");
  //      AliAnalysisAlien::SetupPar("PWGflowTasks");
  //    }
}


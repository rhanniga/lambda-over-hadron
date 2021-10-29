#include "AliAnalysisTaskHadronPhiCorr_current.h"

void LoadLibraries();

void RunMacro(){
//     check the function for asymmetric TPC cut in ConfigHFEemcalMod....the rest is still necessary????

   // Firstly, set some variables
   const char* launch = "grid"; // grid, local (if your data is on your local machine, doesn't connect at all)
   const char*  mode = "terminate"; //test, full, terminate  (test= connect to grid but run locally, full= run on grid, terminate= merge output on grid)
   Bool_t pre_final_stage = kFALSE; //TRUE = merging done on grid, FALSE = merge happens locally   
   Int_t cyclenumber = 1;
   Bool_t debug = kTRUE;
   char* work_dir = "PhiCorrelations_LHC16q";
   //char* output_dir = "output_2020_07_30_hphi_alltrig_0_100";
   char* output_dir = "output_2021_08_26_hh_highesttrig_0_20";
   Int_t ttl = 50000;
   Int_t noffiles = 50;
   //Int_t runcycle[]={0,16,31};
   Int_t runcycle[] ={0,31};
   //Int_t runcycle[]={0,1,12,27};
   //Int_t runcycle[]={0,27};
   Bool_t UseParfiles = kFALSE;

   // load libraries
   LoadLibraries();

// create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGCF/Correlations -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGCF/Correlations/Base -I$ALICE_PHYSICS/include -g");
    
    alienHandler->SetAdditionalLibs("AliAnalysisTaskHadronPhiCorr_current.cxx AliAnalysisTaskHadronPhiCorr_current.h AddTaskHadronPhiCorr_current.C libPWGHFhfe.so libCDB.so libSTEER.so libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libPWGTools.so libPWGCFCorrelationsBase.so");
    
  alienHandler->SetAnalysisSource("AliAnalysisTaskHadronPhiCorr_current.cxx");
  //alienHandler->SetOverwriteMode();
  alienHandler->SetRunMode(mode);
  alienHandler->SetNtestFiles(2);
  //alienHandler->SetAPIVersion("V1.1x");
  //alienHandler->SetAliPhysicsVersion("vAN-20191009_ROOT6-1"); //stable version for analysis  
  //alienHandler->SetAliPhysicsVersion("vAN-20190525-1"); //version for test for MC running
  alienHandler->SetAliPhysicsVersion("v5-09-57g-01_ROOT6-1"); //test for LHC16k5b

  //alienHandler->SetFileForTestMode("File_LHC12dPass1.txt");  //txt file that tells where to look for local files if launch=local
  //alienHandler->SetGridDataDir("/alice/sim/2017/LHC17f2b_cent_woSDD/");
  //alienHandler->SetGridDataDir("/alice/sim/2017/LHC17f2b_fast/"); //for MC
//  alienHandler->SetGridDataDir("/alice/sim/2018/LHC18j2_fast/"); //for pp MC
//  alienHandler->SetDataPattern("AOD209/*/*AOD.root"); //for pp MC
//  alienHandler->SetGridDataDir("/alice/sim/2016/LHC16k5a/"); //for pp MC
//  alienHandler->SetDataPattern("AOD/*/*AOD.root"); //for pp MC
  //alienHandler->SetDataPattern("/AOD202/*/*AOD.root"); //for MC
  //alienHandler->SetDataPattern("*ESDs.root");
  //alienHandler->SetGridDataDir("//alice/data/2017/LHC17p/");
  //alienHandler->SetDataPattern("*/pass1_FAST/AOD208/*/*AOD.root"); // for data
  alienHandler->SetGridDataDir("//alice/sim/2020/LHC20f11c_fast/"); // for data
  alienHandler->SetDataPattern("AOD/*/*AOD.root"); // for data
//  alienHandler->SetDataPattern("pass2_FAST/AOD244/*/*AOD.root"); // for data
  //alienHandler->SetDataPattern("*/AOD/*/*AOD.root");
  // alienHandler->SetRunPrefix("000"); // IMPORTANT! Only need for real data, comment this line out for MC data

   
//LHC16r - 8 TeV pPb data
    //Int_t runArray[] = {266318, 266317, 266316, 266305, 266304, 266300, 266299, 266296, 266208, 266197, 266196, 266193, 266190, 266189, 266187, 266117, 266086, 266085, 266084, 266083, 266081, 266076, 266074, 266034, 265797, 265795, 265789, 265788, 265756, 265754, 265746, 265744, 265742, 265741, 265714, 265713,265709, 265705, 265701, 265700, 265698, 265697, 265696, 265607, 265596, 265594};
   //Int_t runArray[] = {266318, 266317, 266316, 266208, 266197, 266196, 266187, 265754, 265744, 265607, 265596, 265594};

//LHC16q - 5 TeV pPb data
   Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};  

//LHC16k5 a&b (a = Pythia 6, b = Pythia 8) - 27 runs
   //Int_t runArray[] = {244340, 244343, 244351, 244355, 244359, 244364, 244377, 244411, 244416, 244418, 244421, 244453, 244456, 244480, 244481, 244482, 244483, 244484, 244531, 244540, 244542, 244617, 244618, 244619, 244626, 244627, 244628};

//LHC18j2 - 5 TeV pp Pythia anchored to LHC17p 
  //Int_t runArray[] = {282008, 282016, 282021, 282025, 282031, 282050, 282051, 282078, 282098, 282099, 282118, 282119, 282120, 282122, 282123, 282125, 282126, 282127, 282146, 282147, 282189, 282206, 282224, 282227, 282229, 282230, 282247, 282302, 282303, 282304, 282305, 282306, 282307, 282309, 282312, 282313, 282314, 282340, 282341, 282342, 282343, 282365, 282366, 282367};
  //Int_t runArray[] = {282008,282031};
  for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
   {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
   }

   printf("\n\nSetting Up alienHandler.\n\n");
   alienHandler->SetGridWorkingDir(work_dir);
   alienHandler->SetGridOutputDir(output_dir);
   alienHandler->SetDefaultOutputs(kTRUE);
   alienHandler->SetAnalysisMacro("PhiInvMass.C");
   alienHandler->SetSplitMaxInputFileNumber(noffiles);
   alienHandler->SetExecutable("PhiInvMass.sh");
   alienHandler->SetExecutableCommand("aliroot -b -q");
   alienHandler->SetTTL(ttl); //10000
   alienHandler->SetInputFormat("xml-single");
   alienHandler->SetJDLName("PhiInvMass.jdl");
   alienHandler->SetPrice(1);
   alienHandler->SetSplitMode("se");
   alienHandler->SetMasterResubmitThreshold(10);
   alienHandler->SetMergeExcludes("EventStat_temp.root");
   alienHandler->SetOutputToRunNo(kTRUE);
   alienHandler->SetKeepLogs(kTRUE);
   alienHandler->SetMaxMergeFiles(15);
   alienHandler->SetMaxMergeStages(7);
   alienHandler->SetMergeViaJDL(pre_final_stage);
//    alienHandler->SetOneStageMerging(kFALSE);   //???????????????????????????????-------------------
    if (!alienHandler) return;

    
// load the necessary macros
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
// Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
   gROOT->ProcessLine(".include $PWD/.");

   gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/base  -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_PHYSICS/PWGCF/Correlations/Base -I$ALICE_PHYSICS/PWGCF/Correlations -g ");
   // gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/macros -I$ALICE_PHYSICS/include -g");

   //printf("\n!!!!!!!!!!!!!!!!!!!!!!\n AliAnalysis Manager \n\n");
   AliAnalysisManager *mgr = new AliAnalysisManager("PhiAnalysis");
   mgr->SetGridHandler(alienHandler);

   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);
//   AliESDInputHandler* esdH = new AliESDInputHandler();
//   mgr->SetInputEventHandler(esdH);

    //AliMCEventHandler* mcH = new AliMCEventHandler();
    //mgr->SetMCtruthEventHandler(mcH);   
//    mcH->SetReadTR(kFALSE);

   //gROOT->LoadMacro("AddTaskPhiCorr.C");
   //gROOT->LoadMacro("./AliAnalysisTaskHadronPhiCorr_current.cxx++g");
   //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
   //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");


    //switch on aliphysicsselection
    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kTRUE, kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C")))); //"kTRUE, kTRUE" for MC, "kFALSE, kTRUE" for Data 
    
    gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")));
        //Only set true for MC
    Bool_t isMC = kTRUE;
    
    gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")));



    TGrid::Connect("alien://");

 
    //TFile* testfile = TFile::Open("/Users/jtblair/alidock/alirepos/utaustin/efficiency/fits_17f2btriggerefficiency.root");
    //TFile* testfile = TFile::Open("alien:///alice/cern.ch/user/j/jblair/PhiEfficiency/fits_17f2b_secondarytest.root");
    TFile* testfile = TFile::Open("alien:///alice/cern.ch/user/j/jblair/PhiEfficiency/fits_17f2bTPC80efficiency.root");
    //TFile* testfile = TFile::Open("/Users/jtblair/alidock/alirepos/utaustin/efficiency/fits_17f2b_secondarytest.root");
    
    //create a task
    
    //AliAnalysisTaskHadronPhiCorr_current *task1 = reinterpret_cast<AliAnalysisTaskHadronPhiCorr_current*>(gInterpreter->ProcessLine(Form(".x %s(kTRUE, 0.0, 20.0)", gSystem->ExpandPathName("AddTaskHadronPhiCorr_current.C"))));
    //gInterpreter->LoadMacro("AliAnalysisTaskHadronPhiCorr_current.cxx++g");
    AliAnalysisTaskHadronPhiCorr_current *task1 = reinterpret_cast<AliAnalysisTaskHadronPhiCorr_current*>(gInterpreter->ProcessLine(Form(".x %s(kTRUE, 0.0, 20.0)", gSystem->ExpandPathName("AddTaskHadronPhiCorr_current.C"))));

    task1->SetIsMCTrue(kFALSE);
    task1->SetIsMCKaon(kFALSE);
    task1->SetIsMCKTrack(kFALSE);
    task1->SetUseAccpt(kTRUE);
    task1->SetSingleTrigger(kFALSE);
    task1->SetSelectTrigger(kFALSE);
    task1->SetHighestTriggerOnly(kTRUE);
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
    
    //task1->GetOutputSlot(1)->GetContainer()->SetNameTitle("phiCorr_mult_0_20", "phiCorr_mult_0_20");
/*
    AliAnalysisTaskHadronPhiCorr_current *task2 = reinterpret_cast<AliAnalysisTaskHadronPhiCorr_current*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE, 20.0, 50.0)", gSystem->ExpandPathName("AddTaskHadronPhiCorr_current.C"))));
    printf("after init task2\n");
    task2->SetIsMCTrue(kFALSE);
    task2->SetIsMCKaon(kFALSE);
    task2->SetIsMCKTrack(kFALSE);
    task2->SetUseAccpt(kFALSE);
    task2->SetSingleTrigger(kFALSE);
    task2->SetKaonEtaCut(0.8);
    task2->SetKaonTPCCut(3.0);
    task2->SetKaonTOFCut(3.0);
    task2->SetTOFVeto(kFALSE);
    task2->SetKaonTrkBit(1024);
    task2->SetAssocTrkBit(1024);
    task2->SetTrigTrkBit(AliAODTrack::kIsHybridGCG);
    task2->SetZVertexMin(-10.0);
    task2->SetZVertexMax(10.0);
    task2->SetZVertexNbins(10);
    task2->SetCentEstimator("V0A");
    task2->SetEtaPhiRegion(0);
    //task2->GetOutputSlot(1)->GetContainer()->SetNameTitle("phiCorr_mult_20_50", "phiCorr_mult_20_50");
*/
 /*
   
    AliAnalysisTaskHadronPhiCorr_current *task3 = reinterpret_cast<AliAnalysisTaskHadronPhiCorr_current*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE, 50.0, 80.0)", gSystem->ExpandPathName("AddTaskHadronPhiCorr_current.C"))));
    task3->SetIsMCTrue(kFALSE);
    task3->SetIsMCKaon(kFALSE);
    task3->SetIsMCKTrack(kFALSE);
    task3->SetUseAccpt(kFALSE);
    task3->SetSingleTrigger(kFALSE);
    task3->SetKaonEtaCut(0.8);
    task3->SetKaonTPCCut(3.0);
    task3->SetKaonTOFCut(3.0);
    task3->SetTOFVeto(kFALSE);
    task3->SetKaonTrkBit(1024);
    task3->SetAssocTrkBit(1024);
    task3->SetTrigTrkBit(AliAODTrack::kIsHybridGCG);
    task3->SetZVertexMin(-10.0);
    task3->SetZVertexMax(10.0);
    task3->SetZVertexNbins(10);
    task3->SetCentEstimator("V0A");
    task3->SetEtaPhiRegion(0);
    //task3->GetOutputSlot(1)->GetContainer()->SetNameTitle("phiCorr_mult_50_80", "phiCorr_mult_50_80");
*/
   
    TF1* phiEff = (TF1*)(testfile->Get("phiFit")->Clone("phiEff"));
    //if(!phiEff){
    //    AliFatal("No phi Eff found!!");
    //}

    TF1* hEff = (TF1*)(testfile->Get("hFit")->Clone("hEff"));
    //if(!hEff){
    //    AliFatal("No h Eff found!!");
    //}
    
    TF1* trigEff = (TF1*)(testfile->Get("trigFit")->Clone("trigEff"));
    //if(!trigEff){
    //    AliFatal("No trig Eff found!!");
    //}

    //task1->LoadEfficiencies(testfile);
    task1->LoadEfficiencies(phiEff, hEff, trigEff);
    //task1->LoadEfficiencies("alien:///alice/cern.ch/user/j/jblair/PhiEfficiency/fits_17f2bCENTTPC80efficiency.root");
    //task2->LoadEfficiencies(testfile);
    //task3->LoadEfficiencies("alien:///alice/cern.ch/user/j/jblair/PhiEfficiency/fits_17f2bCENTTPC80efficiency.root");
    //task3->LoadEfficiencies(testfile);

    //testfile->Close("R");

    if (!mgr->InitAnalysis())
     return;

   printf("mgr initialized analysis\n");
   mgr->PrintStatus();
   mgr->SetUseProgressBar(1, 25);
   fprintf(stdout, "\n!!!!!!!!!!!!!\nAbout to launch analysis... \n");
   // Start analysis in grid.
   mgr->StartAnalysis(launch);
   printf("\n!!!!!!!!!!!!!\nDone with StartAnalysis(launch)\n");
   fflush(stdout);
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
#include <iostream>

#include "AliAODTrack.h"

#include "AliAnalysisTaskLambdaHadronRatio.h"


void runMacro(bool local=true, bool full=true, bool gridMerge=true){

  float MULT_LOW = 0;
  float MULT_HIGH = 80;

  float TRIG_BIT = AliAODTrack::kIsHybridGCG;
  float ASSOC_BIT =  1024; 

// PID CUTS -------------------------------
// narrow cuts:

  // float NSIGMA_TPC_PROTON = 0.6*2;
  // float NSIGMA_TOF_PROTON =0.6*2;
  // float NSIGMA_TPC_PION = 0.6*3;
  // float NSIGMA_TOF_PION = 0.6*3;
  // bool TOF_VETO = true;

// wide cuts:

  // float NSIGMA_TPC_PROTON = 1.4*2;
  // float NSIGMA_TOF_PROTON =1.4*2;
  // float NSIGMA_TPC_PION = 1.4*3;
  // float NSIGMA_TOF_PION = 1.4*3;
  // bool TOF_VETO = true;


// normal cuts:

  float NSIGMA_TPC_PROTON = 2;
  float NSIGMA_TOF_PROTON = 2;
  float NSIGMA_TPC_PION = 3;
  float NSIGMA_TOF_PION = 3;
  bool TOF_VETO = true;

// require TOF hit:

  // float NSIGMA_TPC_PROTON = 1.0*2;
  // float NSIGMA_TOF_PROTON =1.0*2;
  // float NSIGMA_TPC_PION = 1.0*3;
  // float NSIGMA_TOF_PION = 1.0*3;
  // bool TOF_VETO = false;

// ---------------------------------------



  char *EFF_FILE_PATH = "eff_out_mult_eta_dep.root";
  char *CENT_ESTIMATOR = "V0A";

  //Starting and ending index of the array containing the run numbers, specifies which range to run over

  int startIndex = 0;
  int endIndex = 17;

   /* int startIndex = 18; */
  //  int endIndex = 28;

  // int startIndex = 15;
  // int endIndex = 28;

  TString work_dir = "lambda_hadron_ratio_final";
  TString output_dir = "normal_pid";
  
  //If we want to download test files from grid then run in one swoop (usually just run completely locally):
  bool gridTest = false;
  int numTestFiles = 2;

  // So we can access files from the grid (for eff cor and the like)
  // TGrid::Connect("alien//");

  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

  AliAnalysisManager *manage = new AliAnalysisManager("");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  manage->SetInputEventHandler(aodH);


  //MULT SELECTION:
  gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")));

  //SELECTION TASK:
  AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE, kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"))));

  //PID response:
  gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")));

  // Generating task object
  gInterpreter->LoadMacro("AliAnalysisTaskLambdaHadronRatio.cxx++g");
  AliAnalysisTaskLambdaHadronRatio *task = reinterpret_cast<AliAnalysisTaskLambdaHadronRatio*>(gInterpreter->ProcessLine(Form(".x AddLambdaHadronRatioTask.C(\"%s\", %f, %f, %f, %f, %f, %f, %f, %f, %d, \"%s\", \"%s\")",
  "lambdaHadronRatio",
  MULT_LOW,
  MULT_HIGH,
  TRIG_BIT,
  ASSOC_BIT,
  NSIGMA_TPC_PROTON,
  NSIGMA_TOF_PROTON,
  NSIGMA_TPC_PION,
  NSIGMA_TOF_PION,
  TOF_VETO,
  EFF_FILE_PATH,
  CENT_ESTIMATOR)));

  if(!manage->InitAnalysis()) return;
  manage->SetDebugLevel(2);
  manage->PrintStatus();
  manage->SetUseProgressBar(1, 25);

  if(local) {
    TChain *chain = new TChain("aodTree");

    chain->Add("~/Wonderland/native/data/265525_0001.root");  
    chain->Add("~/Wonderland/native/data/265525_0007.root");
    chain->Add("~/Wonderland/native/data/265525_0012.root");
    chain->Add("~/Wonderland/native/data/265525_0018.root");
    chain->Add("~/Wonderland/native/data/265525_0024.root");
    // chain->Add("~/Wonderland/native/data/265525_0030.root");
    // chain->Add("~/Wonderland/native/data/265525_0036.root");
    // chain->Add("~/Wonderland/native/data/265525_0042.root");
    // chain->Add("~/Wonderland/native/data/265525_0048.root");
    // chain->Add("~/Wonderland/native/data/265525_0002.root");
    // chain->Add("~/Wonderland/native/data/265525_0008.root");
    // chain->Add("~/Wonderland/native/data/265525_0013.root");
    // chain->Add("~/Wonderland/native/data/265525_0019.root");
    // chain->Add("~/Wonderland/native/data/265525_0025.root");
    // chain->Add("~/Wonderland/native/data/265525_0031.root");
    // chain->Add("~/Wonderland/native/data/265525_0037.root");
    // chain->Add("~/Wonderland/native/data/265525_0049.root");
    // chain->Add("~/Wonderland/native/data/265525_0003.root");
    // chain->Add("~/Wonderland/native/data/265525_0009.root");
    // chain->Add("~/Wonderland/native/data/265525_0014.root");
    // chain->Add("~/Wonderland/native/data/265525_0020.root");
    // chain->Add("~/Wonderland/native/data/265525_0026.root");
    // chain->Add("~/Wonderland/native/data/265525_0032.root");
    // chain->Add("~/Wonderland/native/data/265525_0038.root");
    // chain->Add("~/Wonderland/native/data/265525_0044.root");
    // chain->Add("~/Wonderland/native/data/265525_0050.root");
    // chain->Add("~/Wonderland/native/data/265525_0004.root");
    // chain->Add("~/Wonderland/native/data/265525_000.root ");
    // chain->Add("~/Wonderland/native/data/265525_0015.root");
    // chain->Add("~/Wonderland/native/data/265525_0021.root");
    // chain->Add("~/Wonderland/native/data/265525_0027.root");
    // chain->Add("~/Wonderland/native/data/265525_0033.root");
    // chain->Add("~/Wonderland/native/data/265525_0039.root");
    // chain->Add("~/Wonderland/native/data/265525_0045.root");
    // chain->Add("~/Wonderland/native/data/265525_0005.root");
    // chain->Add("~/Wonderland/native/data/265525_0010.root");
    // chain->Add("~/Wonderland/native/data/265525_0016.root");
    // chain->Add("~/Wonderland/native/data/265525_0022.root");
    // chain->Add("~/Wonderland/native/data/265525_0028.root");
    // chain->Add("~/Wonderland/native/data/265525_0034.root");
    // chain->Add("~/Wonderland/native/data/265525_0040.root");
    // chain->Add("~/Wonderland/native/data/265525_0046.root");
    // chain->Add("~/Wonderland/native/data/265525_0006.root");
    // chain->Add("~/Wonderland/native/data/265525_0011.root");
    // chain->Add("~/Wonderland/native/data/265525_0017.root");
    // chain->Add("~/Wonderland/native/data/265525_0023.root");
    // chain->Add("~/Wonderland/native/data/265525_0029.root");
    // chain->Add("~/Wonderland/native/data/265525_0035.root");
    // chain->Add("~/Wonderland/native/data/265525_0041.root");
    // chain->Add("~/Wonderland/native/data/265525_0047.root");

    manage->StartAnalysis("local", chain);
  }

  else {
    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskLambdaHadronRatio.cxx AliAnalysisTaskLambdaHadronRatio.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskLambdaHadronRatio.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20181028_ROOT6-1");
    alienHandler->SetAPIVersion("V1.1x");
    // select the input data
    alienHandler->SetGridDataDir("/alice/data/2016/LHC16q/");
    alienHandler->SetDataPattern("pass1_FAST/AOD190/*/*AOD.root");
    // MC has no prefix, data has prefix 000
    alienHandler->SetRunPrefix("000");

    // addding runs
    int runArray[] = {265525, 265521, 265501, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
    int runArrayLength = (int)(sizeof(runArray)/sizeof(runArray[0]));

    /* int problematicRunArray[] = {265383, 265387, 265388, 265427, 265419}; */
    /* int problematicRunArray[] = { 265387}; */
    /* int problematicRunArrayLength = (int)(sizeof(problematicRunArray)/sizeof(problematicRunArray[0])); */

    if(endIndex > (runArrayLength-1) || endIndex < 0) {
            std::cout << "Your end index is out of bounds!" << std::endl;
            std::abort();
    }

    if(startIndex > (runArrayLength-1) || startIndex < 0) {
            std::cout << "Your start index is out of bounds!" << std::endl;
            std::abort();
    }

    /* for(int i = 0; i < problematicRunArrayLength; i++) { */
    /*  alienHandler->AddRunNumber(problematicRunArray[i]); */
    /* } */

    for(int i = startIndex; i < endIndex + 1; i++) {
     alienHandler->AddRunNumber(runArray[i]);
    }

    // number of files per subjob
    alienHandler->SetSplitMaxInputFileNumber(40);
    alienHandler->SetExecutable("LambdaHadronRatio.sh");
    alienHandler->SetJDLName("LambdaHadronRatio.jdl");
    alienHandler->SetTTL(30000);
    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    // merging: run with kTRUE to merge on grid
    // after re-running the jobs in SetRunMode("terminate")
    // (see below) mode, set SetMergeViaJDL(kFALSE)
    // to collect final results
    alienHandler->SetMaxMergeFiles(15);
    alienHandler->SetMaxMergeStages(7);
    alienHandler->SetMergeViaJDL(gridMerge);

    // define the output folders
    alienHandler->SetGridWorkingDir(work_dir); 
    alienHandler->SetGridOutputDir(output_dir);

    // connect the alien plugin to the manager
    manage->SetGridHandler(alienHandler);

    if(gridTest) {
      alienHandler->SetNtestFiles(numTestFiles);
      alienHandler->SetRunMode("test");
      manage->StartAnalysis("grid");
    }
    else {
      if(full) {
        alienHandler->SetRunMode("full");
      }
     else {
        alienHandler->SetRunMode("terminate");
     } 
      manage->StartAnalysis("grid");
    }
  }
}

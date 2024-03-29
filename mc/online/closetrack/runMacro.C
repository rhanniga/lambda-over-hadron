#include "AliAnalysisTaskLambdaHadronCloseTrackEff.h"

void runMacro(bool local=true, bool full=true, bool gridMerge=true){

  float MULT_LOW = 0;
  float MULT_HIGH = 80;

  float TRIG_BIT = AliAODTrack::kIsHybridGCG;
  char *CENT_ESTIMATOR = "V0A";

  //Starting and ending index of the array containing the run numbers, specifies which range to run over
  int startIndex = 0; 
  int endIndex = 28;

  // int startIndex = 15;
  // int endIndex = 28;

  TString work_dir = "lambda_hadron_closetrackeff";
  TString output_dir = "closure_cent_" + std::to_string(int(MULT_LOW)) + "_" + std::to_string(int(MULT_HIGH)) + "_opposite_charges";
  
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
  gInterpreter->LoadMacro("AliAnalysisTaskLambdaHadronCloseTrackEff.cxx++g");
  AliAnalysisTaskLambdaHadronCloseTrackEff *task = reinterpret_cast<AliAnalysisTaskLambdaHadronCloseTrackEff*>(gInterpreter->ProcessLine(Form(".x AddLambdaHadronCloseTrackEffTask.C(\"%s\", %f, %f, %f, \"%s\")",
  "lambdaHadronCloseTrackEff",
  MULT_LOW,
  MULT_HIGH,
  TRIG_BIT,
  CENT_ESTIMATOR)));

  if(!manage->InitAnalysis()) return;
  // manage->SetDebugLevel(2);
  // manage->PrintStatus();
  // manage->SetUseProgressBar(1, 25);

  if(local) {
    TChain *chain = new TChain("aodTree");
    chain->Add("~/Wonderland/native/sim/265525_0001.root");
    chain->Add("~/Wonderland/native/sim/265525_0002.root");
    chain->Add("~/Wonderland/native/sim/265525_0003.root");
    chain->Add("~/Wonderland/native/sim/265525_0004.root");
    chain->Add("~/Wonderland/native/sim/265525_0005.root");
    chain->Add("~/Wonderland/native/sim/265525_0006.root");
    chain->Add("~/Wonderland/native/sim/265525_0007.root");
    chain->Add("~/Wonderland/native/sim/265525_0008.root");
    chain->Add("~/Wonderland/native/sim/265525_0009.root");
    chain->Add("~/Wonderland/native/sim/265525_0010.root");
    manage->StartAnalysis("local", chain);
  }

  else {
    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskLambdaHadronCloseTrackEff.cxx AliAnalysisTaskLambdaHadronCloseTrackEff.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskLambdaHadronCloseTrackEff.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20201026_ROOT6-1");
    alienHandler->SetAPIVersion("V1.1x");
    // select the input data
    alienHandler->SetGridDataDir("//alice/sim/2017/LHC17f2b_fast/");
    alienHandler->SetDataPattern("/AOD202/*/*AOD.root");
    // MC has no prefix, data has prefix 000
    alienHandler->SetRunPrefix("");


    // addding runs
    int runArray[] = {265525, 265521, 265501, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
    int runArrayLength = (int)(sizeof(runArray)/sizeof(runArray[0]));

    if(endIndex > (runArrayLength-1) || endIndex < 0) {
            std::cout << "Your end index is out of bounds!" << std::endl;
            std::abort();
    }

    if(startIndex > (runArrayLength-1) || startIndex < 0) {
            std::cout << "Your start index is out of bounds!" << std::endl;
            std::abort();
    }

    for(int i = startIndex; i < endIndex + 1; i++) {
     alienHandler->AddRunNumber(runArray[i]);
    }

    // number of files per subjob
    alienHandler->SetSplitMaxInputFileNumber(40);
    alienHandler->SetExecutable("LambdaHadronCloseTrackEff.sh");
    alienHandler->SetJDLName("LambdaHadronCloseTrackEff.jdl");
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

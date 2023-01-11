#include "TRint.h"
#include "AliAnalysisTaskLambdaHadronEfficiency.h"


void runMacroPP(bool local=true, bool full=true, bool gridMerge=true){


   //Starting and ending index of the array containing the run numbers, specifies which range to run over
   // PP has 44 runs (20 should be enough for eff calc)
   int startIndex = 0;
   int endIndex = 20;
   char* work_dir = "lambda_hadron_efficiency_pp";
   char* output_dir = "full_stats";

   bool gridTest = false;
   int numTestFiles = 2;

   gInterpreter->ProcessLine(".include $ROOTSYS/include");
   gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

   // create and customize the alien handler
   AliAnalysisManager *manage = new AliAnalysisManager("");
   AliAODInputHandler *aodH = new AliAODInputHandler();
   manage->SetInputEventHandler(aodH);

  //MULT SELECTION:
  gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")));

  //SELECTION TASK:
  AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kTRUE, kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"))));

  //PID response:
  gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")));

  gInterpreter->LoadMacro("AliAnalysisTaskLambdaHadronEfficiency.cxx++g");
  AliAnalysisTaskLambdaHadronEfficiency *task = reinterpret_cast<AliAnalysisTaskLambdaHadronEfficiency*>(gInterpreter->ExecuteMacro("AddLambdaHadronEfficiencyTask.C"));

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


   else{

    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskLambdaHadronEfficiency.cxx AliAnalysisTaskLambdaHadronEfficiency.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskLambdaHadronEfficiency.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20201026_ROOT6-1");
    alienHandler->SetAPIVersion("V1.1x");
    // select the input data
    alienHandler->SetGridDataDir("//alice/sim/2018/LHC18j2_fast/");
    alienHandler->SetDataPattern("/AOD209/*/*AOD.root");

    // addding runs

    int runArray[] = {
      282367,
      282366,
      282365,
      282343,
      282342,
      282341,
      282340,
      282314,
      282313,
      282312,
      282309,
      282307,
      282306,
      282305,
      282304,
      282303,
      282302,
      282247,
      282230,
      282229,
      282227,
      282224,
      282206,
      282189,
      282147,
      282146,
      282127,
      282126,
      282125,
      282123,
      282122,
      282120,
      282119,
      282118,
      282099,
      282098,
      282078,
      282051,
      282050,
      282031,
      282025,
      282021,
      282016,
      282008
    };
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
    alienHandler->SetExecutable("LambdaHadronEfficiency.sh");
    alienHandler->SetJDLName("LambdaHadronEfficiency.jdl");
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

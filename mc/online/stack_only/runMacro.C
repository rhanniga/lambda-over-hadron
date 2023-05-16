#include "AliAnalysisTaskLambdaHadronMC.h"

void runMacro(bool local=true, bool full=true, bool gridMerge=true){


  int startIndex = 0; 
  int endIndex = 28; 

   /* int startIndex = 18; */ 
   /* int endIndex = 28; */ 


  TString work_dir = "lambda_hadron_mc";
  TString output_dir = "epos";
  
  //If we want to download test files from grid then run in one swoop (usually just run completely locally):
  bool gridTest = false;
  int numTestFiles = 1;

  // So we can access files from the grid (for eff cor and the like)
  // TGrid::Connect("alien//");

  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

  AliAnalysisManager *manage = new AliAnalysisManager("");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  manage->SetInputEventHandler(aodH);


  // Generating task object
  gInterpreter->LoadMacro("AliAnalysisTaskLambdaHadronMC.cxx++g");
  AliAnalysisTaskLambdaHadronMC *task = reinterpret_cast<AliAnalysisTaskLambdaHadronMC*>(gInterpreter->ProcessLine(Form(".x AddLambdaHadronMCTask.C(\"%s\")",
  "lambdaHadronMC")));

  if(!manage->InitAnalysis()) return;
  manage->SetDebugLevel(2);
  manage->PrintStatus();
  manage->SetUseProgressBar(1, 25);

  if(local) {
    TChain *chain = new TChain("aodTree");
    chain->Add("~/Wonderland/native/sim/265525_0001.root");
    chain->Add("~/Wonderland/native/sim/265525_0002.root");
    chain->Add("~/Wonderland/native/sim/265525_0003.root");
    chain->Add("~/Wonderland/native/sim/265525_0004.root");
    // chain->Add("~/Wonderland/native/sim/265525_0005.root");
    // chain->Add("~/Wonderland/native/sim/265525_0006.root");
    // chain->Add("~/Wonderland/native/sim/265525_0007.root");
    // chain->Add("~/Wonderland/native/sim/265525_0008.root");
    // chain->Add("~/Wonderland/native/sim/265525_0009.root");
    // chain->Add("~/Wonderland/native/sim/265525_0010.root");
    manage->StartAnalysis("local", chain);
  }

  else {
    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskLambdaHadronMC.cxx AliAnalysisTaskLambdaHadronMC.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskLambdaHadronMC.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20230210_O2-1");
    alienHandler->SetAPIVersion("V1.1x");
    // select the input data (these are all the runs anchored to 16q

    // EPOS-LHC set:
    alienHandler->SetGridDataDir("//alice/sim/2017/LHC17f2a_fast_fix/");
    alienHandler->SetDataPattern("/AOD228/*/*AOD.root");

    // DPMJet set: (small)
    // alienHandler->SetGridDataDir("//alice/sim/2017/LHC17f2b_fast/");
    // alienHandler->SetDataPattern("/AOD202/*/*AOD.root");

    // alienHandler->SetGridDataDir("//alice/sim/2017/LHC17f2b_cent_woSDD/");
    /* alienHandler->SetGridDataDir("//alice/sim/2018/LHC18f3_fast_2/"); */
    // alienHandler->SetGridDataDir("//alice/sim/2018/LHC18f3_cent_woSDD_2/");
    // alienHandler->SetGridDataDir("//alice/sim/2018/LHC18j2_fast/");
    // alienHandler->SetDataPattern("/AOD209/*/*AOD.root");

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
    alienHandler->SetExecutable("LambdaHadronMC.sh");
    alienHandler->SetJDLName("LambdaHadronMC.jdl");
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
#include <iostream>

#include "AliAODTrack.h"

#include "AliAnalysisTaskLambdaHadronRatio.h"


void runOffline(int numRuns){

  float MULT_LOW = 0;
  float MULT_HIGH = 20;

  float TRIG_BIT = AliAODTrack::kIsHybridGCG;
  float ASSOC_BIT =  1024; 
  char *EFF_FILE_PATH = "eff_out.root";
  char *CENT_ESTIMATOR = "V0A";

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
  AliAnalysisTaskLambdaHadronRatio *task = reinterpret_cast<AliAnalysisTaskLambdaHadronRatio*>(gInterpreter->ProcessLine(Form(".x AddLambdaHadronRatioTask.C(\"%s\", %f, %f, %f, %f, \"%s\", \"%s\")",
  "lambdaHadronRatio",
  MULT_LOW,
  MULT_HIGH,
  TRIG_BIT,
  ASSOC_BIT,
  EFF_FILE_PATH,
  CENT_ESTIMATOR)));

  if(!manage->InitAnalysis()) return;
  manage->SetDebugLevel(0);
  manage->PrintStatus();

  TChain *chain = new TChain("aodTree");
  for(int i = 0; i < numRuns; i++) {
    TString inputFile;
    if (i < 9) {
      inputFile = "~/Wonderland/native/data/265525_000";
      inputFile += std::to_string(i+1);
      inputFile += ".root";
    }
    else {
      inputFile = "~/Wonderland/native/data/265525_00";
      inputFile += std::to_string(i+1);
      inputFile += ".root";
    }
    chain->Add(inputFile);
  }
  manage->StartAnalysis("local", chain);
}
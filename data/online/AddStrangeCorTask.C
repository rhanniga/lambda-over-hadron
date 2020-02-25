#include <iostream>

AliAnalysisTaskStrangeCor* AddStrangeCorTask(TString name = "name") {

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();

  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskStrangeCor* task = new AliAnalysisTaskStrangeCor(name.Data());

  if(!task) return 0x0;

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("h-lambda", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}

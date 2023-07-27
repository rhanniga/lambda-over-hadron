#include <iostream>

AliAnalysisTaskNSigmaElectron* AddNSigmaElectronTask(TString name = "nSigmaElectron") {

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();

  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskNSigmaElectron* task = new AliAnalysisTaskNSigmaElectron(name.Data());

  if(!task) return 0x0;

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("nsigma-electron", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}

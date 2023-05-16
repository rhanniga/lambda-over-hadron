AliAnalysisTaskLambdaHadronMC* AddLambdaHadronMCTask(
  TString name = "lambdaHadronMC"
  ) {

  // NOTE: The default arguments are placeholders ONLY, everything should be set within the run macro before function is called
  // 1024 is primary tracks (tight DCA cut, BIT(10))
  // 1048576 is kIsHybridGCG (BIT(20))

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();
  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskLambdaHadronMC* task = new AliAnalysisTaskLambdaHadronMC(name.Data());

  if(!task) return 0x0;

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("h-lambda", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}
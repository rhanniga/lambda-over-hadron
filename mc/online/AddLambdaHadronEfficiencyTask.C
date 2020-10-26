AliAnalysisTaskLambdaHadronEfficiency *AddLambdaHadronEfficiencyTask(Float_t multLow = 0.0, Float_t multHigh = 20.0){
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    if (!mgr) return 0x0;

    if (!mgr->GetInputEventHandler()) return 0x0;

    TString file_name = AliAnalysisManager::GetCommonFileName();
    TString task_name = "h-lambda_eff";

    AliAnalysisTaskLambdaHadronEfficiency *lambdaEff = new AliAnalysisTaskLambdaHadronEfficiency(task_name.Data(), multLow, multHigh); 

    lambdaEff->SelectCollisionCandidates(AliVEvent::kINT7);


    mgr->AddTask(lambdaEff);

    mgr->ConnectInput(lambdaEff, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(lambdaEff, 1, mgr->CreateContainer("h-lambda_eff", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

    return lambdaEff;

}

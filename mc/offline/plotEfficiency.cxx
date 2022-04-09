void plotMultEff(THnSparse* reco, THnSparse* real, TH1D** recoVar, TH1D** realVar, TH1D** eff, TH1D** ratio, Float_t* mult, Int_t numMultBins, Int_t axis, TCanvas* ceff, TCanvas* cratio, Bool_t isSingle = kFALSE){

    TString var = "";
    TString label = "";
    Int_t rebin = 1;
    Int_t numbins = 10;
    Double_t binedge[11] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0};

    Int_t colors[] = {kGreen+4,kGreen-2,kSpring+5,kRed+1};
    
    Int_t multAxis = 5;
    if(isSingle) multAxis = 4;

    printf("num mult bins: %i\n", numMultBins);
    switch(axis){
        case 0: var = "PT";
                label = "p_{T}";
                //rebin = 5;
                break;
        case 1: var = "Phi";
                label = "#varphi";
                rebin = 4;
                break;
        case 2: var = "Eta";
                label = "#eta";
                rebin = 2;
                break;
        case 3: var = "Z";
                label = "Z Vertex (cm)";
                break;
        default: var = "Var";
                 label = "var";
                 break;
    }

    for(int imult = 0; imult < numMultBins; imult++){
        reco->GetAxis(multAxis)->SetRangeUser(mult[imult], mult[imult+1]);
        real->GetAxis(multAxis)->SetRangeUser(mult[imult], mult[imult+1]);
        recoVar[imult] = reco->Projection(axis);
        recoVar[imult]->SetStats(kFALSE);
        recoVar[imult]->SetName(Form("%sreco%s_mult%i%i", reco->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        realVar[imult] = real->Projection(axis);
        realVar[imult]->SetName(Form("%sreco%s_mult%i%i", real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        if(axis!=0){
            recoVar[imult]->Rebin(rebin);
            realVar[imult]->Rebin(rebin);
        }else{
            recoVar[imult]= (TH1D*)recoVar[imult]->Rebin(numbins, Form("%sreco_rebin%s_mult%i%i", reco->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])), binedge);
            realVar[imult]= (TH1D*)realVar[imult]->Rebin(numbins, Form("%sreal_rebin%s_mult%i%i", real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])), binedge);
        } 
        eff[imult] = (TH1D*)recoVar[imult]->Clone(Form("%s_%s_eff%s_mult%i%i", reco->GetTitle(), real->GetTitle(), var.Data(), int(mult[imult]), int(mult[imult+1])));
        eff[imult]->Divide(eff[imult], realVar[imult], 1., 1., "B");
        eff[imult]->SetTitle(Form("(%s)/(%s) Efficiency vs. %s for Mult. Bins;%s;Efficiency", reco->GetTitle(), real->GetTitle(), label.Data(), label.Data()));
        eff[imult]->SetLineColor(colors[imult]);
        eff[imult]->SetMarkerColor(colors[imult]);
        eff[imult]->SetMarkerSize(1);
        eff[imult]->SetMarkerStyle(22);
        eff[imult]->GetYaxis()->SetRangeUser(0.0, eff[imult]->GetBinContent(eff[imult]->GetMaximumBin())*1.25);
    }
 
    reco->GetAxis(multAxis)->SetRangeUser(0.0, 90.0);
    real->GetAxis(multAxis)->SetRangeUser(0.0, 90.0);
    recoVar[numMultBins] = reco->Projection(axis);
    recoVar[numMultBins]->SetName(Form("%sreco%s_mult090", reco->GetTitle(), var.Data()));
    realVar[numMultBins] = real->Projection(axis);
    realVar[numMultBins]->SetName(Form("%sreal%s_mult090", real->GetTitle(), var.Data()));
    if(axis!=0){
        recoVar[numMultBins]->Rebin(rebin);
        realVar[numMultBins]->Rebin(rebin);
    }else{
        recoVar[numMultBins]= (TH1D*)recoVar[numMultBins]->Rebin(numbins, Form("%sreco_rebin%s_mult090", reco->GetTitle(), var.Data()), binedge);
        realVar[numMultBins]= (TH1D*)realVar[numMultBins]->Rebin(numbins, Form("%sreal_rebin%s_mult090", real->GetTitle(), var.Data()), binedge);

    }
    eff[numMultBins] = (TH1D*)recoVar[numMultBins]->Clone(Form("%s_%s_eff%s_mult090",reco->GetTitle(), real->GetTitle(), var.Data()));
    eff[numMultBins]->Divide(eff[numMultBins], realVar[numMultBins], 1., 1., "B");
    eff[numMultBins]->SetTitle(Form("(%s)/(%s) Efficiency vs. %s for Mult. Bins", reco->GetTitle(), real->GetTitle(), label.Data()));
    eff[numMultBins]->SetLineColor(kRed+1);
    eff[numMultBins]->SetMarkerColor(kRed+1);
    eff[numMultBins]->SetMarkerSize(1);
    eff[numMultBins]->SetMarkerStyle(22);

    cratio->SetLeftMargin(0.13);
    cratio->cd();
    for(int i=0; i<numMultBins; i++){
        ratio[i] = (TH1D*)eff[i]->Clone(Form("%s_%s_ratio%s_mult%i", reco->GetTitle(), real->GetTitle(), var.Data(), i));
        ratio[i]->SetTitle(Form("(%s)/(%s) Efficiency vs %s Ratio to 0-100%% Efficiency;%s;Ratio", reco->GetTitle(), real->GetTitle(), label.Data(), label.Data()));
        ratio[i]->Divide(ratio[i], eff[3], 1., 1., "B");
        ratio[i]->GetYaxis()->SetRangeUser(0.75, 1.25);
        ratio[i]->Draw("P SAME");
    } 
    //cratio->Print(Form("plots/binom_%s_%s_%sratio.pdf", var.Data(), reco->GetName(), real->GetName()));

    ceff->SetLeftMargin(0.13);
    ceff->cd();
    for(int j=0; j<=numMultBins; j++){
        eff[j]->Draw("P SAME");
    }
    //ceff->Print(Form("plots/binom_%s_%s_%seff.pdf", var.Data(), reco->GetName(), real->GetName()));
    TLegend* effleg = new TLegend(0.2, 0.6, 0.6, 0.8);
    effleg->AddEntry(eff[3], "0-100% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[0], "0-20% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[1], "20-50% Mult. Percentile", "lpe");
    effleg->AddEntry(eff[2], "50-80% Mult. Percentile", "lpe");
    effleg->Draw();

}; 

void plotEfficiency(){

    Float_t mult[4] = {0.0, 20.0, 50.0, 80.0};
    
    TFile* infile = new TFile("~/OneDrive/Research/Output/lambda-over-hadron/mc/efficiency_run_etacutonline.root");
    TList* list = (TList*)infile->Get("h-lambda_eff");
    
    float SIG_MIN = 1.08;
    float SIG_MAX = 1.16 - 0.00000001;
    
    //single track histos
    THnSparseF* realTrigger = (THnSparseF*)list->FindObject("fRealChargedDist");
    THnSparseF* recoTrigger = (THnSparseF*)list->FindObject("fRecoChargedTriggerDist");
    THnSparseF* realCharged = (THnSparseF*)list->FindObject("fRealChargedDist");
    THnSparseF* recoCharged = (THnSparseF*)list->FindObject("fRecoChargedDist");

    realTrigger->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoTrigger->GetAxis(0)->SetRangeUser(0.5, 10.0);
    realCharged->GetAxis(0)->SetRangeUser(0.5, 10.0);
    recoCharged->GetAxis(0)->SetRangeUser(0.5, 10.0);

    // offline eta cuts no longer required
    //recoTrigger->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    // realTrigger->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    //recoCharged->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    // realCharged->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    realTrigger->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    recoTrigger->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    recoCharged->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    realCharged->GetAxis(3)->SetRangeUser(-10.0, 10.0);

    realTrigger->Sumw2();
    recoTrigger->Sumw2();
    recoCharged->Sumw2();
    realCharged->Sumw2();
    
    TH1D* realTrigger_PT_mult[4];
    TH1D* recoTrigger_PT_mult[4];
    TH1D* recoCharged_PT_mult[4];
    TH1D* realCharged_PT_mult[4];

    TH1D* effTrigger_PT_mult[4];
    TH1D* effCharged_PT_mult[4];

    TH1D* ratioTrigger_PT_mult[3];
    TH1D* ratioCharged_PT_mult[3];

    TCanvas* ceffTrigger_PT = new TCanvas("ceffTrigger_PT", "effTrigger_PT", 50, 50, 600, 600);
    TCanvas* ceffCharged_PT = new TCanvas("ceffCharged_PT", "effCharged_PT", 50, 50, 600, 600);

    TCanvas* cratioTrigger_PT = new TCanvas("cratioTrigger_PT", "ratioTrigger_PT", 50, 50, 600, 600);
    TCanvas* cratioCharged_PT = new TCanvas("cratioCharged_PT", "ratioCharged_PT", 50, 50, 600, 600);
    plotMultEff(recoTrigger, realTrigger, recoTrigger_PT_mult, realTrigger_PT_mult, effTrigger_PT_mult, ratioTrigger_PT_mult, mult, 3, 0, ceffTrigger_PT, cratioTrigger_PT, kTRUE);
    plotMultEff(recoCharged, realCharged, recoCharged_PT_mult, realCharged_PT_mult, effCharged_PT_mult, ratioCharged_PT_mult, mult, 3, 0, ceffCharged_PT, cratioCharged_PT, kTRUE);

    // real lambda histo
    THnSparseF* realLambda = (THnSparseF*)list->FindObject("fRealTotalLambdaDist");

    //reco lambda histos
    THnSparseF* recoLambda = (THnSparseF*)list->FindObject("fRecoTotalLambdaDist");
    THnSparseF* etaLambda = (THnSparseF*)list->FindObject("fRecoEtaLambdaDist");
    THnSparseF* etaPtLambda = (THnSparseF*)list->FindObject("fRecoEtaPtLambdaDist");
    THnSparseF* etaPtRefitLambda = (THnSparseF*)list->FindObject("fRecoEtaPtRefitLambdaDist");
    THnSparseF* etaPtRefitRowsLambda = (THnSparseF*)list->FindObject("fRecoEtaPtRefitRowsLambdaDist");
    THnSparseF* etaPtRefitRowsRatioLambda = (THnSparseF*)list->FindObject("fRecoEtaPtRefitRowsRatioLambdaDist");

    //reco lambda from v0 finder histos
    THnSparseF* recoLambdaV0 = (THnSparseF*)list->FindObject("fRecoTotalV0LambdaDist");
    THnSparseF* etaLambdaV0 = (THnSparseF*)list->FindObject("fRecoEtaV0LambdaDist");
    THnSparseF* etaPtLambdaV0 = (THnSparseF*)list->FindObject("fRecoEtaPtV0LambdaDist");
    THnSparseF* etaPtRefitLambdaV0 = (THnSparseF*)list->FindObject("fRecoEtaPtRefitV0LambdaDist");
    THnSparseF* etaPtRefitRowsLambdaV0 = (THnSparseF*)list->FindObject("fRecoEtaPtRefitRowsV0LambdaDist");
    THnSparseF* etaPtRefitRowsRatioLambdaV0 = (THnSparseF*)list->FindObject("fRecoEtaPtRefitRowsRatioV0LambdaDist");


    realLambda->SetTitle("real #Lambda");
    realLambda->SetName("realLambda");

    recoLambda->SetTitle("reco #Lambda");
    recoLambda->SetName("recoLambda");
    etaLambda->SetTitle("reco #Lambda (eta <0.8 on daughters)");
    etaLambda->SetName("etaLambda");
    etaPtLambda->SetTitle("reco #Lambda (eta <0.8, pt > 0.15 on daughters)");
    etaPtLambda->SetName("etaPtLambda");
    etaPtRefitLambda->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit on daughters)");
    etaPtRefitLambda->SetName("etaPtRefitLambda");
    etaPtRefitRowsLambda->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit, ncrossedTPC > 70 on daughters)");
    etaPtRefitRowsLambda->SetName("etaPtRefitRowsLambda");
    etaPtRefitRowsRatioLambda->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit, ncrossedTPC > 70 on daughters)");
    etaPtRefitRowsRatioLambda->SetName("etaPtRefitRowsRatioLambda");

    recoLambdaV0->SetTitle("reco #Lambda from v0");
    recoLambdaV0->SetName("recoLambdaV0");
    etaLambdaV0->SetTitle("reco #Lambda (eta <0.8 on daughters) from v0");
    etaLambdaV0->SetName("etaLambdaV0");
    etaPtLambdaV0->SetTitle("reco #Lambda (eta <0.8, pt > 0.15 on daughters) from v0");
    etaPtLambdaV0->SetName("etaPtLambdaV0");
    etaPtRefitLambdaV0->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit on daughters) from v0");
    etaPtRefitLambdaV0->SetName("etaPtRefitLambdaV0");
    etaPtRefitRowsLambdaV0->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit, ncrossedTPC > 70 on daughters) from v0");
    etaPtRefitRowsLambdaV0->SetName("etaPtRefitRowsLambdaV0");
    etaPtRefitRowsRatioLambdaV0->SetTitle("reco #Lambda (eta <0.8, pt > 0.15, tpc refit, ncrossedTPC > 70 on daughters) from v0");
    etaPtRefitRowsRatioLambdaV0->SetName("etaPtRefitRowsRatioLambdaV0");

    // PT AXIS  
    realLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);

    recoLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitRowsLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitRowsRatioLambda->GetAxis(0)->SetRangeUser(0.5, 10.0);

    recoLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitRowsLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
    etaPtRefitRowsRatioLambdaV0->GetAxis(0)->SetRangeUser(0.5, 10.0);
       
    // ETA AXIS  
    realLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    recoLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitRowsLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitRowsRatioLambda->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    recoLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitRowsLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);
    etaPtRefitRowsRatioLambdaV0->GetAxis(2)->SetRangeUser(-0.8, 0.8);

    // Z VTX AXIS  
    realLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);

    recoLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitRowsLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitRowsRatioLambda->GetAxis(3)->SetRangeUser(-10.0, 10.0);

    recoLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitRowsLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);
    etaPtRefitRowsRatioLambdaV0->GetAxis(3)->SetRangeUser(-10.0, 10.0);


    // MASS AXIS  
    realLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);

    recoLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitRowsLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitRowsRatioLambda->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);

    recoLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitRowsLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);
    etaPtRefitRowsRatioLambdaV0->GetAxis(4)->SetRangeUser(SIG_MIN, SIG_MAX);


    realLambda->Sumw2();

    recoLambda->Sumw2();
    etaLambda->Sumw2();
    etaPtLambda->Sumw2();
    etaPtRefitLambda->Sumw2();
    etaPtRefitRowsLambda->Sumw2();
    etaPtRefitRowsRatioLambda->Sumw2();

    recoLambdaV0->Sumw2();
    etaLambdaV0->Sumw2();
    etaPtLambdaV0->Sumw2();
    etaPtRefitLambdaV0->Sumw2();
    etaPtRefitRowsLambdaV0->Sumw2();
    etaPtRefitRowsRatioLambdaV0->Sumw2();

    // pT vs mult
    TH1D* realLambda_PT_mult[4];

    TH1D* recoLambda_PT_mult[4];
    TH1D* etaLambda_PT_mult[4];
    TH1D* etaPtLambda_PT_mult[4];
    TH1D* etaPtRefitLambda_PT_mult[4];
    TH1D* etaPtRefitRowsLambda_PT_mult[4];
    TH1D* etaPtRefitRowsRatioLambda_PT_mult[4];

    TH1D* recoLambdaV0_PT_mult[4];
    TH1D* etaLambdaV0_PT_mult[4];
    TH1D* etaPtLambdaV0_PT_mult[4];
    TH1D* etaPtRefitLambdaV0_PT_mult[4];
    TH1D* etaPtRefitRowsLambdaV0_PT_mult[4];
    TH1D* etaPtRefitRowsRatioLambdaV0_PT_mult[4];

    TH1D* effPT_mult[4];
    TH1D* etaeffPT_mult[4];
    TH1D* etaPteffPT_mult[4];
    TH1D* etaPtRefiteffPT_mult[4];
    TH1D* etaPtRefitRowseffPT_mult[4];
    TH1D* etaPtRefitRowsRatioeffPT_mult[4];

    TH1D* effV0PT_mult[4];
    TH1D* etaeffV0PT_mult[4];
    TH1D* etaPteffV0PT_mult[4];
    TH1D* etaPtRefiteffV0PT_mult[4];
    TH1D* etaPtRefitRowseffV0PT_mult[4];
    TH1D* etaPtRefitRowsRatioeffV0PT_mult[4];

    TH1D* ratioPT_mult[3];
    TH1D* etaratioPT_mult[3];
    TH1D* etaPtratioPT_mult[3];
    TH1D* etaPtRefitratioPT_mult[3];
    TH1D* etaPtRefitRowsratioPT_mult[3];
    TH1D* etaPtRefitRowsRatioratioPT_mult[3];

    TH1D* ratioV0PT_mult[3];
    TH1D* etaratioV0PT_mult[3];
    TH1D* etaPtratioV0PT_mult[3];
    TH1D* etaPtRefitratioV0PT_mult[3];
    TH1D* etaPtRefitRowsratioV0PT_mult[3];
    TH1D* etaPtRefitRowsRatioratioV0PT_mult[3];
     
    TCanvas* cratioPT = new TCanvas("cratioPT", "cratioPT", 55, 55, 600, 600);
    TCanvas* cetaratioPT = new TCanvas("cetaratioPT", "cetaratioPT", 55, 55, 600, 600);
    TCanvas* cetaPtratioPT = new TCanvas("cetaPtratioPT", "cetaPtratioPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitratioPT = new TCanvas("cetaPtRefitratioPT", "cetaPtRefitratioPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsratioPT = new TCanvas("cetaPtRefitRowsratioPT", "cetaPtRefitRowsratioPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsRatioratioPT = new TCanvas("cetaPtRefitRowsRatioratioPT", "cetaPtRefitRowsRatioratioPT", 55, 55, 600, 600);
     
    TCanvas* cratioV0PT = new TCanvas("cratioV0PT", "cratioV0PT", 55, 55, 600, 600);
    TCanvas* cetaratioV0PT = new TCanvas("cetaratioV0PT", "cetaratioV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtratioV0PT = new TCanvas("cetaPtratioV0PT", "cetaPtratioV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitratioV0PT = new TCanvas("cetaPtRefitratioV0PT", "cetaPtRefitratioV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsratioV0PT = new TCanvas("cetaPtRefitRowsratioV0PT", "cetaPtRefitRowsratioV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsRatioratioV0PT = new TCanvas("cetaPtRefitRowsRatioratioV0PT", "cetaPtRefitRowsRatioratioV0PT", 55, 55, 600, 600);

    TCanvas* ceffPT = new TCanvas("ceffPT", "ceffPT", 55, 55, 600, 600);
    TCanvas* cetaeffPT = new TCanvas("cetaeffPT", "cetaeffPT", 55, 55, 600, 600);
    TCanvas* cetaPteffPT = new TCanvas("cetaPteffPT", "cetaPteffPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefiteffPT = new TCanvas("cetaPtRefiteffPT", "cetaPtRefiteffPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowseffPT = new TCanvas("cetaPtRefitRowseffPT", "cetaPtRefitRowseffPT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsRatioeffPT = new TCanvas("cetaPtRefitRowsRatioeffPT", "cetaPtRefitRowsRatioeffPT", 55, 55, 600, 600);

    TCanvas* ceffV0PT = new TCanvas("ceffV0PT", "ceffV0PT", 55, 55, 600, 600);
    TCanvas* cetaeffV0PT = new TCanvas("cetaeffV0PT", "cetaeffV0PT", 55, 55, 600, 600);
    TCanvas* cetaPteffV0PT = new TCanvas("cetaPteffV0PT", "cetaPteffV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefiteffV0PT = new TCanvas("cetaPtRefiteffV0PT", "cetaPtRefiteffV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowseffV0PT = new TCanvas("cetaPtRefitRowseffV0PT", "cetaPtRefitRowseffV0PT", 55, 55, 600, 600);
    TCanvas* cetaPtRefitRowsRatioeffV0PT = new TCanvas("cetaPtRefitRowsRatioeffV0PT", "cetaPtRefitRowsRatioeffV0PT", 55, 55, 600, 600);

    plotMultEff(recoLambda, realLambda, recoLambda_PT_mult, realLambda_PT_mult, effPT_mult, ratioPT_mult, mult, 3, 0, ceffPT, cratioPT, kFALSE);
    plotMultEff(etaLambda, realLambda, etaLambda_PT_mult, realLambda_PT_mult, etaeffPT_mult, etaratioPT_mult, mult, 3, 0, cetaeffPT, cetaratioPT, kFALSE);
    plotMultEff(etaPtLambda, realLambda, etaPtLambda_PT_mult, realLambda_PT_mult, etaPteffPT_mult, etaPtratioPT_mult, mult, 3, 0, cetaPteffPT, cetaPtratioPT, kFALSE);
    plotMultEff(etaPtRefitLambda, realLambda, etaPtRefitLambda_PT_mult, realLambda_PT_mult, etaPtRefiteffPT_mult, etaPtRefitratioPT_mult, mult, 3, 0, cetaPtRefiteffPT, cetaPtRefitratioPT, kFALSE);
    plotMultEff(etaPtRefitRowsLambda, realLambda, etaPtRefitRowsLambda_PT_mult, realLambda_PT_mult, etaPtRefitRowseffPT_mult, etaPtRefitRowsratioPT_mult, mult, 3, 0, cetaPtRefitRowseffPT, cetaPtRefitRowsratioPT, kFALSE);
    plotMultEff(etaPtRefitRowsRatioLambda, realLambda, etaPtRefitRowsRatioLambda_PT_mult, realLambda_PT_mult, etaPtRefitRowsRatioeffPT_mult, etaPtRefitRowsRatioratioPT_mult, mult, 3, 0, cetaPtRefitRowsRatioeffPT, cetaPtRefitRowsRatioratioPT, kFALSE);

    plotMultEff(recoLambdaV0, realLambda, recoLambdaV0_PT_mult, realLambda_PT_mult, effV0PT_mult, ratioV0PT_mult, mult, 3, 0, ceffV0PT, cratioV0PT, kFALSE);
    plotMultEff(etaLambdaV0, realLambda, etaLambdaV0_PT_mult, realLambda_PT_mult, etaeffV0PT_mult, etaratioV0PT_mult, mult, 3, 0, cetaeffV0PT, cetaratioV0PT, kFALSE);
    plotMultEff(etaPtLambdaV0, realLambda, etaPtLambdaV0_PT_mult, realLambda_PT_mult, etaPteffV0PT_mult, etaPtratioV0PT_mult, mult, 3, 0, cetaPteffV0PT, cetaPtratioV0PT, kFALSE);
    plotMultEff(etaPtRefitLambdaV0, realLambda, etaPtRefitLambdaV0_PT_mult, realLambda_PT_mult, etaPtRefiteffV0PT_mult, etaPtRefitratioV0PT_mult, mult, 3, 0, cetaPtRefiteffV0PT, cetaPtRefitratioV0PT, kFALSE);
    plotMultEff(etaPtRefitRowsLambdaV0, realLambda, etaPtRefitRowsLambdaV0_PT_mult, realLambda_PT_mult, etaPtRefitRowseffV0PT_mult, etaPtRefitRowsratioV0PT_mult, mult, 3, 0, cetaPtRefitRowseffV0PT, cetaPtRefitRowsratioV0PT, kFALSE);
    plotMultEff(etaPtRefitRowsRatioLambdaV0, realLambda, etaPtRefitRowsRatioLambdaV0_PT_mult, realLambda_PT_mult, etaPtRefitRowsRatioeffV0PT_mult, etaPtRefitRowsRatioratioV0PT_mult, mult, 3, 0, cetaPtRefitRowsRatioeffV0PT, cetaPtRefitRowsRatioratioV0PT, kFALSE);

    TFile* output = new TFile("eff_out.root", "RECREATE");

    effCharged_PT_mult[3]->SetName("fAssociatedEff");
    effCharged_PT_mult[3]->Write();

    effTrigger_PT_mult[3]->SetName("fTriggerEff");
    effTrigger_PT_mult[3]->Write();

    etaPtRefitRowsRatioeffPT_mult[3]->SetName("fLambdaEff");
    etaPtRefitRowsRatioeffPT_mult[3]->Write();

    etaPtRefitRowsRatioeffV0PT_mult[3]->SetName("fLambdaV0Eff");
    etaPtRefitRowsRatioeffV0PT_mult[3]->Write();

}

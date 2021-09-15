TH2D* makehhCorrections(TH2D* same2D, TH2D* mix2D){
    
    same2D->GetYaxis()->SetRange(1,same2D->GetYaxis()->GetNbins());
    mix2D->GetYaxis()->SetRange(1, mix2D->GetYaxis()->GetNbins());

    Float_t scale = 0.0;

    TH2D* corr2D;

    scale = 0.5*(float)(mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(0.01), mix2D->GetYaxis()->FindBin(0.0)) + mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(-0.01), mix2D->GetYaxis()->FindBin(0.0)));
    printf("scale: %e \n", scale);
    
    corr2D = (TH2D*)same2D->Clone("corr2D");
    corr2D->Divide(mix2D);
    corr2D->Scale(scale);

    return corr2D;
}
//-------------------------------------------------------------------------------------------
TH2D* projectWithEfficiencyCorrections(THnSparseF* sparse, TH1D* eff, float assocPTLow, float assocPTHigh){
    Int_t lowbin = eff->GetXaxis()->FindBin(assocPTLow + 0.00001);
    Int_t highbin = eff->GetXaxis()->FindBin(assocPTHigh - 0.00001);
    Int_t axes[] = {2,3};
    TH2D* corr;
    TH2D* buff;
    for(Int_t i = lowbin; i <= highbin; i++){
        sparse->GetAxis(1)->SetRange(sparse->GetAxis(1)->FindBin(eff->GetBinCenter(i)), sparse->GetAxis(1)->FindBin(eff->GetBinCenter(i)));
        buff = (TH2D*)sparse->Projection(2,3);
        buff->Sumw2();
        buff->Scale(1.0/eff->GetBinContent(i));
        if(i==lowbin){
            corr = (TH2D*)buff->Clone(Form("effcorr%s", sparse->GetName()));
        }else{
            corr->Add(buff);
        }
    }
    return corr;
}
//--------------------------------------------------------------------------------------------
void makeHHMixCorrectionsZVertexEff(string inputName, int multLow, int multHigh, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh, string suffix="_"){
    //TFile *effFile = new TFile("~/utaustin/efficiency/17f2befficiency.root");
    //TH1D* hadronEff = (TH1D*)effFile->Get("hadronPTEff");

    TFile *histoFile = new TFile(inputName.c_str());
    //string mult = inputName.substr(inputName.find("_", inputName.find("_")+1), inputName.find(".") - inputName.find("_", inputName.find("_")+1));
    string mult = "_" + std::to_string(multLow) + "_" + std::to_string(multHigh);
    TList* list = (TList*) histoFile->Get("hhCorr_mult_0_20_");
   
    TH2D *trigHHDist = (TH2D*)list->FindObject("fTrigHHDist");

    Float_t epsilon = 0.0001;

    THnSparseF* trigDist = (THnSparseF*)list->FindObject("fTrigDist");
    TH1D* trigDist1D = (TH1D*)trigDist->Projection(0);
    float totalTrig = (float)trigDist1D->Integral(trigDist1D->GetXaxis()->FindBin(trigPTLow + epsilon), trigDist1D->GetXaxis()->FindBin(trigPTHigh - epsilon));
    
    //float totalTrig = (float)trigHHDist->Integral(trigHHDist->GetXaxis()->FindBin(trigPTLow + epsilon), trigHHDist->GetXaxis()->FindBin(trigPTHigh - epsilon), 1, trigHHDist->GetYaxis()->GetNbins());
      
    TH1D* vtxZmixbins = (TH1D*)list->FindObject("fVtxZmixbins");
    Int_t numbinsZvtx = vtxZmixbins->GetXaxis()->GetNbins();

    //THnSparseF *dphiHH = (THnSparseF*)list->FindObject("fDphiHH");
    //THnSparseF *dphiHHMixed = (THnSparseF*)list->FindObject("fDPhiHHMixed");

    THnSparseF* dphiHH[numbinsZvtx];
    THnSparseF* dphiHHMixed[numbinsZvtx];

    TH2D* hh[numbinsZvtx];
    TH2D* hhMixed[numbinsZvtx];
    TH2D* hh2D[numbinsZvtx];
    
    TH2D* hhTotal;
    TH2D* hhMixedTotal;
    TH2D* uncorrhh2D;    
   
    
    for(int izvtx = 0; izvtx < numbinsZvtx; izvtx++){

        dphiHH[izvtx] = (THnSparseF*)list->FindObject(Form("fDphiHHz%i", izvtx));
        dphiHH[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon); 
        dphiHH[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon, assocPTHigh - epsilon); 
   
        dphiHHMixed[izvtx] = (THnSparseF*)list->FindObject(Form("fDphiHHMixedz%i", izvtx));
        dphiHHMixed[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon); 
        dphiHHMixed[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon,assocPTHigh - epsilon); 

        
        hh[izvtx] = dphiHH[izvtx]->Projection(2,3);
        hh[izvtx]->Sumw2();
        hhMixed[izvtx] = dphiHHMixed[izvtx]->Projection(2,3);
        hhMixed[izvtx]->Sumw2();
    
        //hh[izvtx] = projectWithEfficiencyCorrections(dphiHH[izvtx], hadronEff, 2.0, 4.0);
        //hhMixed[izvtx] = projectWithEfficiencyCorrections(dphiHHMixed[izvtx], hadronEff, 2.0, 4.0);

        hh2D[izvtx] = makehhCorrections(hh[izvtx], hhMixed[izvtx]);
        hh2D[izvtx]->SetName(Form("hh2Dz%i", izvtx));

        if(izvtx == 0){
            hhTotal = (TH2D*)hh2D[izvtx]->Clone("hh2D");
            uncorrhh2D = (TH2D*)hh[izvtx]->Clone("uncorrhh2D");
        }else{
            hhTotal->Add(hh2D[izvtx]);
            uncorrhh2D->Add(hh[izvtx]);
        }
    }

    hhTotal->Scale(1.0/totalTrig);

    TH1D* hhdphi = hhTotal->ProjectionY("hhdphi", hhTotal->GetXaxis()->FindBin(-1.2 + epsilon), hhTotal->GetXaxis()->FindBin(1.2 - epsilon));
    hhdphi->Scale(1.0/(hhdphi->Integral()));

    uncorrhh2D->Scale(1.0/(uncorrhh2D->Integral(uncorrhh2D->GetXaxis()->FindBin(-1.2 + epsilon), uncorrhh2D->GetXaxis()->FindBin(1.2 - epsilon), 1, uncorrhh2D->GetYaxis()->GetNbins())));
    
    TFile* output = new TFile(Form("trig_%.1f_%.1f_assoc_%.1f_%.1f_effcorr_hh%s.root", trigPTLow, trigPTHigh, assocPTLow, assocPTHigh, mult.c_str()), "RECREATE");
    hhTotal->Write();
    hhdphi->Write();
    uncorrhh2D->Write(); 
    trigHHDist->Write();
}

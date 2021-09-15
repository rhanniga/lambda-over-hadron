TH2D* makeCorrections(TH3D* same, TH3D* mixed, Float_t lowmass, Float_t highmass, TH1D* sameEta, TH1D* mixedEta, Int_t zbin){
    same->GetZaxis()->SetRangeUser(lowmass, highmass);
    mixed->GetZaxis()->SetRangeUser(lowmass, highmass);
    TH2D* same2D = (TH2D*)same->Project3D("xye");
    same2D->Sumw2();
    TH2D* mix2D = (TH2D*)mixed->Project3D("xye");
    mix2D->Sumw2();

    //same2D->RebinX(4);
    //mix2D->RebinX(4);
    //same2D->RebinY(4);
    //mix2D->RebinY(4);

    //same3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    //mix3D->GetYaxis()->SetRangeUser(-1.2, 1.2);
    same2D->GetYaxis()->SetRange(1,same2D->GetYaxis()->GetNbins());
    mix2D->GetYaxis()->SetRange(1, mix2D->GetYaxis()->GetNbins());

    Float_t scale = 0.0;
    
    //get d-eta 1D plots for same and mixed distributions for each zbin
    sameEta = same2D->ProjectionX(Form("sameEta_zvtx_%i", zbin), 1, same2D->GetYaxis()->GetNbins());
    mixedEta = mix2D->ProjectionX(Form("mixedEta_zvtx_%i", zbin), 1, mix2D->GetYaxis()->GetNbins());


    scale = 0.5*(float)(mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(0.01), mix2D->GetYaxis()->FindBin(0.0)) + mix2D->GetBinContent(mix2D->GetXaxis()->FindBin(-0.01), mix2D->GetYaxis()->FindBin(0.0)));
    printf("scale: %e \n", scale);
    same2D->Divide(mix2D);
    same2D->Scale(scale);

    same->GetZaxis()->SetRange(0,0);
    mixed->GetZaxis()->SetRange(0,0);
    return same2D;
}
//-------------------------------------------------------------------------------------------
TH3D* projectWithEfficiencyCorrections(THnSparseF* sparse, TH1D* eff, float assocPTLow, float assocPTHigh){
    Int_t lowbin = eff->GetXaxis()->FindBin(assocPTLow + 0.00001);
    Int_t highbin = eff->GetXaxis()->FindBin(assocPTHigh - 0.00001);
    Int_t axes[] = {2,3,4};
    TH3D* corr;
    TH3D* buff;
    //printf("lowbin: %i, highbin:%i \n\n", lowbin, highbin);
    //printf("low range: %i, high range: %i\n\n", sparse->GetAxis(1)->FindBin(eff->GetBinCenter(lowbin)), sparse->GetAxis(1)->FindBin(eff->GetBinCenter(highbin)));
    for(Int_t i = lowbin; i <= highbin; i++){
        sparse->GetAxis(1)->SetRange(sparse->GetAxis(1)->FindBin(eff->GetBinCenter(i)), sparse->GetAxis(1)->FindBin(eff->GetBinCenter(i)));
        buff = (TH3D*)sparse->Projection(2,3,4);
        buff->Sumw2();
        buff->Scale(1.0/eff->GetBinContent(i));
        //buff->Scale(1.0);
        if(i==lowbin){
            corr = (TH3D*)buff->Clone(Form("effcorr%s", sparse->GetName()));
        }else{
            corr->Add(buff);
        }
    }
    return corr;
}
//--------------------------------------------------------------------------------------------
void makeMixCorrectionsZVertexEff(string inputName, string outputStr, int multLow, int multHigh, float trigPTLow, float trigPTHigh, float assocPTLow, float assocPTHigh, string listSuffix = "_", float pLow = 1.014, float pHigh = 1.026, float lsbLow = 0.995, float lsbHigh = 1.005, float rsbLow = 1.040, float rsbHigh = 1.060){
    //TFile *effFile = new TFile("~/repos/utaustin/efficiency/17f2bTPC80efficiency.root");
    //TH1D* phiEff = (TH1D*)effFile->Get("phiPTEff");
    //TH1D* hadronEff = (TH1D*)effFile->Get("hadronPTEff");
    
    TFile *histoFile = new TFile(inputName.c_str());
    //string mult = inputName.substr(inputName.find("_", inputName.find("_")+1), inputName.find(".") - inputName.find("_", inputName.find("_")+1));
    string mult = "_" + std::to_string(multLow) + "_" + std::to_string(multHigh);
    //TList* list = (TList*) histoFile->Get(Form("truePhiCorr_mult%s", mult.c_str()));
    TList* list = (TList*) histoFile->Get("hhCorr_mult_0_20_");
    
    printf("phiCorr_mult%s%s\n", mult.c_str(), listSuffix.c_str());

    /*
    TH2D *trigSameUSDist = (TH2D*)list->FindObject("fTrigSameUSDist");
    TH2D *trigSameLSDist = (TH2D*)list->FindObject("fTrigSameLSDist");

    float totalTrigSameUS = (float)trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(trigPTLow), trigSameUSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameUSDist->GetYaxis()->GetNbins());
    float totalTrigSameLS = (float)trigSameLSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(trigPTLow), trigSameLSDist->GetXaxis()->FindBin(trigPTHigh), 1, trigSameLSDist->GetYaxis()->GetNbins());
   */
    Float_t epsilon = 0.0001;

    THnSparseF* trigDist = (THnSparseF*)list->FindObject("fTrigDist");
    TH1D* trigDist1D = (TH1D*)trigDist->Projection(0);
    float totalTrigSameUS = (float)trigDist1D->Integral(trigDist1D->GetXaxis()->FindBin(trigPTLow + epsilon), trigDist1D->GetXaxis()->FindBin(trigPTHigh - epsilon));

    TH1D* vtxZmixbins = (TH1D*)list->FindObject("fVtxZmixbins");
    Int_t numbinsZvtx = vtxZmixbins->GetXaxis()->GetNbins();
    //Int_t numbinsZvtx = 10;

    printf("first two histos done \n");

    TH1D* sameUSPeakEta[numbinsZvtx];
    TH1D* mixedUSPeakEta[numbinsZvtx];
    TH1D* sameLSPeakEta[numbinsZvtx];
    TH1D* mixedLSPeakEta[numbinsZvtx];
    TH1D* sameUSRsideEta[numbinsZvtx];
    TH1D* mixedUSRsideEta[numbinsZvtx];
    TH1D* sameLSRsideEta[numbinsZvtx];
    TH1D* mixedLSRsideEta[numbinsZvtx];
    TH1D* sameUSLsideEta[numbinsZvtx];
    TH1D* mixedUSLsideEta[numbinsZvtx];
    TH1D* sameLSLsideEta[numbinsZvtx];
    TH1D* mixedLSLsideEta[numbinsZvtx];

    THnSparseF *dphiHPhi[numbinsZvtx];
    THnSparseF *dphiHKK[numbinsZvtx];
    THnSparseF *dphiHPhiMixed[numbinsZvtx];
    THnSparseF *dphiHKKMixed[numbinsZvtx];

    TH3D *hPhi[numbinsZvtx];
    TH3D *hKK[numbinsZvtx];
    TH3D *hPhiMixed[numbinsZvtx];
    TH3D *hKKMixed[numbinsZvtx];

    TH2D* hPhi2Dpeak[numbinsZvtx];
    TH2D* hKK2Dpeak[numbinsZvtx];
    TH2D* hPhi2DRside[numbinsZvtx];
    TH2D* hKK2DRside[numbinsZvtx];
    TH2D* hPhi2DLside[numbinsZvtx];
    TH2D* hKK2DLside[numbinsZvtx];

    TH1D* USinvmass[numbinsZvtx];
    TH1D* LSinvmass[numbinsZvtx];

    TH2D* hPhi2DpeakTotal;
    TH2D* hKK2DpeakTotal;
    TH2D* hPhi2DRsideTotal;
    TH2D* hKK2DRsideTotal;
    TH2D* hPhi2DLsideTotal;
    TH2D* hKK2DLsideTotal;

    TH3D* hPhiTotal;
    TH3D* hPhiMixedTotal;
    TH3D* hKKTotal;
    TH3D* hKKMixedTotal;

    TH1D* USinvmassTotal;
    TH1D* LSinvmassTotal;

    Float_t peakLow = pLow + epsilon;
    Float_t peakHigh = pHigh - epsilon;
    Float_t leftSBLow = lsbLow + epsilon;
    Float_t leftSBHigh = lsbHigh - epsilon;
    Float_t rightSBLow = rsbLow + epsilon;
    Float_t rightSBHigh = rsbHigh - epsilon;

    for(int izvtx = 0; izvtx < numbinsZvtx; izvtx++){
        dphiHPhi[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiz%i", izvtx));
        dphiHKK[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKz%i", izvtx));
        dphiHPhiMixed[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHPhiMixedz%i", izvtx));
        dphiHKKMixed[izvtx] = (THnSparseF *)list->FindObject(Form("fDphiHKKMixedz%i", izvtx));

        //make 4D THnProjections projection to do mixed event corrections    
        dphiHPhi[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon); 
        dphiHPhi[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon,assocPTHigh - epsilon); 
        dphiHKK[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon);
        dphiHKK[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon,assocPTHigh - epsilon);

        dphiHPhiMixed[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon); 
        dphiHPhiMixed[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon,assocPTHigh - epsilon); 
        dphiHKKMixed[izvtx]->GetAxis(0)->SetRangeUser(trigPTLow + epsilon, trigPTHigh - epsilon);
        dphiHKKMixed[izvtx]->GetAxis(1)->SetRangeUser(assocPTLow + epsilon,assocPTHigh - epsilon);

        dphiHPhi[izvtx]->GetAxis(4)->SetRange(1,dphiHPhi[izvtx]->GetAxis(4)->GetNbins()); 
        dphiHKK[izvtx]->GetAxis(4)->SetRange(1,dphiHKK[izvtx]->GetAxis(4)->GetNbins());
        dphiHPhiMixed[izvtx]->GetAxis(4)->SetRange(1,dphiHPhiMixed[izvtx]->GetAxis(4)->GetNbins());
        dphiHKKMixed[izvtx]->GetAxis(4)->SetRange(1,dphiHKKMixed[izvtx]->GetAxis(4)->GetNbins());

        Int_t axes[] = {2,3,4};
        hPhi[izvtx] = (TH3D*)dphiHPhi[izvtx]->Projection(2,3,4);
        hKK[izvtx] = (TH3D*)dphiHKK[izvtx]->Projection(2,3,4);
        hPhiMixed[izvtx] = (TH3D*)dphiHPhiMixed[izvtx]->Projection(2,3,4);
        hKKMixed[izvtx] = (TH3D*)dphiHKKMixed[izvtx]->Projection(2,3,4);
        //hPhi[izvtx] = projectWithEfficiencyCorrections(dphiHPhi[izvtx], phiEff, 2.0, 4.0);
        //hKK[izvtx] = projectWithEfficiencyCorrections(dphiHKK[izvtx], phiEff, 2.0, 4.0);
        //hPhiMixed[izvtx] = projectWithEfficiencyCorrections(dphiHPhiMixed[izvtx], phiEff, 2.0, 4.0);
        //hKKMixed[izvtx] = projectWithEfficiencyCorrections(dphiHKKMixed[izvtx], phiEff, 2.0, 4.0);

        hPhi[izvtx]->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
        hPhi[izvtx]->GetYaxis()->SetRange(1, hPhi[izvtx]->GetYaxis()->GetNbins());
        USinvmass[izvtx] = (TH1D*)hPhi[izvtx]->Project3D("ze");
        hPhi[izvtx]->GetXaxis()->SetRange(1,hPhi[izvtx]->GetXaxis()->GetNbins());
        hKK[izvtx]->GetXaxis()->SetRangeUser(-1.2 + epsilon, 1.2 - epsilon);
        hKK[izvtx]->GetYaxis()->SetRange(1, hKK[izvtx]->GetYaxis()->GetNbins());
        LSinvmass[izvtx] = (TH1D*)hKK[izvtx]->Project3D("ze");
        hKK[izvtx]->GetXaxis()->SetRange(1, hKK[izvtx]->GetXaxis()->GetNbins());
        
        hPhi2Dpeak[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], peakLow, peakHigh, sameUSPeakEta[izvtx], mixedUSPeakEta[izvtx], izvtx);
        hPhi2Dpeak[izvtx]->SetName(Form("hPhi2Dpeakz%i", izvtx));
        //sameUSPeakEta[izvtx]->SetName(Form("sameUSPeakEta_zvtx_%i", izvtx));
        //mixedUSPeakEta[izvtx]->SetName(Form("mixedUSPeakEta_zvtx_%i", izvtx));

        hKK2Dpeak[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], peakLow, peakHigh, sameLSPeakEta[izvtx], mixedLSPeakEta[izvtx], izvtx);
        hKK2Dpeak[izvtx]->SetName(Form("hKK2Dpeakz%i", izvtx));
        //sameLSPeakEta[izvtx]->SetName(Form("sameLSPeakEta_zvtx_%i", izvtx));
        //mixedLSPeakEta[izvtx]->SetName(Form("mixedLSPeakEta_zvtx_%i", izvtx));

        hPhi2DRside[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], rightSBLow, rightSBHigh, sameUSRsideEta[izvtx], mixedUSRsideEta[izvtx], izvtx);
        hPhi2DRside[izvtx]->SetName(Form("hPhi2DRsidez%i", izvtx));
        //sameUSRsideEta[izvtx]->SetName(Form("sameUSRsideEta_zvtx_%i", izvtx));
        //mixedUSRsideEta[izvtx]->SetName(Form("mixedUSRsideEta_zvtx_%i", izvtx));

        hKK2DRside[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], rightSBLow, rightSBHigh, sameLSRsideEta[izvtx], mixedLSRsideEta[izvtx], izvtx);
        hKK2DRside[izvtx]->SetName(Form("hKK2DRsidez%i", izvtx));
        //sameLSRsideEta[izvtx]->SetName(Form("sameLSRsideEta_zvtx_%i", izvtx));
        //mixedLSRsideEta[izvtx]->SetName(Form("mixedLSRsideEta_zvtx_%i", izvtx));

        hPhi2DLside[izvtx] = makeCorrections(hPhi[izvtx], hPhiMixed[izvtx], leftSBLow, leftSBHigh, sameUSLsideEta[izvtx], mixedUSLsideEta[izvtx], izvtx);
        hPhi2DLside[izvtx]->SetName(Form("hPhi2DLsidez%i", izvtx));
        //sameUSLsideEta[izvtx]->SetName(Form("sameUSLsideEta_zvtx_%i", izvtx));
        //mixedUSLsideEta[izvtx]->SetName(Form("mixedUSLsideEta_zvtx_%i", izvtx));

        hKK2DLside[izvtx] = makeCorrections(hKK[izvtx], hKKMixed[izvtx], leftSBLow, leftSBHigh, sameLSLsideEta[izvtx], mixedLSLsideEta[izvtx], izvtx);
        hKK2DLside[izvtx]->SetName(Form("hKK2DLsidez%i", izvtx));
        //sameLSLsideEta[izvtx]->SetName(Form("sameLSLsideEta_zvtx_%i", izvtx));
        //mixedLSLsideEta[izvtx]->SetName(Form("mixedLSLsideEta_zvtx_%i", izvtx));
        if(izvtx==0){
            hPhi2DpeakTotal = (TH2D*)hPhi2Dpeak[izvtx]->Clone("hPhi2Dpeak");
            hPhi2DLsideTotal = (TH2D*)hPhi2DLside[izvtx]->Clone("hPhi2DLside");
            hPhi2DRsideTotal = (TH2D*)hPhi2DRside[izvtx]->Clone("hPhi2DRside");
            hKK2DpeakTotal = (TH2D*)hKK2Dpeak[izvtx]->Clone("hKK2Dpeak");
            hKK2DLsideTotal = (TH2D*)hKK2DLside[izvtx]->Clone("hKK2DLside");
            hKK2DRsideTotal = (TH2D*)hKK2DRside[izvtx]->Clone("hKK2DRside");

            hPhiTotal = (TH3D*)hPhi[izvtx]->Clone("hPhiTotal");
            hPhiMixedTotal = (TH3D*)hPhiMixed[izvtx]->Clone("hPhiMixedTotal");
            hKKTotal = (TH3D*)hKK[izvtx]->Clone("hKKTotal");
            hKKMixedTotal = (TH3D*)hKKMixed[izvtx]->Clone("hKKMixedTotal");

            USinvmassTotal = (TH1D*)USinvmass[izvtx]->Clone("USinvmassTotal");
            LSinvmassTotal = (TH1D*)LSinvmass[izvtx]->Clone("LSinvmassTotal");
       }else{
            hPhi2DpeakTotal->Add(hPhi2Dpeak[izvtx]);
            hPhi2DRsideTotal->Add(hPhi2DRside[izvtx]);
            hPhi2DLsideTotal->Add(hPhi2DLside[izvtx]);
            hKK2DpeakTotal->Add(hKK2Dpeak[izvtx]);
            hKK2DRsideTotal->Add(hKK2DRside[izvtx]);
            hKK2DLsideTotal->Add(hKK2DLside[izvtx]);

            hPhiTotal->Add(hPhi[izvtx]);
            hPhiMixedTotal->Add(hPhiMixed[izvtx]);
            hKKTotal->Add(hKK[izvtx]);
            hKKMixedTotal->Add(hKKMixed[izvtx]);

            USinvmassTotal->Add(USinvmass[izvtx]);
            LSinvmassTotal->Add(LSinvmass[izvtx]);
       }
    }

    hPhi2DpeakTotal->Scale(1.0/totalTrigSameUS);
    hPhi2DRsideTotal->Scale(1.0/totalTrigSameUS);
    hPhi2DLsideTotal->Scale(1.0/totalTrigSameUS);
    hKK2DpeakTotal->Scale(1.0/totalTrigSameUS);
    hKK2DRsideTotal->Scale(1.0/totalTrigSameUS);
    hKK2DLsideTotal->Scale(1.0/totalTrigSameUS);

    //Create some uncorrected same/mixed event 2D histos
    hPhiTotal->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    hPhiMixedTotal->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    hKKTotal->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    hKKMixedTotal->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    TH2D* uncorrhPhi2Dpeak = (TH2D*)hPhiTotal->Project3D("xye");
    uncorrhPhi2Dpeak->Sumw2();
    uncorrhPhi2Dpeak->SetName("uncorrhPhi2Dpeak");
    TH2D* uncorrhKK2Dpeak = (TH2D*)hKKTotal->Project3D("xye");
    uncorrhKK2Dpeak->Sumw2();
    uncorrhKK2Dpeak->SetName("uncorrhKK2Dpeak");
    TH2D* uncorrhPhiMixed2Dpeak = (TH2D*)hPhiMixedTotal->Project3D("xye");
    uncorrhPhiMixed2Dpeak->Sumw2();
    uncorrhPhiMixed2Dpeak->SetName("uncorrhPhiMixed2Dpeak");
    TH2D* uncorrhKKMixed2Dpeak = (TH2D*)hKKMixedTotal->Project3D("xye");
    uncorrhKKMixed2Dpeak->Sumw2();
    uncorrhKKMixed2Dpeak->SetName("uncorrhKKMixed2Dpeak");

    hPhiTotal->GetZaxis()->SetRangeUser(rightSBLow, rightSBHigh);
    hPhiMixedTotal->GetZaxis()->SetRangeUser(rightSBLow, rightSBHigh);
    hKKTotal->GetZaxis()->SetRangeUser(rightSBLow, rightSBHigh);
    hKKMixedTotal->GetZaxis()->SetRangeUser(rightSBLow, rightSBHigh);
    TH2D* uncorrhPhi2DRside = (TH2D*)hPhiTotal->Project3D("xye");
    uncorrhPhi2DRside->Sumw2();
    uncorrhPhi2DRside->SetName("uncorrhPhi2DRside");
    TH2D* uncorrhKK2DRside = (TH2D*)hKKTotal->Project3D("xye");
    uncorrhKK2DRside->Sumw2();
    uncorrhKK2DRside->SetName("uncorrhKK2DRside");
    TH2D* uncorrhPhiMixed2DRside = (TH2D*)hPhiMixedTotal->Project3D("xye");
    uncorrhPhiMixed2DRside->Sumw2();
    uncorrhPhiMixed2DRside->SetName("uncorrhPhiMixed2DRside");
    TH2D* uncorrhKKMixed2DRside = (TH2D*)hKKMixedTotal->Project3D("xye");
    uncorrhKKMixed2DRside->Sumw2();
    uncorrhKKMixed2DRside->SetName("uncorrhKKMixed2DRside");

    hPhiTotal->GetZaxis()->SetRangeUser(leftSBLow, leftSBHigh);
    hPhiMixedTotal->GetZaxis()->SetRangeUser(leftSBLow, leftSBHigh);
    hKKTotal->GetZaxis()->SetRangeUser(leftSBLow, leftSBHigh);
    hKKMixedTotal->GetZaxis()->SetRangeUser(leftSBLow, leftSBHigh);
    TH2D* uncorrhPhi2DLside = (TH2D*)hPhiTotal->Project3D("xye");
    uncorrhPhi2DLside->Sumw2();
    uncorrhPhi2DLside->SetName("uncorrhPhi2DLside");
    TH2D* uncorrhKK2DLside = (TH2D*)hKKTotal->Project3D("xye");
    uncorrhKK2DLside->Sumw2();
    uncorrhKK2DLside->SetName("uncorrhKK2DLside");
    TH2D* uncorrhPhiMixed2DLside = (TH2D*)hPhiMixedTotal->Project3D("xye");
    uncorrhPhiMixed2DLside->Sumw2();
    uncorrhPhiMixed2DLside->SetName("uncorrhPhiMixed2DLside");
    TH2D* uncorrhKKMixed2DLside = (TH2D*)hKKMixedTotal->Project3D("xye");
    uncorrhKKMixed2DLside->Sumw2();
    uncorrhKKMixed2DLside->SetName("uncorrhKKMixed2DLside");

    //Make some ratio plots of mixed event US over LS for sidbeand and peak regions
    //also make projections onto delta eta and delta phi (on restricted delta eta range)
    TH2D* mixedratioRSB = (TH2D*)uncorrhPhiMixed2DRside->Clone("mixedratioRSB");
    mixedratioRSB->Divide(uncorrhKKMixed2DRside);
    TH1D* mixedratioRSBdeta = mixedratioRSB->ProjectionX("mixedratioRSBdeta", 1, mixedratioRSB->GetYaxis()->GetNbins());
    TH1D* mixedratioRSBdphi = mixedratioRSB->ProjectionY("mixedratioRSBdphi", mixedratioRSB->GetXaxis()->FindBin(-1.2+ epsilon), mixedratioRSB->GetXaxis()->FindBin(1.2 - epsilon));
    mixedratioRSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioLSB = (TH2D*)uncorrhPhiMixed2DLside->Clone("mixedratioLSB");
    mixedratioLSB->Divide(uncorrhKKMixed2DLside);
    TH1D* mixedratioLSBdeta = mixedratioLSB->ProjectionX("mixedratioLSBdeta", 1, mixedratioLSB->GetYaxis()->GetNbins());
    TH1D* mixedratioLSBdphi = mixedratioLSB->ProjectionY("mixedratioLSBdphi", mixedratioLSB->GetXaxis()->FindBin(-1.2 + epsilon), mixedratioLSB->GetXaxis()->FindBin(1.2 - epsilon));
    mixedratioLSB->GetXaxis()->SetRange(0,0);

    TH2D* mixedratioPeak = (TH2D*)uncorrhPhiMixed2Dpeak->Clone("mixedratioPeak");
    mixedratioPeak->Divide(uncorrhKKMixed2Dpeak);
    TH1D* mixedratioPeakdeta = mixedratioPeak->ProjectionX("mixedratioPeakdeta", 1, mixedratioPeak->GetYaxis()->GetNbins());
    TH1D* mixedratioPeakdphi = mixedratioPeak->ProjectionY("mixedratioPeakdphi", mixedratioPeak->GetXaxis()->FindBin(-1.2 + epsilon), mixedratioPeak->GetXaxis()->FindBin(1.2 - epsilon));
    mixedratioPeak->GetXaxis()->SetRange(0,0);

    //make ratio plot of just 1 zvtx bin as a check:
    hPhiMixed[6]->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    TH2D* mixedratioPeakZ2 = (TH2D*)hPhiMixed[6]->Project3D("xye");
    hKKMixed[6]->GetZaxis()->SetRangeUser(peakLow, peakHigh);
    TH2D* hist = (TH2D*)hKKMixed[6]->Project3D("xye");
    mixedratioPeakZ2->Divide(hist);
    TH1D* mixedratioPeakZ2deta = mixedratioPeakZ2->ProjectionX("mixedratioPeakZ2deta");


    
    TFile* output = new TFile(Form("trig_%.1f_%.1f_assoc_%.1f_%.1f_effcorr_hPhi%s_%s.root", trigPTLow, trigPTHigh, assocPTLow, assocPTHigh, mult.c_str(), outputStr.c_str()), "RECREATE");
    hPhi2DpeakTotal->Write();
    hKK2DpeakTotal->Write();
    hPhi2DRsideTotal->Write();
    hKK2DRsideTotal->Write();
    hPhi2DLsideTotal->Write();
    hKK2DLsideTotal->Write();
    uncorrhPhi2Dpeak->Write();
    uncorrhKK2Dpeak->Write();
    uncorrhPhiMixed2Dpeak->Write();
    uncorrhKKMixed2Dpeak->Write();
    uncorrhPhi2DRside->Write();
    uncorrhKK2DRside->Write();
    uncorrhPhiMixed2DRside->Write();
    uncorrhKKMixed2DRside->Write();
    uncorrhPhi2DLside->Write();
    uncorrhKK2DLside->Write();
    uncorrhPhiMixed2DLside->Write();
    uncorrhKKMixed2DLside->Write();

    mixedratioRSB->Write();
    mixedratioRSBdeta->Write();
    mixedratioRSBdphi->Write();
    mixedratioLSB->Write();
    mixedratioLSBdeta->Write();
    mixedratioLSBdphi->Write();
    mixedratioPeak->Write();
    mixedratioPeakdeta->Write();
    mixedratioPeakdphi->Write();

    mixedratioPeakZ2->Write();
    mixedratioPeakZ2deta->Write();

    USinvmassTotal->Write();
    LSinvmassTotal->Write();
/*
    for(int i=0; i< numbinsZvtx; i++){
        sameUSPeakEta[i]->Write();
        mixedUSPeakEta[i]->Write();
        sameLSPeakEta[i]->Write();
        mixedLSPeakEta[i]->Write();
        sameUSRsideEta[i]->Write();
        mixedUSRsideEta[i]->Write();
        sameLSRsideEta[i]->Write();
        mixedLSRsideEta[i]->Write();
        sameUSLsideEta[i]->Write();
        mixedUSLsideEta[i]->Write();
        sameLSLsideEta[i]->Write();
        mixedLSLsideEta[i]->Write();
    }
*/
    //trigSameUSDist->Write();
    //trigSameLSDist->Write();

}

TH2D* makeCorrections(THnSparse* same, THnSparse* mixed, Float_t lowmass, Float_t highmass, TH1D** sameEta, TH1D** mixedEta, float* trigMixScale, float totalTrigSame){


    same->GetAxis(2)->SetRangeUser(lowmass, highmass);
    mixed->GetAxis(2)->SetRangeUser(lowmass, highmass);
    TH3D* same3D = same->Projection(0, 1, 3);
    same3D->Sumw2();
    TH3D* mix3D = mixed->Projection(0, 1, 3);
    mix3D->Sumw2();




    same3D->GetYaxis()->SetRange(1,same3D->GetYaxis()->GetNbins());
    mix3D->GetYaxis()->SetRange(1, mix3D->GetYaxis()->GetNbins());

    Float_t scale = 0.0;
    TH2D* same2D[10];
    TH2D* mix2D[10];
    TH2D* same2DTotal;

    for(int zbin = 0; zbin < 10; zbin++){
        same3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        same2D[zbin] = (TH2D*)same3D->Project3D("xye");
        same2D[zbin]->SetName(Form("2dproj_zbin%i", zbin));

        mix3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        mix2D[zbin] = (TH2D*)mix3D->Project3D("xye");
        mix2D[zbin]->SetName(Form("mix2droj_zbin%i", zbin));

        //get d-eta 1D plots for same and mixed distributions for each zbin
        sameEta[zbin] = same2D[zbin]->ProjectionX(Form("sameEta_zvtx_%i", zbin), 1, same2D[zbin]->GetYaxis()->GetNbins());
        mixedEta[zbin] = mix2D[zbin]->ProjectionX(Form("mixedEta_zvtx_%i", zbin), 1, mix2D[zbin]->GetYaxis()->GetNbins());

        //scale mixed by number of triggers in zvtx bin
        //mix2D[zbin]->Scale(1.0/trigMixScale[zbin]);

        scale = 0.5*(float)(mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)) + mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(-0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)));
        printf("scale: %e \n", scale);
        same2D[zbin]->Divide(mix2D[zbin]);
        same2D[zbin]->Scale(scale);
        if(zbin==0){
            same2DTotal = (TH2D*)same2D[zbin]->Clone("2dproj_total");
        }else{
            same2DTotal->Add(same2D[zbin]);
        }

        same3D->GetZaxis()->SetRange(0,0);
        mix3D->GetZaxis()->SetRange(0,0);
    }

    same->GetAxis(3)->SetRange(0,0);
    mixed->GetAxis(3)->SetRange(0,0);
    float totalSame = same2DTotal->Integral();
    //same2DTotal->Scale(1.0/totalTrigSame);
    printf("Total Trig: %f, Total pairs: %f\n, ratio: %f", totalTrigSame, totalSame, totalSame/totalTrigSame);
    return same2DTotal;
}

//---------------------------------------------------------------------------------------------
TH2D* makehhCorrections(TH3D* same3D, TH3D* mix3D){

    same3D->GetYaxis()->SetRange(1,same3D->GetYaxis()->GetNbins());
    mix3D->GetYaxis()->SetRange(1, mix3D->GetYaxis()->GetNbins());

    Float_t scale = 0.0;
    TH2D* same2D[10];
    TH2D* mix2D[10];
    TH2D* same2DTotal;

    for(int zbin = 0; zbin < 10; zbin++){
        same3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        same2D[zbin] = (TH2D*)same3D->Project3D("xye");
        same2D[zbin]->SetName(Form("2dproj_zbin%i", zbin));

        mix3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        mix2D[zbin] = (TH2D*)mix3D->Project3D("xye");
        mix2D[zbin]->SetName(Form("mix2droj_zbin%i", zbin));

        mix3D->GetZaxis()->SetRange(zbin+1, zbin+1);
        mix2D[zbin] = (TH2D*)mix3D->Project3D("xye");
        mix2D[zbin]->SetName(Form("mix2droj_zbin%i", zbin));

        scale = 0.5*(float)(mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)) + mix2D[zbin]->GetBinContent(mix2D[zbin]->GetXaxis()->FindBin(-0.01), mix2D[zbin]->GetYaxis()->FindBin(0.0)));
        printf("scale: %e \n", scale);
        same2D[zbin]->Divide(mix2D[zbin]);
        same2D[zbin]->Scale(scale);
        if(zbin==0){
            same2DTotal = (TH2D*)same2D[zbin]->Clone("2dproj_total");
        }else{
            same2DTotal->Add(same2D[zbin]);
        }

        same3D->GetZaxis()->SetRange(0,0);
        mix3D->GetZaxis()->SetRange(0,0);
    }

    return same2DTotal;
}

//--------------------------------------------------------------------------------------------
void makeMixCorrections(float trigPTLow, float TRIG_PT_HIGH, float ASSOC_PT_LOW, float ASSOC_PT_HIGH){
    TFile *histoFile = new TFile("../online/output/cent_0_20.root");
    TList* list = (TList*) histoFile->Get("h-lambda");
    int centLow = 0;
    int centHigh = 20;

    float EPSILON = 0.00001;
    
    float TRIG_PT_LOW = 4;
    float TRIG_PT_HIGH = 8 - EPSILON;
    float ASSOC_PT_LOW = 2;
    float ASSOC_PT_HIGH = 4 - EPSILON;

    float LSB_MIN = 1.085;
    float LSB_MAX = 1.10 - EPSILON;
    float SIG_MIN = 1.108;
    float SIG_MAX = 1.124 - EPSILON;
    float RSB_MIN = 1.14;
    float RSB_MAX = 1.155 - EPSILON;


    THnSparseF*triggerDist = (THnSparseF*)list->FindObject("fTriggerDist");
    TH2D *trigSameUSDist = (TH2D*)triggerDist->Projection(3, 0);
    TH2D *trigSameLSDist = (TH2D*)triggerDist->Projection(3, 0);

    float trigMixScalesUS[10] = {};
    float trigMixScalesLS[10] = {};

    for(int i = 0; i < 10; i++){
        trigMixScalesUS[i] = (float) trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_LOW), trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_HIGH), i+1, i+1);
        trigMixScalesLS[i] = (float) trigSameLSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_LOW), trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_HIGH), i+1, i+1);
    }

    float totalTrigSameUS = (float)trigSameUSDist->Integral(trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_LOW), trigSameUSDist->GetXaxis()->FindBin(TRIG_PT_HIGH), 1, trigSameUSDist->GetYaxis()->GetNbins());
    float totalTrigSameLS = (float)trigSameLSDist->Integral(trigSameLSDist->GetXaxis()->FindBin(TRIG_PT_LOW), trigSameLSDist->GetXaxis()->FindBin(TRIG_PT_HIGH), 1, trigSameLSDist->GetYaxis()->GetNbins());

    THnSparseF *dphiHLambda = (THnSparseF *)list->FindObject("fDphiHLambdaEff");
    THnSparseF *dphiHLambdaMixed = (THnSparseF *)list->FindObject("fDphiHLambdaMixed");
    THnSparseF *dphiHLambdaLS = (THnSparseF *)list->FindObject("fDphiHLambdaLS");
    THnSparseF *dphiHLambdaLSMixed = (THnSparseF *)list->FindObject("fDphiHLambdaLSMixed");

    THnSparseF *dphiHH = (THnSparseF*)list->FindObject("fDphiHHEff");
    THnSparseF *dphiHHMixed = (THnSparseF*)list->FindObject("fDphiHHMixed");

    //make 4D THnProjections projection to do mixed event corrections
    dphiHLambda->GetAxis(0)->SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH);
    dphiHLambda->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);
    dphiHLambdaLS->GetAxis(0)->SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH);
    dphiHLambdaLS->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);
    dphiHH->GetAxis(0)->SetRangeUser(TRIG_PT_LOW,TRIG_PT_HIGH);
    dphiHH->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);

    dphiHLambdaMixed->GetAxis(0)->SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH);
    dphiHLambdaMixed->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);
    dphiHLambdaLSMixed->GetAxis(0)->SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH);
    dphiHLambdaLSMixed->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);
    dphiHHMixed->GetAxis(0)->SetRangeUser(TRIG_PT_LOW,TRIG_PT_HIGH);
    dphiHHMixed->GetAxis(1)->SetRangeUser(ASSOC_PT_LOW,ASSOC_PT_HIGH);


    dphiHLambda->GetAxis(4)->SetRange(1,dphiHLambda->GetAxis(4)->GetNbins());
    dphiHLambda->GetAxis(5)->SetRangeUser(-10.0, 10.0);

    dphiHLambdaLS->GetAxis(4)->SetRange(1,dphiHLambdaLS->GetAxis(4)->GetNbins());
    dphiHLambdaLS->GetAxis(5)->SetRangeUser(-10.0, 10.0);

    dphiHLambdaMixed->GetAxis(4)->SetRange(1,dphiHLambdaMixed->GetAxis(4)->GetNbins());
    dphiHLambdaMixed->GetAxis(5)->SetRangeUser(-10.0, 10.0);

    dphiHLambdaLSMixed->GetAxis(4)->SetRange(1,dphiHLambdaLSMixed->GetAxis(4)->GetNbins());
    dphiHLambdaLSMixed->GetAxis(5)->SetRangeUser(-10.0, 10.0);

    dphiHH->GetAxis(4)->SetRangeUser(-10.0, 10.0);
    dphiHHMixed->GetAxis(4)->SetRangeUser(-10.0, 10.0);

    Int_t axes[] = {2,3,4,5};

    THnSparseF* hLambda = (THnSparseF*)dphiHLambda->Projection(4, axes);
    THnSparseF* hLambdaLS = (THnSparseF*)dphiHLambdaLS->Projection(4, axes);

    TH3D* hh = (TH3D*)dphiHH->Projection(2,3,4);
    hh->Sumw2();

    THnSparseF* hLambdaMixed = (THnSparseF*)dphiHLambdaMixed->Projection(4, axes);
    THnSparseF* hLambdaLSMixed = (THnSparseF*)dphiHLambdaMixed->Projection(4, axes);

    TH3D* hhMixed = (TH3D*)dphiHHMixed->Projection(2,3,4);
    hhMixed->Sumw2();

    TH1D* sameUSPeakEta[10];
    TH1D* mixedUSPeakEta[10];
    TH1D* sameLSPeakEta[10];
    TH1D* mixedLSPeakEta[10];
    TH1D* sameUSRsideEta[10];
    TH1D* mixedUSRsideEta[10];
    TH1D* sameLSRsideEta[10];
    TH1D* mixedLSRsideEta[10];
    TH1D* sameUSLsideEta[10];
    TH1D* mixedUSLsideEta[10];
    TH1D* sameLSLsideEta[10];
    TH1D* mixedLSLsideEta[10];

    TH2D* hLambda2Dpeak = makeCorrections(hLambda, hLambdaMixed, SIG_MIN, SIG_MAX, sameUSPeakEta, mixedUSPeakEta, trigMixScalesUS, totalTrigSameUS);
    hLambda2Dpeak->SetName("hLambda2Dpeak");
    for(int i = 0; i<10; i++){
        sameUSPeakEta[i]->SetName(Form("sameUSPeakEta_zvtx_%i", i));
        mixedUSPeakEta[i]->SetName(Form("mixedUSPeakEta_zvtx_%i", i));
    }
    TH2D* hLambdaLS2Dpeak = makeCorrections(hLambdaLS, hLambdaLSMixed, SIG_MIN, SIG_MAX,  sameLSPeakEta, mixedLSPeakEta, trigMixScalesLS, totalTrigSameLS);
    hLambdaLS2Dpeak->SetName("hLambdaLS2Dpeak");
    for(int i = 0; i<10; i++){
        sameLSPeakEta[i]->SetName(Form("sameLSPeakEta_zvtx_%i", i));
        mixedLSPeakEta[i]->SetName(Form("mixedLSPeakEta_zvtx_%i", i));
    }
    TH2D* hLambda2DRside = makeCorrections(hLambda, hLambdaMixed, RSB_MIN, RSB_MAX, sameUSRsideEta, mixedUSRsideEta, trigMixScalesUS, totalTrigSameUS);
    hLambda2DRside->SetName("hLambda2DRside");
    for(int i = 0; i<10; i++){
        sameUSRsideEta[i]->SetName(Form("sameUSRsideEta_zvtx_%i", i));
        mixedUSRsideEta[i]->SetName(Form("mixedUSRsideEta_zvtx_%i", i));
    }
    TH2D* hLambdaLS2DRside = makeCorrections(hLambdaLS, hLambdaLSMixed, RSB_MIN, RSB_MAX, sameLSRsideEta, mixedLSRsideEta, trigMixScalesLS, totalTrigSameLS);
    hLambdaLS2DRside->SetName("hLambdaLS2DRside");
    for(int i = 0; i<10; i++){
        sameLSRsideEta[i]->SetName(Form("sameLSRsideEta_zvtx_%i", i));
        mixedLSRsideEta[i]->SetName(Form("mixedLSRsideEta_zvtx_%i", i));
    }
    TH2D* hLambda2DLside = makeCorrections(hLambda, hLambdaMixed, LSB_MIN, LSB_MAX, sameUSLsideEta, mixedUSLsideEta, trigMixScalesUS, totalTrigSameUS);
    hLambda2DLside->SetName("hLambda2DLside");
    for(int i = 0; i<10; i++){
        sameUSLsideEta[i]->SetName(Form("sameUSLsideEta_zvtx_%i", i));
        mixedUSLsideEta[i]->SetName(Form("mixedUSLsideEta_zvtx_%i", i));
    }
    TH2D* hLambdaLS2DLside = makeCorrections(hLambdaLS, hLambdaLSMixed, LSB_MIN, LSB_MAX, sameLSLsideEta, mixedLSLsideEta, trigMixScalesLS, totalTrigSameLS);
    hLambdaLS2DLside->SetName("hLambdaLS2DLside");
    for(int i = 0; i<10; i++){
        sameLSLsideEta[i]->SetName(Form("sameLSLsideEta_zvtx_%i", i));
        mixedLSLsideEta[i]->SetName(Form("mixedLSLsideEta_zvtx_%i", i));
    }

    // TH2D* hLambda2Dpeak = makeCorrections(hLambda, hLambdaMixed, 1.105, 1.125, sameUSPeakEta, mixedUSPeakEta, trigMixScalesUS, totalTrigSameUS);
    // hLambda2Dpeak->SetName("hLambda2Dpeak");
    // for(int i = 0; i<10; i++){
    //     sameUSPeakEta[i]->SetName(Form("sameUSPeakEta_zvtx_%i", i));
    //     mixedUSPeakEta[i]->SetName(Form("mixedUSPeakEta_zvtx_%i", i));
    // }
    // TH2D* hLambdaLS2Dpeak = makeCorrections(hLambdaLS, hLambdaLSMixed, 1.105, 1.125,  sameLSPeakEta, mixedLSPeakEta, trigMixScalesLS, totalTrigSameLS);
    // hLambdaLS2Dpeak->SetName("hLambdaLS2Dpeak");
    // for(int i = 0; i<10; i++){
    //     sameLSPeakEta[i]->SetName(Form("sameLSPeakEta_zvtx_%i", i));
    //     mixedLSPeakEta[i]->SetName(Form("mixedLSPeakEta_zvtx_%i", i));
    // }
    // TH2D* hLambda2DRside = makeCorrections(hLambda, hLambdaMixed, 1.14, 1.155, sameUSRsideEta, mixedUSRsideEta, trigMixScalesUS, totalTrigSameUS);
    // hLambda2DRside->SetName("hLambda2DRside");
    // for(int i = 0; i<10; i++){
    //     sameUSRsideEta[i]->SetName(Form("sameUSRsideEta_zvtx_%i", i));
    //     mixedUSRsideEta[i]->SetName(Form("mixedUSRsideEta_zvtx_%i", i));
    // }
    // TH2D* hLambdaLS2DRside = makeCorrections(hLambdaLS, hLambdaLSMixed, 1.14, 1.155, sameLSRsideEta, mixedLSRsideEta, trigMixScalesLS, totalTrigSameLS);
    // hLambdaLS2DRside->SetName("hLambdaLS2DRside");
    // for(int i = 0; i<10; i++){
    //     sameLSRsideEta[i]->SetName(Form("sameLSRsideEta_zvtx_%i", i));
    //     mixedLSRsideEta[i]->SetName(Form("mixedLSRsideEta_zvtx_%i", i));
    // }
    // TH2D* hLambda2DLside = makeCorrections(hLambda, hLambdaMixed, 1.08, 1.095, sameUSLsideEta, mixedUSLsideEta, trigMixScalesUS, totalTrigSameUS);
    // hLambda2DLside->SetName("hLambda2DLside");
    // for(int i = 0; i<10; i++){
    //     sameUSLsideEta[i]->SetName(Form("sameUSLsideEta_zvtx_%i", i));
    //     mixedUSLsideEta[i]->SetName(Form("mixedUSLsideEta_zvtx_%i", i));
    // }
    // TH2D* hLambdaLS2DLside = makeCorrections(hLambdaLS, hLambdaLSMixed, 1.08, 1.095, sameLSLsideEta, mixedLSLsideEta, trigMixScalesLS, totalTrigSameLS);
    // hLambdaLS2DLside->SetName("hLambdaLS2DLside");
    // for(int i = 0; i<10; i++){
    //     sameLSLsideEta[i]->SetName(Form("sameLSLsideEta_zvtx_%i", i));
    //     mixedLSLsideEta[i]->SetName(Form("mixedLSLsideEta_zvtx_%i", i));
    // }

    TH2D* hh2D = makehhCorrections(hh, hhMixed);
    hh2D->SetName("hh2D");
    hh2D->Scale(1.0/totalTrigSameUS);

	TH1D* hhdphi = hh2D->ProjectionY("hhdphi", hh2D->GetXaxis()->FindBin(-1.2), hh2D->GetXaxis()->FindBin(1.2));
    hhdphi->Scale(1.0/(hhdphi->Integral()));

    hh->GetZaxis()->SetRangeUser(-10.0, 10.0);
    TH2D* uncorrhh2D = (TH2D*)hh->Project3D("xye");
    uncorrhh2D->SetName("uncorrhh2D");
    uncorrhh2D->Scale(1.0/(uncorrhh2D->Integral(uncorrhh2D->GetXaxis()->FindBin(-1.2), uncorrhh2D->GetXaxis()->FindBin(1.2), 1, uncorrhh2D->GetYaxis()->GetNbins())));



    //Create some uncorrected same/mixed event 2D histos
    hLambda->GetAxis(2)->SetRangeUser(SIG_MIN, SIG_MAX);
    hLambdaMixed->GetAxis(2)->SetRangeUser(SIG_MIN, SIG_MAX);
    TH2D* uncorrhLambda2Dpeak = hLambda->Projection(0,1);
    uncorrhLambda2Dpeak->Sumw2();
    uncorrhLambda2Dpeak->SetName("uncorrhLambda2Dpeak");
    TH2D* uncorrhLambdaMixed2Dpeak = hLambdaMixed->Projection(0,1);
    uncorrhLambdaMixed2Dpeak->Sumw2();
    uncorrhLambdaMixed2Dpeak->SetName("uncorrhLambdaMixed2Dpeak");

    hLambda->GetAxis(2)->SetRangeUser(RSB_MIN, RSB_MAX);
    hLambdaMixed->GetAxis(2)->SetRangeUser(RSB_MIN, RSB_MAX);
    TH2D* uncorrhLambda2DRside = hLambda->Projection(0,1);
    uncorrhLambda2DRside->Sumw2();
    uncorrhLambda2DRside->SetName("uncorrhLambda2DRside");
    TH2D* uncorrhLambdaMixed2DRside = hLambdaMixed->Projection(0,1);
    uncorrhLambdaMixed2DRside->Sumw2();
    uncorrhLambdaMixed2DRside->SetName("uncorrhLambdaMixed2DRside");

    hLambda->GetAxis(2)->SetRangeUser(LSB_MIN, LSB_MAX);
    hLambdaMixed->GetAxis(2)->SetRangeUser(LSB_MIN, LSB_MAX);
    TH2D* uncorrhLambda2DLside = hLambda->Projection(0,1);
    uncorrhLambda2DLside->Sumw2();
    uncorrhLambda2DLside->SetName("uncorrhLambda2DLside");
    TH2D* uncorrhLambdaMixed2DLside = hLambdaMixed->Projection(0,1);
    uncorrhLambdaMixed2DLside->Sumw2();


    //make ratio plot of just 1 zvtx bin as a check:
    hLambdaMixed->GetAxis(2)->SetRangeUser(1.01, 1.03);
    hLambdaMixed->GetAxis(3)->SetRange(6,6);
    TH2D* mixedratioPeakZ2 = hLambdaMixed->Projection(0,1);
    TH2D* hist = hLambdaLSMixed->Projection(0,1);
    mixedratioPeakZ2->Divide(hist);
    TH1D* mixedratioPeakZ2deta = mixedratioPeakZ2->ProjectionX("mixedratioPeakZ2deta");


    //project just zvtx distributions for hLambdaMixed points and hLambdaLSMixed points in peak region:
    hLambdaMixed->GetAxis(3)->SetRange(0,0);
    hLambdaLSMixed->GetAxis(3)->SetRange(0,0);
    TH1D* mixedUSzvtx = hLambdaMixed->Projection(3);
    mixedUSzvtx->SetName("mixedUSzvtx");
    TH1D* mixedLSzvtx = hLambdaLSMixed->Projection(3);
    mixedLSzvtx->SetName("mixedLSzvtx");

    TFile* output = new TFile(Form("trig_%i_%i_assoc_%i_%i_cent_%i_%i_mixcorr_hLambda.root", (int)TRIG_PT_LOW, (int)TRIG_PT_HIGH, (int)ASSOC_PT_LOW, (int)ASSOC_PT_HIGH, (int)centLow, (int)centHigh), "RECREATE");
    
    hLambda2Dpeak->Scale(1.0/totalTrigSameUS);
    hLambda2DRside->Scale(1.0/totalTrigSameUS);
    hLambda2DLside->Scale(1.0/totalTrigSameUS);
    hLambdaLS2Dpeak->Scale(1.0/totalTrigSameUS);
    hLambdaLS2DRside->Scale(1.0/totalTrigSameUS);
    hLambdaLS2DLside->Scale(1.0/totalTrigSameUS);

    hLambda2Dpeak->Write();
    hLambda2DRside->Write();
    hLambda2DLside->Write();
    hLambdaLS2Dpeak->Write();
    hLambdaLS2DRside->Write();
    hLambdaLS2DLside->Write();
    uncorrhLambda2Dpeak->Write();
    uncorrhLambdaMixed2Dpeak->Write();
    uncorrhLambda2DRside->Write();
    uncorrhLambdaMixed2DRside->Write();
    uncorrhLambda2DLside->Write();
    uncorrhLambdaMixed2DLside->Write();

    mixedratioPeakZ2->Write();
    mixedratioPeakZ2deta->Write();

    mixedUSzvtx->Write();
    mixedLSzvtx->Write();

    for(int i=0; i< 10; i++){
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

    trigSameUSDist->Write();
    trigSameLSDist->Write();
	hh2D->Write();
    hh2D->Scale(1/totalTrigSameUS);
    hhdphi->Write();
    uncorrhh2D->Write();

}

Double_t bin2val(Int_t bin){
    Double_t val = (-0.5*TMath::Pi())+((Double_t)bin)*(2.0*TMath::Pi()/16.0);
    return val;
}

Double_t flineStd(Double_t *x, Double_t *par){
    if((x[0] > bin2val(1) && x[0] < bin2val(7)) || (x[0] > bin2val(9) && x[0] < bin2val(15))){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}

Double_t fline6bin(Double_t *x, Double_t *par){
    if((x[0] > bin2val(2) && x[0] < bin2val(6)) || (x[0] > bin2val(9) && x[0] < bin2val(15))){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}

Double_t flineNoLast(Double_t *x, Double_t *par){
    if((x[0] > bin2val(1) && x[0] < bin2val(7)) || x[0] > bin2val(9)){
        TF1::RejectPoint();
        return 0;
    }
    return par[0];
}


TH1D* getHisto(TString filename, TString histotype, TString histoname, TString mult, Float_t etamin, Float_t etamax, Int_t color, Int_t markerstyle){
    TFile *hhFile = new TFile(filename.Data());
    TH2D* histo2D = (TH2D*)hhFile->Get(histoname.Data());
    TString newhistoname = histoname +"_"+ mult;
    histo2D->SetName(newhistoname);
    histo2D->Sumw2();
    TString histoname1D = histotype + "dphi_" + mult;
    TH1D* histo1D = (TH1D*)histo2D->ProjectionY(histoname1D.Data(), histo2D->GetXaxis()->FindBin(etamin+ 0.001), histo2D->GetXaxis()->FindBin(etamax - 0.001));
    if(histotype == "hh" && mult == "50_100"){
        //histo1D->Rebin();
    }
    histo1D->SetLineWidth(2);
    histo1D->SetLineColor(color);
    histo1D->SetMarkerColor(color);
    histo1D->SetMarkerStyle(markerstyle);
    histo1D->SetMarkerSize(2);
    histo1D->GetXaxis()->SetTitle("#Delta#varphi");
    histo1D->SetTitle("");
    histo1D->GetYaxis()->SetTitleOffset(1.60);
    histo1D->GetYaxis()->SetMaxDigits(2);
    histo1D->GetXaxis()->SetTitleSize(0.05);
    histo1D->GetXaxis()->SetTitleOffset(0.90);
    return histo1D;
}

TF1* setupFit(TString fitname, TH1D* hist, Int_t color, Int_t linestyle, Int_t bgMethod){
    //Do straight line fits for systematics
    TF1 *linefit = new TF1(Form("line%s", fitname.Data()), flineStd, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use std. 4 points
    //TF1 *linefit = new TF1(Form("line%s", fitname.Data()), fline6bin, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use 6 points
    //TF1 *linefit = new TF1(Form("line%s", fitname.Data()), flineNoLast, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 1); //use 3 points (no last bin)

    hist->Fit(linefit);

    TF1 *basefit = new TF1(fitname, "gaus(0) + gaus(3) + ([3]/([5]))*exp(-(((x - [4] + 2.0*TMath::Pi())^2)/(2*[5]^2)))+ pol0(6)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

    switch(bgMethod){
        case 0: //4 bins for straight line
            basefit->FixParameter(6, 1.0/4.0*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(16)+hist->GetBinContent(1)));
            break;
        case 1: //6 bins for straight line
            basefit->FixParameter(6, (1.0/6.0)*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(16)+hist->GetBinContent(1)+hist->GetBinContent(2)+hist->GetBinContent(7)));
            break;
        case 2: //only around near peak (4)
            basefit->FixParameter(6, (1.0/4.0)*(hist->GetBinContent(8)+hist->GetBinContent(1)+hist->GetBinContent(2)+hist->GetBinContent(7)));
            break;
        case 3: //not last bin (3)
            basefit->FixParameter(6, (1.0/3.0)*(hist->GetBinContent(8)+hist->GetBinContent(9)+hist->GetBinContent(1)));
            break;
        case 4: //free fit
            basefit->SetParameter(0, hist->GetBinContent(hist->GetXaxis()->FindBin(0.0)) - basefit->GetParameter(6));
            basefit->SetParLimits(0, basefit->GetParameter(0)*0.1, basefit->GetParameter(0)*3.0);
            break;
        default: //fix straight line parameter to the line fit from above
            basefit->FixParameter(6, linefit->GetParameter(0));
            break;
    }

    /*if(fitname == "corrFit"){
        basefit->FixParameter(0, 1.68920E-01 - 3.31591E-04);
    }else if(fitname == "corrFit2050"){
        basefit->FixParameter(0, 1.69976E-01 - 3.03244E-04);
    }else if(fitname == "corrFit50100"){
        basefit->FixParameter(0, 1.73948E-01 - 3.76013E-04);
    }*/
    basefit->FixParameter(1, 0.0);
    //basefit->SetParameter(1, 0.0);
    basefit->SetParLimits(1, -0.5, 0.5);
    basefit->SetParameter(2, 0.5);
    basefit->SetParLimits(2, 0.1, 1.0);
    basefit->SetParameter(3, hist->GetBinContent(hist->GetXaxis()->FindBin(3.14)) - basefit->GetParameter(6));
    basefit->SetParLimits(3, basefit->GetParameter(3)*0.1, basefit->GetParameter(3)*3.0);
    basefit->FixParameter(4, 3.14159);
    //basefit->SetParameter(4, 3.14);
    basefit->SetParLimits(4, 3.0, 3.25);
    basefit->SetParameter(5, 0.5);
    basefit->SetParLimits(5, 0.1, 1.0);

    basefit->SetParNames("peak1", "mean1", "sigma1", "peak2", "mean2", "sigma2", "BG");

    basefit->SetLineColor(color);
    basefit->SetLineWidth(6);
    basefit->SetLineStyle(linestyle);
    return basefit;
}

void calculateRatios(TString input_0_20 = "", TString input_20_50 = "", TString input_50_100 = "", TString scaleSuffix = "avgscale", Int_t bgMethod=0){
    //canvas for probing problems

    TString outputstring = "lambda_hadron_ratios";
    
    TCanvas *testCanvas = new TCanvas("test", "TEST CANVAS PLS IGNORE", 0, 10, 1920, 1080);
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0);

    if(input_0_20.EqualTo("")) input_0_20 = "US_syst_trig_4_8_assoc_2_4_cent_0_20_mixcorr_hLambda.root";
    if(input_20_50.EqualTo("")) input_20_50 = "US_syst_trig_4_8_assoc_2_4_cent_20_50_mixcorr_hLambda.root";
    if(input_50_100.EqualTo("")) input_50_100 = "US_syst_trig_4_8_assoc_2_4_cent_50_100_mixcorr_hLambda.root";

    TH1D* hhdphi_0_20 = getHisto("trig_4_8_assoc_2_4_cent_0_20_mixcorr_hLambda.root", "hh", "hh2D", "Eff_0_20", -1.2, 1.2, kBlue+2, 21);
    TH1D* hLambdadphi_0_20 = getHisto(input_0_20, "hLambda", Form("RSUSsubhLambda2Dpeak%s", scaleSuffix.Data()), "Eff_0_20", -1.2, 1.2, kRed+2, 22);


    //scale for inv. mass range
    //hLambdadphi_0_20->Scale(1.0/0.897); //wide mass
    //hLambdadphi_0_20->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit = setupFit("corrFit", hhdphi_0_20, kBlue, 7, bgMethod);
    TF1 *corrFit2 = setupFit("corrFit2", hLambdadphi_0_20, kRed, 7, bgMethod);

    hhdphi_0_20->Fit("corrFit", "R0");
    hLambdadphi_0_20->Fit("corrFit2", "R0");

    TF1 *hhBG = new TF1("hhBG", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG->SetParLimits(0, 0.00001, 10000000.0);
    hhBG->SetParameter(0, 1.0*corrFit->GetParameter(6));
    hhBG->SetLineStyle(2);

    TF1 *hLambdaBG = new TF1("hLambdaBG", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hLambdaBG->SetParLimits(0, 0.00001, 10000000.0);
    hLambdaBG->SetParameter(0, 1.0*corrFit2->GetParameter(6));
    hLambdaBG->SetLineStyle(2);

    Double_t near0_20hLambdaError = 0;
    Double_t near0_20hhError = 0;
    Double_t away0_20hLambdaError = 0;
    Double_t away0_20hhError = 0;
    //Double_t mid0_20hLambdaError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hLambdadphi_0_20->GetBinError(8),2) + TMath::Power(hLambdadphi_0_20->GetBinError(16),2) + TMath::Power(hLambdadphi_0_20->GetBinError(1),2));
    Double_t mid0_20hLambdaError = corrFit2->GetParError(6);
    //Double_t mid0_20hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_0_20->GetBinError(8),2) + TMath::Power(hhdphi_0_20->GetBinError(16),2) + TMath::Power(hhdphi_0_20->GetBinError(1),2));
    Double_t mid0_20hhError = corrFit->GetParError(6);
    Double_t total0_20hLambdaError;
    Double_t total0_20hhError = 0;

    Double_t near0_20hLambdaYield = hLambdadphi_0_20->IntegralAndError(1,8,near0_20hLambdaError, "width") - hLambdaBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(1), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(8));
    //near0_20hLambdaError = TMath::Sqrt(TMath::Power(near0_20hLambdaError, 2) + TMath::Power(6.0*mid0_20hLambdaError, 2));
    Double_t near0_20hhYield = hhdphi_0_20->IntegralAndError(1,8,near0_20hhError, "width") - hhBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(1), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(8));
    //near0_20hhError = TMath::Sqrt(TMath::Power(near0_20hhError, 2) + TMath::Power(6.0*mid0_20hhError, 2));
    Double_t away0_20hLambdaYield = hLambdadphi_0_20->IntegralAndError(9,16,away0_20hLambdaError, "width") - hLambdaBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(9), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(16));
    //away0_20hLambdaError = TMath::Sqrt(TMath::Power(away0_20hLambdaError, 2) + TMath::Power(8.0*mid0_20hLambdaError, 2));
    Double_t away0_20hhYield = hhdphi_0_20->IntegralAndError(9,16,away0_20hhError, "width")- hhBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(9), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(16));
    //away0_20hhError = TMath::Sqrt(TMath::Power(away0_20hhError, 2) + TMath::Power(8.0*mid0_20hhError, 2));
    Double_t total0_20hLambdaYield = hLambdadphi_0_20->IntegralAndError(1,16,total0_20hLambdaError, "width");
    Double_t total0_20hhYield = hhdphi_0_20->IntegralAndError(1,16,total0_20hhError, "width");
    Double_t mid0_20hLambdaYield = hLambdaBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(1), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(16));
    mid0_20hLambdaError = mid0_20hLambdaError*16.0;
    Double_t mid0_20hhYield = hhBG->Integral(hLambdadphi_0_20->GetXaxis()->GetBinLowEdge(1), hLambdadphi_0_20->GetXaxis()->GetBinUpEdge(16));
    mid0_20hhError = mid0_20hhError*16.0;
    Double_t jet0_20hLambda = hLambdadphi_0_20->Integral(1, 16) - mid0_20hLambdaYield;
    Double_t jet0_20hh = hhdphi_0_20->Integral(1, 16) - mid0_20hhYield;

    Double_t near020 = near0_20hLambdaYield/near0_20hhYield;
    Double_t near020Er = near020*TMath::Sqrt(TMath::Power(near0_20hLambdaError/near0_20hLambdaYield, 2) + TMath::Power(near0_20hhError/near0_20hhYield, 2));
    Double_t away020 = away0_20hLambdaYield/away0_20hhYield;
    Double_t away020Er = away020*TMath::Sqrt(TMath::Power(away0_20hLambdaError/away0_20hLambdaYield, 2) + TMath::Power(away0_20hhError/away0_20hhYield, 2));
    Double_t mid020 = mid0_20hLambdaYield/mid0_20hhYield;
    Double_t mid020Er = mid020*TMath::Sqrt(TMath::Power(mid0_20hLambdaError/mid0_20hLambdaYield, 2) + TMath::Power(mid0_20hhError/mid0_20hhYield, 2));
    Double_t total020 = total0_20hLambdaYield/total0_20hhYield;
    Double_t total020Er = total020*TMath::Sqrt(TMath::Power(total0_20hLambdaError/total0_20hLambdaYield, 2) + TMath::Power(total0_20hhError/total0_20hhYield, 2));

    TH1D *ratios020 = new TH1D("ratios020", "(h-#Lambda^{0} / h-h) Ratios", 4, 0, 4);
    ratios020->GetXaxis()->SetBinLabel(1, "near-side");
    ratios020->SetBinContent(1, near020);
    ratios020->SetBinError(1, near020Er);
    ratios020->GetXaxis()->SetBinLabel(2, "mid");
    ratios020->SetBinContent(2, mid020);
    ratios020->SetBinError(2, mid020Er);
    ratios020->GetXaxis()->SetBinLabel(3, "away-side");
    ratios020->SetBinContent(3, away020);
    ratios020->SetBinError(3, away020Er);
    ratios020->GetXaxis()->SetBinLabel(4, "total");
    ratios020->SetBinContent(4, total020);
    ratios020->SetBinError(4, total020Er);
    ratios020->SetMarkerStyle(22);
    ratios020->SetMarkerColor(kRed+2);
    ratios020->SetLineColor(kRed+2);
    ratios020->SetMarkerSize(2);
    ratios020->SetLineWidth(2);
    testCanvas->cd();
    ratios020->Draw();

    TF1* paperfit020 = new TF2("paperfit020", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit020->SetParameter(0, hLambdaBG->GetParameter(0));


   //20-50 section
    TH1D* hhdphi_20_50 = getHisto("trig_4_8_assoc_2_4_cent_20_50_mixcorr_hLambda.root", "hh", "hh2D", "Eff_20_50", -1.2, 1.2, kBlue+2, 21);
    TH1D* hLambdadphi_20_50 = getHisto(input_20_50, "hLambda", Form("RSUSsubhLambda2Dpeak%s", scaleSuffix.Data()), "Eff_20_50", -1.2, 1.2, kRed+2, 22);


    //scale for inv. mass range
    //hLambdadphi_20_50->Scale(1.0/0.897); //wide mass
    //hLambdadphi_20_50->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit2050 = setupFit("corrFit2050", hhdphi_20_50, kBlue, 7, bgMethod);

    TF1 *corrFit2_2050 = setupFit("corrFit2_2050", hLambdadphi_20_50, kRed, 7, bgMethod);

    hhdphi_20_50->Fit("corrFit2050", "R0");
    hLambdadphi_20_50->Fit("corrFit2_2050", "R0");

    TF1 *hhBG_20_50 = new TF1("hhBG_20_50", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_20_50->SetParameter(0, 1.0*corrFit2050->GetParameter(6));
    hhBG_20_50->SetLineStyle(2);

    TF1 *hLambdaBG_20_50 = new TF1("hLambdaBG_20_50", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hLambdaBG_20_50->SetParLimits(0, 0.00001, 10000000.0);
    hLambdaBG_20_50->SetParameter(0, 1.0*corrFit2_2050->GetParameter(6));
    hLambdaBG_20_50->SetLineStyle(2);

    Double_t near20_50hLambdaError = 0;
    Double_t near20_50hhError = 0;
    Double_t away20_50hLambdaError = 0;
    Double_t away20_50hhError = 0;
    //Double_t mid20_50hLambdaError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hLambdadphi_20_50->GetBinError(8),2) + TMath::Power(hLambdadphi_20_50->GetBinError(16),2) + TMath::Power(hLambdadphi_20_50->GetBinError(1),2));
    Double_t mid20_50hLambdaError = corrFit2_2050->GetParError(6);
    //Double_t mid20_50hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_20_50->GetBinError(8),2) + TMath::Power(hhdphi_20_50->GetBinError(16),2) + TMath::Power(hhdphi_20_50->GetBinError(1),2));
    Double_t mid20_50hhError = corrFit2050->GetParError(6);
    Double_t total20_50hLambdaError = 0;
    Double_t total20_50hhError = 0;

    Double_t near20_50hLambdaYield = hLambdadphi_20_50->IntegralAndError(1,8,near20_50hLambdaError, "width") - hLambdaBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(1), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(8));
    //near20_50hLambdaError = TMath::Sqrt(TMath::Power(near20_50hLambdaError, 2) + TMath::Power(6.0*mid20_50hLambdaError, 2));
    Double_t near20_50hhYield = hhdphi_20_50->IntegralAndError(1,8,near20_50hhError, "width") - hhBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(1), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(8));
    //near20_50hhError = TMath::Sqrt(TMath::Power(near20_50hhError, 2) + TMath::Power(6.0*mid20_50hhError, 2));
    Double_t away20_50hLambdaYield = hLambdadphi_20_50->IntegralAndError(9,16,away20_50hLambdaError, "width") - hLambdaBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(9), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(16));
    //away20_50hLambdaError = TMath::Sqrt(TMath::Power(away20_50hLambdaError, 2) + TMath::Power(8.0*mid20_50hLambdaError, 2));
    Double_t away20_50hhYield = hhdphi_20_50->IntegralAndError(9,16,away20_50hhError, "width")- hhBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(9), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(16));
    //away20_50hhError = TMath::Sqrt(TMath::Power(away20_50hhError, 2) + TMath::Power(8.0*mid20_50hhError, 2));
    Double_t mid20_50hLambdaYield = hLambdaBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(1), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(16));
    mid20_50hLambdaError = mid20_50hLambdaError*16.0;
    Double_t mid20_50hhYield = hhBG_20_50->Integral(hLambdadphi_20_50->GetXaxis()->GetBinLowEdge(1), hLambdadphi_20_50->GetXaxis()->GetBinUpEdge(16));
    mid20_50hhError = mid20_50hhError*16.0;
    Double_t total20_50hLambdaYield = hLambdadphi_20_50->IntegralAndError(1, 16,total20_50hLambdaError, "width");
    Double_t total20_50hhYield = hhdphi_20_50->IntegralAndError(1, 16,total20_50hhError, "width");
    Double_t jet20_50hLambda = hLambdadphi_20_50->Integral(1, 16) - mid20_50hLambdaYield;
    Double_t jet20_50hh = hhdphi_20_50->Integral(1, 16) - mid20_50hhYield;


    Double_t near2050 = near20_50hLambdaYield/near20_50hhYield;
    Double_t near2050Er = near2050*TMath::Sqrt(TMath::Power(near20_50hLambdaError/near20_50hLambdaYield, 2) + TMath::Power(near20_50hhError/near20_50hhYield, 2));
    Double_t away2050 = away20_50hLambdaYield/away20_50hhYield;
    Double_t away2050Er = away2050*TMath::Sqrt(TMath::Power(away20_50hLambdaError/away20_50hLambdaYield, 2) + TMath::Power(away20_50hhError/away20_50hhYield, 2));
    Double_t mid2050 = mid20_50hLambdaYield/mid20_50hhYield;
    Double_t mid2050Er = mid2050*TMath::Sqrt(TMath::Power(mid20_50hLambdaError/mid20_50hLambdaYield, 2) + TMath::Power(mid20_50hhError/mid20_50hhYield, 2));
    Double_t total2050 = total20_50hLambdaYield/total20_50hhYield;
    Double_t total2050Er = total2050*TMath::Sqrt(TMath::Power(total20_50hLambdaError/total20_50hLambdaYield, 2) + TMath::Power(total20_50hhError/total20_50hhYield, 2));

    TH1D *ratios2050 = new TH1D("ratios2050", "(h-#Lambda / h-h) Ratios", 4, 0, 4);
    ratios2050->GetXaxis()->SetBinLabel(1, "near-side");
    ratios2050->SetBinContent(1, near2050);
    ratios2050->SetBinError(1, near2050Er);
    ratios2050->GetXaxis()->SetBinLabel(2, "mid");
    ratios2050->SetBinContent(2, mid2050);
    ratios2050->SetBinError(2, mid2050Er);
    ratios2050->GetXaxis()->SetBinLabel(3, "away-side");
    ratios2050->SetBinContent(3, away2050);
    ratios2050->SetBinError(3, away2050Er);
    ratios2050->GetXaxis()->SetBinLabel(4, "total");
    ratios2050->SetBinContent(4, total2050);
    ratios2050->SetBinError(4, total2050Er);
    ratios2050->SetMarkerStyle(21);
    ratios2050->SetMarkerColor(kBlue+2);
    ratios2050->SetLineColor(kBlue+2);
    ratios2050->SetMarkerSize(2);
    ratios2050->SetLineWidth(2);


    TF1* paperfit2050 = new TF2("paperfit2050", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit2050->SetParameter(0, hLambdaBG_20_50->GetParameter(0));

    //50-100 section
    TH1D* hhdphi_50_100 = getHisto("trig_4_8_assoc_2_4_cent_50_100_mixcorr_hLambda.root", "hh", "hh2D", "Eff_50_100", -1.2, 1.2, kBlue+2, 21);
    TH1D* hLambdadphi_50_100 = getHisto(input_50_100, "hLambda", Form("RSUSsubhLambda2Dpeak%s", scaleSuffix.Data()), "Eff_50_100", -1.2, 1.2, kRed+2, 22);



    //scale for inv. mass range
    //hLambdadphi_50_100->Scale(1.0/0.897); //wide mass
    //hLambdadphi_50_100->Scale(1.0/0.803); //narrow mass

    TF1 *corrFit50100 = setupFit("corrFit50100", hhdphi_50_100, kBlue, 7, bgMethod);

    TF1 *corrFit2_50100 = setupFit("corrFit2_50100", hLambdadphi_50_100, kRed, 7, bgMethod);

    hhdphi_50_100->Fit("corrFit50100", "R0");
    hLambdadphi_50_100->Fit("corrFit2_50100", "R0");

    TF1 *hhBG_50_100 = new TF1("hhBG_50_100", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hhBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hhBG_50_100->SetParameter(0, 1.0*corrFit50100->GetParameter(6));
    hhBG_50_100->SetLineStyle(2);

    TF1 *hLambdaBG_50_100 = new TF1("hLambdaBG_50_100", "pol0(0)", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    hLambdaBG_50_100->SetParLimits(0, 0.00001, 10000000.0);
    hLambdaBG_50_100->SetParameter(0, 1.0*corrFit2_50100->GetParameter(6));
    hLambdaBG_50_100->SetLineStyle(2);

    Double_t near50_100hLambdaError = 0;
    Double_t near50_100hhError = 0;
    Double_t away50_100hLambdaError = 0;
    Double_t away50_100hhError = 0;
    //Double_t mid50_100hLambdaError =(1.0/3.0)*TMath::Sqrt(TMath::Power(hLambdadphi_50_100->GetBinError(8),2) + TMath::Power(hLambdadphi_50_100->GetBinError(16),2) + TMath::Power(hLambdadphi_50_100->GetBinError(1),2));
    Double_t mid50_100hLambdaError = corrFit2_50100->GetParError(6);
    //Double_t mid50_100hhError = (1.0/3.0)*TMath::Sqrt(TMath::Power(hhdphi_50_100->GetBinError(8),2) + TMath::Power(hhdphi_50_100->GetBinError(16),2) + TMath::Power(hhdphi_50_100->GetBinError(1),2));
    Double_t mid50_100hhError = corrFit50100->GetParError(6);
    Double_t total50_100hLambdaError = 0;
    Double_t total50_100hhError = 0;

    Double_t near50_100hLambdaYield = hLambdadphi_50_100->IntegralAndError(1,8,near50_100hLambdaError, "width") - hLambdaBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(1), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(8));
    //near50_100hLambdaError = TMath::Sqrt(TMath::Power(near50_100hLambdaError, 2) + TMath::Power(6.0*mid50_100hLambdaError, 2));
    Double_t near50_100hhYield = hhdphi_50_100->IntegralAndError(1,8,near50_100hhError, "width") - hhBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(1), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(8));
    //near50_100hhError = TMath::Sqrt(TMath::Power(near50_100hhError, 2) + TMath::Power(6.0*mid50_100hhError, 2));
    Double_t away50_100hLambdaYield = hLambdadphi_50_100->IntegralAndError(9,16,away50_100hLambdaError, "width") - hLambdaBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(9), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(16));
    //away50_100hLambdaError = TMath::Sqrt(TMath::Power(away50_100hLambdaError, 2) + TMath::Power(8.0*mid50_100hLambdaError, 2));
    Double_t away50_100hhYield = hhdphi_50_100->IntegralAndError(9,16,away50_100hhError, "width")- hhBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(9), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(16));
    //away50_100hhError = TMath::Sqrt(TMath::Power(away50_100hhError, 2) + TMath::Power(8.0*mid50_100hhError, 2));
    Double_t mid50_100hLambdaYield = hLambdaBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(1), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(16));
    mid50_100hLambdaError = mid50_100hLambdaError*16.0;
    Double_t mid50_100hhYield = hhBG_50_100->Integral(hLambdadphi_50_100->GetXaxis()->GetBinLowEdge(1), hLambdadphi_50_100->GetXaxis()->GetBinUpEdge(16));
    mid50_100hhError = mid50_100hhError*16.0;
    Double_t total50_100hLambdaYield = hLambdadphi_50_100->IntegralAndError(1, 16,total50_100hLambdaError, "width");
    Double_t total50_100hhYield = hhdphi_50_100->IntegralAndError(1, 16,total50_100hhError, "width");
    Double_t jet50_100hLambda = hLambdadphi_50_100->Integral(1, 16) - mid50_100hLambdaYield;
    Double_t jet50_100hh = hhdphi_50_100->Integral(1, 16) - mid50_100hhYield;


    Double_t near50100 = near50_100hLambdaYield/near50_100hhYield;
    Double_t near50100Er = near50100*TMath::Sqrt(TMath::Power(near50_100hLambdaError/near50_100hLambdaYield, 2) + TMath::Power(near50_100hhError/near50_100hhYield, 2));
    Double_t away50100 = away50_100hLambdaYield/away50_100hhYield;
    Double_t away50100Er = away50100*TMath::Sqrt(TMath::Power(away50_100hLambdaError/away50_100hLambdaYield, 2) + TMath::Power(away50_100hhError/away50_100hhYield, 2));
    Double_t mid50100 = mid50_100hLambdaYield/mid50_100hhYield;
    Double_t mid50100Er = mid50100*TMath::Sqrt(TMath::Power(mid50_100hLambdaError/mid50_100hLambdaYield, 2) + TMath::Power(mid50_100hhError/mid50_100hhYield, 2));
    Double_t total50100 = total50_100hLambdaYield/total50_100hhYield;
    Double_t total50100Er = total50100*TMath::Sqrt(TMath::Power(total50_100hLambdaError/total50_100hLambdaYield, 2) + TMath::Power(total50_100hhError/total50_100hhYield, 2));

    TH1D *ratios50100 = new TH1D("ratios50100", "(h-#Lambda / h-h) Ratios", 4, 0, 4);
    ratios50100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios50100->SetBinContent(1, near50100);
    ratios50100->SetBinError(1, near50100Er);
    ratios50100->GetXaxis()->SetBinLabel(2, "mid");
    ratios50100->SetBinContent(2, mid50100);
    ratios50100->SetBinError(2, mid50100Er);
    ratios50100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios50100->SetBinContent(3, away50100);
    ratios50100->SetBinError(3, away50100Er);
    ratios50100->GetXaxis()->SetBinLabel(4, "total");
    ratios50100->SetBinContent(4, total50100);
    ratios50100->SetBinError(4, total50100Er);
    ratios50100->SetMarkerStyle(20);
    ratios50100->SetMarkerColor(kGreen+2);
    ratios50100->SetLineColor(kGreen+2);
    ratios50100->SetMarkerSize(2);
    ratios50100->SetLineWidth(2);

    TF1* paperfit50100 = new TF2("paperfit50100", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperfit50100->SetParameter(0, hLambdaBG_50_100->GetParameter(0));

    TF1* paperhhfit50100 = new TF2("paperhhfit50100", "[0]*(1+2*0.12*0.1*cos(2.0*x))", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    paperhhfit50100->SetParameter(0, hhBG_50_100->GetParameter(0));


    Double_t near0_100hLambdaYield = near0_20hLambdaYield + near20_50hLambdaYield + near50_100hLambdaYield;
    Double_t near0_100hLambdaError = TMath::Sqrt((TMath::Power(near0_20hLambdaError,2) + TMath::Power(near20_50hLambdaError,2) + TMath::Power(near50_100hLambdaError,2)));
    Double_t near0_100hhYield = near0_20hhYield + near20_50hhYield + near50_100hhYield;
    Double_t near0_100hhError = TMath::Sqrt((TMath::Power(near0_20hhError,2) + TMath::Power(near20_50hhError,2) + TMath::Power(near50_100hhError,2)));
    Double_t mid0_100hLambdaYield = mid0_20hLambdaYield + mid20_50hLambdaYield + mid50_100hLambdaYield;
    Double_t mid0_100hLambdaError = TMath::Sqrt((TMath::Power(mid0_20hLambdaError,2) + TMath::Power(mid20_50hLambdaError,2) + TMath::Power(mid50_100hLambdaError,2)));
    Double_t mid0_100hhYield = mid0_20hhYield + mid20_50hhYield + mid50_100hhYield;
    Double_t mid0_100hhError = TMath::Sqrt((TMath::Power(mid0_20hhError,2) + TMath::Power(mid20_50hhError,2) + TMath::Power(mid50_100hhError,2)));
    Double_t away0_100hLambdaYield = away0_20hLambdaYield + away20_50hLambdaYield + away50_100hLambdaYield;
    Double_t away0_100hLambdaError = TMath::Sqrt((TMath::Power(away0_20hLambdaError,2) + TMath::Power(away20_50hLambdaError,2) + TMath::Power(away50_100hLambdaError,2)));
    Double_t away0_100hhYield = away0_20hhYield + away20_50hhYield + away50_100hhYield;
    Double_t away0_100hhError = TMath::Sqrt((TMath::Power(away0_20hhError,2) + TMath::Power(away20_50hhError,2) + TMath::Power(away50_100hhError,2)));
    Double_t total0_100hLambdaYield = total0_20hLambdaYield + total20_50hLambdaYield + total50_100hLambdaYield;
    Double_t total0_100hLambdaError = TMath::Sqrt((TMath::Power(total0_20hLambdaError,2) + TMath::Power(total20_50hLambdaError,2) + TMath::Power(total50_100hLambdaError,2)));
    Double_t total0_100hhYield = total0_20hhYield + total20_50hhYield + total50_100hhYield;
    Double_t total0_100hhError = TMath::Sqrt((TMath::Power(total0_20hhError,2) + TMath::Power(total20_50hhError,2) + TMath::Power(total50_100hhError,2)));



    Double_t near0100 = (near0_20hLambdaYield + near20_50hLambdaYield + near50_100hLambdaYield)/(near0_20hhYield + near20_50hhYield + near50_100hhYield);
    Double_t near0100Er = near0100*TMath::Sqrt((TMath::Power(near0_20hLambdaError,2) + TMath::Power(near20_50hLambdaError,2) + TMath::Power(near50_100hLambdaError,2))/TMath::Power((near0_20hLambdaYield + near20_50hLambdaYield + near50_100hLambdaYield), 2) + (TMath::Power(near0_20hhError,2) + TMath::Power(near20_50hhError,2) + TMath::Power(near50_100hhError,2))/TMath::Power((near0_20hhYield + near20_50hhYield + near50_100hhYield), 2));
    Double_t away0100 = (away0_20hLambdaYield + away20_50hLambdaYield + away50_100hLambdaYield)/(away0_20hhYield + away20_50hhYield + away50_100hhYield);
    Double_t away0100Er = away0100*TMath::Sqrt((TMath::Power(away0_20hLambdaError,2) + TMath::Power(away20_50hLambdaError,2) + TMath::Power(away50_100hLambdaError,2))/TMath::Power((away0_20hLambdaYield + away20_50hLambdaYield + away50_100hLambdaYield), 2) + (TMath::Power(away0_20hhError,2) + TMath::Power(away20_50hhError,2) + TMath::Power(away50_100hhError,2))/TMath::Power((away0_20hhYield + away20_50hhYield + away50_100hhYield), 2));
    Double_t mid0100 = (mid0_20hLambdaYield + mid20_50hLambdaYield + mid50_100hLambdaYield)/(mid0_20hhYield + mid20_50hhYield + mid50_100hhYield);
    Double_t mid0100Er = mid0100*TMath::Sqrt((TMath::Power(mid0_20hLambdaError,2) + TMath::Power(mid20_50hLambdaError,2) + TMath::Power(mid50_100hLambdaError,2))/TMath::Power((mid0_20hLambdaYield + mid20_50hLambdaYield + mid50_100hLambdaYield), 2) + (TMath::Power(mid0_20hhError,2) + TMath::Power(mid20_50hhError,2) + TMath::Power(mid50_100hhError,2))/TMath::Power((mid0_20hhYield + mid20_50hhYield + mid50_100hhYield), 2));
    Double_t total0100 = (total0_20hLambdaYield + total20_50hLambdaYield + total50_100hLambdaYield)/(total0_20hhYield + total20_50hhYield + total50_100hhYield);
    Double_t total0100Er = total0100*TMath::Sqrt((TMath::Power(total0_20hLambdaError,2) + TMath::Power(total20_50hLambdaError,2) + TMath::Power(total50_100hLambdaError,2))/TMath::Power((total0_20hLambdaYield + total20_50hLambdaYield + total50_100hLambdaYield), 2) + (TMath::Power(total0_20hhError,2) + TMath::Power(total20_50hhError,2) + TMath::Power(total50_100hhError,2))/TMath::Power((total0_20hhYield + total20_50hhYield + total50_100hhYield), 2));

    printf("near020Er: %E\n", near020Er);
    printf("mid020Er: %E\n", mid020Er);
    printf("total0100Er: %E\n", total0100Er);
    printf("total020Er: %E\n", total020Er);
    printf("total0_20hLambdaError: %E\n", total0_20hLambdaError);
    TH1D *ratios0100 = new TH1D("ratios0100", "(h-#Lambda / h-h) Ratios", 4, 0, 4);
    ratios0100->GetXaxis()->SetBinLabel(1, "near-side");
    ratios0100->SetBinContent(1, near0100);
    ratios0100->SetBinError(1, near0100Er);
    ratios0100->GetXaxis()->SetBinLabel(2, "mid");
    ratios0100->SetBinContent(2, mid0100);
    ratios0100->SetBinError(2, mid0100Er);
    ratios0100->GetXaxis()->SetBinLabel(3, "away-side");
    ratios0100->SetBinContent(3, away0100);
    ratios0100->SetBinError(3, away0100Er);
    ratios0100->GetXaxis()->SetBinLabel(4, "total");
    ratios0100->SetBinContent(4, total0100);
    ratios0100->SetBinError(4, total0100Er);
    ratios0100->SetMarkerStyle(29);
    ratios0100->SetMarkerColor(kViolet-1);
    ratios0100->SetLineColor(kViolet-1);
    ratios0100->SetLineWidth(2);
    ratios0100->SetMarkerSize(2);

    TLegend  *ratioslegend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    ratioslegend->SetMargin(0.35);
    ratioslegend->AddEntry(ratios020, "0-20%", "p");
    ratioslegend->AddEntry(ratios2050, "20-50%", "p");
    ratioslegend->AddEntry(ratios50100, "50-100%", "p");
    ratioslegend->AddEntry(ratios0100, "0-100%", "p");

    TLine *line = new TLine(3.0, 0.0, 3.0, 0.0040);
    line->SetLineStyle(7);
    line->SetLineWidth(2);

    TCanvas *testc = new TCanvas("test", "test",50, 50, 600, 600);
    testc->cd();
    ratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    ratios2050->GetXaxis()->SetLabelSize(0.07);
    ratios2050->Draw("P SAME");
    ratios50100->GetXaxis()->SetLimits(0.05, 4.05);
    ratios50100->Draw("P SAME");
    ratios020->GetXaxis()->SetLimits(-0.05, 3.95);
    ratios020->Draw("P SAME");
    ratios0100->GetXaxis()->SetLimits(-0.1, 3.9);
    ratios0100->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");

    //initialize all sytematic error values
    Double_t near0_20hLambdaSystError = 0.068;
    Double_t near20_50hLambdaSystError = 0.077;
    Double_t near50_100hLambdaSystError = 0.133;
    Double_t away0_20hLambdaSystError = 0.097;
    Double_t away20_50hLambdaSystError = 0.061;
    Double_t away50_100hLambdaSystError = 0.111;
    Double_t mid0_20hLambdaSystError = 0.017;
    Double_t mid20_50hLambdaSystError = 0.013;
    Double_t mid50_100hLambdaSystError = 0.037;
    Double_t total0_20hLambdaSystError = 0.014;
    Double_t total20_50hLambdaSystError = 0.012;
    Double_t total50_100hLambdaSystError = 0.032;

    Double_t near0_20hhSystError = 0;
    Double_t near20_50hhSystError = 0;
    Double_t near50_100hhSystError = 0;
    Double_t away0_20hhSystError = 0;
    Double_t away20_50hhSystError = 0;
    Double_t away50_100hhSystError = 0;
    Double_t mid0_20hhSystError = 0;
    Double_t mid20_50hhSystError = 0;
    Double_t mid50_100hhSystError = 0;
    Double_t total0_20hhSystError = 0;
    Double_t total20_50hhSystError = 0;
    Double_t total50_100hhSystError = 0;


    //Setup single yield arrays for different regions
    Double_t nearhLambdaYieldArray[3] = {near50_100hLambdaYield*300.0, near20_50hLambdaYield*300.0, near0_20hLambdaYield*300.0};
    Double_t nearhLambdaYieldArrayErr[3] = {near50_100hLambdaError*300.0, near20_50hLambdaError*300.0, near0_20hLambdaError*300.0};
    Double_t nearhLambdaYieldArraySystErr[3] = {near50_100hLambdaYield*near50_100hLambdaSystError*300.0, near20_50hLambdaYield*near20_50hLambdaSystError*300.0, near0_20hLambdaYield*near0_20hLambdaSystError*300.0};
    Double_t awayhLambdaYieldArray[3] = {away50_100hLambdaYield*300.0, away20_50hLambdaYield*300.0, away0_20hLambdaYield*300.0};
    Double_t awayhLambdaYieldArrayErr[3] = {away50_100hLambdaError*300.0, away20_50hLambdaError*300.0, away0_20hLambdaError*300.0};
    Double_t awayhLambdaYieldArraySystErr[3] = {away50_100hLambdaYield*away50_100hLambdaSystError*300.0, away20_50hLambdaYield*away20_50hLambdaSystError*300.0, away0_20hLambdaYield*away0_20hLambdaSystError*300.0};
    Double_t bulkhLambdaYieldArray[3] = {mid50_100hLambdaYield*100.0, mid20_50hLambdaYield*100.0, mid0_20hLambdaYield*100.0};
    Double_t bulkhLambdaYieldArrayErr[3] = {mid50_100hLambdaError*100.0, mid20_50hLambdaError*100.0, mid0_20hLambdaError*100.0};
    Double_t bulkhLambdaYieldArraySystErr[3] = {mid50_100hLambdaYield*mid50_100hLambdaSystError*100.0, mid20_50hLambdaYield*mid20_50hLambdaSystError*100.0, mid0_20hLambdaYield*mid0_20hLambdaSystError*100.0};
    Double_t totalhLambdaYieldArray[3] = {total50_100hLambdaYield*100.0, total20_50hLambdaYield*100.0, total0_20hLambdaYield*100.0};
    Double_t totalhLambdaYieldArrayErr[3] = {total50_100hLambdaError*100.0, total20_50hLambdaError*100.0, total0_20hLambdaError*100.0};
    Double_t totalhLambdaYieldArraySystErr[3] = {total50_100hLambdaYield*total50_100hLambdaSystError*100.0, total20_50hLambdaYield*total20_50hLambdaSystError*100.0, total0_20hLambdaYield*total0_20hLambdaSystError*100.0};

    Double_t nearhhYieldArray[3] = {near50_100hhYield, near20_50hhYield, near0_20hhYield};
    Double_t nearhhYieldArrayErr[3] = {near50_100hhError, near20_50hhError, near0_20hhError};
    Double_t nearhhYieldArraySystErr[3] = {near50_100hhSystError, near20_50hhSystError, near0_20hhSystError};
    Double_t awayhhYieldArray[3] = {away50_100hhYield, away20_50hhYield, away0_20hhYield};
    Double_t awayhhYieldArrayErr[3] = {away50_100hhError, away20_50hhError, away0_20hhError};
    Double_t awayhhYieldArraySystErr[3] = {away50_100hhSystError, away20_50hhSystError, away0_20hhSystError};
    Double_t bulkhhYieldArray[3] = {mid50_100hhYield, mid20_50hhYield, mid0_20hhYield};
    Double_t bulkhhYieldArrayErr[3] = {mid50_100hhError, mid20_50hhError, mid0_20hhError};
    Double_t bulkhhYieldArraySystErr[3] = {mid50_100hhSystError, mid20_50hhSystError, mid0_20hhSystError};
    Double_t totalhhYieldArray[3] = {total50_100hhYield, total20_50hhYield, total0_20hhYield};
    Double_t totalhhYieldArrayErr[3] = {total50_100hhError, total20_50hhError, total0_20hhError};
    Double_t totalhhYieldArraySystErr[3] = {total50_100hhSystError, total20_50hhSystError, total0_20hhSystError};


    //Plot ratio as a function of multiplicity for the different angular regions
    Double_t nearArray[3] = {ratios50100->GetBinContent(1), ratios2050->GetBinContent(1), ratios020->GetBinContent(1)};
    Double_t nearArrayErr[3] = {ratios50100->GetBinError(1), ratios2050->GetBinError(1), ratios020->GetBinError(1)};
    Double_t awayArray[3] = {ratios50100->GetBinContent(3), ratios2050->GetBinContent(3), ratios020->GetBinContent(3)};
    Double_t awayArrayErr[3] = {ratios50100->GetBinError(3), ratios2050->GetBinError(3), ratios020->GetBinError(3)};
    Double_t bulkArray[3] = {ratios50100->GetBinContent(2), ratios2050->GetBinContent(2), ratios020->GetBinContent(2)};
    Double_t bulkArrayErr[3] = {ratios50100->GetBinError(2), ratios2050->GetBinError(2), ratios020->GetBinError(2)};
    Double_t totalArray[3] = {ratios50100->GetBinContent(4), ratios2050->GetBinContent(4), ratios020->GetBinContent(4)};
    Double_t totalArrayErr[3] = {ratios50100->GetBinError(4), ratios2050->GetBinError(4), ratios020->GetBinError(4)};

    //ratio of Near Jet to Underlying Event vs. Multiplicity
    //Double_t jet2UEhLambda[3] = {jet50_100hLambda/mid50_100hLambdaYield, jet20_50hLambda/mid20_50hLambdaYield, jet0_20hLambda/mid0_20hLambdaYield};
    //Double_t jet2UEhh[3] = {jet50_100hh/mid50_100hhYield, jet20_50hh/mid20_50hhYield, jet0_20hh/mid0_20hhYield};
    Double_t jet2UEhLambda[3] = {mid50_100hLambdaYield/total50_100hLambdaYield, mid20_50hLambdaYield/total20_50hLambdaYield, mid0_20hLambdaYield/total0_20hLambdaYield};
    Double_t jet2UEhh[3] = {mid50_100hhYield/total50_100hhYield, mid20_50hhYield/total20_50hhYield, mid0_20hhYield/total0_20hhYield};

    //systematic errors from the changing the fitting parameters
    Double_t nearArraySystErr[3] = {ratios50100->GetBinContent(1)*0.068, ratios2050->GetBinContent(1)*0.077, ratios020->GetBinContent(1)*0.133};
    Double_t awayArraySystErr[3] = {ratios50100->GetBinContent(3)*0.097, ratios2050->GetBinContent(3)*0.061, ratios020->GetBinContent(3)*0.111};
    Double_t bulkArraySystErr[3] = {ratios50100->GetBinContent(2)*0.017, ratios2050->GetBinContent(2)*0.013, ratios020->GetBinContent(2)*0.037};
    Double_t totalArraySystErr[3] = {ratios50100->GetBinContent(4)*0.014, ratios2050->GetBinContent(4)*0.012, ratios020->GetBinContent(4)*0.032};


    Double_t multArray[3] = {35.0, 65.0, 90.0};
    Double_t multArrayErr[3] = {15.0, 15.0, 10.0};
    Double_t multArraySystErr[3] = {2.5, 2.5, 2.5};

    Double_t mult2Array[3] = {36.0, 66.0, 91.0};
    Double_t mult2ArrayErr[3] = {15.0, 15.0, 10.0};

    //jet vs UE ratios
    TGraph* jetratioshLambda = new TGraphErrors(3, multArray, jet2UEhLambda);
    jetratioshLambda->SetMarkerStyle(20);
    jetratioshLambda->SetMarkerSize(2);
    jetratioshLambda->SetMarkerColor(kCyan+1);
    jetratioshLambda->SetLineColor(kCyan+2);
    jetratioshLambda->SetLineWidth(2);
    jetratioshLambda->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshLambda->GetXaxis()->SetTitleSize(0.05);
    jetratioshLambda->GetXaxis()->SetLabelSize(0.04);
    jetratioshLambda->GetXaxis()->SetTitleOffset(0.9);
    jetratioshLambda->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshLambda->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near-Side}{Underlying Event}#right)");
    jetratioshLambda->GetYaxis()->SetTitleSize(0.04);
    jetratioshLambda->GetYaxis()->SetTitleOffset(1.5);
    jetratioshLambda->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraph* jetratioshh = new TGraphErrors(3, multArray, jet2UEhh);
    jetratioshh->SetMarkerStyle(20);
    jetratioshh->SetMarkerSize(2);
    jetratioshh->SetMarkerColor(kOrange+1);
    jetratioshh->SetLineColor(kOrange+2);
    jetratioshh->SetLineWidth(2);
    jetratioshh->GetXaxis()->SetTitle("Multiplicity Percentile");
    jetratioshh->GetXaxis()->SetTitleSize(0.05);
    jetratioshh->GetXaxis()->SetLabelSize(0.04);
    jetratioshh->GetXaxis()->SetTitleOffset(0.9);
    jetratioshh->GetXaxis()->SetRangeUser(0.0, 100.0);
    jetratioshh->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Near-Side}{Underlying Event}#right)");
    jetratioshh->GetYaxis()->SetTitleSize(0.04);
    jetratioshh->GetYaxis()->SetTitleOffset(1.5);
    jetratioshh->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    //trying instead with variable sized histograms:
    Double_t binwidths[5] = {0.0, 20.0, 50.0, 80.0, 100.0};
    TH1D* ratioNearHist = new TH1D("ratioNearHist", "", 4, binwidths);
    TH1D* ratioBulkHist = new TH1D("ratioBulkHist", "", 4, binwidths);
    TH1D* ratioJetHist = new TH1D("ratioJetHist", "", 4, binwidths);
    for(int i =0; i<3; i++){
        ratioNearHist->SetBinContent(i+2, nearArray[i]);
        ratioNearHist->SetBinError(i+2, nearArrayErr[i]);
        ratioBulkHist->SetBinContent(i+2, bulkArray[i]);
        ratioBulkHist->SetBinError(i+2, bulkArrayErr[i]);
        ratioJetHist->SetBinContent(i+2, jet2UEhLambda[i]);
    }
    ratioNearHist->SetMarkerStyle(20);
    ratioNearHist->SetMarkerSize(2);
    ratioNearHist->SetMarkerColor(kRed+1);
    ratioNearHist->SetLineColor(kRed+2);
    ratioNearHist->SetLineWidth(2);
    ratioNearHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioNearHist->GetXaxis()->SetTitleSize(0.05);
    ratioNearHist->GetXaxis()->SetLabelSize(0.04);
    ratioNearHist->GetXaxis()->SetTitleOffset(1.2);
    ratioNearHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioNearHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    ratioNearHist->GetYaxis()->SetTitleSize(0.04);
    ratioNearHist->GetYaxis()->SetTitleOffset(1.5);
    ratioNearHist->GetYaxis()->SetRangeUser(0.0, 0.025);

    ratioJetHist->SetMarkerStyle(20);
    ratioJetHist->SetMarkerSize(2);
    ratioJetHist->SetMarkerColor(kCyan+1);
    ratioJetHist->SetLineColor(kCyan+2);
    ratioJetHist->SetLineWidth(2);
    ratioJetHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioJetHist->GetXaxis()->SetTitleSize(0.05);
    ratioJetHist->GetXaxis()->SetLabelSize(0.04);
    ratioJetHist->GetXaxis()->SetTitleOffset(1.2);
    ratioJetHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioJetHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{Jet Yield}{Underlying Event}#right)");
    ratioJetHist->GetYaxis()->SetTitleSize(0.04);
    ratioJetHist->GetYaxis()->SetTitleOffset(1.5);
    ratioJetHist->GetYaxis()->SetRangeUser(0.0, 0.6);


    ratioBulkHist->SetMarkerStyle(22);
    ratioBulkHist->SetMarkerSize(2);
    ratioBulkHist->SetMarkerColor(kGreen+2);
    ratioBulkHist->SetLineColor(kGreen+3);
    ratioBulkHist->SetLineWidth(2);
    ratioBulkHist->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratioBulkHist->GetXaxis()->SetTitleSize(0.05);
    ratioBulkHist->GetXaxis()->SetLabelSize(0.04);
    ratioBulkHist->GetXaxis()->SetTitleOffset(1.2);
    ratioBulkHist->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratioBulkHist->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    ratioBulkHist->GetYaxis()->SetTitleSize(0.04);
    ratioBulkHist->GetYaxis()->SetTitleOffset(1.5);
    ratioBulkHist->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* ratiosNear = new TGraphErrors(3, multArray, nearArray, multArrayErr, nearArrayErr);
    ratiosNear->SetMarkerStyle(20);
    ratiosNear->SetMarkerSize(2);
    ratiosNear->SetMarkerColor(kRed+1);
    ratiosNear->SetLineColor(kRed+2);
    ratiosNear->SetLineWidth(2);
    ratiosNear->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratiosNear->GetXaxis()->SetTitleSize(0.05);
    ratiosNear->GetXaxis()->SetLabelSize(0.04);
    ratiosNear->GetXaxis()->SetTitleOffset(0.9);
    ratiosNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratiosNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    ratiosNear->GetYaxis()->SetTitleSize(0.04);
    ratiosNear->GetYaxis()->SetTitleOffset(1.5);
    ratiosNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* ratiosNearSyst = new TGraphErrors(3, multArray, nearArray, multArraySystErr, nearArraySystErr);
    ratiosNearSyst->SetMarkerStyle(20);
    ratiosNearSyst->SetMarkerSize(1);
    ratiosNearSyst->SetMarkerColor(kRed+1);
    ratiosNearSyst->SetLineColor(kRed+3);
    ratiosNearSyst->SetFillColor(kWhite);
    ratiosNearSyst->SetLineWidth(2);
    ratiosNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    ratiosNearSyst->GetXaxis()->SetTitleSize(0.05);
    ratiosNearSyst->GetXaxis()->SetLabelSize(0.04);
    ratiosNearSyst->GetXaxis()->SetTitleOffset(0.9);
    ratiosNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    ratiosNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    ratiosNearSyst->GetYaxis()->SetTitleSize(0.04);
    ratiosNearSyst->GetYaxis()->SetTitleOffset(1.5);
    ratiosNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* ratiosAway = new TGraphErrors(3, multArray, awayArray, multArrayErr, awayArrayErr);
    ratiosAway->SetMarkerStyle(21);
    ratiosAway->SetMarkerSize(2);
    ratiosAway->SetMarkerColor(kBlue+1);
    ratiosAway->SetLineColor(kBlue+2);
    ratiosAway->SetLineWidth(2);

    TGraphErrors* ratiosAwaySyst = new TGraphErrors(3, multArray, awayArray, multArraySystErr, awayArraySystErr);
    ratiosAwaySyst->SetMarkerStyle(21);
    ratiosAwaySyst->SetMarkerSize(1);
    ratiosAwaySyst->SetMarkerColor(kBlue+1);
    ratiosAwaySyst->SetLineColor(kBlue+3);
    ratiosAwaySyst->SetLineWidth(2);
    //ratiosAwaySyst->SetFillColor(kWhite);

    TGraphErrors* ratiosBulk = new TGraphErrors(3, multArray, bulkArray, multArrayErr, bulkArrayErr);
    ratiosBulk->SetMarkerStyle(22);
    ratiosBulk->SetMarkerSize(2);
    ratiosBulk->SetMarkerColor(kGreen+2);
    ratiosBulk->SetLineColor(kGreen+3);
    ratiosBulk->SetLineWidth(2);

    TGraphErrors* ratiosBulkSyst = new TGraphErrors(3, multArray, bulkArray, multArraySystErr, bulkArraySystErr);
    ratiosBulkSyst->SetMarkerStyle(22);
    ratiosBulkSyst->SetMarkerSize(1);
    ratiosBulkSyst->SetMarkerColor(kGreen+2);
    ratiosBulkSyst->SetLineColor(kGreen+3);
    ratiosBulkSyst->SetLineWidth(2);


    TGraphErrors* ratiosTot = new TGraphErrors(3, multArray, totalArray, multArrayErr, totalArrayErr);
    ratiosTot->SetMarkerStyle(29);
    ratiosTot->SetMarkerSize(3);
    ratiosTot->SetMarkerColor(kMagenta+2);
    ratiosTot->SetLineColor(kMagenta+3);
    ratiosTot->SetLineWidth(2);
    ratiosTot->SetFillColor(kMagenta+1);
    ratiosTot->SetFillStyle(3144);

    TGraphErrors* ratiosTotSyst = new TGraphErrors(3, multArray, totalArray, multArraySystErr, totalArraySystErr);
    ratiosTotSyst->SetMarkerStyle(29);
    ratiosTotSyst->SetMarkerSize(3);
    ratiosTotSyst->SetMarkerColor(kMagenta+2);
    ratiosTotSyst->SetLineColor(kMagenta+3);
    ratiosTotSyst->SetLineWidth(2);
    ratiosTotSyst->SetFillColor(kMagenta+1);
    ratiosTotSyst->SetFillStyle(3144);


    //setting up scaled TGraph's for comparing Efficiency and Non-efficiency corrected
    TH1D* ratioNearHistScaled = (TH1D*)ratioNearHist->Clone("ratioNearHistScaled");
    TGraphErrors* ratiosNearScaled = (TGraphErrors*)ratiosNear->Clone("ratiosNearScaled");
    TGraphErrors* ratiosAwayScaled = (TGraphErrors*)ratiosAway->Clone("ratiosAwayScaled");
    TGraphErrors* ratiosTotScaled = (TGraphErrors*)ratiosTot->Clone("ratiosTotScaled");
    TGraphErrors* ratiosBulkScaled = (TGraphErrors*)ratiosBulk->Clone("ratiosBulkScaled");

    ratioNearHistScaled->Scale(1.0/totalArray[2]);
    Double_t x,y;
    for(int i = 0; i<3; i++){
       ratiosNearScaled->GetPoint(i, x, y);
       ratiosNearScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosNearScaled->SetPointError(i, ratiosNearScaled->GetErrorX(i), ratiosNearScaled->GetErrorY(i)/totalArray[2]);

       ratiosAwayScaled->GetPoint(i, x, y);
       ratiosAwayScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosAwayScaled->SetPointError(i, ratiosAwayScaled->GetErrorX(i), ratiosAwayScaled->GetErrorY(i)/totalArray[2]);

       ratiosBulkScaled->GetPoint(i, x, y);
       ratiosBulkScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosBulkScaled->SetPointError(i, ratiosBulkScaled->GetErrorX(i), ratiosBulkScaled->GetErrorY(i)/totalArray[2]);

       ratiosTotScaled->GetPoint(i, x, y);
       ratiosTotScaled->SetPoint(i, x, y/totalArray[2]);
       ratiosTotScaled->SetPointError(i, ratiosTotScaled->GetErrorX(i), ratiosTotScaled->GetErrorY(i)/totalArray[2]);
    }

    TLegend  *ratiosMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosMultlegend->SetMargin(0.35);
    ratiosMultlegend->AddEntry(ratiosBulk, "Underlying Event", "pl");
    ratiosMultlegend->AddEntry(ratiosAway, "Away-side (Jet)", "pl");
    ratiosMultlegend->AddEntry(ratiosNear, "Near-side (Jet)", "pl");
    ratiosMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "f");
    ratiosMultlegend->SetLineWidth(0);

    TLegend *ratiosUEMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosUEMultlegend->SetMargin(0.35);
    ratiosUEMultlegend->AddEntry(ratiosBulk, "In U.E.", "pl");
    ratiosUEMultlegend->AddEntry(ratiosTot, "Total (Jet + UE)", "f");
    ratiosUEMultlegend->SetLineWidth(0);

    TLegend *ratiosJetMultlegend = new TLegend(0.183, 0.686, 0.461, 0.928);
    ratiosJetMultlegend->SetMargin(0.35);
    ratiosJetMultlegend->AddEntry(jetratioshh, "Jet/U.E. for (h-h)", "pl");
    ratiosJetMultlegend->AddEntry(jetratioshLambda, "Jet/U.E. for (h-#Lambda)", "pl");
    ratiosJetMultlegend->SetLineWidth(0);


    TPaveText *text2 = new TPaveText(0.6687, 0.6585, 0.9282, 0.7979, "NDC");
    text2->AddText("trigger: 4.0 < p_{T}^{h} < 8.0 GeV/c");
    text2->AddText("assoc: 2.0 < p_{T}^{#Lambda} < 4.0 GeV/c");
    text2->SetTextSizePixels(18);
    text2->SetFillColor(kWhite);
    text2->SetBorderSize(0);


    TPaveText *data = new TPaveText(0.6982, 0.8258, 0.8987, 0.9251, "NDC");
    data->AddText("ALICE");
    data->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    data->GetLine(0)->SetTextSizePixels(32);
    data->GetLine(1)->SetTextSizePixels(24);
    data->SetBorderSize(0);
    data->SetFillColor(kWhite);

    TCanvas* vsMultCanvas = new TCanvas("vsMultCanvas", "vsMultCanvas", 55, 55, 900, 600);
    vsMultCanvas->cd();
    vsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioNearHist->Draw("PE");

    ratioNearHist->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHist->GetXaxis()->SetTickSize(0.0);

    gPad->Update();
    TGaxis* newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxis->Draw();
    // ratiosNearSyst->Draw("5");
    ratiosNear->Draw("P");
    // ratiosAwaySyst->Draw("5");
    ratiosAway->Draw("P");
    // ratiosBulkSyst->Draw("5");
    ratiosBulk->Draw("P");
    ratiosTot->Draw("2");
    // ratiosTotSyst->Draw("[]");
    //ratiosTot->Draw("3");
    ratiosMultlegend->Draw();
    data->Draw();
    text2->Draw();
    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();

    //Jet Ratio Canvas
    TCanvas* JetvsMultCanvas = new TCanvas("JetvsMultCanvas", "JetvsMultCanvas", 55, 55, 900, 600);
    JetvsMultCanvas->cd();
    JetvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("P");

    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    jetratioshh->Draw("P L");
    jetratioshLambda->Draw("P L");
    data->Draw();
    ratiosJetMultlegend->Draw("SAME");



    //scaled ratios
    TCanvas* vsMultCanvasScaled = new TCanvas("vsMultCanvasScaled", "vsMultCanvasScaled", 55, 55, 900, 600);
    vsMultCanvasScaled->cd();
    vsMultCanvasScaled->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioNearHistScaled->Draw("PE");

    ratioNearHistScaled->GetXaxis()->SetLabelOffset(999);
    //ratioNearHist->GetXaxis()->SetTitleOffset(999);
    ratioNearHistScaled->GetXaxis()->SetTickSize(0.0);

    //ratiosNear->Draw("P");
    gPad->Update();
    TGaxis *newaxisScaled = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHistScaled->GetXaxis()->GetXmin(),
            ratioNearHistScaled->GetXaxis()->GetXmax(),
            510,"-");
    newaxisScaled->SetLabelOffset(-0.03);
    //newaxis->SetTitle("Multipliciy % (VOA)");
    //newaxis->SetTitleOffset(1.3);
    newaxisScaled->Draw();
    //ratiosNearSyst->Draw("[]");
    ratiosNearScaled->Draw("P");
    //ratiosAwaySyst->Draw("[]");
    ratiosAwayScaled->Draw("P");
    ratiosBulkScaled->Draw("P");
    //ratiosJet->Draw("P");
    ratiosTotScaled->Draw("2");
    //ratiosTot->Draw("3");
    ratiosMultlegend->Draw();
    data->Draw();
    //ratiosNear->Draw("PL");
    //newaxis->Draw();
    //gPad->Update();

    //Just Underlying Event
    TCanvas* vsUEMultCanvas = new TCanvas("vsUEMultCanvas", "vsUEMultCanvas", 55, 55, 900, 600);
    vsUEMultCanvas->cd();
    vsUEMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    //TH1F* hist = ratiosNear->GetHistogram();
    gStyle->SetErrorX(0.5);
    ratioBulkHist->Draw("PE");

    ratioBulkHist->GetXaxis()->SetLabelOffset(999);
    ratioBulkHist->GetXaxis()->SetTitleOffset(999);
    ratioBulkHist->GetXaxis()->SetTickSize(0.0);

    gPad->Update();
    newaxis->Draw();
    //ratiosBulk->Draw("P");
    ratiosTot->Draw("2");
    ratiosUEMultlegend->Draw();
    data->Draw();


    //yield graphs
    TGraphErrors* yieldsNear = new TGraphErrors(3, multArray, nearhLambdaYieldArray, multArrayErr, nearhLambdaYieldArrayErr);
    yieldsNear->SetMarkerStyle(20);
    yieldsNear->SetMarkerSize(2);
    yieldsNear->SetMarkerColor(kRed+1);
    yieldsNear->SetLineColor(kRed+2);
    yieldsNear->SetLineWidth(2);
    yieldsNear->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldsNear->GetXaxis()->SetTitleSize(0.05);
    yieldsNear->GetXaxis()->SetLabelSize(0.04);
    yieldsNear->GetXaxis()->SetTitleOffset(0.9);
    yieldsNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldsNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    yieldsNear->GetYaxis()->SetTitleSize(0.04);
    yieldsNear->GetYaxis()->SetTitleOffset(1.5);
    yieldsNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldsNearSyst = new TGraphErrors(3, multArray, nearhLambdaYieldArray, multArraySystErr, nearhLambdaYieldArraySystErr);
    yieldsNearSyst->SetMarkerStyle(20);
    yieldsNearSyst->SetMarkerSize(1);
    yieldsNearSyst->SetMarkerColor(kRed+1);
    yieldsNearSyst->SetLineColor(kRed+3);
    yieldsNearSyst->SetFillColor(kWhite);
    yieldsNearSyst->SetLineWidth(2);
    yieldsNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldsNearSyst->GetXaxis()->SetTitleSize(0.05);
    yieldsNearSyst->GetXaxis()->SetLabelSize(0.04);
    yieldsNearSyst->GetXaxis()->SetTitleOffset(0.9);
    yieldsNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldsNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    yieldsNearSyst->GetYaxis()->SetTitleSize(0.04);
    yieldsNearSyst->GetYaxis()->SetTitleOffset(1.5);
    yieldsNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* yieldsAway = new TGraphErrors(3, multArray, awayhLambdaYieldArray, multArrayErr, awayhLambdaYieldArrayErr);
    yieldsAway->SetMarkerStyle(21);
    yieldsAway->SetMarkerSize(2);
    yieldsAway->SetMarkerColor(kBlue+1);
    yieldsAway->SetLineColor(kBlue+2);
    yieldsAway->SetLineWidth(2);

    TGraphErrors* yieldsAwaySyst = new TGraphErrors(3, multArray, awayhLambdaYieldArray, multArraySystErr, awayhLambdaYieldArraySystErr);
    yieldsAwaySyst->SetMarkerStyle(21);
    yieldsAwaySyst->SetMarkerSize(1);
    yieldsAwaySyst->SetMarkerColor(kBlue+1);
    yieldsAwaySyst->SetLineColor(kBlue+3);
    yieldsAwaySyst->SetLineWidth(2);
    yieldsAwaySyst->SetFillColor(kWhite);

    TGraphErrors* yieldsBulk = new TGraphErrors(3, multArray, bulkhLambdaYieldArray, multArrayErr, bulkhLambdaYieldArrayErr);
    yieldsBulk->SetMarkerStyle(22);
    yieldsBulk->SetMarkerSize(2);
    yieldsBulk->SetMarkerColor(kGreen+2);
    yieldsBulk->SetLineColor(kGreen+3);
    yieldsBulk->SetLineWidth(2);

    TGraphErrors* yieldsBulkSyst = new TGraphErrors(3, multArray, bulkhLambdaYieldArray, multArraySystErr, bulkhLambdaYieldArraySystErr);
    yieldsBulkSyst->SetMarkerStyle(22);
    yieldsBulkSyst->SetMarkerSize(1);
    yieldsBulkSyst->SetMarkerColor(kGreen+2);
    yieldsBulkSyst->SetLineColor(kGreen+3);
    yieldsBulkSyst->SetLineWidth(2);
    yieldsBulkSyst->SetFillColor(kWhite);


    TGraphErrors* yieldsTot = new TGraphErrors(3, multArray, totalhLambdaYieldArray, multArrayErr, totalhLambdaYieldArrayErr);
    yieldsTot->SetMarkerStyle(29);
    yieldsTot->SetMarkerSize(3);
    yieldsTot->SetMarkerColor(kMagenta+2);
    yieldsTot->SetLineColor(kMagenta+3);
    yieldsTot->SetLineWidth(2);
    yieldsTot->SetFillColor(kMagenta+1);
    yieldsTot->SetFillStyle(3144);

    TGraphErrors* yieldsTotSyst = new TGraphErrors(3, multArray, totalhLambdaYieldArray, multArrayErr, totalhLambdaYieldArraySystErr);
    yieldsTotSyst->SetMarkerStyle(29);
    yieldsTotSyst->SetMarkerSize(1);
    yieldsTotSyst->SetMarkerColor(kMagenta+2);
    yieldsTotSyst->SetLineColor(kMagenta+3);
    yieldsTotSyst->SetLineWidth(2);
    //yieldsTot->SetFillColor(kMagenta+1);
    //yieldsTot->SetFillStyle(3144);


    TGraphErrors* yieldshhNear = new TGraphErrors(3, multArray, nearhhYieldArray, multArrayErr, nearhhYieldArrayErr);
    yieldshhNear->SetMarkerStyle(24);
    yieldshhNear->SetMarkerSize(2);
    yieldshhNear->SetMarkerColor(kRed+2);
    yieldshhNear->SetLineColor(kRed+3);
    yieldshhNear->SetLineWidth(2);
    yieldshhNear->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldshhNear->GetXaxis()->SetTitleSize(0.05);
    yieldshhNear->GetXaxis()->SetLabelSize(0.04);
    yieldshhNear->GetXaxis()->SetTitleOffset(0.9);
    yieldshhNear->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldshhNear->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    yieldshhNear->GetYaxis()->SetTitleSize(0.04);
    yieldshhNear->GetYaxis()->SetTitleOffset(1.5);
    yieldshhNear->GetYaxis()->SetRangeUser(0.0002, 0.0035);

    TGraphErrors* yieldshhNearSyst = new TGraphErrors(3, multArray, nearhhYieldArray, multArrayErr, nearArraySystErr);
    yieldshhNearSyst->SetMarkerStyle(24);
    yieldshhNearSyst->SetMarkerSize(2);
    yieldshhNearSyst->SetMarkerColor(kRed+2);
    yieldshhNearSyst->SetLineColor(kRed+3);
    yieldshhNearSyst->SetFillColor(kWhite);
    yieldshhNearSyst->SetLineWidth(2);
    yieldshhNearSyst->GetXaxis()->SetTitle("Multiplicity Percentile");
    yieldshhNearSyst->GetXaxis()->SetTitleSize(0.05);
    yieldshhNearSyst->GetXaxis()->SetLabelSize(0.04);
    yieldshhNearSyst->GetXaxis()->SetTitleOffset(0.9);
    yieldshhNearSyst->GetXaxis()->SetRangeUser(0.0, 100.0);
    yieldshhNearSyst->GetYaxis()->SetTitle("Yield Ratio #left(#frac{h-#Lambda}{h-h}#right)");
    yieldshhNearSyst->GetYaxis()->SetTitleSize(0.04);
    yieldshhNearSyst->GetYaxis()->SetTitleOffset(1.5);
    yieldshhNearSyst->GetYaxis()->SetRangeUser(0.0002, 0.0035);



    TGraphErrors* yieldshhAway = new TGraphErrors(3, multArray, awayhhYieldArray, multArrayErr, awayhhYieldArrayErr);
    yieldshhAway->SetMarkerStyle(25);
    yieldshhAway->SetMarkerSize(2);
    yieldshhAway->SetMarkerColor(kBlue+2);
    yieldshhAway->SetLineColor(kBlue+3);
    yieldshhAway->SetLineWidth(2);

    TGraphErrors* yieldshhAwaySyst = new TGraphErrors(3, multArray, awayhhYieldArray, multArrayErr, awayhhYieldArraySystErr);
    yieldshhAwaySyst->SetMarkerStyle(25);
    yieldshhAwaySyst->SetMarkerSize(2);
    yieldshhAwaySyst->SetMarkerColor(kBlue+2);
    yieldshhAwaySyst->SetLineColor(kBlue+3);
    yieldshhAwaySyst->SetLineWidth(2);
    yieldshhAwaySyst->SetFillColor(kWhite);

    TGraphErrors* yieldshhBulk = new TGraphErrors(3, multArray, bulkhhYieldArray, multArrayErr, bulkhhYieldArrayErr);
    yieldshhBulk->SetMarkerStyle(26);
    yieldshhBulk->SetMarkerSize(2);
    yieldshhBulk->SetMarkerColor(kGreen+3);
    yieldshhBulk->SetLineColor(kGreen+4);
    yieldshhBulk->SetLineWidth(2);

    TGraphErrors* yieldshhBulkSyst = new TGraphErrors(3, multArray, bulkhhYieldArray, multArrayErr, bulkhhYieldArraySystErr);
    yieldshhBulkSyst->SetMarkerStyle(26);
    yieldshhBulkSyst->SetMarkerSize(2);
    yieldshhBulkSyst->SetMarkerColor(kGreen+3);
    yieldshhBulkSyst->SetLineColor(kGreen+4);
    yieldshhBulkSyst->SetLineWidth(2);
    yieldshhBulkSyst->SetFillColor(kWhite);


    TGraphErrors* yieldshhTot = new TGraphErrors(3, multArray, totalhhYieldArray, multArrayErr, totalhhYieldArrayErr);
    yieldshhTot->SetMarkerStyle(30);
    yieldshhTot->SetMarkerSize(3);
    yieldshhTot->SetMarkerColor(kMagenta+3);
    yieldshhTot->SetLineColor(kMagenta+4);
    yieldshhTot->SetLineWidth(2);
    yieldshhTot->SetFillColor(kMagenta+1);
    yieldshhTot->SetFillStyle(3144);

    TGraphErrors* yieldshhTotSyst = new TGraphErrors(3, multArray, totalhhYieldArray, multArrayErr, totalhhYieldArraySystErr);
    yieldshhTotSyst->SetMarkerStyle(30);
    yieldshhTotSyst->SetMarkerSize(3);
    yieldshhTotSyst->SetMarkerColor(kMagenta+3);
    yieldshhTotSyst->SetLineColor(kMagenta+4);
    yieldshhTotSyst->SetLineWidth(2);
    //yieldshhTot->SetFillColor(kMagenta+1);
    //yieldshhTot->SetFillStyle(3144);

    //set-up jet/total ratio histograms for h-phi and h-h

    TH1D* jet2totalhLambda = (TH1D*)ratioNearHist->Clone("jet2totalhLambda");
    TH1D* jet2totalhh = (TH1D*)ratioNearHist->Clone("jet2totalhh");
    for(int i = 0; i < 3; i++){
        jet2totalhLambda->SetBinContent(i+2, (nearhLambdaYieldArray[i] + awayhLambdaYieldArray[i])/totalhLambdaYieldArray[i]);
        jet2totalhh->SetBinContent(i+2, (nearhhYieldArray[i] + awayhhYieldArray[i])/totalhhYieldArray[i]);
    }

    jet2totalhLambda->SetLineColor(kMagenta+2);
    jet2totalhLambda->SetMarkerColor(kMagenta+2);
    jet2totalhLambda->GetYaxis()->SetTitle("(Jet/Total) pair-yield Ratio");

    jet2totalhh->SetLineColor(kCyan+2);
    jet2totalhh->SetMarkerColor(kCyan+2);

    TCanvas* jet2totalcanvas = new TCanvas("jet2totalcanvas", "jet2totalcanvas", 55, 55, 900, 600);
    jet2totalcanvas->cd();
    jet2totalcanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    jet2totalhLambda->Draw("AXIS");
    jet2totalhLambda->GetXaxis()->SetLabelOffset(999);
    jet2totalhLambda->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    jet2totalhLambda->Draw("P E SAME");
    jet2totalhh->Draw("P E SAME");

    //draw near-side yields
    TLegend *nearYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    nearYieldMultlegend->SetMargin(0.30);
    nearYieldMultlegend->AddEntry(yieldshhNear, "(h-h) in Near Jet", "pl");
    nearYieldMultlegend->AddEntry(yieldsNear, "(h-#Lambda) in Near Jet (x 300)", "pl");
    nearYieldMultlegend->SetLineWidth(0);


    TCanvas* nearYieldvsMultCanvas = new TCanvas("nearYieldvsMultCanvas", "nearYieldvsMultCanvas", 55, 55, 900, 600);
    nearYieldvsMultCanvas->cd();
    nearYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    yieldsNearSyst->Draw("5");
    yieldsNear->Draw("P");
    yieldshhNear->Draw("P");
    data->Draw();
    text2->Draw();
    nearYieldMultlegend->Draw("SAME");

    //draw away-side yields
    TLegend *awayYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    awayYieldMultlegend->SetMargin(0.30);
    awayYieldMultlegend->AddEntry(yieldshhAway, "(h-h) in Away Jet", "pl");
    awayYieldMultlegend->AddEntry(yieldsAway, "(h-#Lambda) in Away Jet (x 300)", "pl");
    awayYieldMultlegend->SetLineWidth(0);


    TCanvas* awayYieldvsMultCanvas = new TCanvas("awayYieldvsMultCanvas", "awayYieldvsMultCanvas", 55, 55, 900, 600);
    awayYieldvsMultCanvas->cd();
    awayYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    yieldsAwaySyst->Draw("5");
    yieldsAway->Draw("P");
    yieldshhAway->Draw("P");
    data->Draw();
    text2->Draw();
    awayYieldMultlegend->Draw("SAME");


    //draw bulk yields
    TLegend *bulkYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    bulkYieldMultlegend->SetMargin(0.30);
    bulkYieldMultlegend->AddEntry(yieldshhBulk, "(h-h) in UE", "pl");
    bulkYieldMultlegend->AddEntry(yieldsBulk, "(h-#Lambda) in UE (x 100)", "pl");
    bulkYieldMultlegend->SetLineWidth(0);


    TCanvas* bulkYieldvsMultCanvas = new TCanvas("bulkYieldvsMultCanvas", "bulkYieldvsMultCanvas", 55, 55, 900, 600);
    bulkYieldvsMultCanvas->cd();
    bulkYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    yieldsBulkSyst->Draw("5");
    yieldsBulk->Draw("P");
    yieldshhBulk->Draw("P");
    data->Draw();
    text2->Draw();
    bulkYieldMultlegend->Draw("SAME");

    //draw total yields
    TLegend *totalYieldMultlegend = new TLegend(0.1737, 0.781, 0.533, 0.919);
    totalYieldMultlegend->SetMargin(0.30);
    totalYieldMultlegend->AddEntry(yieldshhTot, "(h-h) pairs", "pl");
    totalYieldMultlegend->AddEntry(yieldsTot, "(h-#Lambda) pairs (x 100)", "pl");
    totalYieldMultlegend->SetLineWidth(0);


    TCanvas* totalYieldvsMultCanvas = new TCanvas("totalYieldvsMultCanvas", "totalYieldvsMultCanvas", 55, 55, 900, 600);
    totalYieldvsMultCanvas->cd();
    totalYieldvsMultCanvas->SetMargin(0.126, 0.05, 0.125, 0.05);
    gStyle->SetErrorX(0.5);
    ratioJetHist->Draw("AXIS");
    ratioJetHist->GetYaxis()->SetTitle("Per-trigger Pair Yields");
    ratioJetHist->GetXaxis()->SetLabelOffset(999);
    ratioJetHist->GetXaxis()->SetTickSize(0.0);
    gPad->Update();
    newaxis = new TGaxis(gPad->GetUxmax(),
            gPad->GetUymin(),
            gPad->GetUxmin(),
            gPad->GetUymin(),
            ratioNearHist->GetXaxis()->GetXmin(),
            ratioNearHist->GetXaxis()->GetXmax(),
            510,"-");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
    yieldsTotSyst->Draw("[]");
    yieldsTot->Draw("P");
    yieldshhTot->Draw("P");
    data->Draw();
    text2->Draw();
    totalYieldMultlegend->Draw("SAME");


    //Double Ratio plots:
    Double_t doublenear020 = near020/mid020;
    Double_t doublenear020Er = doublenear020*TMath::Sqrt(TMath::Power((near020Er)/(near020), 2) + TMath::Power((mid020Er)/(mid020), 2));
    Double_t doubleaway020 = away020/mid020;
    Double_t doubleaway020Er = doubleaway020*TMath::Sqrt(TMath::Power((away020Er)/(away020), 2) + TMath::Power((mid020Er)/(mid020), 2));

    Double_t doublenear2050 = near2050/mid2050;
    Double_t doublenear2050Er = doublenear2050*TMath::Sqrt(TMath::Power((near2050Er)/(near2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));
    Double_t doubleaway2050 = away2050/mid2050;
    Double_t doubleaway2050Er = doubleaway2050*TMath::Sqrt(TMath::Power((away2050Er)/(away2050), 2) + TMath::Power((mid2050Er)/(mid2050), 2));

    Double_t doublenear50100 = near50100/mid50100;
    Double_t doublenear50100Er = doublenear50100*TMath::Sqrt(TMath::Power((near50100Er)/(near50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));
    Double_t doubleaway50100 = away50100/mid50100;
    Double_t doubleaway50100Er = doubleaway50100*TMath::Sqrt(TMath::Power((away50100Er)/(away50100), 2) + TMath::Power((mid50100Er)/(mid50100), 2));

    Double_t doublenear0100 = near0100/mid0100;
    Double_t doublenear0100Er = doublenear0100*TMath::Sqrt(TMath::Power((near0100Er)/(near0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));
    Double_t doubleaway0100 = away0100/mid0100;
    Double_t doubleaway0100Er = doubleaway0100*TMath::Sqrt(TMath::Power((away0100Er)/(away0100), 2) + TMath::Power((mid0100Er)/(mid0100), 2));

    TH1D *doubleratios020 = new TH1D("doubleratios020", "(h-#Lambda / h-h) Double Ratios", 2, 0, 2);
    doubleratios020->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios020->SetBinContent(1, doublenear020);
    doubleratios020->SetBinError(1, doublenear020Er);
    doubleratios020->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios020->SetBinContent(2, doubleaway020);
    doubleratios020->SetBinError(2, doubleaway020Er);
    doubleratios020->SetMarkerStyle(22);
    doubleratios020->SetMarkerColor(kRed+2);
    doubleratios020->SetLineColor(kRed+2);
    doubleratios020->SetLineWidth(2);
    doubleratios020->SetMarkerSize(2);

    TH1D *doubleratios2050 = new TH1D("doubleratios2050", "(h-#Lambda / h-h) Double Ratios", 2, 0, 2);
    doubleratios2050->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios2050->SetBinContent(1, doublenear2050);
    doubleratios2050->SetBinError(1, doublenear2050Er);
    doubleratios2050->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios2050->SetBinContent(2, doubleaway2050);
    doubleratios2050->SetBinError(2, doubleaway2050Er);
    doubleratios2050->SetMarkerStyle(21);
    doubleratios2050->SetMarkerColor(kBlue+2);
    doubleratios2050->SetLineColor(kBlue+2);
    doubleratios2050->SetLineWidth(2);
    doubleratios2050->SetMarkerSize(2);

    TH1D *doubleratios50100 = new TH1D("doubleratios50100", "(h-#Lambda / h-h) Double Ratios", 2, 0, 2);
    doubleratios50100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios50100->SetBinContent(1, doublenear50100);
    doubleratios50100->SetBinError(1, doublenear50100Er);
    doubleratios50100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios50100->SetBinContent(2, doubleaway50100);
    doubleratios50100->SetBinError(2, doubleaway50100Er);
    doubleratios50100->SetMarkerStyle(20);
    doubleratios50100->SetMarkerColor(kGreen+2);
    doubleratios50100->SetLineColor(kGreen+2);
    doubleratios50100->SetLineWidth(2);
    doubleratios50100->SetMarkerSize(2);

    TH1D *doubleratios0100 = new TH1D("doubleratios0100", "(h-#Lambda / h-h) Double Ratios", 2, 0, 2);
    doubleratios0100->GetXaxis()->SetBinLabel(1, "#frac{near-side}{mid}");
    doubleratios0100->SetBinContent(1, doublenear0100);
    doubleratios0100->SetBinError(1, doublenear0100Er);
    doubleratios0100->GetXaxis()->SetBinLabel(2, "#frac{away-side}{mid}");
    doubleratios0100->SetBinContent(2, doubleaway0100);
    doubleratios0100->SetBinError(1, doubleaway0100Er);
    doubleratios0100->SetMarkerStyle(29);
    doubleratios0100->SetMarkerColor(kViolet-1);
    doubleratios0100->SetLineColor(kViolet-1);
    doubleratios0100->SetLineWidth(2);
    doubleratios0100->SetMarkerSize(2);


    TCanvas *testDoublec = new TCanvas("testdouble", "testdouble",50, 50, 600, 600);
    testDoublec->cd();
    //doubleratios2050->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    doubleratios2050->GetXaxis()->SetLabelSize(0.07);
    doubleratios2050->Draw("P SAME");
    doubleratios50100->GetXaxis()->SetLimits(-0.05, 1.95);
    doubleratios50100->Draw("P SAME");
    doubleratios020->GetXaxis()->SetLimits(0.05, 2.05);
    doubleratios020->Draw("P SAME");
    doubleratios0100->GetXaxis()->SetLimits(-0.1,1.9);
    doubleratios0100->Draw("P SAME");
    //line->Draw("SAME");
    ratioslegend->Draw("SAME");

    printf("\n\n");
    printf(" h-Lambda  YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hLambdaYield, near20_50hLambdaYield, near50_100hLambdaYield, near0_100hLambdaYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hLambdaYield, mid20_50hLambdaYield, mid50_100hLambdaYield, mid0_100hLambdaYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hLambdaYield, away20_50hLambdaYield, away50_100hLambdaYield, away0_100hLambdaYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hLambdaYield, total20_50hLambdaYield, total50_100hLambdaYield, total0_100hLambdaYield);

    printf("\n\n");
    printf(" h-Lambda  ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hLambdaError, near20_50hLambdaError, near50_100hLambdaError, near0_100hLambdaError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hLambdaError, mid20_50hLambdaError, mid50_100hLambdaError, mid0_100hLambdaError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hLambdaError, away20_50hLambdaError, away50_100hLambdaError, away0_100hLambdaError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hLambdaError, total20_50hLambdaError, total50_100hLambdaError, total0_100hLambdaError);

    printf("\n\n");
    printf(" h-Lambda %%ERRORS ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*near0_20hLambdaError/near0_20hLambdaYield, 100*near20_50hLambdaError/near20_50hLambdaYield, 100*near50_100hLambdaError/near50_100hLambdaYield, 100*near0_100hLambdaError/near0_100hLambdaYield);
    printf("      MID      ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*mid0_20hLambdaError/mid0_20hLambdaYield, 100*mid20_50hLambdaError/mid20_50hLambdaYield, 100*mid50_100hLambdaError/mid50_100hLambdaYield, 100*mid0_100hLambdaError/mid0_100hLambdaYield);
    printf("      AWAY     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||     %2.2f%%     ||\n", 100*away0_20hLambdaError/away0_20hLambdaYield, 100*away20_50hLambdaError/away20_50hLambdaYield, 100*away50_100hLambdaError/away50_100hLambdaYield, 100*away0_100hLambdaError/away0_100hLambdaYield);
    printf("      TOTAL    ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||      %2.2f%%     ||\n", 100*total0_20hLambdaError/total0_20hLambdaYield, 100*total20_50hLambdaError/total20_50hLambdaYield, 100*total50_100hLambdaError/total50_100hLambdaYield, 100*total0_100hLambdaError/total0_100hLambdaYield);



    printf("\n\n");
    printf(" h-h   YIELDS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hhYield, near20_50hhYield, near50_100hhYield, near0_100hhYield);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hhYield, mid20_50hhYield, mid50_100hhYield, mid0_100hhYield);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hhYield, away20_50hhYield, away50_100hhYield, away0_100hhYield);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hhYield, total20_50hhYield, total50_100hhYield, total0_100hhYield);

    printf("\n\n");
    printf(" h-h ERRORS  ||     0 - 20     ||    20 - 50     ||    50 - 100    ||     0 - 100    ||\n");
    printf("=========================================================================================\n");
    printf("      NEAR     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", near0_20hhError, near20_50hhError, near50_100hhError, near0_100hhError);
    printf("      MID      ||  %E  ||  %E  ||  %E  ||  %E  ||\n", mid0_20hhError, mid20_50hhError, mid50_100hhError, mid0_100hhError);
    printf("      AWAY     ||  %E  ||  %E  ||  %E  ||  %E  ||\n", away0_20hhError, away20_50hhError, away50_100hhError, away0_100hhError);
    printf("      TOTAL    ||  %E  ||  %E  ||  %E  ||  %E  ||\n", total0_20hhError, total20_50hhError, total50_100hhError, total0_100hhError);
    printf("\n\n");


    TH1D *yields020hLambda = new TH1D("yields020hLambdaEff", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields020hLambda->GetXaxis()->SetBinLabel(1, "near-side");
    yields020hLambda->SetBinContent(1, near0_20hLambdaYield);
    yields020hLambda->SetBinError(1, near0_20hLambdaError);
    yields020hLambda->GetXaxis()->SetBinLabel(2, "mid");
    yields020hLambda->SetBinContent(2, mid0_20hLambdaYield);
    yields020hLambda->SetBinError(2, mid0_20hLambdaError);
    yields020hLambda->GetXaxis()->SetBinLabel(3, "away-side");
    yields020hLambda->SetBinContent(3, away0_20hLambdaYield);
    yields020hLambda->SetBinError(3, away0_20hLambdaError);
    yields020hLambda->GetXaxis()->SetBinLabel(4, "total");
    yields020hLambda->SetBinContent(4, total0_20hLambdaYield);
    yields020hLambda->SetBinError(4, total0_20hLambdaError);
    yields020hLambda->SetMarkerStyle(22);
    yields020hLambda->SetMarkerColor(kRed+2);
    yields020hLambda->SetLineColor(kRed+2);
    yields020hLambda->SetMarkerSize(2);
    yields020hLambda->SetLineWidth(2);

    TH1D *yields2050hLambda = new TH1D("yields2050hLambdaEff", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields2050hLambda->GetXaxis()->SetBinLabel(1, "near-side");
    yields2050hLambda->SetBinContent(1, near20_50hLambdaYield);
    yields2050hLambda->SetBinError(1, near20_50hLambdaError);
    yields2050hLambda->GetXaxis()->SetBinLabel(2, "mid");
    yields2050hLambda->SetBinContent(2, mid20_50hLambdaYield);
    yields2050hLambda->SetBinError(2, mid20_50hLambdaError);
    yields2050hLambda->GetXaxis()->SetBinLabel(3, "away-side");
    yields2050hLambda->SetBinContent(3, away20_50hLambdaYield);
    yields2050hLambda->SetBinError(3, away20_50hLambdaError);
    yields2050hLambda->GetXaxis()->SetBinLabel(4, "total");
    yields2050hLambda->SetBinContent(4, total20_50hLambdaYield);
    yields2050hLambda->SetBinError(4, total20_50hLambdaError);
    yields2050hLambda->SetMarkerStyle(21);
    yields2050hLambda->SetMarkerColor(kBlue+2);
    yields2050hLambda->SetLineColor(kBlue+2);
    yields2050hLambda->SetMarkerSize(2);
    yields2050hLambda->SetLineWidth(2);

    TH1D *yields50100hLambda = new TH1D("yields50100hLambdaEff", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields50100hLambda->GetXaxis()->SetBinLabel(1, "near-side");
    yields50100hLambda->SetBinContent(1, near50_100hLambdaYield);
    yields50100hLambda->SetBinError(1, near50_100hLambdaError);
    yields50100hLambda->GetXaxis()->SetBinLabel(2, "mid");
    yields50100hLambda->SetBinContent(2, mid50_100hLambdaYield);
    yields50100hLambda->SetBinError(2, mid50_100hLambdaError);
    yields50100hLambda->GetXaxis()->SetBinLabel(3, "away-side");
    yields50100hLambda->SetBinContent(3, away50_100hLambdaYield);
    yields50100hLambda->SetBinError(3, away50_100hLambdaError);
    yields50100hLambda->GetXaxis()->SetBinLabel(4, "total");
    yields50100hLambda->SetBinContent(4, total50_100hLambdaYield);
    yields50100hLambda->SetBinError(4, total50_100hLambdaError);
    yields50100hLambda->SetMarkerStyle(20);
    yields50100hLambda->SetMarkerColor(kGreen+2);
    yields50100hLambda->SetLineColor(kGreen+2);
    yields50100hLambda->SetMarkerSize(2);
    yields50100hLambda->SetLineWidth(2);

    TH1D *yields0100hLambda = new TH1D("yields0100hLambda", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields0100hLambda->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hLambda->SetBinContent(1, near0_100hLambdaYield);
    yields0100hLambda->SetBinError(1, near0_100hLambdaError);
    yields0100hLambda->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hLambda->SetBinContent(2, mid0_100hLambdaYield);
    yields0100hLambda->SetBinError(2, mid0_100hLambdaError);
    yields0100hLambda->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hLambda->SetBinContent(3, away0_100hLambdaYield);
    yields0100hLambda->SetBinError(3, away0_100hLambdaError);
    yields0100hLambda->GetXaxis()->SetBinLabel(4, "total");
    yields0100hLambda->SetBinContent(4, total0_100hLambdaYield);
    yields0100hLambda->SetBinError(4, total0_100hLambdaError);
    yields0100hLambda->SetMarkerStyle(29);
    yields0100hLambda->SetMarkerColor(kViolet-1);
    yields0100hLambda->SetLineColor(kViolet-1);
    yields0100hLambda->SetMarkerSize(2);
    yields0100hLambda->SetLineWidth(2);


    TCanvas *yieldhLambdac = new TCanvas("yieldhLambda", "yieldhLambda",50, 50, 600, 600);
    yieldhLambdac->cd();
    yields2050hLambda->GetYaxis()->SetRangeUser(0.0000, 0.0040);
    yields2050hLambda->GetXaxis()->SetLabelSize(0.07);
    yields2050hLambda->Draw("P SAME");
    yields50100hLambda->GetXaxis()->SetLimits(0.05, 4.05);
    yields50100hLambda->Draw("P SAME");
    yields020hLambda->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hLambda->Draw("P SAME");
    yields0100hLambda->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hLambda->Draw("P SAME");
    line->Draw("SAME");
    ratioslegend->Draw("SAME");


    TH1D *yields020hh = new TH1D("yields020hhEff", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields020hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields020hh->SetBinContent(1, near0_20hhYield);
    yields020hh->SetBinError(1, near0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(2, "mid");
    yields020hh->SetBinContent(2, mid0_20hhYield);
    yields020hh->SetBinError(2, mid0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields020hh->SetBinContent(3, away0_20hhYield);
    yields020hh->SetBinError(3, away0_20hhError);
    yields020hh->GetXaxis()->SetBinLabel(4, "total");
    yields020hh->SetBinContent(4, total0_20hhYield);
    yields020hh->SetBinError(4, total0_20hhError);
    yields020hh->SetMarkerStyle(22);
    yields020hh->SetMarkerColor(kRed+2);
    yields020hh->SetLineColor(kRed+2);
    yields020hh->SetMarkerSize(2);
    yields020hh->SetLineWidth(2);

    TH1D *yields2050hh = new TH1D("yields2050hhEff", "h-h Per Trigger Yields", 4, 0, 4);
    yields2050hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields2050hh->SetBinContent(1, near20_50hhYield);
    yields2050hh->SetBinError(1, near20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(2, "mid");
    yields2050hh->SetBinContent(2, mid20_50hhYield);
    yields2050hh->SetBinError(2, mid20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields2050hh->SetBinContent(3, away20_50hhYield);
    yields2050hh->SetBinError(3, away20_50hhError);
    yields2050hh->GetXaxis()->SetBinLabel(4, "total");
    yields2050hh->SetBinContent(4, total20_50hhYield);
    yields2050hh->SetBinError(4, total20_50hhError);
    yields2050hh->SetMarkerStyle(21);
    yields2050hh->SetMarkerColor(kBlue+2);
    yields2050hh->SetLineColor(kBlue+2);
    yields2050hh->SetMarkerSize(2);
    yields2050hh->SetLineWidth(2);

    TH1D *yields50100hh = new TH1D("yields50100hhEff", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields50100hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields50100hh->SetBinContent(1, near50_100hhYield);
    yields50100hh->SetBinError(1, near50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(2, "mid");
    yields50100hh->SetBinContent(2, mid50_100hhYield);
    yields50100hh->SetBinError(2, mid50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields50100hh->SetBinContent(3, away50_100hhYield);
    yields50100hh->SetBinError(3, away50_100hhError);
    yields50100hh->GetXaxis()->SetBinLabel(4, "total");
    yields50100hh->SetBinContent(4, total50_100hhYield);
    yields50100hh->SetBinError(4, total50_100hhError);
    yields50100hh->SetMarkerStyle(20);
    yields50100hh->SetMarkerColor(kGreen+2);
    yields50100hh->SetLineColor(kGreen+2);
    yields50100hh->SetMarkerSize(2);
    yields50100hh->SetLineWidth(2);

    TH1D *yields0100hh = new TH1D("yields0100hh", "h-#Lambda Per Trigger Yields", 4, 0, 4);
    yields0100hh->GetXaxis()->SetBinLabel(1, "near-side");
    yields0100hh->SetBinContent(1, near0_100hhYield);
    yields0100hh->SetBinError(1, near0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(2, "mid");
    yields0100hh->SetBinContent(2, mid0_100hhYield);
    yields0100hh->SetBinError(2, mid0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(3, "away-side");
    yields0100hh->SetBinContent(3, away0_100hhYield);
    yields0100hh->SetBinError(3, away0_100hhError);
    yields0100hh->GetXaxis()->SetBinLabel(4, "total");
    yields0100hh->SetBinContent(4, total0_100hhYield);
    yields0100hh->SetBinError(4, total0_100hhError);
    yields0100hh->SetMarkerStyle(29);
    yields0100hh->SetMarkerColor(kViolet-1);
    yields0100hh->SetLineColor(kViolet-1);
    yields0100hh->SetMarkerSize(2);
    yields0100hh->SetLineWidth(2);

    TLine *linehh = new TLine(3.0, 0.0, 3.0, 10.000);
    linehh->SetLineStyle(7);
    linehh->SetLineWidth(2);

    TCanvas *yieldhhc = new TCanvas("yieldhh", "yieldhh",50, 50, 600, 600);
    yieldhhc->cd();
    yields2050hh->GetYaxis()->SetRangeUser(0.0000, 10.000);
    yields2050hh->GetXaxis()->SetLabelSize(0.07);
    yields2050hh->Draw("P SAME");
    yields50100hh->GetXaxis()->SetLimits(0.05, 4.05);
    yields50100hh->Draw("P SAME");
    yields020hh->GetXaxis()->SetLimits(-0.05, 3.95);
    yields020hh->Draw("P SAME");
    yields0100hh->GetXaxis()->SetLimits(-0.1, 3.9);
    yields0100hh->Draw("P SAME");
    linehh->Draw("SAME");
    ratioslegend->Draw("SAME");




    TH1D *ratio = (TH1D*)hLambdadphi_0_20->Clone("ratio");
    ratio->Divide(hhdphi_0_20);


    TLegend  *legend = new TLegend(0.3791, 0.1518, 0.8772, 0.2688);
    legend->SetMargin(0.15);
    legend->AddEntry(corrFit2, "Hadron-#Lambda(1020) Correlation", "l");
    legend->AddEntry(corrFit, "Hadron-hadron Correlations", "l");

    TPaveText *text = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text->AddText("ALICE");
    text->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text->AddText("0%-20% Multiplicity");
    text->SetTextSizePixels(20);
    text->SetBorderSize(0);
    text->SetFillColor(kWhite);

    TPaveText *text2050 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text2050->AddText("ALICE");
    text2050->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text2050->AddText("20%-50% Multiplicity");
    text2050->SetTextSizePixels(20);
    text2050->SetBorderSize(0);
    text2050->SetFillColor(kWhite);

    TPaveText *text50100 = new TPaveText(0.4815, 0.7056, 0.8658, 0.8551, "NDC");
    text50100->AddText("ALICE");
    text50100->AddText("p-Pb #sqrt{s_{NN}} = 5 TeV");
    text50100->AddText("50%-80% Multiplicity");
    text50100->SetBorderSize(0);
    text50100->SetTextSizePixels(20);
    text50100->SetFillColor(kWhite);

    Double_t dphi020syst = 0.056;
    Double_t dphi2050syst = 0.056;
    Double_t dphi5080syst = 0.073;

    Double_t hhdphi020syst = 0.042;
    Double_t hhdphi2050syst = 0.042;
    Double_t hhdphi5080syst = 0.042;


    TH1D* hLambdadphi_0_20_syst = (TH1D*)hLambdadphi_0_20->Clone("hLambdadphi_0_20_syst");
    hLambdadphi_0_20_syst->SetFillColor(kGray);
    TH1D* hLambdadphi_20_50_syst = (TH1D*)hLambdadphi_20_50->Clone("hLambdadphi_20_50_syst");
    hLambdadphi_20_50_syst->SetFillColor(kGray);
    TH1D* hLambdadphi_50_100_syst = (TH1D*)hLambdadphi_50_100->Clone("hLambdadphi_50_100_syst");
    hLambdadphi_50_100_syst->SetFillColor(kGray);

    TH1D* hhdphi_0_20_syst = (TH1D*)hhdphi_0_20->Clone("hhdphi_0_20_syst");
    hhdphi_0_20_syst->SetFillColor(kGray);
    TH1D* hhdphi_20_50_syst = (TH1D*)hhdphi_20_50->Clone("hhdphi_20_50_syst");
    hhdphi_20_50_syst->SetFillColor(kGray);
    TH1D* hhdphi_50_100_syst = (TH1D*)hhdphi_50_100->Clone("hhdphi_50_100_syst");
    hhdphi_50_100_syst->SetFillColor(kGray);

    for(int i = 1; i <= 16; i++){
        hLambdadphi_0_20_syst->SetBinError(i, hLambdadphi_0_20_syst->GetBinContent(i)*dphi020syst);
        hLambdadphi_20_50_syst->SetBinError(i, hLambdadphi_20_50_syst->GetBinContent(i)*dphi2050syst);
        hLambdadphi_50_100_syst->SetBinError(i, hLambdadphi_50_100_syst->GetBinContent(i)*dphi5080syst);
        hhdphi_0_20_syst->SetBinError(i, hhdphi_0_20_syst->GetBinContent(i)*hhdphi020syst);
        hhdphi_20_50_syst->SetBinError(i, hhdphi_20_50_syst->GetBinContent(i)*hhdphi2050syst);
        hhdphi_50_100_syst->SetBinError(i, hhdphi_50_100_syst->GetBinContent(i)*hhdphi5080syst);

    }

    TCanvas *c0_20 = new TCanvas("c0_20", "c0_20", 50, 50, 550, 600);
    c0_20->cd();
    c0_20->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    //hhdphi_0_20->Draw("E0 X0");
    hLambdadphi_0_20->GetYaxis()->SetTitle("Per Trigger (h-#Lambda) Pairs");
    hLambdadphi_0_20_syst->Draw("E2 SAME");
    hLambdadphi_0_20->Draw("E0 X0 HIST SAME");
    //corrFit->Draw("SAME");
    //corrFit2->Draw("SAME");
    //hhBG->Draw("SAME");
    hLambdaBG->Draw("SAME");
    text->Draw();
    text2->Draw();
    //legend->Draw();



    TCanvas *c0_20pp = new TCanvas("c0_20pp", "c0_20pp", 50, 50, 550, 600);
    c0_20pp->cd();
    c0_20pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.001, 0.25);
    hhdphi_0_20->GetYaxis()->SetTitle("Per Trigger (h-h) Pairs");
    hhdphi_0_20_syst->Draw("E2 SAME");
    hhdphi_0_20->Draw("E0 X0 HIST SAME");
    //hLambdadphi_0_20->Draw("E0 X0 SAME");
    // corrFit->Draw("SAME");
    corrFit2->Draw("SAME");
    //hhBG->Draw("SAME");
    //hLambdaBG->Draw("SAME");
    text->Draw();
    text2->Draw();
    //legend->Draw();

    //"cartoon" plot showing different yield regions
    TH1D* regionHist = (TH1D*)hhdphi_0_20->Clone("regionHist");
    regionHist->SetLineWidth(3);
    regionHist->SetLineColor(kBlack);

    TH1D* nearHist = (TH1D*)regionHist->Clone("nearHist");
    nearHist->SetFillColor(kRed+2);
    TH1D* awayHist = (TH1D*)regionHist->Clone("awayHist");
    awayHist->SetFillColor(kBlue+1);
    TH1D* bulkHist = (TH1D*)regionHist->Clone("bulkHist");
    bulkHist->SetLineColor(kGreen+2);
    bulkHist->SetLineWidth(1);
    bulkHist->SetFillColor(kGreen+2);

    for(int ibin = 1; ibin <= nearHist->GetXaxis()->GetNbins(); ibin++){
        if(ibin <= nearHist->GetXaxis()->GetNbins()/2.0){
            awayHist->SetBinContent(ibin, 0.);
        }else{
            nearHist->SetBinContent(ibin, 0.);
        }
        bulkHist->SetBinContent(ibin, hhBG->GetParameter(0));
    }

    TLegend* regionLeg = new TLegend(0.5, 0.6, 0.8, 0.8);
    regionLeg->AddEntry(nearHist, "Near-side Jet", "f");
    regionLeg->AddEntry(awayHist, "Away-side Jet", "f");
    regionLeg->AddEntry(bulkHist, "Underlying Event", "f");
    TCanvas* cregions = new TCanvas("cregions", "cregions", 50, 50, 600, 600);
    cregions->cd();
    regionHist->GetYaxis()->SetRangeUser(0.0, 0.8);
    regionHist->Draw("HIST");
    nearHist->Draw("HIST SAME");
    awayHist->Draw("HIST SAME");
    bulkHist->Draw("HIST SAME");
    regionHist->Draw("HIST SAME");
    regionLeg->Draw();


    TCanvas *c20_50 = new TCanvas("c20_50", "c20_50", 50, 50, 550, 600);
    c20_50->cd();
    c20_50->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    //hhdphi_20_50->Draw("E0 X0");
    hLambdadphi_20_50_syst->Draw("E2 SAME");
    hLambdadphi_20_50->Draw("E0 X0 HIST SAME");
    //corrFit2050->Draw("SAME");
    //corrFit2_2050->Draw("SAME");
    //hhBG_20_50->Draw("SAME");
    hLambdaBG_20_50->Draw("SAME");
    text2050->Draw();
    text2->Draw();
    //legend->Draw();

    TCanvas *c20_50pp = new TCanvas("c20_50pp", "c20_50pp", 50, 50, 550, 600);
    c20_50pp->cd();
    c20_50pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_0_20->GetYaxis()->SetRangeUser(0.03, 0.11);
    hhdphi_20_50_syst->Draw("E2 SAME");
    hhdphi_20_50->Draw("E0 X0 HIST SAME");
    //hLambdadphi_20_50->Draw("E0 X0 SAME");
    //corrFit2050->Draw("SAME");
    //corrFit2_2050->Draw("SAME");
    hhBG_20_50->Draw("SAME");
    //hLambdaBG_20_50->Draw("SAME");
    text2050->Draw();
    text2->Draw();
    //legend->Draw();

    TCanvas *c50_100 = new TCanvas("c50_100", "c50_100", 50, 50, 550, 600);
    c50_100->cd();
    c50_100->SetMargin(0.12, 0.05, 0.1, 0.05);
    //hhdphi_50_100->Draw("E0 X0");
    hLambdadphi_50_100_syst->Draw("E2 SAME");
    hLambdadphi_50_100->Draw("E0 X0 HIST SAME");
    //corrFit50100->Draw("SAME");
    //corrFit2_50100->Draw("SAME");
    //hhBG_50_100->Draw("SAME");
    hLambdaBG_50_100->Draw("SAME");
    text50100->Draw();
    text2->Draw();
    //legend->Draw();

    TCanvas *c50_100pp = new TCanvas("c50_100pp", "c50_100pp", 50, 50, 550, 600);
    c50_100pp->cd();
    c50_100pp->SetMargin(0.12, 0.05, 0.1, 0.05);
    hhdphi_50_100_syst->Draw("E2 SAME");
    hhdphi_50_100->Draw("E0 X0 HIST SAME");
    //hLambdadphi_50_100->Draw("E0 X0 SAME");
    //corrFit50100->Draw("SAME");
    //corrFit2_50100->Draw("SAME");
    hhBG_50_100->Draw("SAME");
    //hLambdaBG_50_100->Draw("SAME");
    text50100->Draw();
    text2->Draw();
    //legend->Draw();

    //Set-up and Draw the hh and hlambda correlations on single canvases

    TH1D* hlambda020 = (TH1D*)hLambdadphi_0_20->Clone("hlambda020");
    hlambda020->SetMarkerColor(kRed);
    hlambda020->SetMarkerSize(1.3);
    hlambda020->SetMarkerStyle(33);
    hlambda020->SetLineColor(kRed);
    hlambda020->Add(hLambdaBG, -1.0);
    hlambda020->Scale(1.0/(2.4*hlambda020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hlambda2050 = (TH1D*)hLambdadphi_20_50->Clone("hlambda2050");
    hlambda2050->SetMarkerColor(kOrange);
    hlambda2050->SetLineColor(kOrange);
    hlambda2050->SetMarkerSize(1.3);
    hlambda2050->SetMarkerStyle(22);
    hlambda2050->Add(hLambdaBG_20_50, -1.0);
    hlambda2050->Scale(1.0/(2.4*hlambda2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hlambda5080 = (TH1D*)hLambdadphi_50_100->Clone("hlambda5080");
    hlambda5080->SetMarkerColor(kAzure);
    hlambda5080->SetMarkerSize(1.3);
    hlambda5080->SetMarkerStyle(23);
    hlambda5080->SetLineColor(kAzure);
    hlambda5080->Add(hLambdaBG_50_100, -1.0);
    hlambda5080->Scale(1.0/(2.4*hlambda5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hh020 = (TH1D*)hhdphi_0_20->Clone("hh020");
    hh020->SetMarkerColor(kRed);
    hh020->SetMarkerSize(1.3);
    hh020->SetMarkerStyle(34);
    hh020->SetLineColor(kRed);
    hh020->Add(hhBG, -1.0);
    hh020->Scale(1.0/(2.4*hh020->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hh2050 = (TH1D*)hhdphi_20_50->Clone("hh2050");
    hh2050->SetMarkerColor(kOrange);
    hh2050->SetMarkerSize(1.3);
    hh2050->SetMarkerStyle(21);
    hh2050->SetLineColor(kOrange);
    hh2050->Add(hhBG_20_50, -1.0);
    hh2050->Scale(1.0/(2.4*hh2050->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TH1D* hh5080 = (TH1D*)hhdphi_50_100->Clone("hh5080");
    hh5080->SetMarkerColor(kAzure);
    hh5080->SetMarkerSize(1.3);
    hh5080->SetMarkerStyle(20);
    hh5080->SetLineColor(kAzure);
    hh5080->Add(hhBG_50_100, -1.0);
    hh5080->Scale(1.0/(2.4*hh5080->GetXaxis()->GetBinWidth(1))); //scale by delta-eta and bin width

    TF1* fline = new TF1("fline", "pol0", -0.5*TMath::Pi(), 1.5*TMath::Pi());
    fline->SetParameter(0, 0.0);
    fline->SetLineColor(kBlack);
    fline->SetLineWidth(2);
    fline->SetLineStyle(7);

    TLegend* hlambdaleg = new TLegend(0.456, 0.641, 0.613, 0.937);
    hlambdaleg->SetBorderSize(0);
    hlambdaleg->AddEntry(hlambda020, "0-20%", "lep");
    hlambdaleg->AddEntry(hlambda2050, "20-50%", "lep");
    hlambdaleg->AddEntry(hlambda5080, "50-80%", "lep");

    TLegend* hhleg = new TLegend(0.456, 0.641, 0.613, 0.937);
    hhleg->SetBorderSize(0);
    hhleg->AddEntry(hh020, "0-20%", "lep");
    hhleg->AddEntry(hh2050, "20-50%", "lep");
    hhleg->AddEntry(hh5080, "50-80%", "lep");

    TPaveText* hlambdatext = new TPaveText(0.181, 0.827, 0.356, 0.925, "NDC");
    hlambdatext->SetBorderSize(0);
    hlambdatext->SetFillColor(kWhite);
    hlambdatext->AddText("h-#Lambda");

    TPaveText* hhtext = new TPaveText(0.704, 0.604, 0.861, 0.700, "NDC");
    hhtext->SetBorderSize(0);
    hhtext->SetFillColor(kWhite);
    hhtext->AddText("h-h");

    TCanvas* chlambda = new TCanvas("chlambda", "chlambda", 50, 50, 550, 600);
    chlambda->cd();
    chlambda->SetMargin(0.12, 0.05, 0.1, 0.05);
    hlambda020->GetYaxis()->SetTitle("1/N_{trig} dN/d#Delta#varphi per #Delta#eta - constant (rad^{-1})");
    hlambda020->Draw("E0 X0 P SAME");
    hlambda2050->Draw("E0 X0 P SAME");
    hlambda5080->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    hlambdatext->Draw();
    text2->Draw();
    hlambdaleg->Draw();

    TCanvas* chh = new TCanvas("chhh", "chhh", 50, 50, 550, 600);
    chh->cd();
    chh->SetMargin(0.12, 0.05, 0.1, 0.05);
    hh020->Draw("E0 X0 P SAME");
    hh2050->Draw("E0 X0 P SAME");
    hh5080->Draw("E0 X0 P SAME");
    fline->Draw("SAME");
    hhtext->Draw();
    text2->Draw();
    hhleg->Draw();





    /*
    TCanvas *cratio = new TCanvas("cratio", "cratio", 50, 50, 550, 600);
    cratio->cd();
    cratio->SetMargin(0.12, 0.05, 0.1, 0.05);
    ratio->Draw("E0 X0");
*/

    //Do plots of the yields for the systematics:
  /*TH1F* hhyieldSyst_0_20 = new TH1F("hhyieldSyst_0_20", "hhyieldSyst_0_20", 1000, 0.0, 10.0);
    TH1F* hhyieldSyst_20_50 = new TH1F("hhyieldSyst_20_50", "hhyieldSyst_20_50", 1000, 0.0, 10.0);
    TH1F* hhyieldSyst_50_100 = new TH1F("hhyieldSyst_50_100", "hhyieldSyst_50_100", 1000, 0.0, 10.0);
    TH1F* hlambdayieldSyst_0_20 = new TH1F("hlambdayieldSyst_0_20", "hlambdayieldSyst_0_20", 700, 3.0E-02, 10.0E-02);
    TH1F* hlambdayieldSyst_20_50 = new TH1F("hlambdayieldSyst_20_50", "hlambdayieldSyst_20_50", 600, 1.0E-02, 8.0E-02);
    TH1F* hlambdayieldSyst_50_100 = new TH1F("hlambdayieldSyst_50_100", "hlambdayieldSyst_50_100", 700, 0.5E-02, 7.5E-02);

    TH1F* hhnearyieldSyst_0_20 = new TH1F("hhnearyieldSyst_0_20", "hhnearyieldSyst_0_20", 1000, 0.0, 2.0);
    TH1F* hhnearyieldSyst_20_50 = new TH1F("hhnearyieldSyst_20_50", "hhnearyieldSyst_20_50", 1000, 0.0, 2.0);
    TH1F* hhnearyieldSyst_50_100 = new TH1F("hhnearyieldSyst_50_100", "hhnearyieldSyst_50_100", 1000, 0.0, 2.0);
    TH1F* hlambdanearyieldSyst_0_20 = new TH1F("hlambdanearyieldSyst_0_20", "hlambdanearyieldSyst_0_20", 800, 0.5E-03, 4.5E-03);
    TH1F* hlambdanearyieldSyst_20_50 = new TH1F("hlambdanearyieldSyst_20_50", "hlambdanearyieldSyst_20_50", 800, 0.5E-03, 4.5E-03);
    TH1F* hlambdanearyieldSyst_50_100 = new TH1F("hlambdanearyieldSyst_50_100", "hlambdanearyieldSyst_50_100", 800, 0.5E-03, 4.5E-03);

    TH1F* hhawayyieldSyst_0_20 = new TH1F("hhawayyieldSyst_0_20", "hhawayyieldSyst_0_20", 1000, 0.0, 2.0);
    TH1F* hhawayyieldSyst_20_50 = new TH1F("hhawayyieldSyst_20_50", "hhawayyieldSyst_20_50", 1000, 0.0, 2.0);
    TH1F* hhawayyieldSyst_50_100 = new TH1F("hhawayyieldSyst_50_100", "hhawayyieldSyst_50_100", 1000, 0.0, 2.0);
    TH1F* hlambdaawayyieldSyst_0_20 = new TH1F("hlambdaawayyieldSyst_0_20", "hlambdaawayyieldSyst_0_20", 800, 0.5E-03, 4.5E-03);
    TH1F* hlambdaawayyieldSyst_20_50 = new TH1F("hlambdaawayyieldSyst_20_50", "hlambdaawayyieldSyst_20_50", 800, 0.5E-03, 4.5E-03);
    TH1F* hlambdaawayyieldSyst_50_100 = new TH1F("hlambdaawayyieldSyst_50_100", "hlambdaawayyieldSyst_50_100", 800, 0.5E-03, 4.5E-03);

    TH1F* hhbulkyieldSyst_0_20 = new TH1F("hhbulkyieldSyst_0_20", "hhbulkyieldSyst_0_20", 1000, 0.0, 10.0);
    TH1F* hhbulkyieldSyst_20_50 = new TH1F("hhbulkyieldSyst_20_50", "hhbulkyieldSyst_20_50", 1000, 0.0, 10.0);
    TH1F* hhbulkyieldSyst_50_100 = new TH1F("hhbulkyieldSyst_50_100", "hhbulkyieldSyst_50_100", 1000, 0.0, 10.0);
    TH1F* hlambdabulkyieldSyst_0_20 = new TH1F("hlambdabulkyieldSyst_0_20", "hlambdabulkyieldSyst_0_20", 700, 3.0E-02, 10.0E-02);
    TH1F* hlambdabulkyieldSyst_20_50 = new TH1F("hlambdabulkyieldSyst_20_50", "hlambdabulkyieldSyst_20_50", 600, 1.0E-02, 8.0E-02);
    TH1F* hlambdabulkyieldSyst_50_100 = new TH1F("hlambdabulkyieldSyst_50_100", "hlambdabulkyieldSyst_50_100", 700, 0.5E-02, 7.5E-02);*/


    TH1F* hhyieldSyst_0_20 = new TH1F("hhyieldSyst_0_20", "hhyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhyieldSyst_20_50 = new TH1F("hhyieldSyst_20_50", "hhyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhyieldSyst_50_100 = new TH1F("hhyieldSyst_50_100", "hhyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hlambdayieldSyst_0_20 = new TH1F("hlambdayieldSyst_0_20", "hlambdayieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hlambdayieldSyst_20_50 = new TH1F("hlambdayieldSyst_20_50", "hlambdayieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hlambdayieldSyst_50_100 = new TH1F("hlambdayieldSyst_50_100", "hlambdayieldSyst_50_100", 500, 0.75, 1.25);

    TH1F* hhnearyieldSyst_0_20 = new TH1F("hhnearyieldSyst_0_20", "hhnearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhnearyieldSyst_20_50 = new TH1F("hhnearyieldSyst_20_50", "hhnearyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhnearyieldSyst_50_100 = new TH1F("hhnearyieldSyst_50_100", "hhnearyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hlambdanearyieldSyst_0_20 = new TH1F("hlambdanearyieldSyst_0_20", "hlambdanearyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hlambdanearyieldSyst_20_50 = new TH1F("hlambdanearyieldSyst_20_50", "hlambdanearyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hlambdanearyieldSyst_50_100 = new TH1F("hlambdanearyieldSyst_50_100", "hlambdanearyieldSyst_50_100", 500, 0.75, 1.25);

    TH1F* hhawayyieldSyst_0_20 = new TH1F("hhawayyieldSyst_0_20", "hhawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhawayyieldSyst_20_50 = new TH1F("hhawayyieldSyst_20_50", "hhawayyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhawayyieldSyst_50_100 = new TH1F("hhawayyieldSyst_50_100", "hhawayyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hlambdaawayyieldSyst_0_20 = new TH1F("hlambdaawayyieldSyst_0_20", "hlambdaawayyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hlambdaawayyieldSyst_20_50 = new TH1F("hlambdaawayyieldSyst_20_50", "hlambdaawayyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hlambdaawayyieldSyst_50_100 = new TH1F("hlambdaawayyieldSyst_50_100", "hlambdaawayyieldSyst_50_100", 500, 0.75, 1.25);

    TH1F* hhbulkyieldSyst_0_20 = new TH1F("hhbulkyieldSyst_0_20", "hhbulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hhbulkyieldSyst_20_50 = new TH1F("hhbulkyieldSyst_20_50", "hhbulkyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hhbulkyieldSyst_50_100 = new TH1F("hhbulkyieldSyst_50_100", "hhbulkyieldSyst_50_100", 500, 0.75, 1.25);
    TH1F* hlambdabulkyieldSyst_0_20 = new TH1F("hlambdabulkyieldSyst_0_20", "hlambdabulkyieldSyst_0_20", 500, 0.75, 1.25);
    TH1F* hlambdabulkyieldSyst_20_50 = new TH1F("hlambdabulkyieldSyst_20_50", "hlambdabulkyieldSyst_20_50", 500, 0.75, 1.25);
    TH1F* hlambdabulkyieldSyst_50_100 = new TH1F("hlambdabulkyieldSyst_50_100", "hlambdabulkyieldSyst_50_100", 500, 0.75, 1.25);

    Float_t stdtotal0_20 = 5.850084E-02;
    Float_t stdtotal20_50 = 3.489754E-02;
    Float_t stdtotal50_100 = 1.881542E-02;
    Float_t stdnear0_20 = 1.509740E-03;
    Float_t stdnear20_50 = 1.450373E-03;
    Float_t stdnear50_100 = 8.015203E-04;
    Float_t stdaway0_20 = 2.234562E-03;
    Float_t stdaway20_50 = 2.274811E-03;
    Float_t stdaway50_100 = 1.119789E-03;
    Float_t stdmid0_20 = 5.476692E-02;
    Float_t stdmid20_50 = 3.111197E-02;
    Float_t stdmid50_100 = 1.686431E-02;



    hhyieldSyst_0_20->Fill(total0_20hhYield);
    hhyieldSyst_20_50->Fill(total20_50hhYield);
    hhyieldSyst_50_100->Fill(total50_100hhYield);
    hhnearyieldSyst_0_20->Fill(near0_20hhYield);
    hhnearyieldSyst_20_50->Fill(near20_50hhYield);
    hhnearyieldSyst_50_100->Fill(near50_100hhYield);
    hhawayyieldSyst_0_20->Fill(away0_20hhYield);
    hhawayyieldSyst_20_50->Fill(away20_50hhYield);
    hhawayyieldSyst_50_100->Fill(away50_100hhYield);
    hhbulkyieldSyst_0_20->Fill(mid0_20hhYield);
    hhbulkyieldSyst_20_50->Fill(mid20_50hhYield);
    hhbulkyieldSyst_50_100->Fill(mid50_100hhYield);

    hlambdayieldSyst_0_20->Fill(total0_20hLambdaYield/stdtotal0_20);
    hlambdayieldSyst_20_50->Fill(total20_50hLambdaYield/stdtotal20_50);
    hlambdayieldSyst_50_100->Fill(total50_100hLambdaYield/stdtotal50_100);
    hlambdanearyieldSyst_0_20->Fill(near0_20hLambdaYield/stdnear0_20);
    hlambdanearyieldSyst_20_50->Fill(near20_50hLambdaYield/stdnear20_50);
    hlambdanearyieldSyst_50_100->Fill(near50_100hLambdaYield/stdnear50_100);
    hlambdaawayyieldSyst_0_20->Fill(away0_20hLambdaYield/stdaway0_20);
    hlambdaawayyieldSyst_20_50->Fill(away20_50hLambdaYield/stdaway20_50);
    hlambdaawayyieldSyst_50_100->Fill(away50_100hLambdaYield/stdaway50_100);
    hlambdabulkyieldSyst_0_20->Fill(mid0_20hLambdaYield/stdmid0_20);
    hlambdabulkyieldSyst_20_50->Fill(mid20_50hLambdaYield/stdmid20_50);
    hlambdabulkyieldSyst_50_100->Fill(mid50_100hLambdaYield/stdmid50_100);


/*
    TCanvas *csyst = new TCanvas("csyst", "csyst", 50, 50, 1000, 500);
    csyst->Divide(1,3);
    csyst->cd(1);
    hhyieldSyst_0_20->Draw();
    csyst->cd(2);
    hhyieldSyst_20_50->Draw();
    csyst->cd(3);
    hhyieldSyst_50_100->Draw();
*/
    TFile* output = new TFile(Form("fitsyst_%s.root", outputstring.Data()), "RECREATE");

    hhyieldSyst_0_20->Write();
    hhyieldSyst_20_50->Write();
    hhyieldSyst_50_100->Write();
    hhnearyieldSyst_0_20->Write();
    hhnearyieldSyst_20_50->Write();
    hhnearyieldSyst_50_100->Write();
    hhawayyieldSyst_0_20->Write();
    hhawayyieldSyst_20_50->Write();
    hhawayyieldSyst_50_100->Write();
    hhbulkyieldSyst_0_20->Write();
    hhbulkyieldSyst_20_50->Write();
    hhbulkyieldSyst_50_100->Write();

    hhdphi_0_20->Write();
    hLambdadphi_0_20->Write();
    hhdphi_20_50->Write();
    hLambdadphi_20_50->Write();
    hhdphi_50_100->Write();
    hLambdadphi_50_100->Write();


    hlambdayieldSyst_0_20->Write();
    hlambdayieldSyst_20_50->Write();
    hlambdayieldSyst_50_100->Write();
    hlambdanearyieldSyst_0_20->Write();
    hlambdanearyieldSyst_20_50->Write();
    hlambdanearyieldSyst_50_100->Write();
    hlambdaawayyieldSyst_0_20->Write();
    hlambdaawayyieldSyst_20_50->Write();
    hlambdaawayyieldSyst_50_100->Write();
    hlambdabulkyieldSyst_0_20->Write();
    hlambdabulkyieldSyst_20_50->Write();
    hlambdabulkyieldSyst_50_100->Write();

    ratioNearHist->Write();
    ratiosNear->Write("ratiosNear");
    ratiosAway->Write("ratiosAway");
    ratiosBulk->Write("ratiosBulk");
    ratiosTot->Write("ratiosTot");


 }


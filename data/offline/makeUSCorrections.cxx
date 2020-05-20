void makeUSCorrections(string inputFile){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hLambda2Dpeak = (TH2D*)input->Get("hLambda2Dpeak");
    TH2D* hLambda2DLside = (TH2D*)input->Get("hLambda2DLside");
    TH2D* hLambda2DRside = (TH2D*)input->Get("hLambda2DRside");
    TH2D* hLambdaLS2Dpeak = (TH2D*)input->Get("hLambdaLS2Dpeak");
    TH2D* hLambdaLS2DLside = (TH2D*)input->Get("hLambdaLS2DLside");
    TH2D* hLambdaLS2DRside = (TH2D*)input->Get("hLambdaLS2DRside");

    TH2D* trigDistSameUS = (TH2D*)input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = (TH2D*)input->Get("fTrigSameLSDist");

    hLambda2Dpeak->SetName("uncorrectedhLambda2Dpeak");

    TH2D* hLambdaBGPeakRegionL = (TH2D*)hLambda2DLside->Clone("hLambdaBGPeakRegionL");
    hLambdaBGPeakRegionL->Scale(1.0/(hLambda2DLside->Integral(hLambda2DLside->GetXaxis()->FindBin(-1.2), hLambda2DLside->GetXaxis()->FindBin(1.2), 1, hLambda2DLside->GetYaxis()->GetNbins())));
    hLambdaBGPeakRegionL->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hLambdaBGPeakRegionL_deta = (TH1D*)hLambdaBGPeakRegionL->ProjectionX("hLambdaBGPeakRegionL_deta", 1, hLambdaBGPeakRegionL->GetYaxis()->GetNbins());
    TH1D* hLambdaBGPeakRegionL_dphi = (TH1D*)hLambdaBGPeakRegionL->ProjectionY("hLambdaBGPeakRegionL_dphi", hLambdaBGPeakRegionL->GetXaxis()->FindBin(-1.2), hLambdaBGPeakRegionL->GetXaxis()->FindBin(1.2));

    TH2D* hLambdaBGPeakRegionR = (TH2D*)hLambda2DRside->Clone("hLambdaBGPeakRegionR");
    hLambdaBGPeakRegionR->Scale(1.0/(hLambda2DRside->Integral(hLambda2DRside->GetXaxis()->FindBin(-1.2), hLambda2DRside->GetXaxis()->FindBin(1.2), 1, hLambda2DRside->GetYaxis()->GetNbins())));
    hLambdaBGPeakRegionR->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hLambdaBGPeakRegionR_deta = (TH1D*)hLambdaBGPeakRegionR->ProjectionX("hLambdaBGPeakRegionR_deta", 1, hLambdaBGPeakRegionR->GetYaxis()->GetNbins());
    TH1D* hLambdaBGPeakRegionR_dphi = (TH1D*)hLambdaBGPeakRegionR->ProjectionY("hLambdaBGPeakRegionR_dphi", hLambdaBGPeakRegionR->GetXaxis()->FindBin(-1.2), hLambdaBGPeakRegionR->GetXaxis()->FindBin(1.2));

    TH2D* hLambdaBGPeakRegion = (TH2D*)hLambdaBGPeakRegionL->Clone("hLambdaBGPeakregion");
    hLambdaBGPeakRegion->Add(hLambdaBGPeakRegionR);
    hLambdaBGPeakRegion->Scale(0.5);
    hLambdaBGPeakRegion->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hLambdaBGPeakRegion_deta = (TH1D*)hLambdaBGPeakRegion->ProjectionX("hLambdaBGPeakRegion_deta", 1, hLambdaBGPeakRegion->GetYaxis()->GetNbins());
    TH1D* hLambdaBGPeakRegion_dphi = (TH1D*)hLambdaBGPeakRegion->ProjectionY("hLambdaBGPeakRegion_dphi", hLambdaBGPeakRegion->GetXaxis()->FindBin(-1.2), hLambdaBGPeakRegion->GetXaxis()->FindBin(1.2));


    //US residual checks between SB average and the Left and Right separately
    TH2D* resLeftVsAvg = (TH2D*)hLambdaBGPeakRegionL->Clone("resLeftVsAvg");
    resLeftVsAvg->Add(hLambdaBGPeakRegion, -1.0);
    resLeftVsAvg->Divide(hLambdaBGPeakRegionL);
    resLeftVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_deta = (TH1D*)hLambdaBGPeakRegionL_deta->Clone("resLeftVsAvg_deta");
    resLeftVsAvg_deta->Add(hLambdaBGPeakRegion_deta, -1.0);
    resLeftVsAvg_deta->Divide(hLambdaBGPeakRegionL_deta);
    resLeftVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_dphi = (TH1D*)hLambdaBGPeakRegionL_dphi->Clone("resLeftVsAvg_dphi");
    resLeftVsAvg_dphi->Add(hLambdaBGPeakRegion_dphi, -1.0);
    resLeftVsAvg_dphi->Divide(hLambdaBGPeakRegionL_dphi);

    TH2D* resRightVsAvg = (TH2D*)hLambdaBGPeakRegionR->Clone("resRightVsAbg");
    resRightVsAvg->Add(hLambdaBGPeakRegion, -1.0);
    resRightVsAvg->Divide(hLambdaBGPeakRegionR);
    resRightVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_deta = (TH1D*)hLambdaBGPeakRegionR_deta->Clone("resRightVsAvg_deta");
    resRightVsAvg_deta->Add(hLambdaBGPeakRegion_deta, -1.0);
    resRightVsAvg_deta->Divide(hLambdaBGPeakRegionR_deta);
    resRightVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_dphi = (TH1D*)hLambdaBGPeakRegionR_dphi->Clone("resRightVsAvg_dphi");
    resRightVsAvg_dphi->Add(hLambdaBGPeakRegion_dphi, -1.0);
    resRightVsAvg_dphi->Divide(hLambdaBGPeakRegionR_dphi);



    Float_t leftscale = hLambda2DLside->Integral(hLambda2DLside->GetXaxis()->FindBin(-1.2), hLambda2DLside->GetXaxis()->FindBin(1.2), 1, hLambda2DLside->GetYaxis()->GetNbins())/hLambdaLS2DLside->Integral(hLambda2DLside->GetXaxis()->FindBin(-1.2), hLambda2DLside->GetXaxis()->FindBin(1.2), 1, hLambda2DLside->GetYaxis()->GetNbins());

    TH2D* LLSsubhLambda2DLside = (TH2D*)hLambda2DLside->Clone("LLSsubhLambda2DLside");
    TH2D* LLSsubhLambda2Dpeak = (TH2D*)hLambda2Dpeak->Clone("LLSsubhLambda2Dpeak");
    LLSsubhLambda2DLside->Add(hLambdaLS2DLside, -1.0*leftscale);
    //LLSsubhLambda2DLside->Divide(hLambdaLS2DLside);
    //LLSsubhLambda2DLside->Scale(1.0/leftscale);
    LLSsubhLambda2Dpeak->Add(hLambdaLS2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhLambda2DLside_deta = LLSsubhLambda2DLside->ProjectionX("LLSsubhLambda2DLside_deta", 1, LLSsubhLambda2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhLambda2DLside_dphi = LLSsubhLambda2DLside->ProjectionY("LLSsubhLambda2DLside_dphi", LLSsubhLambda2DLside->GetXaxis()->FindBin(-1.2), LLSsubhLambda2DLside->GetXaxis()->FindBin(1.2));

    //Float_t rightscale = hLambda2DRside->Integral(1, hLambda2DRside->GetXaxis()->GetNbins(), 1, hLambda2DRside->GetYaxis()->GetNbins())/hLambdaLS2DRside->Integral(1, hLambda2DRside->GetXaxis()->GetNbins(), 1, hLambda2DRside->GetYaxis()->GetNbins());
    gStyle->SetOptStat(0);
    Float_t rightscale = hLambda2DRside->Integral(hLambda2DRside->GetXaxis()->FindBin(-1.2), hLambda2DRside->GetXaxis()->FindBin(1.2), 1, hLambda2DRside->GetYaxis()->GetNbins())/hLambdaLS2DRside->Integral(hLambda2DRside->GetXaxis()->FindBin(-1.2), hLambda2DRside->GetXaxis()->FindBin(1.2), 1, hLambda2DRside->GetYaxis()->GetNbins());
    TH2D* RLSsubhLambda2DRside = (TH2D*)hLambda2DRside->Clone("RLSsubhLambda2DRside");
    TH2D* RLSsubhLambda2Dpeak = (TH2D*)hLambda2Dpeak->Clone("RLSsubhLambda2Dpeak");
    TH2D* RLSsubhLambda2DRsideScaled = (TH2D*)hLambda2DRside->Clone("RLSsubhLambda2DRsideScaled");
    RLSsubhLambda2DRsideScaled->GetXaxis()->SetRangeUser(-1.2, 1.2);
    RLSsubhLambda2DRsideScaled->Scale(rightscale);
    TH1D* RLSsubhLambdaDPhiRsideScaled = RLSsubhLambda2DRsideScaled->ProjectionY();
    RLSsubhLambdaDPhiRsideScaled->SetTitle("h-#Lambda^{0} scaled R-sideband #Delta#varphi distribution");
    RLSsubhLambdaDPhiRsideScaled->SetLineColor(6);
    RLSsubhLambdaDPhiRsideScaled->SetLineWidth(3);
    RLSsubhLambdaDPhiRsideScaled->SetMarkerColor(6);
    TCanvas *presCanvas = new TCanvas("presCanvas", "Presentation Canvas", 0, 10, 1600, 1200);
    presCanvas->cd();
    RLSsubhLambdaDPhiRsideScaled->Draw();

    RLSsubhLambda2DRside->Add(hLambdaLS2DRside, -1.0*rightscale);
    //RLSsubhLambda2DRside->Divide(hLambdaLS2DRside);
    //RLSsubhLambda2DRside->Scale(1.0/rightscale);
    RLSsubhLambda2Dpeak->Add(hLambdaLS2Dpeak, -1.0*rightscale);
    RLSsubhLambda2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);

    TH1D* RLSsubhLambda2DRside_deta = RLSsubhLambda2DRside->ProjectionX("RLSsubhLambda2DRside_deta", 1, RLSsubhLambda2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhLambda2DRside_dphi = RLSsubhLambda2DRside->ProjectionY("RLSsubhLambda2DRside_dphi", RLSsubhLambda2DRside->GetXaxis()->FindBin(-1.2), RLSsubhLambda2DRside->GetXaxis()->FindBin(1.2));

    TH1D* RLSsubhLambda2Dpeak_deta = RLSsubhLambda2Dpeak->ProjectionX("RLSsubhLambda2Dpeak_deta", 1, RLSsubhLambda2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhLambda2Dpeak_dphi = RLSsubhLambda2Dpeak->ProjectionY("RLSsubhLambda2Dpeak_dphi", RLSsubhLambda2Dpeak->GetXaxis()->FindBin(-1.2), RLSsubhLambda2Dpeak->GetXaxis()->FindBin(1.2));


    TH1D* scales = new TH1D("scales", "scales", 2, -1, 1);
    scales->SetBinContent(1, leftscale);
    scales->SetBinContent(2, rightscale);

    TH2D* rebinRLSsubhLambda2Dpeak = (TH2D*)RLSsubhLambda2Dpeak->Clone("rebinRLSsubhLambda2Dpeak");
    rebinRLSsubhLambda2Dpeak->Rebin2D(2, 2);

    //Using US estimate for BG to subtract off the from the peak region:

    Float_t scaleUS = (rightscale)*hLambdaLS2Dpeak->Integral(hLambdaLS2Dpeak->GetXaxis()->FindBin(-1.2), hLambdaLS2Dpeak->GetXaxis()->FindBin(1.2), 1, hLambdaLS2Dpeak->GetYaxis()->GetNbins());
    Float_t scaletest = (rightscale)*hLambda2Dpeak->Integral(hLambda2Dpeak->GetXaxis()->FindBin(-1.2), hLambda2Dpeak->GetXaxis()->FindBin(1.2), 1, hLambda2Dpeak->GetYaxis()->GetNbins());


    printf("\n\nscaleUS = %e\n\ntestscale = %e \n\n", scaleUS, scaletest);

    //avg of right and left US sideband tests
    TH2D* AvgUSsubhLambda2Dpeak = (TH2D*)hLambda2Dpeak->Clone("AvgUSsubhLambda2Dpeak");
    AvgUSsubhLambda2Dpeak->Add(hLambdaBGPeakRegion, -1.0*scaleUS);
    AvgUSsubhLambda2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhLambda2Dpeak_deta = (TH1D*)AvgUSsubhLambda2Dpeak->ProjectionX("AvgUSsubhLambda2Dpeak_deta", 1, AvgUSsubhLambda2Dpeak->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhLambda2Dpeak_dphi = (TH1D*)AvgUSsubhLambda2Dpeak->ProjectionY("AvgUSsubhLambda2Dpeak_dphi", AvgUSsubhLambda2Dpeak->GetXaxis()->FindBin(-1.2), AvgUSsubhLambda2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhLambda2Dpeakleftscale = (TH2D*)hLambda2Dpeak->Clone("AvgUSsubhLambda2Dpeakleftscale");
    AvgUSsubhLambda2Dpeakleftscale->Add(hLambdaBGPeakRegion, -1.0*scaleUS*leftscale/rightscale);
    AvgUSsubhLambda2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhLambda2Dpeakleftscale_deta = (TH1D*)AvgUSsubhLambda2Dpeakleftscale->ProjectionX("AvgUSsubhLambda2Dpeakleftscale_deta", 1, AvgUSsubhLambda2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhLambda2Dpeakleftscale_dphi = (TH1D*)AvgUSsubhLambda2Dpeakleftscale->ProjectionY("AvgUSsubhLambda2Dpeakleftscale_dphi", AvgUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(-1.2), AvgUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhLambda2Dpeakavgscale = (TH2D*)hLambda2Dpeak->Clone("AvgUSsubhLambda2Dpeakavgscale");
    AvgUSsubhLambda2Dpeakavgscale->Add(hLambdaBGPeakRegion, -1.0*scaleUS*(leftscale + rightscale)/(2.0*rightscale));
    AvgUSsubhLambda2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhLambda2Dpeakavgscale_deta = (TH1D*)AvgUSsubhLambda2Dpeakavgscale->ProjectionX("AvgUSsubhLambda2Dpeakavgscale_deta", 1, AvgUSsubhLambda2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhLambda2Dpeakavgscale_dphi = (TH1D*)AvgUSsubhLambda2Dpeakavgscale->ProjectionY("AvgUSsubhLambda2Dpeakavgscale_dphi", AvgUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(-1.2), AvgUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //right side US sideband tests
    TH2D* RSUSsubhLambda2Dpeak = (TH2D*)hLambda2Dpeak->Clone("RSUSsubhLambda2Dpeak");
    RSUSsubhLambda2Dpeak->Add(hLambdaBGPeakRegionR, -1.0*scaleUS);
    RSUSsubhLambda2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhLambda2Dpeak_deta = (TH1D*)RSUSsubhLambda2Dpeak->ProjectionX("RSUSsubhLambda2Dpeak_deta", 1, RSUSsubhLambda2Dpeak->GetYaxis()->GetNbins());
    TH1D* RSUSsubhLambda2Dpeak_dphi = (TH1D*)RSUSsubhLambda2Dpeak->ProjectionY("RSUSsubhLambda2Dpeak_dphi", RSUSsubhLambda2Dpeak->GetXaxis()->FindBin(-1.2), RSUSsubhLambda2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* RSUSsubhLambda2Dpeakleftscale = (TH2D*)hLambda2Dpeak->Clone("RSUSsubhLambda2Dpeakleftscale");
    RSUSsubhLambda2Dpeakleftscale->Add(hLambdaBGPeakRegionR, -1.0*scaleUS*leftscale/rightscale);
    RSUSsubhLambda2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhLambda2Dpeakleftscale_deta = (TH1D*)RSUSsubhLambda2Dpeakleftscale->ProjectionX("RSUSsubhLambda2Dpeakleftscale_deta", 1, RSUSsubhLambda2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhLambda2Dpeakleftscale_dphi = (TH1D*)RSUSsubhLambda2Dpeakleftscale->ProjectionY("RSUSsubhLambda2Dpeakleftscale_dphi", RSUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(-1.2), RSUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    // //GOTO HERE
    // TH2D* RSUSsubhLambda2Dpeakavgscale = (TH2D*)hLambda2Dpeak->Clone("RSUSsubhLambda2Dpeakavgscale");
    // RSUSsubhLambda2Dpeakavgscale->Add(hLambdaBGPeakRegionR, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    // RSUSsubhLambda2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    // TH1D* RSUSsubhLambda2Dpeakavgscale_deta = (TH1D*)RSUSsubhLambda2Dpeakavgscale->ProjectionX("RSUSsubhLambda2Dpeakavgscale_deta", 1, RSUSsubhLambda2Dpeakavgscale->GetYaxis()->GetNbins());
    // TH1D* RSUSsubhLambda2Dpeakavgscale_dphi = (TH1D*)RSUSsubhLambda2Dpeakavgscale->ProjectionY("RSUSsubhLambda2Dpeakavgscale_dphi", RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(1.2));
    // //END GOTO

    double signalOverTotalArray[3] = {0.3620, 0.4584, 0.5214}; // have to change this for each multiplicity bin
    double scaleFactorArray[3] = {0.731, 0.706, 0.735}; // have to change this for each multiplicity bin and for each method...

    TH2D* RSUSsubhLambda2Dpeakavgscale = (TH2D*)hLambda2Dpeak->Clone("RSUSsubhLambda2Dpeakavgscale");

    // Using BG = TOTAL - SIGNAL = TOTAL(1-S/TOTAL)
    double bgIntegral = RSUSsubhLambda2Dpeakavgscale->Integral(RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(1.2), 1, RSUSsubhLambda2Dpeakavgscale->GetYaxis()->GetNbins())*(1 - signalOverTotalArray[2]);

    RSUSsubhLambda2Dpeakavgscale->Add(hLambdaBGPeakRegionR, -1.0*bgIntegral);
    RSUSsubhLambda2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhLambda2Dpeakavgscale_deta = (TH1D*)RSUSsubhLambda2Dpeakavgscale->ProjectionX("RSUSsubhLambda2Dpeakavgscale_deta", 1, RSUSsubhLambda2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhLambda2Dpeakavgscale_dphi = (TH1D*)RSUSsubhLambda2Dpeakavgscale->ProjectionY("RSUSsubhLambda2Dpeakavgscale_dphi", RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //left side US sideband tests
    TH2D* LSUSsubhLambda2Dpeak = (TH2D*)hLambda2Dpeak->Clone("LSUSsubhLambda2Dpeak");
    LSUSsubhLambda2Dpeak->Add(hLambdaBGPeakRegionL, -1.0*scaleUS);
    LSUSsubhLambda2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhLambda2Dpeak_deta = (TH1D*)LSUSsubhLambda2Dpeak->ProjectionX("LSUSsubhLambda2Dpeak_deta", 1, LSUSsubhLambda2Dpeak->GetYaxis()->GetNbins());
    TH1D* LSUSsubhLambda2Dpeak_dphi = (TH1D*)LSUSsubhLambda2Dpeak->ProjectionY("LSUSsubhLambda2Dpeak_dphi", LSUSsubhLambda2Dpeak->GetXaxis()->FindBin(-1.2), LSUSsubhLambda2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhLambda2Dpeakleftscale = (TH2D*)hLambda2Dpeak->Clone("LSUSsubhLambda2Dpeakleftscale");
    LSUSsubhLambda2Dpeakleftscale->Add(hLambdaBGPeakRegionL, -1.0*scaleUS*leftscale/rightscale);
    LSUSsubhLambda2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhLambda2Dpeakleftscale_deta = (TH1D*)LSUSsubhLambda2Dpeakleftscale->ProjectionX("LSUSsubhLambda2Dpeakleftscale_deta", 1, LSUSsubhLambda2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhLambda2Dpeakleftscale_dphi = (TH1D*)LSUSsubhLambda2Dpeakleftscale->ProjectionY("LSUSsubhLambda2Dpeakleftscale_dphi", LSUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(-1.2), LSUSsubhLambda2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhLambda2Dpeakavgscale = (TH2D*)hLambda2Dpeak->Clone("LSUSsubhLambda2Dpeakavgscale");
    LSUSsubhLambda2Dpeakavgscale->Add(hLambdaBGPeakRegionL, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    LSUSsubhLambda2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhLambda2Dpeakavgscale_deta = (TH1D*)LSUSsubhLambda2Dpeakavgscale->ProjectionX("LSUSsubhLambda2Dpeakavgscale_deta", 1, LSUSsubhLambda2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhLambda2Dpeakavgscale_dphi = (TH1D*)LSUSsubhLambda2Dpeakavgscale->ProjectionY("LSUSsubhLambda2Dpeakavgscale_dphi", LSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(-1.2), LSUSsubhLambda2Dpeakavgscale->GetXaxis()->FindBin(1.2));


    TH2D* resUSvsLS = (TH2D*)AvgUSsubhLambda2Dpeak->Clone("resUSvsLS");
    resUSvsLS->Add(RLSsubhLambda2Dpeak, -1.0);
    resUSvsLS->Divide(AvgUSsubhLambda2Dpeak);
    resUSvsLS->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_deta = (TH1D*)AvgUSsubhLambda2Dpeak_deta->Clone("resUSvsLS_deta");
    resUSvsLS_deta->Add(RLSsubhLambda2Dpeak_deta, -1.0);
    resUSvsLS_deta->Divide(AvgUSsubhLambda2Dpeak_deta);
    resUSvsLS_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_dphi = (TH1D*)AvgUSsubhLambda2Dpeak_dphi->Clone("resUSvsLS_dphi");
    resUSvsLS_dphi->Add(RLSsubhLambda2Dpeak_dphi, -1.0);
    resUSvsLS_dphi->Divide(AvgUSsubhLambda2Dpeak_dphi);



    TFile* output = new TFile(Form("US_syst_%s", inputFile.c_str()), "RECREATE");
    LLSsubhLambda2DLside->Write();
    LLSsubhLambda2DLside_deta->Write();
    LLSsubhLambda2DLside_dphi->Write();
    LLSsubhLambda2Dpeak->Write();
    RLSsubhLambda2DRside->Write();
    RLSsubhLambda2DRside_deta->Write();
    RLSsubhLambda2DRside_dphi->Write();
    RLSsubhLambda2Dpeak->Write();
    RLSsubhLambda2Dpeak_deta->Write();
    RLSsubhLambda2Dpeak_dphi->Write();
    rebinRLSsubhLambda2Dpeak->Write();
    scales->Write();
    hLambda2Dpeak->Write();
    hLambdaBGPeakRegionL->Write();
    hLambdaBGPeakRegionL_deta->Write();
    hLambdaBGPeakRegionL_dphi->Write();
    hLambdaBGPeakRegionR->Write();
    hLambdaBGPeakRegionR_deta->Write();
    hLambdaBGPeakRegionR_dphi->Write();
    hLambdaBGPeakRegion->Write();
    hLambdaBGPeakRegion_deta->Write();
    hLambdaBGPeakRegion_dphi->Write();
    resLeftVsAvg->Write();
    resLeftVsAvg_deta->Write();
    resLeftVsAvg_dphi->Write();
    resRightVsAvg->Write();
    resRightVsAvg_deta->Write();
    resRightVsAvg_dphi->Write();
    AvgUSsubhLambda2Dpeak->Write();
    AvgUSsubhLambda2Dpeak_deta->Write();
    AvgUSsubhLambda2Dpeak_dphi->Write();
    AvgUSsubhLambda2Dpeakleftscale->Write();
    AvgUSsubhLambda2Dpeakleftscale_deta->Write();
    AvgUSsubhLambda2Dpeakleftscale_dphi->Write();
    AvgUSsubhLambda2Dpeakavgscale->Write();
    AvgUSsubhLambda2Dpeakavgscale_deta->Write();
    AvgUSsubhLambda2Dpeakavgscale_dphi->Write();
    RSUSsubhLambda2Dpeak->Write();
    RSUSsubhLambda2Dpeak_deta->Write();
    RSUSsubhLambda2Dpeak_dphi->Write();
    RSUSsubhLambda2Dpeakleftscale->Write();
    RSUSsubhLambda2Dpeakleftscale_deta->Write();
    RSUSsubhLambda2Dpeakleftscale_dphi->Write();
    RSUSsubhLambda2Dpeakavgscale->Write();
    RSUSsubhLambda2Dpeakavgscale_deta->Write();
    RSUSsubhLambda2Dpeakavgscale_dphi->Write();
    LSUSsubhLambda2Dpeak->Write();
    LSUSsubhLambda2Dpeak_deta->Write();
    LSUSsubhLambda2Dpeak_dphi->Write();
    LSUSsubhLambda2Dpeakleftscale->Write();
    LSUSsubhLambda2Dpeakleftscale_deta->Write();
    LSUSsubhLambda2Dpeakleftscale_dphi->Write();
    LSUSsubhLambda2Dpeakavgscale->Write();
    LSUSsubhLambda2Dpeakavgscale_deta->Write();
    LSUSsubhLambda2Dpeakavgscale_dphi->Write();
    resUSvsLS->Write();
    resUSvsLS_deta->Write();
    resUSvsLS_dphi->Write();
}



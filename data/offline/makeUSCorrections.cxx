void makeUSCorrections(string inputFile){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hStrangePart2Dpeak = (TH2D*)input->Get("hStrangePart2Dpeak");
    TH2D* hStrangePart2DLside = (TH2D*)input->Get("hStrangePart2DLside");
    TH2D* hStrangePart2DRside = (TH2D*)input->Get("hStrangePart2DRside");
    TH2D* hStrangePartLS2Dpeak = (TH2D*)input->Get("hStrangePartLS2Dpeak");
    TH2D* hStrangePartLS2DLside = (TH2D*)input->Get("hStrangePartLS2DLside");
    TH2D* hStrangePartLS2DRside = (TH2D*)input->Get("hStrangePartLS2DRside");

    TH2D* trigDistSameUS = (TH2D*)input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = (TH2D*)input->Get("fTrigSameLSDist");

    hStrangePart2Dpeak->SetName("uncorrectedhStrangePart2Dpeak");

    TH2D* hStrangePartBGPeakRegionL = (TH2D*)hStrangePart2DLside->Clone("hStrangePartBGPeakRegionL");
    hStrangePartBGPeakRegionL->Scale(1.0/(hStrangePart2DLside->Integral(hStrangePart2DLside->GetXaxis()->FindBin(-1.2), hStrangePart2DLside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DLside->GetYaxis()->GetNbins())));
    hStrangePartBGPeakRegionL->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hStrangePartBGPeakRegionL_deta = (TH1D*)hStrangePartBGPeakRegionL->ProjectionX("hStrangePartBGPeakRegionL_deta", 1, hStrangePartBGPeakRegionL->GetYaxis()->GetNbins());
    TH1D* hStrangePartBGPeakRegionL_dphi = (TH1D*)hStrangePartBGPeakRegionL->ProjectionY("hStrangePartBGPeakRegionL_dphi", hStrangePartBGPeakRegionL->GetXaxis()->FindBin(-1.2), hStrangePartBGPeakRegionL->GetXaxis()->FindBin(1.2));

    TH2D* hStrangePartBGPeakRegionR = (TH2D*)hStrangePart2DRside->Clone("hStrangePartBGPeakRegionR");
    hStrangePartBGPeakRegionR->Scale(1.0/(hStrangePart2DRside->Integral(hStrangePart2DRside->GetXaxis()->FindBin(-1.2), hStrangePart2DRside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DRside->GetYaxis()->GetNbins())));
    hStrangePartBGPeakRegionR->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hStrangePartBGPeakRegionR_deta = (TH1D*)hStrangePartBGPeakRegionR->ProjectionX("hStrangePartBGPeakRegionR_deta", 1, hStrangePartBGPeakRegionR->GetYaxis()->GetNbins());
    TH1D* hStrangePartBGPeakRegionR_dphi = (TH1D*)hStrangePartBGPeakRegionR->ProjectionY("hStrangePartBGPeakRegionR_dphi", hStrangePartBGPeakRegionR->GetXaxis()->FindBin(-1.2), hStrangePartBGPeakRegionR->GetXaxis()->FindBin(1.2));

    TH2D* hStrangePartBGPeakRegion = (TH2D*)hStrangePartBGPeakRegionL->Clone("hStrangePartBGPeakregion");
    hStrangePartBGPeakRegion->Add(hStrangePartBGPeakRegionR);
    hStrangePartBGPeakRegion->Scale(0.5);
    hStrangePartBGPeakRegion->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hStrangePartBGPeakRegion_deta = (TH1D*)hStrangePartBGPeakRegion->ProjectionX("hStrangePartBGPeakRegion_deta", 1, hStrangePartBGPeakRegion->GetYaxis()->GetNbins());
    TH1D* hStrangePartBGPeakRegion_dphi = (TH1D*)hStrangePartBGPeakRegion->ProjectionY("hStrangePartBGPeakRegion_dphi", hStrangePartBGPeakRegion->GetXaxis()->FindBin(-1.2), hStrangePartBGPeakRegion->GetXaxis()->FindBin(1.2));


    //US residual checks between SB average and the Left and Right separately
    TH2D* resLeftVsAvg = (TH2D*)hStrangePartBGPeakRegionL->Clone("resLeftVsAbg");
    resLeftVsAvg->Add(hStrangePartBGPeakRegion, -1.0);
    resLeftVsAvg->Divide(hStrangePartBGPeakRegionL);
    resLeftVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_deta = (TH1D*)hStrangePartBGPeakRegionL_deta->Clone("resLeftVsAvg_deta");
    resLeftVsAvg_deta->Add(hStrangePartBGPeakRegion_deta, -1.0);
    resLeftVsAvg_deta->Divide(hStrangePartBGPeakRegionL_deta);
    resLeftVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_dphi = (TH1D*)hStrangePartBGPeakRegionL_dphi->Clone("resLeftVsAvg_dphi");
    resLeftVsAvg_dphi->Add(hStrangePartBGPeakRegion_dphi, -1.0);
    resLeftVsAvg_dphi->Divide(hStrangePartBGPeakRegionL_dphi);

    TH2D* resRightVsAvg = (TH2D*)hStrangePartBGPeakRegionR->Clone("resRightVsAbg");
    resRightVsAvg->Add(hStrangePartBGPeakRegion, -1.0);
    resRightVsAvg->Divide(hStrangePartBGPeakRegionR);
    resRightVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_deta = (TH1D*)hStrangePartBGPeakRegionR_deta->Clone("resRightVsAvg_deta");
    resRightVsAvg_deta->Add(hStrangePartBGPeakRegion_deta, -1.0);
    resRightVsAvg_deta->Divide(hStrangePartBGPeakRegionR_deta);
    resRightVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_dphi = (TH1D*)hStrangePartBGPeakRegionR_dphi->Clone("resRightVsAvg_dphi");
    resRightVsAvg_dphi->Add(hStrangePartBGPeakRegion_dphi, -1.0);
    resRightVsAvg_dphi->Divide(hStrangePartBGPeakRegionR_dphi);



    Float_t leftscale = hStrangePart2DLside->Integral(hStrangePart2DLside->GetXaxis()->FindBin(-1.2), hStrangePart2DLside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DLside->GetYaxis()->GetNbins())/hStrangePartLS2DLside->Integral(hStrangePart2DLside->GetXaxis()->FindBin(-1.2), hStrangePart2DLside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DLside->GetYaxis()->GetNbins());

    TH2D* LLSsubhStrangePart2DLside = (TH2D*)hStrangePart2DLside->Clone("LLSsubhStrangePart2DLside");
    TH2D* LLSsubhStrangePart2Dpeak = (TH2D*)hStrangePart2Dpeak->Clone("LLSsubhStrangePart2Dpeak");
    LLSsubhStrangePart2DLside->Add(hStrangePartLS2DLside, -1.0*leftscale);
    //LLSsubhStrangePart2DLside->Divide(hStrangePartLS2DLside);
    //LLSsubhStrangePart2DLside->Scale(1.0/leftscale);
    LLSsubhStrangePart2Dpeak->Add(hStrangePartLS2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhStrangePart2DLside_deta = LLSsubhStrangePart2DLside->ProjectionX("LLSsubhStrangePart2DLside_deta", 1, LLSsubhStrangePart2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhStrangePart2DLside_dphi = LLSsubhStrangePart2DLside->ProjectionY("LLSsubhStrangePart2DLside_dphi", LLSsubhStrangePart2DLside->GetXaxis()->FindBin(-1.2), LLSsubhStrangePart2DLside->GetXaxis()->FindBin(1.2));

    //Float_t rightscale = hStrangePart2DRside->Integral(1, hStrangePart2DRside->GetXaxis()->GetNbins(), 1, hStrangePart2DRside->GetYaxis()->GetNbins())/hStrangePartLS2DRside->Integral(1, hStrangePart2DRside->GetXaxis()->GetNbins(), 1, hStrangePart2DRside->GetYaxis()->GetNbins());
    Float_t rightscale = hStrangePart2DRside->Integral(hStrangePart2DRside->GetXaxis()->FindBin(-1.2), hStrangePart2DRside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DRside->GetYaxis()->GetNbins())/hStrangePartLS2DRside->Integral(hStrangePart2DRside->GetXaxis()->FindBin(-1.2), hStrangePart2DRside->GetXaxis()->FindBin(1.2), 1, hStrangePart2DRside->GetYaxis()->GetNbins());
    TH2D* RLSsubhStrangePart2DRside = (TH2D*)hStrangePart2DRside->Clone("RLSsubhStrangePart2DRside");
    TH2D* RLSsubhStrangePart2Dpeak = (TH2D*)hStrangePart2Dpeak->Clone("RLSsubhStrangePart2Dpeak");
    RLSsubhStrangePart2DRside->Add(hStrangePartLS2DRside, -1.0*rightscale);
    //RLSsubhStrangePart2DRside->Divide(hStrangePartLS2DRside);
    //RLSsubhStrangePart2DRside->Scale(1.0/rightscale);
    RLSsubhStrangePart2Dpeak->Add(hStrangePartLS2Dpeak, -1.0*rightscale);
    RLSsubhStrangePart2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);

    TH1D* RLSsubhStrangePart2DRside_deta = RLSsubhStrangePart2DRside->ProjectionX("RLSsubhStrangePart2DRside_deta", 1, RLSsubhStrangePart2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhStrangePart2DRside_dphi = RLSsubhStrangePart2DRside->ProjectionY("RLSsubhStrangePart2DRside_dphi", RLSsubhStrangePart2DRside->GetXaxis()->FindBin(-1.2), RLSsubhStrangePart2DRside->GetXaxis()->FindBin(1.2));

    TH1D* RLSsubhStrangePart2Dpeak_deta = RLSsubhStrangePart2Dpeak->ProjectionX("RLSsubhStrangePart2Dpeak_deta", 1, RLSsubhStrangePart2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhStrangePart2Dpeak_dphi = RLSsubhStrangePart2Dpeak->ProjectionY("RLSsubhStrangePart2Dpeak_dphi", RLSsubhStrangePart2Dpeak->GetXaxis()->FindBin(-1.2), RLSsubhStrangePart2Dpeak->GetXaxis()->FindBin(1.2));


    TH1D* scales = new TH1D("scales", "scales", 2, -1, 1);
    scales->SetBinContent(1, leftscale);
    scales->SetBinContent(2, rightscale);

    TH2D* rebinRLSsubhStrangePart2Dpeak = (TH2D*)RLSsubhStrangePart2Dpeak->Clone("rebinRLSsubhStrangePart2Dpeak");
    rebinRLSsubhStrangePart2Dpeak->Rebin2D(2, 2);

    //Using US estimate for BG to subtract off the from the peak region:

    Float_t scaleUS = (rightscale)*hStrangePartLS2Dpeak->Integral(hStrangePartLS2Dpeak->GetXaxis()->FindBin(-1.2), hStrangePartLS2Dpeak->GetXaxis()->FindBin(1.2), 1, hStrangePartLS2Dpeak->GetYaxis()->GetNbins());
    Float_t scaletest = (rightscale)*hStrangePart2Dpeak->Integral(hStrangePart2Dpeak->GetXaxis()->FindBin(-1.2), hStrangePart2Dpeak->GetXaxis()->FindBin(1.2), 1, hStrangePart2Dpeak->GetYaxis()->GetNbins());


    printf("\n\nscaleUS = %e\n\ntestscale = %e \n\n", scaleUS, scaletest);

    //avg of right and left US sideband tests
    TH2D* AvgUSsubhStrangePart2Dpeak = (TH2D*)hStrangePart2Dpeak->Clone("AvgUSsubhStrangePart2Dpeak");
    AvgUSsubhStrangePart2Dpeak->Add(hStrangePartBGPeakRegion, -1.0*scaleUS);
    AvgUSsubhStrangePart2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhStrangePart2Dpeak_deta = (TH1D*)AvgUSsubhStrangePart2Dpeak->ProjectionX("AvgUSsubhStrangePart2Dpeak_deta", 1, AvgUSsubhStrangePart2Dpeak->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhStrangePart2Dpeak_dphi = (TH1D*)AvgUSsubhStrangePart2Dpeak->ProjectionY("AvgUSsubhStrangePart2Dpeak_dphi", AvgUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(-1.2), AvgUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhStrangePart2Dpeakleftscale = (TH2D*)hStrangePart2Dpeak->Clone("AvgUSsubhStrangePart2Dpeakleftscale");
    AvgUSsubhStrangePart2Dpeakleftscale->Add(hStrangePartBGPeakRegion, -1.0*scaleUS*leftscale/rightscale);
    AvgUSsubhStrangePart2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhStrangePart2Dpeakleftscale_deta = (TH1D*)AvgUSsubhStrangePart2Dpeakleftscale->ProjectionX("AvgUSsubhStrangePart2Dpeakleftscale_deta", 1, AvgUSsubhStrangePart2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhStrangePart2Dpeakleftscale_dphi = (TH1D*)AvgUSsubhStrangePart2Dpeakleftscale->ProjectionY("AvgUSsubhStrangePart2Dpeakleftscale_dphi", AvgUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(-1.2), AvgUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhStrangePart2Dpeakavgscale = (TH2D*)hStrangePart2Dpeak->Clone("AvgUSsubhStrangePart2Dpeakavgscale");
    AvgUSsubhStrangePart2Dpeakavgscale->Add(hStrangePartBGPeakRegion, -1.0*scaleUS*(leftscale + rightscale)/(2.0*rightscale));
    AvgUSsubhStrangePart2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhStrangePart2Dpeakavgscale_deta = (TH1D*)AvgUSsubhStrangePart2Dpeakavgscale->ProjectionX("AvgUSsubhStrangePart2Dpeakavgscale_deta", 1, AvgUSsubhStrangePart2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhStrangePart2Dpeakavgscale_dphi = (TH1D*)AvgUSsubhStrangePart2Dpeakavgscale->ProjectionY("AvgUSsubhStrangePart2Dpeakavgscale_dphi", AvgUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(-1.2), AvgUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //right side US sideband tests
    TH2D* RSUSsubhStrangePart2Dpeak = (TH2D*)hStrangePart2Dpeak->Clone("RSUSsubhStrangePart2Dpeak");
    RSUSsubhStrangePart2Dpeak->Add(hStrangePartBGPeakRegionR, -1.0*scaleUS);
    RSUSsubhStrangePart2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhStrangePart2Dpeak_deta = (TH1D*)RSUSsubhStrangePart2Dpeak->ProjectionX("RSUSsubhStrangePart2Dpeak_deta", 1, RSUSsubhStrangePart2Dpeak->GetYaxis()->GetNbins());
    TH1D* RSUSsubhStrangePart2Dpeak_dphi = (TH1D*)RSUSsubhStrangePart2Dpeak->ProjectionY("RSUSsubhStrangePart2Dpeak_dphi", RSUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(-1.2), RSUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* RSUSsubhStrangePart2Dpeakleftscale = (TH2D*)hStrangePart2Dpeak->Clone("RSUSsubhStrangePart2Dpeakleftscale");
    RSUSsubhStrangePart2Dpeakleftscale->Add(hStrangePartBGPeakRegionR, -1.0*scaleUS*leftscale/rightscale);
    RSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhStrangePart2Dpeakleftscale_deta = (TH1D*)RSUSsubhStrangePart2Dpeakleftscale->ProjectionX("RSUSsubhStrangePart2Dpeakleftscale_deta", 1, RSUSsubhStrangePart2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhStrangePart2Dpeakleftscale_dphi = (TH1D*)RSUSsubhStrangePart2Dpeakleftscale->ProjectionY("RSUSsubhStrangePart2Dpeakleftscale_dphi", RSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(-1.2), RSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* RSUSsubhStrangePart2Dpeakavgscale = (TH2D*)hStrangePart2Dpeak->Clone("RSUSsubhStrangePart2Dpeakavgscale");
    RSUSsubhStrangePart2Dpeakavgscale->Add(hStrangePartBGPeakRegionR, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    RSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhStrangePart2Dpeakavgscale_deta = (TH1D*)RSUSsubhStrangePart2Dpeakavgscale->ProjectionX("RSUSsubhStrangePart2Dpeakavgscale_deta", 1, RSUSsubhStrangePart2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhStrangePart2Dpeakavgscale_dphi = (TH1D*)RSUSsubhStrangePart2Dpeakavgscale->ProjectionY("RSUSsubhStrangePart2Dpeakavgscale_dphi", RSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //left side US sideband tests
    TH2D* LSUSsubhStrangePart2Dpeak = (TH2D*)hStrangePart2Dpeak->Clone("LSUSsubhStrangePart2Dpeak");
    LSUSsubhStrangePart2Dpeak->Add(hStrangePartBGPeakRegionL, -1.0*scaleUS);
    LSUSsubhStrangePart2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhStrangePart2Dpeak_deta = (TH1D*)LSUSsubhStrangePart2Dpeak->ProjectionX("LSUSsubhStrangePart2Dpeak_deta", 1, LSUSsubhStrangePart2Dpeak->GetYaxis()->GetNbins());
    TH1D* LSUSsubhStrangePart2Dpeak_dphi = (TH1D*)LSUSsubhStrangePart2Dpeak->ProjectionY("LSUSsubhStrangePart2Dpeak_dphi", LSUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(-1.2), LSUSsubhStrangePart2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhStrangePart2Dpeakleftscale = (TH2D*)hStrangePart2Dpeak->Clone("LSUSsubhStrangePart2Dpeakleftscale");
    LSUSsubhStrangePart2Dpeakleftscale->Add(hStrangePartBGPeakRegionL, -1.0*scaleUS*leftscale/rightscale);
    LSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhStrangePart2Dpeakleftscale_deta = (TH1D*)LSUSsubhStrangePart2Dpeakleftscale->ProjectionX("LSUSsubhStrangePart2Dpeakleftscale_deta", 1, LSUSsubhStrangePart2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhStrangePart2Dpeakleftscale_dphi = (TH1D*)LSUSsubhStrangePart2Dpeakleftscale->ProjectionY("LSUSsubhStrangePart2Dpeakleftscale_dphi", LSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(-1.2), LSUSsubhStrangePart2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhStrangePart2Dpeakavgscale = (TH2D*)hStrangePart2Dpeak->Clone("LSUSsubhStrangePart2Dpeakavgscale");
    LSUSsubhStrangePart2Dpeakavgscale->Add(hStrangePartBGPeakRegionL, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    LSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhStrangePart2Dpeakavgscale_deta = (TH1D*)LSUSsubhStrangePart2Dpeakavgscale->ProjectionX("LSUSsubhStrangePart2Dpeakavgscale_deta", 1, LSUSsubhStrangePart2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhStrangePart2Dpeakavgscale_dphi = (TH1D*)LSUSsubhStrangePart2Dpeakavgscale->ProjectionY("LSUSsubhStrangePart2Dpeakavgscale_dphi", LSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(-1.2), LSUSsubhStrangePart2Dpeakavgscale->GetXaxis()->FindBin(1.2));


    TH2D* resUSvsLS = (TH2D*)AvgUSsubhStrangePart2Dpeak->Clone("resUSvsLS");
    resUSvsLS->Add(RLSsubhStrangePart2Dpeak, -1.0);
    resUSvsLS->Divide(AvgUSsubhStrangePart2Dpeak);
    resUSvsLS->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_deta = (TH1D*)AvgUSsubhStrangePart2Dpeak_deta->Clone("resUSvsLS_deta");
    resUSvsLS_deta->Add(RLSsubhStrangePart2Dpeak_deta, -1.0);
    resUSvsLS_deta->Divide(AvgUSsubhStrangePart2Dpeak_deta);
    resUSvsLS_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_dphi = (TH1D*)AvgUSsubhStrangePart2Dpeak_dphi->Clone("resUSvsLS_dphi");
    resUSvsLS_dphi->Add(RLSsubhStrangePart2Dpeak_dphi, -1.0);
    resUSvsLS_dphi->Divide(AvgUSsubhStrangePart2Dpeak_dphi);



    // TFile* output = new TFile(Form("US_syst_%s", inputFile.c_str()), "RECREATE");
    TFile* output = new TFile("hello.root", "RECREATE");
    LLSsubhStrangePart2DLside->Write();
    LLSsubhStrangePart2DLside_deta->Write();
    LLSsubhStrangePart2DLside_dphi->Write();
    LLSsubhStrangePart2Dpeak->Write();
    RLSsubhStrangePart2DRside->Write();
    RLSsubhStrangePart2DRside_deta->Write();
    RLSsubhStrangePart2DRside_dphi->Write();
    RLSsubhStrangePart2Dpeak->Write();
    RLSsubhStrangePart2Dpeak_deta->Write();
    RLSsubhStrangePart2Dpeak_dphi->Write();
    rebinRLSsubhStrangePart2Dpeak->Write();
    scales->Write();
    hStrangePart2Dpeak->Write();
    hStrangePartBGPeakRegionL->Write();
    hStrangePartBGPeakRegionL_deta->Write();
    hStrangePartBGPeakRegionL_dphi->Write();
    hStrangePartBGPeakRegionR->Write();
    hStrangePartBGPeakRegionR_deta->Write();
    hStrangePartBGPeakRegionR_dphi->Write();
    hStrangePartBGPeakRegion->Write();
    hStrangePartBGPeakRegion_deta->Write();
    hStrangePartBGPeakRegion_dphi->Write();
    resLeftVsAvg->Write();
    resLeftVsAvg_deta->Write();
    resLeftVsAvg_dphi->Write();
    resRightVsAvg->Write();
    resRightVsAvg_deta->Write();
    resRightVsAvg_dphi->Write();
    AvgUSsubhStrangePart2Dpeak->Write();
    AvgUSsubhStrangePart2Dpeak_deta->Write();
    AvgUSsubhStrangePart2Dpeak_dphi->Write();
    AvgUSsubhStrangePart2Dpeakleftscale->Write();
    AvgUSsubhStrangePart2Dpeakleftscale_deta->Write();
    AvgUSsubhStrangePart2Dpeakleftscale_dphi->Write();
    AvgUSsubhStrangePart2Dpeakavgscale->Write();
    AvgUSsubhStrangePart2Dpeakavgscale_deta->Write();
    AvgUSsubhStrangePart2Dpeakavgscale_dphi->Write();
    RSUSsubhStrangePart2Dpeak->Write();
    RSUSsubhStrangePart2Dpeak_deta->Write();
    RSUSsubhStrangePart2Dpeak_dphi->Write();
    RSUSsubhStrangePart2Dpeakleftscale->Write();
    RSUSsubhStrangePart2Dpeakleftscale_deta->Write();
    RSUSsubhStrangePart2Dpeakleftscale_dphi->Write();
    RSUSsubhStrangePart2Dpeakavgscale->Write();
    RSUSsubhStrangePart2Dpeakavgscale_deta->Write();
    RSUSsubhStrangePart2Dpeakavgscale_dphi->Write();
    LSUSsubhStrangePart2Dpeak->Write();
    LSUSsubhStrangePart2Dpeak_deta->Write();
    LSUSsubhStrangePart2Dpeak_dphi->Write();
    LSUSsubhStrangePart2Dpeakleftscale->Write();
    LSUSsubhStrangePart2Dpeakleftscale_deta->Write();
    LSUSsubhStrangePart2Dpeakleftscale_dphi->Write();
    LSUSsubhStrangePart2Dpeakavgscale->Write();
    LSUSsubhStrangePart2Dpeakavgscale_deta->Write();
    LSUSsubhStrangePart2Dpeakavgscale_dphi->Write();
    resUSvsLS->Write();
    resUSvsLS_deta->Write();
    resUSvsLS_dphi->Write();
}



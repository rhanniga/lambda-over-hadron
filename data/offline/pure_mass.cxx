void pure_mass(){
    gStyle->SetOptStat(0);

    TFile* file = new TFile("../online/output/lambda_0_20.root");
    TList* list = (TList*) file->Get("h-lambda");
    THnSparseF* lambdaUSDist= (THnSparseF*)list->FindObject("fDphiHLambda");
    THnSparseF* lambdaLSDist= (THnSparseF*)list->FindObject("fDphiHLambdaLS");

    lambdaUSDist->GetAxis(0)->SetRangeUser(4, 8);
    lambdaLSDist->GetAxis(0)->SetRangeUser(4, 8);
    lambdaUSDist->GetAxis(1)->SetRangeUser(2.0, 4.0);
    lambdaLSDist->GetAxis(1)->SetRangeUser(2.0, 4.0);

    TH1D* USInvMass = lambdaUSDist->Projection(4);
    TH1D* LSInvMass = lambdaLSDist->Projection(4);

    USInvMass->Sumw2();
    LSInvMass->Sumw2();
    USInvMass->Rebin(3);
    LSInvMass->Rebin(3);

    float sidebandUS = (float)USInvMass->Integral(USInvMass->GetXaxis()->FindBin(1.132), USInvMass->GetXaxis()->FindBin(1.14));
    float sidebandLS = (float)LSInvMass->Integral(LSInvMass->GetXaxis()->FindBin(1.132), LSInvMass->GetXaxis()->FindBin(1.14));
    float scale = sidebandUS/sidebandLS;

    LSInvMass->Scale(scale);

    int left_bin = LSInvMass->FindBin(1.13);
    int right_bin = LSInvMass->FindBin(1.14);

    double signal = USInvMass->Integral(left_bin, right_bin) - LSInvMass->Integral(left_bin, right_bin);
    double bg = LSInvMass->Integral(left_bin, right_bin);

    std::cout << "signal: " << signal << " bg: " << bg << " sig/bg: " << signal/bg << std::endl;

    LSInvMass->SetLineColor(kMagenta);
    LSInvMass->SetLineWidth(1);
    USInvMass->SetLineWidth(1);

    USInvMass->SetTitle("");
    USInvMass->GetXaxis()->SetTitle("m_{p-Pi} (GeV/c^{2})");
    USInvMass->GetXaxis()->SetTitleSize(0.05);
    USInvMass->SetLineColor(kBlue);

    TH1D* corrected = (TH1D*)USInvMass->Clone("corrected");
    corrected->Add(LSInvMass, -1.0);

    TF1* fit = new TF1("fit", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)",1.09, 1.14);
    fit->SetParameter(0, 1.31792e+02);
    fit->SetParLimits(0, 1.31792e+02*0.9, 1.31792e+02*1.1);
	fit->SetParameter(1, 1.11584e+00);
    fit->SetParLimits(1, 1.11584e+00*0.9, 1.11584e+00*1.1);
	fit->SetParameter(2, 5.81684e-04);
    fit->SetParLimits(2, 5.81684e-04*0.9, 5.81684e-04*1.1);
    fit->SetParameter(3, 0.006);
    fit->SetParLimits(3, 0.006*0.9, 0.006*1.1);
    fit->SetParameter(4, 0);
    fit->SetParameter(5, 0);
	fit->SetLineColor(kBlue + 2);
	fit->SetLineStyle(7);
	fit->SetLineWidth(3);
	fit->SetNpx(500);
	corrected->Fit(fit, "R");


    int ndf = fit->GetNDF();
    float chi_square = fit->GetChisquare();

    std::cout << ndf << " is the NDF \n" << chi_square << " is the chi ssquare \n" << chi_square/ndf << " is chi_square/ndf" << std::endl;

	TF1* voigtFit = new TF1("voigtFit", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.09, 1.14);
	voigtFit->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3));

	TF1* bgFit = new TF1("bgFit", "pol1(0)", 1.09, 1.4);
	bgFit->SetParameters(fit->GetParameter(4), fit->GetParameter(5));

    corrected->SetMarkerSize(1);
    corrected->SetMarkerStyle(20);
    corrected->SetMarkerColor(kGreen + 3);
    corrected->SetLineColor(kGreen + 3);
    corrected->GetXaxis()->SetRangeUser(1.09, 1.138);

    bgFit->SetLineColor(kMagenta + 1);

    TCanvas *c = new TCanvas("c", "c", 50, 50, 1920, 1080);
    c->cd();
    USInvMass->Draw();
    LSInvMass->Draw("SAME");

    TLegend *corrleg = new TLegend(0.46, 0.39, 0.88, 0.57);
    corrleg->AddEntry(corrected, "Corrected US Inv. Mass", "p");
    corrleg->AddEntry(fit, "Inv. Mass Fit", "l");
    corrleg->AddEntry(bgFit, "Residual BG Fit", "l");

    TCanvas *c2 = new TCanvas("c2", "c2", 50, 50, 1920, 1080);
    c2->cd();
    corrected->Draw();
    bgFit->Draw("SAME");
    corrleg->Draw();

}

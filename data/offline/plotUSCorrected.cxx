void plotUSCorrected(string inputname){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile* eta20File = new TFile(inputname.c_str());

    TH2D* eta20peak = (TH2D*)eta20File->Get("AvgUSsubhKStar2Dpeak");
    TH2D* eta20RSB = (TH2D*)eta20File->Get("RSUSsubhKStar2Dpeak");
    TH2D* eta20LSB = (TH2D*)eta20File->Get("LSUSsubhKStar2Dpeak");


    eta20peak->GetXaxis()->SetTitle("#Delta#eta");
    eta20peak->GetXaxis()->SetTitleSize(0.05);
    eta20peak->GetXaxis()->SetTitleOffset(1.3);
    eta20peak->GetYaxis()->SetTitle("#Delta#varphi");
    eta20peak->GetYaxis()->SetTitleSize(0.05);
    eta20peak->GetYaxis()->SetTitleOffset(1.3);
    eta20peak->SetTitle("");
    //eta20peak->SetStats(kFALSE);
    //eta20peak->Scale(1.0/(eta20peak->Integral(eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2), 1, eta20peak->GetYaxis()->GetNbins())));

    eta20RSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20RSB->GetXaxis()->SetTitleSize(0.05);
    eta20RSB->GetXaxis()->SetTitleOffset(1.3);
    eta20RSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20RSB->GetYaxis()->SetTitleSize(0.05);
    eta20RSB->GetYaxis()->SetTitleOffset(1.3);
    eta20RSB->SetTitle("");
    eta20RSB->SetStats(kFALSE);

    eta20LSB->GetXaxis()->SetTitle("#Delta#eta");
    eta20LSB->GetXaxis()->SetTitleSize(0.05);
    eta20LSB->GetXaxis()->SetTitleOffset(1.3);
    eta20LSB->GetYaxis()->SetTitle("#Delta#varphi");
    eta20LSB->GetYaxis()->SetTitleSize(0.05);
    eta20LSB->GetYaxis()->SetTitleOffset(1.3);
    eta20LSB->SetTitle("");
    eta20LSB->SetStats(kFALSE);

    TH1D* eta20peakEta = eta20peak->ProjectionX("eta20peakEta", 1, eta20peak->GetYaxis()->GetNbins());
    eta20peakEta->GetXaxis()->SetTitleOffset(1.0);
    eta20peakEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20peakEta->SetStats(kFALSE);
    TH1D* eta20peakPhi = eta20peak->ProjectionY("eta20peakPhi", eta20peak->GetXaxis()->FindBin(-1.5), eta20peak->GetXaxis()->FindBin(1.5));
    eta20peakPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20peakPhi->SetStats(kFALSE);
    eta20peakPhi->SetLineColor(kBlue);
    TH1D* eta20peakPhiNarrow = eta20peak->ProjectionY("eta20peakPhiNarrow", eta20peak->GetXaxis()->FindBin(-1.2), eta20peak->GetXaxis()->FindBin(1.2));
    eta20peakPhiNarrow->SetLineColor(kViolet);
    eta20peakPhiNarrow->SetStats(kFALSE);
    eta20peakPhiNarrow->GetXaxis()->SetTitleOffset(1.0);
    TH1D* eta20peakPhiNarrowest = eta20peak->ProjectionY("eta20peakPhiNarrowest", eta20peak->GetXaxis()->FindBin(-1.0), eta20peak->GetXaxis()->FindBin(1.0));
    eta20peakPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20RSBEta = eta20RSB->ProjectionX("eta20RSBEta", 1, eta20RSB->GetYaxis()->GetNbins());
    eta20RSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSBEta->SetStats(kFALSE);
    TH1D* eta20RSBPhi = eta20RSB->ProjectionY("eta20RSBPhi", eta20RSB->GetXaxis()->FindBin(-1.5), eta20RSB->GetXaxis()->FindBin(1.5));
    eta20RSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20RSBPhi->SetLineColor(kBlue);
    eta20RSBPhi->SetStats(kFALSE);
    TH1D* eta20RSBPhiNarrow = eta20RSB->ProjectionY("eta20RSBPhiNarrow", eta20RSB->GetXaxis()->FindBin(-1.2), eta20RSB->GetXaxis()->FindBin(1.2));
    eta20RSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20RSBPhiNarrowest = eta20RSB->ProjectionY("eta20RSBPhiNarrowest", eta20RSB->GetXaxis()->FindBin(-1.0), eta20RSB->GetXaxis()->FindBin(1.0));
    eta20RSBPhiNarrowest->SetLineColor(kRed);

    TH1D* eta20LSBEta = eta20LSB->ProjectionX("eta20LSBEta", 1, eta20LSB->GetYaxis()->GetNbins());
    eta20LSBEta->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBEta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSBEta->SetStats(kFALSE);
    TH1D* eta20LSBPhi = eta20LSB->ProjectionY("eta20LSBPhi", eta20LSB->GetXaxis()->FindBin(-1.5), eta20LSB->GetXaxis()->FindBin(1.5));
    eta20LSBPhi->GetXaxis()->SetTitleOffset(1.0);
    eta20LSBPhi->SetLineColor(kBlue);
    eta20LSBPhi->SetStats(kFALSE);
    TH1D* eta20LSBPhiNarrow = eta20LSB->ProjectionY("eta20LSBPhiNarrow", eta20LSB->GetXaxis()->FindBin(-1.2), eta20LSB->GetXaxis()->FindBin(1.2));
    eta20LSBPhiNarrow->SetLineColor(kViolet);
    TH1D* eta20LSBPhiNarrowest = eta20LSB->ProjectionY("eta20LSBPhiNarrowest", eta20LSB->GetXaxis()->FindBin(-1.0), eta20LSB->GetXaxis()->FindBin(1.0));
    eta20LSBPhiNarrowest->SetLineColor(kRed);


    //reset eta range to narrow view for 2D plotting and rebin
//    eta20peak->Rebin2D(2,2);
//    eta20RSB->Rebin2D(2,2);
//    eta20LSB->Rebin2D(2,2);
    eta20peak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20RSB->GetXaxis()->SetRangeUser(-1.2, 1.2);
    eta20LSB->GetXaxis()->SetRangeUser(-1.2, 1.2);

    //set up legend for delta-phi projection
    TLegend* legend = new TLegend(0.5612, 0.6716, 0.89, 0.8712);
    legend->SetMargin(0.15);
    legend->AddEntry(eta20peakPhi, "-1.5 < #Delta#eta < 1.5", "le");
    legend->AddEntry(eta20peakPhiNarrow, "-1.2 < #Delta#eta < 1.2", "le");
    legend->AddEntry(eta20peakPhiNarrowest, "-1.0 < #Delta#eta < 1.0", "le");
    legend->SetBorderSize(0);

    TCanvas* ceta20peak = new TCanvas("ceta20peak", "ceta20peak", 50, 50, 800, 800);
    ceta20peak->Divide(2,2);
    ceta20peak->cd(1)->SetTheta(50);
    ceta20peak->cd(1)->SetPhi(50);
    eta20peak->Draw("SURF1");
    ceta20peak->cd(2);
    eta20peakEta->Draw("H");
    ceta20peak->cd(3);
    eta20peakPhiNarrow->GetYaxis()->SetRangeUser(0.5*eta20peakPhiNarrowest->GetMinimum(), 1.2*eta20peakPhiNarrow->GetMaximum());
    eta20peakPhi->Draw("H");
    eta20peakPhiNarrow->Draw("H SAME");
    eta20peakPhiNarrowest->Draw("H SAME");
    legend->Draw();

    TCanvas* ceta20RSB = new TCanvas("ceta20RSB", "ceta20RSB", 60, 60, 800, 800);
    ceta20RSB->Divide(2,2);
    ceta20RSB->cd(1)->SetTheta(50);
    ceta20RSB->cd(1)->SetPhi(50);
    eta20RSB->Draw("SURF1");
    ceta20RSB->cd(2);
    eta20RSBEta->Draw("H");
    ceta20RSB->cd(3);
    //eta20RSBPhi->GetYaxis()->SetRangeUser(1.5*eta20RSBPhiNarrowest->GetMinimum(), 1.5*eta20RSBPhi->GetMaximum());
    eta20RSBPhi->Draw("H");
    eta20RSBPhiNarrow->Draw("H SAME");
    eta20RSBPhiNarrowest->Draw("H SAME");
    legend->Draw();

    TCanvas* ceta20LSB = new TCanvas("ceta20LSB", "ceta20LSB", 70, 70, 800, 800);
    ceta20LSB->Divide(2,2);
    ceta20LSB->cd(1)->SetTheta(50);
    ceta20LSB->cd(1)->SetPhi(50);
    eta20LSB->Draw("SURF1");
    ceta20LSB->cd(2);
    eta20LSBEta->Draw("H");
    ceta20LSB->cd(3);
    //eta20LSBPhi->GetYaxis()->SetRangeUser(1.5*eta20LSBPhiNarrowest->GetMinimum(), 1.5* eta20LSBPhi->GetMaximum());
    eta20LSBPhi->Draw("H");
    eta20LSBPhiNarrow->Draw("H SAME");
    eta20LSBPhiNarrowest->Draw("H SAME");
    legend->Draw();
}

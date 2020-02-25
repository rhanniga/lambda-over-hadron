void raw_mass() {

    gStyle->SetOptStat(0); //because fuck statistics

    TFile *histo_file = new TFile("lambda.root");
    TList* list = (TList*) histo_file->Get("h-lambda");

    THnSparse *k_star_cor = (THnSparse*)list->FindObject("fDphiHKStar");
    k_star_cor->Sumw2();


    THnSparse *ls_k_star_cor = (THnSparse*)list->FindObject("fDphiHKStarLS");
    ls_k_star_cor->Sumw2();


    //Some raw mass plots with different pt ranges:

    TCanvas *raw_mass = new TCanvas("raw_mass", "Raw K* Invariant Mass At Different Pt", 0, 10, 1920, 1080);
    raw_mass->Divide(2, 2);

    TCanvas *raw_mass_2_4 = new TCanvas("raw_mass_2_4", "Raw K* Invariant Mass at 2 GeV to 4 GeV", 0, 10, 1920, 1080);

    int pt_ranges[5] = {2, 4, 6, 8, 12};
    int colors[4] = {4, 1, 6, 2};

    for(int i = 0; i < 4; i++) {

        gStyle->SetHistLineColor(colors[i]);
        gStyle->SetHistFillColor(colors[i]);

        raw_mass->cd(i+1);
        k_star_cor->GetAxis(1)->SetRangeUser(pt_ranges[i], pt_ranges[i+1]);
        auto* mass_proj = k_star_cor->Projection(4);

        mass_proj->SetName(Form("mass_proj_%d_%d", pt_ranges[i], pt_ranges[i+1]));
        mass_proj->SetTitle(Form("%d GeV < p_{t} < %d GeV", pt_ranges[i], pt_ranges[i+1]));

        mass_proj->GetXaxis()->SetTitle("Mass (GeV)");
        mass_proj->GetYaxis()->SetTitle("Count");

        if(i == 0) {
            raw_mass_2_4->cd();
            mass_proj->DrawClone("PMC");
            raw_mass->cd(i+1);
        }

        mass_proj->Draw("PMC");
    }

    //Doing like-sign BG subtraction for the sake of doing it (only in 2-4 GeV)

    TCanvas *ls_subtraction = new TCanvas("ls_subtraction", "LS Subtraction", 0, 10, 1920, 1080);

    ls_subtraction->cd();

    auto* k_star = (THnSparseF*)k_star_cor->Clone("k_star");
    auto* ls_k_star = (THnSparseF*)ls_k_star_cor->Clone("ls_k_star_cor");
    k_star->GetAxis(1)->SetRangeUser(2, 4);
    ls_k_star->GetAxis(1)->SetRangeUser(2, 4);

    auto* mass_proj_2_4 = k_star->Projection(4);
    auto* ls_mass_proj_2_4 = ls_k_star->Projection(4);

    mass_proj_2_4->SetTitle("Raw US and LS Inv Mass (no scaling)");

    mass_proj_2_4->SetLineColor(4);
    ls_mass_proj_2_4->SetLineColor(6);

    mass_proj_2_4->SetMarkerColor(4);
    ls_mass_proj_2_4->SetMarkerColor(6);

    mass_proj_2_4->DrawClone("PMC");
    ls_mass_proj_2_4->DrawClone("PMC SAME");

    //Looking into the D0 region
/*    TCanvas* d_mass = new TCanvas("d_mass", "Invariant Mass in Different D Regions", 0, 10, 1920, 1080);
    d_mass->Divide(2, 2);

    full_inv_mass->GetXaxis()->SetRangeUser(1.75, 2.1); //Pb-p D0 yield paper
    for(int i = 0; i < 4; i++) {

        d_mass->cd(i+1);

        gStyle->SetHistLineColor(5);
        gStyle->SetHistFillColor(5);

        full_inv_mass->GetYaxis()->SetRangeUser(pt_ranges[i], pt_ranges[i+1]);
        auto* d_mass_proj = full_inv_mass->ProjectionX();

        d_mass_proj->SetName(Form("d_mass_proj_%d_%d", pt_ranges[i], pt_ranges[i+1]));
        d_mass_proj->SetTitle(Form("%d GeV < p_{t} < %d GeV", pt_ranges[i], pt_ranges[i+1]));

        d_mass_proj->DrawClone("PMC");
    }



    //Looking into the full invariant mass range

    TCanvas* full_mass = new TCanvas("full_mass", "Full Invariant Mass Range (over all Pt), 0, 10, 1920, 1080");

    full_mass->cd();

    full_inv_mass->GetYaxis()->SetRangeUser(0, 0);
    full_inv_mass->GetXaxis()->SetRangeUser(0, 0);
    auto* full_mass_proj = full_inv_mass->ProjectionX();
    full_mass_proj->SetTitle("Full K-pi invariant mass range");
    full_mass_proj->Draw("PMC");
*/

    //Saving everything
    ls_subtraction->SaveAs("lambda/ls_raw_FS.png");
    raw_mass->SaveAs("lambda/raw_mass_divided_FS.png");
    raw_mass_2_4->SaveAs("lambda/raw_mass_2_4_FS.png");
   // d_mass->SaveAs("lambda/d_mass_divided_FS.png");
  //  full_mass->SaveAs("lambda/full_mass_FS.png");
}

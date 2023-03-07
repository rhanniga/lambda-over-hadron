import ROOT as rt
rt.gStyle.SetOptStat(0)

c = rt.TCanvas("c", "c", 800, 600)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.12)
c.SetTopMargin(0.05)

for PT_MODE in [0, 1, 2]:

    if PT_MODE == 0:
        avg6_infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        avg4_infile = rt.TFile("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        zyam_infile = rt.TFile("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        gaus_infile = rt.TFile("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        von_infile = rt.TFile("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        v2_infile = rt.TFile("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        avg6_infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        avg4_infile = rt.TFile("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        zyam_infile = rt.TFile("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        gaus_infile = rt.TFile("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        von_infile = rt.TFile("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        v2_infile = rt.TFile("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        avg6_infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        avg4_infile = rt.TFile("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        zyam_infile = rt.TFile("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        gaus_infile = rt.TFile("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        von_infile = rt.TFile("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        v2_infile = rt.TFile("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    file_dict = {"avg6": avg6_infile, "avg4": avg4_infile, "zyam": zyam_infile, "gaus": gaus_infile, "von": von_infile, "v2": v2_infile}

    for key, file in file_dict.items():

        mult_strings = ["_0_20", "_20_50", "_50_80"]
        h_h_colors = [rt.kRed + 2, rt.kBlue + 2, rt.kGreen + 2]
        h_lambda_colors = [rt.kRed - 3, rt.kAzure - 3, rt.kGreen - 2]


        color_index = 0
        for mult_string in mult_strings:
            h_lambda_dphi = file.Get("h_lambda_dphi_subtracted" + mult_string)
            h_h_dphi = file.Get("h_h_dphi" + mult_string)
            if key in ["avg6", "avg4", "zyam"]:
                h_lambda_ue_fit = file.Get("ue_line" + mult_string)
                h_h_ue_fit = file.Get("hh_ue_line" + mult_string)
            elif key in ["v2"]:
                h_lambda_ue_fit = file.Get("v2_fit" + mult_string)
                h_h_ue_fit = file.Get("hh_v2_fit" + mult_string)
            elif key in ["gaus"]:
                h_lambda_ue_fit = file.Get("ue_avg_fit" + mult_string)
                h_h_ue_fit = file.Get("hh_ue_avg_fit" + mult_string)
                h_lambda_total_fit = file.Get("fit_function" + mult_string)
                h_h_total_fit = file.Get("hh_fit_function" + mult_string)
            elif key in ["von"]:
                h_lambda_ue_fit = file.Get("v2_fit" + mult_string)
                h_h_ue_fit = file.Get("hh_v2_fit" + mult_string)
                h_lambda_total_fit = file.Get("von_fit" + mult_string)
                h_h_total_fit = file.Get("hh_von_fit" + mult_string)
            
            h_lambda_dphi.SetTitle("")
            h_lambda_dphi.GetYaxis().SetTitle("#frac{1}{N_{trig}}(h-#Lambda pairs per #Delta#varphi bin with |#Delta#eta| < 1.2)")
            h_lambda_dphi.GetYaxis().SetRangeUser(0.8*h_lambda_dphi.GetMinimum(), 1.2*h_lambda_dphi.GetMaximum())
            h_lambda_dphi.GetYaxis().SetMaxDigits(3)
            h_lambda_dphi.GetXaxis().SetTitle("#Delta#varphi_{h-#Lambda}")
            h_lambda_dphi.GetXaxis().SetTitleSize(0.05)
            h_lambda_dphi.GetXaxis().SetTitleOffset(1)
            h_lambda_dphi.SetLineColor(h_lambda_colors[color_index])
            h_lambda_dphi.SetLineWidth(2)
            h_lambda_dphi.SetMarkerColor(h_lambda_colors[color_index])
            h_lambda_dphi.SetMarkerStyle(43)
            h_lambda_dphi.SetMarkerSize(2)

            h_lambda_ue_fit.SetLineColor(rt.kGray + 3)
            h_lambda_ue_fit.SetLineWidth(2)
            h_lambda_ue_fit.SetLineStyle(rt.kDashed)

            if key in ["gaus", "von"]:
                h_lambda_total_fit.SetLineColor(rt.kGray + 3)
                h_lambda_total_fit.SetLineWidth(2)
                h_lambda_total_fit.SetLineStyle(rt.kSolid)
            
            h_h_dphi.SetTitle("")
            h_h_dphi.SetLineColor(h_h_colors[color_index])
            h_h_dphi.GetYaxis().SetTitle("#frac{1}{N_{trig}}(h-h pairs per #Delta#varphi bin with |#Delta#eta| < 1.2)")
            h_h_dphi.GetYaxis().SetRangeUser(0.8*h_h_dphi.GetMinimum(), 1.3*h_h_dphi.GetMaximum())
            h_h_dphi.GetYaxis().SetMaxDigits(3)
            h_h_dphi.GetXaxis().SetTitle("#Delta#varphi_{h-h}")
            h_h_dphi.GetXaxis().SetTitleOffset(1)
            h_h_dphi.GetXaxis().SetTitleSize(0.05)
            h_h_dphi.SetLineWidth(2)
            h_h_dphi.SetMarkerColor(h_h_colors[color_index])
            h_h_dphi.SetMarkerStyle(34)
            h_h_dphi.SetMarkerSize(1.5)

            h_h_ue_fit.SetLineColor(rt.kGray + 3)
            h_h_ue_fit.SetLineWidth(2)
            h_h_ue_fit.SetLineStyle(rt.kDashed)

            if key in ["gaus", "von"]:
                h_h_total_fit.SetLineColor(rt.kGray + 3)
                h_h_total_fit.SetLineWidth(2)
                h_h_total_fit.SetLineStyle(rt.kSolid)

            
            h_lambda_dphi.Draw()
            h_lambda_ue_fit.Draw("same")
            if key in ["gaus", "von"]:
                h_lambda_total_fit.Draw("same")
            leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            if PT_MODE == 0:
                if mult_string == "_0_20":
                    leg.AddEntry(h_lambda_dphi, "Data, 0-20% V0A mult. (main p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_lambda_dphi, "Data, 20-50% V0A mult. (main p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_lambda_dphi, "Data, 50-80% V0A mult. (main p_{T} bin)", "lep")
            elif PT_MODE == 1:
                if mult_string == "_0_20":
                    leg.AddEntry(h_lambda_dphi, "Data, 0-20% V0A mult. (low p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_lambda_dphi, "Data, 20-50% V0A mult. (low p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_lambda_dphi, "Data, 50-80% V0A mult. (low p_{T} bin)", "lep")
            elif PT_MODE == 2:
                if mult_string == "_0_20":
                    leg.AddEntry(h_lambda_dphi, "Data, 0-20% V0A mult. (high p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_lambda_dphi, "Data, 20-50% V0A mult. (high p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_lambda_dphi, "Data, 50-80% V0A mult. (high p_{T} bin)", "lep")
            if key == "avg6":
                leg.AddEntry(h_lambda_ue_fit, "6-bin avg.", "l")
            elif key == "avg4":
                leg.AddEntry(h_lambda_ue_fit, "4-bin avg.", "l")
            elif key == "zyam":
                leg.AddEntry(h_lambda_ue_fit, "ZYAM", "l")
            elif key == "v2":
                leg.AddEntry(h_lambda_ue_fit, "v_{2} assumption", "l")
            elif key == "gaus":
                leg.AddEntry(h_lambda_total_fit, "Gaus + pol0", "l")
                leg.AddEntry(h_lambda_ue_fit, "6-bin avg.", "l")
            elif key == "von":
                leg.AddEntry(h_lambda_total_fit, "Von Mises + pol0", "l")
                leg.AddEntry(h_lambda_ue_fit, "v_{2} assumption", "l")
            leg.Draw("SAME")
            c.Draw()
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_" + key + mult_string + ".pdf")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_" + key + mult_string + "_lowpt.pdf")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_" + key + mult_string + "_highpt.pdf")
            
            h_h_dphi.Draw()
            h_h_ue_fit.Draw("same")
            if key in ["gaus", "von"]:
                h_h_total_fit.Draw("same")
            leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            if PT_MODE == 0:
                if mult_string == "_0_20":
                    leg.AddEntry(h_h_dphi, "Data, 0-20% V0A mult. (main p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_h_dphi, "Data, 20-50% V0A mult. (main p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_h_dphi, "Data, 50-80% V0A mult. (main p_{T} bin)", "lep")
            elif PT_MODE == 1:
                if mult_string == "_0_20":
                    leg.AddEntry(h_h_dphi, "Data, 0-20% V0A mult. (low p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_h_dphi, "Data, 20-50% V0A mult. (low p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_h_dphi, "Data, 50-80% V0A mult. (low p_{T} bin)", "lep")
            elif PT_MODE == 2:
                if mult_string == "_0_20":
                    leg.AddEntry(h_h_dphi, "Data, 0-20% V0A mult. (high p_{T} bin)", "lep")
                elif mult_string == "_20_50":
                    leg.AddEntry(h_h_dphi, "Data, 20-50% V0A mult. (high p_{T} bin)", "lep")
                elif mult_string == "_50_80":
                    leg.AddEntry(h_h_dphi, "Data, 50-80% V0A mult. (high p_{T} bin)", "lep")
            if key == "avg6":
                leg.AddEntry(h_h_ue_fit, "6-bin avg.", "l")
            elif key == "avg4":
                leg.AddEntry(h_h_ue_fit, "4-bin avg.", "l")
            elif key == "zyam":
                leg.AddEntry(h_h_ue_fit, "ZYAM", "l")
            elif key == "v2":
                leg.AddEntry(h_h_ue_fit, "v_{2} assumption", "l")
            elif key == "gaus":
                leg.AddEntry(h_h_total_fit, "Gaus + pol0", "l")
                leg.AddEntry(h_h_ue_fit, "6-bin avg.", "l")
            elif key == "von":
                leg.AddEntry(h_h_total_fit, "Von Mises + pol0", "l")
                leg.AddEntry(h_h_ue_fit, "v_{2} assumption", "l")
            leg.Draw("SAME")
            c.Draw()
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_" + key + mult_string + ".pdf")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_" + key + mult_string + "_lowpt.pdf")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_" + key + mult_string + "_highpt.pdf")
            
            color_index += 1
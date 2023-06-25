import ROOT as rt
rt.gStyle.SetOptStat(0)

from array import array as arr


def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))
def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error


c = rt.TCanvas("c","c", 800, 600)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.15)


low_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
high_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
normal_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")


for i, f in enumerate([low_pt_infile, high_pt_infile, normal_pt_infile]):

    h_h_dphi_0_20 = f.Get("h_h_dphi_0_20")
    h_h_dphi_20_50 = f.Get("h_h_dphi_20_50")
    h_h_dphi_50_80 = f.Get("h_h_dphi_50_80")

    h_h_von_fit_0_20 = f.Get("hh_von_fit_0_20")
    h_h_von_fit_20_50 = f.Get("hh_von_fit_20_50")
    h_h_von_fit_50_80 = f.Get("hh_von_fit_50_80")


    h_lambda_dphi_0_20 = f.Get("h_lambda_dphi_subtracted_0_20")
    h_lambda_dphi_20_50 = f.Get("h_lambda_dphi_subtracted_20_50")
    h_lambda_dphi_50_80 = f.Get("h_lambda_dphi_subtracted_50_80")

    h_lambda_von_fit_0_20 = f.Get("von_fit_0_20")
    h_lambda_von_fit_20_50 = f.Get("von_fit_20_50")
    h_lambda_von_fit_50_80 = f.Get("von_fit_50_80")
    
    h_h_dphi_0_20.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-h} per #Delta#it{#varphi} bin)")
    h_h_dphi_20_50.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-h} per #Delta#it{#varphi} bin)")
    h_h_dphi_50_80.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-h} per #Delta#it{#varphi} bin)")

    h_lambda_dphi_0_20.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-#Lambda} per #Delta#it{#varphi} bin)")
    h_lambda_dphi_20_50.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-#Lambda} per #Delta#it{#varphi} bin)")
    h_lambda_dphi_50_80.SetTitle(";#Delta#it{#varphi};#frac{1}{N_{trig}}(N_{h-#Lambda} per #Delta#it{#varphi} bin)")

    h_h_dphi_0_20.SetLineColor(rt.kRed + 2)
    h_h_dphi_20_50.SetLineColor(rt.kBlue + 2)
    h_h_dphi_50_80.SetLineColor(rt.kMagenta + 2)
    h_h_dphi_0_20.SetMarkerColor(rt.kRed + 2)
    h_h_dphi_20_50.SetMarkerColor(rt.kBlue + 2)
    h_h_dphi_50_80.SetMarkerColor(rt.kMagenta + 2)
    h_h_dphi_0_20.SetLineWidth(2)
    h_h_dphi_20_50.SetLineWidth(2)
    h_h_dphi_50_80.SetLineWidth(2)
    h_h_dphi_0_20.SetMarkerStyle(43)
    h_h_dphi_20_50.SetMarkerStyle(43)
    h_h_dphi_50_80.SetMarkerStyle(43)
    h_h_dphi_0_20.GetYaxis().SetRangeUser(h_h_dphi_0_20.GetMinimum()*0.8, h_h_dphi_0_20.GetMaximum()*1.3)
    h_h_dphi_20_50.GetYaxis().SetRangeUser(h_h_dphi_20_50.GetMinimum()*0.8, h_h_dphi_20_50.GetMaximum()*1.3)
    h_h_dphi_50_80.GetYaxis().SetRangeUser(h_h_dphi_50_80.GetMinimum()*0.8, h_h_dphi_50_80.GetMaximum()*1.3)


    h_h_dphi_0_20.SetMarkerSize(1.5)
    h_h_dphi_20_50.SetMarkerSize(1.5)
    h_h_dphi_50_80.SetMarkerSize(1.5)
    h_h_von_fit_0_20.SetLineColor(rt.kRed + 2)
    h_h_von_fit_20_50.SetLineColor(rt.kBlue + 2)
    h_h_von_fit_50_80.SetLineColor(rt.kMagenta + 2)
    h_h_von_fit_0_20.SetLineWidth(2)
    h_h_von_fit_20_50.SetLineWidth(2)
    h_h_von_fit_50_80.SetLineWidth(2)
    h_h_von_fit_0_20.SetLineStyle(2)
    h_h_von_fit_20_50.SetLineStyle(2)
    h_h_von_fit_50_80.SetLineStyle(2)


    h_lambda_dphi_0_20.SetLineColor(rt.kRed - 4)
    h_lambda_dphi_20_50.SetLineColor(rt.kBlue - 4)
    h_lambda_dphi_50_80.SetLineColor(rt.kMagenta - 4)
    h_lambda_dphi_0_20.SetMarkerColor(rt.kRed - 4)
    h_lambda_dphi_20_50.SetMarkerColor(rt.kBlue - 4)
    h_lambda_dphi_50_80.SetMarkerColor(rt.kMagenta - 4)
    h_lambda_dphi_0_20.SetLineWidth(2)
    h_lambda_dphi_20_50.SetLineWidth(2)
    h_lambda_dphi_50_80.SetLineWidth(2)
    h_lambda_dphi_0_20.SetMarkerStyle(43)
    h_lambda_dphi_20_50.SetMarkerStyle(43)
    h_lambda_dphi_50_80.SetMarkerStyle(43)
    h_lambda_dphi_0_20.SetMarkerSize(1.5)
    h_lambda_dphi_20_50.SetMarkerSize(1.5)
    h_lambda_dphi_50_80.SetMarkerSize(1.5)
    h_lambda_dphi_0_20.GetYaxis().SetRangeUser(h_lambda_dphi_0_20.GetMinimum()*0.8, h_lambda_dphi_0_20.GetMaximum()*1.3)
    h_lambda_dphi_20_50.GetYaxis().SetRangeUser(h_lambda_dphi_20_50.GetMinimum()*0.8, h_lambda_dphi_20_50.GetMaximum()*1.3)
    h_lambda_dphi_50_80.GetYaxis().SetRangeUser(h_lambda_dphi_50_80.GetMinimum()*0.8, h_lambda_dphi_50_80.GetMaximum()*1.3)
    h_lambda_von_fit_0_20.SetLineColor(rt.kRed - 4)
    h_lambda_von_fit_20_50.SetLineColor(rt.kBlue - 4)
    h_lambda_von_fit_50_80.SetLineColor(rt.kMagenta - 4)
    h_lambda_von_fit_0_20.SetLineWidth(2)
    h_lambda_von_fit_20_50.SetLineWidth(2)
    h_lambda_von_fit_50_80.SetLineWidth(2)
    h_lambda_von_fit_0_20.SetLineStyle(2)
    h_lambda_von_fit_20_50.SetLineStyle(2)
    h_lambda_von_fit_50_80.SetLineStyle(2)


    h_h_0_20_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_h_0_20_leg.AddEntry(h_h_dphi_0_20, "h-h data, 0-20%", "lep")
    h_h_0_20_leg.AddEntry(h_h_von_fit_0_20, "Von mises fit", "l")
    h_h_0_20_leg.SetLineWidth(0)
    h_h_0_20_leg.SetBorderSize(0)
    h_h_20_50_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_h_20_50_leg.AddEntry(h_h_dphi_20_50, "h-h data, 20-50%", "lep")
    h_h_20_50_leg.AddEntry(h_h_von_fit_20_50, "Von mises fit", "l")
    h_h_20_50_leg.SetLineWidth(0)
    h_h_20_50_leg.SetBorderSize(0)
    h_h_50_80_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_h_50_80_leg.AddEntry(h_h_dphi_50_80, "h-h data, 50-80%", "lep")
    h_h_50_80_leg.AddEntry(h_h_von_fit_50_80, "Von mises fit", "l")
    h_h_50_80_leg.SetLineWidth(0)
    h_h_50_80_leg.SetBorderSize(0)

    if i == 0:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_0_20_lowpt.pdf")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_20_50_lowpt.pdf")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_50_80_lowpt.pdf")
    elif i == 1:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_0_20_highpt.pdf")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_20_50_highpt.pdf")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_50_80_highpt.pdf")
    else:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_0_20_midpt.pdf")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_20_50_midpt.pdf")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_h_dphi_with_von_50_80_midpt.pdf")

    if i == 0:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_0_20_lowpt.png")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_20_50_lowpt.png")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_50_80_lowpt.png")
    elif i == 1:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_0_20_highpt.png")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_20_50_highpt.png")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_50_80_highpt.png")
    else:
        h_h_dphi_0_20.Draw("PE")
        h_h_von_fit_0_20.Draw("SAME")
        h_h_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_0_20_midpt.png")
        h_h_dphi_20_50.Draw("PE")
        h_h_von_fit_20_50.Draw("SAME")
        h_h_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_20_50_midpt.png")
        h_h_dphi_50_80.Draw("PE")
        h_h_von_fit_50_80.Draw("SAME")
        h_h_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_h_dphi_with_von_50_80_midpt.png")




    h_lambda_0_20_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_lambda_0_20_leg.AddEntry(h_lambda_dphi_0_20, "h-#Lambda data, 0-20%", "lep")
    h_lambda_0_20_leg.AddEntry(h_lambda_von_fit_0_20, "Von mises fit", "l")
    h_lambda_0_20_leg.SetLineWidth(0)
    h_lambda_0_20_leg.SetBorderSize(0)
    h_lambda_20_50_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_lambda_20_50_leg.AddEntry(h_lambda_dphi_20_50, "h-#Lambda data, 20-50%", "lep")
    h_lambda_20_50_leg.AddEntry(h_lambda_von_fit_20_50, "Von mises fit", "l")
    h_lambda_20_50_leg.SetLineWidth(0)
    h_lambda_20_50_leg.SetBorderSize(0)
    h_lambda_50_80_leg = rt.TLegend(0.6, 0.75, 0.9, 0.9)
    h_lambda_50_80_leg.AddEntry(h_lambda_dphi_50_80, "h-#Lambda data, 50-80%", "lep")
    h_lambda_50_80_leg.AddEntry(h_lambda_von_fit_50_80, "Von mises fit", "l")
    h_lambda_50_80_leg.SetLineWidth(0)
    h_lambda_50_80_leg.SetBorderSize(0)

    if i == 0:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_0_20_lowpt.pdf")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_20_50_lowpt.pdf")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_50_80_lowpt.pdf")
    elif i == 1:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_0_20_highpt.pdf")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_20_50_highpt.pdf")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_50_80_highpt.pdf")
    else:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_0_20_midpt.pdf")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_20_50_midpt.pdf")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/h_lambda_dphi_with_von_50_80_midpt.pdf")

    if i == 0:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_0_20_lowpt.png")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_20_50_lowpt.png")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_50_80_lowpt.png")
    elif i == 1:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_0_20_highpt.png")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_20_50_highpt.png")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_50_80_highpt.png")
    else:
        h_lambda_dphi_0_20.Draw("PE")
        h_lambda_von_fit_0_20.Draw("SAME")
        h_lambda_0_20_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_0_20_midpt.png")
        h_lambda_dphi_20_50.Draw("PE")
        h_lambda_von_fit_20_50.Draw("SAME")
        h_lambda_20_50_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_20_50_midpt.png")
        h_lambda_dphi_50_80.Draw("PE")
        h_lambda_von_fit_50_80.Draw("SAME")
        h_lambda_50_80_leg.Draw("SAME")
        c.Draw()
        c.SaveAs("figures/width_tmp/png/h_lambda_dphi_with_von_50_80_midpt.png")



    # parameter 4 == near side, parameter 6 == away side

    # # only doing width comp with phi for 2-4 GeV bin
    # if i != 2:
    #     continue

    syst_infile = rt.TFile("../systematics_scripts/width_syst_out.root", "READ")

    if i == 0:
        dpmjet_infile = rt.TFile("../output/dpmjet_graphs_assoc_2_4.root")

        h_lambda_ns_syst_graph = syst_infile.Get("ns_total_syst")
        h_lambda_as_syst_graph = syst_infile.Get("as_total_syst")
        h_h_ns_syst_graph = syst_infile.Get("hh_ns_total_syst")
        h_h_as_syst_graph = syst_infile.Get("hh_as_total_syst")
    elif i == 1:
        dpmjet_infile = rt.TFile("../output/dpmjet_graphs_assoc_15_25.root")

        h_lambda_ns_syst_graph = syst_infile.Get("ns_total_syst_lowpt")
        h_lambda_as_syst_graph = syst_infile.Get("as_total_syst_lowpt")
        h_h_ns_syst_graph = syst_infile.Get("hh_ns_total_syst_lowpt")
        h_h_as_syst_graph = syst_infile.Get("hh_as_total_syst_lowpt")
    elif i == 2:
        dpmjet_infile = rt.TFile("../output/dpmjet_graphs_assoc_25_4.root")

        h_lambda_ns_syst_graph = syst_infile.Get("ns_total_syst_highpt")
        h_lambda_as_syst_graph = syst_infile.Get("as_total_syst_highpt")
        h_h_ns_syst_graph = syst_infile.Get("hh_ns_total_syst_highpt")
        h_h_as_syst_graph = syst_infile.Get("hh_as_total_syst_highpt")

    h_lambda_near_width_graph_dpmjet = dpmjet_infile.Get("h_lambda_near_width_graph")
    h_lambda_away_width_graph_dpmjet = dpmjet_infile.Get("h_lambda_away_width_graph")
    h_h_near_width_graph_dpmjet = dpmjet_infile.Get("h_h_near_width_graph")
    h_h_away_width_graph_dpmjet = dpmjet_infile.Get("h_h_away_width_graph")

    mult_array = arr('d', [35, 65, 90])
    mult_array_err = arr('d', [15, 15, 10])
    mult_array_err_syst = arr('d', [3, 3, 3])
    tmp_width_array_err = arr('d', [0, 0, 0])

    h_h_kappa_near_0_20 = h_h_von_fit_0_20.GetParameter(4)
    h_h_kappa_near_20_50 = h_h_von_fit_20_50.GetParameter(4)
    h_h_kappa_near_50_80 = h_h_von_fit_50_80.GetParameter(4)
    h_h_kappa_near_error_0_20 = h_h_von_fit_0_20.GetParError(4)
    h_h_kappa_near_error_20_50 = h_h_von_fit_20_50.GetParError(4)
    h_h_kappa_near_error_50_80 = h_h_von_fit_50_80.GetParError(4)
    h_h_width_near_0_20 = get_width_from_kappa(h_h_kappa_near_0_20)
    h_h_width_near_20_50 = get_width_from_kappa(h_h_kappa_near_20_50)
    h_h_width_near_50_80 = get_width_from_kappa(h_h_kappa_near_50_80)
    h_h_width_error_near_0_20 = get_width_error_from_kappa(h_h_kappa_near_0_20, h_h_kappa_near_error_0_20)
    h_h_width_error_near_20_50 = get_width_error_from_kappa(h_h_kappa_near_20_50, h_h_kappa_near_error_20_50)
    h_h_width_error_near_50_80 = get_width_error_from_kappa(h_h_kappa_near_50_80, h_h_kappa_near_error_50_80)
    h_h_near_width_array = arr('d', [h_h_width_near_50_80, h_h_width_near_20_50, h_h_width_near_0_20])
    h_h_near_width_error_array = arr('d', [h_h_width_error_near_50_80, h_h_width_error_near_20_50, h_h_width_error_near_0_20])

    h_h_near_width_syst_error_array = arr('d', [
        h_h_ns_syst_graph.GetY()[0]*h_h_width_near_50_80,
        h_h_ns_syst_graph.GetY()[1]*h_h_width_near_20_50,
        h_h_ns_syst_graph.GetY()[2]*h_h_width_near_0_20
    ])

    h_h_kappa_away_0_20 = h_h_von_fit_0_20.GetParameter(6)
    h_h_kappa_away_20_50 = h_h_von_fit_20_50.GetParameter(6)
    h_h_kappa_away_50_80 = h_h_von_fit_50_80.GetParameter(6)
    h_h_kappa_away_error_0_20 = h_h_von_fit_0_20.GetParError(6)
    h_h_kappa_away_error_20_50 = h_h_von_fit_20_50.GetParError(6)
    h_h_kappa_away_error_50_80 = h_h_von_fit_50_80.GetParError(6)
    h_h_width_away_0_20 = get_width_from_kappa(h_h_kappa_away_0_20)
    h_h_width_away_20_50 = get_width_from_kappa(h_h_kappa_away_20_50)
    h_h_width_away_50_80 = get_width_from_kappa(h_h_kappa_away_50_80)
    h_h_width_error_away_0_20 = get_width_error_from_kappa(h_h_kappa_away_0_20, h_h_kappa_away_error_0_20)
    h_h_width_error_away_20_50 = get_width_error_from_kappa(h_h_kappa_away_20_50, h_h_kappa_away_error_20_50)
    h_h_width_error_away_50_80 = get_width_error_from_kappa(h_h_kappa_away_50_80, h_h_kappa_away_error_50_80)
    h_h_away_width_array = arr('d', [h_h_width_away_50_80, h_h_width_away_20_50, h_h_width_away_0_20])
    h_h_away_width_error_array = arr('d', [h_h_width_error_away_50_80, h_h_width_error_away_20_50, h_h_width_error_away_0_20])

    h_h_away_width_syst_error_array = arr('d', [
        h_h_as_syst_graph.GetY()[0]*h_h_width_away_50_80,
        h_h_as_syst_graph.GetY()[1]*h_h_width_away_20_50,
        h_h_as_syst_graph.GetY()[2]*h_h_width_away_0_20
    ])


    h_lambda_kappa_near_0_20 = h_lambda_von_fit_0_20.GetParameter(4)
    h_lambda_kappa_near_20_50 = h_lambda_von_fit_20_50.GetParameter(4)
    h_lambda_kappa_near_50_80 = h_lambda_von_fit_50_80.GetParameter(4)
    h_lambda_kappa_near_error_0_20 = h_lambda_von_fit_0_20.GetParError(4)
    h_lambda_kappa_near_error_20_50 = h_lambda_von_fit_20_50.GetParError(4)
    h_lambda_kappa_near_error_50_80 = h_lambda_von_fit_50_80.GetParError(4)
    h_lambda_width_near_0_20 = get_width_from_kappa(h_lambda_kappa_near_0_20)
    h_lambda_width_near_20_50 = get_width_from_kappa(h_lambda_kappa_near_20_50)
    h_lambda_width_near_50_80 = get_width_from_kappa(h_lambda_kappa_near_50_80)
    h_lambda_width_error_near_0_20 = get_width_error_from_kappa(h_lambda_kappa_near_0_20, h_lambda_kappa_near_error_0_20)
    h_lambda_width_error_near_20_50 = get_width_error_from_kappa(h_lambda_kappa_near_20_50, h_lambda_kappa_near_error_20_50)
    h_lambda_width_error_near_50_80 = get_width_error_from_kappa(h_lambda_kappa_near_50_80, h_lambda_kappa_near_error_50_80)
    h_lambda_near_width_array = arr('d', [h_lambda_width_near_50_80, h_lambda_width_near_20_50, h_lambda_width_near_0_20])
    h_lambda_near_width_error_array = arr('d', [h_lambda_width_error_near_50_80, h_lambda_width_error_near_20_50, h_lambda_width_error_near_0_20])

    h_lambda_near_width_syst_error_array = arr('d', [
        h_lambda_ns_syst_graph.GetY()[0]*h_lambda_width_near_50_80,
        h_lambda_ns_syst_graph.GetY()[1]*h_lambda_width_near_20_50,
        h_lambda_ns_syst_graph.GetY()[2]*h_lambda_width_near_0_20
    ])

    h_lambda_kappa_away_0_20 = h_lambda_von_fit_0_20.GetParameter(6)
    h_lambda_kappa_away_20_50 = h_lambda_von_fit_20_50.GetParameter(6)
    h_lambda_kappa_away_50_80 = h_lambda_von_fit_50_80.GetParameter(6)
    h_lambda_kappa_away_error_0_20 = h_lambda_von_fit_0_20.GetParError(6)
    h_lambda_kappa_away_error_20_50 = h_lambda_von_fit_20_50.GetParError(6)
    h_lambda_kappa_away_error_50_80 = h_lambda_von_fit_50_80.GetParError(6)
    h_lambda_width_away_0_20 = get_width_from_kappa(h_lambda_kappa_away_0_20)
    h_lambda_width_away_20_50 = get_width_from_kappa(h_lambda_kappa_away_20_50)
    h_lambda_width_away_50_80 = get_width_from_kappa(h_lambda_kappa_away_50_80)
    h_lambda_width_error_away_0_20 = get_width_error_from_kappa(h_lambda_kappa_away_0_20, h_lambda_kappa_away_error_0_20)
    h_lambda_width_error_away_20_50 = get_width_error_from_kappa(h_lambda_kappa_away_20_50, h_lambda_kappa_away_error_20_50)
    h_lambda_width_error_away_50_80 = get_width_error_from_kappa(h_lambda_kappa_away_50_80, h_lambda_kappa_away_error_50_80)
    h_lambda_away_width_array = arr('d', [h_lambda_width_away_50_80, h_lambda_width_away_20_50, h_lambda_width_away_0_20])
    h_lambda_away_width_error_array = arr('d', [h_lambda_width_error_away_50_80, h_lambda_width_error_away_20_50, h_lambda_width_error_away_0_20])

    h_lambda_away_width_syst_error_array = arr('d', [
        h_lambda_as_syst_graph.GetY()[0]*h_lambda_width_away_50_80,
        h_lambda_as_syst_graph.GetY()[1]*h_lambda_width_away_20_50,
        h_lambda_as_syst_graph.GetY()[2]*h_lambda_width_away_0_20
    ])

    h_phi_fit_file = rt.TFile("../output/h_phi_fits.root", "READ")
    h_phi_von_fit_0_20 = h_phi_fit_file.Get("h_phi_von_fit_0_20")
    h_phi_von_fit_20_50 = h_phi_fit_file.Get("h_phi_von_fit_20_50")
    h_phi_von_fit_50_80 = h_phi_fit_file.Get("h_phi_von_fit_50_80")
    h_phi_fit_file.Close()



    h_phi_kappa_near_0_20 = h_phi_von_fit_0_20.GetParameter(4)
    h_phi_kappa_near_20_50 = h_phi_von_fit_20_50.GetParameter(4)
    h_phi_kappa_near_50_80 = h_phi_von_fit_50_80.GetParameter(4)
    h_phi_kappa_near_error_0_20 = h_phi_von_fit_0_20.GetParError(4)
    h_phi_kappa_near_error_20_50 = h_phi_von_fit_20_50.GetParError(4)
    h_phi_kappa_near_error_50_80 = h_phi_von_fit_50_80.GetParError(4)
    h_phi_width_near_0_20 = get_width_from_kappa(h_phi_kappa_near_0_20)
    h_phi_width_near_20_50 = get_width_from_kappa(h_phi_kappa_near_20_50)
    h_phi_width_near_50_80 = get_width_from_kappa(h_phi_kappa_near_50_80)
    h_phi_width_error_near_0_20 = get_width_error_from_kappa(h_phi_kappa_near_0_20, h_phi_kappa_near_error_0_20)
    h_phi_width_error_near_20_50 = get_width_error_from_kappa(h_phi_kappa_near_20_50, h_phi_kappa_near_error_20_50)
    h_phi_width_error_near_50_80 = get_width_error_from_kappa(h_phi_kappa_near_50_80, h_phi_kappa_near_error_50_80)
    h_phi_near_width_array = arr('d', [h_phi_width_near_50_80, h_phi_width_near_20_50, h_phi_width_near_0_20])
    h_phi_near_width_error_array = arr('d', [h_phi_width_error_near_50_80, h_phi_width_error_near_20_50, h_phi_width_error_near_0_20])

    h_phi_kappa_away_0_20 = h_phi_von_fit_0_20.GetParameter(6)
    h_phi_kappa_away_20_50 = h_phi_von_fit_20_50.GetParameter(6)
    h_phi_kappa_away_50_80 = h_phi_von_fit_50_80.GetParameter(6)
    h_phi_kappa_away_error_0_20 = h_phi_von_fit_0_20.GetParError(6)
    h_phi_kappa_away_error_20_50 = h_phi_von_fit_20_50.GetParError(6)
    h_phi_kappa_away_error_50_80 = h_phi_von_fit_50_80.GetParError(6)
    h_phi_width_away_0_20 = get_width_from_kappa(h_phi_kappa_away_0_20)
    h_phi_width_away_20_50 = get_width_from_kappa(h_phi_kappa_away_20_50)
    h_phi_width_away_50_80 = get_width_from_kappa(h_phi_kappa_away_50_80)
    h_phi_width_error_away_0_20 = get_width_error_from_kappa(h_phi_kappa_away_0_20, h_phi_kappa_away_error_0_20)
    h_phi_width_error_away_20_50 = get_width_error_from_kappa(h_phi_kappa_away_20_50, h_phi_kappa_away_error_20_50)
    h_phi_width_error_away_50_80 = get_width_error_from_kappa(h_phi_kappa_away_50_80, h_phi_kappa_away_error_50_80)
    h_phi_away_width_array = arr('d', [h_phi_width_away_50_80, h_phi_width_away_20_50, h_phi_width_away_0_20])
    h_phi_away_width_error_array = arr('d', [h_phi_width_error_away_50_80, h_phi_width_error_away_20_50, h_phi_width_error_away_0_20])

    print("THE 0-20 FIT Chi2/NDF: ", h_phi_von_fit_0_20.GetChisquare()/h_phi_von_fit_0_20.GetNDF())
    print("THE 20-50 FIT Chi2/NDF: ", h_phi_von_fit_20_50.GetChisquare()/h_phi_von_fit_20_50.GetNDF())
    print("THE 50-80 FIT Chi2/NDF: ", h_phi_von_fit_50_80.GetChisquare()/h_phi_von_fit_50_80.GetNDF())

    print("THE 0-20 WIDTH ERROR: ", h_phi_width_error_away_0_20)
    print("THE 20-50 WIDTH ERROR: ", h_phi_width_error_away_20_50)
    print("THE 50-80 WIDTH ERROR: ", h_phi_width_error_away_50_80)


    
    # h_h_near_width_graph = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err, tmp_width_array_err)
    h_h_near_width_graph = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err, h_h_near_width_error_array)
    h_h_near_width_graph_syst = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err_syst, h_h_near_width_syst_error_array)

    h_h_near_width_graph.SetLineColor(rt.kRed+3)
    h_h_near_width_graph.SetLineWidth(2)
    h_h_near_width_graph.SetMarkerStyle(20)
    h_h_near_width_graph.SetMarkerSize(2)
    h_h_near_width_graph.SetMarkerColor(rt.kRed+3)

    h_h_near_width_graph_dpmjet.SetLineColor(rt.kRed+3)
    h_h_near_width_graph_dpmjet.SetLineWidth(1)
    h_h_near_width_graph_dpmjet.SetMarkerStyle(0)
    h_h_near_width_graph_dpmjet.SetMarkerSize(0)
    h_h_near_width_graph_dpmjet.SetFillColor(rt.kRed+3)
    h_h_near_width_graph_dpmjet.SetFillStyle(3144)

    h_h_near_width_graph_syst.SetMarkerStyle(21)
    h_h_near_width_graph_syst.SetMarkerSize(0)
    h_h_near_width_graph_syst.SetMarkerColor(rt.kRed+3)
    h_h_near_width_graph_syst.SetLineColor(rt.kRed+3)
    h_h_near_width_graph_syst.SetLineWidth(2)
    h_h_near_width_graph_syst.SetFillStyle(0)

    # h_h_away_width_graph = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err, tmp_width_array_err)
    h_h_away_width_graph = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err, h_h_away_width_error_array)
    h_h_away_width_graph_syst = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err_syst, h_h_away_width_syst_error_array)

    h_h_away_width_graph.SetLineColor(rt.kViolet - 6)
    h_h_away_width_graph.SetLineWidth(2)
    h_h_away_width_graph.SetMarkerStyle(20)
    h_h_away_width_graph.SetMarkerSize(2)
    h_h_away_width_graph.SetMarkerColor(rt.kViolet - 6)

    h_h_away_width_graph_dpmjet.SetLineColor(rt.kViolet-6)
    h_h_away_width_graph_dpmjet.SetLineWidth(1)
    h_h_away_width_graph_dpmjet.SetMarkerStyle(0)
    h_h_away_width_graph_dpmjet.SetMarkerSize(0)
    h_h_away_width_graph_dpmjet.SetFillColor(rt.kViolet-6)
    h_h_away_width_graph_dpmjet.SetFillStyle(3144)

    h_h_away_width_graph_syst.SetMarkerStyle(21)
    h_h_away_width_graph_syst.SetMarkerSize(0)
    h_h_away_width_graph_syst.SetMarkerColor(rt.kViolet - 6)
    h_h_away_width_graph_syst.SetLineColor(rt.kViolet - 6)
    h_h_away_width_graph_syst.SetLineWidth(2)
    h_h_away_width_graph_syst.SetFillStyle(0)


    # h_lambda_near_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err, tmp_width_array_err)
    h_lambda_near_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err, h_lambda_near_width_error_array)
    h_lambda_near_width_graph_syst = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err_syst, h_lambda_near_width_syst_error_array)

    h_lambda_near_width_graph.SetLineColor(rt.kOrange-1)
    h_lambda_near_width_graph.SetLineWidth(2)
    h_lambda_near_width_graph.SetMarkerStyle(43)
    h_lambda_near_width_graph.SetMarkerSize(2.5)
    h_lambda_near_width_graph.SetMarkerColor(rt.kOrange-1)

    h_lambda_near_width_graph_dpmjet.SetLineColor(rt.kOrange-1)
    h_lambda_near_width_graph_dpmjet.SetLineWidth(1)
    h_lambda_near_width_graph_dpmjet.SetMarkerStyle(0)
    h_lambda_near_width_graph_dpmjet.SetMarkerSize(0)
    h_lambda_near_width_graph_dpmjet.SetFillColor(rt.kOrange-1)
    h_lambda_near_width_graph_dpmjet.SetFillStyle(3144)

    h_lambda_near_width_graph_syst.SetMarkerStyle(21)
    h_lambda_near_width_graph_syst.SetMarkerSize(0)
    h_lambda_near_width_graph_syst.SetMarkerColor(rt.kOrange-1)
    h_lambda_near_width_graph_syst.SetLineColor(rt.kOrange-1)
    h_lambda_near_width_graph_syst.SetLineWidth(2)
    h_lambda_near_width_graph_syst.SetFillStyle(0)

    h_lambda_away_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_away_width_array, mult_array_err, h_lambda_away_width_error_array)
    h_lambda_away_width_graph_syst = rt.TGraphErrors(3, mult_array, h_lambda_away_width_array, mult_array_err_syst, h_lambda_away_width_syst_error_array)

    h_lambda_away_width_graph.SetLineColor(rt.kAzure-1)
    h_lambda_away_width_graph.SetLineWidth(2)
    h_lambda_away_width_graph.SetMarkerStyle(43)
    h_lambda_away_width_graph.SetMarkerSize(2.5)
    h_lambda_away_width_graph.SetMarkerColor(rt.kAzure-1)

    h_lambda_away_width_graph_dpmjet.SetLineColor(rt.kAzure-1)
    h_lambda_away_width_graph_dpmjet.SetLineWidth(1)
    h_lambda_away_width_graph_dpmjet.SetMarkerStyle(0)
    h_lambda_away_width_graph_dpmjet.SetMarkerSize(0)
    h_lambda_away_width_graph_dpmjet.SetFillColor(rt.kAzure-1)
    h_lambda_away_width_graph_dpmjet.SetFillStyle(3144)

    h_lambda_away_width_graph_syst.SetMarkerStyle(21)
    h_lambda_away_width_graph_syst.SetMarkerSize(0)
    h_lambda_away_width_graph_syst.SetMarkerColor(rt.kAzure-1)
    h_lambda_away_width_graph_syst.SetLineColor(rt.kAzure-1)
    h_lambda_away_width_graph_syst.SetLineWidth(2)
    h_lambda_away_width_graph_syst.SetFillStyle(0)

    # h_phi_near_width_graph = rt.TGraphErrors(3, mult_array, h_phi_near_width_array, mult_array_err, tmp_width_array_err)
    h_phi_near_width_graph = rt.TGraphErrors(3, mult_array, h_phi_near_width_array, mult_array_err, h_phi_near_width_error_array)
    h_phi_near_width_graph.SetLineColor(rt.kRed)
    h_phi_near_width_graph.SetLineWidth(2)
    h_phi_near_width_graph.SetMarkerStyle(23)
    h_phi_near_width_graph.SetMarkerSize(2)
    h_phi_near_width_graph.SetMarkerColor(rt.kRed)
    # h_phi_away_width_graph = rt.TGraphErrors(3, mult_array, h_phi_away_width_array, mult_array_err, tmp_width_array_err)
    h_phi_away_width_graph = rt.TGraphErrors(3, mult_array, h_phi_away_width_array, mult_array_err, h_phi_away_width_error_array)
    h_phi_away_width_graph.SetLineColor(rt.kBlue)
    h_phi_away_width_graph.SetLineWidth(2)
    h_phi_away_width_graph.SetMarkerStyle(23)
    h_phi_away_width_graph.SetMarkerSize(2)
    h_phi_away_width_graph.SetMarkerColor(rt.kBlue)


    mult_bin_widths = arr('d', [0.0, 20.0, 50.0, 80.0, 100.0])
    plotting_hist = rt.TH1D("plotting_hist", "", 4, mult_bin_widths)


    plotting_hist.SetMarkerStyle(20)
    plotting_hist.SetMarkerSize(1)
    plotting_hist.SetMarkerColor(rt.kRed+1)
    plotting_hist.SetLineColor(rt.kRed+2)
    plotting_hist.SetLineWidth(2)
    plotting_hist.GetXaxis().SetTitle("Multiplicity (%)")
    plotting_hist.GetYaxis().SetTitle("Von Mises width")
    plotting_hist.GetXaxis().SetTitleSize(0.05)
    plotting_hist.GetXaxis().SetLabelSize(0.045)
    plotting_hist.GetYaxis().SetLabelSize(0.045)
    plotting_hist.GetXaxis().SetTitleOffset(1.2)
    plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
    plotting_hist.SetStats(0)
    plotting_hist.Draw("PE")

    plotting_hist.GetXaxis().SetLabelOffset(999)
    plotting_hist.GetXaxis().SetTickSize(0)

    if i == 0:
        plotting_hist.GetYaxis().SetRangeUser(0.1, 1.1)
    else:
        plotting_hist.GetYaxis().SetRangeUser(0.1, 1.1)
    rt.gPad.Update()
    new_axis = rt.TGaxis(rt.gPad.GetUxmax(),
            rt.gPad.GetUymin(),
            rt.gPad.GetUxmin(),
            rt.gPad.GetUymin(),
            plotting_hist.GetXaxis().GetXmin(),
            plotting_hist.GetXaxis().GetXmax(),
            510,"-")
    new_axis.SetLabelSize(0.045)
    new_axis.SetLabelFont(plotting_hist.GetXaxis().GetLabelFont())
    new_axis.SetLabelOffset(-0.03)
    new_axis.Draw()

    h_h_width_legend = rt.TLegend(0.17, 0.7, 0.33, 0.9)
    h_h_width_legend.AddEntry(h_h_near_width_graph, "h-h, near-side", "p")
    h_h_width_legend.AddEntry(h_h_away_width_graph, "h-h, away-side", "p")
    h_h_width_legend.SetBorderSize(0)
    h_h_width_legend.Draw()
    h_h_near_width_graph.Draw("sameP")
    h_h_away_width_graph.Draw("sameP")
    h_h_near_width_graph_syst.Draw("same e2")
    h_h_away_width_graph_syst.Draw("same e2")
    h_h_near_width_graph_dpmjet.Draw("same e3")
    h_h_away_width_graph_dpmjet.Draw("same e3")

    h_lambda_width_legend = rt.TLegend(0.35, 0.7, 0.51, 0.9)
    h_lambda_width_legend.AddEntry(h_lambda_near_width_graph, "h-#Lambda, near-side", "p")
    h_lambda_width_legend.AddEntry(h_lambda_away_width_graph, "h-#Lambda, away-side", "p")
    h_lambda_width_legend.SetBorderSize(0)
    h_lambda_width_legend.Draw()

    h_lambda_near_width_graph.Draw("sameP")
    h_lambda_away_width_graph.Draw("sameP")
    h_lambda_near_width_graph_syst.Draw("same e2")
    h_lambda_away_width_graph_syst.Draw("same e2")
    h_lambda_near_width_graph_dpmjet.Draw("same e3")
    h_lambda_away_width_graph_dpmjet.Draw("same e3")

    dpmjet_width_legend = rt.TLegend(0.55, 0.75, 0.69, 0.85)
    tmp_graph = h_lambda_near_width_graph_dpmjet.Clone("tmp_graph")
    tmp_graph.SetLineColor(rt.kGray+1)
    tmp_graph.SetFillColor(rt.kGray+1)
    dpmjet_width_legend.AddEntry(tmp_graph, "DPMJet", "f")
    dpmjet_width_legend.SetBorderSize(0)
    dpmjet_width_legend.Draw()

    # h_phi_width_legend = rt.TLegend(0.51, 0.7, 0.68, 0.9)
    # h_phi_width_legend.AddEntry(h_phi_near_width_graph, "h-#phi, near-side", "p")
    # h_phi_width_legend.AddEntry(h_phi_away_width_graph, "h-#phi, away-side", "p")
    # h_phi_width_legend.SetBorderSize(0)
    # h_phi_width_legend.Draw()
    # h_phi_near_width_graph.Draw("sameP")
    # h_phi_away_width_graph.Draw("sameP")

    c.Draw()
    if i == 0:
        c.SaveAs("figures/width_tmp/von_mises_widths_lowpt.pdf")
    elif i == 1:
        c.SaveAs("figures/width_tmp/von_mises_widths_highpt.pdf")
    else:
        c.SaveAs("figures/width_tmp/von_mises_widths.pdf")


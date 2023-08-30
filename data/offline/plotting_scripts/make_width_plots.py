import ROOT as rt
rt.gStyle.SetOptStat(0)

import math

from array import array as arr

DIHADRON_NEAR_COLOR = rt.kRed + 2
DIHADRON_AWAY_COLOR = rt.kBlue + 2

LAMBDA_NEAR_COLOR = rt.kPink + 7
LAMBDA_AWAY_COLOR = rt.kAzure + 7

def get_ratio_error(numerator, denominator, numerator_error, denominator_error):
    return math.sqrt((numerator_error/denominator)**2 + (numerator*denominator_error/denominator**2)**2)

def get_ratio_graph(numerator_graph, denominator_graph, name):

    ratio_graph = rt.TGraphAsymmErrors()

    flat_error = numerator_graph.GetErrorY(1)

    for i in range(denominator_graph.GetN()):

        x = denominator_graph.GetX()[i]
        y = denominator_graph.GetY()[i]
        x_error = denominator_graph.GetErrorX(i)
        y_error = denominator_graph.GetErrorY(i)

        numerator_y = numerator_graph.Eval(x)

        ratio = numerator_y/y
        ratio_error = get_ratio_error(numerator_y, y, flat_error, y_error)

        ratio_graph.SetPoint(i, x, ratio)
        ratio_graph.SetPointError(i, x_error, x_error, ratio_error, ratio_error)

    ratio_graph.SetName(name)

    return ratio_graph.Clone()

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

PAD_WIDTH = 0.43
LEFT_MARGIN = 0.23
TOP_MARGIN = 0.05
BOTTOM_MARGIN = 0.30
FIRST_PAD_WIDTH = PAD_WIDTH + LEFT_MARGIN*PAD_WIDTH
LABEL_SIZE = 0.045
TITLE_SIZE = 0.05
TITLE_OFFSET = 1.35
LABEL_OFFSET = 0.01
RATIO_PAD_HEIGHT = 0.25
FIRST_RATIO_PAD_HEIGHT = RATIO_PAD_HEIGHT + BOTTOM_MARGIN*RATIO_PAD_HEIGHT

unity_line = rt.TLine(3, 1, 53, 1)
unity_line.SetLineStyle(2)
unity_line.SetLineColor(rt.kBlack)


width_final_canvas = rt.TCanvas("width_final_canvas", "", 1500, 800)
width_final_canvas.SetMargin(0,0,0,0)
width_final_canvas_new_x_axis = rt.TCanvas("width_final_canvas_new_x_axis", "", 1500, 800)
width_final_canvas_new_x_axis.SetMargin(0,0,0,0)

low_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
high_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
normal_pt_infile = rt.TFile("../output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")


c.cd()
for i, f in enumerate([normal_pt_infile, low_pt_infile, high_pt_infile]):

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

    if i == 1:
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
    elif i == 2:
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

    if i == 1:
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
    elif i == 2:
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

    if i == 1:
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
    elif i == 2:
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

    if i == 1:
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
    elif i == 2:
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

    h_lambda_near_width_graph_new_x_axis_dpmjet = dpmjet_infile.Get("h_lambda_near_width_graph_new_x_axis")
    h_lambda_away_width_graph_new_x_axis_dpmjet = dpmjet_infile.Get("h_lambda_away_width_graph_new_x_axis")
    h_h_near_width_graph_new_x_axis_dpmjet = dpmjet_infile.Get("h_h_near_width_graph_new_x_axis")
    h_h_away_width_graph_new_x_axis_dpmjet = dpmjet_infile.Get("h_h_away_width_graph_new_x_axis")

    mult_array = arr('d', [35, 65, 90])
    mult_array_err = arr('d', [15, 15, 10])
    mult_array_err_syst = arr('d', [3, 3, 3])
    tmp_width_array_err = arr('d', [0, 0, 0])

    nch_0_20 = 42.4
    nch_0_20_err = 0.9
    nch_20_50 = 27.6
    nch_20_50_err = 0.5
    nch_50_80 = 17.7
    nch_50_80_err = 0.4
    nch_array = arr("d", [nch_50_80, nch_20_50, nch_0_20])
    nch_array_err = arr("d", [nch_50_80_err, nch_20_50_err, nch_0_20_err])

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



    print(f"----------------------PT MODE: {i}----------------------")
    print("---------------h lambda width values----------------")
    print(f"0-20\% & {h_lambda_width_near_0_20:.2e} $\\pm$ {h_lambda_width_error_near_0_20:.2e} &  {h_lambda_width_away_0_20:.2e} $\\pm$ {h_lambda_width_error_away_0_20:.2e} \\\\")
    print(f"20-50\% & {h_lambda_width_near_20_50:.2e} $\\pm$ {h_lambda_width_error_near_20_50:.2e} &  {h_lambda_width_away_20_50:.2e} $\\pm$ {h_lambda_width_error_away_20_50:.2e} \\\\")
    print(f"50-80\% & {h_lambda_width_near_50_80:.2e} $\\pm$ {h_lambda_width_error_near_50_80:.2e} &  {h_lambda_width_away_50_80:.2e} $\\pm$ {h_lambda_width_error_away_50_80:.2e} \\\\")
    print("---------------h h width values----------------")
    print(f"0-20\% & {h_h_width_near_0_20:.2e} $\\pm$ {h_h_width_error_near_0_20:.2e} &  {h_h_width_away_0_20:.2e} $\\pm$ {h_h_width_error_away_0_20:.2e} \\\\")
    print(f"20-50\% & {h_h_width_near_20_50:.2e} $\\pm$ {h_h_width_error_near_20_50:.2e} &  {h_h_width_away_20_50:.2e} $\\pm$ {h_h_width_error_away_20_50:.2e} \\\\")
    print(f"50-80\% & {h_h_width_near_50_80:.2e} $\\pm$ {h_h_width_error_near_50_80:.2e} &  {h_h_width_away_50_80:.2e} $\\pm$ {h_h_width_error_away_50_80:.2e} \\\\")


    
    # h_h_near_width_graph = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err, tmp_width_array_err)
    h_h_near_width_graph = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err, h_h_near_width_error_array).Clone()
    h_h_near_width_graph_syst = rt.TGraphErrors(3, mult_array, h_h_near_width_array, mult_array_err_syst, h_h_near_width_syst_error_array).Clone()

    h_h_near_width_graph_new_x_axis = rt.TGraphErrors(3, nch_array, h_h_near_width_array, nch_array_err, h_h_near_width_error_array).Clone()
    h_h_near_width_graph_new_x_axis_syst = rt.TGraphErrors(3, nch_array, h_h_near_width_array, nch_array_err, h_h_near_width_syst_error_array).Clone()

    h_h_near_width_graph.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph.SetLineWidth(2)
    h_h_near_width_graph.SetMarkerStyle(20)
    h_h_near_width_graph.SetMarkerSize(2)
    h_h_near_width_graph.SetMarkerColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_dpmjet.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_dpmjet.SetLineWidth(1)
    h_h_near_width_graph_dpmjet.SetMarkerStyle(0)
    h_h_near_width_graph_dpmjet.SetMarkerSize(0)
    h_h_near_width_graph_dpmjet.SetFillColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_dpmjet.SetFillStyle(3144)
    h_h_near_width_graph_syst.SetMarkerStyle(21)
    h_h_near_width_graph_syst.SetMarkerSize(0)
    h_h_near_width_graph_syst.SetMarkerColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_syst.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_syst.SetLineWidth(2)
    h_h_near_width_graph_syst.SetFillStyle(0)

    h_h_near_width_graph_new_x_axis.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis.SetLineWidth(2)
    h_h_near_width_graph_new_x_axis.SetMarkerStyle(20)
    h_h_near_width_graph_new_x_axis.SetMarkerSize(2)
    h_h_near_width_graph_new_x_axis.SetMarkerColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis_dpmjet.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis_dpmjet.SetLineWidth(1)
    h_h_near_width_graph_new_x_axis_dpmjet.SetMarkerStyle(0)
    h_h_near_width_graph_new_x_axis_dpmjet.SetMarkerSize(0)
    h_h_near_width_graph_new_x_axis_dpmjet.SetFillColorAlpha(DIHADRON_NEAR_COLOR, 0.5)
    h_h_near_width_graph_new_x_axis_dpmjet.SetFillStyle(3144)
    h_h_near_width_graph_new_x_axis_syst.SetMarkerStyle(21)
    h_h_near_width_graph_new_x_axis_syst.SetMarkerSize(0)
    h_h_near_width_graph_new_x_axis_syst.SetMarkerColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis_syst.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis_syst.SetLineWidth(2)
    h_h_near_width_graph_new_x_axis_syst.SetFillStyle(0)

    # h_h_away_width_graph = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err, tmp_width_array_err)
    h_h_away_width_graph = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err, h_h_away_width_error_array).Clone()
    h_h_away_width_graph_syst = rt.TGraphErrors(3, mult_array, h_h_away_width_array, mult_array_err_syst, h_h_away_width_syst_error_array).Clone()
    h_h_away_width_graph_new_x_axis = rt.TGraphErrors(3, nch_array, h_h_away_width_array, nch_array_err, h_h_away_width_error_array).Clone()
    h_h_away_width_graph_new_x_axis_syst = rt.TGraphErrors(3, nch_array, h_h_away_width_array, nch_array_err, h_h_away_width_syst_error_array).Clone()

    h_h_away_width_graph.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph.SetLineWidth(2)
    h_h_away_width_graph.SetMarkerStyle(20)
    h_h_away_width_graph.SetMarkerSize(2)
    h_h_away_width_graph.SetMarkerColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_dpmjet.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_dpmjet.SetLineWidth(1)
    h_h_away_width_graph_dpmjet.SetMarkerStyle(0)
    h_h_away_width_graph_dpmjet.SetMarkerSize(0)
    h_h_away_width_graph_dpmjet.SetFillColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_dpmjet.SetFillStyle(3144)
    h_h_away_width_graph_syst.SetMarkerStyle(21)
    h_h_away_width_graph_syst.SetMarkerSize(0)
    h_h_away_width_graph_syst.SetMarkerColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_syst.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_syst.SetLineWidth(2)
    h_h_away_width_graph_syst.SetFillStyle(0)

    h_h_away_width_graph_new_x_axis.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis.SetLineWidth(2)
    h_h_away_width_graph_new_x_axis.SetMarkerStyle(20)
    h_h_away_width_graph_new_x_axis.SetMarkerSize(2)
    h_h_away_width_graph_new_x_axis.SetMarkerColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis_dpmjet.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis_dpmjet.SetLineWidth(1)
    h_h_away_width_graph_new_x_axis_dpmjet.SetMarkerStyle(0)
    h_h_away_width_graph_new_x_axis_dpmjet.SetMarkerSize(0)
    h_h_away_width_graph_new_x_axis_dpmjet.SetFillColorAlpha(DIHADRON_AWAY_COLOR, 0.5)
    h_h_away_width_graph_new_x_axis_dpmjet.SetFillStyle(3144)
    h_h_away_width_graph_new_x_axis_syst.SetMarkerStyle(21)
    h_h_away_width_graph_new_x_axis_syst.SetMarkerSize(0)
    h_h_away_width_graph_new_x_axis_syst.SetMarkerColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis_syst.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis_syst.SetLineWidth(2)
    h_h_away_width_graph_new_x_axis_syst.SetFillStyle(0)


    # h_lambda_near_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err, tmp_width_array_err)
    h_lambda_near_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err, h_lambda_near_width_error_array).Clone()
    h_lambda_near_width_graph_syst = rt.TGraphErrors(3, mult_array, h_lambda_near_width_array, mult_array_err_syst, h_lambda_near_width_syst_error_array).Clone()

    h_lambda_near_width_graph_new_x_axis = rt.TGraphErrors(3, nch_array, h_lambda_near_width_array, nch_array_err, h_lambda_near_width_error_array).Clone()
    h_lambda_near_width_graph_new_x_axis_syst = rt.TGraphErrors(3, nch_array, h_lambda_near_width_array, nch_array_err, h_lambda_near_width_syst_error_array).Clone()

    h_lambda_near_width_graph.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph.SetLineWidth(2)
    h_lambda_near_width_graph.SetMarkerStyle(43)
    h_lambda_near_width_graph.SetMarkerSize(2.5)
    h_lambda_near_width_graph.SetMarkerColor(LAMBDA_NEAR_COLOR)

    h_lambda_near_width_graph_dpmjet.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_dpmjet.SetLineWidth(1)
    h_lambda_near_width_graph_dpmjet.SetMarkerStyle(0)
    h_lambda_near_width_graph_dpmjet.SetMarkerSize(0)
    h_lambda_near_width_graph_dpmjet.SetFillColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_dpmjet.SetFillStyle(3144)

    h_lambda_near_width_graph_syst.SetMarkerStyle(21)
    h_lambda_near_width_graph_syst.SetMarkerSize(0)
    h_lambda_near_width_graph_syst.SetMarkerColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_syst.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_syst.SetLineWidth(2)
    h_lambda_near_width_graph_syst.SetFillStyle(0)

    h_lambda_near_width_graph_new_x_axis.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis.SetLineWidth(2)
    h_lambda_near_width_graph_new_x_axis.SetMarkerStyle(21)
    h_lambda_near_width_graph_new_x_axis.SetMarkerSize(2)
    h_lambda_near_width_graph_new_x_axis.SetMarkerColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetLineWidth(1)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetMarkerStyle(0)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetMarkerSize(0)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetFillColorAlpha(LAMBDA_NEAR_COLOR, 0.5)
    h_lambda_near_width_graph_new_x_axis_dpmjet.SetFillStyle(3144)
    h_lambda_near_width_graph_new_x_axis_syst.SetMarkerStyle(21)
    h_lambda_near_width_graph_new_x_axis_syst.SetMarkerSize(0)
    h_lambda_near_width_graph_new_x_axis_syst.SetMarkerColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis_syst.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis_syst.SetLineWidth(2)
    h_lambda_near_width_graph_new_x_axis_syst.SetFillStyle(0)

    h_lambda_away_width_graph = rt.TGraphErrors(3, mult_array, h_lambda_away_width_array, mult_array_err, h_lambda_away_width_error_array).Clone()
    h_lambda_away_width_graph_syst = rt.TGraphErrors(3, mult_array, h_lambda_away_width_array, mult_array_err_syst, h_lambda_away_width_syst_error_array).Clone()

    h_lambda_away_width_graph_new_x_axis = rt.TGraphErrors(3, nch_array, h_lambda_away_width_array, nch_array_err, h_lambda_away_width_error_array).Clone()
    h_lambda_away_width_graph_new_x_axis_syst = rt.TGraphErrors(3, nch_array, h_lambda_away_width_array, nch_array_err, h_lambda_away_width_syst_error_array).Clone()

    h_lambda_away_width_graph.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph.SetLineWidth(2)
    h_lambda_away_width_graph.SetMarkerStyle(22)
    h_lambda_away_width_graph.SetMarkerSize(2.5)
    h_lambda_away_width_graph.SetMarkerColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_dpmjet.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_dpmjet.SetLineWidth(1)
    h_lambda_away_width_graph_dpmjet.SetMarkerStyle(0)
    h_lambda_away_width_graph_dpmjet.SetMarkerSize(0)
    h_lambda_away_width_graph_dpmjet.SetFillColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_dpmjet.SetFillStyle(3144)
    h_lambda_away_width_graph_syst.SetMarkerStyle(21)
    h_lambda_away_width_graph_syst.SetMarkerSize(0)
    h_lambda_away_width_graph_syst.SetMarkerColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_syst.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_syst.SetLineWidth(2)
    h_lambda_away_width_graph_syst.SetFillStyle(0)

    h_lambda_away_width_graph_new_x_axis.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis.SetLineWidth(2)
    h_lambda_away_width_graph_new_x_axis.SetMarkerStyle(21)
    h_lambda_away_width_graph_new_x_axis.SetMarkerSize(2)
    h_lambda_away_width_graph_new_x_axis.SetMarkerColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetLineWidth(1)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetMarkerStyle(0)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetMarkerSize(0)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetFillColorAlpha(LAMBDA_AWAY_COLOR, 0.5)
    h_lambda_away_width_graph_new_x_axis_dpmjet.SetFillStyle(3144)
    h_lambda_away_width_graph_new_x_axis_syst.SetMarkerStyle(21)
    h_lambda_away_width_graph_new_x_axis_syst.SetMarkerSize(0)
    h_lambda_away_width_graph_new_x_axis_syst.SetMarkerColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis_syst.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis_syst.SetLineWidth(2)
    h_lambda_away_width_graph_new_x_axis_syst.SetFillStyle(0)

    h_lambda_near_width_graph_comparison_dpmjet = get_ratio_graph(h_lambda_near_width_graph_dpmjet, h_lambda_near_width_graph, "h_lambda_near_width_graph_comparison_dpmjet")
    h_lambda_near_width_graph_comparison_dpmjet.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_comparison_dpmjet.SetLineWidth(1)
    h_lambda_near_width_graph_comparison_dpmjet.SetMarkerStyle(0)
    h_lambda_near_width_graph_comparison_dpmjet.SetMarkerSize(0)
    h_lambda_near_width_graph_comparison_dpmjet.SetFillColorAlpha(LAMBDA_NEAR_COLOR, 0.4)
    h_lambda_near_width_graph_comparison_dpmjet.SetFillStyle(3144)
    h_lambda_away_width_graph_comparison_dpmjet = get_ratio_graph(h_lambda_away_width_graph_dpmjet, h_lambda_away_width_graph, "h_lambda_away_width_graph_comparison_dpmjet")
    h_lambda_away_width_graph_comparison_dpmjet.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_comparison_dpmjet.SetLineWidth(1)
    h_lambda_away_width_graph_comparison_dpmjet.SetMarkerStyle(0)
    h_lambda_away_width_graph_comparison_dpmjet.SetMarkerSize(0)
    h_lambda_away_width_graph_comparison_dpmjet.SetFillColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_comparison_dpmjet.SetFillStyle(3144)
    h_h_near_width_graph_comparison_dpmjet = get_ratio_graph(h_h_near_width_graph_dpmjet, h_h_near_width_graph, "h_h_near_width_graph_comparison_dpmjet")
    h_h_near_width_graph_comparison_dpmjet.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_comparison_dpmjet.SetLineWidth(1)
    h_h_near_width_graph_comparison_dpmjet.SetMarkerStyle(0)
    h_h_near_width_graph_comparison_dpmjet.SetMarkerSize(0)
    h_h_near_width_graph_comparison_dpmjet.SetFillColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_comparison_dpmjet.SetFillStyle(3144)
    h_h_away_width_graph_comparison_dpmjet = get_ratio_graph(h_h_away_width_graph_dpmjet, h_h_away_width_graph, "h_h_away_width_graph_comparison_dpmjet")
    h_h_away_width_graph_comparison_dpmjet.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_comparison_dpmjet.SetLineWidth(1)
    h_h_away_width_graph_comparison_dpmjet.SetMarkerStyle(0)
    h_h_away_width_graph_comparison_dpmjet.SetMarkerSize(0)
    h_h_away_width_graph_comparison_dpmjet.SetFillColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_comparison_dpmjet.SetFillStyle(3144)

    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet = get_ratio_graph(h_lambda_near_width_graph_new_x_axis_dpmjet, h_lambda_near_width_graph_new_x_axis, "h_lambda_near_width_graph_new_x_axis_comparison_dpmjet")
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetLineWidth(1)
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetMarkerStyle(0)
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetMarkerSize(0)
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetFillColorAlpha(LAMBDA_NEAR_COLOR, 0.5)
    h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.SetFillStyle(3144)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet = get_ratio_graph(h_lambda_away_width_graph_new_x_axis_dpmjet, h_lambda_away_width_graph_new_x_axis, "h_lambda_away_width_graph_new_x_axis_comparison_dpmjet")
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetLineWidth(1)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetMarkerStyle(0)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetMarkerSize(0)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetFillColorAlpha(LAMBDA_AWAY_COLOR, 0.5)
    h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.SetFillStyle(3144)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet = get_ratio_graph(h_h_near_width_graph_new_x_axis_dpmjet, h_h_near_width_graph_new_x_axis, "h_h_near_width_graph_new_x_axis_comparison_dpmjet")
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetLineWidth(1)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetMarkerStyle(0)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetMarkerSize(0)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetFillColorAlpha(DIHADRON_NEAR_COLOR, 0.5)
    h_h_near_width_graph_new_x_axis_comparison_dpmjet.SetFillStyle(3144)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet = get_ratio_graph(h_h_away_width_graph_new_x_axis_dpmjet, h_h_away_width_graph_new_x_axis, "h_h_away_width_graph_new_x_axis_comparison_dpmjet")
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetLineWidth(1)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetMarkerStyle(0)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetMarkerSize(0)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetFillColorAlpha(DIHADRON_AWAY_COLOR, 0.5)
    h_h_away_width_graph_new_x_axis_comparison_dpmjet.SetFillStyle(3144)

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
    plotting_hist.GetYaxis().SetTitle("Extracted width")
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
    dpmjet_width_legend.AddEntry(tmp_graph, "DPMJET", "f")
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
    if i == 1:
        c.SaveAs("figures/width_tmp/von_mises_widths_lowpt_with_dpmjet.pdf")
    elif i == 2:
        c.SaveAs("figures/width_tmp/von_mises_widths_highpt_with_dpmjet.pdf")
    else:
        c.SaveAs("figures/width_tmp/von_mises_widths_with_dpmjet.pdf")


    plotting_hist = rt.TH1D("plotting_hist_new_x_axis", "", 65, 0, 65)
    plotting_hist.GetXaxis().SetTitle("#LT#frac{d#it{N}_{ch}}{d#it{#eta}}#GT_{|#it{#eta}_{lab}| < 0.5}")
    plotting_hist.GetXaxis().SetTitleSize(0.05)
    plotting_hist.GetXaxis().SetLabelSize(0.045)
    plotting_hist.GetYaxis().SetLabelSize(0.045)
    plotting_hist.GetXaxis().SetTitleOffset(1.35)
    plotting_hist.GetXaxis().SetRangeUser(3, 53)

    plotting_hist.GetYaxis().SetTitle("Extracted width")
    plotting_hist.GetXaxis().SetTitleSize(0.05)
    plotting_hist.GetXaxis().SetLabelSize(0.045)
    plotting_hist.GetYaxis().SetLabelSize(0.045)
    plotting_hist.GetXaxis().SetTitleOffset(1.2)
    plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
    plotting_hist.SetStats(0)
    if i == 0:
        plotting_hist.GetYaxis().SetRangeUser(0.13, 1.1)
    else:
        plotting_hist.GetYaxis().SetRangeUser(0.13, 1.1)
    plotting_hist.Draw("PE")
    left_plotting_hist = plotting_hist.Clone("left_plotting_hist")
    right_plotting_hist = plotting_hist.Clone("right_plotting_hist")


    h_h_width_legend = rt.TLegend(0.17, 0.7, 0.33, 0.9)
    h_h_width_legend.AddEntry(h_h_near_width_graph, "h-h, near-side", "p")
    h_h_width_legend.AddEntry(h_h_away_width_graph, "h-h, away-side", "p")
    h_h_width_legend.SetBorderSize(0)
    h_h_width_legend.Draw()

    h_h_near_width_graph_new_x_axis.Draw("sameP")
    h_h_away_width_graph_new_x_axis.Draw("sameP")
    h_h_near_width_graph_new_x_axis_syst.Draw("same e2")
    h_h_away_width_graph_new_x_axis_syst.Draw("same e2")
    h_h_near_width_graph_new_x_axis_dpmjet.Draw("same e3")
    h_h_away_width_graph_new_x_axis_dpmjet.Draw("same e3")

    h_lambda_width_legend = rt.TLegend(0.35, 0.7, 0.51, 0.9)
    h_lambda_width_legend.AddEntry(h_lambda_near_width_graph, "h-#Lambda, near-side", "p")
    h_lambda_width_legend.AddEntry(h_lambda_away_width_graph, "h-#Lambda, away-side", "p")
    h_lambda_width_legend.SetBorderSize(0)
    h_lambda_width_legend.Draw()

    h_lambda_near_width_graph_new_x_axis.Draw("sameP")
    h_lambda_away_width_graph_new_x_axis.Draw("sameP")
    h_lambda_near_width_graph_new_x_axis_syst.Draw("same e2")
    h_lambda_away_width_graph_new_x_axis_syst.Draw("same e2")
    h_lambda_near_width_graph_new_x_axis_dpmjet.Draw("same e3")
    h_lambda_away_width_graph_new_x_axis_dpmjet.Draw("same e3")

    h_lambda_near_width_fit = rt.TF1("h_lambda_near_width_fit", "pol1", 7, 51)
    h_lambda_near_width_fit.SetLineColor(LAMBDA_NEAR_COLOR)
    h_lambda_near_width_fit.SetLineStyle(2)
    h_lambda_near_width_fit.SetLineWidth(2)
    h_lambda_near_width_graph_new_x_axis_syst.Fit(h_lambda_near_width_fit, "R")
    h_lambda_away_width_fit = rt.TF1("h_lambda_away_width_fit", "pol1", 7, 51)
    h_lambda_away_width_fit.SetLineColor(LAMBDA_AWAY_COLOR)
    h_lambda_away_width_fit.SetLineStyle(2)
    h_lambda_away_width_fit.SetLineWidth(2)
    h_lambda_away_width_graph_new_x_axis_syst.Fit(h_lambda_away_width_fit, "R")

    h_h_near_width_fit = rt.TF1("h_h_near_width_fit", "pol1", 7, 51)
    h_h_near_width_fit.SetLineColor(DIHADRON_NEAR_COLOR)
    h_h_near_width_fit.SetLineStyle(2)
    h_h_near_width_fit.SetLineWidth(2)
    h_h_near_width_graph_new_x_axis_syst.Fit(h_h_near_width_fit, "R")
    h_h_away_width_fit = rt.TF1("h_h_away_width_fit", "pol1", 7, 51)
    h_h_away_width_fit.SetLineColor(DIHADRON_AWAY_COLOR)
    h_h_away_width_fit.SetLineStyle(2)
    h_h_away_width_fit.SetLineWidth(2)
    h_h_away_width_graph_new_x_axis_syst.Fit(h_h_away_width_fit, "R")

    dpmjet_width_legend = rt.TLegend(0.55, 0.75, 0.69, 0.85)
    tmp_graph = h_lambda_near_width_graph_dpmjet.Clone("tmp_graph")
    tmp_graph.SetLineColor(rt.kGray+1)
    tmp_graph.SetFillColor(rt.kGray+1)
    dpmjet_width_legend.AddEntry(tmp_graph, "DPMJET", "f")
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
    if i == 1:
        c.SaveAs("figures/width_tmp/von_mises_widths_lowpt_with_dpmjet_new_x_axis.pdf")
    elif i == 2:
        c.SaveAs("figures/width_tmp/von_mises_widths_highpt_with_dpmjet_new_x_axis.pdf")
    else:
        c.SaveAs("figures/width_tmp/von_mises_widths_with_dpmjet_new_x_axis.pdf")

    
    if i == 0:
        continue

    width_final_canvas_new_x_axis.cd()
    cur_c = width_final_canvas_new_x_axis

    if i == 1:
        left_pad = rt.TPad("lpad", "", 0, FIRST_RATIO_PAD_HEIGHT, FIRST_PAD_WIDTH, 1)
        left_pad.SetMargin(LEFT_MARGIN, 0, 0, TOP_MARGIN)
        left_pad.Draw()
        left_pad.cd()
        left_plotting_hist.GetYaxis().SetTitleOffset(TITLE_OFFSET - 0.2)
        left_plotting_hist.GetYaxis().SetTitleSize(TITLE_SIZE)
        left_plotting_hist.GetXaxis().SetRangeUser(3, 53)
        left_plotting_hist.DrawCopy("PE")
        label_x_start = 0.70
        label_y_start = 0.99
        label_text_space = 0.07
        alice_data_label = rt.TLatex()
        alice_data_label.SetNDC()
        alice_data_label.SetTextSize(0.06)
        alice_data_label.SetTextAlign(13)
        # alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
        alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{ALICE p#minusPb}")
        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{#sqrt{#it{s}_{NN}} = 5.02 TeV}")
        pt_label_x_start = 0.28
        pt_label_y_start = 0.93
        pt_label_text_space = 0.07
        pt_label = rt.TLatex()
        pt_label.SetNDC()
        pt_label.SetTextSize(0.05)
        pt_label.SetTextAlign(13)
        pt_label.DrawLatex(pt_label_x_start, pt_label_y_start - 1*pt_label_text_space, "#bf{4.0 < #it{p}^{h}_{T,trig}    <  8.0 GeV/#it{c}}")
        pt_label.DrawLatex(pt_label_x_start, pt_label_y_start - 2*pt_label_text_space, "#bf{|#Delta#it{#eta}| < 1.2}")
        assoc_pt_label_x_start = 0.28
        assoc_pt_label_y_start = 0.93
        assoc_pt_label_text_space = 0.07
        assoc_pt_label = rt.TLatex()
        assoc_pt_label.SetNDC()
        assoc_pt_label.SetTextSize(0.05)
        assoc_pt_label.SetTextAlign(13)
        assoc_pt_label.DrawLatex(assoc_pt_label_x_start, assoc_pt_label_y_start, "#bf{1.5 < #it{p}^{#Lambda, h}_{T,assoc} < 2.5 GeV/#it{c}}")


    else:
        right_pad = rt.TPad("rpad", "", FIRST_PAD_WIDTH, FIRST_RATIO_PAD_HEIGHT, FIRST_PAD_WIDTH + PAD_WIDTH, 1)
        right_pad.SetMargin(0, 0, 0, TOP_MARGIN)
        scale = FIRST_PAD_WIDTH / PAD_WIDTH
        right_pad.Draw()
        right_pad.cd()
        right_plotting_hist.GetXaxis().SetLabelSize(scale*LABEL_SIZE)
        right_plotting_hist.GetXaxis().SetTitleSize(scale*TITLE_SIZE)
        right_plotting_hist.GetXaxis().SetLabelOffset(LABEL_OFFSET - 0.014)
        right_plotting_hist.GetXaxis().SetTitleOffset(TITLE_OFFSET - 0.25)
        right_plotting_hist.GetXaxis().SetRangeUser(3, 53)
        right_plotting_hist.DrawCopy("PE")
        assoc_pt_label_x_start = 0.04
        assoc_pt_label_y_start = 0.93
        assoc_pt_label_text_space = 0.07
        assoc_pt_label = rt.TLatex()
        assoc_pt_label.SetNDC()
        assoc_pt_label.SetTextSize(scale*0.041)
        assoc_pt_label.SetTextAlign(13)
        assoc_pt_label.DrawLatex(assoc_pt_label_x_start, assoc_pt_label_y_start, "#bf{2.5 < #it{p}^{#Lambda, h}_{T,assoc} < 4.0 GeV/#it{c}}")

        as_legend = rt.TLegend(0.6, 0.7, 1.0, 0.9)
        as_legend.SetBorderSize(0)
        as_legend.SetFillStyle(0)
        as_legend.AddEntry(h_lambda_away_width_graph_new_x_axis, "h-#Lambda, away-side", "PL")
        as_legend.AddEntry(h_h_away_width_graph_new_x_axis, "h-h, away-side", "PL")
        as_legend.Draw("SAME")

        ns_legend = rt.TLegend(0.6, 0.25, 1.0, 0.45)
        ns_legend.SetBorderSize(0)
        ns_legend.SetFillStyle(0)
        ns_legend.AddEntry(h_lambda_near_width_graph_new_x_axis, "h-#Lambda, near-side", "PL")
        ns_legend.AddEntry(h_h_near_width_graph_new_x_axis, "h-h, near-side", "PL")
        ns_legend.Draw("SAME")

        dpmjet_legend = rt.TLegend(0.1, 0.65, 0.4, 0.78)
        dpmjet_legend.SetBorderSize(0)
        dpmjet_legend.SetFillStyle(0)
        tmp_graph.SetFillStyle(3244)
        fit_line = rt.TLine()
        fit_line.SetLineWidth(2)
        fit_line.SetLineColor(rt.kBlack)
        fit_line.SetLineStyle(2)
        dpmjet_legend.AddEntry(tmp_graph, "DPMJET", "F")
        dpmjet_legend.AddEntry(fit_line, "Fit (data)", "L")
        dpmjet_legend.Draw("SAME")



    h_lambda_near_width_graph_new_x_axis_dpmjet.Draw("E3 SAME")
    h_lambda_away_width_graph_new_x_axis_dpmjet.Draw("E3 SAME")

    h_h_near_width_graph_new_x_axis_dpmjet.Draw("E3 SAME")
    h_h_away_width_graph_new_x_axis_dpmjet.Draw("E3 SAME")

    h_lambda_near_width_graph_new_x_axis.Draw("PE SAME")
    h_lambda_near_width_graph_new_x_axis_syst.Draw("E2 SAME")
    h_lambda_away_width_graph_new_x_axis.Draw("PE SAME")
    h_lambda_away_width_graph_new_x_axis_syst.Draw("E2 SAME")
    h_lambda_near_width_fit.Draw("SAME")
    h_lambda_away_width_fit.Draw("SAME")
    h_h_near_width_graph_new_x_axis.Draw("PE SAME")
    h_h_near_width_graph_new_x_axis_syst.Draw("E2 SAME")
    h_h_near_width_fit.Draw("SAME")
    h_h_away_width_graph_new_x_axis.Draw("PE SAME")
    h_h_away_width_graph_new_x_axis_syst.Draw("E2 SAME")
    h_h_away_width_fit.Draw("SAME")

    print("######################################################")
    print(f"pt mode: {i}")
    print(f"lambda near-side slope: {h_lambda_near_width_fit.GetParameter(1)} +/- {h_lambda_near_width_fit.GetParError(1)}")
    print(f"lambda away-side slope: {h_lambda_away_width_fit.GetParameter(1)} +/- {h_lambda_away_width_fit.GetParError(1)}")
    print(f"h near-side slope: {h_h_near_width_fit.GetParameter(1)} +/- {h_h_near_width_fit.GetParError(1)}")
    print(f"h away-side slope: {h_h_away_width_fit.GetParameter(1)} +/- {h_h_away_width_fit.GetParError(1)}")
    print("######################################################")
    width_final_canvas_new_x_axis.cd()
    cur_c = width_final_canvas_new_x_axis
    if i == 1:
        bottom_left_pad = rt.TPad("blpad", "", 0, 0, FIRST_PAD_WIDTH, FIRST_RATIO_PAD_HEIGHT)
        bottom_left_pad.SetMargin(LEFT_MARGIN, 0, BOTTOM_MARGIN, 0)
        bottom_left_pad.SetTicks(1, 1)
        bottom_left_pad.Draw()
        bottom_left_pad.cd()
        test_plotting_hist = rt.TH1D("test_plotting_left_hist", "", 65, 0, 65)
        for i in range(1, 66):
            test_plotting_hist.SetBinContent(i, -9999999)
            test_plotting_hist.SetBinError(i, 0)

        test_plotting_hist.GetXaxis().SetRangeUser(3, 53)
        test_plotting_hist.GetXaxis().SetLabelSize(0.1)
        test_plotting_hist.GetXaxis().SetTitle("<d#it{N}_{ch}/d#it{#eta}>")
        test_plotting_hist.GetXaxis().SetTitleSize(0.13)
        test_plotting_hist.GetXaxis().SetTitleOffset(0.93)

        test_plotting_hist.GetYaxis().SetLabelSize(0.07)
        test_plotting_hist.GetYaxis().SetTitle("Model/Data")
        test_plotting_hist.GetYaxis().SetTitleSize(0.1)
        test_plotting_hist.GetYaxis().SetTitleOffset(0.5)
        test_plotting_hist.GetYaxis().SetRangeUser(0.22, 1.32)
        test_plotting_hist.DrawCopy("PE")

        h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_h_near_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_h_away_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        unity_line.Draw("SAME")

    else:
        bottom_right_pad = rt.TPad("brpad", "", FIRST_PAD_WIDTH, 0, FIRST_PAD_WIDTH + PAD_WIDTH, FIRST_RATIO_PAD_HEIGHT)
        bottom_right_pad.SetMargin(0, 0, BOTTOM_MARGIN, 0)
        bottom_right_pad.SetTicks(1, 1)
        bottom_right_pad.Draw()
        bottom_right_pad.cd()
        test_plotting_hist = rt.TH1D("test_plotting_right_hist", "", 65, 0, 65)
        for i in range(1, 66):
            test_plotting_hist.SetBinContent(i, -9999999)
            test_plotting_hist.SetBinError(i, 0)
        test_plotting_hist.GetXaxis().SetRangeUser(3, 53)
        test_plotting_hist.GetXaxis().SetLabelSize(0.1)
        test_plotting_hist.GetXaxis().SetTitle("<d#it{N}_{ch}/d#it{#eta}>")
        test_plotting_hist.GetXaxis().SetTitleSize(0.13)
        test_plotting_hist.GetXaxis().SetTitleOffset(0.93)
        test_plotting_hist.GetYaxis().SetRangeUser(0.22, 1.32)
        test_plotting_hist.DrawCopy("PE")

        h_lambda_near_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_lambda_away_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_h_near_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        h_h_away_width_graph_new_x_axis_comparison_dpmjet.Draw("E3 SAME")
        unity_line.Draw("SAME")

        cur_c.SaveAs("final_width_plot_new_x_axis_model_ratio.pdf")
    c.cd()
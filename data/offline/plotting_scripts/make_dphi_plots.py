import ROOT as rt
import math
import array as arr
import time
rt.gStyle.SetOptStat(0)


colors = [rt.kGreen - 2, rt.kRed, rt.kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow, rt.kGray, rt.kViolet, rt.kTeal, rt.kSpring, rt.kAzure, rt.kPink, rt.kCopper, rt.kOrange+7, rt.kSpring+9, rt.kTeal-7, rt.kAzure+2, rt.kPink+10, rt.kCopper+3, rt.kOrange-3, rt.kSpring-5, rt.kTeal+3, rt.kAzure-6, rt.kPink-7, rt.kCopper-9]

total_dphi_sys_0_20 = 0.031
total_dphi_sys_20_50 = 0.032
total_dphi_sys_50_80 = 0.036
hh_total_dphi_sys_0_20 = 0.035
hh_total_dphi_sys_20_50 = 0.035
hh_total_dphi_sys_50_80 = 0.035


TOTAL_PAD_LENGTH = 1200
TOTAL_PAD_HEIGHT = 1050
dphi_all_for_real = rt.TCanvas(f"dphi_all_for_real", f"dphi_all_for_real", 50, 50, TOTAL_PAD_LENGTH, TOTAL_PAD_HEIGHT)
dphi_all_for_real.SetMargin(0, 0, 0, 0)

for PT_MODE in [0, 1, 2]:

    c = rt.TCanvas(f"c_{PT_MODE}",f"c_{PT_MODE}",800,600)
    c.SetLeftMargin(0.15)

    if PT_MODE == 0:
        central_file = rt.TFile("../output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        central_file = rt.TFile("../output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        central_file = rt.TFile("../output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    file_dict = {
        "central value": central_file
    }



    dphi_dict_0_20 = {}
    dphi_dict_20_50 = {}
    dphi_dict_50_80 = {}

    hh_dphi_dict_0_20 = {}
    hh_dphi_dict_20_50 = {}
    hh_dphi_dict_50_80 = {}

    for key in file_dict:
        dphi_dict_0_20[key] = file_dict[key].Get("h_lambda_dphi_subtracted_0_20")
        dphi_dict_20_50[key] = file_dict[key].Get("h_lambda_dphi_subtracted_20_50")
        dphi_dict_50_80[key] = file_dict[key].Get("h_lambda_dphi_subtracted_50_80")

        bin_width = dphi_dict_0_20[key].GetBinWidth(1)

        dphi_dict_0_20[key].Scale(1/bin_width)
        dphi_dict_20_50[key].Scale(1/bin_width)
        dphi_dict_50_80[key].Scale(1/bin_width)

        dphi_dict_0_20[key].GetYaxis().SetRangeUser(0.8*dphi_dict_0_20[key].GetMinimum(), 1.2*dphi_dict_0_20[key].GetMaximum())
        dphi_dict_20_50[key].GetYaxis().SetRangeUser(0.8*dphi_dict_20_50[key].GetMinimum(), 1.2*dphi_dict_20_50[key].GetMaximum())
        dphi_dict_50_80[key].GetYaxis().SetRangeUser(0.8*dphi_dict_50_80[key].GetMinimum(), 1.2*dphi_dict_50_80[key].GetMaximum())

        dphi_dict_0_20[key].SetTitle("")
        dphi_dict_20_50[key].SetTitle("")
        dphi_dict_50_80[key].SetTitle("")

        dphi_dict_0_20[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        dphi_dict_0_20[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")
        dphi_dict_20_50[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        dphi_dict_20_50[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")
        dphi_dict_50_80[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        dphi_dict_50_80[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")

        hh_dphi_dict_0_20[key] = file_dict[key].Get("h_h_dphi_0_20")
        hh_dphi_dict_20_50[key] = file_dict[key].Get("h_h_dphi_20_50")
        hh_dphi_dict_50_80[key] = file_dict[key].Get("h_h_dphi_50_80")

        hh_bin_width = hh_dphi_dict_0_20[key].GetBinWidth(1)

        hh_dphi_dict_0_20[key].Scale(1/hh_bin_width)
        hh_dphi_dict_20_50[key].Scale(1/hh_bin_width)
        hh_dphi_dict_50_80[key].Scale(1/hh_bin_width)

        hh_dphi_dict_0_20[key].GetYaxis().SetRangeUser(0.8*hh_dphi_dict_0_20[key].GetMinimum(), 1.2*hh_dphi_dict_0_20[key].GetMaximum())
        hh_dphi_dict_20_50[key].GetYaxis().SetRangeUser(0.8*hh_dphi_dict_20_50[key].GetMinimum(), 1.2*hh_dphi_dict_20_50[key].GetMaximum())
        hh_dphi_dict_50_80[key].GetYaxis().SetRangeUser(0.8*hh_dphi_dict_50_80[key].GetMinimum(), 1.2*hh_dphi_dict_50_80[key].GetMaximum())

        hh_dphi_dict_0_20[key].SetTitle("")
        hh_dphi_dict_20_50[key].SetTitle("")
        hh_dphi_dict_50_80[key].SetTitle("")

        hh_dphi_dict_0_20[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        hh_dphi_dict_0_20[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")
        hh_dphi_dict_20_50[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        hh_dphi_dict_20_50[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")
        hh_dphi_dict_50_80[key].GetXaxis().SetTitle("#Delta#it{#varphi}")
        hh_dphi_dict_50_80[key].GetYaxis().SetTitle("#frac{1}{N_{trig}} #frac{dN_{pair}}{d#Delta#it{#varphi}}")




    central_dphi_0_20_just_syst = dphi_dict_0_20["central value"].Clone("central_dphi_0_20")
    central_dphi_20_50_just_syst = dphi_dict_20_50["central value"].Clone("central_dphi_20_50")
    central_dphi_50_80_just_syst = dphi_dict_50_80["central value"].Clone("central_dphi_50_80")

    hh_central_dphi_0_20_just_syst = hh_dphi_dict_0_20["central value"].Clone("hh_central_dphi_0_20")
    hh_central_dphi_20_50_just_syst = hh_dphi_dict_20_50["central value"].Clone("hh_central_dphi_20_50")
    hh_central_dphi_50_80_just_syst = hh_dphi_dict_50_80["central value"].Clone("hh_central_dphi_50_80")


    for bin in range(1, 17):
        central_dphi_0_20_just_syst.SetBinError(bin, total_dphi_sys_0_20*central_dphi_0_20_just_syst.GetBinContent(bin))
        central_dphi_20_50_just_syst.SetBinError(bin, total_dphi_sys_20_50*central_dphi_20_50_just_syst.GetBinContent(bin))
        central_dphi_50_80_just_syst.SetBinError(bin, total_dphi_sys_50_80*central_dphi_50_80_just_syst.GetBinContent(bin))

        hh_central_dphi_0_20_just_syst.SetBinError(bin, hh_total_dphi_sys_0_20*hh_central_dphi_0_20_just_syst.GetBinContent(bin))
        hh_central_dphi_20_50_just_syst.SetBinError(bin, hh_total_dphi_sys_20_50*hh_central_dphi_20_50_just_syst.GetBinContent(bin))
        hh_central_dphi_50_80_just_syst.SetBinError(bin, hh_total_dphi_sys_50_80*hh_central_dphi_50_80_just_syst.GetBinContent(bin))



    for ONLY_SYST in [True, False]:
        dphi_dict_0_20["central value"].SetLineColor(rt.kRed)
        if ONLY_SYST:
            dphi_dict_0_20["central value"].Draw("E0")
        else:
            dphi_dict_0_20["central value"].Draw()
        central_dphi_0_20_just_syst.SetMarkerStyle(43)
        central_dphi_0_20_just_syst.SetMarkerSize(2)
        central_dphi_0_20_just_syst.SetMarkerColor(rt.kRed)
        central_dphi_0_20_just_syst.SetLineColor(rt.kRed)
        central_dphi_0_20_just_syst.SetFillColor(rt.kRed)
        central_dphi_0_20_just_syst.SetFillStyle(3004)
        central_dphi_0_20_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(central_dphi_0_20_just_syst, "0-20% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(central_dphi_0_20_just_syst, "0-20% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_0_20_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_0_20_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_0_20_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_0_20_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_0_20_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_0_20_stats_and_syst_highpt.eps")

        dphi_dict_20_50["central value"].SetLineColor(rt.kOrange - 3)
        if ONLY_SYST:
            dphi_dict_20_50["central value"].Draw("E0")
        else:
            dphi_dict_20_50["central value"].Draw()
        central_dphi_20_50_just_syst.SetMarkerStyle(43)
        central_dphi_20_50_just_syst.SetMarkerSize(2)
        central_dphi_20_50_just_syst.SetMarkerColor(rt.kOrange - 3)
        central_dphi_20_50_just_syst.SetLineColor(rt.kOrange - 3)
        central_dphi_20_50_just_syst.SetFillColor(rt.kOrange - 3)
        central_dphi_20_50_just_syst.SetFillStyle(3004)
        central_dphi_20_50_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(central_dphi_20_50_just_syst, "20-50% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(central_dphi_20_50_just_syst, "20-50% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_20_50_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_20_50_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_20_50_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_20_50_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_20_50_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_20_50_stats_and_syst_highpt.eps")
        dphi_dict_50_80["central value"].SetLineColor(rt.kBlue - 2)
        if ONLY_SYST:
            dphi_dict_50_80["central value"].Draw("E0")
        else:
            dphi_dict_50_80["central value"].Draw()
        central_dphi_50_80_just_syst.SetMarkerStyle(43)
        central_dphi_50_80_just_syst.SetMarkerSize(2)
        central_dphi_50_80_just_syst.SetMarkerColor(rt.kBlue - 2)
        central_dphi_50_80_just_syst.SetLineColor(rt.kBlue - 2)
        central_dphi_50_80_just_syst.SetFillColor(rt.kBlue - 2)
        central_dphi_50_80_just_syst.SetFillStyle(3004)
        central_dphi_50_80_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(central_dphi_50_80_just_syst, "50-80% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(central_dphi_50_80_just_syst, "50-80% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_50_80_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_50_80_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_50_80_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_lambda_dphi_50_80_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_lambda_dphi_50_80_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_lambda_dphi_50_80_stats_and_syst_highpt.eps")

        hh_dphi_dict_0_20["central value"].SetLineColor(rt.kRed + 2)
        if ONLY_SYST:
            hh_dphi_dict_0_20["central value"].Draw("E0")
        else:
            hh_dphi_dict_0_20["central value"].Draw()
        hh_central_dphi_0_20_just_syst.SetMarkerStyle(43)
        hh_central_dphi_0_20_just_syst.SetMarkerSize(2)
        hh_central_dphi_0_20_just_syst.SetMarkerColor(rt.kRed + 2) 
        hh_central_dphi_0_20_just_syst.SetLineColor(rt.kRed + 2)
        hh_central_dphi_0_20_just_syst.SetFillColor(rt.kRed + 2)
        hh_central_dphi_0_20_just_syst.SetFillStyle(3004)
        hh_central_dphi_0_20_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(hh_central_dphi_0_20_just_syst, "0-20% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(hh_central_dphi_0_20_just_syst, "0-20% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_0_20_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_0_20_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_0_20_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_0_20_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_0_20_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_0_20_stats_and_syst_highpt.eps")
        hh_dphi_dict_20_50["central value"].SetLineColor(rt.kOrange - 3)
        if ONLY_SYST:
            hh_dphi_dict_20_50["central value"].Draw("E0")
        else:
            hh_dphi_dict_20_50["central value"].Draw()
        hh_central_dphi_20_50_just_syst.SetMarkerStyle(43)
        hh_central_dphi_20_50_just_syst.SetMarkerSize(2)
        hh_central_dphi_20_50_just_syst.SetMarkerColor(rt.kOrange + 3)
        hh_central_dphi_20_50_just_syst.SetLineColor(rt.kOrange + 4)
        hh_central_dphi_20_50_just_syst.SetFillColor(rt.kOrange + 3)
        hh_central_dphi_20_50_just_syst.SetFillStyle(3004)
        hh_central_dphi_20_50_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(hh_central_dphi_20_50_just_syst, "20-50% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(hh_central_dphi_20_50_just_syst, "20-50% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_20_50_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_20_50_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_20_50_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_20_50_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_20_50_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_20_50_stats_and_syst_highpt.eps")
        hh_dphi_dict_50_80["central value"].SetLineColor(rt.kBlue - 2)
        if ONLY_SYST:
            hh_dphi_dict_50_80["central value"].Draw("E0")
        else:
            hh_dphi_dict_50_80["central value"].Draw()
        hh_central_dphi_50_80_just_syst.SetMarkerStyle(43)
        hh_central_dphi_50_80_just_syst.SetMarkerSize(2)
        hh_central_dphi_50_80_just_syst.SetMarkerColor(rt.kBlue + 2)
        hh_central_dphi_50_80_just_syst.SetLineColor(rt.kBlue + 2)
        hh_central_dphi_50_80_just_syst.SetFillColor(rt.kBlue + 2)
        hh_central_dphi_50_80_just_syst.SetFillStyle(3004)
        hh_central_dphi_50_80_just_syst.Draw("SAME E2")
        leg = rt.TLegend(0.6, 0.8, 0.9, 0.9)
        if ONLY_SYST:
            leg.AddEntry(hh_central_dphi_50_80_just_syst, "50-80% data (only syst. errors)", "lf")
        else:
            leg.AddEntry(hh_central_dphi_50_80_just_syst, "50-80% data (final)", "lepf")
        leg.Draw("same")
        c.Draw()
        if ONLY_SYST:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_50_80_onlysyst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_50_80_onlysyst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_50_80_onlysyst_highpt.eps")
        else:
            if PT_MODE == 0:
                c.SaveAs("figures/h_h_dphi_50_80_stats_and_syst.eps")
            elif PT_MODE == 1:
                c.SaveAs("figures/h_h_dphi_50_80_stats_and_syst_lowpt.eps")
            elif PT_MODE == 2:
                c.SaveAs("figures/h_h_dphi_50_80_stats_and_syst_highpt.eps")



    v0_dphi_dist_0_20 = dphi_dict_0_20["central value"].Clone()
    v0_dphi_dist_20_50 = dphi_dict_20_50["central value"].Clone()
    v0_dphi_dist_50_80 = dphi_dict_50_80["central value"].Clone()

    v0_ue_line_0_20 = file_dict["central value"].Get("ue_line_0_20")
    v0_ue_line_20_50 = file_dict["central value"].Get("ue_line_20_50")
    v0_ue_line_50_80 = file_dict["central value"].Get("ue_line_50_80")

    v0_ue_line_0_20.SetParameter(0, v0_ue_line_0_20.GetParameter(0)/bin_width)
    v0_ue_line_20_50.SetParameter(0, v0_ue_line_20_50.GetParameter(0)/bin_width)
    v0_ue_line_50_80.SetParameter(0, v0_ue_line_50_80.GetParameter(0)/bin_width)

    v0_ue_line_0_20.SetLineColor(rt.kPink - 2)
    v0_ue_line_20_50.SetLineColor(rt.kOrange - 4)
    v0_ue_line_50_80.SetLineColor(rt.kBlue - 4)

    central_dphi_0_20_just_syst.SetFillColor(rt.kPink - 2)
    central_dphi_20_50_just_syst.SetFillColor(rt.kOrange - 4)
    central_dphi_50_80_just_syst.SetFillColor(rt.kBlue - 4)

    central_dphi_0_20_just_syst.SetFillStyle(3244)
    central_dphi_20_50_just_syst.SetFillStyle(3244)
    central_dphi_50_80_just_syst.SetFillStyle(3244)

    central_dphi_0_20_just_syst.SetFillColorAlpha(rt.kPink - 2, 0.5)
    central_dphi_20_50_just_syst.SetFillColorAlpha(rt.kOrange - 4, 0.5)
    central_dphi_50_80_just_syst.SetFillColorAlpha(rt.kBlue - 4, 0.5)

    central_dphi_0_20_just_syst.SetMarkerSize(0)
    central_dphi_20_50_just_syst.SetMarkerSize(0)
    central_dphi_50_80_just_syst.SetMarkerSize(0)

    v0_ue_line_0_20.SetLineStyle(3)
    v0_ue_line_20_50.SetLineStyle(3)
    v0_ue_line_50_80.SetLineStyle(3)

    v0_ue_line_0_20.SetLineWidth(1)
    v0_ue_line_20_50.SetLineWidth(1)
    v0_ue_line_50_80.SetLineWidth(1)

    v0_dphi_dist_0_20.SetLineColor(rt.kPink - 1)
    v0_dphi_dist_20_50.SetLineColor(rt.kOrange - 3)
    v0_dphi_dist_50_80.SetLineColor(rt.kBlue - 3)


    v0_dphi_dist_0_20.SetMarkerStyle(43)
    v0_dphi_dist_20_50.SetMarkerStyle(43)
    v0_dphi_dist_50_80.SetMarkerStyle(43)

    v0_dphi_dist_0_20.SetMarkerColor(rt.kPink - 1)
    v0_dphi_dist_20_50.SetMarkerColor(rt.kOrange - 3)
    v0_dphi_dist_50_80.SetMarkerColor(rt.kBlue - 3)


    v0_dphi_dist_0_20.SetMarkerSize(1.5)
    v0_dphi_dist_20_50.SetMarkerSize(1.5)
    v0_dphi_dist_50_80.SetMarkerSize(1.5)


    central_dphi_0_20_just_syst.SetTitle("")
    central_dphi_0_20_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")
    central_dphi_0_20_just_syst.GetXaxis().CenterTitle()
    central_dphi_0_20_just_syst.GetYaxis().SetTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{pair}}{d#Delta#it{#varphi}}")


    central_dphi_20_50_just_syst.SetTitle("")
    central_dphi_20_50_just_syst.GetXaxis().CenterTitle()
    central_dphi_20_50_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")

    central_dphi_50_80_just_syst.SetTitle("")
    central_dphi_50_80_just_syst.GetXaxis().CenterTitle()
    central_dphi_50_80_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")

    if PT_MODE == 0:
        continue

    TOP_MARGIN = 0.1
    BOTTOM_MARGIN = 0.1
    LEFT_MARGIN = 0.25
    FIRST_PAD_WIDTH = 0.375
    PAD_WIDTH = 0.3

    PAD_HEIGHT = 0.5
    BUFFER = 0.087
    SCALE = PAD_WIDTH/(PAD_WIDTH + BUFFER)

    EPSILON = 0.0001

    X_LABEL_SIZE = 0.06
    Y_LABEL_SIZE = 0.06
    X_TITLE_SIZE = 0.065
    Y_TITLE_SIZE = 0.065
    X_TITLE_OFFSET = 0.63
    Y_TITLE_OFFSET = 2.1

    MAX_Y_RANGE = v0_dphi_dist_0_20.GetMaximum()
    MIN_Y_RANGE = -0.01


    # dphi_all = rt.TCanvas(f"dphi_all_{PT_MODE}", f"dphi_all_{PT_MODE}", 50, 50, 1200, 525)
    # dphi_all.SetMargin(0, 0, 0, 0)

    # three pads, each completely next to each other
    # but the first pad needs a left margin, but the actual pads need to be the same size
    # so the first pad needs to be a little bigger

    dphi_all_for_real.cd()
    hl020pad = rt.TPad("hl020pad", "", 0, PAD_HEIGHT, FIRST_PAD_WIDTH, 1.0)
    hl020pad.SetMargin(LEFT_MARGIN, 0.0, EPSILON, TOP_MARGIN)
    hl020pad.SetTicks(1, 1)
    hl020pad.Draw()
    hl020pad.cd()

    central_dphi_0_20_just_syst.GetYaxis().SetMaxDigits(3)
    central_dphi_0_20_just_syst.GetYaxis().SetTitleSize(SCALE*Y_TITLE_SIZE)
    central_dphi_0_20_just_syst.GetYaxis().SetLabelSize(SCALE*Y_LABEL_SIZE)
    central_dphi_0_20_just_syst.GetYaxis().SetTitleOffset(Y_TITLE_OFFSET)
    central_dphi_0_20_just_syst.GetYaxis().SetRangeUser(0, MAX_Y_RANGE)
    central_dphi_0_20_just_syst.GetXaxis().SetLabelSize(SCALE*X_LABEL_SIZE)
    central_dphi_0_20_just_syst.GetXaxis().SetTitleSize(SCALE*X_TITLE_SIZE*1.06)
    central_dphi_0_20_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET + 0.15)

    central_dphi_0_20_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    central_dphi_0_20_just_syst.Draw("E2")
    v0_dphi_dist_0_20.Draw("SAME")
    v0_ue_line_0_20.Draw("SAME")


    label_x_start = 0.37
    label_y_start = 0.30
    label_text_space = 0.07
    alice_data_label = rt.TLatex()
    alice_data_label.SetNDC()
    alice_data_label.SetTextSize(0.06)
    alice_data_label.SetTextAlign(13)
    # alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
    alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{ALICE p#minusPb}")
    alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{#sqrt{#it{s}_{NN}} = 5.02 TeV}")



    dphi_all_for_real.cd()

    hl2050pad = rt.TPad("hl2050pad", "", FIRST_PAD_WIDTH, PAD_HEIGHT,  FIRST_PAD_WIDTH + PAD_WIDTH, 1.0)
    hl2050pad.SetMargin(0.0, 0.0, EPSILON, TOP_MARGIN)
    hl2050pad.SetTicks(1, 1)
    hl2050pad.Draw()
    hl2050pad.cd()


    central_dphi_20_50_just_syst.GetXaxis().SetLabelOffset(-0.002)
    central_dphi_20_50_just_syst.GetXaxis().SetTitleOffset(0.83)
    central_dphi_20_50_just_syst.GetXaxis().SetLabelSize(X_LABEL_SIZE)
    central_dphi_20_50_just_syst.GetXaxis().SetTitleSize(X_TITLE_SIZE)
    central_dphi_20_50_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET)
    central_dphi_20_50_just_syst.GetYaxis().SetLabelSize(0.0)
    central_dphi_20_50_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    v0_dphi_dist_20_50.GetYaxis().SetLabelSize(0.0)
    v0_dphi_dist_20_50.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    central_dphi_20_50_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)
    central_dphi_20_50_just_syst.Draw("E2")
    v0_dphi_dist_20_50.Draw("SAME")
    v0_ue_line_20_50.Draw("SAME")


    pt_label_x_start = 0.16
    pt_label_y_start = 0.11
    pt_label_text_space = 0.07

    pt_data_label = rt.TLatex()
    pt_data_label.SetNDC()
    pt_data_label.SetTextSize(0.05/SCALE)
    pt_data_label.SetTextAlign(13)
    pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start, "#bf{4.0 <   #it{p}^{h}_{T,trig}    < 8.0 GeV/#it{c}}")
    # if PT_MODE == 0:
    #     pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{2.0 <   #it{p}^{#Lambda, h}_{T,assoc} < 4.0 GeV/#it{c}}")
    # elif PT_MODE == 1:
    #     pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{1.5 <   #it{p}^{#Lambda, h}_{T,assoc} < 2.5 GeV/#it{c}}")
    # elif PT_MODE == 2:
    #     pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{2.5 <   #it{p}^{#Lambda, h}_{T,assoc} < 4.0 GeV/#it{c}}")



    dphi_all_for_real.cd()

    hl5080pad = rt.TPad("hl5080pad", "", FIRST_PAD_WIDTH + PAD_WIDTH, PAD_HEIGHT, FIRST_PAD_WIDTH + 2*PAD_WIDTH, 1.0)
    hl5080pad.SetMargin(0.0, EPSILON, EPSILON, TOP_MARGIN)
    hl5080pad.SetTicks(1, 1)
    hl5080pad.Draw()
    hl5080pad.cd()


    central_dphi_50_80_just_syst.GetXaxis().SetLabelOffset(-0.002)
    central_dphi_50_80_just_syst.GetXaxis().SetTitleOffset(0.83)
    central_dphi_50_80_just_syst.GetXaxis().SetLabelSize(X_LABEL_SIZE)
    central_dphi_50_80_just_syst.GetXaxis().SetTitleSize(X_TITLE_SIZE)
    central_dphi_50_80_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET)
    central_dphi_50_80_just_syst.GetYaxis().SetLabelSize(0.0)
    central_dphi_50_80_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    central_dphi_50_80_just_syst.Draw("E2")
    v0_dphi_dist_50_80.Draw("SAME")
    v0_ue_line_50_80.Draw("SAME")

    color_legend = rt.TLegend(0.06, 0.65, 0.30, 0.86)
    color_legend.SetTextSize(0.05/SCALE)
    color_legend.SetBorderSize(0)
    color_legend.SetFillStyle(0)
    legend_box_0_20 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    legend_box_0_20.SetFillStyle(1001)
    legend_box_0_20.SetFillColor(rt.kPink - 1)
    legend_box_0_20.SetLineWidth(0)
    legend_box_20_50 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    legend_box_20_50.SetFillStyle(1001)
    legend_box_20_50.SetFillColor(rt.kOrange - 3)
    legend_box_20_50.SetLineWidth(0)
    legend_box_50_80 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    legend_box_50_80.SetFillStyle(1001)
    legend_box_50_80.SetFillColor(rt.kBlue - 3)
    legend_box_50_80.SetLineWidth(0)
    color_legend.AddEntry(legend_box_0_20, "V0A 0#minus20%", "f")
    color_legend.AddEntry(legend_box_20_50, "V0A 20#minus50%", "f")
    color_legend.AddEntry(legend_box_50_80, "V0A 50#minus80%", "f")
    color_legend.Draw()


    marker_legend = rt.TLegend(0.5, 0.67, 0.93, 0.87)
    marker_legend.SetTextSize(0.07/SCALE)
    marker_legend.SetBorderSize(0)
    marker_legend.SetFillStyle(0)

    legend_marker = rt.TMarker(0.55, 0.65, 20)
    legend_marker.SetMarkerColor(rt.kBlack)
    legend_marker.SetMarkerStyle(43)
    legend_marker.SetMarkerSize(3)

    legend_line = rt.TLine(0.7, 0.65, 0.95, 0.65)
    legend_line.SetLineColor(rt.kBlack)
    legend_line.SetLineWidth(2)
    legend_line.SetLineStyle(3)

    marker_legend.AddEntry(legend_marker, " h#minus#Lambda data", "lp")
    marker_legend.AddEntry(legend_line, " UE fit", "l")

    marker_legend.Draw()

    # quit()

    # dphi_all.Draw()

    # quit()

    # if PT_MODE == 0:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all.eps")
    # elif PT_MODE == 1:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all_lowpt.eps")
    # elif PT_MODE == 2:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all_highpt.eps")

    # if PT_MODE == 0:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all.pdf")
    # elif PT_MODE == 1:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all_lowpt.pdf")
    # elif PT_MODE == 2:
    #     dphi_all.SaveAs("figures/h_lambda_dphi_all_highpt.pdf")

    # dphi_all.Delete()

    hh_dphi_dist_0_20 = hh_dphi_dict_0_20["central value"].Clone()
    hh_dphi_dist_20_50 = hh_dphi_dict_20_50["central value"].Clone()
    hh_dphi_dist_50_80 = hh_dphi_dict_50_80["central value"].Clone()

    hh_ue_line_0_20 = file_dict["central value"].Get("hh_ue_line_0_20")
    hh_ue_line_20_50 = file_dict["central value"].Get("hh_ue_line_20_50")
    hh_ue_line_50_80 = file_dict["central value"].Get("hh_ue_line_50_80")

    hh_ue_line_0_20.SetParameter(0, hh_ue_line_0_20.GetParameter(0)/hh_bin_width)
    hh_ue_line_20_50.SetParameter(0, hh_ue_line_20_50.GetParameter(0)/hh_bin_width)
    hh_ue_line_50_80.SetParameter(0, hh_ue_line_50_80.GetParameter(0)/hh_bin_width)

    hh_ue_line_0_20.SetLineColor(rt.kPink - 6)
    hh_ue_line_20_50.SetLineColor(rt.kOrange - 6)
    hh_ue_line_50_80.SetLineColor(rt.kBlue - 6)

    hh_ue_line_0_20.SetLineStyle(3)
    hh_ue_line_20_50.SetLineStyle(3)
    hh_ue_line_50_80.SetLineStyle(3)

    hh_ue_line_0_20.SetLineWidth(1)
    hh_ue_line_20_50.SetLineWidth(1)
    hh_ue_line_50_80.SetLineWidth(1)

    hh_dphi_dist_0_20.SetLineColor(rt.kPink - 6)
    hh_dphi_dist_20_50.SetLineColor(rt.kOrange - 6)
    hh_dphi_dist_50_80.SetLineColor(rt.kBlue - 6)

    hh_central_dphi_0_20_just_syst.SetFillColor(rt.kPink - 6)
    hh_central_dphi_20_50_just_syst.SetFillColor(rt.kOrange - 6)
    hh_central_dphi_50_80_just_syst.SetFillColor(rt.kBlue - 6)

    hh_central_dphi_0_20_just_syst.SetMarkerSize(0)
    hh_central_dphi_20_50_just_syst.SetMarkerSize(0)
    hh_central_dphi_50_80_just_syst.SetMarkerSize(0)

    hh_central_dphi_0_20_just_syst.SetFillStyle(3244)
    hh_central_dphi_20_50_just_syst.SetFillStyle(3244)
    hh_central_dphi_50_80_just_syst.SetFillStyle(3244)

    hh_central_dphi_0_20_just_syst.SetFillColorAlpha(rt.kPink - 6, 0.5)
    hh_central_dphi_20_50_just_syst.SetFillColorAlpha(rt.kOrange - 6, 0.5)
    hh_central_dphi_50_80_just_syst.SetFillColorAlpha(rt.kBlue - 6, 0.5)

    hh_dphi_dist_0_20.SetMarkerStyle(34)
    hh_dphi_dist_20_50.SetMarkerStyle(34)
    hh_dphi_dist_50_80.SetMarkerStyle(34)

    hh_dphi_dist_0_20.SetMarkerColor(rt.kPink - 6)
    hh_dphi_dist_20_50.SetMarkerColor(rt.kOrange - 6)
    hh_dphi_dist_50_80.SetMarkerColor(rt.kBlue - 6)

    hh_dphi_dist_0_20.SetMarkerSize(1.5)
    hh_dphi_dist_20_50.SetMarkerSize(1.5)
    hh_dphi_dist_50_80.SetMarkerSize(1.5)

    hh_central_dphi_0_20_just_syst.SetTitle("")
    hh_central_dphi_0_20_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")
    hh_central_dphi_0_20_just_syst.GetXaxis().CenterTitle()
    hh_central_dphi_0_20_just_syst.GetYaxis().SetTitle("#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{pair}}{d#Delta#it{#varphi}}")


    hh_central_dphi_20_50_just_syst.SetTitle("")
    hh_central_dphi_20_50_just_syst.GetXaxis().CenterTitle()
    hh_central_dphi_20_50_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")

    hh_central_dphi_50_80_just_syst.SetTitle("")
    hh_central_dphi_50_80_just_syst.GetXaxis().CenterTitle()
    hh_central_dphi_50_80_just_syst.GetXaxis().SetTitle("#Delta#it{#varphi}")

    # TOP_MARGIN = 0.1
    # BOTTOM_MARGIN = 0.1
    # LEFT_MARGIN = 0.24
    # PAD_WIDTH = 0.3
    # BUFFER = 0.087
    # SCALE = PAD_WIDTH/(PAD_WIDTH + BUFFER)

    X_LABEL_SIZE = 0.06
    Y_LABEL_SIZE = 0.06
    X_TITLE_SIZE = 0.065
    Y_TITLE_SIZE = 0.065
    X_TITLE_OFFSET = 0.63
    Y_TITLE_OFFSET = 2.1

    MAX_Y_RANGE = hh_dphi_dist_0_20.GetMaximum()
    MIN_Y_RANGE = -0.05

    # dphi_all = rt.TCanvas(f"hh_dphi_all_{PT_MODE}", f"hh_dphi_all_{PT_MODE}", 50, 50, 1200, 525)
    # dphi_all.SetMargin(0, 0, 0, 0)
    dphi_all_for_real.cd()

    hl020pad = rt.TPad("hl020pad", "", 0, 0, FIRST_PAD_WIDTH, PAD_HEIGHT)
    hl020pad.SetMargin(LEFT_MARGIN, 0.0, BOTTOM_MARGIN, 0)
    hl020pad.SetTicks(1, 1)
    hl020pad.Draw()
    hl020pad.cd()

    hh_central_dphi_0_20_just_syst.GetYaxis().SetMaxDigits(3)
    hh_central_dphi_0_20_just_syst.GetYaxis().SetTitleSize(SCALE*Y_TITLE_SIZE)
    hh_central_dphi_0_20_just_syst.GetYaxis().SetLabelSize(SCALE*Y_LABEL_SIZE)
    hh_central_dphi_0_20_just_syst.GetYaxis().SetTitleOffset(Y_TITLE_OFFSET)
    hh_central_dphi_0_20_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)
    hh_central_dphi_0_20_just_syst.GetXaxis().SetLabelSize(SCALE*X_LABEL_SIZE)
    hh_central_dphi_0_20_just_syst.GetXaxis().SetTitleSize(SCALE*X_TITLE_SIZE*1.06)
    hh_central_dphi_0_20_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET + 0.15)

    hh_central_dphi_0_20_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)
    hh_central_dphi_0_20_just_syst.Draw("E2")
    hh_dphi_dist_0_20.Draw("SAME")
    hh_ue_line_0_20.Draw("SAME")

    label_x_start = 0.29
    label_y_start = 0.38
    label_text_space = 0.07

    # alice_data_label = rt.TLatex()
    # alice_data_label.SetNDC()
    # alice_data_label.SetTextSize(0.06)
    # alice_data_label.SetTextAlign(13)
    # alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
    # alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}")
    # alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{|#Delta#it{#eta}| < 1.2}")

    dphi_all_for_real.cd()

    hl2050pad = rt.TPad("hl2050pad", "", FIRST_PAD_WIDTH, 0, FIRST_PAD_WIDTH + PAD_WIDTH, PAD_HEIGHT)
    hl2050pad.SetMargin(0.0, 0.0, BOTTOM_MARGIN, 0)
    hl2050pad.SetTicks(1, 1)
    hl2050pad.Draw()
    hl2050pad.cd()


    hh_central_dphi_20_50_just_syst.GetXaxis().SetLabelOffset(-0.002)
    hh_central_dphi_20_50_just_syst.GetXaxis().SetTitleOffset(0.83)
    hh_central_dphi_20_50_just_syst.GetXaxis().SetLabelSize(X_LABEL_SIZE)
    hh_central_dphi_20_50_just_syst.GetXaxis().SetTitleSize(X_TITLE_SIZE)
    hh_central_dphi_20_50_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET)
    hh_central_dphi_20_50_just_syst.GetYaxis().SetLabelSize(0.0)
    hh_central_dphi_20_50_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    hh_dphi_dist_20_50.GetYaxis().SetLabelSize(0.0)
    hh_dphi_dist_20_50.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    hh_central_dphi_20_50_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)
    hh_central_dphi_20_50_just_syst.Draw("E2")
    hh_dphi_dist_20_50.Draw("SAME")
    hh_ue_line_20_50.Draw("SAME")

    pt_label_x_start = 0.16
    pt_label_y_start = 1.02
    pt_label_text_space = 0.073

    pt_data_label = rt.TLatex()
    pt_data_label.SetNDC()
    pt_data_label.SetTextSize(0.05/SCALE)
    pt_data_label.SetTextAlign(13)
    if PT_MODE == 0:
        pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{2.0 <   #it{p}^{#Lambda, h}_{T,assoc} < 4.0 GeV/#it{c}}")
    elif PT_MODE == 1:
        pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{1.5 <   #it{p}^{#Lambda, h}_{T,assoc} < 2.5 GeV/#it{c}}")
    elif PT_MODE == 2:
        pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - pt_label_text_space, "#bf{2.5 <   #it{p}^{#Lambda, h}_{T,assoc} < 4.0 GeV/#it{c}}")
    pt_data_label.DrawLatex(pt_label_x_start, pt_label_y_start - 2*pt_label_text_space - 0.02, "#bf{|#Delta#it{#eta}| < 1.2}")



    dphi_all_for_real.cd()

    hl5080pad = rt.TPad("hl5080pad", "", FIRST_PAD_WIDTH + PAD_WIDTH, 0, 2*PAD_WIDTH + FIRST_PAD_WIDTH, PAD_HEIGHT)
    hl5080pad.SetMargin(0.0, 0.0, BOTTOM_MARGIN, 0)
    hl5080pad.SetTicks(1, 1)
    hl5080pad.Draw()
    hl5080pad.cd()


    hh_central_dphi_50_80_just_syst.GetXaxis().SetLabelOffset(-0.002)
    hh_central_dphi_50_80_just_syst.GetXaxis().SetTitleOffset(0.83)
    hh_central_dphi_50_80_just_syst.GetXaxis().SetLabelSize(X_LABEL_SIZE)
    hh_central_dphi_50_80_just_syst.GetXaxis().SetTitleSize(X_TITLE_SIZE)
    hh_central_dphi_50_80_just_syst.GetXaxis().SetTitleOffset(X_TITLE_OFFSET)
    hh_central_dphi_50_80_just_syst.GetYaxis().SetLabelSize(0.0)
    hh_central_dphi_50_80_just_syst.GetYaxis().SetRangeUser(MIN_Y_RANGE, MAX_Y_RANGE)

    hh_central_dphi_50_80_just_syst.Draw("E2")
    hh_dphi_dist_50_80.Draw("SAME")
    hh_ue_line_50_80.Draw("SAME")

    color_hh_legend = rt.TLegend(0.06, 0.75, 0.3, 0.95)
    color_hh_legend.SetTextSize(0.05/SCALE)
    color_hh_legend.SetBorderSize(0)
    color_hh_legend.SetFillStyle(0)
    hh_legend_box_0_20 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    hh_legend_box_0_20.SetFillStyle(1001)
    hh_legend_box_0_20.SetFillColor(rt.kPink - 6)
    hh_legend_box_0_20.SetLineWidth(0)
    hh_legend_box_20_50 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    hh_legend_box_20_50.SetFillStyle(1001)
    hh_legend_box_20_50.SetFillColor(rt.kOrange - 6)
    hh_legend_box_20_50.SetLineWidth(0)
    hh_legend_box_50_80 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    hh_legend_box_50_80.SetFillStyle(1001)
    hh_legend_box_50_80.SetFillColor(rt.kBlue - 6)
    hh_legend_box_50_80.SetLineWidth(0)
    color_hh_legend.AddEntry(hh_legend_box_0_20, "V0A 0#minus20%", "f")
    color_hh_legend.AddEntry(hh_legend_box_20_50, "V0A 20#minus50%", "f")
    color_hh_legend.AddEntry(hh_legend_box_50_80, "V0A 50#minus80%", "f")
    color_hh_legend.Draw()

    # color_legend = rt.TLegend(0.15, 0.65, 0.47, 0.86)
    # color_legend.SetTextSize(0.05/SCALE)
    # color_legend.SetBorderSize(0)
    # color_legend.SetFillStyle(0)

    # legend_box_0_20 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    # legend_box_0_20.SetFillStyle(1001)
    # legend_box_0_20.SetFillColor(rt.kPink - 6)
    # legend_box_0_20.SetLineWidth(0)

    # legend_box_20_50 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    # legend_box_20_50.SetFillStyle(1001)
    # legend_box_20_50.SetFillColor(rt.kOrange - 6)
    # legend_box_20_50.SetLineWidth(0)

    # legend_box_50_80 = rt.TBox(0.55, 0.65, 0.85, 0.86)
    # legend_box_50_80.SetFillStyle(1001)
    # legend_box_50_80.SetFillColor(rt.kBlue - 6)
    # legend_box_50_80.SetLineWidth(0)

    # color_legend.AddEntry(legend_box_0_20, "V0A 0#minus20%", "f")
    # color_legend.AddEntry(legend_box_20_50, "V0A 20#minus50%", "f")
    # color_legend.AddEntry(legend_box_50_80, "V0A 50#minus80%", "f")

    # color_legend.Draw()

    hh_marker_legend = rt.TLegend(0.5, 0.75, 0.93, 0.95)
    hh_marker_legend.SetTextSize(0.07/SCALE)
    hh_marker_legend.SetBorderSize(0)
    hh_marker_legend.SetFillStyle(0)

    legend_hh_marker = rt.TMarker(0.55, 0.65, 20)
    legend_hh_marker.SetMarkerColor(rt.kBlack)
    legend_hh_marker.SetMarkerStyle(34)
    legend_hh_marker.SetMarkerSize(2)

    legend_hh_line = rt.TLine(0.7, 0.65, 0.95, 0.65)
    legend_hh_line.SetLineColor(rt.kBlack)
    legend_hh_line.SetLineWidth(2)
    legend_hh_line.SetLineStyle(3)

    hh_marker_legend.AddEntry(legend_hh_marker, " h#minush data", "lp")
    hh_marker_legend.AddEntry(legend_hh_line, " UE fit", "l")

    hh_marker_legend.Draw()


    if PT_MODE == 1:
        dphi_all_for_real.SaveAs("figures/dphi_final_lowpt.eps")
    elif PT_MODE == 2:
        dphi_all_for_real.SaveAs("figures/dphi_final_highpt.eps")

    if PT_MODE == 1:
        dphi_all_for_real.SaveAs("figures/dphi_final_lowpt.pdf")
    elif PT_MODE == 2:
        dphi_all_for_real.SaveAs("figures/dphi_final_highpt.pdf")


    # dphi_all.Draw()

    # if PT_MODE == 0:
    #     dphi_all.SaveAs("figures/h_h_dphi_all.eps")
    # elif PT_MODE == 1:
    #     dphi_all.SaveAs("figures/h_h_dphi_all_lowpt.eps")
    # elif PT_MODE == 2:
    #     dphi_all.SaveAs("figures/h_h_dphi_all_highpt.eps")

    # if PT_MODE == 0:
    #     dphi_all.SaveAs("figures/h_h_dphi_all.pdf")
    # elif PT_MODE == 1:
    #     dphi_all.SaveAs("figures/h_h_dphi_all_lowpt.pdf")
    # elif PT_MODE == 2:
    #     dphi_all.SaveAs("figures/h_h_dphi_all_highpt.pdf")
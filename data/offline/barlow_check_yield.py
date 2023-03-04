import ROOT as rt

import math
import array as arr

from time import sleep


rt.gStyle.SetOptStat(0)

c = rt.TCanvas("c","c",800,600)
c.SetLeftMargin(0.10)
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)

colors = [rt.kGreen - 2, rt.kRed, rt.kBlue, rt.kOrange, rt.kMagenta, rt.kCyan, rt.kYellow, rt.kGray, rt.kViolet, rt.kTeal, rt.kSpring, rt.kAzure, rt.kPink, rt.kCopper, rt.kOrange+7, rt.kSpring+9, rt.kTeal-7, rt.kAzure+2, rt.kPink+10, rt.kCopper+3, rt.kOrange-3, rt.kSpring-5, rt.kTeal+3, rt.kAzure-6, rt.kPink-7, rt.kCopper-9]

pt_mode = 0 # 0 = central, 1 = low, 2 = high 

def do_barlow(yield_dict, name):

    color_index = 0

    central_yield = yield_dict["central value"]
    barlow_legend = rt.TLegend(0.15, 0.15, 0.32, 0.32)
    barlow_legend.SetBorderSize(0)
    for key in yield_dict:
        if key == "central value":
            base_hist = rt.TH1D("barlow_base","", 3, 0, 3)
            base_hist.GetXaxis().SetTitle("Kinematic Region")
            base_hist.GetXaxis().SetBinLabel(1, "Near Side")
            base_hist.GetXaxis().SetBinLabel(2, "Away Side")
            base_hist.GetXaxis().SetBinLabel(3, "Underlying Event")
            base_hist.GetYaxis().SetTitle("(Yield_{var}-Yield_{central})/#sqrt{|#sigma_{var}^{2} - #sigma_{central}^{2}|}")
            base_hist.GetYaxis().SetTitleOffset(1.3)
            base_hist.GetYaxis().SetRangeUser(-10, 10)
            base_hist.SetLineColor(rt.kBlack)
            base_hist.SetLineStyle(0)
            base_hist.SetLineWidth(0)
            base_hist.Draw()
        else:
            cur_yield = yield_dict[key]
            barlow_hist = rt.TH1D("barlow_"+key,"",3, 0, 3).Clone("barlow_hist_" + key)
            for i in range(3):
                # avg6 (non-negative) shares the same error as avg6 sometimes, as does gaus (for the UE, which is the same as our central)
                if (key == "avg6 (non-negative)" or key == "gaus") and cur_yield[i][1] == central_yield[i][1]:
                    cur_yield[i][1] = 0

                barlow_hist.SetBinContent(i+1, (cur_yield[i][0] - central_yield[i][0])/math.sqrt(abs(cur_yield[i][1]**2 - central_yield[i][1]**2)))
            barlow_legend.AddEntry(barlow_hist.Clone(), key, "P")
            barlow_hist.SetLineColor(colors[color_index])
            barlow_hist.SetMarkerColor(colors[color_index])
            barlow_hist.SetMarkerStyle(43)
            barlow_hist.SetMarkerSize(2)
            barlow_hist.DrawCopy("SAME P")
        color_index += 1

    barlow_min = rt.TLine(0, -1, 3, -1)
    barlow_max = rt.TLine(0, 1, 3, 1)
    barlow_min.SetLineColor(rt.kRed)
    barlow_max.SetLineColor(rt.kRed)
    barlow_min.SetLineStyle(rt.kDashed)
    barlow_max.SetLineStyle(rt.kDashed)
    barlow_min.SetLineWidth(2)
    barlow_max.SetLineWidth(2)
    barlow_min.Draw("SAME")
    barlow_max.Draw("SAME")
    barlow_legend.Draw("SAME")
    sleep(0.1)
    c.SaveAs(name)



if pt_mode == 0: # CENTRAL 2 - 4 BIN

    yield_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_zyam_file = rt.TFile.Open("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_avg4_file = rt.TFile.Open("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_avg6nonneg_file = rt.TFile.Open("output/yield_variation/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_gaus_file = rt.TFile.Open("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_v2_file = rt.TFile.Open("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    yield_von_file = rt.TFile.Open("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")


elif pt_mode == 1: # LOW 1.5 - 2.5 BIN

    yield_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_zyam_file = rt.TFile.Open("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_avg4_file = rt.TFile.Open("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_avg6nonneg_file = rt.TFile.Open("output/yield_variation/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_gaus_file = rt.TFile.Open("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_v2_file = rt.TFile.Open("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")
    yield_von_file = rt.TFile.Open("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_lowpt.root")


elif pt_mode == 2: # HIGH 2.5 - 4 BIN

    yield_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_zyam_file = rt.TFile.Open("output/yield_variation/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_avg4_file = rt.TFile.Open("output/yield_variation/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_avg6nonneg_file = rt.TFile.Open("output/yield_variation/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_gaus_file = rt.TFile.Open("output/yield_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_v2_file = rt.TFile.Open("output/yield_variation/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")
    yield_von_file = rt.TFile.Open("output/yield_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal_highpt.root")

signal_file_dict = {
    "central value": yield_central_file,
    "zyam": yield_zyam_file,
    "avg4": yield_avg4_file,
    "avg6 (non-negative)": yield_avg6nonneg_file,
    "v2": yield_v2_file,
    "von": yield_von_file,
    "gaus": yield_gaus_file
}


h_lambda_yield_dict_0_20 = {}
h_lambda_yield_dict_20_50 = {}
h_lambda_yield_dict_50_80 = {}

h_h_yield_dict_0_20 = {}
h_h_yield_dict_20_50 = {}
h_h_yield_dict_50_80 = {}

for key, file in signal_file_dict.items():
    # total is obviously not considered, as it is mostly the same for all variations
    h_lambda_near_graph = file.Get("h_lambda_near_graph")
    h_lambda_away_graph = file.Get("h_lambda_away_graph")
    h_lambda_ue_graph = file.Get("h_lambda_ue_graph")

    h_h_near_graph = file.Get("h_h_near_graph")
    h_h_away_graph = file.Get("h_h_away_graph")
    h_h_ue_graph = file.Get("h_h_ue_graph")

    # ARRAY INDEXING: [0] = near, [1] = away, [2] = ue
    # SUBARRAY INDEXING: [0] = yield, [1] = error
    h_lambda_yield_dict_0_20[key] = [[h_lambda_near_graph.GetY()[2], h_lambda_near_graph.GetEY()[2]], [h_lambda_away_graph.GetY()[2], h_lambda_away_graph.GetEY()[2]], [h_lambda_ue_graph.GetY()[2], h_lambda_ue_graph.GetEY()[2]]]
    h_lambda_yield_dict_20_50[key] = [[h_lambda_near_graph.GetY()[1], h_lambda_near_graph.GetEY()[1]], [h_lambda_away_graph.GetY()[1], h_lambda_away_graph.GetEY()[1]], [h_lambda_ue_graph.GetY()[1], h_lambda_ue_graph.GetEY()[1]]]
    h_lambda_yield_dict_50_80[key] = [[h_lambda_near_graph.GetY()[0], h_lambda_near_graph.GetEY()[0]], [h_lambda_away_graph.GetY()[0], h_lambda_away_graph.GetEY()[0]], [h_lambda_ue_graph.GetY()[0], h_lambda_ue_graph.GetEY()[0]]]

    h_h_yield_dict_0_20[key] = [[h_h_near_graph.GetY()[2], h_h_near_graph.GetEY()[2]], [h_h_away_graph.GetY()[2], h_h_away_graph.GetEY()[2]], [h_h_ue_graph.GetY()[2], h_h_ue_graph.GetEY()[2]]]
    h_h_yield_dict_20_50[key] = [[h_h_near_graph.GetY()[1], h_h_near_graph.GetEY()[1]], [h_h_away_graph.GetY()[1], h_h_away_graph.GetEY()[1]], [h_h_ue_graph.GetY()[1], h_h_ue_graph.GetEY()[1]]]
    h_h_yield_dict_50_80[key] = [[h_h_near_graph.GetY()[0], h_h_near_graph.GetEY()[0]], [h_h_away_graph.GetY()[0], h_h_away_graph.GetEY()[0]], [h_h_ue_graph.GetY()[0], h_h_ue_graph.GetEY()[0]]]


if pt_mode == 0:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80.pdf")
elif pt_mode == 1:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20_lowpt.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50_lowpt.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80_lowpt.pdf")
elif pt_mode == 2:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20_highpt.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50_highpt.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80_highpt.pdf")
if pt_mode == 0:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80.pdf")
elif pt_mode == 1:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20_lowpt.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50_lowpt.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80_lowpt.pdf")
elif pt_mode == 2:
    do_barlow(h_lambda_yield_dict_0_20, "figures/h_lambda_yield_barlow_0_20_highpt.pdf")
    do_barlow(h_lambda_yield_dict_20_50, "figures/h_lambda_yield_barlow_20_50_highpt.pdf")
    do_barlow(h_lambda_yield_dict_50_80, "figures/h_lambda_yield_barlow_50_80_highpt.pdf")

if pt_mode == 0:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80.pdf")
elif pt_mode == 1:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20_lowpt.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50_lowpt.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80_lowpt.pdf")
elif pt_mode == 2:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20_highpt.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50_highpt.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80_highpt.pdf")
if pt_mode == 0:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80.pdf")
elif pt_mode == 1:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20_lowpt.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50_lowpt.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80_lowpt.pdf")
elif pt_mode == 2:
    do_barlow(h_h_yield_dict_0_20, "figures/h_h_yield_barlow_0_20_highpt.pdf")
    do_barlow(h_h_yield_dict_20_50, "figures/h_h_yield_barlow_20_50_highpt.pdf")
    do_barlow(h_h_yield_dict_50_80, "figures/h_h_yield_barlow_50_80_highpt.pdf")
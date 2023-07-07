import ROOT as rt

import math
import array as arr

from time import sleep

rt.gStyle.SetOptStat(0)

c = rt.TCanvas("c","c",800,600)
c.SetLeftMargin(0.10)
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)

colors = [rt.kGreen + 2, rt.kRed + 2, rt.kBlue + 2,  rt.kMagenta + 2,  rt.kYellow + 2,  rt.kViolet + 2, rt.kOrange + 2]

def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error

def do_barlow(von_dict, name):

    color_index = 0
    central_fit = von_dict["central value"]
    central_ns_width = get_width_from_kappa(central_fit.GetParameter(4))
    central_ns_width_error = get_width_error_from_kappa(central_fit.GetParameter(4), central_fit.GetParError(4))
    central_as_width = get_width_from_kappa(central_fit.GetParameter(6))
    central_as_width_error = get_width_error_from_kappa(central_fit.GetParameter(6), central_fit.GetParError(6))
    barlow_legend = rt.TLegend(0.15, 0.15, 0.32, 0.32)
    barlow_legend.SetBorderSize(0)
    for key in von_dict:
        if key == "central value":
            base_hist = rt.TH1D("barlow_base","", 2, 0, 2)
            base_hist.GetXaxis().SetTitle("Region")
            base_hist.GetXaxis().SetBinLabel(1, "Near-side")
            base_hist.GetXaxis().SetBinLabel(2, "Away-side")
            base_hist.GetYaxis().SetTitle("(y_{var}-y_{central})/#sqrt{|#sigma_{var}^{2} - #sigma_{central}^{2}|}")
            base_hist.GetYaxis().SetTitleOffset(1.3)
            base_hist.GetYaxis().SetRangeUser(-10, 10)
            base_hist.SetLineColor(rt.kBlack)
            base_hist.SetLineStyle(0)
            base_hist.SetLineWidth(0)
            base_hist.Draw()
        else:
            current_fit = von_dict[key]
            current_ns_width = get_width_from_kappa(current_fit.GetParameter(4))
            current_ns_width_error = get_width_error_from_kappa(current_fit.GetParameter(4), current_fit.GetParError(4))
            current_as_width = get_width_from_kappa(current_fit.GetParameter(6))
            current_as_width_error = get_width_error_from_kappa(current_fit.GetParameter(6), current_fit.GetParError(6))
            barlow_hist = rt.TH1D("barlow_"+key,"", 2, 0, 2)
            barlow_hist.SetBinContent(1, (current_ns_width - central_ns_width)/math.sqrt(abs(current_ns_width_error**2 - central_ns_width_error**2)))
            barlow_hist.SetBinContent(2, (current_as_width - central_as_width)/math.sqrt(abs(current_as_width_error**2 - central_as_width_error**2)))
            barlow_legend.AddEntry(barlow_hist.Clone(), key, "P")
            barlow_hist.SetLineColor(colors[color_index])
            barlow_hist.SetMarkerColor(colors[color_index])
            barlow_hist.SetMarkerStyle(34)
            barlow_hist.SetMarkerSize(1.5)
            barlow_hist.DrawCopy("SAME P")
        color_index += 1

    barlow_min = rt.TLine(0, -1, 2, -1)
    barlow_max = rt.TLine(0, 1, 2, 1)
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

for pt_mode in range(3):

    # SIGNAL VARIATIONS
    if pt_mode == 0: # CENTRAL 2 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_wide_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_wider_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif pt_mode == 1: # LOW 1.5 - 2.5 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_wide_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_wider_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif pt_mode == 2: # HIGH 2.5 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_wide_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_wider_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
    signal_file_dict = {
        "central value": signal_central_file,
        "narrow": signal_narrow_file,
        "wide": signal_wide_file,
        "narrower": signal_narrower_file,
        "wider": signal_wider_file
    }

    signal_von_dict_0_20 = {}
    signal_von_dict_20_50 = {}
    signal_von_dict_50_80 = {}
    for key, file in signal_file_dict.items():
        signal_von_dict_0_20[key] = file.Get("von_fit_0_20")
        signal_von_dict_20_50[key] = file.Get("von_fit_20_50")
        signal_von_dict_50_80[key] = file.Get("von_fit_50_80")
    if pt_mode == 0:
        do_barlow(signal_von_dict_0_20, "figures/width_signal_barlow_0_20.pdf")
        do_barlow(signal_von_dict_20_50, "figures/width_signal_barlow_20_50.pdf")
        do_barlow(signal_von_dict_50_80, "figures/width_signal_barlow_50_80.pdf")
    elif pt_mode == 1:
        do_barlow(signal_von_dict_0_20, "figures/width_signal_barlow_0_20_lowpt.pdf")
        do_barlow(signal_von_dict_20_50, "figures/width_signal_barlow_20_50_lowpt.pdf")
        do_barlow(signal_von_dict_50_80, "figures/width_signal_barlow_50_80_lowpt.pdf")
    elif pt_mode == 2:
        do_barlow(signal_von_dict_0_20, "figures/width_signal_barlow_0_20_highpt.pdf")
        do_barlow(signal_von_dict_20_50, "figures/width_signal_barlow_20_50_highpt.pdf")
        do_barlow(signal_von_dict_50_80, "figures/width_signal_barlow_50_80_highpt.pdf")

    if pt_mode == 0:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_wide_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_narrow_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif pt_mode == 1:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_wide_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_narrow_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif pt_mode == 2:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_wide_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_narrow_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
    sideband_file_dict = {
        "central value": sideband_central_file,
        "narrow": sideband_narrow_file,
        "wide": sideband_wide_file,
        "shifted right": sideband_shifted_right_file,
        "shifted left": sideband_shifted_left_file
    }

    sideband_von_dict_0_20 = {}
    sideband_von_dict_20_50 = {}
    sideband_von_dict_50_80 = {}
    for key, file in sideband_file_dict.items():
        sideband_von_dict_0_20[key] = file.Get("von_fit_0_20")
        sideband_von_dict_20_50[key] = file.Get("von_fit_20_50")
        sideband_von_dict_50_80[key] = file.Get("von_fit_50_80")
    if pt_mode == 0:
        do_barlow(sideband_von_dict_0_20, "figures/width_sideband_barlow_0_20.pdf")
        do_barlow(sideband_von_dict_20_50, "figures/width_sideband_barlow_20_50.pdf")
        do_barlow(sideband_von_dict_50_80, "figures/width_sideband_barlow_50_80.pdf")
    elif pt_mode == 1:
        do_barlow(sideband_von_dict_0_20, "figures/width_sideband_barlow_0_20_lowpt.pdf")
        do_barlow(sideband_von_dict_20_50, "figures/width_sideband_barlow_20_50_lowpt.pdf")
        do_barlow(sideband_von_dict_50_80, "figures/width_sideband_barlow_50_80_lowpt.pdf")
    elif pt_mode == 2:
        do_barlow(sideband_von_dict_0_20, "figures/width_sideband_barlow_0_20_highpt.pdf")
        do_barlow(sideband_von_dict_20_50, "figures/width_sideband_barlow_20_50_highpt.pdf")
        do_barlow(sideband_von_dict_50_80, "figures/width_sideband_barlow_50_80_highpt.pdf")

    if pt_mode == 0:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_narrow.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_wide.root")
        pid_tof_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_tof.root")
    elif pt_mode == 1:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_narrow_lowpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_wide_lowpt.root")
        pid_tof_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_tof_lowpt.root")
    elif pt_mode == 2:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_narrow_highpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_wide_highpt.root")
        pid_tof_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_tof_highpt.root")
    pid_file_dict = {
        "central value": pid_central_file,
        "narrow": pid_narrow_file,
        "wide": pid_wide_file,
        "require tof": pid_tof_file
    }

    pid_von_dict_0_20 = {}
    pid_von_dict_20_50 = {}
    pid_von_dict_50_80 = {}
    for key, file in pid_file_dict.items():
        pid_von_dict_0_20[key] = file.Get("von_fit_0_20")
        pid_von_dict_20_50[key] = file.Get("von_fit_20_50")
        pid_von_dict_50_80[key] = file.Get("von_fit_50_80")
    if pt_mode == 0:
        do_barlow(pid_von_dict_0_20, "figures/width_pid_barlow_0_20.pdf")
        do_barlow(pid_von_dict_20_50, "figures/width_pid_barlow_20_50.pdf")
        do_barlow(pid_von_dict_50_80, "figures/width_pid_barlow_50_80.pdf")
    elif pt_mode == 1:
        do_barlow(pid_von_dict_0_20, "figures/width_pid_barlow_0_20_lowpt.pdf")
        do_barlow(pid_von_dict_20_50, "figures/width_pid_barlow_20_50_lowpt.pdf")
        do_barlow(pid_von_dict_50_80, "figures/width_pid_barlow_50_80_lowpt.pdf")
    elif pt_mode == 2:
        do_barlow(pid_von_dict_0_20, "figures/width_pid_barlow_0_20_highpt.pdf")
        do_barlow(pid_von_dict_20_50, "figures/width_pid_barlow_20_50_highpt.pdf")
        do_barlow(pid_von_dict_50_80, "figures/width_pid_barlow_50_80_highpt.pdf")
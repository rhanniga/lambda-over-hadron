import ROOT as rt
rt.gStyle.SetOptStat(0)

from array import array
import math

colors = [rt.kGreen + 2, rt.kRed + 2, rt.kBlue + 2,  rt.kMagenta + 2,  rt.kOrange + 2]

def get_von_fit_avgbg(dphi_dist, pt_mode, avg6=True):

    dphi_dist = dphi_dist.Clone()

    if avg6:
        avg = (dphi_dist.GetBinContent(1) 
                + dphi_dist.GetBinContent(2)
                + dphi_dist.GetBinContent(7)
                + dphi_dist.GetBinContent(8)
                + dphi_dist.GetBinContent(9)
                + dphi_dist.GetBinContent(16))/6
    else:
        avg = (dphi_dist.GetBinContent(1)
                + dphi_dist.GetBinContent(8)
                + dphi_dist.GetBinContent(9)
                + dphi_dist.GetBinContent(16))/4
    
    von_fit_string = "[0]"
    von_fit_string += " + [1]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x- TMath::Pi()))"
    von_fit = rt.TF1("von_fit", von_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

    von_fit.FixParameter(0, avg)
    von_fit.SetParameter(1, 0.01)
    von_fit.SetParameter(2, 7)
    von_fit.SetParameter(3, 0.01)
    von_fit.SetParameter(4, 2)

    if pt_mode == 1:
        dphi_dist.Fit(von_fit, "R", "", -rt.TMath.Pi()/2, rt.TMath.Pi() - 0.0001)
    else:
        dphi_dist.Fit(von_fit, "R")

    return von_fit.Clone()

def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error

def get_rms(von_dict, dphi_dict, output_name, pid=False):

    c1 = rt.TCanvas("c1", "c1", 800, 600)
    c1.SetLeftMargin(0.15)

    c2 = rt.TCanvas("c2", "c2", 800, 600)
    c2.SetLeftMargin(0.15)

    color_index = 0

    leg = rt.TLegend(0.6, 0.67, 0.8, 0.8)
    leg.SetBorderSize(0)
    width_leg = rt.TLegend(0.2, 0.67, 0.45, 0.8)
    width_leg.SetBorderSize(0)

    ns_rms = 0
    as_rms = 0
    ns_n = 0
    as_n = 0


    for key in von_dict:
        if key == "central value":

            c1.cd()
            central_fit = von_dict[key].Clone(f"central_fit_{key}")
            central_fit.SetNpx(1000)
            central_fit.GetXaxis().SetRangeUser(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            central_fit.SetLineColor(rt.kBlack)
            central_fit.SetLineWidth(2)
            central_fit.SetTitle("")
            central_fit.GetYaxis().SetRangeUser(0.7*central_fit.GetMinimum(), 1.3*central_fit.GetMaximum())
            central_fit.GetYaxis().SetTitle("Counts")
            central_fit.GetXaxis().SetTitle("#Delta#varphi")
            von_fit_string = "[0]/(2*TMath::Pi()*TMath::BesselI0([1]))*TMath::Exp([1]*TMath::Cos(x))"
            von_fit_string += " + [2]/(2*TMath::Pi()*TMath::BesselI0([3]))*TMath::Exp([3]*TMath::Cos(x- TMath::Pi()))"
            just_von_fit = rt.TF1("just_von_fit", von_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            just_von_fit.SetNpx(1000)
            just_von_fit.SetParameter(0, central_fit.GetParameter(3))
            just_von_fit.SetParameter(1, central_fit.GetParameter(4))
            just_von_fit.SetParameter(2, central_fit.GetParameter(5))
            just_von_fit.SetParameter(3, central_fit.GetParameter(6))
            just_von_fit.SetLineColor(rt.kBlack)
            just_von_fit.SetLineWidth(2)
            just_von_fit.GetYaxis().SetRangeUser(0.7*just_von_fit.GetMinimum(), 1.3*just_von_fit.GetMaximum())
            # just_von_fit.Draw()
            central_fit.Draw()

            central_plot = dphi_dict[key].Clone(f"central_plot_{key}")
            central_plot.GetXaxis().SetRangeUser(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            central_plot.SetLineColor(rt.kBlack)
            central_plot.SetLineWidth(2)
            central_plot.GetYaxis().SetRangeUser(0.7*central_plot.GetMinimum(), 1.3*central_plot.GetMaximum())
            central_plot.GetYaxis().SetTitle("Counts")
            central_plot.GetXaxis().SetTitle("#Delta#varphi")
            central_plot.DrawCopy("SAME")


            c2.cd()
            central_ns_width = get_width_from_kappa(central_fit.GetParameter(4))
            central_ns_width_error = get_width_error_from_kappa(central_fit.GetParameter(4), central_fit.GetParError(4))
            central_as_width = get_width_from_kappa(central_fit.GetParameter(6))
            central_as_width_error = get_width_error_from_kappa(central_fit.GetParameter(6), central_fit.GetParError(6))
            central_width_plot = rt.TH1D("central_width_plot", "", 2, 0, 2)
            central_width_plot.SetBinContent(1, central_ns_width)
            central_width_plot.SetBinError(1, central_ns_width_error)
            central_width_plot.SetBinContent(2, central_as_width)
            central_width_plot.SetBinError(2, central_as_width_error)
            central_width_plot.SetMarkerStyle(20)
            central_width_plot.SetMarkerSize(2)
            central_width_plot.SetMarkerColor(rt.kBlack)
            central_width_plot.SetLineColor(rt.kBlack)
            central_width_plot.SetLineWidth(2)
            central_width_plot.GetYaxis().SetRangeUser(0.7*central_width_plot.GetMinimum(), 1.3*central_width_plot.GetMaximum())
            central_width_plot.GetYaxis().SetTitle("Width")
            central_width_plot.GetXaxis().SetTitle("Region")
            central_width_plot.GetXaxis().SetBinLabel(1, "NS")
            central_width_plot.GetXaxis().SetBinLabel(2, "AS")
            central_width_plot.DrawCopy("P")

            width_leg.AddEntry(central_width_plot, key, "l")
            leg.AddEntry(central_width_plot, key, "l")

        else:
            c1.cd()
            varied_fit = von_dict[key].Clone(f"varied_fit_{key}")
            varied_fit.SetNpx(1000)
            if pid:
                l = -rt.TMath.Pi()/2
                r = 3*rt.TMath.Pi()/2
                scale_factor = central_fit.Integral(l, r)/varied_fit.Integral(l, r)
                varied_fit.SetParameter(0, varied_fit.GetParameter(0)*scale_factor)
                varied_fit.SetParameter(3, varied_fit.GetParameter(3)*scale_factor)
                varied_fit.SetParameter(5, varied_fit.GetParameter(5)*scale_factor)
            varied_fit.SetLineColor(colors[color_index])
            varied_fit.SetLineWidth(2)

            # if key == "gaus":
            #     just_gaus_fit = rt.TF1("just_gaus_fit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            #     just_gaus_fit.SetNpx(1000)
            #     just_gaus_fit.SetParameter(0, varied_fit.GetParameter(3))
            #     just_gaus_fit.SetParameter(1, varied_fit.GetParameter(4))
            #     just_gaus_fit.SetParameter(2, varied_fit.GetParameter(5))
            #     just_gaus_fit.SetParameter(3, varied_fit.GetParameter(6))
            #     just_gaus_fit.SetParameter(4, varied_fit.GetParameter(7))
            #     just_gaus_fit.SetParameter(5, varied_fit.GetParameter(8))
            #     just_gaus_fit.SetParameter(6, varied_fit.GetParameter(9))
            #     just_gaus_fit.SetParameter(7, varied_fit.GetParameter(10))
            #     just_gaus_fit.SetParameter(8, varied_fit.GetParameter(11))
            #     just_gaus_fit.SetParameter(9, varied_fit.GetParameter(12))
            #     just_gaus_fit.SetParameter(10, varied_fit.GetParameter(13))
            #     just_gaus_fit.SetParameter(11, varied_fit.GetParameter(14))
            #     just_gaus_fit.GetYaxis().SetRangeUser(0.7*just_gaus_fit.GetMinimum(), 1.3*just_gaus_fit.GetMaximum())
            #     just_gaus_fit.Draw("same")

            # else:
            # varied_fit.DrawCopy("same")
            varied_fit.DrawCopy("same")

            varied_plot = dphi_dict[key].Clone(f"varied_plot_{key}")
            if pid:
                scale_factor = central_plot.Integral()/varied_plot.Integral()
                print(f"THE PID SCALE FACTOR IS: {scale_factor}\n KEY: {key}")
                varied_plot.Scale(central_plot.Integral()/varied_plot.Integral())
            varied_plot.SetLineColor(colors[color_index])
            varied_plot.SetLineWidth(2)
            if key not in ["gaus", "avg6", "avg4"]:
                varied_plot.DrawCopy("same")

            c2.cd()
            if key == "gaus":
                varied_ns_width = varied_fit.GetParameter(5)
                varied_ns_width_error = varied_fit.GetParError(5)
                varied_as_width = varied_fit.GetParameter(11)
                varied_as_width_error = varied_fit.GetParError(11)
            else:

                if key == "avg6" or key == "avg4":
                    varied_ns_width = get_width_from_kappa(varied_fit.GetParameter(2))
                    varied_ns_width_error = get_width_error_from_kappa(varied_fit.GetParameter(2), varied_fit.GetParError(2))
                    varied_as_width = get_width_from_kappa(varied_fit.GetParameter(4))
                    varied_as_width_error = get_width_error_from_kappa(varied_fit.GetParameter(4), varied_fit.GetParError(4))
                else:
                    varied_ns_width = get_width_from_kappa(varied_fit.GetParameter(4))
                    varied_ns_width_error = get_width_error_from_kappa(varied_fit.GetParameter(4), varied_fit.GetParError(4))
                    varied_as_width = get_width_from_kappa(varied_fit.GetParameter(6))
                    varied_as_width_error = get_width_error_from_kappa(varied_fit.GetParameter(6), varied_fit.GetParError(6))

            varied_width_plot = rt.TH1D("varied_width_plot", "varied_width_plot", 2, 0, 2)
            varied_width_plot.SetBinContent(1, varied_ns_width)
            varied_width_plot.SetBinError(1, varied_ns_width_error)
            varied_width_plot.SetBinContent(2, varied_as_width)
            varied_width_plot.SetBinError(2, varied_as_width_error)
            varied_width_plot.SetMarkerStyle(20)
            varied_width_plot.SetMarkerSize(2)
            varied_width_plot.SetMarkerColor(colors[color_index])
            varied_width_plot.SetLineColor(colors[color_index])
            varied_width_plot.SetLineWidth(2)
            varied_width_plot.GetYaxis().SetRangeUser(0.7*varied_width_plot.GetMinimum(), 1.3*varied_width_plot.GetMaximum())
            varied_width_plot.GetYaxis().SetTitle("Width")
            varied_width_plot.GetXaxis().SetTitle("Region")
            varied_width_plot.GetXaxis().SetBinLabel(1, "NS")
            varied_width_plot.GetXaxis().SetBinLabel(2, "AS")
            width_leg.AddEntry(varied_width_plot.Clone(), key, "l")
            leg.AddEntry(varied_width_plot.Clone(), key, "l")
            varied_width_plot.DrawCopy("P SAME")

            ns_rms += (1 - (varied_ns_width/central_ns_width))**2
            as_rms += (1 - (varied_as_width/central_as_width))**2

            ns_n += 1
            as_n += 1

        color_index += 1


    
    c1.cd()
    leg.Draw("same")
    c1.SaveAs(output_name + ".pdf")

    c2.cd()
    width_leg.Draw("same")
    c2.SaveAs(output_name + "_widths.pdf")

    return math.sqrt(ns_rms/ns_n), math.sqrt(as_rms/as_n)

outfile = rt.TFile("width_syst_out.root", "RECREATE")
for PT_MODE in range(3):

    print("PT MODE: " + str(PT_MODE))
    # SIGNAL VARIATIONS
    if PT_MODE == 0: # CENTRAL 2 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1: # LOW 1.5 - 2.5 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2: # HIGH 2.5 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    signal_file_dict = {
        "central value": signal_central_file,
        "narrow": signal_narrow_file,
        "narrower": signal_narrower_file
    }

    signal_von_dict_0_20 = {}
    signal_von_dict_20_50 = {}
    signal_von_dict_50_80 = {}
    signal_von_dict_0_80 = {}

    signal_dphi_dict_0_20 = {}
    signal_dphi_dict_20_50 = {}
    signal_dphi_dict_50_80 = {}
    signal_dphi_dict_0_80 = {}

    for key, file in signal_file_dict.items():

        signal_von_dict_0_20[key] = file.Get("von_fit_0_20")
        signal_von_dict_20_50[key] = file.Get("von_fit_20_50")
        signal_von_dict_50_80[key] = file.Get("von_fit_50_80")
        signal_von_dict_0_80[key] = file.Get("von_fit_0_80")

        signal_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        signal_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        signal_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        signal_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80")
    elif PT_MODE == 1:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20_lowpt")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50_lowpt")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80_lowpt")
    elif PT_MODE == 2:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20_highpt")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50_highpt")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80_highpt")
    


    if PT_MODE == 0:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    sideband_file_dict = {
        "central value": sideband_central_file,
        "shifted right": sideband_shifted_right_file,
        "shifted left": sideband_shifted_left_file
    }

    sideband_von_dict_0_20 = {}
    sideband_von_dict_20_50 = {}
    sideband_von_dict_50_80 = {}
    sideband_von_dict_0_80 = {}

    sideband_dphi_dict_0_20 = {}
    sideband_dphi_dict_20_50 = {}
    sideband_dphi_dict_50_80 = {}
    sideband_dphi_dict_0_80 = {}

    for key, file in sideband_file_dict.items():

        sideband_von_dict_0_20[key] = file.Get("von_fit_0_20")
        sideband_von_dict_20_50[key] = file.Get("von_fit_20_50")
        sideband_von_dict_50_80[key] = file.Get("von_fit_50_80")
        sideband_von_dict_0_80[key] = file.Get("von_fit_0_80")

        sideband_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        sideband_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        sideband_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        sideband_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80")
    elif PT_MODE == 1:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20_lowpt")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50_lowpt")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80_lowpt")
    elif PT_MODE == 2:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20_highpt")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50_highpt")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80_highpt")


    if PT_MODE == 0:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_narrow.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_wide.root")
    elif PT_MODE == 1:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_narrow_lowpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_wide_lowpt.root")
    elif PT_MODE == 2:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_narrow_highpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_wide_highpt.root")


    pid_file_dict = {
        "central value": pid_central_file,
        "narrow": pid_narrow_file,
        "wide": pid_wide_file
    }

    pid_von_dict_0_20 = {}
    pid_von_dict_20_50 = {}
    pid_von_dict_50_80 = {}
    pid_von_dict_0_80 = {}

    pid_dphi_dict_0_20 = {}
    pid_dphi_dict_20_50 = {}
    pid_dphi_dict_50_80 = {}
    pid_dphi_dict_0_80 = {}

    for key, file in pid_file_dict.items():

        pid_von_dict_0_20[key] = file.Get("von_fit_0_20")
        pid_von_dict_20_50[key] = file.Get("von_fit_20_50")
        pid_von_dict_50_80[key] = file.Get("von_fit_50_80")
        pid_von_dict_0_80[key] = file.Get("von_fit_0_80")

        pid_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        pid_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        pid_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        pid_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80", True)
    elif PT_MODE == 1:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20_lowpt", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50_lowpt", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80_lowpt", True)
    elif PT_MODE == 2:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20_highpt", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50_highpt", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80_highpt", True)


    if PT_MODE == 0:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")


    technique_file_dict = {
        "central value": technique_central_file,
        "gaus": technique_gaus_file,
        "avg6" : technique_central_file,
        "avg4" : technique_central_file
    }

    technique_von_dict_0_20 = {}
    technique_von_dict_20_50 = {}
    technique_von_dict_50_80 = {}
    technique_von_dict_0_80 = {}
    technique_dphi_dict_0_20 = {}
    technique_dphi_dict_20_50 = {}
    technique_dphi_dict_50_80 = {}
    technique_dphi_dict_0_80 = {}

    hh_technique_von_dict_0_20 = {}
    hh_technique_von_dict_20_50 = {}
    hh_technique_von_dict_50_80 = {}
    hh_technique_von_dict_0_80 = {}

    hh_technique_dphi_dict_0_20 = {}
    hh_technique_dphi_dict_20_50 = {}
    hh_technique_dphi_dict_50_80 = {}
    hh_technique_dphi_dict_0_80 = {}

    for key, file in technique_file_dict.items():

        technique_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        technique_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        technique_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        technique_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

        hh_technique_dphi_dict_0_20[key] = file.Get("h_h_dphi_0_20")
        hh_technique_dphi_dict_20_50[key] = file.Get("h_h_dphi_20_50")
        hh_technique_dphi_dict_50_80[key] = file.Get("h_h_dphi_50_80")
        hh_technique_dphi_dict_0_80[key] = file.Get("h_h_dphi_0_80")

        if key == "gaus":

            technique_von_dict_0_20[key] = file.Get("fit_function_0_20")
            technique_von_dict_20_50[key] = file.Get("fit_function_20_50")
            technique_von_dict_50_80[key] = file.Get("fit_function_50_80")
            technique_von_dict_0_80[key] = file.Get("fit_function_0_80")
            
            hh_technique_von_dict_0_20[key] = file.Get("hh_fit_function_0_20")
            hh_technique_von_dict_20_50[key] = file.Get("hh_fit_function_20_50")
            hh_technique_von_dict_50_80[key] = file.Get("hh_fit_function_50_80")
            hh_technique_von_dict_0_80[key] = file.Get("hh_fit_function_0_80")

        elif key == "central value":

            technique_von_dict_0_20[key] = file.Get("von_fit_0_20")
            technique_von_dict_20_50[key] = file.Get("von_fit_20_50")
            technique_von_dict_50_80[key] = file.Get("von_fit_50_80")
            technique_von_dict_0_80[key] = file.Get("von_fit_0_80")

            hh_technique_von_dict_0_20[key] = file.Get("hh_von_fit_0_20")
            hh_technique_von_dict_20_50[key] = file.Get("hh_von_fit_20_50")
            hh_technique_von_dict_50_80[key] = file.Get("hh_von_fit_50_80")
            hh_technique_von_dict_0_80[key] = file.Get("hh_von_fit_0_80")

        elif key == "avg6":

            technique_von_dict_0_20[key] = get_von_fit_avgbg(technique_dphi_dict_0_20[key], PT_MODE)
            technique_von_dict_20_50[key] = get_von_fit_avgbg(technique_dphi_dict_20_50[key], PT_MODE)
            technique_von_dict_50_80[key] = get_von_fit_avgbg(technique_dphi_dict_50_80[key], PT_MODE)
            technique_von_dict_0_80[key] = get_von_fit_avgbg(technique_dphi_dict_0_80[key], PT_MODE)

            hh_technique_von_dict_0_20[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_20[key], PT_MODE)
            hh_technique_von_dict_20_50[key] = get_von_fit_avgbg(hh_technique_dphi_dict_20_50[key], PT_MODE)
            hh_technique_von_dict_50_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_50_80[key], PT_MODE)
            hh_technique_von_dict_0_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_80[key], PT_MODE)

        elif key == "avg4":

            technique_von_dict_0_20[key] = get_von_fit_avgbg(technique_dphi_dict_0_20[key], PT_MODE, False)
            technique_von_dict_20_50[key] = get_von_fit_avgbg(technique_dphi_dict_20_50[key], PT_MODE, False)
            technique_von_dict_50_80[key] = get_von_fit_avgbg(technique_dphi_dict_50_80[key], PT_MODE, False)
            technique_von_dict_0_80[key] = get_von_fit_avgbg(technique_dphi_dict_0_80[key], PT_MODE, False)

            hh_technique_von_dict_0_20[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_20[key], PT_MODE, False)
            hh_technique_von_dict_20_50[key] = get_von_fit_avgbg(hh_technique_dphi_dict_20_50[key], PT_MODE, False)
            hh_technique_von_dict_50_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_50_80[key], PT_MODE, False)
            hh_technique_von_dict_0_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_80[key], PT_MODE, False)



    if PT_MODE == 0:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80")

        
    elif PT_MODE == 1:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20_lowpt")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50_lowpt")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80_lowpt")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20_lowpt")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50_lowpt")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80_lowpt")

    elif PT_MODE == 2:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20_highpt")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50_highpt")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80_highpt")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20_highpt")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50_highpt")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80_highpt")



    ns_total_rms_0_20 = math.sqrt(ns_signal_rms_0_20**2 + ns_sideband_rms_0_20**2 + ns_pid_rms_0_20**2 + ns_technique_rms_0_20**2)
    print("ns_total_rms_0_20: " + str(ns_total_rms_0_20))
    as_total_rms_0_20 = math.sqrt(as_signal_rms_0_20**2 + as_sideband_rms_0_20**2 + as_pid_rms_0_20**2 + as_technique_rms_0_20**2)
    print("as_total_rms_0_20: " + str(as_total_rms_0_20))

    c = rt.TCanvas("c", "c", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    ns_sideband_syst_hist = rt.TH1D("ns_sideband_syst_hist", "", 3, 0, 3)
    ns_sideband_syst_hist.SetBinContent(1, ns_sideband_rms_0_20)
    ns_sideband_syst_hist.SetBinContent(2, ns_sideband_rms_20_50)
    ns_sideband_syst_hist.SetBinContent(3, ns_sideband_rms_50_80)
    ns_sideband_syst_hist.SetBinError(1, 0)
    ns_sideband_syst_hist.SetBinError(2, 0)
    ns_sideband_syst_hist.SetBinError(3, 0)
    ns_sideband_syst_hist.SetLineColor(rt.kGreen - 2)
    ns_sideband_syst_hist.SetLineWidth(2)
    ns_signal_syst_hist = rt.TH1D("ns_signal_syst_hist", "", 3, 0, 3)
    ns_signal_syst_hist.SetBinContent(1, ns_signal_rms_0_20)
    ns_signal_syst_hist.SetBinContent(2, ns_signal_rms_20_50)
    ns_signal_syst_hist.SetBinContent(3, ns_signal_rms_50_80)
    ns_signal_syst_hist.SetBinError(1, 0)
    ns_signal_syst_hist.SetBinError(2, 0)
    ns_signal_syst_hist.SetBinError(3, 0)
    ns_signal_syst_hist.SetLineColor(rt.kMagenta - 2)
    ns_signal_syst_hist.SetLineWidth(2)
    ns_pid_syst_hist = rt.TH1D("ns_pid_syst_hist", "", 3, 0, 3)
    ns_pid_syst_hist.SetBinContent(1, ns_pid_rms_0_20)
    ns_pid_syst_hist.SetBinContent(2, ns_pid_rms_20_50)
    ns_pid_syst_hist.SetBinContent(3, ns_pid_rms_50_80)
    ns_pid_syst_hist.SetBinError(1, 0)
    ns_pid_syst_hist.SetBinError(2, 0)
    ns_pid_syst_hist.SetBinError(3, 0)
    ns_pid_syst_hist.SetLineColor(rt.kAzure - 2)
    ns_pid_syst_hist.SetLineWidth(2)
    ns_technique_syst_hist = rt.TH1D("ns_technique_syst_hist", "", 3, 0, 3)
    ns_technique_syst_hist.SetBinContent(1, ns_technique_rms_0_20)
    ns_technique_syst_hist.SetBinContent(2, ns_technique_rms_20_50)
    ns_technique_syst_hist.SetBinContent(3, ns_technique_rms_50_80)
    ns_technique_syst_hist.SetBinError(1, 0)
    ns_technique_syst_hist.SetBinError(2, 0)
    ns_technique_syst_hist.SetBinError(3, 0)
    ns_technique_syst_hist.SetLineColor(rt.kViolet + 2)
    ns_technique_syst_hist.SetLineWidth(2)


    as_sideband_syst_hist = rt.TH1D("as_sideband_syst_hist", "", 3, 0, 3)
    as_sideband_syst_hist.SetBinContent(1, as_sideband_rms_0_20)
    as_sideband_syst_hist.SetBinContent(2, as_sideband_rms_20_50)
    as_sideband_syst_hist.SetBinContent(3, as_sideband_rms_50_80)
    as_sideband_syst_hist.SetBinError(1, 0)
    as_sideband_syst_hist.SetBinError(2, 0)
    as_sideband_syst_hist.SetBinError(3, 0)
    as_sideband_syst_hist.SetLineColor(rt.kGreen - 2)
    as_sideband_syst_hist.SetLineWidth(2)
    as_signal_syst_hist = rt.TH1D("as_signal_syst_hist", "", 3, 0, 3)
    as_signal_syst_hist.SetBinContent(1, as_signal_rms_0_20)
    as_signal_syst_hist.SetBinContent(2, as_signal_rms_20_50)
    as_signal_syst_hist.SetBinContent(3, as_signal_rms_50_80)
    as_signal_syst_hist.SetBinError(1, 0)
    as_signal_syst_hist.SetBinError(2, 0)
    as_signal_syst_hist.SetBinError(3, 0)
    as_signal_syst_hist.SetLineColor(rt.kMagenta - 2)
    as_signal_syst_hist.SetLineWidth(2)
    as_pid_syst_hist = rt.TH1D("as_pid_syst_hist", "", 3, 0, 3)
    as_pid_syst_hist.SetBinContent(1, as_pid_rms_0_20)
    as_pid_syst_hist.SetBinContent(2, as_pid_rms_20_50)
    as_pid_syst_hist.SetBinContent(3, as_pid_rms_50_80)
    as_pid_syst_hist.SetBinError(1, 0)
    as_pid_syst_hist.SetBinError(2, 0)
    as_pid_syst_hist.SetBinError(3, 0)
    as_pid_syst_hist.SetLineColor(rt.kAzure - 2)
    as_pid_syst_hist.SetLineWidth(2)
    as_technique_syst_hist = rt.TH1D("as_technique_syst_hist", "", 3, 0, 3)
    as_technique_syst_hist.SetBinContent(1, as_technique_rms_0_20)
    as_technique_syst_hist.SetBinContent(2, as_technique_rms_20_50)
    as_technique_syst_hist.SetBinContent(3, as_technique_rms_50_80)
    as_technique_syst_hist.SetBinError(1, 0)
    as_technique_syst_hist.SetBinError(2, 0)
    as_technique_syst_hist.SetBinError(3, 0)
    as_technique_syst_hist.SetLineColor(rt.kViolet + 2)
    as_technique_syst_hist.SetLineWidth(2)

    # the topo syst is a LITTLE higher at low pT
    if PT_MODE == 1:
        topo_rms = 0.032
        topo_nch_dep_rms = 0.005
        topo_syst_hist = rt.TH1D("topo_syst_hist", "", 3, 0, 3)
        topo_syst_hist.SetBinContent(1, topo_rms)
        topo_syst_hist.SetBinContent(2, topo_rms)
        topo_syst_hist.SetBinContent(3, topo_rms)
        topo_syst_hist.SetBinError(1, 0)
        topo_syst_hist.SetBinError(2, 0)
        topo_syst_hist.SetBinError(3, 0)
        topo_syst_hist.SetLineColor(rt.kRed - 2)
        topo_syst_hist.SetLineWidth(2)
    else:
        topo_rms = 0.030
        topo_nch_dep_rms = 0.005
        topo_syst_hist = rt.TH1D("topo_syst_hist", "", 3, 0, 3)
        topo_syst_hist.SetBinContent(1, topo_rms)
        topo_syst_hist.SetBinContent(2, topo_rms)
        topo_syst_hist.SetBinContent(3, topo_rms)
        topo_syst_hist.SetBinError(1, 0)
        topo_syst_hist.SetBinError(2, 0)
        topo_syst_hist.SetBinError(3, 0)
        topo_syst_hist.SetLineColor(rt.kRed - 2)
        topo_syst_hist.SetLineWidth(2)

    # the matbud syst is a LITTLE higher at low pT
    if PT_MODE == 1:
        matbud_rms = 0.011
        matbud_nch_dep_rms = 0.0
        matbud_syst_hist = rt.TH1D("matbud_syst_hist", "", 3, 0, 3)
        matbud_syst_hist.SetBinContent(1, matbud_rms)
        matbud_syst_hist.SetBinContent(2, matbud_rms)
        matbud_syst_hist.SetBinContent(3, matbud_rms)
        matbud_syst_hist.SetBinError(1, 0)
        matbud_syst_hist.SetBinError(2, 0)
        matbud_syst_hist.SetBinError(3, 0)
        matbud_syst_hist.SetLineColor(rt.kOrange - 2)
        matbud_syst_hist.SetLineWidth(2)
    else:
        matbud_rms = 0.006
        matbud_nch_dep_rms = 0.000
        matbud_syst_hist = rt.TH1D("matbud_syst_hist", "", 3, 0, 3)
        matbud_syst_hist.SetBinContent(1, matbud_rms)
        matbud_syst_hist.SetBinContent(2, matbud_rms)
        matbud_syst_hist.SetBinContent(3, matbud_rms)
        matbud_syst_hist.SetBinError(1, 0)
        matbud_syst_hist.SetBinError(2, 0)
        matbud_syst_hist.SetBinError(3, 0)
        matbud_syst_hist.SetLineColor(rt.kOrange - 2)
        matbud_syst_hist.SetLineWidth(2)

    ns_total_syst_hist = rt.TH1D("ns_total_syst_hist", "", 3, 0, 3)
    ns_total_syst_hist.SetBinContent(1, math.sqrt(ns_sideband_rms_0_20**2 + ns_signal_rms_0_20**2 + ns_pid_rms_0_20**2 + ns_technique_rms_0_20**2 + topo_rms**2 + matbud_rms**2 ))
    ns_total_syst_hist.SetBinContent(2, math.sqrt(ns_sideband_rms_20_50**2 + ns_signal_rms_20_50**2 + ns_pid_rms_20_50**2 + ns_technique_rms_20_50**2 + topo_rms**2 + matbud_rms**2))
    ns_total_syst_hist.SetBinContent(3, math.sqrt(ns_sideband_rms_50_80**2 + ns_signal_rms_50_80**2 + ns_pid_rms_50_80**2 + ns_technique_rms_50_80**2 + topo_rms**2 + matbud_rms**2))
    ns_total_syst_hist.SetBinError(1, 0)
    ns_total_syst_hist.SetBinError(2, 0)
    ns_total_syst_hist.SetBinError(3, 0)
    ns_total_syst_hist.SetLineColor(rt.kBlack)
    ns_total_syst_hist.SetLineWidth(2)
    ns_total_syst_hist.SetTitle("")
    ns_total_syst_hist.GetXaxis().SetBinLabel(1, "0-20%")
    ns_total_syst_hist.GetXaxis().SetBinLabel(2, "20-50%")
    ns_total_syst_hist.GetXaxis().SetBinLabel(3, "50-80%")
    ns_total_syst_hist.GetXaxis().SetTitle("Mult. Percentile")
    ns_total_syst_hist.GetYaxis().SetTitle("Systematic Uncertainty (%)")

    as_total_syst_hist = rt.TH1D("as_total_syst_hist", "", 3, 0, 3)
    as_total_syst_hist.SetBinContent(1, math.sqrt(as_sideband_rms_0_20**2 + as_signal_rms_0_20**2 + as_pid_rms_0_20**2 + as_technique_rms_0_20**2 + topo_rms**2 + matbud_rms**2 ))
    as_total_syst_hist.SetBinContent(2, math.sqrt(as_sideband_rms_20_50**2 + as_signal_rms_20_50**2 + as_pid_rms_20_50**2 + as_technique_rms_20_50**2 + topo_rms**2 + matbud_rms**2))
    as_total_syst_hist.SetBinContent(3, math.sqrt(as_sideband_rms_50_80**2 + as_signal_rms_50_80**2 + as_pid_rms_50_80**2 + as_technique_rms_50_80**2 + topo_rms**2 + matbud_rms**2))
    as_total_syst_hist.SetBinError(1, 0)
    as_total_syst_hist.SetBinError(2, 0)
    as_total_syst_hist.SetBinError(3, 0)
    as_total_syst_hist.SetLineColor(rt.kBlack)
    as_total_syst_hist.SetLineWidth(2)
    as_total_syst_hist.SetTitle("")
    as_total_syst_hist.GetXaxis().SetBinLabel(1, "0-20%")
    as_total_syst_hist.GetXaxis().SetBinLabel(2, "20-50%")
    as_total_syst_hist.GetXaxis().SetBinLabel(3, "50-80%")
    as_total_syst_hist.GetXaxis().SetTitle("Mult. Percentile")
    as_total_syst_hist.GetYaxis().SetTitle("Systematic Uncertainty (%)")

    leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(ns_sideband_syst_hist, "Sideband", "l")
    leg.AddEntry(ns_signal_syst_hist, "Signal", "l")
    leg.AddEntry(ns_pid_syst_hist, "PID", "l")
    leg.AddEntry(ns_technique_syst_hist, "Technique", "l")
    leg.AddEntry(topo_syst_hist, "Topo", "l")
    leg.AddEntry(matbud_syst_hist, "Mat. budget", "l")
    leg.AddEntry(ns_total_syst_hist, "Total", "l")

    ns_total_syst_hist.Scale(100)
    ns_sideband_syst_hist.Scale(100)
    ns_signal_syst_hist.Scale(100)
    ns_pid_syst_hist.Scale(100)
    ns_technique_syst_hist.Scale(100)

    as_total_syst_hist.Scale(100)
    as_sideband_syst_hist.Scale(100)
    as_signal_syst_hist.Scale(100)
    as_pid_syst_hist.Scale(100)
    as_technique_syst_hist.Scale(100)

    topo_syst_hist.Scale(100)
    matbud_syst_hist.Scale(100)

    ns_total_syst_hist.GetYaxis().SetRangeUser(0, ns_total_syst_hist.GetMaximum()*1.8)
    as_total_syst_hist.GetYaxis().SetRangeUser(0, as_total_syst_hist.GetMaximum()*1.8)


    ns_total_syst_hist.Draw("hist")
    ns_sideband_syst_hist.Draw("hist same")
    ns_signal_syst_hist.Draw("hist same")
    ns_pid_syst_hist.Draw("hist same")
    ns_technique_syst_hist.Draw("hist same")
    topo_syst_hist.Draw("hist same")
    matbud_syst_hist.Draw("hist same")

    leg.Draw("SAME")
    if PT_MODE == 0:
        c.SaveAs("../figures/systematics_ns_width_postbarlow.pdf")
    elif PT_MODE == 1:
        c.SaveAs("../figures/systematics_ns_width_postbarlow_lowpt.pdf")
    elif PT_MODE == 2:
        c.SaveAs("../figures/systematics_ns_width_postbarlow_highpt.pdf")

    as_total_syst_hist.Draw("hist")
    as_sideband_syst_hist.Draw("hist same")
    as_signal_syst_hist.Draw("hist same")
    as_pid_syst_hist.Draw("hist same")
    as_technique_syst_hist.Draw("hist same")
    topo_syst_hist.Draw("hist same")
    matbud_syst_hist.Draw("hist same")

    leg.Draw("SAME")
    if PT_MODE == 0:
        c.SaveAs("../figures/systematics_as_width_postbarlow.pdf")
    elif PT_MODE == 1:
        c.SaveAs("../figures/systematics_as_width_postbarlow_lowpt.pdf")
    elif PT_MODE == 2:
        c.SaveAs("../figures/systematics_as_width_postbarlow_highpt.pdf")

    
    as_total_syst_0_20 = as_total_syst_hist.GetBinContent(1)
    as_total_syst_20_50 = as_total_syst_hist.GetBinContent(2)
    as_total_syst_50_80 = as_total_syst_hist.GetBinContent(3)

    ns_total_syst_0_20 = ns_total_syst_hist.GetBinContent(1)
    ns_total_syst_20_50 = ns_total_syst_hist.GetBinContent(2)
    ns_total_syst_50_80 = ns_total_syst_hist.GetBinContent(3)
    print("PT MODE: " + str(PT_MODE))
    # SIGNAL VARIATIONS
    if PT_MODE == 0: # CENTRAL 2 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1: # LOW 1.5 - 2.5 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2: # HIGH 2.5 - 4 BIN
        signal_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrow_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrower_file = rt.TFile.Open("../output/signal_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    signal_file_dict = {
        "central value": signal_central_file,
        "narrow": signal_narrow_file,
        "narrower": signal_narrower_file
    }

    signal_von_dict_0_20 = {}
    signal_von_dict_20_50 = {}
    signal_von_dict_50_80 = {}
    signal_von_dict_0_80 = {}

    signal_dphi_dict_0_20 = {}
    signal_dphi_dict_20_50 = {}
    signal_dphi_dict_50_80 = {}
    signal_dphi_dict_0_80 = {}

    for key, file in signal_file_dict.items():

        signal_von_dict_0_20[key] = file.Get("von_fit_0_20")
        signal_von_dict_20_50[key] = file.Get("von_fit_20_50")
        signal_von_dict_50_80[key] = file.Get("von_fit_50_80")
        signal_von_dict_0_80[key] = file.Get("von_fit_0_80")

        signal_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        signal_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        signal_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        signal_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80")
    elif PT_MODE == 1:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20_lowpt")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50_lowpt")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80_lowpt")
    elif PT_MODE == 2:
        ns_signal_rms_0_20, as_signal_rms_0_20 = get_rms(signal_von_dict_0_20, signal_dphi_dict_0_20, "../figures/signal_variations_width_0_20_highpt")
        ns_signal_rms_20_50, as_signal_rms_20_50 = get_rms(signal_von_dict_20_50, signal_dphi_dict_20_50, "../figures/signal_variations_width_20_50_highpt")
        ns_signal_rms_50_80, as_signal_rms_50_80 = get_rms(signal_von_dict_50_80, signal_dphi_dict_50_80, "../figures/signal_variations_width_50_80_highpt")
    


    if PT_MODE == 0:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        sideband_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_right_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_left_file = rt.TFile.Open("../output/sideband_variation/v0_von_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    sideband_file_dict = {
        "central value": sideband_central_file,
        "shifted right": sideband_shifted_right_file,
        "shifted left": sideband_shifted_left_file
    }

    sideband_von_dict_0_20 = {}
    sideband_von_dict_20_50 = {}
    sideband_von_dict_50_80 = {}
    sideband_von_dict_0_80 = {}

    sideband_dphi_dict_0_20 = {}
    sideband_dphi_dict_20_50 = {}
    sideband_dphi_dict_50_80 = {}
    sideband_dphi_dict_0_80 = {}

    for key, file in sideband_file_dict.items():

        sideband_von_dict_0_20[key] = file.Get("von_fit_0_20")
        sideband_von_dict_20_50[key] = file.Get("von_fit_20_50")
        sideband_von_dict_50_80[key] = file.Get("von_fit_50_80")
        sideband_von_dict_0_80[key] = file.Get("von_fit_0_80")

        sideband_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        sideband_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        sideband_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        sideband_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80")
    elif PT_MODE == 1:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20_lowpt")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50_lowpt")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80_lowpt")
    elif PT_MODE == 2:
        ns_sideband_rms_0_20, as_sideband_rms_0_20 = get_rms(sideband_von_dict_0_20, sideband_dphi_dict_0_20, "../figures/sideband_variations_width_0_20_highpt")
        ns_sideband_rms_20_50, as_sideband_rms_20_50 = get_rms(sideband_von_dict_20_50, sideband_dphi_dict_20_50, "../figures/sideband_variations_width_20_50_highpt")
        ns_sideband_rms_50_80, as_sideband_rms_50_80 = get_rms(sideband_von_dict_50_80, sideband_dphi_dict_50_80, "../figures/sideband_variations_width_50_80_highpt")


    if PT_MODE == 0:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_narrow.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_wide.root")
    elif PT_MODE == 1:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_narrow_lowpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_wide_lowpt.root")
    elif PT_MODE == 2:
        pid_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        pid_narrow_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_narrow_highpt.root")
        pid_wide_file = rt.TFile.Open("../output/pid_variation/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_wide_highpt.root")


    pid_file_dict = {
        "central value": pid_central_file,
        "narrow": pid_narrow_file,
        "wide": pid_wide_file
    }

    pid_von_dict_0_20 = {}
    pid_von_dict_20_50 = {}
    pid_von_dict_50_80 = {}
    pid_von_dict_0_80 = {}

    pid_dphi_dict_0_20 = {}
    pid_dphi_dict_20_50 = {}
    pid_dphi_dict_50_80 = {}
    pid_dphi_dict_0_80 = {}

    for key, file in pid_file_dict.items():

        pid_von_dict_0_20[key] = file.Get("von_fit_0_20")
        pid_von_dict_20_50[key] = file.Get("von_fit_20_50")
        pid_von_dict_50_80[key] = file.Get("von_fit_50_80")
        pid_von_dict_0_80[key] = file.Get("von_fit_0_80")

        pid_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        pid_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        pid_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        pid_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

    if PT_MODE == 0:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80", True)
    elif PT_MODE == 1:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20_lowpt", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50_lowpt", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80_lowpt", True)
    elif PT_MODE == 2:
        ns_pid_rms_0_20, as_pid_rms_0_20 = get_rms(pid_von_dict_0_20, pid_dphi_dict_0_20, "../figures/pid_variations_width_0_20_highpt", True)
        ns_pid_rms_20_50, as_pid_rms_20_50 = get_rms(pid_von_dict_20_50, pid_dphi_dict_20_50, "../figures/pid_variations_width_20_50_highpt", True)
        ns_pid_rms_50_80, as_pid_rms_50_80 = get_rms(pid_von_dict_50_80, pid_dphi_dict_50_80, "../figures/pid_variations_width_50_80_highpt", True)


    if PT_MODE == 0:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        technique_central_file = rt.TFile.Open("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        technique_gaus_file = rt.TFile.Open("../output/technique_variation/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")


    technique_file_dict = {
        "central value": technique_central_file,
        "gaus": technique_gaus_file,
        "avg6" : technique_central_file,
        "avg4" : technique_central_file
    }

    technique_von_dict_0_20 = {}
    technique_von_dict_20_50 = {}
    technique_von_dict_50_80 = {}
    technique_von_dict_0_80 = {}
    technique_dphi_dict_0_20 = {}
    technique_dphi_dict_20_50 = {}
    technique_dphi_dict_50_80 = {}
    technique_dphi_dict_0_80 = {}

    hh_technique_von_dict_0_20 = {}
    hh_technique_von_dict_20_50 = {}
    hh_technique_von_dict_50_80 = {}
    hh_technique_von_dict_0_80 = {}

    hh_technique_dphi_dict_0_20 = {}
    hh_technique_dphi_dict_20_50 = {}
    hh_technique_dphi_dict_50_80 = {}
    hh_technique_dphi_dict_0_80 = {}

    for key, file in technique_file_dict.items():

        technique_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        technique_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        technique_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        technique_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")

        hh_technique_dphi_dict_0_20[key] = file.Get("h_h_dphi_0_20")
        hh_technique_dphi_dict_20_50[key] = file.Get("h_h_dphi_20_50")
        hh_technique_dphi_dict_50_80[key] = file.Get("h_h_dphi_50_80")
        hh_technique_dphi_dict_0_80[key] = file.Get("h_h_dphi_0_80")

        if key == "gaus":

            technique_von_dict_0_20[key] = file.Get("fit_function_0_20")
            technique_von_dict_20_50[key] = file.Get("fit_function_20_50")
            technique_von_dict_50_80[key] = file.Get("fit_function_50_80")
            technique_von_dict_0_80[key] = file.Get("fit_function_0_80")
            
            hh_technique_von_dict_0_20[key] = file.Get("hh_fit_function_0_20")
            hh_technique_von_dict_20_50[key] = file.Get("hh_fit_function_20_50")
            hh_technique_von_dict_50_80[key] = file.Get("hh_fit_function_50_80")
            hh_technique_von_dict_0_80[key] = file.Get("hh_fit_function_0_80")

        elif key == "central value":

            technique_von_dict_0_20[key] = file.Get("von_fit_0_20")
            technique_von_dict_20_50[key] = file.Get("von_fit_20_50")
            technique_von_dict_50_80[key] = file.Get("von_fit_50_80")
            technique_von_dict_0_80[key] = file.Get("von_fit_0_80")

            hh_technique_von_dict_0_20[key] = file.Get("hh_von_fit_0_20")
            hh_technique_von_dict_20_50[key] = file.Get("hh_von_fit_20_50")
            hh_technique_von_dict_50_80[key] = file.Get("hh_von_fit_50_80")
            hh_technique_von_dict_0_80[key] = file.Get("hh_von_fit_0_80")

        elif key == "avg6":

            technique_von_dict_0_20[key] = get_von_fit_avgbg(technique_dphi_dict_0_20[key], PT_MODE)
            technique_von_dict_20_50[key] = get_von_fit_avgbg(technique_dphi_dict_20_50[key], PT_MODE)
            technique_von_dict_50_80[key] = get_von_fit_avgbg(technique_dphi_dict_50_80[key], PT_MODE)
            technique_von_dict_0_80[key] = get_von_fit_avgbg(technique_dphi_dict_0_80[key], PT_MODE)

            hh_technique_von_dict_0_20[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_20[key], PT_MODE)
            hh_technique_von_dict_20_50[key] = get_von_fit_avgbg(hh_technique_dphi_dict_20_50[key], PT_MODE)
            hh_technique_von_dict_50_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_50_80[key], PT_MODE)
            hh_technique_von_dict_0_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_80[key], PT_MODE)

        elif key == "avg4":

            technique_von_dict_0_20[key] = get_von_fit_avgbg(technique_dphi_dict_0_20[key], PT_MODE, False)
            technique_von_dict_20_50[key] = get_von_fit_avgbg(technique_dphi_dict_20_50[key], PT_MODE, False)
            technique_von_dict_50_80[key] = get_von_fit_avgbg(technique_dphi_dict_50_80[key], PT_MODE, False)
            technique_von_dict_0_80[key] = get_von_fit_avgbg(technique_dphi_dict_0_80[key], PT_MODE, False)

            hh_technique_von_dict_0_20[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_20[key], PT_MODE, False)
            hh_technique_von_dict_20_50[key] = get_von_fit_avgbg(hh_technique_dphi_dict_20_50[key], PT_MODE, False)
            hh_technique_von_dict_50_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_50_80[key], PT_MODE, False)
            hh_technique_von_dict_0_80[key] = get_von_fit_avgbg(hh_technique_dphi_dict_0_80[key], PT_MODE, False)



    if PT_MODE == 0:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80")

        
    elif PT_MODE == 1:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20_lowpt")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50_lowpt")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80_lowpt")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20_lowpt")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50_lowpt")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80_lowpt")

    elif PT_MODE == 2:

        ns_technique_rms_0_20, as_technique_rms_0_20 = get_rms(technique_von_dict_0_20, technique_dphi_dict_0_20, "../figures/technique_variations_width_0_20_highpt")
        ns_technique_rms_20_50, as_technique_rms_20_50 = get_rms(technique_von_dict_20_50, technique_dphi_dict_20_50, "../figures/technique_variations_width_20_50_highpt")
        ns_technique_rms_50_80, as_technique_rms_50_80 = get_rms(technique_von_dict_50_80, technique_dphi_dict_50_80, "../figures/technique_variations_width_50_80_highpt")

        ns_hh_technique_rms_0_20, as_hh_technique_rms_0_20 = get_rms(hh_technique_von_dict_0_20, hh_technique_dphi_dict_0_20, "../figures/hh_technique_variations_width_0_20_highpt")
        ns_hh_technique_rms_20_50, as_hh_technique_rms_20_50 = get_rms(hh_technique_von_dict_20_50, hh_technique_dphi_dict_20_50, "../figures/hh_technique_variations_width_20_50_highpt")
        ns_hh_technique_rms_50_80, as_hh_technique_rms_50_80 = get_rms(hh_technique_von_dict_50_80, hh_technique_dphi_dict_50_80, "../figures/hh_technique_variations_width_50_80_highpt")



    ns_total_rms_0_20 = math.sqrt(ns_signal_rms_0_20**2 + ns_sideband_rms_0_20**2 + ns_pid_rms_0_20**2 + ns_technique_rms_0_20**2)
    print("ns_total_rms_0_20: " + str(ns_total_rms_0_20))
    as_total_rms_0_20 = math.sqrt(as_signal_rms_0_20**2 + as_sideband_rms_0_20**2 + as_pid_rms_0_20**2 + as_technique_rms_0_20**2)
    print("as_total_rms_0_20: " + str(as_total_rms_0_20))

    c = rt.TCanvas("c", "c", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    ns_sideband_syst_hist = rt.TH1D("ns_sideband_syst_hist", "", 3, 0, 3)
    ns_sideband_syst_hist.SetBinContent(1, ns_sideband_rms_0_20)
    ns_sideband_syst_hist.SetBinContent(2, ns_sideband_rms_20_50)
    ns_sideband_syst_hist.SetBinContent(3, ns_sideband_rms_50_80)
    ns_sideband_syst_hist.SetBinError(1, 0)
    ns_sideband_syst_hist.SetBinError(2, 0)
    ns_sideband_syst_hist.SetBinError(3, 0)
    ns_sideband_syst_hist.SetLineColor(rt.kGreen - 2)
    ns_sideband_syst_hist.SetLineWidth(2)
    ns_signal_syst_hist = rt.TH1D("ns_signal_syst_hist", "", 3, 0, 3)
    ns_signal_syst_hist.SetBinContent(1, ns_signal_rms_0_20)
    ns_signal_syst_hist.SetBinContent(2, ns_signal_rms_20_50)
    ns_signal_syst_hist.SetBinContent(3, ns_signal_rms_50_80)
    ns_signal_syst_hist.SetBinError(1, 0)
    ns_signal_syst_hist.SetBinError(2, 0)
    ns_signal_syst_hist.SetBinError(3, 0)
    ns_signal_syst_hist.SetLineColor(rt.kMagenta - 2)
    ns_signal_syst_hist.SetLineWidth(2)
    ns_pid_syst_hist = rt.TH1D("ns_pid_syst_hist", "", 3, 0, 3)
    ns_pid_syst_hist.SetBinContent(1, ns_pid_rms_0_20)
    ns_pid_syst_hist.SetBinContent(2, ns_pid_rms_20_50)
    ns_pid_syst_hist.SetBinContent(3, ns_pid_rms_50_80)
    ns_pid_syst_hist.SetBinError(1, 0)
    ns_pid_syst_hist.SetBinError(2, 0)
    ns_pid_syst_hist.SetBinError(3, 0)
    ns_pid_syst_hist.SetLineColor(rt.kAzure - 2)
    ns_pid_syst_hist.SetLineWidth(2)
    ns_technique_syst_hist = rt.TH1D("ns_technique_syst_hist", "", 3, 0, 3)
    ns_technique_syst_hist.SetBinContent(1, ns_technique_rms_0_20)
    ns_technique_syst_hist.SetBinContent(2, ns_technique_rms_20_50)
    ns_technique_syst_hist.SetBinContent(3, ns_technique_rms_50_80)
    ns_technique_syst_hist.SetBinError(1, 0)
    ns_technique_syst_hist.SetBinError(2, 0)
    ns_technique_syst_hist.SetBinError(3, 0)
    ns_technique_syst_hist.SetLineColor(rt.kViolet + 2)
    ns_technique_syst_hist.SetLineWidth(2)


    as_sideband_syst_hist = rt.TH1D("as_sideband_syst_hist", "", 3, 0, 3)
    as_sideband_syst_hist.SetBinContent(1, as_sideband_rms_0_20)
    as_sideband_syst_hist.SetBinContent(2, as_sideband_rms_20_50)
    as_sideband_syst_hist.SetBinContent(3, as_sideband_rms_50_80)
    as_sideband_syst_hist.SetBinError(1, 0)
    as_sideband_syst_hist.SetBinError(2, 0)
    as_sideband_syst_hist.SetBinError(3, 0)
    as_sideband_syst_hist.SetLineColor(rt.kGreen - 2)
    as_sideband_syst_hist.SetLineWidth(2)
    as_signal_syst_hist = rt.TH1D("as_signal_syst_hist", "", 3, 0, 3)
    as_signal_syst_hist.SetBinContent(1, as_signal_rms_0_20)
    as_signal_syst_hist.SetBinContent(2, as_signal_rms_20_50)
    as_signal_syst_hist.SetBinContent(3, as_signal_rms_50_80)
    as_signal_syst_hist.SetBinError(1, 0)
    as_signal_syst_hist.SetBinError(2, 0)
    as_signal_syst_hist.SetBinError(3, 0)
    as_signal_syst_hist.SetLineColor(rt.kMagenta - 2)
    as_signal_syst_hist.SetLineWidth(2)
    as_pid_syst_hist = rt.TH1D("as_pid_syst_hist", "", 3, 0, 3)
    as_pid_syst_hist.SetBinContent(1, as_pid_rms_0_20)
    as_pid_syst_hist.SetBinContent(2, as_pid_rms_20_50)
    as_pid_syst_hist.SetBinContent(3, as_pid_rms_50_80)
    as_pid_syst_hist.SetBinError(1, 0)
    as_pid_syst_hist.SetBinError(2, 0)
    as_pid_syst_hist.SetBinError(3, 0)
    as_pid_syst_hist.SetLineColor(rt.kAzure - 2)
    as_pid_syst_hist.SetLineWidth(2)
    as_technique_syst_hist = rt.TH1D("as_technique_syst_hist", "", 3, 0, 3)
    as_technique_syst_hist.SetBinContent(1, as_technique_rms_0_20)
    as_technique_syst_hist.SetBinContent(2, as_technique_rms_20_50)
    as_technique_syst_hist.SetBinContent(3, as_technique_rms_50_80)
    as_technique_syst_hist.SetBinError(1, 0)
    as_technique_syst_hist.SetBinError(2, 0)
    as_technique_syst_hist.SetBinError(3, 0)
    as_technique_syst_hist.SetLineColor(rt.kViolet + 2)
    as_technique_syst_hist.SetLineWidth(2)

    # the topo syst is a LITTLE higher at low pT
    if PT_MODE == 1:
        topo_rms = 0.032
        topo_nch_dep_rms = 0.005
        topo_syst_hist = rt.TH1D("topo_syst_hist", "", 3, 0, 3)
        topo_syst_hist.SetBinContent(1, topo_rms)
        topo_syst_hist.SetBinContent(2, topo_rms)
        topo_syst_hist.SetBinContent(3, topo_rms)
        topo_syst_hist.SetBinError(1, 0)
        topo_syst_hist.SetBinError(2, 0)
        topo_syst_hist.SetBinError(3, 0)
        topo_syst_hist.SetLineColor(rt.kRed - 2)
        topo_syst_hist.SetLineWidth(2)
    else:
        topo_rms = 0.030
        topo_nch_dep_rms = 0.005
        topo_syst_hist = rt.TH1D("topo_syst_hist", "", 3, 0, 3)
        topo_syst_hist.SetBinContent(1, topo_rms)
        topo_syst_hist.SetBinContent(2, topo_rms)
        topo_syst_hist.SetBinContent(3, topo_rms)
        topo_syst_hist.SetBinError(1, 0)
        topo_syst_hist.SetBinError(2, 0)
        topo_syst_hist.SetBinError(3, 0)
        topo_syst_hist.SetLineColor(rt.kRed - 2)
        topo_syst_hist.SetLineWidth(2)

    # the matbud syst is a LITTLE higher at low pT
    if PT_MODE == 1:
        matbud_rms = 0.011
        matbud_nch_dep_rms = 0.0
        matbud_syst_hist = rt.TH1D("matbud_syst_hist", "", 3, 0, 3)
        matbud_syst_hist.SetBinContent(1, matbud_rms)
        matbud_syst_hist.SetBinContent(2, matbud_rms)
        matbud_syst_hist.SetBinContent(3, matbud_rms)
        matbud_syst_hist.SetBinError(1, 0)
        matbud_syst_hist.SetBinError(2, 0)
        matbud_syst_hist.SetBinError(3, 0)
        matbud_syst_hist.SetLineColor(rt.kOrange - 2)
        matbud_syst_hist.SetLineWidth(2)
    else:
        matbud_rms = 0.006
        matbud_nch_dep_rms = 0.000
        matbud_syst_hist = rt.TH1D("matbud_syst_hist", "", 3, 0, 3)
        matbud_syst_hist.SetBinContent(1, matbud_rms)
        matbud_syst_hist.SetBinContent(2, matbud_rms)
        matbud_syst_hist.SetBinContent(3, matbud_rms)
        matbud_syst_hist.SetBinError(1, 0)
        matbud_syst_hist.SetBinError(2, 0)
        matbud_syst_hist.SetBinError(3, 0)
        matbud_syst_hist.SetLineColor(rt.kOrange - 2)
        matbud_syst_hist.SetLineWidth(2)

    ns_total_syst_hist = rt.TH1D("ns_total_syst_hist", "", 3, 0, 3)
    ns_total_syst_hist.SetBinContent(1, math.sqrt(ns_sideband_rms_0_20**2 + ns_signal_rms_0_20**2 + ns_pid_rms_0_20**2 + ns_technique_rms_0_20**2 + topo_rms**2 + matbud_rms**2 ))
    ns_total_syst_hist.SetBinContent(2, math.sqrt(ns_sideband_rms_20_50**2 + ns_signal_rms_20_50**2 + ns_pid_rms_20_50**2 + ns_technique_rms_20_50**2 + topo_rms**2 + matbud_rms**2))
    ns_total_syst_hist.SetBinContent(3, math.sqrt(ns_sideband_rms_50_80**2 + ns_signal_rms_50_80**2 + ns_pid_rms_50_80**2 + ns_technique_rms_50_80**2 + topo_rms**2 + matbud_rms**2))
    ns_total_syst_hist.SetBinError(1, 0)
    ns_total_syst_hist.SetBinError(2, 0)
    ns_total_syst_hist.SetBinError(3, 0)
    ns_total_syst_hist.SetLineColor(rt.kBlack)
    ns_total_syst_hist.SetLineWidth(2)
    ns_total_syst_hist.SetTitle("")
    ns_total_syst_hist.GetXaxis().SetBinLabel(1, "0-20%")
    ns_total_syst_hist.GetXaxis().SetBinLabel(2, "20-50%")
    ns_total_syst_hist.GetXaxis().SetBinLabel(3, "50-80%")
    ns_total_syst_hist.GetXaxis().SetTitle("Mult. Percentile")
    ns_total_syst_hist.GetYaxis().SetTitle("Systematic Uncertainty (%)")

    as_total_syst_hist = rt.TH1D("as_total_syst_hist", "", 3, 0, 3)
    as_total_syst_hist.SetBinContent(1, math.sqrt(as_sideband_rms_0_20**2 + as_signal_rms_0_20**2 + as_pid_rms_0_20**2 + as_technique_rms_0_20**2 + topo_rms**2 + matbud_rms**2 ))
    as_total_syst_hist.SetBinContent(2, math.sqrt(as_sideband_rms_20_50**2 + as_signal_rms_20_50**2 + as_pid_rms_20_50**2 + as_technique_rms_20_50**2 + topo_rms**2 + matbud_rms**2))
    as_total_syst_hist.SetBinContent(3, math.sqrt(as_sideband_rms_50_80**2 + as_signal_rms_50_80**2 + as_pid_rms_50_80**2 + as_technique_rms_50_80**2 + topo_rms**2 + matbud_rms**2))
    as_total_syst_hist.SetBinError(1, 0)
    as_total_syst_hist.SetBinError(2, 0)
    as_total_syst_hist.SetBinError(3, 0)
    as_total_syst_hist.SetLineColor(rt.kBlack)
    as_total_syst_hist.SetLineWidth(2)
    as_total_syst_hist.SetTitle("")
    as_total_syst_hist.GetXaxis().SetBinLabel(1, "0-20%")
    as_total_syst_hist.GetXaxis().SetBinLabel(2, "20-50%")
    as_total_syst_hist.GetXaxis().SetBinLabel(3, "50-80%")
    as_total_syst_hist.GetXaxis().SetTitle("Mult. Percentile")
    as_total_syst_hist.GetYaxis().SetTitle("Systematic Uncertainty (%)")

    leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(ns_sideband_syst_hist, "Sideband", "l")
    leg.AddEntry(ns_signal_syst_hist, "Signal", "l")
    leg.AddEntry(ns_pid_syst_hist, "PID", "l")
    leg.AddEntry(ns_technique_syst_hist, "Technique", "l")
    leg.AddEntry(topo_syst_hist, "Topo", "l")
    leg.AddEntry(matbud_syst_hist, "Mat. budget", "l")
    leg.AddEntry(ns_total_syst_hist, "Total", "l")

    ns_total_syst_hist.Scale(100)
    ns_sideband_syst_hist.Scale(100)
    ns_signal_syst_hist.Scale(100)
    ns_pid_syst_hist.Scale(100)
    ns_technique_syst_hist.Scale(100)

    as_total_syst_hist.Scale(100)
    as_sideband_syst_hist.Scale(100)
    as_signal_syst_hist.Scale(100)
    as_pid_syst_hist.Scale(100)
    as_technique_syst_hist.Scale(100)

    topo_syst_hist.Scale(100)
    matbud_syst_hist.Scale(100)

    ns_total_syst_hist.GetYaxis().SetRangeUser(0, ns_total_syst_hist.GetMaximum()*1.8)
    as_total_syst_hist.GetYaxis().SetRangeUser(0, as_total_syst_hist.GetMaximum()*1.8)


    ns_total_syst_hist.Draw("hist")
    ns_sideband_syst_hist.Draw("hist same")
    ns_signal_syst_hist.Draw("hist same")
    ns_pid_syst_hist.Draw("hist same")
    ns_technique_syst_hist.Draw("hist same")
    topo_syst_hist.Draw("hist same")
    matbud_syst_hist.Draw("hist same")

    leg.Draw("SAME")
    if PT_MODE == 0:
        c.SaveAs("../figures/systematics_ns_width_postbarlow.pdf")
    elif PT_MODE == 1:
        c.SaveAs("../figures/systematics_ns_width_postbarlow_lowpt.pdf")
    elif PT_MODE == 2:
        c.SaveAs("../figures/systematics_ns_width_postbarlow_highpt.pdf")

    as_total_syst_hist.Draw("hist")
    as_sideband_syst_hist.Draw("hist same")
    as_signal_syst_hist.Draw("hist same")
    as_pid_syst_hist.Draw("hist same")
    as_technique_syst_hist.Draw("hist same")
    topo_syst_hist.Draw("hist same")
    matbud_syst_hist.Draw("hist same")

    leg.Draw("SAME")
    if PT_MODE == 0:
        c.SaveAs("../figures/systematics_as_width_postbarlow.pdf")
    elif PT_MODE == 1:
        c.SaveAs("../figures/systematics_as_width_postbarlow_lowpt.pdf")
    elif PT_MODE == 2:
        c.SaveAs("../figures/systematics_as_width_postbarlow_highpt.pdf")

    
    ns_total_syst_0_20 = ns_total_syst_hist.GetBinContent(1)/100
    ns_total_syst_20_50 = ns_total_syst_hist.GetBinContent(2)/100
    ns_total_syst_50_80 = ns_total_syst_hist.GetBinContent(3)/100

    as_total_syst_0_20 = as_total_syst_hist.GetBinContent(1)/100
    as_total_syst_20_50 = as_total_syst_hist.GetBinContent(2)/100
    as_total_syst_50_80 = as_total_syst_hist.GetBinContent(3)/100

    hh_ns_total_syst_0_20 = math.sqrt(ns_hh_technique_rms_0_20**2 + 0.035**2)
    hh_ns_total_syst_20_50 = math.sqrt(ns_hh_technique_rms_20_50**2 + 0.035**2)
    hh_ns_total_syst_50_80 = math.sqrt(ns_hh_technique_rms_50_80**2 + 0.035**2)

    hh_as_total_syst_0_20 = math.sqrt(as_hh_technique_rms_0_20**2 + 0.035**2)
    hh_as_total_syst_20_50 = math.sqrt(as_hh_technique_rms_20_50**2 + 0.035**2)
    hh_as_total_syst_50_80 = math.sqrt(as_hh_technique_rms_50_80**2 + 0.035**2)

    ns_outgraph = rt.TGraph(3, array('d', [0, 1, 2]), array('d', [ns_total_syst_50_80, ns_total_syst_20_50, ns_total_syst_0_20]))
    as_outgraph = rt.TGraph(3, array('d', [0, 1, 2]), array('d', [as_total_syst_50_80, as_total_syst_20_50, as_total_syst_0_20]))

    hh_ns_outgraph = rt.TGraph(3, array('d', [0, 1, 2]), array('d', [hh_ns_total_syst_50_80, hh_ns_total_syst_20_50, hh_ns_total_syst_0_20]))
    hh_as_outgraph = rt.TGraph(3, array('d', [0, 1, 2]), array('d', [hh_as_total_syst_50_80, hh_as_total_syst_20_50, hh_as_total_syst_0_20]))

    outfile.cd()
    if PT_MODE == 0:
        ns_outgraph.Write("ns_total_syst")
        as_outgraph.Write("as_total_syst")
        hh_ns_outgraph.Write("hh_ns_total_syst")
        hh_as_outgraph.Write("hh_as_total_syst")
    elif PT_MODE == 1:
        ns_outgraph.Write("ns_total_syst_lowpt")
        as_outgraph.Write("as_total_syst_lowpt")
        hh_ns_outgraph.Write("hh_ns_total_syst_lowpt")
        hh_as_outgraph.Write("hh_as_total_syst_lowpt")
    elif PT_MODE == 2:
        ns_outgraph.Write("ns_total_syst_highpt")
        as_outgraph.Write("as_total_syst_highpt")
        hh_ns_outgraph.Write("hh_ns_total_syst_highpt")
        hh_as_outgraph.Write("hh_as_total_syst_highpt")
import ROOT as rt

colors = [rt.kGreen + 2, rt.kRed + 2, rt.kBlue + 2,  rt.kMagenta + 2,  rt.kYellow + 2,  rt.kViolet + 2, rt.kOrange + 2]

def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error

def get_rms(von_dict, dphi_dict, output_name, pid=False):

    c1 = rt.TCanvas("c1", "c1", 800, 600)
    c1.SetLeftMargin(0.15)

    color_index = 0
    leg = rt.TLegend(0.67, 0.67, 0.8, 0.8)
    leg.SetBorderSize(0)
    for key in von_dict:
        if key == "central value":

            central_fit = von_dict[key]
            central_fit.SetNpx(1000)
            central_fit.GetXaxis().SetRangeUser(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            central_fit.SetLineColor(rt.kBlack)
            central_fit.SetLineWidth(2)
            central_fit.SetTitle("")
            central_fit.GetYaxis().SetRangeUser(0.7*central_fit.GetMinimum(), 1.3*central_fit.GetMaximum())
            central_fit.GetYaxis().SetTitle("Counts")
            central_fit.GetXaxis().SetTitle("#Delta#varphi")
            central_fit.Draw()

            central_ns_width = get_width_from_kappa(central_fit.GetParameter(4))
            centra_ns_width_error = get_width_error_from_kappa(central_fit.GetParameter(4), central_fit.GetParError(4))
            central_as_width = get_width_from_kappa(central_fit.GetParameter(5))
            central_as_width_error = get_width_error_from_kappa(central_fit.GetParameter(5), central_fit.GetParError(5))

            central_plot = dphi_dict[key]
            central_plot.GetXaxis().SetRangeUser(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
            central_plot.SetLineColor(rt.kBlack)
            central_plot.SetLineWidth(2)
            central_plot.GetYaxis().SetRangeUser(0.7*central_plot.GetMinimum(), 1.3*central_plot.GetMaximum())
            central_plot.GetYaxis().SetTitle("Counts")
            central_plot.GetXaxis().SetTitle("#Delta#varphi")
            central_plot.Draw("SAME")

            leg.AddEntry(central_plot, key, "l")
        else:
            varied_fit = von_dict[key]
            varied_fit.SetNpx(1000)
            if pid:
                varied_fit.Scale(central_fit.Integral()/varied_fit.Integral())
            varied_fit.SetLineColor(colors[color_index])
            varied_fit.SetLineWidth(2)
            varied_fit.Draw("same")

            varied_ns_width = get_width_from_kappa(varied_fit.GetParameter(4))
            varied_ns_width_error = get_width_error_from_kappa(varied_fit.GetParameter(4), varied_fit.GetParError(4))
            varied_as_width = get_width_from_kappa(varied_fit.GetParameter(5))
            varied_as_width_error = get_width_error_from_kappa(varied_fit.GetParameter(5), varied_fit.GetParError(5))

            varied_plot = dphi_dict[key]
            if pid:
                varied_plot.Scale(central_plot.Integral()/varied_plot.Integral())
            varied_plot.SetLineColor(colors[color_index])
            varied_plot.SetLineWidth(2)
            varied_plot.Draw("same")

            leg.AddEntry(varied_plot, key, "l")
        color_index += 1
    
    leg.Draw("SAME")
    c1.SaveAs(output_name + ".pdf")

    return 1, 1



for PT_MODE in range(3):

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
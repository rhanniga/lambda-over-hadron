import ROOT as rt
rt.gStyle.SetOptStat(0)

import math

for PT_MODE in range(3):
    BARLOW = False

    colors = [rt.kGreen + 2, rt.kRed + 2, rt.kBlue + 2,  rt.kMagenta + 2,  rt.kYellow + 2,  rt.kViolet + 2, rt.kOrange + 2]


    def get_nch_dep_rms(dphi_dict, dphi_dict_minbias, output_string, pid=False):

        c1 = rt.TCanvas("c1","c1",800,600)
        c1.SetLeftMargin(0.15)

        nch_dep_rms = 0
        n = 0

        leg = rt.TLegend(0.67, 0.67, 0.8, 0.8)
        leg.SetBorderSize(0)

        color_index = 0

        pointer_list = []

        for key in dphi_dict:
            if key == "central value":

                central_plot = dphi_dict[key].Clone(f"central_plot{color_index}_{key}")
                central_ratio = central_plot.Clone(f"central_ratio{color_index}_{key}")
                central_ratio.SetTitle("")
                central_ratio.GetYaxis().SetTitle("R")
                central_ratio.GetXaxis().SetTitle("#Delta#varphi")
                central_ratio.GetYaxis().SetRangeUser(0.7, 1.3)
                central_ratio.Divide(central_plot)
                central_ratio.SetLineColor(rt.kBlack)
                central_ratio.SetLineWidth(2)

                central_plot_minbias = dphi_dict_minbias[key].Clone(f"central_plot_minbias{color_index}")
                central_ratio_minbias = central_plot_minbias.Clone(f"central_ratio_minbias{color_index}")
                central_ratio_minbias.SetTitle("")
                central_ratio_minbias.GetYaxis().SetTitle("R")
                central_ratio_minbias.GetXaxis().SetTitle("#Delta#varphi")
                central_ratio_minbias.GetYaxis().SetRangeUser(0.7, 1.3)
                central_ratio_minbias.Divide(central_ratio_minbias)
                central_ratio_minbias.SetLineColor(rt.kBlack)
                central_ratio_minbias.SetLineWidth(2)

                central_double_ratio = central_ratio/central_ratio_minbias
                main_clone = central_double_ratio.Clone(f"main_clone{color_index}")
                pointer_list.append(main_clone)
                main_clone.Draw("SAME")


            else:
                varied_plot = dphi_dict[key].Clone(f"varied_plot{color_index}")
                ratio = varied_plot.Clone(f"ratio{color_index}")
                if pid:
                    ratio.Scale(central_plot.Integral()/ratio.Integral())
                ratio.Divide(central_plot)
                ratio.SetLineColor(colors[color_index])
                ratio.SetLineWidth(2)

                varied_plot_minbias = dphi_dict_minbias[key].Clone(f"varied_plot_minbias{color_index}")
                ratio_minbias = varied_plot_minbias.Clone(f"ratio_minbias{color_index}")
                if pid:
                    ratio_minbias.Scale(central_plot_minbias.Integral()/ratio_minbias.Integral())
                ratio_minbias.Divide(central_plot_minbias)
                ratio_minbias.SetLineColor(colors[color_index])
                ratio_minbias.SetLineWidth(2)

                double_ratio = ratio/ratio_minbias
                double_ratio.SetLineColor(colors[color_index])
                double_ratio.SetLineWidth(2)
                baby_clone = double_ratio.Clone(f"baby_clone{color_index}")
                pointer_list.append(baby_clone)
                baby_clone.Draw("SAME")
                leg.AddEntry(baby_clone, key, "l")

                for i in range(1, double_ratio.GetNbinsX()+1):
                    nch_dep_rms += (double_ratio.GetBinContent(i) - 1)**2
                    n += 1


            color_index += 1
        

        leg.Draw("SAME")
        c1.Draw()
        c1.SaveAs(output_string + "_ratio.pdf")

        nch_dep_rms = math.sqrt(nch_dep_rms/n)
        
        return nch_dep_rms



    def get_rms(dphi_dict, output_string, pid=False):

        c1 = rt.TCanvas("c1","c1",800,600)
        c1.SetLeftMargin(0.15)

        color_index = 0
        leg = rt.TLegend(0.67, 0.67, 0.8, 0.8)
        leg.SetBorderSize(0)
        for key in dphi_dict:
            if key == "central value":
                central_plot = dphi_dict[key]
                central_plot.SetLineColor(rt.kBlack)
                central_plot.SetLineWidth(2)
                central_plot.GetYaxis().SetRangeUser(0.7*central_plot.GetMinimum(), 1.3*central_plot.GetMaximum())
                central_plot.SetTitle("")
                central_plot.GetYaxis().SetTitle("Counts")
                central_plot.GetXaxis().SetTitle("#Delta#varphi")
                central_plot.Draw()
                leg.AddEntry(central_plot, key, "l")
            else:
                varied_plot = dphi_dict[key]
                if pid:
                    varied_plot.Scale(central_plot.Integral()/varied_plot.Integral()) 
                varied_plot.SetLineColor(colors[color_index])
                varied_plot.SetLineWidth(2)
                varied_plot.Draw("same")
                leg.AddEntry(varied_plot, key, "l")
            color_index += 1

        leg.Draw("SAME")
        c1.Draw()
        c1.SaveAs(output_string + ".pdf")

        c2 = rt.TCanvas("c2","c2",800,600)
        c2.SetLeftMargin(0.15)

        color_index = 0
        rms = 0
        n = 0

        pointer_list = []
        for key in dphi_dict:
            if key == "central value":
                central_plot = dphi_dict[key].Clone(f"central_plot_2{color_index}_{key}")
                central_ratio = central_plot.Clone(f"central_ratio_2{color_index}_{key}")
                central_ratio.SetTitle("")
                central_ratio.GetYaxis().SetTitle("Ratio")
                central_ratio.GetXaxis().SetTitle("#Delta#varphi")
                central_ratio.GetYaxis().SetRangeUser(0.7, 1.3)
                central_ratio.Divide(central_plot)
                central_ratio.SetLineColor(rt.kBlack)
                central_ratio.SetLineWidth(2)
                central_ratio.Draw("SAME")
            else:
                varied_plot = dphi_dict[key].Clone(f"varied_plot_2{color_index}")
                if pid:
                    varied_plot.Scale(central_plot.Integral()/varied_plot.Integral())
                ratio = varied_plot/central_plot
                ratio.SetName(f"ratio_2{color_index}")
                ratio.SetLineColor(colors[color_index])
                ratio.SetLineWidth(2)
                clone = ratio.Clone(f"fuck_this_{color_index}")
                pointer_list.append(clone)
                clone.Draw("SAME")
                
                for i in range(1, ratio.GetNbinsX()+1):
                    rms += (ratio.GetBinContent(i) - 1)**2
                    n += 1
            color_index += 1
        
        rms = math.sqrt(rms/n)

        leg.Draw("SAME")
        c2.Draw()
        c2.SaveAs(output_string + "_ratio.pdf")

        return rms



    if PT_MODE == 0: # CENTRAL 2 - 4 BIN
        signal_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrow_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_narrower_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_wide_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        signal_wider_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1: # LOW 1.5 - 2.5 BIN
        signal_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrow_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_narrower_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_wide_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        signal_wider_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2: # HIGH 2.5 - 4 BIN
        signal_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrow_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1108_1124_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_narrower_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1112_112_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_wide_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_11_1132_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        signal_wider_file = rt.TFile.Open("output/signal_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1096_1136_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    signal_file_dict = {
        "central value": signal_central_file,
        "narrow": signal_narrow_file,
        "wide": signal_wide_file,
        "narrower": signal_narrower_file,
        "wider": signal_wider_file
    }

    signal_dphi_dict_0_20 = {}
    signal_dphi_dict_20_50 = {}
    signal_dphi_dict_50_80 = {}
    signal_dphi_dict_0_80 = {}
    for key, file in signal_file_dict.items():
        signal_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        signal_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        signal_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        signal_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")
    if PT_MODE == 0:
        signal_rms_0_20 = get_rms(signal_dphi_dict_0_20, "figures/signal_variations_dphi_0_20")
        signal_rms_20_50 = get_rms(signal_dphi_dict_20_50, "figures/signal_variations_dphi_20_50")
        signal_rms_50_80 = get_rms(signal_dphi_dict_50_80, "figures/signal_variations_dphi_50_80")
        signal_nch_dep_rms_0_20 = get_nch_dep_rms(signal_dphi_dict_0_20, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_0_20")
        signal_nch_dep_rms_20_50 = get_nch_dep_rms(signal_dphi_dict_20_50, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_20_50")
        signal_nch_dep_rms_50_80 = get_nch_dep_rms(signal_dphi_dict_50_80, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_50_80")
    elif PT_MODE == 1:
        signal_rms_0_20 = get_rms(signal_dphi_dict_0_20, "figures/signal_variations_dphi_0_20_lowpt")
        signal_rms_20_50 = get_rms(signal_dphi_dict_20_50, "figures/signal_variations_dphi_20_50_lowpt")
        signal_rms_50_80 = get_rms(signal_dphi_dict_50_80, "figures/signal_variations_dphi_50_80_lowpt")
        signal_nch_dep_rms_0_20 = get_nch_dep_rms(signal_dphi_dict_0_20, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_0_20_lowpt")
        signal_nch_dep_rms_20_50 = get_nch_dep_rms(signal_dphi_dict_20_50, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_20_50_lowpt")
        signal_nch_dep_rms_50_80 = get_nch_dep_rms(signal_dphi_dict_50_80, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_50_80_lowpt")
    elif PT_MODE == 2:
        signal_rms_0_20 = get_rms(signal_dphi_dict_0_20, "figures/signal_variations_dphi_0_20_highpt")
        signal_rms_20_50 = get_rms(signal_dphi_dict_20_50, "figures/signal_variations_dphi_20_50_highpt")
        signal_rms_50_80 = get_rms(signal_dphi_dict_50_80, "figures/signal_variations_dphi_50_80_highpt")
        signal_nch_dep_rms_0_20 = get_nch_dep_rms(signal_dphi_dict_0_20, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_0_20_highpt")
        signal_nch_dep_rms_20_50 = get_nch_dep_rms(signal_dphi_dict_20_50, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_20_50_highpt")
        signal_nch_dep_rms_50_80 = get_nch_dep_rms(signal_dphi_dict_50_80, signal_dphi_dict_0_80, "figures/signal_nch_dep_rms_dphi_50_80_highpt")



    if PT_MODE == 0:
        sideband_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_right_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_shifted_left_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_wide_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        sideband_narrow_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    elif PT_MODE == 1:
        sideband_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_right_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_shifted_left_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_wide_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        sideband_narrow_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
    elif PT_MODE == 2:
        sideband_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_right_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_114_1155_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_shifted_left_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1084_1096_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_wide_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_116_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        sideband_narrow_file = rt.TFile.Open("output/sideband_variation/v0_avg6_sideband_subtraction_rsb_1135_1145_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

    sideband_file_dict = {
        "central value": sideband_central_file,
        "narrow": sideband_narrow_file,
        "wide": sideband_wide_file,
        "shifted right": sideband_shifted_right_file,
        "shifted left": sideband_shifted_left_file
    }

    sideband_dphi_dict_0_20 = {}
    sideband_dphi_dict_20_50 = {}
    sideband_dphi_dict_50_80 = {}
    sideband_dphi_dict_0_80 = {}
    for key, file in sideband_file_dict.items():
        sideband_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        sideband_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        sideband_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        sideband_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")
    if PT_MODE == 0:
        sideband_rms_0_20 = get_rms(sideband_dphi_dict_0_20, "figures/sideband_variations_dphi_0_20")
        sideband_rms_20_50 = get_rms(sideband_dphi_dict_20_50, "figures/sideband_variations_dphi_20_50")
        sideband_rms_50_80 = get_rms(sideband_dphi_dict_50_80, "figures/sideband_variations_dphi_50_80")
        sideband_nch_dep_rms_0_20 = get_nch_dep_rms(sideband_dphi_dict_0_20, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_0_20")
        sideband_nch_dep_rms_20_50 = get_nch_dep_rms(sideband_dphi_dict_20_50, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_20_50")
        sideband_nch_dep_rms_50_80 = get_nch_dep_rms(sideband_dphi_dict_50_80, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_50_80")
    elif PT_MODE == 1:
        sideband_rms_0_20 = get_rms(sideband_dphi_dict_0_20, "figures/sideband_variations_dphi_0_20_lowpt")
        sideband_rms_20_50 = get_rms(sideband_dphi_dict_20_50, "figures/sideband_variations_dphi_20_50_lowpt")
        sideband_rms_50_80 = get_rms(sideband_dphi_dict_50_80, "figures/sideband_variations_dphi_50_80_lowpt")
        sideband_nch_dep_rms_0_20 = get_nch_dep_rms(sideband_dphi_dict_0_20, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_0_20_lowpt")
        sideband_nch_dep_rms_20_50 = get_nch_dep_rms(sideband_dphi_dict_20_50, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_20_50_lowpt")
        sideband_nch_dep_rms_50_80 = get_nch_dep_rms(sideband_dphi_dict_50_80, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_50_80_lowpt")
    elif PT_MODE == 2:
        sideband_rms_0_20 = get_rms(sideband_dphi_dict_0_20, "figures/sideband_variations_dphi_0_20_highpt")
        sideband_rms_20_50 = get_rms(sideband_dphi_dict_20_50, "figures/sideband_variations_dphi_20_50_highpt")
        sideband_rms_50_80 = get_rms(sideband_dphi_dict_50_80, "figures/sideband_variations_dphi_50_80_highpt")
        sideband_nch_dep_rms_0_20 = get_nch_dep_rms(sideband_dphi_dict_0_20, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_0_20_highpt")
        sideband_nch_dep_rms_20_50 = get_nch_dep_rms(sideband_dphi_dict_20_50, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_20_50_highpt")
        sideband_nch_dep_rms_50_80 = get_nch_dep_rms(sideband_dphi_dict_50_80, sideband_dphi_dict_0_80, "figures/sideband_nch_dep_rms_dphi_50_80_highpt")


    if PT_MODE == 0:
        pid_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        pid_narrow_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_narrow.root")
        pid_wide_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_wide.root")
        pid_tof_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_tof.root")

    elif PT_MODE == 1:
        pid_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        pid_narrow_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_narrow_lowpt.root")
        pid_wide_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_wide_lowpt.root")
        pid_tof_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_tof_lowpt.root")
    elif PT_MODE == 2:
        pid_central_file = rt.TFile.Open("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        pid_narrow_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_narrow_highpt.root")
        pid_wide_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_wide_highpt.root")
        pid_tof_file = rt.TFile.Open("output/pid_variation/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_tof_highpt.root")

    pid_file_dict = {
        "central value": pid_central_file,
        "narrow": pid_narrow_file,
        "wide": pid_wide_file,
        "require tof": pid_tof_file
    }

    pid_dphi_dict_0_20 = {}
    pid_dphi_dict_20_50 = {}
    pid_dphi_dict_50_80 = {}
    pid_dphi_dict_0_80 = {}
    for key, file in pid_file_dict.items():
        pid_dphi_dict_0_20[key] = file.Get("h_lambda_dphi_subtracted_0_20")
        pid_dphi_dict_20_50[key] = file.Get("h_lambda_dphi_subtracted_20_50")
        pid_dphi_dict_50_80[key] = file.Get("h_lambda_dphi_subtracted_50_80")
        pid_dphi_dict_0_80[key] = file.Get("h_lambda_dphi_subtracted_0_80")
    if PT_MODE == 0:
        pid_rms_0_20 = get_rms(pid_dphi_dict_0_20, "figures/pid_variations_dphi_0_20", True)
        pid_rms_20_50 = get_rms(pid_dphi_dict_20_50, "figures/pid_variations_dphi_20_50", True)
        pid_rms_50_80 = get_rms(pid_dphi_dict_50_80, "figures/pid_variations_dphi_50_80", True)
        pid_nch_dep_rms_0_20 = get_nch_dep_rms(pid_dphi_dict_0_20, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_0_20", True)
        pid_nch_dep_rms_20_50 = get_nch_dep_rms(pid_dphi_dict_20_50, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_20_50", True)
        pid_nch_dep_rms_50_80 = get_nch_dep_rms(pid_dphi_dict_50_80, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_50_80", True)
    elif PT_MODE == 1:
        pid_rms_0_20 = get_rms(pid_dphi_dict_0_20, "figures/pid_variations_dphi_0_20_lowpt", True)
        pid_rms_20_50 = get_rms(pid_dphi_dict_20_50, "figures/pid_variations_dphi_20_50_lowpt", True)
        pid_rms_50_80 = get_rms(pid_dphi_dict_50_80, "figures/pid_variations_dphi_50_80_lowpt", True)
        pid_nch_dep_rms_0_20 = get_nch_dep_rms(pid_dphi_dict_0_20, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_0_20_lowpt", True)
        pid_nch_dep_rms_20_50 = get_nch_dep_rms(pid_dphi_dict_20_50, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_20_50_lowpt", True)
        pid_nch_dep_rms_50_80 = get_nch_dep_rms(pid_dphi_dict_50_80, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_50_80_lowpt", True)
    elif PT_MODE == 2:
        pid_rms_0_20 = get_rms(pid_dphi_dict_0_20, "figures/pid_variations_dphi_0_20_highpt", True)
        pid_rms_20_50 = get_rms(pid_dphi_dict_20_50, "figures/pid_variations_dphi_20_50_highpt", True)
        pid_rms_50_80 = get_rms(pid_dphi_dict_50_80, "figures/pid_variations_dphi_50_80_highpt", True)
        pid_nch_dep_rms_0_20 = get_nch_dep_rms(pid_dphi_dict_0_20, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_0_20_highpt", True)
        pid_nch_dep_rms_20_50 = get_nch_dep_rms(pid_dphi_dict_20_50, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_20_50_highpt", True)
        pid_nch_dep_rms_50_80 = get_nch_dep_rms(pid_dphi_dict_50_80, pid_dphi_dict_0_80, "figures/pid_nch_dep_rms_dphi_50_80_highpt", True)


    c = rt.TCanvas("c", "c", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)

    sideband_syst_hist = rt.TH1D("sideband_syst_hist", "", 3, 0, 3)
    sideband_syst_hist.SetBinContent(1, sideband_rms_0_20)
    sideband_syst_hist.SetBinContent(2, sideband_rms_20_50)
    sideband_syst_hist.SetBinContent(3, sideband_rms_50_80)
    sideband_syst_hist.SetBinError(1, 0)
    sideband_syst_hist.SetBinError(2, 0)
    sideband_syst_hist.SetBinError(3, 0)
    sideband_syst_hist.SetLineColor(rt.kGreen - 2)
    sideband_syst_hist.SetLineWidth(2)

    signal_syst_hist = rt.TH1D("signal_syst_hist", "", 3, 0, 3)
    signal_syst_hist.SetBinContent(1, signal_rms_0_20)
    signal_syst_hist.SetBinContent(2, signal_rms_20_50)
    signal_syst_hist.SetBinContent(3, signal_rms_50_80)
    signal_syst_hist.SetBinError(1, 0)
    signal_syst_hist.SetBinError(2, 0)
    signal_syst_hist.SetBinError(3, 0)
    signal_syst_hist.SetLineColor(rt.kMagenta - 2)
    signal_syst_hist.SetLineWidth(2)

    pid_syst_hist = rt.TH1D("pid_syst_hist", "", 3, 0, 3)
    pid_syst_hist.SetBinContent(1, pid_rms_0_20)
    pid_syst_hist.SetBinContent(2, pid_rms_20_50)
    pid_syst_hist.SetBinContent(3, pid_rms_50_80)
    pid_syst_hist.SetBinError(1, 0)
    pid_syst_hist.SetBinError(2, 0)
    pid_syst_hist.SetBinError(3, 0)
    pid_syst_hist.SetLineColor(rt.kAzure - 2)
    pid_syst_hist.SetLineWidth(2)

    # the topo syst is a LITTLE higher at low pT
    if PT_MODE == 1:
        topo_rms = 0.032
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
        topo_syst_hist = rt.TH1D("topo_syst_hist", "", 3, 0, 3)
        topo_syst_hist.SetBinContent(1, topo_rms)
        topo_syst_hist.SetBinContent(2, topo_rms)
        topo_syst_hist.SetBinContent(3, topo_rms)
        topo_syst_hist.SetBinError(1, 0)
        topo_syst_hist.SetBinError(2, 0)
        topo_syst_hist.SetBinError(3, 0)
        topo_syst_hist.SetLineColor(rt.kRed - 2)
        topo_syst_hist.SetLineWidth(2)


    total_syst_hist = rt.TH1D("total_syst_hist", "", 3, 0, 3)
    total_syst_hist.SetBinContent(1, math.sqrt(sideband_rms_0_20**2 + signal_rms_0_20**2 + pid_rms_0_20**2 + topo_rms**2))
    total_syst_hist.SetBinContent(2, math.sqrt(sideband_rms_20_50**2 + signal_rms_20_50**2 + pid_rms_20_50**2 + topo_rms**2))
    total_syst_hist.SetBinContent(3, math.sqrt(sideband_rms_50_80**2 + signal_rms_50_80**2 + pid_rms_50_80**2 + topo_rms**2))
    total_syst_hist.SetBinError(1, 0)
    total_syst_hist.SetBinError(2, 0)
    total_syst_hist.SetBinError(3, 0)
    total_syst_hist.SetLineColor(rt.kBlack)
    total_syst_hist.SetLineWidth(2)

    total_syst_hist.SetTitle("")
    total_syst_hist.GetXaxis().SetBinLabel(1, "0-20%")
    total_syst_hist.GetXaxis().SetBinLabel(2, "20-50%")
    total_syst_hist.GetXaxis().SetBinLabel(3, "50-80%")
    total_syst_hist.GetXaxis().SetTitle("Mult. Percentile")
    total_syst_hist.GetYaxis().SetTitle("Systematic Uncertainty (%)")

    leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(sideband_syst_hist, "Sideband", "l")
    leg.AddEntry(signal_syst_hist, "Signal", "l")
    leg.AddEntry(pid_syst_hist, "PID", "l")
    leg.AddEntry(topo_syst_hist, "Topo", "l")
    leg.AddEntry(total_syst_hist, "Total", "l")

    total_syst_hist.Scale(100)
    sideband_syst_hist.Scale(100)
    signal_syst_hist.Scale(100)
    pid_syst_hist.Scale(100)
    topo_syst_hist.Scale(100)

    total_syst_hist.GetYaxis().SetRangeUser(0, 9)


    total_syst_hist.Draw("hist")
    sideband_syst_hist.Draw("hist same")
    signal_syst_hist.Draw("hist same")
    pid_syst_hist.Draw("hist same")
    topo_syst_hist.Draw("hist same")
    leg.Draw("SAME")
    if PT_MODE == 0:
        c.SaveAs("figures/systematics_dphi_prebarlow.pdf")
    elif PT_MODE == 1:
        c.SaveAs("figures/systematics_dphi_prebarlow_lowpt.pdf")
    elif PT_MODE == 2:
        c.SaveAs("figures/systematics_dphi_prebarlow_highpt.pdf")


    # print a latex table with the systematics (signal, sideband, pid, topo, total)
    if PT_MODE == 0:
        print("---------------------CENTRAL PT---------------------")
    elif PT_MODE == 1:
        print("---------------------LOW PT---------------------")
    elif PT_MODE == 2:
        print("---------------------HIGH PT---------------------")
    print(f"0-20\% & {signal_rms_0_20*100:.2e} $\pm$ {sideband_rms_0_20*100:.2e} & {pid_rms_0_20*100:.2e} & {topo_rms*100:.2e} & {math.sqrt(sideband_rms_0_20**2 + signal_rms_0_20**2 + pid_rms_0_20**2 + topo_rms**2)*100:.2e} \\\\")
    print(f"20-50\% & {signal_rms_20_50*100:.2e} $\pm$ {sideband_rms_20_50*100:.2e} & {pid_rms_20_50*100:.2e} & {topo_rms*100:.2e} & {math.sqrt(sideband_rms_20_50**2 + signal_rms_20_50**2 + pid_rms_20_50**2 + topo_rms**2)*100:.2e} \\\\")
    print(f"50-80\% & {signal_rms_50_80*100:.2e} $\pm$ {sideband_rms_50_80*100:.2e} & {pid_rms_50_80*100:.2e} & {topo_rms*100:.2e} & {math.sqrt(sideband_rms_50_80**2 + signal_rms_50_80**2 + pid_rms_50_80**2 + topo_rms**2)*100:.2e} \\\\")

import ROOT as rt
rt.gStyle.SetOptStat(0)

import math

colors = [rt.kGreen + 2, rt.kRed + 2, rt.kBlue + 2,  rt.kMagenta + 2,  rt.kOrange + 2]

TRIGGER_V2 = 0.092
LAMBDA_V2 = 0.111
ASSOCIATED_V2 = 0.112

LAMBDA_V2_LOWPT = 0.075
ASSOCIATED_V2_LOWPT = 0.100

LAMBDA_V2_HIGHPT = 0.137
ASSOCIATED_V2_HIGHPT = 0.119

def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error

def get_fit(dphi_dist, bin_string):

    fit_dphi_dist = dphi_dist.Clone(f"{bin_string}_fit_dist")

    v2_fitpars_file_string = "../v2_fitpars.txt"
    # the exact ordering is: key: pedestal, pedestal error, trigger v2 (fixed to weighted avg), associated v2 (fixed to weighted avg)

    v2_fitpars_dict = {
        line[0] : [float(line[1]), float(line[2]), float(line[3]), float(line[4])] for line in [line.split() for line in open(v2_fitpars_file_string).readlines()]
    }

    left_v2_fit = rt.TF1("left_v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, -1)
    right_v2_fit = rt.TF1("right_v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", 1, rt.TMath.Pi()/2)
    total_v2_fit = rt.TF1("v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

    if "h_h_" in bin_string:
        if "cent_0_20" in bin_string:
            left_v2_fit.FixParameter(1, TRIGGER_V2)
            right_v2_fit.FixParameter(1, TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2)
            
        elif "cent_20_50" in bin_string:
            left_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2)

        elif "cent_50_80" in bin_string:
            left_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2)

        elif "cent_0_80" in bin_string:
            left_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2)

    elif "h_lambda_" in bin_string:

        if "cent_0_20" in bin_string:
            left_v2_fit.FixParameter(1, TRIGGER_V2)
            right_v2_fit.FixParameter(1, TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, LAMBDA_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, LAMBDA_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, LAMBDA_V2)
                right_v2_fit.FixParameter(2, LAMBDA_V2)


        elif "cent_20_50" in bin_string:
            left_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2)

        elif "cent_50_80" in bin_string:
            left_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2)

        elif "cent_0_80" in bin_string:
            left_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)

            if "assoc_15_25" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_LOWPT)
            elif "assoc_25_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_HIGHPT)
            elif "assoc_2_4" in bin_string:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2)


    left_v2_fit.SetParameter(0, 0.5*dphi_dist.GetMaximum())
    right_v2_fit.SetParameter(0, 0.5*dphi_dist.GetMaximum())

    fit_dphi_dist.Fit(left_v2_fit, "R")
    fit_dphi_dist.Fit(right_v2_fit, "R+")


    total_v2_fit.SetParameter(0, (left_v2_fit.GetParameter(0)+right_v2_fit.GetParameter(0))/2)
    total_v2_fit.SetParameter(1, (left_v2_fit.GetParameter(1)+right_v2_fit.GetParameter(1))/2)
    total_v2_fit.SetParameter(2, (left_v2_fit.GetParameter(2)+right_v2_fit.GetParameter(2))/2)

    von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
    von_fit = rt.TF1("von_fit", von_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

    von_fit.FixParameter(0, total_v2_fit.GetParameter(0))
    von_fit.FixParameter(1, total_v2_fit.GetParameter(1))
    von_fit.FixParameter(2, total_v2_fit.GetParameter(2))

    if "assoc_2_4" in bin_string:
        # von_fit.SetParLimits(3, 0.0001, 1.1)
        von_fit.SetParameter(3, 0.015)
        # von_fit.SetParLimits(4, 0.8, 20)
        von_fit.SetParameter(4, 7.6)
        # von_fit.SetParLimits(5, 0.0001, 1.1)
        von_fit.SetParameter(5, 0.015)
        # von_fit.SetParLimits(6, 0.8, 20)
        von_fit.SetParameter(6, 2.2)
    elif "assoc_15_25" in bin_string:
        # von_fit.SetParLimits(3, 0.02/3, 0.02*3)
        von_fit.SetParameter(3, 0.0137579)
        # von_fit.SetParLimits(4, 4.6/3, 4.6*3)
        von_fit.SetParameter(4, 8.0)
        # von_fit.SetParLimits(5, 0.057/5, 0.057*5)
        von_fit.SetParameter(5, 1.33787e-2)
        # von_fit.SetParLimits(6, 2/4, 8)
        von_fit.SetParameter(6, 2.01496)
    elif "assoc_25_4" in bin_string:
        # von_fit.SetParLimits(3, 0.009/3, 0.009*3)
        von_fit.SetParameter(3, 0.009)
        # von_fit.SetParLimits(4, 8.88/3, 8.88*3)
        von_fit.SetParameter(4, 8.88)
        # von_fit.SetParLimits(5, 0.0076/5, 0.0076*5)
        von_fit.SetParameter(5, 0.0076)
        # von_fit.SetParLimits(6, 2.55/3, 2.55*3)
        von_fit.SetParameter(6, 2.55)

    if "assoc_15_25" in bin_string:
        fit_dphi_dist.Fit(von_fit, "R", "", -rt.TMath.Pi()/2, rt.TMath.Pi() - 0.0001)
    else:
        fit_dphi_dist.Fit(von_fit, "R")

    return von_fit

def get_random_variation(dphi_dist, sys_error):

    random_dist = dphi_dist.Clone(dphi_dist.GetName() + "_randomized")

    for bin in range(1, dphi_dist.GetNbinsX() + 1):
        low_bound = dphi_dist.GetBinContent(bin) - (sys_error*dphi_dist.GetBinContent(bin))/3
        high_bound = dphi_dist.GetBinContent(bin) + (sys_error*dphi_dist.GetBinContent(bin))/3
        random_dist.SetBinContent(bin, rt.gRandom.Uniform(low_bound, high_bound))
    
    return random_dist

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

PT_BINS = ["assoc_15_25", "assoc_25_4"]

CENT_BINS = ["_20_50"]


for PT_BIN in PT_BINS:

    if PT_BIN == "assoc_2_4":
        infile = rt.TFile("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
        h_lambda_sys_error = 0.035
    elif PT_BIN == "assoc_15_25":
        infile = rt.TFile("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
        h_lambda_sys_error = 0.035
    elif PT_BIN == "assoc_25_4":
        infile = rt.TFile("../output/central_value/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")
        h_lambda_sys_error = 0.035

    
    for CENT_BIN in CENT_BINS:


        h_lambda_dphi_dict = {}
        h_lambda_von_dict = {}

        h_lambda_dphi_nominal = infile.Get("h_lambda_dphi_subtracted" + CENT_BIN).Clone()
        h_lambda_von_nominal = infile.Get("von_fit" + CENT_BIN).Clone()

        # h_lambda_dphi_nominal = infile.Get("h_h_dphi" + CENT_BIN).Clone()
        # h_lambda_von_nominal = infile.Get("hh_von_fit" + CENT_BIN).Clone()

        h_lambda_dphi_dict["central value"] = h_lambda_dphi_nominal
        h_lambda_von_dict["central value"] = h_lambda_von_nominal

        for trial in range(1, 5):

            h_lambda_dphi_var = get_random_variation(h_lambda_dphi_nominal, h_lambda_sys_error)
            h_lambda_von_var = get_fit(h_lambda_dphi_var, "h_lambda_cent" + CENT_BIN + "_trigger_4_8_" + PT_BIN)
            h_lambda_dphi_dict[f"variation {trial}"] = h_lambda_dphi_var
            h_lambda_von_dict[f"variation {trial}"] = h_lambda_von_var

        ns_rms, as_rms = get_rms(h_lambda_von_dict, h_lambda_dphi_dict, "h_lambda_cent" + CENT_BIN + "_trigger_4_8_" + PT_BIN + "_topo_variations")

        print(ns_rms, as_rms)
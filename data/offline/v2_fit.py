import ROOT as rt

rt.gStyle.SetOptStat(0)

# V2 VALUES (obtained from weighted average across pt bins, using v2 vals from: https://alice-notes.web.cern.ch/system/files/notes/analysis/1228/2022-09-15-AN_PIDflowInSmallSystems_v4.pdf 
# and yield vals from: https://www.hepdata.net/record/ins1244523?version=2)

TRIGGER_V2 = 0.092
LAMBDA_V2 = 0.111
ASSOCIATED_V2 = 0.112

LAMBDA_V2_LOWPT = 0.075
ASSOCIATED_V2_LOWPT = 0.100

LAMBDA_V2_HIGHPT = 0.137
ASSOCIATED_V2_HIGHPT = 0.119

def find_v2_baseline(h, name, outfile):

    c = rt.TCanvas("c", "c", 800, 600)
    h.GetXaxis().SetRangeUser(-1.2, 1.2 - 0.001)

    h_dphi = h.ProjectionY().Clone("h_dphi")

    left_v2_fit = rt.TF1("left_v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, -1)
    right_v2_fit = rt.TF1("right_v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", 1, rt.TMath.Pi()/2)
    total_v2_fit = rt.TF1("v2_fit", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

    if "h_h_" in name:
        if "_0_20" in name:
            left_v2_fit.FixParameter(1, TRIGGER_V2)
            right_v2_fit.FixParameter(1, TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, ASSOCIATED_V2)
            
        elif "_20_50" in name:
            left_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.85*ASSOCIATED_V2)

        elif "_50_80" in name:
            left_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.5*ASSOCIATED_V2)

        elif "_0_80" in name:
            left_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2)
                right_v2_fit.FixParameter(2, 0.95*ASSOCIATED_V2)

    elif "h_lambda_" in name:

        if "_0_20" in name:
            left_v2_fit.FixParameter(1, TRIGGER_V2)
            right_v2_fit.FixParameter(1, TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, LAMBDA_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, LAMBDA_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, LAMBDA_V2)
                right_v2_fit.FixParameter(2, LAMBDA_V2)


        elif "_20_50" in name:
            left_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.85*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.85*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.85*LAMBDA_V2)

        elif "_50_80" in name:
            left_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.5*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.5*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.5*LAMBDA_V2)

        elif "_0_80" in name:
            left_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)
            right_v2_fit.FixParameter(1, 0.95*TRIGGER_V2)

            if "lowpt" in name:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_LOWPT)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_LOWPT)
            elif "highpt" in name:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_HIGHPT)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2_HIGHPT)
            elif "central" in name:
                left_v2_fit.FixParameter(2, 0.95*LAMBDA_V2)
                right_v2_fit.FixParameter(2, 0.95*LAMBDA_V2)


    left_v2_fit.SetParameter(0, 0.5*h_dphi.GetMaximum())
    right_v2_fit.SetParameter(0, 0.5*h_dphi.GetMaximum())

    h_dphi.Fit(left_v2_fit, "R")
    h_dphi.Fit(right_v2_fit, "R+")


    total_v2_fit.SetParameter(0, (left_v2_fit.GetParameter(0)+right_v2_fit.GetParameter(0))/2)
    total_v2_fit.SetParameter(1, (left_v2_fit.GetParameter(1)+right_v2_fit.GetParameter(1))/2)
    total_v2_fit.SetParameter(2, (left_v2_fit.GetParameter(2)+right_v2_fit.GetParameter(2))/2)

    h_dphi.SetTitle("")
    h_dphi.GetXaxis().SetTitle("#Delta#varphi")
    h_dphi.GetYaxis().SetTitle("Counts")
    h_dphi.SetLineColor(rt.kBlack)
    h_dphi.SetLineWidth(2)
    h_dphi.SetMarkerColor(rt.kBlack)
    h_dphi.SetMarkerStyle(34)
    h_dphi.SetMarkerSize(2)

    total_v2_fit.SetLineColor(rt.kBlue -2)
    total_v2_fit.SetLineWidth(2)
    leg = rt.TLegend(0.16, 0.65, 0.44, 0.87)
    leg.AddEntry(h_dphi, name + " data (|#Delta#eta| < 1.2)", "lep")
    leg.AddEntry(total_v2_fit, "v_{2} fit", "l")
    leg.SetBorderSize(0)
    h_dphi.GetYaxis().SetRangeUser(0.8*h_dphi.GetMinimum(), 1.2*h_dphi.GetMaximum())
    h_dphi.Draw()
    total_v2_fit.Draw("SAME")
    leg.Draw("SAME")
    c.SaveAs("figures/v2fit_" + name + ".pdf")
    baseline = total_v2_fit.GetParameter(0)
    baseline_error = total_v2_fit.GetParError(0)
    trigger_v2 = total_v2_fit.GetParameter(1)
    associated_v2 = total_v2_fit.GetParameter(2)
    outfile.write(name + " " + str(baseline) + " " + str(baseline_error) + " " + str(trigger_v2) + " " + str(associated_v2) + "\n")



infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
infile_lowpt = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
infile_highpt = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

file_dict = {
    "central": infile,
    "lowpt": infile_lowpt,
    "highpt": infile_highpt,
}

h_lambda_base_hist_name = "h_lambda_2d_subtracted"
h_h_base_hist_name = "h_h_2d_mixcor"
mult_strings = ["_0_20", "_20_50", "_50_80", "_0_80"]

with open("v2_fitpars.txt", "w") as outfile:
    for key, f in file_dict.items():
        for mult_string in mult_strings:
            h_lambda = f.Get(h_lambda_base_hist_name + mult_string)
            h_h = f.Get(h_h_base_hist_name + mult_string)
            find_v2_baseline(h_h, "h_h_" + key + mult_string, outfile)
            find_v2_baseline(h_lambda,"h_lambda_" + key + mult_string, outfile)

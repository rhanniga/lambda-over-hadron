import ROOT as rt
rt.gStyle.SetOptStat(0)
import array as arr

import math

import strangehelper as sh

from systematics_helper import AnalysisBin, AnalysisBins
from yield_extractor import YieldExtractor, FitType

EPSILON = 0.0001

# define the bins that will *likely* never change

CENTRALITY_BINS = AnalysisBins("centrality_bins",
                                [AnalysisBin("cent_50_80", 1.0, 2 - EPSILON), 
                                AnalysisBin("cent_20_50", 2.0, 3.0 - EPSILON),
                                AnalysisBin("cent_0_20", 3.0, 4.0 - EPSILON)
                                # [AnalysisBin("cent_0_100", 0, 0)
                                ])

TRIGGER_PT_BINS = AnalysisBins("trigger_pt_bins", 
                               [AnalysisBin("trigger_4_8", 4.0, 8.0 - EPSILON)])

ASSOCIATED_PT_BINS = AnalysisBins("associated_pt_bins",
                                [AnalysisBin("assoc_2_4", 2.0, 4.0 - EPSILON),
                                AnalysisBin("assoc_15_25", 1.5, 2.5 - EPSILON),
                                AnalysisBin("assoc_25_4", 2.5, 4.0 - EPSILON)])

DELTA_ETA_BINS = AnalysisBins("delta_eta_bins",
                                [AnalysisBin("delta_eta_12", -1.2, 1.2 - EPSILON)])

ETA_BINS = AnalysisBins("eta_bins",
                        [AnalysisBin("eta_08", -0.8, 0.8 - EPSILON)])
    
MODEL_BINS = AnalysisBins("model_bins",
                        [AnalysisBin("dpmjet", 1, 1, "output/AnalysisResults.root")])
                        # AnalysisBin("dpmjet_noeta", 1, 1, "output/dpmjet_17f2b_fast_noeta.root"),
                        # [AnalysisBin("epos", 2, 2, "output/epos_17f2a_fast.root")])
                        # [AnalysisBin("phsd", 3, 3, "output/phsd_out_200mil.root")])


def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
    return deriv * kappa_error

def get_von_fit(dphi_dist, avg4=False):
    dphi_dist = dphi_dist.Clone()
    if avg4:
        avg = (dphi_dist.GetBinContent(1) 
                + dphi_dist.GetBinContent(8)
                + dphi_dist.GetBinContent(9)
                + dphi_dist.GetBinContent(16))/4
    else:
        avg = (dphi_dist.GetBinContent(1) 
                + dphi_dist.GetBinContent(2)
                + dphi_dist.GetBinContent(7)
                + dphi_dist.GetBinContent(8)
                + dphi_dist.GetBinContent(9)
                + dphi_dist.GetBinContent(16))/6
    von_fit_string = "[0]"
    von_fit_string += " + [1]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x- TMath::Pi()))"
    von_fit = rt.TF1("von_fit", von_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
    von_fit.FixParameter(0, avg)
    von_fit.SetParLimits(1, 0.0001, 2)
    von_fit.SetParameter(1, 0.01)
    von_fit.SetParLimits(2, 0.1, 20)
    von_fit.SetParameter(2, 7)
    von_fit.SetParLimits(3, 0.0001, 2)
    von_fit.SetParameter(3, 0.01)
    von_fit.SetParLimits(4, 0.1, 20)
    von_fit.SetParameter(4, 2)
    dphi_dist.Fit(von_fit, "R")
    return von_fit.Clone()

def get_gaus_fit(dphi_dist):

    dphi_dist = dphi_dist.Clone()
    avg = (dphi_dist.GetBinContent(1) 
            + dphi_dist.GetBinContent(2)
            + dphi_dist.GetBinContent(7)
            + dphi_dist.GetBinContent(8)
            + dphi_dist.GetBinContent(9)
            + dphi_dist.GetBinContent(16))/6
    gaus_fit_string = "[0]"
    gaus_fit_string += " + gaus(1) + gaus(4) + gaus(7) + gaus(10)"
    gaus_fit = rt.TF1("gaus_fit", gaus_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
    gaus_fit.FixParameter(0, avg)

    gaus_fit.SetParLimits(1, 0.0001, 2)
    gaus_fit.SetParameter(1, 0.01)
    gaus_fit.FixParameter(2, 0)
    gaus_fit.SetParLimits(3, 0.05, 2)
    gaus_fit.SetParameter(3, 0.3)

    gaus_fit.SetParLimits(4, 0, 1)
    gaus_fit.SetParameter(4, 0)
    gaus_fit.FixParameter(5, 2*rt.TMath.Pi())
    gaus_fit.SetParLimits(6, 0.05, 2)
    gaus_fit.SetParameter(6, 0.6)

    gaus_fit.SetParLimits(7, 0.0001, 2)
    gaus_fit.SetParameter(7, 0.01)
    gaus_fit.FixParameter(8, rt.TMath.Pi())
    gaus_fit.SetParLimits(9, 0.05, 2)
    gaus_fit.SetParameter(9, 0.3)

    gaus_fit.SetParLimits(10, 0.0001, 2)
    gaus_fit.SetParameter(10, 0.01)
    gaus_fit.FixParameter(11, -rt.TMath.Pi())
    gaus_fit.SetParLimits(12, 0.05, 2)
    gaus_fit.SetParameter(12, 0.1)

    dphi_dist.Fit(gaus_fit, "R")
    return gaus_fit.Clone()


def get_dphi_dists(model, 
                   centrality_bin,
                   trigger_pt_bin,
                   associated_pt_bin,
                   delta_eta_bin,
                   eta_bin,
                   output_file=None):

    # get unique string for cloning to avoid name conflicts
    unique_clone_string = model.name + "_" + \
                        centrality_bin.name + "_" + \
                        trigger_pt_bin.name + "_" + \
                        associated_pt_bin.name + "_" + \
                        delta_eta_bin.name 

    # projecting trigger distribution within the specified centrality and trigger pt bins down to 1d (pt)

    if model.name == "phsd":
        trigger_dist = model.input_file.Get("fTriggerDist_MC") 

        h_h_dist = model.input_file.Get("fDphiHH_MC")
        h_lambda_dist = model.input_file.Get("fDphiHLambda_MC")
        h_phi_dist = model.input_file.Get("fDphiHPhi_MC")

        h_h_dist_no_eta_cut = model.input_file.Get("fDphiHH_MC").Clone("h_h_dist_no_eta_cut")
        h_lambda_dist_no_eta_cut = model.input_file.Get("fDphiHLambda_MC").Clone("h_lambda_dist_no_eta_cut")
        h_phi_dist_no_eta_cut = model.input_file.Get("fDphiHPhi_MC").Clone("h_phi_dist_no_eta_cut")
    else:
        trigger_dist = model.input_file.Get("h-lambda").FindObject("fTriggerDist_MC")
        h_h_dist = model.input_file.Get("h-lambda").FindObject("fDphiHH_MC")
        h_h_mixed_dist = model.input_file.Get("h-lambda").FindObject("fDphiHHMixed_MC")
        h_lambda_dist = model.input_file.Get("h-lambda").FindObject("fDphiHLambda_MC")
        h_lambda_mixed_dist = model.input_file.Get("h-lambda").FindObject("fDphiHLambdaMixed_MC")
        h_phi_dist = model.input_file.Get("h-lambda").FindObject("fDphiHPhi_MC")
        h_phi_mixed_dist = model.input_file.Get("h-lambda").FindObject("fDphiHPhiMixed_MC")

        # h_h_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHH_MC_no_eta_cut")
        # h_lambda_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHLambda_MC_no_eta_cut")
        # h_phi_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHPhi_MC_no_eta_cut")

        h_h_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHH_MC").Clone("h_h_dist_no_eta_cut")
        h_lambda_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHLambda_MC").Clone("h_lambda_dist_no_eta_cut")
        h_phi_dist_no_eta_cut = model.input_file.Get("h-lambda").FindObject("fDphiHPhi_MC").Clone("h_phi_dist_no_eta_cut")

    trigger_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    trigger_dist.GetAxis(2).SetRangeUser(eta_bin.lower_bound, eta_bin.upper_bound)
    trigger_dist.GetAxis(4).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    trigger_pt_dist = trigger_dist.Projection(0).Clone("trigger_pt_dist_" + unique_clone_string)

    num_triggers = trigger_pt_dist.Integral()

    h_h_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_h_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_h_dist.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_h_dist = h_h_dist.Projection(2, 3, 4).Clone("h_h_dist_" + unique_clone_string)

    h_h_dist_no_eta_cut.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_h_dist_no_eta_cut.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_h_dist_no_eta_cut.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_h_dist_no_eta_cut = h_h_dist_no_eta_cut.Projection(2, 3, 4).Clone("h_h_dist_no_eta_cut_" + unique_clone_string)

    if model.name != "phsd":
        h_h_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
        h_h_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
        h_h_mixed_dist = h_h_mixed_dist.Projection(2, 3, 4).Clone("h_h_mixed_dist_" + unique_clone_string)

    if model.name == "phsd":
        h_h_2d_mixcor = h_h_dist.Project3D("xye").Clone("h_h_2d_mixcor_" + unique_clone_string)
        h_h_2d_mixcor_no_eta_cut = h_h_dist_no_eta_cut.Project3D("xye").Clone("h_h_2d_mixcor_no_eta_cut_" + unique_clone_string)
    else:  
        h_h_2d_mixcor_no_eta_cut = h_h_dist_no_eta_cut.Project3D("xye").Clone("h_h_2d_mixcor_no_eta_cut_" + unique_clone_string)
        h_h_2d_mixcor = sh.make_h_h_mixed_corrections(h_h_dist, h_h_mixed_dist, unique_clone_string, output_file)
    

    # do per-trigger scaling
    h_h_2d_mixcor.Scale(1.0/num_triggers)
    h_h_2d_mixcor_no_eta_cut.Scale(1.0/num_triggers)

    # project down to dPhi in given dEta range
    h_h_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
    h_h_dphi = h_h_2d_mixcor.ProjectionY("h_h_dphi_" + unique_clone_string)

    h_lambda_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_lambda_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_lambda_dist.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_lambda_dist = h_lambda_dist.Projection(2, 3, 4).Clone("h_lambda_dist_" + unique_clone_string)

    h_lambda_dist_no_eta_cut.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_lambda_dist_no_eta_cut.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_lambda_dist_no_eta_cut.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_lambda_dist_no_eta_cut = h_lambda_dist_no_eta_cut.Projection(2, 3, 4).Clone("h_lambda_dist_no_eta_cut_" + unique_clone_string)

    if model.name != "phsd":
        h_lambda_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
        h_lambda_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
        h_lambda_mixed_dist = h_lambda_mixed_dist.Projection(2, 3, 4).Clone("h_lambda_mixed_dist_" + unique_clone_string) 

    # do h-lambda mixed event correction
    if model.name == "phsd":
        h_lambda_2d_mixcor = h_lambda_dist.Project3D("xye").Clone("h_lambda_2d_mixcor_" + unique_clone_string)
        h_lambda_2d_mixcor_no_eta_cut = h_lambda_dist_no_eta_cut.Project3D("xye").Clone("h_lambda_2d_mixcor_no_eta_cut_" + unique_clone_string)
    else:
        h_lambda_2d_mixcor_no_eta_cut = h_lambda_dist_no_eta_cut.Project3D("xye").Clone("h_lambda_2d_mixcor_no_eta_cut_" + unique_clone_string)
        h_lambda_2d_mixcor = sh.make_h_h_mixed_corrections(h_lambda_dist, h_lambda_mixed_dist, unique_clone_string, output_file)

    # do per-trigger scaling 
    h_lambda_2d_mixcor_no_eta_cut.Scale(1.0/num_triggers)
    h_lambda_2d_mixcor.Scale(1.0/num_triggers)

    # correct for branching ratio
    lambda_br = 0.639
    h_lambda_2d_mixcor_no_eta_cut.Scale(1/lambda_br)
    h_lambda_2d_mixcor.Scale(1/lambda_br)



    test_c = rt.TCanvas("test_c", "", 800, 600)
    test_c.SetLeftMargin(0.15)
    test_c.SetRightMargin(0.05)
    test_c.SetBottomMargin(0.12)
    test_c.SetTopMargin(0.05)

    # test_c.SetPhi(50)
    # test_c.SetTheta(20)
    # h_lambda_2d_mixcor.SetTitle("")
    # h_lambda_2d_mixcor.GetXaxis().SetTitle("#Delta#eta")
    # h_lambda_2d_mixcor.GetYaxis().SetTitle("#Delta#varphi")
    # h_lambda_2d_mixcor.GetXaxis().SetTitleSize(0.05)
    # h_lambda_2d_mixcor.GetYaxis().SetTitleSize(0.05)
    # h_lambda_2d_mixcor.GetXaxis().SetTitleOffset(1.2)
    # h_lambda_2d_mixcor.GetYaxis().SetTitleOffset(1.2)
    # h_lambda_2d_mixcor.GetXaxis().SetRangeUser(-1.2, 1.2)
    # h_lambda_2d_mixcor.Draw("SURF1")






    h_lambda_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
    h_lambda_dphi = h_lambda_2d_mixcor.ProjectionY("h_lambda_dphi_" + unique_clone_string)

    h_phi_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_phi_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_phi_dist.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_phi_dist = h_phi_dist.Projection(2, 3, 4).Clone("h_phi_dist_" + unique_clone_string)

    h_phi_dist_no_eta_cut.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_phi_dist_no_eta_cut.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_phi_dist_no_eta_cut.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_phi_dist_no_eta_cut = h_phi_dist_no_eta_cut.Projection(2, 3, 4).Clone("h_phi_dist_no_eta_cut_" + unique_clone_string)

    if model.name != "phsd":
        h_phi_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
        h_phi_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
        h_phi_mixed_dist = h_phi_mixed_dist.Projection(2, 3, 4).Clone("h_phi_mixed_dist_" + unique_clone_string)

    # do h-phi mixed event correction
    if model.name == "phsd":
        h_phi_2d_mixcor = h_phi_dist.Project3D("xye").Clone("h_phi_2d_mixcor_" + unique_clone_string)
        h_phi_2d_mixcor_no_eta_cut = h_phi_dist_no_eta_cut.Project3D("xye").Clone("h_phi_2d_mixcor_no_eta_cut_" + unique_clone_string)
    else:
        h_phi_2d_mixcor_no_eta_cut = h_phi_dist_no_eta_cut.Project3D("xye").Clone("h_phi_2d_mixcor_no_eta_cut_" + unique_clone_string)
        h_phi_2d_mixcor = sh.make_h_h_mixed_corrections(h_phi_dist, h_phi_mixed_dist, unique_clone_string, output_file)

    # do per-trigger scaling
    h_phi_2d_mixcor.Scale(1.0/num_triggers)
    h_phi_2d_mixcor_no_eta_cut.Scale(1.0/num_triggers)

    # correct for branching ratio
    phi_br = 0.492
    h_phi_2d_mixcor.Scale(1/phi_br)
    h_phi_2d_mixcor_no_eta_cut.Scale(1/phi_br)


    # h_phi_2d_mixcor.SetTitle("")
    # h_phi_2d_mixcor.GetXaxis().SetTitle("#Delta#eta")
    # h_phi_2d_mixcor.GetYaxis().SetTitle("#Delta#varphi")
    # h_phi_2d_mixcor.GetXaxis().SetTitleSize(0.05)
    # h_phi_2d_mixcor.GetYaxis().SetTitleSize(0.05)
    # h_phi_2d_mixcor.GetXaxis().SetTitleOffset(1.2)
    # h_phi_2d_mixcor.GetYaxis().SetTitleOffset(1.2)
    # h_phi_2d_mixcor.GetXaxis().SetRangeUser(-1.2, 1.2)
    # h_phi_2d_mixcor.Draw("SURF1")
    # test_c.SaveAs("h_phi_2d_phsd.pdf")


    h_phi_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
    h_phi_dphi = h_phi_2d_mixcor.ProjectionY("h_phi_dphi_" + unique_clone_string)

    h_phi_dphi.SetTitle("")
    h_phi_dphi.GetXaxis().SetTitle("#Delta#varphi")
    h_phi_dphi.GetYaxis().SetTitle("1/N_{trig} dN/d#Delta#varphi")
    h_phi_dphi.SetMarkerStyle(20)
    h_phi_dphi.SetMarkerSize(1.0)
    h_phi_dphi.SetLineWidth(2)
    h_phi_dphi.SetLineColor(rt.kBlue+2)
    h_phi_dphi.SetMarkerColor(rt.kBlue+2)

    # h_phi_dphi.Draw()
    # test_c.SaveAs("h_phi_1d_epos.pdf")
    # quit()


    return h_lambda_dphi, h_phi_dphi, h_h_dphi

c = rt.TCanvas("c", "", 800, 600)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.12)


for model in MODEL_BINS.analysis_bins:
    for trigger_pt_bin in TRIGGER_PT_BINS.analysis_bins:
        for assoc_pt_bin in ASSOCIATED_PT_BINS.analysis_bins:
            if model.name == "dpmjet":
                outfile = rt.TFile(f"dpmjet_graphs_{assoc_pt_bin.name}.root", "recreate")
            elif model.name == "phsd":
                outfile = rt.TFile(f"phsd_graphs_{assoc_pt_bin.name}.root", "recreate")
            for delta_eta_bin in DELTA_ETA_BINS.analysis_bins:

                h_lambda_near_width_cent = []
                h_lambda_near_width_cent_err = []
                h_lambda_away_width_cent = []
                h_lambda_away_width_cent_err = []

                h_lambda_near_yield_cent = []
                h_lambda_near_yield_cent_err = []
                h_lambda_away_yield_cent = []
                h_lambda_away_yield_cent_err = []

                h_phi_near_width_cent = []
                h_phi_near_width_cent_err = []
                h_phi_away_width_cent = []
                h_phi_away_width_cent_err = []

                h_h_near_width_cent = []
                h_h_near_width_cent_err = []
                h_h_away_width_cent = []
                h_h_away_width_cent_err = []

                h_h_near_yield_cent = []
                h_h_near_yield_cent_err = []
                h_h_away_yield_cent = []
                h_h_away_yield_cent_err = []

                lambda_hadron_near_ratio_cent = []
                lambda_hadron_near_ratio_cent_err = []
                lambda_hadron_away_ratio_cent = []
                lambda_hadron_away_ratio_cent_err = []
                lambda_hadron_ue_ratio_cent = []
                lambda_hadron_ue_ratio_cent_err = []
                lambda_hadron_total_ratio_cent = []
                lambda_hadron_total_ratio_cent_err = []

                phi_hadron_near_ratio_cent = []
                phi_hadron_near_ratio_cent_err = []
                phi_hadron_away_ratio_cent = []
                phi_hadron_away_ratio_cent_err = []
                phi_hadron_ue_ratio_cent = []
                phi_hadron_ue_ratio_cent_err = []
                phi_hadron_total_ratio_cent = []
                phi_hadron_total_ratio_cent_err = []

                lambda_phi_near_ratio_cent = []
                lambda_phi_near_ratio_cent_err = []
                lambda_phi_away_ratio_cent = []
                lambda_phi_away_ratio_cent_err = []
                lambda_phi_ue_ratio_cent = []
                lambda_phi_ue_ratio_cent_err = []
                lambda_phi_total_ratio_cent = []
                lambda_phi_total_ratio_cent_err = []

                for centrality_bin in CENTRALITY_BINS.analysis_bins:
                    h_lambda_dphi, h_phi_dphi, h_h_dphi = get_dphi_dists(model, centrality_bin, trigger_pt_bin, assoc_pt_bin, delta_eta_bin, ETA_BINS.analysis_bins[0], None)

                    h_lambda_yield_extractor = YieldExtractor(h_lambda_dphi)
                    h_phi_yield_extractor = YieldExtractor(h_phi_dphi)
                    h_h_yield_extractor = YieldExtractor(h_h_dphi)

                    h_lambda_von_fit = get_von_fit(h_lambda_dphi)
                    h_phi_von_fit = get_von_fit(h_phi_dphi)
                    h_h_von_fit = get_von_fit(h_h_dphi)

                    h_lambda_von_fit_avg4 = get_von_fit(h_lambda_dphi, True)
                    h_phi_von_fit_avg4 = get_von_fit(h_phi_dphi, True)
                    h_h_von_fit_avg4 = get_von_fit(h_h_dphi, True)

                    h_lambda_gaus_fit = get_gaus_fit(h_lambda_dphi)
                    h_phi_gaus_fit = get_gaus_fit(h_phi_dphi)
                    h_h_gaus_fit = get_gaus_fit(h_h_dphi)

                    h_lambda_near_kappa = h_lambda_von_fit.GetParameter(2)
                    h_lambda_near_kappa_err = h_lambda_von_fit.GetParError(2)
                    h_lambda_near_width = get_width_from_kappa(h_lambda_near_kappa)
                    h_lambda_near_width_err = get_width_error_from_kappa(h_lambda_near_kappa, h_lambda_near_kappa_err)
                    h_lambda_near_kappa_avg4 = h_lambda_von_fit_avg4.GetParameter(2)
                    h_lambda_near_kappa_avg4_err = h_lambda_von_fit_avg4.GetParError(2)
                    h_lambda_near_width_avg4 = get_width_from_kappa(h_lambda_near_kappa_avg4)
                    h_lambda_near_width_avg4_err = get_width_error_from_kappa(h_lambda_near_kappa_avg4, h_lambda_near_kappa_avg4_err)
                    h_lambda_near_width_gaus = h_lambda_gaus_fit.GetParameter(3)
                    h_lambda_near_width_gaus_err = h_lambda_gaus_fit.GetParError(3)
                    h_lambda_near_width_syserr = ((h_lambda_near_width_avg4/h_lambda_near_width) - 1)**2 + ((h_lambda_near_width_gaus/h_lambda_near_width) - 1)**2
                    h_lambda_near_width_syserr = h_lambda_near_width * math.sqrt(h_lambda_near_width_syserr/2)
                    h_lambda_near_width_err = math.sqrt(h_lambda_near_width_err**2 + h_lambda_near_width_syserr**2)
                    h_lambda_near_width_cent.append(h_lambda_near_width)
                    h_lambda_near_width_cent_err.append(h_lambda_near_width_err)
                    h_lambda_away_kappa = h_lambda_von_fit.GetParameter(4)
                    h_lambda_away_kappa_err = h_lambda_von_fit.GetParError(4)
                    h_lambda_away_width = get_width_from_kappa(h_lambda_away_kappa)
                    h_lambda_away_width_err = get_width_error_from_kappa(h_lambda_away_kappa, h_lambda_away_kappa_err)
                    h_lambda_away_kappa_avg4 = h_lambda_von_fit_avg4.GetParameter(4)
                    h_lambda_away_kappa_avg4_err = h_lambda_von_fit_avg4.GetParError(4)
                    h_lambda_away_width_avg4 = get_width_from_kappa(h_lambda_away_kappa_avg4)
                    h_lambda_away_width_avg4_err = get_width_error_from_kappa(h_lambda_away_kappa_avg4, h_lambda_away_kappa_avg4_err)
                    h_lambda_away_width_gaus = h_lambda_gaus_fit.GetParameter(9)
                    h_lambda_away_width_gaus_err = h_lambda_gaus_fit.GetParError(9)
                    h_lambda_away_width_syserr = ((h_lambda_away_width_avg4/h_lambda_away_width) - 1)**2 + ((h_lambda_away_width_gaus/h_lambda_away_width) - 1)**2
                    h_lambda_away_width_syserr = h_lambda_away_width * math.sqrt(h_lambda_away_width_syserr/2)
                    h_lambda_away_width_err = math.sqrt(h_lambda_away_width_err**2 + h_lambda_away_width_syserr**2)
                    h_lambda_away_width_cent.append(h_lambda_away_width)
                    h_lambda_away_width_cent_err.append(h_lambda_away_width_err)


                    

                    h_phi_near_kappa = h_phi_von_fit.GetParameter(2)
                    h_phi_near_kappa_err = h_phi_von_fit.GetParError(2)
                    h_phi_near_width = get_width_from_kappa(h_phi_near_kappa)
                    h_phi_near_width_err = get_width_error_from_kappa(h_phi_near_kappa, h_phi_near_kappa_err)
                    h_phi_near_kappa_avg4 = h_phi_von_fit_avg4.GetParameter(2)
                    h_phi_near_kappa_avg4_err = h_phi_von_fit_avg4.GetParError(2)
                    h_phi_near_width_avg4 = get_width_from_kappa(h_phi_near_kappa_avg4)
                    h_phi_near_width_avg4_err = get_width_error_from_kappa(h_phi_near_kappa_avg4, h_phi_near_kappa_avg4_err)
                    h_phi_near_width_gaus = h_phi_gaus_fit.GetParameter(3)
                    h_phi_near_width_gaus_err = h_phi_gaus_fit.GetParError(3)
                    h_phi_near_width_syserr = ((h_phi_near_width_avg4/h_phi_near_width) - 1)**2 + ((h_phi_near_width_gaus/h_phi_near_width) - 1)**2
                    h_phi_near_width_syserr = h_phi_near_width * math.sqrt(h_phi_near_width_syserr/2)
                    h_phi_near_width_err = math.sqrt(h_phi_near_width_err**2 + h_phi_near_width_syserr**2)
                    h_phi_near_width_cent.append(h_phi_near_width)
                    h_phi_near_width_cent_err.append(h_phi_near_width_err)
                    h_phi_away_kappa = h_phi_von_fit.GetParameter(4)
                    h_phi_away_kappa_err = h_phi_von_fit.GetParError(4)
                    h_phi_away_width = get_width_from_kappa(h_phi_away_kappa)
                    h_phi_away_width_err = get_width_error_from_kappa(h_phi_away_kappa, h_phi_away_kappa_err)
                    h_phi_away_kappa_avg4 = h_phi_von_fit_avg4.GetParameter(4)
                    h_phi_away_kappa_avg4_err = h_phi_von_fit_avg4.GetParError(4)
                    h_phi_away_width_avg4 = get_width_from_kappa(h_phi_away_kappa_avg4)
                    h_phi_away_width_avg4_err = get_width_error_from_kappa(h_phi_away_kappa_avg4, h_phi_away_kappa_avg4_err)
                    h_phi_away_width_gaus = h_phi_gaus_fit.GetParameter(9)
                    h_phi_away_width_gaus_err = h_phi_gaus_fit.GetParError(9)
                    h_phi_away_width_syserr = ((h_phi_away_width_avg4/h_phi_away_width) - 1)**2 + ((h_phi_away_width_gaus/h_phi_away_width) - 1)**2
                    h_phi_away_width_syserr = h_phi_away_width * math.sqrt(h_phi_away_width_syserr/2)
                    h_phi_away_width_err = math.sqrt(h_phi_away_width_err**2 + h_phi_away_width_syserr**2)
                    h_phi_away_width_cent.append(h_phi_away_width)
                    h_phi_away_width_cent_err.append(h_phi_away_width_err)

                    h_h_near_kappa = h_h_von_fit.GetParameter(2)
                    h_h_near_kappa_err = h_h_von_fit.GetParError(2)
                    h_h_near_width = get_width_from_kappa(h_h_near_kappa)
                    h_h_near_width_err = get_width_error_from_kappa(h_h_near_kappa, h_h_near_kappa_err)
                    h_h_near_kappa_avg4 = h_h_von_fit_avg4.GetParameter(2)
                    h_h_near_kappa_avg4_err = h_h_von_fit_avg4.GetParError(2)
                    h_h_near_width_avg4 = get_width_from_kappa(h_h_near_kappa_avg4)
                    h_h_near_width_avg4_err = get_width_error_from_kappa(h_h_near_kappa_avg4, h_h_near_kappa_avg4_err)
                    h_h_near_width_gaus = h_h_gaus_fit.GetParameter(3)
                    h_h_near_width_gaus_err = h_h_gaus_fit.GetParError(3)
                    h_h_near_width_syserr = ((h_h_near_width_avg4/h_h_near_width) - 1)**2 + ((h_h_near_width_gaus/h_h_near_width) - 1)**2
                    h_h_near_width_syserr = h_h_near_width * math.sqrt(h_h_near_width_syserr/2)
                    h_h_near_width_err = math.sqrt(h_h_near_width_err**2 + h_h_near_width_syserr**2)
                    h_h_near_width_cent.append(h_h_near_width)
                    h_h_near_width_cent_err.append(h_h_near_width_err)
                    h_h_away_kappa = h_h_von_fit.GetParameter(4)
                    h_h_away_kappa_err = h_h_von_fit.GetParError(4)
                    h_h_away_width = get_width_from_kappa(h_h_away_kappa)
                    h_h_away_width_err = get_width_error_from_kappa(h_h_away_kappa, h_h_away_kappa_err)
                    h_h_away_kappa_avg4 = h_h_von_fit_avg4.GetParameter(4)
                    h_h_away_kappa_avg4_err = h_h_von_fit_avg4.GetParError(4)
                    h_h_away_width_avg4 = get_width_from_kappa(h_h_away_kappa_avg4)
                    h_h_away_width_avg4_err = get_width_error_from_kappa(h_h_away_kappa_avg4, h_h_away_kappa_avg4_err)
                    h_h_away_width_gaus = h_h_gaus_fit.GetParameter(9)
                    h_h_away_width_gaus_err = h_h_gaus_fit.GetParError(9)
                    h_h_away_width_syserr = ((h_h_away_width_avg4/h_h_away_width) - 1)**2 + ((h_h_away_width_gaus/h_h_away_width) - 1)**2
                    h_h_away_width_syserr = h_h_away_width * math.sqrt(h_h_away_width_syserr/2)
                    h_h_away_width_err = math.sqrt(h_h_away_width_err**2 + h_h_away_width_syserr**2)
                    h_h_away_width_cent.append(h_h_away_width)
                    h_h_away_width_cent_err.append(h_h_away_width_err)




                    h_lambda_yield_extractor.extract_all_yields()
                    h_phi_yield_extractor.extract_all_yields()
                    h_h_yield_extractor.extract_all_yields()

                    lambda_near_yield = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["ns"][0]
                    lambda_near_yield_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["ns"][1]
                    lambda_away_yield = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["as"][0]
                    lambda_away_yield_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["as"][1]
                    lambda_ue_yield = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["ue"][0]
                    lambda_ue_yield_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["ue"][1]
                    lambda_total_yield = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["total"][0]
                    lambda_total_yield_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX]["total"][1]
                    lambda_near_yield_avg4 = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["ns"][0]
                    lambda_near_yield_avg4_err = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["ns"][1]
                    lambda_away_yield_avg4 = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["as"][0]
                    lambda_away_yield_avg4_err = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["as"][1]
                    lambda_ue_yield_avg4 = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["ue"][0]
                    lambda_ue_yield_avg4_err = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["ue"][1]
                    lambda_total_yield_avg4 = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["total"][0]
                    lambda_total_yield_avg4_err = h_lambda_yield_extractor.yields[FitType.AVG_FOUR]["total"][1]
                    lambda_near_yield_avg6nonneg = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][0]
                    lambda_near_yield_avg6nonneg_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][1]
                    lambda_away_yield_avg6nonneg = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][0]
                    lambda_away_yield_avg6nonneg_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][1]
                    lambda_ue_yield_avg6nonneg = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][0]
                    lambda_ue_yield_avg6nonneg_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][1]
                    lambda_total_yield_avg6nonneg = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][0]
                    lambda_total_yield_avg6nonneg_err = h_lambda_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][1]
                    lambda_near_yield_syserror = ((lambda_near_yield_avg4/lambda_near_yield) - 1)**2 + ((lambda_near_yield_avg6nonneg/lambda_near_yield) - 1)**2
                    lambda_near_yield_syserror = lambda_near_yield * math.sqrt(lambda_near_yield_syserror/2)
                    lambda_away_yield_syserror = ((lambda_away_yield_avg4/lambda_away_yield) - 1)**2 + ((lambda_away_yield_avg6nonneg/lambda_away_yield) - 1)**2
                    lambda_away_yield_syserror = lambda_away_yield * math.sqrt(lambda_away_yield_syserror/2)
                    lambda_ue_yield_syserror = ((lambda_ue_yield_avg4/lambda_ue_yield) - 1)**2 + ((lambda_ue_yield_avg6nonneg/lambda_ue_yield) - 1)**2
                    lambda_ue_yield_syserror = lambda_ue_yield * math.sqrt(lambda_ue_yield_syserror/2)
                    lambda_total_yield_syserror = ((lambda_total_yield_avg4/lambda_total_yield) - 1)**2 + ((lambda_total_yield_avg6nonneg/lambda_total_yield) - 1)**2
                    lambda_total_yield_syserror = lambda_total_yield * math.sqrt(lambda_total_yield_syserror/2)

                    lambda_near_yield_err = math.sqrt(lambda_near_yield_err**2 + lambda_near_yield_syserror**2)
                    lambda_away_yield_err = math.sqrt(lambda_away_yield_err**2 + lambda_away_yield_syserror**2)
                    lambda_ue_yield_err = math.sqrt(lambda_ue_yield_err**2 + lambda_ue_yield_syserror**2)
                    lambda_total_yield_err = math.sqrt(lambda_total_yield_err**2 + lambda_total_yield_syserror**2)

                    phi_near_yield = h_phi_yield_extractor.yields[FitType.AVG_SIX]["ns"][0]
                    phi_near_yield_err = h_phi_yield_extractor.yields[FitType.AVG_SIX]["ns"][1]
                    phi_away_yield = h_phi_yield_extractor.yields[FitType.AVG_SIX]["as"][0]
                    phi_away_yield_err = h_phi_yield_extractor.yields[FitType.AVG_SIX]["as"][1]
                    phi_ue_yield = h_phi_yield_extractor.yields[FitType.AVG_SIX]["ue"][0]
                    phi_ue_yield_err = h_phi_yield_extractor.yields[FitType.AVG_SIX]["ue"][1]
                    phi_total_yield = h_phi_yield_extractor.yields[FitType.AVG_SIX]["total"][0]
                    phi_total_yield_err = h_phi_yield_extractor.yields[FitType.AVG_SIX]["total"][1]
                    phi_near_yield_avg4 = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["ns"][0]
                    phi_near_yield_avg4_err = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["ns"][1]
                    phi_away_yield_avg4 = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["as"][0]
                    phi_away_yield_avg4_err = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["as"][1]
                    phi_ue_yield_avg4 = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["ue"][0]
                    phi_ue_yield_avg4_err = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["ue"][1]
                    phi_total_yield_avg4 = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["total"][0]
                    phi_total_yield_avg4_err = h_phi_yield_extractor.yields[FitType.AVG_FOUR]["total"][1]
                    phi_near_yield_avg6nonneg = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][0]
                    phi_near_yield_avg6nonneg_err = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][1]
                    phi_away_yield_avg6nonneg = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][0]
                    phi_away_yield_avg6nonneg_err = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][1]
                    phi_ue_yield_avg6nonneg = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][0]
                    phi_ue_yield_avg6nonneg_err = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][1]
                    phi_total_yield_avg6nonneg = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][0]
                    phi_total_yield_avg6nonneg_err = h_phi_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][1]
                    phi_near_yield_syserror = ((phi_near_yield_avg4/phi_near_yield) - 1)**2 + ((phi_near_yield_avg6nonneg/phi_near_yield) - 1)**2
                    phi_near_yield_syserror = phi_near_yield * math.sqrt(phi_near_yield_syserror/2)
                    phi_away_yield_syserror = ((phi_away_yield_avg4/phi_away_yield) - 1)**2 + ((phi_away_yield_avg6nonneg/phi_away_yield) - 1)**2
                    phi_away_yield_syserror = phi_away_yield * math.sqrt(phi_away_yield_syserror/2)
                    phi_ue_yield_syserror = ((phi_ue_yield_avg4/phi_ue_yield) - 1)**2 + ((phi_ue_yield_avg6nonneg/phi_ue_yield) - 1)**2
                    phi_ue_yield_syserror = phi_ue_yield * math.sqrt(phi_ue_yield_syserror/2)
                    phi_total_yield_syserror = ((phi_total_yield_avg4/phi_total_yield) - 1)**2 + ((phi_total_yield_avg6nonneg/phi_total_yield) - 1)**2
                    phi_total_yield_syserror = phi_total_yield * math.sqrt(phi_total_yield_syserror/2)

                    phi_near_yield_err = math.sqrt(phi_near_yield_err**2 + phi_near_yield_syserror**2)
                    phi_away_yield_err = math.sqrt(phi_away_yield_err**2 + phi_away_yield_syserror**2)
                    phi_ue_yield_err = math.sqrt(phi_ue_yield_err**2 + phi_ue_yield_syserror**2)
                    phi_total_yield_err = math.sqrt(phi_total_yield_err**2 + phi_total_yield_syserror**2)

                    h_near_yield = h_h_yield_extractor.yields[FitType.AVG_SIX]["ns"][0]
                    h_near_yield_err = h_h_yield_extractor.yields[FitType.AVG_SIX]["ns"][1]
                    h_away_yield = h_h_yield_extractor.yields[FitType.AVG_SIX]["as"][0]
                    h_away_yield_err = h_h_yield_extractor.yields[FitType.AVG_SIX]["as"][1]
                    h_ue_yield = h_h_yield_extractor.yields[FitType.AVG_SIX]["ue"][0]
                    h_ue_yield_err = h_h_yield_extractor.yields[FitType.AVG_SIX]["ue"][1]
                    h_total_yield = h_h_yield_extractor.yields[FitType.AVG_SIX]["total"][0]
                    h_total_yield_err = h_h_yield_extractor.yields[FitType.AVG_SIX]["total"][1]
                    h_near_yield_avg4 = h_h_yield_extractor.yields[FitType.AVG_FOUR]["ns"][0]
                    h_near_yield_avg4_err = h_h_yield_extractor.yields[FitType.AVG_FOUR]["ns"][1]
                    h_away_yield_avg4 = h_h_yield_extractor.yields[FitType.AVG_FOUR]["as"][0]
                    h_away_yield_avg4_err = h_h_yield_extractor.yields[FitType.AVG_FOUR]["as"][1]
                    h_ue_yield_avg4 = h_h_yield_extractor.yields[FitType.AVG_FOUR]["ue"][0]
                    h_ue_yield_avg4_err = h_h_yield_extractor.yields[FitType.AVG_FOUR]["ue"][1]
                    h_total_yield_avg4 = h_h_yield_extractor.yields[FitType.AVG_FOUR]["total"][0]
                    h_total_yield_avg4_err = h_h_yield_extractor.yields[FitType.AVG_FOUR]["total"][1]
                    h_near_yield_avg6nonneg = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][0]
                    h_near_yield_avg6nonneg_err = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ns"][1]
                    h_away_yield_avg6nonneg = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][0]
                    h_away_yield_avg6nonneg_err = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["as"][1]
                    h_ue_yield_avg6nonneg = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][0]
                    h_ue_yield_avg6nonneg_err = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["ue"][1]
                    h_total_yield_avg6nonneg = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][0]
                    h_total_yield_avg6nonneg_err = h_h_yield_extractor.yields[FitType.AVG_SIX_NONNEG]["total"][1]
                    h_near_yield_syserror = ((h_near_yield_avg4/h_near_yield) - 1)**2 + ((h_near_yield_avg6nonneg/h_near_yield) - 1)**2
                    h_near_yield_syserror = h_near_yield * math.sqrt(h_near_yield_syserror/2)
                    h_away_yield_syserror = ((h_away_yield_avg4/h_away_yield) - 1)**2 + ((h_away_yield_avg6nonneg/h_away_yield) - 1)**2
                    h_away_yield_syserror = h_away_yield * math.sqrt(h_away_yield_syserror/2)
                    h_ue_yield_syserror = ((h_ue_yield_avg4/h_ue_yield) - 1)**2 + ((h_ue_yield_avg6nonneg/h_ue_yield) - 1)**2
                    h_ue_yield_syserror = h_ue_yield * math.sqrt(h_ue_yield_syserror/2)
                    h_total_yield_syserror = ((h_total_yield_avg4/h_total_yield) - 1)**2 + ((h_total_yield_avg6nonneg/h_total_yield) - 1)**2
                    h_total_yield_syserror = h_total_yield * math.sqrt(h_total_yield_syserror/2)

                    print(f"----------CENT BIN {centrality_bin.name}-------ASSOCIATED PT BIN {assoc_pt_bin.name}----------")
                    print("h-lambda near-side width systematic error: " + str(h_lambda_near_width_syserr/h_lambda_near_width))
                    print("h-lambda away-side width systematic error: " + str(h_lambda_away_width_syserr/h_lambda_away_width))
                    print("h-h near-side width systematic error: " + str(h_h_near_width_syserr/h_h_near_width))
                    print("h-h away-side width systematic error: " + str(h_h_away_width_syserr/h_h_away_width))

                    print("h-lambda near-side yield systematic error: " + str(lambda_near_yield_syserror/lambda_near_yield))
                    print("h-lambda away-side yield systematic error: " + str(lambda_away_yield_syserror/lambda_away_yield))
                    print("h-lambda ue yield systematic error: " + str(lambda_ue_yield_syserror/lambda_ue_yield))
                    print("h-lambda total yield systematic error: " + str(lambda_total_yield_syserror/lambda_total_yield))

                    print("h-h near-side yield systematic error: " + str(h_near_yield_syserror/h_near_yield))
                    print("h-h away-side yield systematic error: " + str(h_away_yield_syserror/h_away_yield))
                    print("h-h ue yield systematic error: " + str(h_ue_yield_syserror/h_ue_yield))
                    print("h-h total yield systematic error: " + str(h_total_yield_syserror/h_total_yield))

                    print("h-phi near-side yield systematic error: " + str(phi_near_yield_syserror/phi_near_yield))
                    print("h-phi away-side yield systematic error: " + str(phi_away_yield_syserror/phi_away_yield))
                    print("h-phi ue yield systematic error: " + str(phi_ue_yield_syserror/phi_ue_yield))
                    print("h-phi total yield systematic error: " + str(phi_total_yield_syserror/phi_total_yield))


                    h_near_yield_err = math.sqrt(h_near_yield_err**2 + h_near_yield_syserror**2)
                    h_away_yield_err = math.sqrt(h_away_yield_err**2 + h_away_yield_syserror**2)
                    h_ue_yield_err = math.sqrt(h_ue_yield_err**2 + h_ue_yield_syserror**2)
                    h_total_yield_err = math.sqrt(h_total_yield_err**2 + h_total_yield_syserror**2)


                    h_lambda_near_yield_cent.append(lambda_near_yield)
                    h_lambda_near_yield_cent_err.append(lambda_near_yield_err)
                    h_lambda_away_yield_cent.append(lambda_away_yield)
                    h_lambda_away_yield_cent_err.append(lambda_away_yield_err)

                    h_h_near_yield_cent.append(h_near_yield)
                    h_h_near_yield_cent_err.append(h_near_yield_err)
                    h_h_away_yield_cent.append(h_away_yield)
                    h_h_away_yield_cent_err.append(h_away_yield_err)



                    lambda_hadron_near_ratio = lambda_near_yield/h_near_yield
                    lambda_hadron_near_ratio_err = lambda_hadron_near_ratio * math.sqrt((lambda_near_yield_err/lambda_near_yield)**2 + (h_near_yield_err/h_near_yield)**2)
                    lambda_hadron_near_ratio_cent.append(lambda_hadron_near_ratio)
                    lambda_hadron_near_ratio_cent_err.append(lambda_hadron_near_ratio_err)
                    lambda_hadron_away_ratio = lambda_away_yield/h_away_yield
                    lambda_hadron_away_ratio_err = lambda_hadron_away_ratio * math.sqrt((lambda_away_yield_err/lambda_away_yield)**2 + (h_away_yield_err/h_away_yield)**2)
                    lambda_hadron_away_ratio_cent.append(lambda_hadron_away_ratio)
                    lambda_hadron_away_ratio_cent_err.append(lambda_hadron_away_ratio_err)
                    lambda_hadron_ue_ratio = lambda_ue_yield/h_ue_yield
                    lambda_hadron_ue_ratio_err = lambda_hadron_ue_ratio * math.sqrt((lambda_ue_yield_err/lambda_ue_yield)**2 + (h_ue_yield_err/h_ue_yield)**2)
                    lambda_hadron_ue_ratio_cent.append(lambda_hadron_ue_ratio)
                    lambda_hadron_ue_ratio_cent_err.append(lambda_hadron_ue_ratio_err)
                    lambda_hadron_total_ratio = lambda_total_yield/h_total_yield
                    lambda_hadron_total_ratio_err = lambda_hadron_total_ratio * math.sqrt((lambda_total_yield_err/lambda_total_yield)**2 + (h_total_yield_err/h_total_yield)**2)
                    lambda_hadron_total_ratio_cent.append(lambda_hadron_total_ratio)
                    lambda_hadron_total_ratio_cent_err.append(lambda_hadron_total_ratio_err)

                    phi_hadron_near_ratio = phi_near_yield/h_near_yield
                    phi_hadron_near_ratio_err = phi_hadron_near_ratio * math.sqrt((phi_near_yield_err/phi_near_yield)**2 + (h_near_yield_err/h_near_yield)**2)
                    phi_hadron_near_ratio_cent.append(phi_hadron_near_ratio)
                    phi_hadron_near_ratio_cent_err.append(phi_hadron_near_ratio_err)
                    phi_hadron_away_ratio = phi_away_yield/h_away_yield
                    phi_hadron_away_ratio_err = phi_hadron_away_ratio * math.sqrt((phi_away_yield_err/phi_away_yield)**2 + (h_away_yield_err/h_away_yield)**2)
                    phi_hadron_away_ratio_cent.append(phi_hadron_away_ratio)
                    phi_hadron_away_ratio_cent_err.append(phi_hadron_away_ratio_err)
                    phi_hadron_ue_ratio = phi_ue_yield/h_ue_yield
                    phi_hadron_ue_ratio_err = phi_hadron_ue_ratio * math.sqrt((phi_ue_yield_err/phi_ue_yield)**2 + (h_ue_yield_err/h_ue_yield)**2)
                    phi_hadron_ue_ratio_cent.append(phi_hadron_ue_ratio)
                    phi_hadron_ue_ratio_cent_err.append(phi_hadron_ue_ratio_err)
                    phi_hadron_total_ratio = phi_total_yield/h_total_yield
                    phi_hadron_total_ratio_err = phi_hadron_total_ratio * math.sqrt((phi_total_yield_err/phi_total_yield)**2 + (h_total_yield_err/h_total_yield)**2)
                    phi_hadron_total_ratio_cent.append(phi_hadron_total_ratio)
                    phi_hadron_total_ratio_cent_err.append(phi_hadron_total_ratio_err)

                    lambda_phi_near_ratio = lambda_near_yield/phi_near_yield
                    lambda_phi_near_ratio_err = lambda_phi_near_ratio * math.sqrt((lambda_near_yield_err/lambda_near_yield)**2 + (phi_near_yield_err/phi_near_yield)**2)
                    lambda_phi_near_ratio_cent.append(lambda_phi_near_ratio)
                    lambda_phi_near_ratio_cent_err.append(lambda_phi_near_ratio_err)
                    lambda_phi_away_ratio = lambda_away_yield/phi_away_yield
                    lambda_phi_away_ratio_err = lambda_phi_away_ratio * math.sqrt((lambda_away_yield_err/lambda_away_yield)**2 + (phi_away_yield_err/phi_away_yield)**2)
                    lambda_phi_away_ratio_cent.append(lambda_phi_away_ratio)
                    lambda_phi_away_ratio_cent_err.append(lambda_phi_away_ratio_err)
                    lambda_phi_ue_ratio = lambda_ue_yield/phi_ue_yield
                    lambda_phi_ue_ratio_err = lambda_phi_ue_ratio * math.sqrt((lambda_ue_yield_err/lambda_ue_yield)**2 + (phi_ue_yield_err/phi_ue_yield)**2)
                    lambda_phi_ue_ratio_cent.append(lambda_phi_ue_ratio)
                    lambda_phi_ue_ratio_cent_err.append(lambda_phi_ue_ratio_err)
                    lambda_phi_total_ratio = lambda_total_yield/phi_total_yield
                    lambda_phi_total_ratio_err = lambda_phi_total_ratio * math.sqrt((lambda_total_yield_err/lambda_total_yield)**2 + (phi_total_yield_err/phi_total_yield)**2)
                    lambda_phi_total_ratio_cent.append(lambda_phi_total_ratio)
                    lambda_phi_total_ratio_cent_err.append(lambda_phi_total_ratio_err)



                mult_list = arr.array('d', [35, 65, 90])
                mult_error_list = arr.array('d', [15, 15, 10])
                nch_list = arr.array('d', [7, 20, 45])
                nch_error_list = arr.array('d', [0, 0, 0])

                h_lambda_near_width_cent = arr.array('d', h_lambda_near_width_cent)
                h_lambda_near_width_cent_err = arr.array('d', h_lambda_near_width_cent_err)
                h_lambda_away_width_cent = arr.array('d', h_lambda_away_width_cent)
                h_lambda_away_width_cent_err = arr.array('d', h_lambda_away_width_cent_err)

                h_phi_near_width_cent = arr.array('d', h_phi_near_width_cent)
                h_phi_near_width_cent_err = arr.array('d', h_phi_near_width_cent_err)
                h_phi_away_width_cent = arr.array('d', h_phi_away_width_cent)
                h_phi_away_width_cent_err = arr.array('d', h_phi_away_width_cent_err)

                h_h_near_width_cent = arr.array('d', h_h_near_width_cent)
                h_h_near_width_cent_err = arr.array('d', h_h_near_width_cent_err)
                h_h_away_width_cent = arr.array('d', h_h_away_width_cent)
                h_h_away_width_cent_err = arr.array('d', h_h_away_width_cent_err)

                h_lambda_near_yield_cent = arr.array('d', h_lambda_near_yield_cent)
                h_lambda_near_yield_cent_err = arr.array('d', h_lambda_near_yield_cent_err)
                h_lambda_away_yield_cent = arr.array('d', h_lambda_away_yield_cent)
                h_lambda_away_yield_cent_err = arr.array('d', h_lambda_away_yield_cent_err)

                h_h_near_yield_cent = arr.array('d', h_h_near_yield_cent)
                h_h_near_yield_cent_err = arr.array('d', h_h_near_yield_cent_err)
                h_h_away_yield_cent = arr.array('d', h_h_away_yield_cent)
                h_h_away_yield_cent_err = arr.array('d', h_h_away_yield_cent_err)


                lambda_hadron_near_ratio_cent = arr.array('d', lambda_hadron_near_ratio_cent)
                lambda_hadron_near_ratio_cent_err = arr.array('d', lambda_hadron_near_ratio_cent_err)
                lambda_hadron_away_ratio_cent = arr.array('d', lambda_hadron_away_ratio_cent)
                lambda_hadron_away_ratio_cent_err = arr.array('d', lambda_hadron_away_ratio_cent_err)
                lambda_hadron_ue_ratio_cent = arr.array('d', lambda_hadron_ue_ratio_cent)
                lambda_hadron_ue_ratio_cent_err = arr.array('d', lambda_hadron_ue_ratio_cent_err)
                lambda_hadron_total_ratio_cent = arr.array('d', lambda_hadron_total_ratio_cent)
                lambda_hadron_total_ratio_cent_err = arr.array('d', lambda_hadron_total_ratio_cent_err)

                phi_hadron_near_ratio_cent = arr.array('d', phi_hadron_near_ratio_cent)
                phi_hadron_near_ratio_cent_err = arr.array('d', phi_hadron_near_ratio_cent_err)
                phi_hadron_away_ratio_cent = arr.array('d', phi_hadron_away_ratio_cent)
                phi_hadron_away_ratio_cent_err = arr.array('d', phi_hadron_away_ratio_cent_err)
                phi_hadron_ue_ratio_cent = arr.array('d', phi_hadron_ue_ratio_cent)
                phi_hadron_ue_ratio_cent_err = arr.array('d', phi_hadron_ue_ratio_cent_err)
                phi_hadron_total_ratio_cent = arr.array('d', phi_hadron_total_ratio_cent)
                phi_hadron_total_ratio_cent_err = arr.array('d', phi_hadron_total_ratio_cent_err)

                lambda_phi_near_ratio_cent = arr.array('d', lambda_phi_near_ratio_cent)
                lambda_phi_near_ratio_cent_err = arr.array('d', lambda_phi_near_ratio_cent_err)
                lambda_phi_away_ratio_cent = arr.array('d', lambda_phi_away_ratio_cent)
                lambda_phi_away_ratio_cent_err = arr.array('d', lambda_phi_away_ratio_cent_err)
                lambda_phi_ue_ratio_cent = arr.array('d', lambda_phi_ue_ratio_cent)
                lambda_phi_ue_ratio_cent_err = arr.array('d', lambda_phi_ue_ratio_cent_err)
                lambda_phi_total_ratio_cent = arr.array('d', lambda_phi_total_ratio_cent)
                lambda_phi_total_ratio_cent_err = arr.array('d', lambda_phi_total_ratio_cent_err)

                h_lambda_near_width_graph = rt.TGraphErrors(3, mult_list, h_lambda_near_width_cent, mult_error_list, h_lambda_near_width_cent_err)
                h_lambda_near_width_graph.SetMarkerStyle(20)
                h_lambda_near_width_graph.SetMarkerSize(1)
                h_lambda_near_width_graph.SetMarkerColor(rt.kRed+1)
                h_lambda_near_width_graph.SetLineColor(rt.kRed+2)
                h_lambda_near_width_graph.SetLineWidth(2)
                h_lambda_away_width_graph = rt.TGraphErrors(3, mult_list, h_lambda_away_width_cent, mult_error_list, h_lambda_away_width_cent_err)
                h_lambda_away_width_graph.SetMarkerStyle(20)
                h_lambda_away_width_graph.SetMarkerSize(1)
                h_lambda_away_width_graph.SetMarkerColor(rt.kBlue+1)
                h_lambda_away_width_graph.SetLineColor(rt.kBlue+2)
                h_lambda_away_width_graph.SetLineWidth(2)

                h_lambda_near_yield_graph = rt.TGraphErrors(3, mult_list, h_lambda_near_yield_cent, mult_error_list, h_lambda_near_yield_cent_err)
                h_lambda_near_yield_graph.SetMarkerStyle(20)
                h_lambda_near_yield_graph.SetMarkerSize(1)
                h_lambda_near_yield_graph.SetMarkerColor(rt.kRed+1)
                h_lambda_near_yield_graph.SetLineColor(rt.kRed+2)
                h_lambda_near_yield_graph.SetLineWidth(2)
                h_lambda_away_yield_graph = rt.TGraphErrors(3, mult_list, h_lambda_away_yield_cent, mult_error_list, h_lambda_away_yield_cent_err)
                h_lambda_away_yield_graph.SetMarkerStyle(20)
                h_lambda_away_yield_graph.SetMarkerSize(1)
                h_lambda_away_yield_graph.SetMarkerColor(rt.kBlue+1)
                h_lambda_away_yield_graph.SetLineColor(rt.kBlue+2)
                h_lambda_away_yield_graph.SetLineWidth(2)

                h_lambda_near_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_lambda_near_width_cent, nch_error_list, h_lambda_near_width_cent_err)
                h_lambda_near_width_graph_new_x_axis.SetMarkerStyle(20)
                h_lambda_near_width_graph_new_x_axis.SetMarkerSize(1)
                h_lambda_near_width_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                h_lambda_near_width_graph_new_x_axis.SetLineColor(rt.kRed+2)
                h_lambda_near_width_graph_new_x_axis.SetLineWidth(2)
                h_lambda_away_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_lambda_away_width_cent, nch_error_list, h_lambda_away_width_cent_err)
                h_lambda_away_width_graph_new_x_axis.SetMarkerStyle(20)
                h_lambda_away_width_graph_new_x_axis.SetMarkerSize(1)
                h_lambda_away_width_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                h_lambda_away_width_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                h_lambda_away_width_graph_new_x_axis.SetLineWidth(2)

                h_lambda_near_yield_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_lambda_near_yield_cent, nch_error_list, h_lambda_near_yield_cent_err)
                h_lambda_near_yield_graph_new_x_axis.SetMarkerStyle(20)
                h_lambda_near_yield_graph_new_x_axis.SetMarkerSize(1)
                h_lambda_near_yield_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                h_lambda_near_yield_graph_new_x_axis.SetLineColor(rt.kRed+2)
                h_lambda_near_yield_graph_new_x_axis.SetLineWidth(2)
                h_lambda_away_yield_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_lambda_away_yield_cent, nch_error_list, h_lambda_away_yield_cent_err)
                h_lambda_away_yield_graph_new_x_axis.SetMarkerStyle(20)
                h_lambda_away_yield_graph_new_x_axis.SetMarkerSize(1)
                h_lambda_away_yield_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                h_lambda_away_yield_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                h_lambda_away_yield_graph_new_x_axis.SetLineWidth(2)

                h_phi_near_width_graph = rt.TGraphErrors(3, mult_list, h_phi_near_width_cent, mult_error_list, h_phi_near_width_cent_err)
                h_phi_near_width_graph.SetMarkerStyle(20)
                h_phi_near_width_graph.SetMarkerSize(1)
                h_phi_near_width_graph.SetMarkerColor(rt.kRed+1)
                h_phi_near_width_graph.SetLineColor(rt.kRed+2)
                h_phi_near_width_graph.SetLineWidth(2)
                h_phi_away_width_graph = rt.TGraphErrors(3, mult_list, h_phi_away_width_cent, mult_error_list, h_phi_away_width_cent_err)
                h_phi_away_width_graph.SetMarkerStyle(20)
                h_phi_away_width_graph.SetMarkerSize(1)
                h_phi_away_width_graph.SetMarkerColor(rt.kBlue+1)
                h_phi_away_width_graph.SetLineColor(rt.kBlue+2)
                h_phi_away_width_graph.SetLineWidth(2)

                h_phi_near_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_phi_near_width_cent, nch_error_list, h_phi_near_width_cent_err)
                h_phi_near_width_graph_new_x_axis.SetMarkerStyle(20)
                h_phi_near_width_graph_new_x_axis.SetMarkerSize(1)
                h_phi_near_width_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                h_phi_near_width_graph_new_x_axis.SetLineColor(rt.kRed+2)
                h_phi_near_width_graph_new_x_axis.SetLineWidth(2)
                h_phi_away_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_phi_away_width_cent, nch_error_list, h_phi_away_width_cent_err)
                h_phi_away_width_graph_new_x_axis.SetMarkerStyle(20)
                h_phi_away_width_graph_new_x_axis.SetMarkerSize(1)
                h_phi_away_width_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                h_phi_away_width_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                h_phi_away_width_graph_new_x_axis.SetLineWidth(2)

                h_h_near_width_graph = rt.TGraphErrors(3, mult_list, h_h_near_width_cent, mult_error_list, h_h_near_width_cent_err)
                h_h_near_width_graph.SetMarkerStyle(20)
                h_h_near_width_graph.SetMarkerSize(1)
                h_h_near_width_graph.SetMarkerColor(rt.kRed+1)
                h_h_near_width_graph.SetLineColor(rt.kRed+2)
                h_h_near_width_graph.SetLineWidth(2)
                h_h_away_width_graph = rt.TGraphErrors(3, mult_list, h_h_away_width_cent, mult_error_list, h_h_away_width_cent_err)
                h_h_away_width_graph.SetMarkerStyle(20)
                h_h_away_width_graph.SetMarkerSize(1)
                h_h_away_width_graph.SetMarkerColor(rt.kBlue+1)
                h_h_away_width_graph.SetLineColor(rt.kBlue+2)
                h_h_away_width_graph.SetLineWidth(2)

                h_h_near_yield_graph = rt.TGraphErrors(3, mult_list, h_h_near_yield_cent, mult_error_list, h_h_near_yield_cent_err)
                h_h_near_yield_graph.SetMarkerStyle(20)
                h_h_near_yield_graph.SetMarkerSize(1)
                h_h_near_yield_graph.SetMarkerColor(rt.kRed+1)
                h_h_near_yield_graph.SetLineColor(rt.kRed+2)
                h_h_near_yield_graph.SetLineWidth(2)
                h_h_away_yield_graph = rt.TGraphErrors(3, mult_list, h_h_away_yield_cent, mult_error_list, h_h_away_yield_cent_err)
                h_h_away_yield_graph.SetMarkerStyle(20)
                h_h_away_yield_graph.SetMarkerSize(1)
                h_h_away_yield_graph.SetMarkerColor(rt.kBlue+1)
                h_h_away_yield_graph.SetLineColor(rt.kBlue+2)
                h_h_away_yield_graph.SetLineWidth(2)

                h_h_near_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_h_near_width_cent, nch_error_list, h_h_near_width_cent_err)
                h_h_near_width_graph_new_x_axis.SetMarkerStyle(20)
                h_h_near_width_graph_new_x_axis.SetMarkerSize(1)
                h_h_near_width_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                h_h_near_width_graph_new_x_axis.SetLineColor(rt.kRed+2)
                h_h_near_width_graph_new_x_axis.SetLineWidth(2)
                h_h_away_width_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_h_away_width_cent, nch_error_list, h_h_away_width_cent_err)
                h_h_away_width_graph_new_x_axis.SetMarkerStyle(20)
                h_h_away_width_graph_new_x_axis.SetMarkerSize(1)
                h_h_away_width_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                h_h_away_width_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                h_h_away_width_graph_new_x_axis.SetLineWidth(2)

                h_h_near_yield_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_h_near_yield_cent, nch_error_list, h_h_near_yield_cent_err)
                h_h_near_yield_graph_new_x_axis.SetMarkerStyle(20)
                h_h_near_yield_graph_new_x_axis.SetMarkerSize(1)
                h_h_near_yield_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                h_h_near_yield_graph_new_x_axis.SetLineColor(rt.kRed+2)
                h_h_near_yield_graph_new_x_axis.SetLineWidth(2)
                h_h_away_yield_graph_new_x_axis = rt.TGraphErrors(3, nch_list, h_h_away_yield_cent, nch_error_list, h_h_away_yield_cent_err)
                h_h_away_yield_graph_new_x_axis.SetMarkerStyle(20)
                h_h_away_yield_graph_new_x_axis.SetMarkerSize(1)
                h_h_away_yield_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                h_h_away_yield_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                h_h_away_yield_graph_new_x_axis.SetLineWidth(2)

                lambda_hadron_near_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_hadron_near_ratio_cent, mult_error_list, lambda_hadron_near_ratio_cent_err)
                lambda_hadron_near_ratio_graph.SetMarkerStyle(20)
                lambda_hadron_near_ratio_graph.SetMarkerSize(1)
                lambda_hadron_near_ratio_graph.SetMarkerColor(rt.kRed+1)
                lambda_hadron_near_ratio_graph.SetLineColor(rt.kRed+2)
                lambda_hadron_near_ratio_graph.SetLineWidth(2)
                lambda_hadron_away_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_hadron_away_ratio_cent, mult_error_list, lambda_hadron_away_ratio_cent_err)
                lambda_hadron_away_ratio_graph.SetMarkerStyle(20)
                lambda_hadron_away_ratio_graph.SetMarkerSize(1)
                lambda_hadron_away_ratio_graph.SetMarkerColor(rt.kBlue+1)
                lambda_hadron_away_ratio_graph.SetLineColor(rt.kBlue+2)
                lambda_hadron_away_ratio_graph.SetLineWidth(2)
                lambda_hadron_ue_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_hadron_ue_ratio_cent, mult_error_list, lambda_hadron_ue_ratio_cent_err)
                lambda_hadron_ue_ratio_graph.SetMarkerStyle(20)
                lambda_hadron_ue_ratio_graph.SetMarkerSize(1)
                lambda_hadron_ue_ratio_graph.SetMarkerColor(rt.kGreen+1)
                lambda_hadron_ue_ratio_graph.SetLineColor(rt.kGreen+2)
                lambda_hadron_ue_ratio_graph.SetLineWidth(2)
                lambda_hadron_total_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_hadron_total_ratio_cent, mult_error_list, lambda_hadron_total_ratio_cent_err)
                lambda_hadron_total_ratio_graph.SetMarkerStyle(20)
                lambda_hadron_total_ratio_graph.SetMarkerSize(1)
                lambda_hadron_total_ratio_graph.SetMarkerColor(rt.kMagenta+1)
                lambda_hadron_total_ratio_graph.SetLineColor(rt.kMagenta+2)
                lambda_hadron_total_ratio_graph.SetLineWidth(2)

                lambda_hadron_near_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_hadron_near_ratio_cent, nch_error_list, lambda_hadron_near_ratio_cent_err)
                lambda_hadron_near_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_hadron_near_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_hadron_near_ratio_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                lambda_hadron_near_ratio_graph_new_x_axis.SetLineColor(rt.kRed+2)
                lambda_hadron_near_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_hadron_away_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_hadron_away_ratio_cent, nch_error_list, lambda_hadron_away_ratio_cent_err)
                lambda_hadron_away_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_hadron_away_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_hadron_away_ratio_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                lambda_hadron_away_ratio_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                lambda_hadron_away_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_hadron_ue_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_hadron_ue_ratio_cent, nch_error_list, lambda_hadron_ue_ratio_cent_err)
                lambda_hadron_ue_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_hadron_ue_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_hadron_ue_ratio_graph_new_x_axis.SetMarkerColor(rt.kGreen+1)
                lambda_hadron_ue_ratio_graph_new_x_axis.SetLineColor(rt.kGreen+2)
                lambda_hadron_ue_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_hadron_total_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_hadron_total_ratio_cent, nch_error_list, lambda_hadron_total_ratio_cent_err)
                lambda_hadron_total_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_hadron_total_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_hadron_total_ratio_graph_new_x_axis.SetMarkerColor(rt.kMagenta+1)
                lambda_hadron_total_ratio_graph_new_x_axis.SetLineColor(rt.kMagenta+2)
                lambda_hadron_total_ratio_graph_new_x_axis.SetLineWidth(2)

                phi_hadron_near_ratio_graph = rt.TGraphErrors(3, mult_list, phi_hadron_near_ratio_cent, mult_error_list, phi_hadron_near_ratio_cent_err)
                phi_hadron_near_ratio_graph.SetMarkerStyle(20)
                phi_hadron_near_ratio_graph.SetMarkerSize(1)
                phi_hadron_near_ratio_graph.SetMarkerColor(rt.kRed+1)
                phi_hadron_near_ratio_graph.SetLineColor(rt.kRed+2)
                phi_hadron_near_ratio_graph.SetLineWidth(2)
                phi_hadron_away_ratio_graph = rt.TGraphErrors(3, mult_list, phi_hadron_away_ratio_cent, mult_error_list, phi_hadron_away_ratio_cent_err)
                phi_hadron_away_ratio_graph.SetMarkerStyle(20)
                phi_hadron_away_ratio_graph.SetMarkerSize(1)
                phi_hadron_away_ratio_graph.SetMarkerColor(rt.kBlue+1)
                phi_hadron_away_ratio_graph.SetLineColor(rt.kBlue+2)
                phi_hadron_away_ratio_graph.SetLineWidth(2)
                phi_hadron_ue_ratio_graph = rt.TGraphErrors(3, mult_list, phi_hadron_ue_ratio_cent, mult_error_list, phi_hadron_ue_ratio_cent_err)
                phi_hadron_ue_ratio_graph.SetMarkerStyle(20)
                phi_hadron_ue_ratio_graph.SetMarkerSize(1)
                phi_hadron_ue_ratio_graph.SetMarkerColor(rt.kGreen+1)
                phi_hadron_ue_ratio_graph.SetLineColor(rt.kGreen+2)
                phi_hadron_ue_ratio_graph.SetLineWidth(2)
                phi_hadron_total_ratio_graph = rt.TGraphErrors(3, mult_list, phi_hadron_total_ratio_cent, mult_error_list, phi_hadron_total_ratio_cent_err)
                phi_hadron_total_ratio_graph.SetMarkerStyle(20)
                phi_hadron_total_ratio_graph.SetMarkerSize(1)
                phi_hadron_total_ratio_graph.SetMarkerColor(rt.kMagenta+1)
                phi_hadron_total_ratio_graph.SetLineColor(rt.kMagenta+2)
                phi_hadron_total_ratio_graph.SetLineWidth(2)

                phi_hadron_near_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, phi_hadron_near_ratio_cent, nch_error_list, phi_hadron_near_ratio_cent_err)
                phi_hadron_near_ratio_graph_new_x_axis.SetMarkerStyle(20)
                phi_hadron_near_ratio_graph_new_x_axis.SetMarkerSize(1)
                phi_hadron_near_ratio_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                phi_hadron_near_ratio_graph_new_x_axis.SetLineColor(rt.kRed+2)
                phi_hadron_near_ratio_graph_new_x_axis.SetLineWidth(2)
                phi_hadron_away_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, phi_hadron_away_ratio_cent, nch_error_list, phi_hadron_away_ratio_cent_err)
                phi_hadron_away_ratio_graph_new_x_axis.SetMarkerStyle(20)
                phi_hadron_away_ratio_graph_new_x_axis.SetMarkerSize(1)
                phi_hadron_away_ratio_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                phi_hadron_away_ratio_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                phi_hadron_away_ratio_graph_new_x_axis.SetLineWidth(2)
                phi_hadron_ue_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, phi_hadron_ue_ratio_cent, nch_error_list, phi_hadron_ue_ratio_cent_err)
                phi_hadron_ue_ratio_graph_new_x_axis.SetMarkerStyle(20)
                phi_hadron_ue_ratio_graph_new_x_axis.SetMarkerSize(1)
                phi_hadron_ue_ratio_graph_new_x_axis.SetMarkerColor(rt.kGreen+1)
                phi_hadron_ue_ratio_graph_new_x_axis.SetLineColor(rt.kGreen+2)
                phi_hadron_ue_ratio_graph_new_x_axis.SetLineWidth(2)
                phi_hadron_total_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, phi_hadron_total_ratio_cent, nch_error_list, phi_hadron_total_ratio_cent_err)
                phi_hadron_total_ratio_graph_new_x_axis.SetMarkerStyle(20)
                phi_hadron_total_ratio_graph_new_x_axis.SetMarkerSize(1)
                phi_hadron_total_ratio_graph_new_x_axis.SetMarkerColor(rt.kMagenta+1)
                phi_hadron_total_ratio_graph_new_x_axis.SetLineColor(rt.kMagenta+2)
                phi_hadron_total_ratio_graph_new_x_axis.SetLineWidth(2)

                lambda_phi_near_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_phi_near_ratio_cent, mult_error_list, lambda_phi_near_ratio_cent_err)
                lambda_phi_near_ratio_graph.SetMarkerStyle(20)
                lambda_phi_near_ratio_graph.SetMarkerSize(1)
                lambda_phi_near_ratio_graph.SetMarkerColor(rt.kRed+1)
                lambda_phi_near_ratio_graph.SetLineColor(rt.kRed+2)
                lambda_phi_near_ratio_graph.SetLineWidth(2)
                lambda_phi_away_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_phi_away_ratio_cent, mult_error_list, lambda_phi_away_ratio_cent_err)
                lambda_phi_away_ratio_graph.SetMarkerStyle(20)
                lambda_phi_away_ratio_graph.SetMarkerSize(1)
                lambda_phi_away_ratio_graph.SetMarkerColor(rt.kBlue+1)
                lambda_phi_away_ratio_graph.SetLineColor(rt.kBlue+2)
                lambda_phi_away_ratio_graph.SetLineWidth(2)
                lambda_phi_ue_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_phi_ue_ratio_cent, mult_error_list, lambda_phi_ue_ratio_cent_err)
                lambda_phi_ue_ratio_graph.SetMarkerStyle(20)
                lambda_phi_ue_ratio_graph.SetMarkerSize(1)
                lambda_phi_ue_ratio_graph.SetMarkerColor(rt.kGreen+1)
                lambda_phi_ue_ratio_graph.SetLineColor(rt.kGreen+2)
                lambda_phi_ue_ratio_graph.SetLineWidth(2)
                lambda_phi_total_ratio_graph = rt.TGraphErrors(3, mult_list, lambda_phi_total_ratio_cent, mult_error_list, lambda_phi_total_ratio_cent_err)
                lambda_phi_total_ratio_graph.SetMarkerStyle(20)
                lambda_phi_total_ratio_graph.SetMarkerSize(1)
                lambda_phi_total_ratio_graph.SetMarkerColor(rt.kMagenta+1)
                lambda_phi_total_ratio_graph.SetLineColor(rt.kMagenta+2)
                lambda_phi_total_ratio_graph.SetLineWidth(2)

                lambda_phi_near_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_phi_near_ratio_cent, nch_error_list, lambda_phi_near_ratio_cent_err)
                lambda_phi_near_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_phi_near_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_phi_near_ratio_graph_new_x_axis.SetMarkerColor(rt.kRed+1)
                lambda_phi_near_ratio_graph_new_x_axis.SetLineColor(rt.kRed+2)
                lambda_phi_near_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_phi_away_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_phi_away_ratio_cent, nch_error_list, lambda_phi_away_ratio_cent_err)
                lambda_phi_away_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_phi_away_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_phi_away_ratio_graph_new_x_axis.SetMarkerColor(rt.kBlue+1)
                lambda_phi_away_ratio_graph_new_x_axis.SetLineColor(rt.kBlue+2)
                lambda_phi_away_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_phi_ue_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_phi_ue_ratio_cent, nch_error_list, lambda_phi_ue_ratio_cent_err)
                lambda_phi_ue_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_phi_ue_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_phi_ue_ratio_graph_new_x_axis.SetMarkerColor(rt.kGreen+1)
                lambda_phi_ue_ratio_graph_new_x_axis.SetLineColor(rt.kGreen+2)
                lambda_phi_ue_ratio_graph_new_x_axis.SetLineWidth(2)
                lambda_phi_total_ratio_graph_new_x_axis = rt.TGraphErrors(3, nch_list, lambda_phi_total_ratio_cent, nch_error_list, lambda_phi_total_ratio_cent_err)
                lambda_phi_total_ratio_graph_new_x_axis.SetMarkerStyle(20)
                lambda_phi_total_ratio_graph_new_x_axis.SetMarkerSize(1)
                lambda_phi_total_ratio_graph_new_x_axis.SetMarkerColor(rt.kMagenta+1)
                lambda_phi_total_ratio_graph_new_x_axis.SetLineColor(rt.kMagenta+2)
                lambda_phi_total_ratio_graph_new_x_axis.SetLineWidth(2)

                outfile.cd()

                h_lambda_near_width_graph.Write(f"h_lambda_near_width_graph")
                h_lambda_away_width_graph.Write(f"h_lambda_away_width_graph")
                h_lambda_near_width_graph_new_x_axis.Write(f"h_lambda_near_width_graph_new_x_axis")
                h_lambda_away_width_graph_new_x_axis.Write(f"h_lambda_away_width_graph_new_x_axis")

                h_lambda_near_yield_graph.Write(f"h_lambda_near_yield_graph")
                h_lambda_away_yield_graph.Write(f"h_lambda_away_yield_graph")
                h_lambda_near_yield_graph_new_x_axis.Write(f"h_lambda_near_yield_graph_new_x_axis")
                h_lambda_away_yield_graph_new_x_axis.Write(f"h_lambda_away_yield_graph_new_x_axis")

                h_phi_near_width_graph.Write(f"h_phi_near_width_graph")
                h_phi_away_width_graph.Write(f"h_phi_away_width_graph")
                h_phi_near_width_graph_new_x_axis.Write(f"h_phi_near_width_graph_new_x_axis")
                h_phi_away_width_graph_new_x_axis.Write(f"h_phi_away_width_graph_new_x_axis")

                h_h_near_width_graph.Write(f"h_h_near_width_graph")
                h_h_away_width_graph.Write(f"h_h_away_width_graph")
                h_h_near_width_graph_new_x_axis.Write(f"h_h_near_width_graph_new_x_axis")
                h_h_away_width_graph_new_x_axis.Write(f"h_h_away_width_graph_new_x_axis")

                h_h_near_yield_graph.Write(f"h_h_near_yield_graph")
                h_h_away_yield_graph.Write(f"h_h_away_yield_graph")
                h_h_near_yield_graph_new_x_axis.Write(f"h_h_near_yield_graph_new_x_axis")
                h_h_away_yield_graph_new_x_axis.Write(f"h_h_away_yield_graph_new_x_axis")

                lambda_hadron_near_ratio_graph.Write(f"lambda_hadron_near_ratio_graph")
                lambda_hadron_away_ratio_graph.Write(f"lambda_hadron_away_ratio_graph")
                lambda_hadron_ue_ratio_graph.Write(f"lambda_hadron_ue_ratio_graph")
                lambda_hadron_total_ratio_graph.Write(f"lambda_hadron_total_ratio_graph")
                phi_hadron_near_ratio_graph.Write(f"phi_hadron_near_ratio_graph")
                phi_hadron_away_ratio_graph.Write(f"phi_hadron_away_ratio_graph")
                phi_hadron_ue_ratio_graph.Write(f"phi_hadron_ue_ratio_graph")
                phi_hadron_total_ratio_graph.Write(f"phi_hadron_total_ratio_graph")
                lambda_phi_near_ratio_graph.Write(f"lambda_phi_near_ratio_graph")
                lambda_phi_away_ratio_graph.Write(f"lambda_phi_away_ratio_graph")
                lambda_phi_ue_ratio_graph.Write(f"lambda_phi_ue_ratio_graph")
                lambda_phi_total_ratio_graph.Write(f"lambda_phi_total_ratio_graph")

                lambda_hadron_near_ratio_graph_new_x_axis.Write(f"lambda_hadron_near_ratio_graph_new_x_axis")
                lambda_hadron_away_ratio_graph_new_x_axis.Write(f"lambda_hadron_away_ratio_graph_new_x_axis")
                lambda_hadron_ue_ratio_graph_new_x_axis.Write(f"lambda_hadron_ue_ratio_graph_new_x_axis")
                lambda_hadron_total_ratio_graph_new_x_axis.Write(f"lambda_hadron_total_ratio_graph_new_x_axis")
                phi_hadron_near_ratio_graph_new_x_axis.Write(f"phi_hadron_near_ratio_graph_new_x_axis")
                phi_hadron_away_ratio_graph_new_x_axis.Write(f"phi_hadron_away_ratio_graph_new_x_axis")
                phi_hadron_ue_ratio_graph_new_x_axis.Write(f"phi_hadron_ue_ratio_graph_new_x_axis")
                phi_hadron_total_ratio_graph_new_x_axis.Write(f"phi_hadron_total_ratio_graph_new_x_axis")
                lambda_phi_near_ratio_graph_new_x_axis.Write(f"lambda_phi_near_ratio_graph_new_x_axis")
                lambda_phi_away_ratio_graph_new_x_axis.Write(f"lambda_phi_away_ratio_graph_new_x_axis")
                lambda_phi_ue_ratio_graph_new_x_axis.Write(f"lambda_phi_ue_ratio_graph_new_x_axis")
                lambda_phi_total_ratio_graph_new_x_axis.Write(f"lambda_phi_total_ratio_graph_new_x_axis")


            #     legend = rt.TLegend(0.65, 0.65, 0.85, 0.85)
            #     legend.SetBorderSize(0)
            #     legend.SetFillStyle(0)
            #     legend.AddEntry(lambda_hadron_near_ratio_graph, "Near-side (jet)", "lep")
            #     legend.AddEntry(lambda_hadron_away_ratio_graph, "Away-side (jet)", "lep")
            #     legend.AddEntry(lambda_hadron_ue_ratio_graph, "Underlying Event", "lep")
            #     # legend.AddEntry(lambda_hadron_total_ratio_graph, "Total", "lep")

            #     for lambda_phi_ratio in [True, False]:
            #         mult_bin_widths = arr.array('d', [0.0, 20.0, 50.0, 80.0, 100.0])
            #         plotting_hist = rt.TH1D("plotting_hist", "", 4, mult_bin_widths)
            #         plotting_hist.SetMarkerStyle(20)
            #         plotting_hist.SetMarkerSize(1)
            #         plotting_hist.SetMarkerColor(rt.kRed+1)
            #         plotting_hist.SetLineColor(rt.kRed+2)
            #         plotting_hist.SetLineWidth(2)
            #         plotting_hist.GetXaxis().SetTitle("Multiplicity (%)")
            #         plotting_hist.GetXaxis().SetTitleSize(0.05)
            #         plotting_hist.GetXaxis().SetLabelSize(0.045)
            #         plotting_hist.GetYaxis().SetLabelSize(0.045)
            #         plotting_hist.GetXaxis().SetTitleOffset(1.2)
            #         plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
            #         plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
            #         plotting_hist.SetStats(0)
            #         if lambda_phi_ratio:
            #             plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minus#phi}#right)")
            #             plotting_hist.GetYaxis().SetRangeUser(0, 100)
            #         else:
            #             plotting_hist.GetYaxis().SetRangeUser(0.0, 0.05)
            #             plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#phi}{h#minush}#right)")
            #         plotting_hist.GetYaxis().SetTitleSize(0.05)
            #         plotting_hist.GetYaxis().SetTitleOffset(1.3)
            #         plotting_hist.SetStats(0)
            #         plotting_hist.GetXaxis().SetLabelOffset(999)
            #         plotting_hist.Draw("PE")
            #         plotting_hist.GetXaxis().SetTickSize(0)
            #         rt.gPad.Update()
            #         new_axis = rt.TGaxis(rt.gPad.GetUxmax(),
            #                 rt.gPad.GetUymin(),
            #                 rt.gPad.GetUxmin(),
            #                 rt.gPad.GetUymin(),
            #                 plotting_hist.GetXaxis().GetXmin(),
            #                 plotting_hist.GetXaxis().GetXmax(),
            #                 510,"-")
            #         new_axis.SetLabelSize(0.045)
            #         new_axis.SetLabelFont(plotting_hist.GetXaxis().GetLabelFont())
            #         new_axis.SetLabelOffset(-0.03)
            #         new_axis.Draw()

            #         if lambda_phi_ratio:
            #             lambda_phi_near_ratio_graph.Draw("PZ SAME")
            #             lambda_phi_away_ratio_graph.Draw("PZ SAME")
            #             lambda_phi_ue_ratio_graph.Draw("PZ SAME")
            #             legend.Draw("SAME")
            #             c.SaveAs("figures/lambda_phi_ratio_" + model.name + "_" + assoc_pt_bin.name + ".pdf")
            #         else:
            #             plotting_hist.GetYaxis().SetRangeUser(0.0, 0.08)
            #             plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#phi}{h#minush}#right)")
            #             phi_hadron_near_ratio_graph.Draw("PZ SAME")
            #             phi_hadron_away_ratio_graph.Draw("PZ SAME")
            #             phi_hadron_ue_ratio_graph.Draw("PZ SAME")
            #             legend.Draw("SAME")
            #             c.SaveAs("figures/phi_hadron_ratio_" + model.name + "_" + assoc_pt_bin.name + ".pdf")

            #             plotting_hist.GetYaxis().SetRangeUser(0.0, 0.1)
            #             plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minush}#right)")
            #             lambda_hadron_near_ratio_graph.Draw("PZ SAME")
            #             lambda_hadron_away_ratio_graph.Draw("PZ SAME")
            #             lambda_hadron_ue_ratio_graph.Draw("PZ SAME")
            #             legend.Draw("SAME")
            #             c.SaveAs("figures/lambda_hadron_ratio_" + model.name + "_" + assoc_pt_bin.name + ".pdf")
            outfile.Close()
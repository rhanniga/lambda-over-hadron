import ROOT as rt
import array as arr
import math

from strangehelper import make_mixed_corrections

rt.gStyle.SetOptStat(0)

EPSILON = 0.00001

# decide whether to do sideband subtraction or not
DO_SIDEBAND_SUBTRACTION = False
USE_MC_KINEMATICS = True

# ETA CUTS
ETA_MIN = -0.8
ETA_MAX = 0.8 - EPSILON
DELTA_ETA_MIN = -1.2
DELTA_ETA_MAX = 1.2 - EPSILON 

TRIG_PT_LOW = 4
TRIG_PT_HIGH = 8 - EPSILON

for PT_MODE in [0, 1, 2]:
    if PT_MODE == 0:
        ASSOC_PT_LOW = 1.0
        ASSOC_PT_HIGH = 4.0 - EPSILON
    elif PT_MODE == 1:
        ASSOC_PT_LOW = 1.5
        ASSOC_PT_HIGH = 2.5 - EPSILON
    elif PT_MODE == 2:
        ASSOC_PT_LOW = 2.5
        ASSOC_PT_HIGH = 4.0 - EPSILON

    # SIGNAL AND SIDEBAND REGION CUTS
    SIG_MIN = 1.102
    SIG_MAX = 1.130 - EPSILON
    RSB_MIN = 1.135
    RSB_MAX = 1.50 - EPSILON

    c = rt.TCanvas("main_canvas", "Main Canvas", 55, 55, 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)

    # input_file = rt.TFile("../online/closure/v0/output/closure_full_stats.root")
    # input_file = rt.TFile("../online/closure/v0/output/closure_18f3_FAST_pt1.root")
    input_file = rt.TFile("../online/closure/v0/output/closure_20f11c2_FAST.root")
    input_list = input_file.Get("h-lambda")

    trig_dist = input_list.FindObject("fTriggerDistEff")
    lambda_dist = input_list.FindObject("fTriggeredLambdaDist")
    trig_dist_mc = input_list.FindObject("fTriggerDist_MC")
    lambda_dist_mc = input_list.FindObject("fTriggeredLambdaDist_MC")

    h_h = input_list.FindObject("fDphiHHEff")
    h_h_mixed = input_list.FindObject("fDphiHHMixed")
    h_h_mc = input_list.FindObject("fDphiHH_MC")
    h_h_mixed_mc = input_list.FindObject("fDphiHHMixed_MC")

    h_lambda = input_list.FindObject("fDphiHLambdaEff_MCKin")
    h_lambda_mixed = input_list.FindObject("fDphiHLambdaMixed_MCKin")
    h_lambda_mc = input_list.FindObject("fDphiHLambda_MC")
    h_lambda_mixed_mc = input_list.FindObject("fDphiHLambdaMixed_MC")

    trig_dist.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_h.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_h_mixed.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)

    h_lambda.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_lambda_mixed.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)

    trig_dist_mc.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_h_mc.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_h_mixed_mc.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)

    h_lambda_mc.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
    h_lambda_mixed_mc.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)

    lambda_dist.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_h.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_h_mixed.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_lambda.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_lambda_mixed.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    lambda_dist_mc.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_h_mc.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_h_mixed_mc.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_lambda_mc.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
    h_lambda_mixed_mc.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)

    trig_dist.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
    trig_dist_mc.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
    lambda_dist.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
    lambda_dist_mc.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)

    trig_pt_dist = trig_dist.Projection(0).Clone("trig_pt_dist")
    trig_phi_dist = trig_dist.Projection(1).Clone("trig_phi_dist")
    trig_eta_dist = trig_dist.Projection(2).Clone("trig_eta_dist")
    trig_2d_dist = trig_dist.Projection(0, 3).Clone("trig_2d_dist")

    trig_pt_dist.SetTitle("Trigger #font[12]{p}_{T} Distribution (red: reconstruced, blue: MC truth)")
    trig_pt_dist.Sumw2()

    trig_phi_dist.SetTitle("Trigger #varphi Distribution (red: reconstructed, blue: MC truth")
    trig_phi_dist.Sumw2()

    trig_eta_dist.SetTitle("Trigger #eta Distribution (red: reconstructed, blue: MC truth")
    trig_eta_dist.Sumw2()

    trig_pt_dist_mc = trig_dist_mc.Projection(0).Clone("trig_pt_dist_mc")
    trig_phi_dist_mc = trig_dist_mc.Projection(1).Clone("trig_phi_dist_mc")
    trig_eta_dist_mc = trig_dist_mc.Projection(2).Clone("trig_eta_dist_mc")
    trig_2d_dist_mc = trig_dist_mc.Projection(0, 3).Clone("trig_2d_dist_mc")


    trig_pt_dist_mc.SetTitle("Trigger #font[12]{p}_{T} Distribution (red: reconstruced, blue: MC truth)")
    trig_pt_dist_mc.Sumw2()

    trig_phi_dist_mc.SetTitle("Trigger #varphi Distribution (red: reconstructed, blue: MC truth")
    trig_phi_dist_mc.Sumw2()

    trig_eta_dist_mc.SetTitle("Trigger #eta Distribution (red: reconstructed, blue: MC truth")
    trig_eta_dist_mc.Sumw2()

    num_trigs = trig_2d_dist.Integral()
    num_trigs_mc = trig_2d_dist_mc.Integral()

    axes = arr.array('i', [2, 3, 4, 5])
    h_lambda = h_lambda.Projection(4, axes)
    h_lambda_mc = h_lambda_mc.Projection(4, axes)
    h_lambda_mixed = h_lambda_mixed.Projection(4, axes)
    h_lambda_mixed_mc = h_lambda_mixed_mc.Projection(4, axes)

    h_h = h_h.Projection(2, 3, 4)
    h_h_mc = h_h_mc.Projection(2, 3, 4)
    h_h_mixed = h_h_mixed.Projection(2, 3, 4)
    h_h_mixed_mc = h_h_mixed_mc.Projection(2, 3, 4)

    h_lambda_2d_nomixcor = h_lambda.Projection(0, 1).Clone("h_lambda_2d_nomixcor")
    h_lambda_2d_nomixcor_mc = h_lambda_mc.Projection(0, 1).Clone("h_lambda_2d_nomixcor_mc")
    h_lambda_mixed_2d = h_lambda_mixed.Projection(0, 1).Clone("h_lambda_mixed_2d")
    h_lambda_mixed_2d_mc = h_lambda_mixed_mc.Projection(0, 1).Clone("h_lambda_mixed_2d_mc")

    h_h_2d_nomixcor = h_h.Project3D("xye").Clone("h_h_2d_nomixcor")
    h_h_2d_nomixcor_mc = h_h_mc.Project3D("xye").Clone("h_h_2d_nomixcor_mc")
    h_h_mixed_2d = h_h_mixed.Project3D("xye").Clone("h_h_mixed_2d")
    h_h_mixed_2d_mc = h_h_mixed_mc.Project3D("xye").Clone("h_h_mixed_2d_mc")


    h_lambda_2d_mixcor_sig = make_mixed_corrections(h_lambda, h_lambda_mixed, SIG_MIN, SIG_MAX)
    h_lambda_2d_mixcor_sig_mc = make_mixed_corrections(h_lambda_mc, h_lambda_mixed_mc, SIG_MIN, SIG_MAX)

    h_lambda_2d_mixcor_sig_mc_ratio = h_lambda_2d_mixcor_sig.Clone("h_lambda_2d_mixcor_sig_mc_ratio")
    h_lambda_2d_mixcor_sig_mc_ratio.Divide(h_lambda_2d_mixcor_sig_mc)

    if PT_MODE == 0:
        template_file = rt.TFile("twotrack_template.root", "RECREATE")
    elif PT_MODE == 1:
        template_file = rt.TFile("twotrack_template_lowpt.root", "RECREATE")
    elif PT_MODE == 2:
        template_file = rt.TFile("twotrack_template_highpt.root", "RECREATE")
    template_file.cd()
    h_lambda_2d_mixcor_sig_mc_ratio.Write("twotrack_template")
    template_file.Close()
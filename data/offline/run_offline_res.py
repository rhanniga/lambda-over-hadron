#! /usr/bin/env python
import math
import configparser

import array as arr
import ROOT as rt

from strangehelper import make_mixed_corrections

# epsilon used to avoid bin edge nightmares (if you pick a value that lies on bin edge, it defaults to right bin)
EPSILON = 0.00001

# read in the parameters from the config file
config = configparser.ConfigParser()
config.read("config/offline_res_config.ini")


# decide whether to do sideband subtraction or not
DO_SIDEBAND_SUBTRACTION = config.getboolean("GENERAL", "DO_SIDEBAND_SUBTRACTION")

# decide whether to do analysis with single highest pt trigger
DO_HIGHEST_PT = config.getboolean("GENERAL", "DO_HIGHEST_PT")

# UE line method
USE_AVG_4 = config.getboolean("GENERAL", "USE_AVG_4")
USE_AVG_6 = config.getboolean("GENERAL", "USE_AVG_6")
USE_ZYAM = config.getboolean("GENERAL", "USE_ZYAM")
assert sum([USE_AVG_4, USE_AVG_6, USE_ZYAM]) == 1, "Only select 1 method for UE line please"

# ETA CUTS 
ETA_MIN = config.getfloat("ETA_CUTS", "ETA_MIN")
ETA_MAX = config.getfloat("ETA_CUTS", "ETA_MAX") - EPSILON
DELTA_ETA_MIN = config.getfloat("ETA_CUTS", "DELTA_ETA_MIN")
DELTA_ETA_MAX = config.getfloat("ETA_CUTS", "DELTA_ETA_MAX") - EPSILON


# PT CUTS
TRIG_PT_LOW = config.getfloat("PT_CUTS", "TRIG_PT_LOW")
TRIG_PT_HIGH = config.getfloat("PT_CUTS", "TRIG_PT_HIGH") - EPSILON
ASSOC_PT_LOW = config.getfloat("PT_CUTS", "ASSOC_PT_LOW")
ASSOC_PT_HIGH = config.getfloat("PT_CUTS", "ASSOC_PT_HIGH") - EPSILON


# SIGNAL AND SIDEBAND REGION CUTS
SIG_MIN = config.getfloat("REGION_CUTS", "SIG_MIN")
SIG_MAX = config.getfloat("REGION_CUTS", "SIG_MAX") - EPSILON
RSB_MIN_0_20 = config.getfloat("REGION_CUTS", "RSB_MIN_0_20")
RSB_MAX_0_20 = config.getfloat("REGION_CUTS", "RSB_MAX_0_20") - EPSILON
RSB_MIN_20_50 = config.getfloat("REGION_CUTS", "RSB_MIN_20_50")
RSB_MAX_20_50 = config.getfloat("REGION_CUTS", "RSB_MAX_20_50") - EPSILON
RSB_MIN_50_80 = config.getfloat("REGION_CUTS", "RSB_MIN_50_80")
RSB_MAX_50_80 = config.getfloat("REGION_CUTS", "RSB_MAX_50_80") - EPSILON


# Output file containing all of the relevant results (long name )
output_file_string = "output/res_" 
output_file_string += ("highest_pt_" if DO_HIGHEST_PT else "") 
output_file_string += ("avg4_" if USE_AVG_4 else "") 
output_file_string += ("avg6_" if USE_AVG_6 else "") 
output_file_string += ("zyam_" if USE_ZYAM else "") 
output_file_string += ("sideband_subtraction_" if DO_SIDEBAND_SUBTRACTION else "")
output_file_string += "sig_" + str(SIG_MIN).replace(".", "") + "_" + str(SIG_MAX + EPSILON).replace(".", "") + "_"
output_file_string += "trig_" + str(TRIG_PT_LOW).replace(".", "") + "_" + str(TRIG_PT_HIGH + EPSILON).replace(".", "") + "_"
output_file_string += "assoc_" + str(ASSOC_PT_LOW).replace(".", "") + "_" + str(ASSOC_PT_HIGH + EPSILON).replace(".", "") + "_"
output_file_string += "delta_eta_" + str(DELTA_ETA_MAX + EPSILON).replace(".", "") + ".root"

PID_CORRECTION = (1.025/(1.71e-1)) # found from fitting dphi ratios of final dists with and without TOF veto

output_file = rt.TFile(output_file_string, "RECREATE")

############################################################################################################
############################################################################################################
############################################# 0-20 CENTRALITY ##############################################
############################################################################################################
############################################################################################################


input_file_0_20 = rt.TFile("~/OneDrive/Research/Output/lambda-over-hadron/data/res_cent_0_20.root")
input_list_0_20 = input_file_0_20.Get("h-lambda")
input_file_0_20.Close()


trig_dist_0_20 = input_list_0_20.FindObject("fTriggerDist_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fTriggerDist")
assoc_dist_0_20 = input_list_0_20.FindObject("fAssociatedHDist")
lambda_dist_0_20 = input_list_0_20.FindObject("fTriggeredLambdaDist")
lambda_ls_dist_0_20 = input_list_0_20.FindObject("fTriggeredLambdaLSDist")


h_h_0_20 = input_list_0_20.FindObject("fDphiHH_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHH")
h_h_mixed_0_20 = input_list_0_20.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHHMixed")

h_lambda_0_20 = input_list_0_20.FindObject("fDphiHLambda_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHLambda")
h_lambda_mixed_0_20 = input_list_0_20.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
assoc_dist_0_20.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
lambda_dist_0_20.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
lambda_ls_dist_0_20.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_ls_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


# Getting the single-particle trigger distributions
trig_pt_dist_0_20 = trig_dist_0_20.Projection(0).Clone("trig_pt_dist_0_20")
trig_phi_dist_0_20 = trig_dist_0_20.Projection(1).Clone("trig_phi_dist_0_20")
trig_eta_dist_0_20 = trig_dist_0_20.Projection(2).Clone("trig_eta_dist_0_20")
trig_2d_dist_0_20 = trig_dist_0_20.Projection(0, 3).Clone("trig_2d_dist_0_20")

# Get total number of triggers (used for per-trigger normalization)
num_trigs_0_20 = trig_2d_dist_0_20.Integral()

output_file.cd()
trig_pt_dist_0_20.Write()
trig_phi_dist_0_20.Write()
trig_eta_dist_0_20.Write()
trig_2d_dist_0_20.Write()

### ASSOCIATED HADRON SECTION ### 

assoc_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


# Getting the single-particle trigger distributions

assoc_pt_dist_0_20 = assoc_dist_0_20.Projection(0).Clone("assoc_pt_dist_0_20")
assoc_phi_dist_0_20 = assoc_dist_0_20.Projection(1).Clone("assoc_phi_dist_0_20")
assoc_eta_dist_0_20 = assoc_dist_0_20.Projection(2).Clone("assoc_eta_dist_0_20")
assoc_2d_dist_0_20 = assoc_dist_0_20.Projection(0, 3).Clone("assoc_2d_dist_0_20")

output_file.cd()
assoc_pt_dist_0_20.Write()
assoc_phi_dist_0_20.Write()
assoc_eta_dist_0_20.Write()
assoc_2d_dist_0_20.Write()

### SIGNAL ANALYSIS SECTION ###

lambda_mass_dist_0_20 = lambda_dist_0_20.Projection(3).Clone("lambda_mass_dist_0_20")
lambda_mass_ls_dist_0_20 = lambda_ls_dist_0_20.Projection(3).Clone("lambda_mass_ls_dist_0_20")


# scale LS to match lambda dist in RSB
left_rsb_bin_0_20 = lambda_mass_dist_0_20.FindBin(RSB_MIN_0_20)
right_rsb_bin_0_20 = lambda_mass_dist_0_20.FindBin(RSB_MAX_0_20)
lambda_mass_ls_dist_0_20.Scale(lambda_mass_dist_0_20.Integral(left_rsb_bin_0_20, right_rsb_bin_0_20)/lambda_mass_ls_dist_0_20.Integral(left_rsb_bin_0_20, right_rsb_bin_0_20))

residual_0_20 = lambda_mass_dist_0_20.Clone("residual_0_20")
residual_0_20.Add(lambda_mass_ls_dist_0_20, -1)

RSB_region_0_20 = rt.TBox(RSB_MIN_0_20, 0, RSB_MAX_0_20, lambda_mass_dist_0_20.GetMaximum()*1.055)
RSB_region_0_20.SetLineColor(rt.kRed)
RSB_region_0_20.SetFillColor(rt.kRed)
RSB_region_0_20.SetFillStyle(3003)


RSB_min_line_0_20 = rt.TLine(RSB_MIN_0_20, 0, RSB_MIN_0_20, lambda_mass_dist_0_20.GetMaximum()*1.05)
RSB_min_line_0_20.SetLineColor(rt.kRed)
RSB_min_line_0_20.SetLineWidth(2)
RSB_min_line_0_20.SetLineStyle(2)


RSB_max_line_0_20 = rt.TLine(RSB_MAX_0_20, 0, RSB_MAX_0_20, lambda_mass_dist_0_20.GetMaximum()*1.05)
RSB_max_line_0_20.SetLineColor(rt.kRed)
RSB_max_line_0_20.SetLineWidth(2)
RSB_max_line_0_20.SetLineStyle(2)


left_signal_bin_0_20 = lambda_mass_dist_0_20.FindBin(SIG_MIN)
right_signal_bin_0_20 = lambda_mass_dist_0_20.FindBin(SIG_MAX)

lambda_bg_0_20 = 0
lambda_total_0_20 = 0
for bin_num in range(left_signal_bin_0_20, right_signal_bin_0_20 + 1):
    lambda_total_0_20 += lambda_mass_dist_0_20.GetBinContent(bin_num)
    lambda_bg_0_20 += lambda_mass_ls_dist_0_20.GetBinContent(bin_num)

lambda_signal_0_20 = lambda_total_0_20 - lambda_bg_0_20
lambda_signal_total_ratio_0_20 = lambda_signal_0_20/lambda_total_0_20

scale_factor_0_20 = residual_0_20.Integral(1, left_rsb_bin_0_20)/residual_0_20.Integral(left_signal_bin_0_20, right_signal_bin_0_20)

output_file.cd()
RSB_region_0_20.Write("RSB_region_0_20")
RSB_min_line_0_20.Write("RSB_min_line_0_20")
RSB_max_line_0_20.Write("RSB_max_line_0_20")

lambda_mass_dist_0_20.Write()
lambda_mass_ls_dist_0_20.Write()
residual_0_20.Write()
RSB_region_0_20.Write()
RSB_min_line_0_20.Write()
RSB_max_line_0_20.Write()

### MIXED EVENT CORRECTION SECTION ###

axes = arr.array('i', [2, 3, 4, 5])
h_lambda_0_20 = h_lambda_0_20.Projection(4, axes)
h_lambda_mixed_0_20 = h_lambda_mixed_0_20.Projection(4, axes)

h_h_0_20 = h_h_0_20.Projection(2, 3, 4)
h_h_mixed_0_20 = h_h_mixed_0_20.Projection(2, 3, 4)


# Setting up 2-d correlation plots before the mixed event correction
h_lambda_2d_nomixcor_0_20 = h_lambda_0_20.Projection(0, 1).Clone("h_lambda_2d_nomixcor_0_20")
h_lambda_mixed_2d_0_20 = h_lambda_mixed_0_20.Projection(0, 1).Clone("h_lambda_mixed_2d_0_20")

h_h_2d_nomixcor_0_20 = h_h_0_20.Project3D("xye").Clone("h_h_2d_nomixcor_0_20")
h_h_mixed_2d_0_20 = h_h_mixed_0_20.Project3D("xye").Clone("h_h_mixed_2d_0_20")

h_lambda_2d_mixcor_sig_0_20 = make_mixed_corrections(h_lambda_0_20, h_lambda_mixed_0_20, SIG_MIN, SIG_MAX)
h_lambda_2d_mixcor_rsb_0_20 = make_mixed_corrections(h_lambda_0_20, h_lambda_mixed_0_20, RSB_MIN_0_20, RSB_MAX_0_20)

h_h_2d_mixcor_0_20 = make_mixed_corrections(h_h_0_20, h_h_mixed_0_20, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_0_20.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_0_20.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_h_2d_mixcor_0_20.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_0_20.Scale(1.0/num_trigs_0_20)
h_lambda_2d_mixcor_rsb_0_20.Scale(1.0/num_trigs_0_20)
h_h_2d_mixcor_0_20.Scale(1.0/num_trigs_0_20)

# scaling by total signal/signal region done here
h_lambda_2d_mixcor_sig_0_20.Scale(scale_factor_0_20)
h_lambda_2d_mixcor_rsb_0_20.Scale(scale_factor_0_20)
h_lambda_2d_mixcor_sig_0_20.Scale(PID_CORRECTION)
h_lambda_2d_mixcor_rsb_0_20.Scale(PID_CORRECTION)


output_file.cd()
h_lambda_2d_nomixcor_0_20.Write()
h_lambda_mixed_2d_0_20.Write()
h_h_2d_nomixcor_0_20.Write()
h_h_mixed_2d_0_20.Write()
h_lambda_2d_mixcor_sig_0_20.Write()
h_lambda_2d_mixcor_rsb_0_20.Write()
h_h_2d_mixcor_0_20.Write()

# SIDEBAND SUBTRACTION SECTION


# Normalize side band to 1
h_lambda_2d_mixcor_rsb_0_20.Scale(1/h_lambda_2d_mixcor_rsb_0_20.Integral())


# using RSB for sideband subtraction
h_lambda_2d_subtracted_0_20 = h_lambda_2d_mixcor_sig_0_20.Clone("h_lambda_2d_subtracted_0_20")
bg_integral_0_20 = (1 - lambda_signal_total_ratio_0_20)*h_lambda_2d_subtracted_0_20.Integral()
h_lambda_2d_subtracted_0_20.Add(h_lambda_2d_mixcor_rsb_0_20, -bg_integral_0_20)

h_lambda_dphi_rsb_0_20 = h_lambda_2d_mixcor_rsb_0_20.ProjectionY("h_lambda_dphi_rsb_0_20")
output_file.cd()
h_lambda_dphi_rsb_0_20.Write()


# INTEGRAL AND RATIO SECTION

h_lambda_dphi_subtracted_0_20 = h_lambda_2d_subtracted_0_20.ProjectionY("h_lambda_dphi_subtracted_0_20")

if USE_AVG_4:
    ue_line_0_20 = rt.TF1("ue_line_0_20", "pol0", -2, 6)
    ue_upper_line_0_20 = rt.TF1("ue_upper_line_0_20", "pol0", -2, 6)
    ue_lower_line_0_20 = rt.TF1("ue_lower_line_0_20", "pol0", -2, 6)
    zero_line_0_20 = rt.TF1("zero_line_0_20", "pol0", -2, 6)
    zero_upper_line_0_20 = rt.TF1("zero_upper_line_0_20", "pol0", -2, 6)
    zero_lower_line_0_20 = rt.TF1("zero_lower_line_0_20", "pol0", -2, 6)
    ue_avg_0_20 = (h_lambda_dphi_subtracted_0_20.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(16))/4

    ue_avg_error_0_20 = (1/4)*(math.sqrt(h_lambda_dphi_subtracted_0_20.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_0_20.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(16)**2))


    ue_line_0_20.SetParameter(0, ue_avg_0_20)
    ue_upper_line_0_20.SetParameter(0, ue_avg_0_20 + ue_avg_error_0_20)
    ue_lower_line_0_20.SetParameter(0, ue_avg_0_20 - ue_avg_error_0_20)

    zero_line_0_20.SetParameter(0, 0)
    zero_upper_line_0_20.SetParameter(0, ue_avg_error_0_20)
    zero_lower_line_0_20.SetParameter(0, -ue_avg_error_0_20)

elif USE_AVG_6:
    ue_line_0_20 = rt.TF1("ue_line_0_20", "pol0", -2, 6)
    ue_upper_line_0_20 = rt.TF1("ue_upper_line_0_20", "pol0", -2, 6)
    ue_lower_line_0_20 = rt.TF1("ue_lower_line_0_20", "pol0", -2, 6)
    zero_line_0_20 = rt.TF1("zero_line_0_20", "pol0", -2, 6)
    zero_upper_line_0_20 = rt.TF1("zero_upper_line_0_20", "pol0", -2, 6)
    zero_lower_line_0_20 = rt.TF1("zero_lower_line_0_20", "pol0", -2, 6)
    ue_avg_0_20 = (h_lambda_dphi_subtracted_0_20.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(16))/6

    ue_avg_error_0_20 = (1/6)*(math.sqrt(h_lambda_dphi_subtracted_0_20.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_0_20.GetBinError(2)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(7)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_0_20.GetBinError(16)**2))


    ue_line_0_20.SetParameter(0, ue_avg_0_20)
    ue_upper_line_0_20.SetParameter(0, ue_avg_0_20 + ue_avg_error_0_20)
    ue_lower_line_0_20.SetParameter(0, ue_avg_0_20 - ue_avg_error_0_20)

    zero_line_0_20.SetParameter(0, 0)
    zero_upper_line_0_20.SetParameter(0, ue_avg_error_0_20)
    zero_lower_line_0_20.SetParameter(0, -ue_avg_error_0_20)

elif USE_ZYAM:
    ue_line_0_20 = rt.TF1("ue_line_0_20", "pol0", -2, 6)
    ue_upper_line_0_20 = rt.TF1("ue_upper_line_0_20", "pol0", -2, 6)
    ue_lower_line_0_20 = rt.TF1("ue_lower_line_0_20", "pol0", -2, 6)
    zero_line_0_20 = rt.TF1("zero_line_0_20", "pol0", -2, 6)
    zero_upper_line_0_20 = rt.TF1("zero_upper_line_0_20", "pol0", -2, 6)
    zero_lower_line_0_20 = rt.TF1("zero_lower_line_0_20", "pol0", -2, 6)
    min_bin = h_lambda_dphi_subtracted_0_20.GetMinimumBin()
    ue_avg_0_20 = h_lambda_dphi_subtracted_0_20.GetBinContent(min_bin)
    ue_avg_error_0_20 = h_lambda_dphi_subtracted_0_20.GetBinError(min_bin)


    ue_line_0_20.SetParameter(0, ue_avg_0_20)
    ue_upper_line_0_20.SetParameter(0, ue_avg_0_20 + ue_avg_error_0_20)
    ue_lower_line_0_20.SetParameter(0, ue_avg_0_20 - ue_avg_error_0_20)

    zero_line_0_20.SetParameter(0, 0)
    zero_upper_line_0_20.SetParameter(0, ue_avg_error_0_20)
    zero_lower_line_0_20.SetParameter(0, -ue_avg_error_0_20)
else:
    raise NotImplementedError("UE line mode not supported")


h_lambda_dphi_subtracted_0_20_zeroed = h_lambda_dphi_subtracted_0_20.Clone("h_lambda_dphi_subtracted_0_20_zeroed")
h_lambda_dphi_subtracted_0_20_zeroed.Add(ue_line_0_20, -1)


DPHI_BINS = h_lambda_dphi_subtracted_0_20.GetNbinsX()

h_lambda_total_integral_0_20 = 0
h_lambda_near_integral_0_20 = 0
h_lambda_away_integral_0_20 = 0
h_lambda_ue_integral_0_20 = ue_avg_0_20*DPHI_BINS

h_lambda_total_integral_error_0_20 = 0
h_lambda_near_integral_error_0_20 = 0
h_lambda_away_integral_error_0_20 = 0
h_lambda_ue_integral_error_0_20 = ue_avg_error_0_20*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_lambda_total_integral_0_20 += h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num)
    h_lambda_total_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2
    part = h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num) - ue_avg_0_20
    if part < 0:
        continue
    if bin_num < 9:
        h_lambda_near_integral_0_20 += part
        h_lambda_near_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2
        h_lambda_near_integral_error_0_20 += ue_avg_error_0_20**2
    else:
        h_lambda_away_integral_0_20 += part
        h_lambda_away_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2
        h_lambda_away_integral_error_0_20 += ue_avg_error_0_20**2

h_lambda_total_integral_error_0_20 = math.sqrt(h_lambda_total_integral_error_0_20)
h_lambda_near_integral_error_0_20 = math.sqrt(h_lambda_near_integral_error_0_20)
h_lambda_away_integral_error_0_20 = math.sqrt(h_lambda_away_integral_error_0_20)



h_h_dphi_0_20 = h_h_2d_mixcor_0_20.ProjectionY("h_h_dphi_0_20")

if USE_AVG_4:
    hh_ue_line_0_20 = rt.TF1("hh_ue_line_0_20", "pol0", -2, 6)
    hh_ue_upper_line_0_20 = rt.TF1("hh_ue_upper_line_0_20", "pol0", -2, 6)
    hh_ue_lower_line_0_20 = rt.TF1("hh_ue_lower_line_0_20", "pol0", -2, 6)
    hh_zero_line_0_20 = rt.TF1("hh_zero_line_0_20", "pol0", -2, 6)
    hh_zero_upper_line_0_20 = rt.TF1("hh_zero_upper_line_0_20", "pol0", -2, 6)
    hh_zero_lower_line_0_20 = rt.TF1("hh_zero_lower_line_0_20", "pol0", -2, 6)
    hh_ue_avg_0_20 = (h_h_dphi_0_20.GetBinContent(1) 
                   + h_h_dphi_0_20.GetBinContent(8)
                   + h_h_dphi_0_20.GetBinContent(9)
                   + h_h_dphi_0_20.GetBinContent(16))/4

    hh_ue_avg_error_0_20 = (1/4)*(math.sqrt(h_h_dphi_0_20.GetBinError(1)**2 
                   + h_h_dphi_0_20.GetBinError(8)**2
                   + h_h_dphi_0_20.GetBinError(9)**2
                   + h_h_dphi_0_20.GetBinError(16)**2))


    hh_ue_line_0_20.SetParameter(0, hh_ue_avg_0_20)
    hh_ue_upper_line_0_20.SetParameter(0, hh_ue_avg_0_20 + hh_ue_avg_error_0_20)
    hh_ue_lower_line_0_20.SetParameter(0, hh_ue_avg_0_20 - hh_ue_avg_error_0_20)

    hh_zero_line_0_20.SetParameter(0, 0)
    hh_zero_upper_line_0_20.SetParameter(0, hh_ue_avg_error_0_20)
    hh_zero_lower_line_0_20.SetParameter(0, -hh_ue_avg_error_0_20)

elif USE_AVG_6:
    hh_ue_line_0_20 = rt.TF1("hh_ue_line_0_20", "pol0", -2, 6)
    hh_ue_upper_line_0_20 = rt.TF1("hh_ue_upper_line_0_20", "pol0", -2, 6)
    hh_ue_lower_line_0_20 = rt.TF1("hh_ue_lower_line_0_20", "pol0", -2, 6)
    hh_zero_line_0_20 = rt.TF1("hh_zero_line_0_20", "pol0", -2, 6)
    hh_zero_upper_line_0_20 = rt.TF1("hh_zero_upper_line_0_20", "pol0", -2, 6)
    hh_zero_lower_line_0_20 = rt.TF1("hh_zero_lower_line_0_20", "pol0", -2, 6)
    hh_ue_avg_0_20 = (h_h_dphi_0_20.GetBinContent(1) 
                   + h_h_dphi_0_20.GetBinContent(2)
                   + h_h_dphi_0_20.GetBinContent(7)
                   + h_h_dphi_0_20.GetBinContent(8)
                   + h_h_dphi_0_20.GetBinContent(9)
                   + h_h_dphi_0_20.GetBinContent(16))/6

    hh_ue_avg_error_0_20 = (1/6)*(math.sqrt(h_h_dphi_0_20.GetBinError(1)**2 
                   + h_h_dphi_0_20.GetBinError(2)**2
                   + h_h_dphi_0_20.GetBinError(7)**2
                   + h_h_dphi_0_20.GetBinError(8)**2
                   + h_h_dphi_0_20.GetBinError(9)**2
                   + h_h_dphi_0_20.GetBinError(16)**2))


    hh_ue_line_0_20.SetParameter(0, hh_ue_avg_0_20)
    hh_ue_upper_line_0_20.SetParameter(0, hh_ue_avg_0_20 + hh_ue_avg_error_0_20)
    hh_ue_lower_line_0_20.SetParameter(0, hh_ue_avg_0_20 - hh_ue_avg_error_0_20)

    hh_zero_line_0_20.SetParameter(0, 0)
    hh_zero_upper_line_0_20.SetParameter(0, hh_ue_avg_error_0_20)
    hh_zero_lower_line_0_20.SetParameter(0, -hh_ue_avg_error_0_20)
elif USE_ZYAM:
    hh_ue_line_0_20 = rt.TF1("hh_ue_line_0_20", "pol0", -2, 6)
    hh_ue_upper_line_0_20 = rt.TF1("hh_ue_upper_line_0_20", "pol0", -2, 6)
    hh_ue_lower_line_0_20 = rt.TF1("hh_ue_lower_line_0_20", "pol0", -2, 6)
    hh_zero_line_0_20 = rt.TF1("hh_zero_line_0_20", "pol0", -2, 6)
    hh_zero_upper_line_0_20 = rt.TF1("hh_zero_upper_line_0_20", "pol0", -2, 6)
    hh_zero_lower_line_0_20 = rt.TF1("hh_zero_lower_line_0_20", "pol0", -2, 6)
    
    min_bin = h_h_dphi_0_20.GetMinimumBin()
    hh_ue_avg_0_20 = h_h_dphi_0_20.GetBinContent(min_bin)
    hh_ue_avg_error_0_20 = h_h_dphi_0_20.GetBinError(min_bin)

    hh_ue_line_0_20.SetParameter(0, hh_ue_avg_0_20)
    hh_ue_upper_line_0_20.SetParameter(0, hh_ue_avg_0_20 + hh_ue_avg_error_0_20)
    hh_ue_lower_line_0_20.SetParameter(0, hh_ue_avg_0_20 - hh_ue_avg_error_0_20)

    hh_zero_line_0_20.SetParameter(0, 0)
    hh_zero_upper_line_0_20.SetParameter(0, hh_ue_avg_error_0_20)
    hh_zero_lower_line_0_20.SetParameter(0, -hh_ue_avg_error_0_20)
else:
    raise NotImplementedError("UE line mode not supported")


h_h_dphi_0_20_zeroed = h_h_dphi_0_20.Clone("h_h_dphi_0_20_zeroed")
h_h_dphi_0_20_zeroed.Add(hh_ue_line_0_20, -1)

h_h_total_integral_0_20 = 0
h_h_near_integral_0_20 = 0
h_h_away_integral_0_20 = 0
h_h_ue_integral_0_20 = hh_ue_avg_0_20*DPHI_BINS

h_h_total_integral_error_0_20 = 0
h_h_near_integral_error_0_20 = 0
h_h_away_integral_error_0_20 = 0
h_h_ue_integral_error_0_20 = hh_ue_avg_error_0_20*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_h_total_integral_0_20 += h_h_dphi_0_20.GetBinContent(bin_num)
    h_h_total_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
    part = h_h_dphi_0_20.GetBinContent(bin_num) - hh_ue_avg_0_20
    if part < 0:
        continue
    if bin_num < 9:
        h_h_near_integral_0_20 += part
        h_h_near_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
        h_h_near_integral_error_0_20 += hh_ue_avg_error_0_20**2
    else:
        h_h_away_integral_0_20 += part
        h_h_away_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
        h_h_away_integral_error_0_20 += hh_ue_avg_error_0_20**2
h_h_total_integral_error_0_20 = math.sqrt(h_h_total_integral_error_0_20)
h_h_near_integral_error_0_20 = math.sqrt(h_h_near_integral_error_0_20)
h_h_away_integral_error_0_20 = math.sqrt(h_h_away_integral_error_0_20)

output_file.cd()
h_lambda_2d_subtracted_0_20.Write()

h_lambda_dphi_subtracted_0_20.Write()
h_lambda_dphi_subtracted_0_20_zeroed.Write()
h_h_dphi_0_20.Write()
h_h_dphi_0_20_zeroed.Write()

ue_upper_line_0_20.Write()
ue_line_0_20.Write()
ue_lower_line_0_20.Write()
zero_upper_line_0_20.Write()
zero_line_0_20.Write()
zero_lower_line_0_20.Write()

hh_ue_upper_line_0_20.Write()
hh_ue_line_0_20.Write()
hh_ue_lower_line_0_20.Write()
hh_zero_upper_line_0_20.Write()
hh_zero_line_0_20.Write()
hh_zero_lower_line_0_20.Write()

near_ratio_0_20 = h_lambda_near_integral_0_20/h_h_near_integral_0_20
away_ratio_0_20 = h_lambda_away_integral_0_20/h_h_away_integral_0_20
ue_ratio_0_20 = h_lambda_ue_integral_0_20/h_h_ue_integral_0_20
total_ratio_0_20 = h_lambda_total_integral_0_20/h_h_total_integral_0_20

near_ratio_error_0_20 = near_ratio_0_20*math.sqrt((h_lambda_near_integral_error_0_20/h_lambda_near_integral_0_20)**2
                                                 + (h_h_near_integral_error_0_20/h_h_near_integral_0_20)**2)
away_ratio_error_0_20 = away_ratio_0_20*math.sqrt((h_lambda_away_integral_error_0_20/h_lambda_away_integral_0_20)**2
                                                 + (h_h_away_integral_error_0_20/h_h_away_integral_0_20)**2)
ue_ratio_error_0_20 = ue_ratio_0_20*math.sqrt((h_lambda_ue_integral_error_0_20/h_lambda_ue_integral_0_20)**2
                                                 + (h_h_ue_integral_error_0_20/h_h_ue_integral_0_20)**2)
total_ratio_error_0_20 = total_ratio_0_20*math.sqrt((h_lambda_total_integral_error_0_20/h_lambda_total_integral_0_20)**2
                                                 + (h_h_total_integral_error_0_20/h_h_total_integral_0_20)**2)



############################################################################################################
############################################################################################################
############################################ 20-50 CENTRALITY ##############################################
############################################################################################################
############################################################################################################

input_file_20_50 = rt.TFile("~/OneDrive/Research/Output/lambda-over-hadron/data/res_cent_20_50.root")
input_list_20_50 = input_file_20_50.Get("h-lambda")
input_file_20_50.Close()


trig_dist_20_50 = input_list_20_50.FindObject("fTriggerDist_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fTriggerDist")
lambda_dist_20_50 = input_list_20_50.FindObject("fTriggeredLambdaDist")
lambda_ls_dist_20_50 = input_list_20_50.FindObject("fTriggeredLambdaLSDist")


h_h_20_50 = input_list_20_50.FindObject("fDphiHH_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHH")
h_h_mixed_20_50 = input_list_20_50.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHHMixed")

h_lambda_20_50 = input_list_20_50.FindObject("fDphiHLambda_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHLambda")
h_lambda_mixed_20_50 = input_list_20_50.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_20_50.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
lambda_ls_dist_20_50.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_20_50.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_20_50.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_ls_dist_20_50.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


# Getting the single-particle trigger distributions
trig_pt_dist_20_50 = trig_dist_20_50.Projection(0).Clone("trig_pt_dist_20_50")
trig_phi_dist_20_50 = trig_dist_20_50.Projection(1).Clone("trig_phi_dist_20_50")
trig_eta_dist_20_50 = trig_dist_20_50.Projection(2).Clone("trig_eta_dist_20_50")
trig_2d_dist_20_50 = trig_dist_20_50.Projection(0, 3).Clone("trig_2d_dist_20_50")

# Get total number of triggers (used for per-trigger normalization)
num_trigs_20_50 = trig_2d_dist_20_50.Integral()

output_file.cd()
trig_pt_dist_20_50.Write()
trig_phi_dist_20_50.Write()
trig_eta_dist_20_50.Write()
trig_2d_dist_20_50.Write()

### SIGNAL ANALYSIS SECTION ###

lambda_mass_dist_20_50 = lambda_dist_20_50.Projection(3).Clone("lambda_mass_dist_20_50")
lambda_mass_ls_dist_20_50 = lambda_ls_dist_20_50.Projection(3).Clone("lambda_mass_ls_dist_20_50")


# scale LS to match lambda dist in RSB
left_rsb_bin_20_50 = lambda_mass_dist_20_50.FindBin(RSB_MIN_20_50)
right_rsb_bin_20_50 = lambda_mass_dist_20_50.FindBin(RSB_MAX_20_50)
lambda_mass_ls_dist_20_50.Scale(lambda_mass_dist_20_50.Integral(left_rsb_bin_20_50, right_rsb_bin_20_50)/lambda_mass_ls_dist_20_50.Integral(left_rsb_bin_20_50, right_rsb_bin_20_50))

residual_20_50 = lambda_mass_dist_20_50.Clone("residual_20_50")
residual_20_50.Add(lambda_mass_ls_dist_20_50, -1)

RSB_region_20_50 = rt.TBox(RSB_MIN_20_50, 0, RSB_MAX_20_50, lambda_mass_dist_20_50.GetMaximum()*1.055)
RSB_region_20_50.SetLineColor(rt.kRed)
RSB_region_20_50.SetFillColor(rt.kRed)
RSB_region_20_50.SetFillStyle(3003)


RSB_min_line_20_50 = rt.TLine(RSB_MIN_20_50, 0, RSB_MIN_20_50, lambda_mass_dist_20_50.GetMaximum()*1.05)
RSB_min_line_20_50.SetLineColor(rt.kRed)
RSB_min_line_20_50.SetLineWidth(2)
RSB_min_line_20_50.SetLineStyle(2)


RSB_max_line_20_50 = rt.TLine(RSB_MAX_20_50, 0, RSB_MAX_20_50, lambda_mass_dist_20_50.GetMaximum()*1.05)
RSB_max_line_20_50.SetLineColor(rt.kRed)
RSB_max_line_20_50.SetLineWidth(2)
RSB_max_line_20_50.SetLineStyle(2)


left_signal_bin_20_50 = lambda_mass_dist_20_50.FindBin(SIG_MIN)
right_signal_bin_20_50 = lambda_mass_dist_20_50.FindBin(SIG_MAX)

lambda_bg_20_50 = 0
lambda_total_20_50 = 0
for bin_num in range(left_signal_bin_20_50, right_signal_bin_20_50 + 1):
    lambda_total_20_50 += lambda_mass_dist_20_50.GetBinContent(bin_num)
    lambda_bg_20_50 += lambda_mass_ls_dist_20_50.GetBinContent(bin_num)

lambda_signal_20_50 = lambda_total_20_50 - lambda_bg_20_50
lambda_signal_total_ratio_20_50 = lambda_signal_20_50/lambda_total_20_50

scale_factor_20_50 = residual_20_50.Integral(1, left_rsb_bin_20_50)/residual_20_50.Integral(left_signal_bin_20_50, right_signal_bin_20_50)

output_file.cd()
RSB_region_20_50.Write("RSB_region_20_50")
RSB_min_line_20_50.Write("RSB_min_line_20_50")
RSB_max_line_20_50.Write("RSB_max_line_20_50")
lambda_mass_dist_20_50.Write()
lambda_mass_ls_dist_20_50.Write()
residual_20_50.Write()
RSB_region_20_50.Write()
RSB_min_line_20_50.Write()
RSB_max_line_20_50.Write()

### MIXED EVENT CORRECTION SECTION ###

axes = arr.array('i', [2, 3, 4, 5])
h_lambda_20_50 = h_lambda_20_50.Projection(4, axes)
h_lambda_mixed_20_50 = h_lambda_mixed_20_50.Projection(4, axes)

h_h_20_50 = h_h_20_50.Projection(2, 3, 4)
h_h_mixed_20_50 = h_h_mixed_20_50.Projection(2, 3, 4)


# Setting up 2-d correlation plots before the mixed event correction
h_lambda_2d_nomixcor_20_50 = h_lambda_20_50.Projection(0, 1).Clone("h_lambda_2d_nomixcor_20_50")
h_lambda_mixed_2d_20_50 = h_lambda_mixed_20_50.Projection(0, 1).Clone("h_lambda_mixed_2d_20_50")

h_h_2d_nomixcor_20_50 = h_h_20_50.Project3D("xye").Clone("h_h_2d_nomixcor_20_50")
h_h_mixed_2d_20_50 = h_h_mixed_20_50.Project3D("xye").Clone("h_h_mixed_2d_20_50")

h_lambda_2d_mixcor_sig_20_50 = make_mixed_corrections(h_lambda_20_50, h_lambda_mixed_20_50, SIG_MIN, SIG_MAX)
h_lambda_2d_mixcor_rsb_20_50 = make_mixed_corrections(h_lambda_20_50, h_lambda_mixed_20_50, RSB_MIN_20_50, RSB_MAX_20_50)

h_h_2d_mixcor_20_50 = make_mixed_corrections(h_h_20_50, h_h_mixed_20_50, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_20_50.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_20_50.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_h_2d_mixcor_20_50.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_20_50.Scale(1.0/num_trigs_20_50)
h_lambda_2d_mixcor_rsb_20_50.Scale(1.0/num_trigs_20_50)
h_h_2d_mixcor_20_50.Scale(1.0/num_trigs_20_50)

# scaling by total signal/signal region done here
h_lambda_2d_mixcor_sig_20_50.Scale(scale_factor_20_50)
h_lambda_2d_mixcor_rsb_20_50.Scale(scale_factor_20_50)
h_lambda_2d_mixcor_sig_20_50.Scale(PID_CORRECTION)
h_lambda_2d_mixcor_rsb_20_50.Scale(PID_CORRECTION)

output_file.cd()
h_lambda_2d_nomixcor_20_50.Write()
h_lambda_mixed_2d_20_50.Write()
h_h_2d_nomixcor_20_50.Write()
h_h_mixed_2d_20_50.Write()
h_lambda_2d_mixcor_sig_20_50.Write()
h_lambda_2d_mixcor_rsb_20_50.Write()
h_h_2d_mixcor_20_50.Write()

# SIDEBAND SUBTRACTION SECTION

# Normalize side band to 1
h_lambda_2d_mixcor_rsb_20_50.Scale(1/h_lambda_2d_mixcor_rsb_20_50.Integral())


# using RSB for sideband subtraction
h_lambda_2d_subtracted_20_50 = h_lambda_2d_mixcor_sig_20_50.Clone("h_lambda_2d_subtracted_20_50")
bg_integral_20_50 = (1 - lambda_signal_total_ratio_20_50)*h_lambda_2d_subtracted_20_50.Integral()
h_lambda_2d_subtracted_20_50.Add(h_lambda_2d_mixcor_rsb_20_50, -bg_integral_20_50)


# INTEGRAL AND RATIO SECTION

h_lambda_dphi_subtracted_20_50 = h_lambda_2d_subtracted_20_50.ProjectionY("h_lambda_dphi_subtracted_20_50")

if USE_AVG_4:
    ue_line_20_50 = rt.TF1("ue_line_20_50", "pol0", -2, 6)
    ue_upper_line_20_50 = rt.TF1("ue_upper_line_20_50", "pol0", -2, 6)
    ue_lower_line_20_50 = rt.TF1("ue_lower_line_20_50", "pol0", -2, 6)
    zero_line_20_50 = rt.TF1("zero_line_20_50", "pol0", -2, 6)
    zero_upper_line_20_50 = rt.TF1("zero_upper_line_20_50", "pol0", -2, 6)
    zero_lower_line_20_50 = rt.TF1("zero_lower_line_20_50", "pol0", -2, 6)
    ue_avg_20_50 = (h_lambda_dphi_subtracted_20_50.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(8)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(9)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(16))/4

    ue_avg_error_20_50 = (1/4)*(math.sqrt(h_lambda_dphi_subtracted_20_50.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_20_50.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(16)**2))


    ue_line_20_50.SetParameter(0, ue_avg_20_50)
    ue_upper_line_20_50.SetParameter(0, ue_avg_20_50 + ue_avg_error_20_50)
    ue_lower_line_20_50.SetParameter(0, ue_avg_20_50 - ue_avg_error_20_50)

    zero_line_20_50.SetParameter(0, 0)
    zero_upper_line_20_50.SetParameter(0, ue_avg_error_20_50)
    zero_lower_line_20_50.SetParameter(0, -ue_avg_error_20_50)

elif USE_AVG_6:
    ue_line_20_50 = rt.TF1("ue_line_20_50", "pol0", -2, 6)
    ue_upper_line_20_50 = rt.TF1("ue_upper_line_20_50", "pol0", -2, 6)
    ue_lower_line_20_50 = rt.TF1("ue_lower_line_20_50", "pol0", -2, 6)
    zero_line_20_50 = rt.TF1("zero_line_20_50", "pol0", -2, 6)
    zero_upper_line_20_50 = rt.TF1("zero_upper_line_20_50", "pol0", -2, 6)
    zero_lower_line_20_50 = rt.TF1("zero_lower_line_20_50", "pol0", -2, 6)
    ue_avg_20_50 = (h_lambda_dphi_subtracted_20_50.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(2)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(7)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(8)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(9)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(16))/6

    ue_avg_error_20_50 = (1/6)*(math.sqrt(h_lambda_dphi_subtracted_20_50.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_20_50.GetBinError(2)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(7)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_20_50.GetBinError(16)**2))


    ue_line_20_50.SetParameter(0, ue_avg_20_50)
    ue_upper_line_20_50.SetParameter(0, ue_avg_20_50 + ue_avg_error_20_50)
    ue_lower_line_20_50.SetParameter(0, ue_avg_20_50 - ue_avg_error_20_50)

    zero_line_20_50.SetParameter(0, 0)
    zero_upper_line_20_50.SetParameter(0, ue_avg_error_20_50)
    zero_lower_line_20_50.SetParameter(0, -ue_avg_error_20_50)

elif USE_ZYAM:
    ue_line_20_50 = rt.TF1("ue_line_20_50", "pol0", -2, 6)
    ue_upper_line_20_50 = rt.TF1("ue_upper_line_20_50", "pol0", -2, 6)
    ue_lower_line_20_50 = rt.TF1("ue_lower_line_20_50", "pol0", -2, 6)
    zero_line_20_50 = rt.TF1("zero_line_20_50", "pol0", -2, 6)
    zero_upper_line_20_50 = rt.TF1("zero_upper_line_20_50", "pol0", -2, 6)
    zero_lower_line_20_50 = rt.TF1("zero_lower_line_20_50", "pol0", -2, 6)
    min_bin = h_lambda_dphi_subtracted_20_50.GetMinimumBin()
    ue_avg_20_50 = h_lambda_dphi_subtracted_20_50.GetBinContent(min_bin)
    ue_avg_error_20_50 = h_lambda_dphi_subtracted_20_50.GetBinError(min_bin)


    ue_line_20_50.SetParameter(0, ue_avg_20_50)
    ue_upper_line_20_50.SetParameter(0, ue_avg_20_50 + ue_avg_error_20_50)
    ue_lower_line_20_50.SetParameter(0, ue_avg_20_50 - ue_avg_error_20_50)

    zero_line_20_50.SetParameter(0, 0)
    zero_upper_line_20_50.SetParameter(0, ue_avg_error_20_50)
    zero_lower_line_20_50.SetParameter(0, -ue_avg_error_20_50)
else:
    raise NotImplementedError("UE line mode not supported")


h_lambda_dphi_subtracted_20_50_zeroed = h_lambda_dphi_subtracted_20_50.Clone("h_lambda_dphi_subtracted_20_50_zeroed")
h_lambda_dphi_subtracted_20_50_zeroed.Add(ue_line_20_50, -1)


DPHI_BINS = h_lambda_dphi_subtracted_20_50.GetNbinsX()

h_lambda_total_integral_20_50 = 0
h_lambda_near_integral_20_50 = 0
h_lambda_away_integral_20_50 = 0
h_lambda_ue_integral_20_50 = ue_avg_20_50*DPHI_BINS

h_lambda_total_integral_error_20_50 = 0
h_lambda_near_integral_error_20_50 = 0
h_lambda_away_integral_error_20_50 = 0
h_lambda_ue_integral_error_20_50 = ue_avg_error_20_50*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_lambda_total_integral_20_50 += h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num)
    h_lambda_total_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2
    part = h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num) - ue_avg_20_50
    if part < 0:
        continue
    if bin_num < 9:
        h_lambda_near_integral_20_50 += part
        h_lambda_near_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2
        h_lambda_near_integral_error_20_50 += ue_avg_error_20_50**2
    else:
        h_lambda_away_integral_20_50 += part
        h_lambda_away_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2
        h_lambda_away_integral_error_20_50 += ue_avg_error_20_50**2

h_lambda_total_integral_error_20_50 = math.sqrt(h_lambda_total_integral_error_20_50)
h_lambda_near_integral_error_20_50 = math.sqrt(h_lambda_near_integral_error_20_50)
h_lambda_away_integral_error_20_50 = math.sqrt(h_lambda_away_integral_error_20_50)



h_h_dphi_20_50 = h_h_2d_mixcor_20_50.ProjectionY("h_h_dphi_20_50")

if USE_AVG_4:
    hh_ue_line_20_50 = rt.TF1("hh_ue_line_20_50", "pol0", -2, 6)
    hh_ue_upper_line_20_50 = rt.TF1("hh_ue_upper_line_20_50", "pol0", -2, 6)
    hh_ue_lower_line_20_50 = rt.TF1("hh_ue_lower_line_20_50", "pol0", -2, 6)
    hh_zero_line_20_50 = rt.TF1("hh_zero_line_20_50", "pol0", -2, 6)
    hh_zero_upper_line_20_50 = rt.TF1("hh_zero_upper_line_20_50", "pol0", -2, 6)
    hh_zero_lower_line_20_50 = rt.TF1("hh_zero_lower_line_20_50", "pol0", -2, 6)
    hh_ue_avg_20_50 = (h_h_dphi_20_50.GetBinContent(1) 
                   + h_h_dphi_20_50.GetBinContent(8)
                   + h_h_dphi_20_50.GetBinContent(9)
                   + h_h_dphi_20_50.GetBinContent(16))/4

    hh_ue_avg_error_20_50 = (1/4)*(math.sqrt(h_h_dphi_20_50.GetBinError(1)**2 
                   + h_h_dphi_20_50.GetBinError(8)**2
                   + h_h_dphi_20_50.GetBinError(9)**2
                   + h_h_dphi_20_50.GetBinError(16)**2))


    hh_ue_line_20_50.SetParameter(0, hh_ue_avg_20_50)
    hh_ue_upper_line_20_50.SetParameter(0, hh_ue_avg_20_50 + hh_ue_avg_error_20_50)
    hh_ue_lower_line_20_50.SetParameter(0, hh_ue_avg_20_50 - hh_ue_avg_error_20_50)

    hh_zero_line_20_50.SetParameter(0, 0)
    hh_zero_upper_line_20_50.SetParameter(0, hh_ue_avg_error_20_50)
    hh_zero_lower_line_20_50.SetParameter(0, -hh_ue_avg_error_20_50)

elif USE_AVG_6:
    hh_ue_line_20_50 = rt.TF1("hh_ue_line_20_50", "pol0", -2, 6)
    hh_ue_upper_line_20_50 = rt.TF1("hh_ue_upper_line_20_50", "pol0", -2, 6)
    hh_ue_lower_line_20_50 = rt.TF1("hh_ue_lower_line_20_50", "pol0", -2, 6)
    hh_zero_line_20_50 = rt.TF1("hh_zero_line_20_50", "pol0", -2, 6)
    hh_zero_upper_line_20_50 = rt.TF1("hh_zero_upper_line_20_50", "pol0", -2, 6)
    hh_zero_lower_line_20_50 = rt.TF1("hh_zero_lower_line_20_50", "pol0", -2, 6)
    hh_ue_avg_20_50 = (h_h_dphi_20_50.GetBinContent(1) 
                   + h_h_dphi_20_50.GetBinContent(2)
                   + h_h_dphi_20_50.GetBinContent(7)
                   + h_h_dphi_20_50.GetBinContent(8)
                   + h_h_dphi_20_50.GetBinContent(9)
                   + h_h_dphi_20_50.GetBinContent(16))/6

    hh_ue_avg_error_20_50 = (1/6)*(math.sqrt(h_h_dphi_20_50.GetBinError(1)**2 
                   + h_h_dphi_20_50.GetBinError(2)**2
                   + h_h_dphi_20_50.GetBinError(7)**2
                   + h_h_dphi_20_50.GetBinError(8)**2
                   + h_h_dphi_20_50.GetBinError(9)**2
                   + h_h_dphi_20_50.GetBinError(16)**2))


    hh_ue_line_20_50.SetParameter(0, hh_ue_avg_20_50)
    hh_ue_upper_line_20_50.SetParameter(0, hh_ue_avg_20_50 + hh_ue_avg_error_20_50)
    hh_ue_lower_line_20_50.SetParameter(0, hh_ue_avg_20_50 - hh_ue_avg_error_20_50)

    hh_zero_line_20_50.SetParameter(0, 0)
    hh_zero_upper_line_20_50.SetParameter(0, hh_ue_avg_error_20_50)
    hh_zero_lower_line_20_50.SetParameter(0, -hh_ue_avg_error_20_50)
elif USE_ZYAM:
    hh_ue_line_20_50 = rt.TF1("hh_ue_line_20_50", "pol0", -2, 6)
    hh_ue_upper_line_20_50 = rt.TF1("hh_ue_upper_line_20_50", "pol0", -2, 6)
    hh_ue_lower_line_20_50 = rt.TF1("hh_ue_lower_line_20_50", "pol0", -2, 6)
    hh_zero_line_20_50 = rt.TF1("hh_zero_line_20_50", "pol0", -2, 6)
    hh_zero_upper_line_20_50 = rt.TF1("hh_zero_upper_line_20_50", "pol0", -2, 6)
    hh_zero_lower_line_20_50 = rt.TF1("hh_zero_lower_line_20_50", "pol0", -2, 6)
    
    min_bin = h_h_dphi_20_50.GetMinimumBin()
    hh_ue_avg_20_50 = h_h_dphi_20_50.GetBinContent(min_bin)
    hh_ue_avg_error_20_50 = h_h_dphi_20_50.GetBinError(min_bin)

    hh_ue_line_20_50.SetParameter(0, hh_ue_avg_20_50)
    hh_ue_upper_line_20_50.SetParameter(0, hh_ue_avg_20_50 + hh_ue_avg_error_20_50)
    hh_ue_lower_line_20_50.SetParameter(0, hh_ue_avg_20_50 - hh_ue_avg_error_20_50)

    hh_zero_line_20_50.SetParameter(0, 0)
    hh_zero_upper_line_20_50.SetParameter(0, hh_ue_avg_error_20_50)
    hh_zero_lower_line_20_50.SetParameter(0, -hh_ue_avg_error_20_50)
else:
    raise NotImplementedError("UE line mode not supported")

h_h_dphi_20_50_zeroed = h_h_dphi_20_50.Clone("h_h_dphi_20_50_zeroed")
h_h_dphi_20_50_zeroed.Add(hh_ue_line_20_50, -1)

h_h_total_integral_20_50 = 0
h_h_near_integral_20_50 = 0
h_h_away_integral_20_50 = 0
h_h_ue_integral_20_50 = hh_ue_avg_20_50*DPHI_BINS

h_h_total_integral_error_20_50 = 0
h_h_near_integral_error_20_50 = 0
h_h_away_integral_error_20_50 = 0
h_h_ue_integral_error_20_50 = hh_ue_avg_error_20_50*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_h_total_integral_20_50 += h_h_dphi_20_50.GetBinContent(bin_num)
    h_h_total_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
    part = h_h_dphi_20_50.GetBinContent(bin_num) - hh_ue_avg_20_50
    if part < 0:
        continue
    if bin_num < 9:
        h_h_near_integral_20_50 += part
        h_h_near_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
        h_h_near_integral_error_20_50 += hh_ue_avg_error_20_50**2
    else:
        h_h_away_integral_20_50 += part
        h_h_away_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
        h_h_away_integral_error_20_50 += hh_ue_avg_error_20_50**2
h_h_total_integral_error_20_50 = math.sqrt(h_h_total_integral_error_20_50)
h_h_near_integral_error_20_50 = math.sqrt(h_h_near_integral_error_20_50)
h_h_away_integral_error_20_50 = math.sqrt(h_h_away_integral_error_20_50)

output_file.cd()
h_lambda_2d_subtracted_20_50.Write()

h_lambda_dphi_subtracted_20_50.Write()
h_lambda_dphi_subtracted_20_50_zeroed.Write()
h_h_dphi_20_50.Write()
h_h_dphi_20_50_zeroed.Write()

ue_upper_line_20_50.Write()
ue_line_20_50.Write()
ue_lower_line_20_50.Write()
zero_upper_line_20_50.Write()
zero_line_20_50.Write()
zero_lower_line_20_50.Write()

hh_ue_upper_line_20_50.Write()
hh_ue_line_20_50.Write()
hh_ue_lower_line_20_50.Write()
hh_zero_upper_line_20_50.Write()
hh_zero_line_20_50.Write()
hh_zero_lower_line_20_50.Write()


near_ratio_20_50 = h_lambda_near_integral_20_50/h_h_near_integral_20_50
away_ratio_20_50 = h_lambda_away_integral_20_50/h_h_away_integral_20_50
ue_ratio_20_50 = h_lambda_ue_integral_20_50/h_h_ue_integral_20_50
total_ratio_20_50 = h_lambda_total_integral_20_50/h_h_total_integral_20_50

near_ratio_error_20_50 = near_ratio_20_50*math.sqrt((h_lambda_near_integral_error_20_50/h_lambda_near_integral_20_50)**2
                                                 + (h_h_near_integral_error_20_50/h_h_near_integral_20_50)**2)
away_ratio_error_20_50 = away_ratio_20_50*math.sqrt((h_lambda_away_integral_error_20_50/h_lambda_away_integral_20_50)**2
                                                 + (h_h_away_integral_error_20_50/h_h_away_integral_20_50)**2)
ue_ratio_error_20_50 = ue_ratio_20_50*math.sqrt((h_lambda_ue_integral_error_20_50/h_lambda_ue_integral_20_50)**2
                                                 + (h_h_ue_integral_error_20_50/h_h_ue_integral_20_50)**2)
total_ratio_error_20_50 = total_ratio_20_50*math.sqrt((h_lambda_total_integral_error_20_50/h_lambda_total_integral_20_50)**2
                                                 + (h_h_total_integral_error_20_50/h_h_total_integral_20_50)**2)


############################################################################################################
############################################################################################################
########################################## 50 - 80 CENTRALITY ##############################################
############################################################################################################
############################################################################################################

input_file_50_80 = rt.TFile("~/OneDrive/Research/Output/lambda-over-hadron/data/res_cent_50_80.root")
input_list_50_80 = input_file_50_80.Get("h-lambda")
input_file_50_80.Close()


trig_dist_50_80 = input_list_50_80.FindObject("fTriggerDist_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fTriggerDist")
lambda_dist_50_80 = input_list_50_80.FindObject("fTriggeredLambdaDist")
lambda_ls_dist_50_80 = input_list_50_80.FindObject("fTriggeredLambdaLSDist")


h_h_50_80 = input_list_50_80.FindObject("fDphiHH_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHH")
h_h_mixed_50_80 = input_list_50_80.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHHMixed")

h_lambda_50_80 = input_list_50_80.FindObject("fDphiHLambda_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHLambda")
h_lambda_mixed_50_80 = input_list_50_80.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_50_80.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
lambda_ls_dist_50_80.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_50_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_50_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_ls_dist_50_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


# Getting the single-particle trigger distributions
trig_pt_dist_50_80 = trig_dist_50_80.Projection(0).Clone("trig_pt_dist_50_80")
trig_phi_dist_50_80 = trig_dist_50_80.Projection(1).Clone("trig_phi_dist_50_80")
trig_eta_dist_50_80 = trig_dist_50_80.Projection(2).Clone("trig_eta_dist_50_80")
trig_2d_dist_50_80 = trig_dist_50_80.Projection(0, 3).Clone("trig_2d_dist_50_80")

# Get total number of triggers (used for per-trigger normalization)
num_trigs_50_80 = trig_2d_dist_50_80.Integral()

output_file.cd()
trig_pt_dist_50_80.Write()
trig_phi_dist_50_80.Write()
trig_eta_dist_50_80.Write()
trig_2d_dist_50_80.Write()

### SIGNAL ANALYSIS SECTION ###

lambda_mass_dist_50_80 = lambda_dist_50_80.Projection(3).Clone("lambda_mass_dist_50_80")
lambda_mass_ls_dist_50_80 = lambda_ls_dist_50_80.Projection(3).Clone("lambda_mass_ls_dist_50_80")


# scale LS to match lambda dist in RSB
left_rsb_bin_50_80 = lambda_mass_dist_50_80.FindBin(RSB_MIN_50_80)
right_rsb_bin_50_80 = lambda_mass_dist_50_80.FindBin(RSB_MAX_50_80)
lambda_mass_ls_dist_50_80.Scale(lambda_mass_dist_50_80.Integral(left_rsb_bin_50_80, right_rsb_bin_50_80)/lambda_mass_ls_dist_50_80.Integral(left_rsb_bin_50_80, right_rsb_bin_50_80))

residual_50_80 = lambda_mass_dist_50_80.Clone("residual_50_80")
residual_50_80.Add(lambda_mass_ls_dist_50_80, -1)

RSB_region_50_80 = rt.TBox(RSB_MIN_50_80, 0, RSB_MAX_50_80, lambda_mass_dist_50_80.GetMaximum()*1.055)
RSB_region_50_80.SetLineColor(rt.kRed)
RSB_region_50_80.SetFillColor(rt.kRed)
RSB_region_50_80.SetFillStyle(3003)


RSB_min_line_50_80 = rt.TLine(RSB_MIN_50_80, 0, RSB_MIN_50_80, lambda_mass_dist_50_80.GetMaximum()*1.05)
RSB_min_line_50_80.SetLineColor(rt.kRed)
RSB_min_line_50_80.SetLineWidth(2)
RSB_min_line_50_80.SetLineStyle(2)


RSB_max_line_50_80 = rt.TLine(RSB_MAX_50_80, 0, RSB_MAX_50_80, lambda_mass_dist_50_80.GetMaximum()*1.05)
RSB_max_line_50_80.SetLineColor(rt.kRed)
RSB_max_line_50_80.SetLineWidth(2)
RSB_max_line_50_80.SetLineStyle(2)


left_signal_bin_50_80 = lambda_mass_dist_50_80.FindBin(SIG_MIN)
right_signal_bin_50_80 = lambda_mass_dist_50_80.FindBin(SIG_MAX)

lambda_bg_50_80 = 0
lambda_total_50_80 = 0
for bin_num in range(left_signal_bin_50_80, right_signal_bin_50_80 + 1):
    lambda_total_50_80 += lambda_mass_dist_50_80.GetBinContent(bin_num)
    lambda_bg_50_80 += lambda_mass_ls_dist_50_80.GetBinContent(bin_num)

lambda_signal_50_80 = lambda_total_50_80 - lambda_bg_50_80
lambda_signal_total_ratio_50_80 = lambda_signal_50_80/lambda_total_50_80

scale_factor_50_80 = residual_50_80.Integral(1, left_rsb_bin_50_80)/residual_50_80.Integral(left_signal_bin_50_80, right_signal_bin_50_80)

output_file.cd()
RSB_region_50_80.Write("RSB_region_50_80")
RSB_min_line_50_80.Write("RSB_min_line_50_80")
RSB_max_line_50_80.Write("RSB_max_line_50_80")
lambda_mass_dist_50_80.Write()
lambda_mass_ls_dist_50_80.Write()
residual_50_80.Write()
RSB_region_50_80.Write()
RSB_min_line_50_80.Write()
RSB_max_line_50_80.Write()

### MIXED EVENT CORRECTION SECTION ###

axes = arr.array('i', [2, 3, 4, 5])
h_lambda_50_80 = h_lambda_50_80.Projection(4, axes)
h_lambda_mixed_50_80 = h_lambda_mixed_50_80.Projection(4, axes)

h_h_50_80 = h_h_50_80.Projection(2, 3, 4)
h_h_mixed_50_80 = h_h_mixed_50_80.Projection(2, 3, 4)


# Setting up 2-d correlation plots before the mixed event correction
h_lambda_2d_nomixcor_50_80 = h_lambda_50_80.Projection(0, 1).Clone("h_lambda_2d_nomixcor_50_80")
h_lambda_mixed_2d_50_80 = h_lambda_mixed_50_80.Projection(0, 1).Clone("h_lambda_mixed_2d_50_80")

h_h_2d_nomixcor_50_80 = h_h_50_80.Project3D("xye").Clone("h_h_2d_nomixcor_50_80")
h_h_mixed_2d_50_80 = h_h_mixed_50_80.Project3D("xye").Clone("h_h_mixed_2d_50_80")

h_lambda_2d_mixcor_sig_50_80 = make_mixed_corrections(h_lambda_50_80, h_lambda_mixed_50_80, SIG_MIN, SIG_MAX)
h_lambda_2d_mixcor_rsb_50_80 = make_mixed_corrections(h_lambda_50_80, h_lambda_mixed_50_80, RSB_MIN_50_80, RSB_MAX_50_80)

h_h_2d_mixcor_50_80 = make_mixed_corrections(h_h_50_80, h_h_mixed_50_80, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_50_80.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_50_80.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)
h_h_2d_mixcor_50_80.GetXaxis().SetRangeUser(-DELTA_ETA_MAX, DELTA_ETA_MAX)




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_50_80.Scale(1.0/num_trigs_50_80)
h_lambda_2d_mixcor_rsb_50_80.Scale(1.0/num_trigs_50_80)
h_h_2d_mixcor_50_80.Scale(1.0/num_trigs_50_80)

# scaling by total signal/signal region done here
h_lambda_2d_mixcor_sig_50_80.Scale(scale_factor_50_80)
h_lambda_2d_mixcor_rsb_50_80.Scale(scale_factor_50_80)
h_lambda_2d_mixcor_sig_50_80.Scale(PID_CORRECTION)
h_lambda_2d_mixcor_rsb_50_80.Scale(PID_CORRECTION)

output_file.cd()
h_lambda_2d_nomixcor_50_80.Write()
h_lambda_mixed_2d_50_80.Write()
h_h_2d_nomixcor_50_80.Write()
h_h_mixed_2d_50_80.Write()
h_lambda_2d_mixcor_sig_50_80.Write()
h_lambda_2d_mixcor_rsb_50_80.Write()
h_h_2d_mixcor_50_80.Write()


# SIDEBAND SUBTRACTION SECTION

# Normalize side band to 1
h_lambda_2d_mixcor_rsb_50_80.Scale(1/h_lambda_2d_mixcor_rsb_50_80.Integral())


# using RSB for sideband subtraction
h_lambda_2d_subtracted_50_80 = h_lambda_2d_mixcor_sig_50_80.Clone("h_lambda_2d_subtracted_50_80")
bg_integral_50_80 = (1 - lambda_signal_total_ratio_50_80)*h_lambda_2d_subtracted_50_80.Integral()
h_lambda_2d_subtracted_50_80.Add(h_lambda_2d_mixcor_rsb_50_80, -bg_integral_50_80)


# INTEGRAL AND RATIO SECTION

h_lambda_dphi_subtracted_50_80 = h_lambda_2d_subtracted_50_80.ProjectionY("h_lambda_dphi_subtracted_50_80")

if USE_AVG_4:
    ue_line_50_80 = rt.TF1("ue_line_50_80", "pol0", -2, 6)
    ue_upper_line_50_80 = rt.TF1("ue_upper_line_50_80", "pol0", -2, 6)
    ue_lower_line_50_80 = rt.TF1("ue_lower_line_50_80", "pol0", -2, 6)
    zero_line_50_80 = rt.TF1("zero_line_50_80", "pol0", -2, 6)
    zero_upper_line_50_80 = rt.TF1("zero_upper_line_50_80", "pol0", -2, 6)
    zero_lower_line_50_80 = rt.TF1("zero_lower_line_50_80", "pol0", -2, 6)
    ue_avg_50_80 = (h_lambda_dphi_subtracted_50_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(16))/4

    ue_avg_error_50_80 = (1/4)*(math.sqrt(h_lambda_dphi_subtracted_50_80.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_50_80.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(16)**2))


    ue_line_50_80.SetParameter(0, ue_avg_50_80)
    ue_upper_line_50_80.SetParameter(0, ue_avg_50_80 + ue_avg_error_50_80)
    ue_lower_line_50_80.SetParameter(0, ue_avg_50_80 - ue_avg_error_50_80)

    zero_line_50_80.SetParameter(0, 0)
    zero_upper_line_50_80.SetParameter(0, ue_avg_error_50_80)
    zero_lower_line_50_80.SetParameter(0, -ue_avg_error_50_80)

elif USE_AVG_6:
    ue_line_50_80 = rt.TF1("ue_line_50_80", "pol0", -2, 6)
    ue_upper_line_50_80 = rt.TF1("ue_upper_line_50_80", "pol0", -2, 6)
    ue_lower_line_50_80 = rt.TF1("ue_lower_line_50_80", "pol0", -2, 6)
    zero_line_50_80 = rt.TF1("zero_line_50_80", "pol0", -2, 6)
    zero_upper_line_50_80 = rt.TF1("zero_upper_line_50_80", "pol0", -2, 6)
    zero_lower_line_50_80 = rt.TF1("zero_lower_line_50_80", "pol0", -2, 6)
    ue_avg_50_80 = (h_lambda_dphi_subtracted_50_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(16))/6

    ue_avg_error_50_80 = (1/6)*(math.sqrt(h_lambda_dphi_subtracted_50_80.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_50_80.GetBinError(2)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(7)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_50_80.GetBinError(16)**2))


    ue_line_50_80.SetParameter(0, ue_avg_50_80)
    ue_upper_line_50_80.SetParameter(0, ue_avg_50_80 + ue_avg_error_50_80)
    ue_lower_line_50_80.SetParameter(0, ue_avg_50_80 - ue_avg_error_50_80)

    zero_line_50_80.SetParameter(0, 0)
    zero_upper_line_50_80.SetParameter(0, ue_avg_error_50_80)
    zero_lower_line_50_80.SetParameter(0, -ue_avg_error_50_80)

elif USE_ZYAM:
    ue_line_50_80 = rt.TF1("ue_line_50_80", "pol0", -2, 6)
    ue_upper_line_50_80 = rt.TF1("ue_upper_line_50_80", "pol0", -2, 6)
    ue_lower_line_50_80 = rt.TF1("ue_lower_line_50_80", "pol0", -2, 6)
    zero_line_50_80 = rt.TF1("zero_line_50_80", "pol0", -2, 6)
    zero_upper_line_50_80 = rt.TF1("zero_upper_line_50_80", "pol0", -2, 6)
    zero_lower_line_50_80 = rt.TF1("zero_lower_line_50_80", "pol0", -2, 6)
    min_bin = h_lambda_dphi_subtracted_50_80.GetMinimumBin()
    ue_avg_50_80 = h_lambda_dphi_subtracted_50_80.GetBinContent(min_bin)
    ue_avg_error_50_80 = h_lambda_dphi_subtracted_50_80.GetBinError(min_bin)


    ue_line_50_80.SetParameter(0, ue_avg_50_80)
    ue_upper_line_50_80.SetParameter(0, ue_avg_50_80 + ue_avg_error_50_80)
    ue_lower_line_50_80.SetParameter(0, ue_avg_50_80 - ue_avg_error_50_80)

    zero_line_50_80.SetParameter(0, 0)
    zero_upper_line_50_80.SetParameter(0, ue_avg_error_50_80)
    zero_lower_line_50_80.SetParameter(0, -ue_avg_error_50_80)
else:
    raise NotImplementedError("UE line mode not supported")


h_lambda_dphi_subtracted_50_80_zeroed = h_lambda_dphi_subtracted_50_80.Clone("h_lambda_dphi_subtracted_50_80_zeroed")
h_lambda_dphi_subtracted_50_80_zeroed.Add(ue_line_50_80, -1)


DPHI_BINS = h_lambda_dphi_subtracted_50_80.GetNbinsX()

h_lambda_total_integral_50_80 = 0
h_lambda_near_integral_50_80 = 0
h_lambda_away_integral_50_80 = 0
h_lambda_ue_integral_50_80 = ue_avg_50_80*DPHI_BINS

h_lambda_total_integral_error_50_80 = 0
h_lambda_near_integral_error_50_80 = 0
h_lambda_away_integral_error_50_80 = 0
h_lambda_ue_integral_error_50_80 = ue_avg_error_50_80*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_lambda_total_integral_50_80 += h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num)
    h_lambda_total_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2
    part = h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num) - ue_avg_50_80
    if part < 0:
        continue
    if bin_num < 9:
        h_lambda_near_integral_50_80 += part
        h_lambda_near_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2
        h_lambda_near_integral_error_50_80 += ue_avg_error_50_80**2
    else:
        h_lambda_away_integral_50_80 += part
        h_lambda_away_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2
        h_lambda_away_integral_error_50_80 += ue_avg_error_50_80**2

h_lambda_total_integral_error_50_80 = math.sqrt(h_lambda_total_integral_error_50_80)
h_lambda_near_integral_error_50_80 = math.sqrt(h_lambda_near_integral_error_50_80)
h_lambda_away_integral_error_50_80 = math.sqrt(h_lambda_away_integral_error_50_80)

h_h_dphi_50_80 = h_h_2d_mixcor_50_80.ProjectionY("h_h_dphi_50_80")

if USE_AVG_4:
    hh_ue_line_50_80 = rt.TF1("hh_ue_line_50_80", "pol0", -2, 6)
    hh_ue_upper_line_50_80 = rt.TF1("hh_ue_upper_line_50_80", "pol0", -2, 6)
    hh_ue_lower_line_50_80 = rt.TF1("hh_ue_lower_line_50_80", "pol0", -2, 6)
    hh_zero_line_50_80 = rt.TF1("hh_zero_line_50_80", "pol0", -2, 6)
    hh_zero_upper_line_50_80 = rt.TF1("hh_zero_upper_line_50_80", "pol0", -2, 6)
    hh_zero_lower_line_50_80 = rt.TF1("hh_zero_lower_line_50_80", "pol0", -2, 6)
    hh_ue_avg_50_80 = (h_h_dphi_50_80.GetBinContent(1) 
                   + h_h_dphi_50_80.GetBinContent(8)
                   + h_h_dphi_50_80.GetBinContent(9)
                   + h_h_dphi_50_80.GetBinContent(16))/4

    hh_ue_avg_error_50_80 = (1/4)*(math.sqrt(h_h_dphi_50_80.GetBinError(1)**2 
                   + h_h_dphi_50_80.GetBinError(8)**2
                   + h_h_dphi_50_80.GetBinError(9)**2
                   + h_h_dphi_50_80.GetBinError(16)**2))


    hh_ue_line_50_80.SetParameter(0, hh_ue_avg_50_80)
    hh_ue_upper_line_50_80.SetParameter(0, hh_ue_avg_50_80 + hh_ue_avg_error_50_80)
    hh_ue_lower_line_50_80.SetParameter(0, hh_ue_avg_50_80 - hh_ue_avg_error_50_80)

    hh_zero_line_50_80.SetParameter(0, 0)
    hh_zero_upper_line_50_80.SetParameter(0, hh_ue_avg_error_50_80)
    hh_zero_lower_line_50_80.SetParameter(0, -hh_ue_avg_error_50_80)

elif USE_AVG_6:
    hh_ue_line_50_80 = rt.TF1("hh_ue_line_50_80", "pol0", -2, 6)
    hh_ue_upper_line_50_80 = rt.TF1("hh_ue_upper_line_50_80", "pol0", -2, 6)
    hh_ue_lower_line_50_80 = rt.TF1("hh_ue_lower_line_50_80", "pol0", -2, 6)
    hh_zero_line_50_80 = rt.TF1("hh_zero_line_50_80", "pol0", -2, 6)
    hh_zero_upper_line_50_80 = rt.TF1("hh_zero_upper_line_50_80", "pol0", -2, 6)
    hh_zero_lower_line_50_80 = rt.TF1("hh_zero_lower_line_50_80", "pol0", -2, 6)
    hh_ue_avg_50_80 = (h_h_dphi_50_80.GetBinContent(1) 
                   + h_h_dphi_50_80.GetBinContent(2)
                   + h_h_dphi_50_80.GetBinContent(7)
                   + h_h_dphi_50_80.GetBinContent(8)
                   + h_h_dphi_50_80.GetBinContent(9)
                   + h_h_dphi_50_80.GetBinContent(16))/6

    hh_ue_avg_error_50_80 = (1/6)*(math.sqrt(h_h_dphi_50_80.GetBinError(1)**2 
                   + h_h_dphi_50_80.GetBinError(2)**2
                   + h_h_dphi_50_80.GetBinError(7)**2
                   + h_h_dphi_50_80.GetBinError(8)**2
                   + h_h_dphi_50_80.GetBinError(9)**2
                   + h_h_dphi_50_80.GetBinError(16)**2))


    hh_ue_line_50_80.SetParameter(0, hh_ue_avg_50_80)
    hh_ue_upper_line_50_80.SetParameter(0, hh_ue_avg_50_80 + hh_ue_avg_error_50_80)
    hh_ue_lower_line_50_80.SetParameter(0, hh_ue_avg_50_80 - hh_ue_avg_error_50_80)

    hh_zero_line_50_80.SetParameter(0, 0)
    hh_zero_upper_line_50_80.SetParameter(0, hh_ue_avg_error_50_80)
    hh_zero_lower_line_50_80.SetParameter(0, -hh_ue_avg_error_50_80)
elif USE_ZYAM:
    hh_ue_line_50_80 = rt.TF1("hh_ue_line_50_80", "pol0", -2, 6)
    hh_ue_upper_line_50_80 = rt.TF1("hh_ue_upper_line_50_80", "pol0", -2, 6)
    hh_ue_lower_line_50_80 = rt.TF1("hh_ue_lower_line_50_80", "pol0", -2, 6)
    hh_zero_line_50_80 = rt.TF1("hh_zero_line_50_80", "pol0", -2, 6)
    hh_zero_upper_line_50_80 = rt.TF1("hh_zero_upper_line_50_80", "pol0", -2, 6)
    hh_zero_lower_line_50_80 = rt.TF1("hh_zero_lower_line_50_80", "pol0", -2, 6)
    
    min_bin = h_h_dphi_50_80.GetMinimumBin()
    hh_ue_avg_50_80 = h_h_dphi_50_80.GetBinContent(min_bin)
    hh_ue_avg_error_50_80 = h_h_dphi_50_80.GetBinError(min_bin)

    hh_ue_line_50_80.SetParameter(0, hh_ue_avg_50_80)
    hh_ue_upper_line_50_80.SetParameter(0, hh_ue_avg_50_80 + hh_ue_avg_error_50_80)
    hh_ue_lower_line_50_80.SetParameter(0, hh_ue_avg_50_80 - hh_ue_avg_error_50_80)

    hh_zero_line_50_80.SetParameter(0, 0)
    hh_zero_upper_line_50_80.SetParameter(0, hh_ue_avg_error_50_80)
    hh_zero_lower_line_50_80.SetParameter(0, -hh_ue_avg_error_50_80)
else:
    raise NotImplementedError("UE line mode not supported")

h_h_dphi_50_80_zeroed = h_h_dphi_50_80.Clone("h_h_dphi_50_80_zeroed")
h_h_dphi_50_80_zeroed.Add(hh_ue_line_50_80, -1)

h_h_total_integral_50_80 = 0
h_h_near_integral_50_80 = 0
h_h_away_integral_50_80 = 0
h_h_ue_integral_50_80 = hh_ue_avg_50_80*DPHI_BINS

h_h_total_integral_error_50_80 = 0
h_h_near_integral_error_50_80 = 0
h_h_away_integral_error_50_80 = 0
h_h_ue_integral_error_50_80 = hh_ue_avg_error_50_80*DPHI_BINS

for bin_num in range(1, DPHI_BINS + 1):
    h_h_total_integral_50_80 += h_h_dphi_50_80.GetBinContent(bin_num)
    h_h_total_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
    part = h_h_dphi_50_80.GetBinContent(bin_num) - hh_ue_avg_50_80
    if part < 0:
        continue
    if bin_num < 9:
        h_h_near_integral_50_80 += part
        h_h_near_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
        h_h_near_integral_error_50_80 += hh_ue_avg_error_50_80**2
    else:
        h_h_away_integral_50_80 += part
        h_h_away_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
        h_h_away_integral_error_50_80 += hh_ue_avg_error_50_80**2
h_h_total_integral_error_50_80 = math.sqrt(h_h_total_integral_error_50_80)
h_h_near_integral_error_50_80 = math.sqrt(h_h_near_integral_error_50_80)
h_h_away_integral_error_50_80 = math.sqrt(h_h_away_integral_error_50_80)

output_file.cd()
h_lambda_2d_subtracted_50_80.Write()

h_lambda_dphi_subtracted_50_80.Write()
h_lambda_dphi_subtracted_50_80_zeroed.Write()
h_h_dphi_50_80.Write()
h_h_dphi_50_80_zeroed.Write()

ue_upper_line_50_80.Write()
ue_line_50_80.Write()
ue_lower_line_50_80.Write()
zero_upper_line_50_80.Write()
zero_line_50_80.Write()
zero_lower_line_50_80.Write()

hh_ue_upper_line_50_80.Write()
hh_ue_line_50_80.Write()
hh_ue_lower_line_50_80.Write()
hh_zero_upper_line_50_80.Write()
hh_zero_line_50_80.Write()
hh_zero_lower_line_50_80.Write()

near_ratio_50_80 = h_lambda_near_integral_50_80/h_h_near_integral_50_80
away_ratio_50_80 = h_lambda_away_integral_50_80/h_h_away_integral_50_80
ue_ratio_50_80 = h_lambda_ue_integral_50_80/h_h_ue_integral_50_80
total_ratio_50_80 = h_lambda_total_integral_50_80/h_h_total_integral_50_80

near_ratio_error_50_80 = near_ratio_50_80*math.sqrt((h_lambda_near_integral_error_50_80/h_lambda_near_integral_50_80)**2
                                                 + (h_h_near_integral_error_50_80/h_h_near_integral_50_80)**2)
away_ratio_error_50_80 = away_ratio_50_80*math.sqrt((h_lambda_away_integral_error_50_80/h_lambda_away_integral_50_80)**2
                                                 + (h_h_away_integral_error_50_80/h_h_away_integral_50_80)**2)
ue_ratio_error_50_80 = ue_ratio_50_80*math.sqrt((h_lambda_ue_integral_error_50_80/h_lambda_ue_integral_50_80)**2
                                                 + (h_h_ue_integral_error_50_80/h_h_ue_integral_50_80)**2)
total_ratio_error_50_80 = total_ratio_50_80*math.sqrt((h_lambda_total_integral_error_50_80/h_lambda_total_integral_50_80)**2
                                                 + (h_h_total_integral_error_50_80/h_h_total_integral_50_80)**2)


############################################################################################################
############################################################################################################
############################################ FINAL RATIOS ##################################################
############################################################################################################
############################################################################################################

mult_list = arr.array('d', [35, 65, 90])
mult_error_list = arr.array('d', [15, 15, 10])


h_lambda_near_list = arr.array('d', [h_lambda_near_integral_50_80, h_lambda_near_integral_20_50, h_lambda_near_integral_0_20])
h_lambda_near_error_list = arr.array('d', [h_lambda_near_integral_error_50_80, h_lambda_near_integral_error_20_50, h_lambda_near_integral_error_0_20])

h_lambda_away_list = arr.array('d', [h_lambda_away_integral_50_80, h_lambda_away_integral_20_50, h_lambda_away_integral_0_20])
h_lambda_away_error_list = arr.array('d', [h_lambda_away_integral_error_50_80, h_lambda_away_integral_error_20_50, h_lambda_away_integral_error_0_20])

h_lambda_ue_list = arr.array('d', [h_lambda_ue_integral_50_80, h_lambda_ue_integral_20_50, h_lambda_ue_integral_0_20])
h_lambda_ue_error_list = arr.array('d', [h_lambda_ue_integral_error_50_80, h_lambda_ue_integral_error_20_50, h_lambda_ue_integral_error_0_20])

h_lambda_total_list = arr.array('d', [h_lambda_total_integral_50_80, h_lambda_total_integral_20_50, h_lambda_total_integral_0_20])
h_lambda_total_error_list = arr.array('d', [h_lambda_total_integral_error_50_80, h_lambda_total_integral_error_20_50, h_lambda_total_integral_error_0_20])

h_h_near_list = arr.array('d', [h_h_near_integral_50_80, h_h_near_integral_20_50, h_h_near_integral_0_20])
h_h_near_error_list = arr.array('d', [h_h_near_integral_error_50_80, h_h_near_integral_error_20_50, h_h_near_integral_error_0_20])

h_h_away_list = arr.array('d', [h_h_away_integral_50_80, h_h_away_integral_20_50, h_h_away_integral_0_20])
h_h_away_error_list = arr.array('d', [h_h_away_integral_error_50_80, h_h_away_integral_error_20_50, h_h_away_integral_error_0_20])

h_h_ue_list = arr.array('d', [h_h_ue_integral_50_80, h_h_ue_integral_20_50, h_h_ue_integral_0_20])
h_h_ue_error_list = arr.array('d', [h_h_ue_integral_error_50_80, h_h_ue_integral_error_20_50, h_h_ue_integral_error_0_20])

h_h_total_list = arr.array('d', [h_h_total_integral_50_80, h_h_total_integral_20_50, h_h_total_integral_0_20])
h_h_total_error_list = arr.array('d', [h_h_total_integral_error_50_80, h_h_total_integral_error_20_50, h_h_total_integral_error_0_20])


near_ratio_list = arr.array('d', [near_ratio_50_80, near_ratio_20_50, near_ratio_0_20])
near_ratio_error_list = arr.array('d', [near_ratio_error_0_20, near_ratio_error_20_50, near_ratio_error_50_80])

away_ratio_list = arr.array('d', [away_ratio_50_80, away_ratio_20_50, away_ratio_0_20])
away_ratio_error_list = arr.array('d', [away_ratio_error_50_80, away_ratio_error_20_50, away_ratio_error_0_20])

ue_ratio_list = arr.array('d', [ue_ratio_50_80, ue_ratio_20_50, ue_ratio_0_20])
ue_ratio_error_list = arr.array('d', [ue_ratio_error_50_80, ue_ratio_error_20_50, ue_ratio_error_0_20])

total_ratio_list = arr.array('d', [total_ratio_50_80, total_ratio_20_50, total_ratio_0_20])
total_ratio_error_list = arr.array('d', [total_ratio_error_50_80, total_ratio_error_20_50, total_ratio_error_0_20])


h_lambda_near_graph = rt.TGraphErrors(3, mult_list, h_lambda_near_list, mult_error_list, h_lambda_near_error_list)
h_lambda_away_graph = rt.TGraphErrors(3, mult_list, h_lambda_away_list, mult_error_list, h_lambda_away_error_list)
h_lambda_ue_graph = rt.TGraphErrors(3, mult_list, h_lambda_ue_list, mult_error_list, h_lambda_ue_error_list)
h_lambda_total_graph = rt.TGraphErrors(3, mult_list, h_lambda_total_list, mult_error_list, h_lambda_total_error_list)

h_h_near_graph = rt.TGraphErrors(3, mult_list, h_h_near_list, mult_error_list, h_h_near_error_list)
h_h_away_graph = rt.TGraphErrors(3, mult_list, h_h_away_list, mult_error_list, h_h_away_error_list)
h_h_ue_graph = rt.TGraphErrors(3, mult_list, h_h_ue_list, mult_error_list, h_h_ue_error_list)
h_h_total_graph = rt.TGraphErrors(3, mult_list, h_h_total_list, mult_error_list, h_h_total_error_list)

near_ratio_graph = rt.TGraphErrors(3, mult_list, near_ratio_list, mult_error_list, near_ratio_error_list)
away_ratio_graph = rt.TGraphErrors(3, mult_list, away_ratio_list, mult_error_list, away_ratio_error_list)
ue_ratio_graph = rt.TGraphErrors(3, mult_list, ue_ratio_list, mult_error_list, ue_ratio_error_list)
total_ratio_graph = rt.TGraphErrors(3, mult_list, total_ratio_list, mult_error_list, total_ratio_error_list)


mult_bin_widths = arr.array('d', [0.0, 20.0, 50.0, 80.0, 100.0])
near_ratio_hist = rt.TH1D("near_ratio_hist", "", 4, mult_bin_widths)
for i in range(3):
    near_ratio_hist.SetBinContent(i+2, near_ratio_list[i])
    near_ratio_hist.SetBinError(i+2, near_ratio_error_list[i])

output_file.cd()
h_lambda_near_graph.Write("h_lambda_near_graph")
h_lambda_away_graph.Write("h_lambda_away_graph")
h_lambda_ue_graph.Write("h_lambda_ue_graph")
h_lambda_total_graph.Write("h_lambda_total_graph")
h_h_near_graph.Write("h_h_near_graph")
h_h_away_graph.Write("h_h_away_graph")
h_h_ue_graph.Write("h_h_ue_graph")
h_h_total_graph.Write("h_h_total_graph")
near_ratio_graph.Write("near_ratio_graph")
away_ratio_graph.Write("away_ratio_graph")
ue_ratio_graph.Write("ue_ratio_graph")
total_ratio_graph.Write("total_ratio_graph")
near_ratio_hist.Write("near_ratio_hist")
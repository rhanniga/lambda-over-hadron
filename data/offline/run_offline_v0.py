#! /usr/bin/env python
import math
import configparser

import array as arr
import ROOT as rt


from sys import argv
from strangehelper_old import make_mixed_corrections,  get_straight_line, apply_two_track_correction

# epsilon used to avoid bin edge nightmares (if you pick a value that lies on bin edge, it defaults to right bin)
EPSILON = 0.00001

# read in the parameters from the config file
config = configparser.ConfigParser()
config.read(argv[1])


# decide whether to do sideband subtraction or not
DO_SIDEBAND_SUBTRACTION = config.getboolean("GENERAL", "DO_SIDEBAND_SUBTRACTION")

# decide whether to do analysis with single highest pt trigger
DO_HIGHEST_PT = config.getboolean("GENERAL", "DO_HIGHEST_PT")

# UE line method

USE_AVG_4 = config.getboolean("GENERAL", "USE_AVG_4")
USE_AVG_6 = config.getboolean("GENERAL", "USE_AVG_6")
USE_AVG_6_NONNEGATIVE = config.getboolean("GENERAL", "USE_AVG_6_NONNEGATIVE")
USE_ZYAM = config.getboolean("GENERAL", "USE_ZYAM")
USE_FIT = config.getboolean("GENERAL", "USE_FIT")
USE_V2 = config.getboolean("GENERAL", "USE_V2")
USE_VON = config.getboolean("GENERAL", "USE_VON")

assert sum([USE_AVG_4, USE_AVG_6, USE_AVG_6_NONNEGATIVE, USE_ZYAM, USE_FIT, USE_V2, USE_VON]) == 1, "Only select 1 method for UE line please"

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
RSB_MIN = config.getfloat("REGION_CUTS", "RSB_MIN")
RSB_MAX = config.getfloat("REGION_CUTS", "RSB_MAX") - EPSILON


# PID CUTS (already applied online, need to correct for % loss)
IS_NORMAL_PID = config.getboolean("PID_CUTS", "IS_NORMAL")
IS_WIDE_PID = config.getboolean("PID_CUTS", "IS_WIDE")
IS_NARROW_PID = config.getboolean("PID_CUTS", "IS_NARROW")
IS_TOF_PID = config.getboolean("PID_CUTS", "IS_TOF")

assert sum([IS_NORMAL_PID, IS_NARROW_PID, IS_WIDE_PID, IS_TOF_PID]) == 1, "Only 1 PID cut is applied online"

two_sigma_normal = 0.9544
three_sigma_normal = 0.9974

two_sigma_wide = 0.9948 # 2 * 1.4 sigma
three_sigma_wide = 1.0 # 3 * 1.4 sigma

two_sigma_narrow = 0.7700 # 2 * 0.6
three_sigma_narrow = 0.9282 # 3 * 0.6

if IS_NORMAL_PID:
    PID_CORRECTION = 1/(two_sigma_normal*three_sigma_normal)
if IS_WIDE_PID:
    PID_CORRECTION = 1/(two_sigma_wide*three_sigma_wide)
if IS_NARROW_PID:
    PID_CORRECTION = 1/(two_sigma_narrow*three_sigma_narrow)
if IS_TOF_PID:
    PID_CORRECTION = (1.02/(1.71e-1)) # found from fitting dphi ratios of final dists with and without TOF veto

FULL_REGION_MIN = 1.10
FULL_REGION_MAX = 1.132 - EPSILON

# need bin width scaling for yields extracted from integration
N_DPHI_BINS = 16
BIN_WIDTH = 2*rt.TMath.Pi()/N_DPHI_BINS


v2_fitpars_file_string = "v2_fitpars.txt"
# funny and inefficienct dict comprehension to practice my python prowess
# the exact ordering is: key: pedestal, pedestal error, trigger v2 (fixed to weighted avg), associated v2 (fixed to weighted avg)
v2_fitpars_dict = {
    line[0] : [float(line[1]), float(line[2]), float(line[3]), float(line[4])] for line in [line.split() for line in open(v2_fitpars_file_string).readlines()]
}

# go ahead and load the two-track correction templates

if "lowpt" in argv[1]:
    template_file = rt.TFile("templates/twotrack_template_lowpt.root")
elif "highpt" in argv[1]:
    template_file = rt.TFile("templates/twotrack_template_highpt.root")
else:
    template_file = rt.TFile("templates/twotrack_template.root")
TWOTRACK_TEMPLATE = template_file.Get("twotrack_template")


BRANCHING_RATIO = 0.639 # from PDG

# Output file containing all of the relevant results (long name )
if "signal" in argv[1]:
    output_file_string = "output/signal_variation/v0_" 
elif "sideband" in argv[1]:
    output_file_string = "output/sideband_variation/v0_"
elif "pid" in argv[1]:
    output_file_string = "output/pid_variation/v0_"
elif "yield" in argv[1]:
    output_file_string = "output/yield_variation/v0_"
elif "central" in argv[1]:
    output_file_string = "output/central_value/v0_"
else:
    print("Unrecognized string in input file name")
    quit()

output_file_string += ("highest_pt_" if DO_HIGHEST_PT else "") 
output_file_string += ("avg4_" if USE_AVG_4 else "") 
output_file_string += ("avg6_" if USE_AVG_6 else "") 
output_file_string += ("avg6nonneg_" if USE_AVG_6_NONNEGATIVE else "") 
output_file_string += ("fullfit_" if USE_FIT else "") 
output_file_string += ("v2_" if USE_V2 else "") 
output_file_string += ("von_" if USE_VON else "") 
output_file_string += ("zyam_" if USE_ZYAM else "") 
output_file_string += ("sideband_subtraction_" if DO_SIDEBAND_SUBTRACTION else "")
output_file_string += "rsb_" + str(RSB_MIN).replace(".", "") + "_" + str(RSB_MAX + EPSILON).replace(".", "") + "_"
output_file_string += "sig_" + str(SIG_MIN).replace(".", "") + "_" + str(SIG_MAX + EPSILON).replace(".", "") + "_"
output_file_string += "trig_" + str(TRIG_PT_LOW).replace(".", "") + "_" + str(TRIG_PT_HIGH + EPSILON).replace(".", "") + "_"
output_file_string += "assoc_" + str(ASSOC_PT_LOW).replace(".", "") + "_" + str(ASSOC_PT_HIGH + EPSILON).replace(".", "") + "_"
output_file_string += "delta_eta_" + str(DELTA_ETA_MAX + EPSILON).replace(".", "") + "_"
output_file_string += ("normal" if IS_NORMAL_PID else "" )
output_file_string += ("wide" if IS_WIDE_PID else "" )
output_file_string += ("narrow" if IS_NARROW_PID else "" )
output_file_string += ("tof" if IS_TOF_PID else "" )

PT_MODE = 0 # 0 = central, 1 = low, 2 = high

if "lowpt" in argv[1]:
    output_file_string += "_lowpt.root"
    PT_MODE = 1
elif "highpt" in argv[1]:
    output_file_string += "_highpt.root"
    PT_MODE = 2
else:
    output_file_string += ".root"

output_file = rt.TFile(output_file_string, "RECREATE")


############################################################################################################
############################################################################################################
############################################# 0-20 CENTRALITY ##############################################
############################################################################################################
############################################################################################################

if IS_NORMAL_PID:
    input_file_0_20 = rt.TFile("../online/v0/output/v0_normal_0_20.root")
elif IS_WIDE_PID:
    input_file_0_20 = rt.TFile("../online/v0/output/v0_wide_0_20.root")
elif IS_NARROW_PID:
    input_file_0_20 = rt.TFile("../online/v0/output/v0_narrow_0_20.root")
elif IS_TOF_PID:
    input_file_0_20 = rt.TFile("../online/v0/output/v0_tof_0_20.root")
else:
    print("ERROR: No PID cut specified")
    exit(1)

input_list_0_20 = input_file_0_20.Get("h-lambda")
input_file_0_20.Close()

trig_dist_0_20 = input_list_0_20.FindObject("fTriggerDistEff_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fTriggerDistEff")
lambda_dist_0_20 = input_list_0_20.FindObject("fTriggeredLambdaDist")


h_h_0_20 = input_list_0_20.FindObject("fDphiHHEff_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHHEff")
h_h_mixed_0_20 = input_list_0_20.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHHMixed")

h_lambda_0_20 = input_list_0_20.FindObject("fDphiHLambdaEff_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHLambdaEff")
h_lambda_mixed_0_20 = input_list_0_20.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_0_20.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_0_20.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_0_20.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_mixed_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_0_20.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_0_20.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


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

### SIGNAL ANALYSIS SECTION ###

lambda_mass_dist_0_20 = lambda_dist_0_20.Projection(3).Clone("lambda_mass_dist_0_20")

bin_1 = lambda_mass_dist_0_20.FindBin(1.102)
bin_2 = lambda_mass_dist_0_20.FindBin(1.132)
point_one = [1.102, lambda_mass_dist_0_20.GetBinContent(bin_1)]
point_two = [1.132, lambda_mass_dist_0_20.GetBinContent(bin_2)]
bg_starting_params_0_20 = get_straight_line(point_one, point_two)
lambda_mass_fit_0_20 = rt.TF1("lambda_mass_fit_0_20", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)", 1.102, 1.127)
lambda_mass_fit_0_20.SetNpx(1000)
lambda_mass_fit_0_20.SetParameter(0, 1.74e03)
lambda_mass_fit_0_20.SetParameter(1, 1.116)
lambda_mass_fit_0_20.SetParameter(2, 1.30576e-03 )
lambda_mass_fit_0_20.SetParameter(3, 1.804166e-03)
lambda_mass_fit_0_20.SetParameter(4, bg_starting_params_0_20[0])
lambda_mass_fit_0_20.SetParameter(5, bg_starting_params_0_20[1])
lambda_mass_dist_fit_0_20 = lambda_mass_dist_0_20.Clone("lambda_mass_dist_fit_0_20")


lambda_mass_bg_fit_0_20 = rt.TF1("bg_fit_0_20", "pol1", 1.096, 1.136)
lambda_mass_bg_fit_0_20.SetNpx(1000)
lambda_mass_bg_fit_0_20.SetParameter(0, bg_starting_params_0_20[0])
lambda_mass_bg_fit_0_20.SetParameter(1, bg_starting_params_0_20[1])
lambda_mass_dist_fit_0_20.Fit(lambda_mass_fit_0_20, "RS")

voigt_fit_0_20 = rt.TF1("voigt_fit_0_20", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.102, 1.127)
voigt_fit_0_20.SetNpx(1000)
voigt_fit_0_20.SetParameter(0, lambda_mass_fit_0_20.GetParameter(0))
voigt_fit_0_20.SetParameter(1, lambda_mass_fit_0_20.GetParameter(1))
voigt_fit_0_20.SetParameter(2, lambda_mass_fit_0_20.GetParameter(2))
voigt_fit_0_20.SetParameter(3, lambda_mass_fit_0_20.GetParameter(3))


residual_0_20 = lambda_mass_dist_0_20.Clone("residual_0_20")
residual_0_20.GetXaxis().SetRangeUser(1.094, 1.142)
residual_0_20.Add(lambda_mass_bg_fit_0_20, -1)


RSB_region_0_20 = rt.TBox(RSB_MIN, 0, RSB_MAX, lambda_mass_dist_0_20.GetMaximum()*1.055)
RSB_region_0_20.SetLineColor(rt.kRed)
RSB_region_0_20.SetFillColor(rt.kRed)
RSB_region_0_20.SetFillStyle(3003)


RSB_min_line_0_20 = rt.TLine(RSB_MIN, 0, RSB_MIN, lambda_mass_dist_0_20.GetMaximum()*1.05)
RSB_min_line_0_20.SetLineColor(rt.kRed)
RSB_min_line_0_20.SetLineWidth(2)
RSB_min_line_0_20.SetLineStyle(2)


RSB_max_line_0_20 = rt.TLine(RSB_MAX, 0, RSB_MAX, lambda_mass_dist_0_20.GetMaximum()*1.05)
RSB_max_line_0_20.SetLineColor(rt.kRed)
RSB_max_line_0_20.SetLineWidth(2)
RSB_max_line_0_20.SetLineStyle(2)

left_signal_bin_0_20 = lambda_mass_dist_0_20.FindBin(SIG_MIN)
right_signal_bin_0_20 = lambda_mass_dist_0_20.FindBin(SIG_MAX)

full_region_bin_1 = residual_0_20.FindBin(FULL_REGION_MIN)
full_region_bin_2 = residual_0_20.FindBin(FULL_REGION_MAX)

lambda_bg_0_20 = 0
lambda_total_0_20 = 0
for bin_num in range(left_signal_bin_0_20, right_signal_bin_0_20 + 1):
    lambda_total_0_20 += lambda_mass_dist_0_20.GetBinContent(bin_num)
    lambda_bg_0_20 += lambda_mass_bg_fit_0_20.Eval(lambda_mass_dist_0_20.GetBinCenter(bin_num))


lambda_signal_0_20 = lambda_total_0_20 - lambda_bg_0_20
lambda_signal_total_ratio_0_20 = lambda_signal_0_20/lambda_total_0_20

signal_scale_0_20 = residual_0_20.Integral(full_region_bin_1, full_region_bin_2)/residual_0_20.Integral(left_signal_bin_0_20, right_signal_bin_0_20)

output_file.cd()
RSB_region_0_20.Write("RSB_region_0_20")
RSB_min_line_0_20.Write("RSB_min_line_0_20")
RSB_max_line_0_20.Write("RSB_max_line_0_20")

lambda_mass_dist_0_20.Write()
lambda_mass_dist_fit_0_20.Write()
lambda_mass_fit_0_20.Write()
lambda_mass_bg_fit_0_20.Write()
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
h_lambda_2d_mixcor_rsb_0_20 = make_mixed_corrections(h_lambda_0_20, h_lambda_mixed_0_20, RSB_MIN, RSB_MAX)

h_h_2d_mixcor_0_20 = make_mixed_corrections(h_h_0_20, h_h_mixed_0_20, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_0_20.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_0_20.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_h_2d_mixcor_0_20.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)

h_lambda_2d_mixcor_sig_0_20.SetName("h_lambda_2d_mixcor_sig_0_20")
h_lambda_2d_mixcor_rsb_0_20.SetName("h_lambda_2d_mixcor_rsb_0_20")
h_h_2d_mixcor_0_20.SetName("h_h_2d_mixcor_0_20")




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_0_20.Scale(1.0/num_trigs_0_20)
h_lambda_2d_mixcor_rsb_0_20.Scale(1.0/num_trigs_0_20)
h_h_2d_mixcor_0_20.Scale(1.0/num_trigs_0_20)

output_file.cd()
h_lambda_2d_nomixcor_0_20.Write()
h_lambda_mixed_2d_0_20.Write()
h_h_2d_nomixcor_0_20.Write()
h_h_mixed_2d_0_20.Write()
h_lambda_2d_mixcor_sig_0_20.Write()
h_lambda_2d_mixcor_rsb_0_20.Write()
h_h_2d_mixcor_0_20.Write()

# SIDEBAND SUBTRACTION SECTION
if DO_SIDEBAND_SUBTRACTION:

    # Normalize side band to 1
    h_lambda_2d_mixcor_rsb_0_20.Scale(1/h_lambda_2d_mixcor_rsb_0_20.Integral())


    # using RSB for sideband subtraction
    h_lambda_2d_subtracted_0_20 = h_lambda_2d_mixcor_sig_0_20.Clone("h_lambda_2d_subtracted_0_20")
    bg_integral_0_20 = (1 - lambda_signal_total_ratio_0_20)*h_lambda_2d_subtracted_0_20.Integral()
    h_lambda_2d_subtracted_0_20.Add(h_lambda_2d_mixcor_rsb_0_20, -bg_integral_0_20)

    # save the RSB deltaphi distribution
    h_lambda_dphi_rsb_0_20 = h_lambda_2d_mixcor_rsb_0_20.ProjectionY("h_lambda_dphi_rsb_0_20")
    output_file.cd()
    h_lambda_dphi_rsb_0_20.Write()

else:
    h_lambda_2d_subtracted_0_20 = h_lambda_2d_mixcor_sig_0_20.Clone("h_lambda_2d_subtracted_0_20")


# Correct for signal scale
h_lambda_2d_subtracted_0_20.Scale(signal_scale_0_20)
# Correct for PID loss
h_lambda_2d_subtracted_0_20.Scale(PID_CORRECTION)
# Correct for branching ratio
h_lambda_2d_subtracted_0_20.Scale(1/BRANCHING_RATIO)

# Correct for two-track efficiency
h_lambda_2d_subtracted_0_20 = apply_two_track_correction(h_lambda_2d_subtracted_0_20, TWOTRACK_TEMPLATE)



# INTEGRAL AND RATIO SECTION
h_lambda_dphi_subtracted_0_20 = h_lambda_2d_subtracted_0_20.ProjectionY("h_lambda_dphi_subtracted_0_20")

# this is completely insane but i need to scale the PID stuff exactly prior to fitting
if IS_WIDE_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_0_20.Scale(0.9223032488724654)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_0_20.Scale(0.9168476238506074)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_0_20.Scale(0.9293847448569807)
elif IS_NARROW_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_0_20.Scale(1.1940715128854646)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_0_20.Scale(1.1788960330443934)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_0_20.Scale(1.1970913624774955)
elif IS_TOF_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_0_20.Scale(0.9925605504303531)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_0_20.Scale(1.768552008127845)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_0_20.Scale(0.8240094208251248)





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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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

elif USE_FIT:

    # fitting to v2 component + four gaussians
    fit_string_0_20 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    fit_string_0_20 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    fit_function_0_20 = rt.TF1("fit_function_0_20", fit_string_0_20, -2, 6)
    fit_function_0_20.SetNpx(1000)
    if PT_MODE == 0:
        fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][0])
        fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][2])
        fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][0])
        fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][2])
        fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][0])
        fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][2])
        fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    fit_function_0_20.SetParLimits(3, 5e-4, 0.6)
    fit_function_0_20.SetParameter(3, 0.009)
    fit_function_0_20.FixParameter(4, 0)
    fit_function_0_20.SetParLimits(5, 0.1, 2)
    fit_function_0_20.SetParameter(5, 0.8)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    fit_function_0_20.SetParLimits(6, 0,  0.6)
    fit_function_0_20.SetParameter(6, 0)
    fit_function_0_20.FixParameter(7, 2*rt.TMath.Pi())
    fit_function_0_20.SetParLimits(8, 0.1, 2)
    fit_function_0_20.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    fit_function_0_20.SetParLimits(9, 5e-4, 0.6)
    fit_function_0_20.SetParameter(9, 0.003)
    fit_function_0_20.FixParameter(10, rt.TMath.Pi())
    fit_function_0_20.SetParLimits(11, 0.1, 2)
    fit_function_0_20.SetParameter(11, 1)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    fit_function_0_20.SetParLimits(12, 0, 0.6)
    fit_function_0_20.SetParameter(12, 0)
    fit_function_0_20.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    fit_function_0_20.SetParLimits(14, 0.1, 2)
    fit_function_0_20.SetParameter(14, 0.1)


    h_lambda_dphi_subtracted_with_fit_0_20 = h_lambda_dphi_subtracted_0_20.Clone("h_lambda_dphi_subtracted_with_fit_0_20")
    fit_result_0_20 = h_lambda_dphi_subtracted_with_fit_0_20.Fit(fit_function_0_20, "RS")

    v2_fit_0_20 = rt.TF1("v2_fit_0_20", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][3])

elif USE_V2:

    ue_avg_0_20 = (h_lambda_dphi_subtracted_0_20.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(16))/6

    v2_fit_0_20 = rt.TF1("v2_fit_0_20", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    ue_avg_0_20 = (h_lambda_dphi_subtracted_0_20.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_20.GetBinContent(16))/6

    v2_fit_0_20 = rt.TF1("v2_fit_0_20", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
    von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"

    von_fit_0_20 = rt.TF1("von_fit_0_20", von_fit_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

    if PT_MODE == 0:
        von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][0])
        von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][2])
        von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][3])
        # von_fit_0_20.SetParLimits(3, 0.015/3, 0.015*3)
        von_fit_0_20.SetParameter(3, 0.015)
        # von_fit_0_20.SetParLimits(4, 7.6/3, 7.6*3)
        von_fit_0_20.SetParameter(4, 7.6)
        # von_fit_0_20.SetParLimits(5, 0.014/3, 0.014*3)
        von_fit_0_20.SetParameter(5, 0.014)
        # von_fit_0_20.SetParLimits(6, 2.2/3, 2.2*3)
        von_fit_0_20.SetParameter(6, 2.2)
    elif PT_MODE == 1:
        von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][0])
        von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][2])
        von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][3])
        # von_fit_0_20.SetParLimits(3, 0.02/3, 0.02*3)
        von_fit_0_20.SetParameter(3, 0.0137579)
        # von_fit_0_20.SetParLimits(4, 4.6/3, 4.6*3)
        von_fit_0_20.SetParameter(4, 8.0)
        # von_fit_0_20.SetParLimits(5, 0.057/5, 0.057*5)
        von_fit_0_20.SetParameter(5, 1.33787e-2)
        # von_fit_0_20.SetParLimits(6, 2/4, 8)
        von_fit_0_20.SetParameter(6, 2.01496)
    elif PT_MODE == 2:
        von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][0])
        von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][2])
        von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][3])
        # von_fit_0_20.SetParLimits(3, 0.009/3, 0.009*3)
        von_fit_0_20.SetParameter(3, 0.009)
        # von_fit_0_20.SetParLimits(4, 8.88/3, 8.88*3)
        von_fit_0_20.SetParameter(4, 8.88)
        # von_fit_0_20.SetParLimits(5, 0.0076/5, 0.0076*5)
        von_fit_0_20.SetParameter(5, 0.0076)
        # von_fit_0_20.SetParLimits(6, 2.55/3, 2.55*3)
        von_fit_0_20.SetParameter(6, 2.55)


    h_lambda_dphi_subtracted_with_fit_0_20 = h_lambda_dphi_subtracted_0_20.Clone("h_lambda_dphi_subtracted_with_fit_0_20")
    h_lambda_dphi_subtracted_with_fit_0_20.GetYaxis().SetRangeUser(0.9*h_lambda_dphi_subtracted_with_fit_0_20.GetMinimum(), 1.2*h_lambda_dphi_subtracted_with_fit_0_20.GetMaximum())

    if PT_MODE == 1: # low pt has away-side width too large for von mises fit, only works if we symmetrize the peaks
        fit_result_0_20 = h_lambda_dphi_subtracted_with_fit_0_20.Fit(von_fit_0_20, "R", "", -rt.TMath.Pi()/2, rt.TMath.Pi() - EPSILON)
    else:
        fit_result_0_20 = h_lambda_dphi_subtracted_with_fit_0_20.Fit(von_fit_0_20, "R")

    if PT_MODE == 0:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][0])
        v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][1])
        v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][2])
        v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_20_trigger_4_8_assoc_25_4"][3])

else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_lambda_dphi_subtracted_0_20_zeroed = h_lambda_dphi_subtracted_0_20.Clone("h_lambda_dphi_subtracted_0_20_zeroed")
    h_lambda_dphi_subtracted_0_20_zeroed.Add(v2_fit_0_20, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_lambda_total_integral_error_0_20 = fit_function_0_20.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_0_20 = fit_function_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_0_20 = v2_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_0_20 = v2_fit_0_20.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_0_20_halved = h_lambda_ue_integral_error_0_20/2
    h_lambda_near_integral_error_0_20 = math.sqrt((fit_function_0_20.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_near_integral_0_20 = fit_function_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_0_20 = math.sqrt((fit_function_0_20.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_away_integral_0_20 = fit_function_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_lambda_total_integral_error_0_20 = von_fit_0_20.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_0_20 = von_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_0_20 = v2_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_0_20 = v2_fit_0_20.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_20_halved = h_lambda_ue_integral_error_0_20/2
    h_lambda_near_integral_error_0_20 = math.sqrt((von_fit_0_20.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_near_integral_0_20 = von_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_0_20 = math.sqrt((von_fit_0_20.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_away_integral_0_20 = von_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_dphi_subtracted_0_20_zeroed = h_lambda_dphi_subtracted_0_20.Clone("h_lambda_dphi_subtracted_0_20_zeroed")
    h_lambda_dphi_subtracted_0_20_zeroed.Add(v2_fit_0_20, -1)

elif USE_V2:

    DPHI_BINS = h_lambda_dphi_subtracted_0_20.GetNbinsX()
    h_lambda_total_integral_0_20 = 0
    h_lambda_near_integral_0_20 = 0
    h_lambda_away_integral_0_20 = 0
    h_lambda_ue_integral_0_20 = v2_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

    h_lambda_total_integral_error_0_20 = 0
    h_lambda_near_integral_error_0_20 = 0
    h_lambda_away_integral_error_0_20 = 0
    h_lambda_ue_integral_error_0_20 = v2_fit_0_20.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_20_halved = h_lambda_ue_integral_error_0_20/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_0_20 += h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num)
        h_lambda_total_integral_error_0_20 += (h_lambda_dphi_subtracted_0_20.GetBinError(bin_num))**2
        part = h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num) - v2_fit_0_20.Eval(h_lambda_dphi_subtracted_0_20.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_0_20 += part
            h_lambda_near_integral_error_0_20 += (h_lambda_dphi_subtracted_0_20.GetBinError(bin_num))**2
        else:
            h_lambda_away_integral_0_20 += part
            h_lambda_away_integral_error_0_20 += (h_lambda_dphi_subtracted_0_20.GetBinError(bin_num))**2
        
    h_lambda_total_integral_error_0_20 = math.sqrt(h_lambda_total_integral_error_0_20)
    h_lambda_near_integral_error_0_20 = math.sqrt(h_lambda_near_integral_error_0_20)
    h_lambda_near_integral_error_0_20 = math.sqrt(h_lambda_near_integral_error_0_20**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_away_integral_error_0_20 = math.sqrt(h_lambda_away_integral_error_0_20)
    h_lambda_away_integral_error_0_20 = math.sqrt(h_lambda_away_integral_error_0_20**2 + h_lambda_ue_integral_error_0_20_halved**2)

else: 

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
    h_lambda_ue_integral_error_0_20 = ue_avg_error_0_20*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_20_halved = h_lambda_ue_integral_error_0_20/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_0_20 += h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num)
        h_lambda_total_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2
        part = h_lambda_dphi_subtracted_0_20.GetBinContent(bin_num) - ue_avg_0_20
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_0_20 += part
            h_lambda_near_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2
        else:
            h_lambda_away_integral_0_20 += part
            h_lambda_away_integral_error_0_20 += h_lambda_dphi_subtracted_0_20.GetBinError(bin_num)**2

    h_lambda_total_integral_error_0_20 = math.sqrt(h_lambda_total_integral_error_0_20)
    h_lambda_near_integral_error_0_20 = math.sqrt(h_lambda_near_integral_error_0_20)
    h_lambda_near_integral_error_0_20 = math.sqrt(h_lambda_near_integral_error_0_20**2 + h_lambda_ue_integral_error_0_20_halved**2)
    h_lambda_away_integral_error_0_20 = math.sqrt(h_lambda_away_integral_error_0_20)
    h_lambda_away_integral_error_0_20 = math.sqrt(h_lambda_away_integral_error_0_20**2 + h_lambda_ue_integral_error_0_20_halved**2)


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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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
elif USE_FIT:

    # fitting to v2 component + four gaussians
    hh_fit_string_0_20 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_fit_string_0_20 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    hh_fit_function_0_20 = rt.TF1("hh_fit_function_0_20", hh_fit_string_0_20, -2, 6)
    hh_fit_function_0_20.SetNpx(1000)
    if PT_MODE == 0:
        hh_fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][0])
        hh_fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][2])
        hh_fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][0])
        hh_fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][2])
        hh_fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_fit_function_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][0])
        hh_fit_function_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][2])
        hh_fit_function_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][3])
    # setting parameters for first gaussian (centered at 0)
    hh_fit_function_0_20.SetParLimits(3, 0.02, 1)
    hh_fit_function_0_20.SetParameter(3, 0.16)
    hh_fit_function_0_20.FixParameter(4, 0)
    hh_fit_function_0_20.SetParLimits(5, 0.1, 1)
    hh_fit_function_0_20.SetParameter(5, 0.3)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    hh_fit_function_0_20.SetParLimits(6, 0, 1)
    hh_fit_function_0_20.SetParameter(6, 0)
    hh_fit_function_0_20.FixParameter(7, 2*rt.TMath.Pi())
    hh_fit_function_0_20.SetParLimits(8, 0.1, 2)
    hh_fit_function_0_20.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    hh_fit_function_0_20.SetParLimits(9, 0.02, 1)
    hh_fit_function_0_20.SetParameter(9, 0.3)
    hh_fit_function_0_20.FixParameter(10, rt.TMath.Pi())
    hh_fit_function_0_20.SetParLimits(11, 0.3, 2)
    hh_fit_function_0_20.SetParameter(11, 0.7)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    hh_fit_function_0_20.SetParLimits(12, 0, 1)
    hh_fit_function_0_20.SetParameter(12, 0)
    hh_fit_function_0_20.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    hh_fit_function_0_20.SetParLimits(14, 0.1, 3)
    hh_fit_function_0_20.SetParameter(14, 0.1)


    # setting parameters for constant background
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

    hh_fit_function_0_20.SetParameter(12, hh_ue_avg_0_20)


    h_h_dphi_with_fit_0_20 = h_h_dphi_0_20.Clone("h_h_dphi_with_fit_0_20")
    h_h_dphi_with_fit_0_20.Fit(hh_fit_function_0_20, "R")

    hh_ue_avg_fit_0_20 = rt.TF1("hh_ue_avg_fit_0_20", "pol0", -2, 6)
    hh_ue_avg_fit_0_20.SetParameter(0, hh_ue_avg_0_20)
    hh_ue_avg_fit_0_20.SetParError(0, hh_ue_avg_error_0_20)

elif USE_V2:
    hh_ue_avg_0_20 = (h_h_dphi_0_20.GetBinContent(1) 
                   + h_h_dphi_0_20.GetBinContent(2)
                   + h_h_dphi_0_20.GetBinContent(7)
                   + h_h_dphi_0_20.GetBinContent(8)
                   + h_h_dphi_0_20.GetBinContent(9)
                   + h_h_dphi_0_20.GetBinContent(16))/6

    hh_v2_fit_0_20 = rt.TF1("hh_v2_fit_0_20", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    hh_ue_avg_0_20 = (h_h_dphi_0_20.GetBinContent(1) 
                   + h_h_dphi_0_20.GetBinContent(2)
                   + h_h_dphi_0_20.GetBinContent(7)
                   + h_h_dphi_0_20.GetBinContent(8)
                   + h_h_dphi_0_20.GetBinContent(9)
                   + h_h_dphi_0_20.GetBinContent(16))/6

    hh_v2_fit_0_20 = rt.TF1("hh_v2_fit_0_20", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    hh_von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    hh_von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
    hh_von_fit_0_20 = rt.TF1("hh_von_fit_0_20", hh_von_fit_string, -2, 6)

    if PT_MODE == 0:
        hh_von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][0])
        hh_von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][2])
        hh_von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][0])
        hh_von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][2])
        hh_von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_von_fit_0_20.FixParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][0])
        hh_von_fit_0_20.FixParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][2])
        hh_von_fit_0_20.FixParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][3])

    hh_von_fit_0_20.SetParLimits(3, 0, 1)
    hh_von_fit_0_20.SetParameter(3, 0.02)
    hh_von_fit_0_20.SetParLimits(4, 0, 100)
    hh_von_fit_0_20.SetParameter(4, 1)
    hh_von_fit_0_20.SetParLimits(5, 0, 1)
    hh_von_fit_0_20.SetParameter(5, 0.01)
    hh_von_fit_0_20.SetParLimits(6, 0, 100)
    hh_von_fit_0_20.SetParameter(6, 1)

    h_h_dphi_with_fit_0_20 = h_h_dphi_0_20.Clone("h_h_dphi_with_fit_0_20")
    hh_fit_result_0_20 = h_h_dphi_with_fit_0_20.Fit(hh_von_fit_0_20, "R")

    if PT_MODE == 0:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_0_20.SetParameter(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_0_20.SetParError(0, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_0_20.SetParameter(1, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_0_20.SetParameter(2, v2_fitpars_dict["h_h_cent_0_20_trigger_4_8_assoc_25_4"][3])

else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_h_dphi_0_20_zeroed = h_h_dphi_0_20.Clone("h_h_dphi_0_20_zeroed")
    h_h_dphi_0_20_zeroed.Add(hh_ue_avg_fit_0_20, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_h_total_integral_error_0_20 = hh_fit_function_0_20.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_0_20 = hh_fit_function_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_0_20 = hh_ue_avg_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_0_20 = hh_ue_avg_error_0_20*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_20_halved = h_h_ue_integral_error_0_20/2
    h_h_near_integral_error_0_20 = math.sqrt((hh_fit_function_0_20.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_near_integral_0_20 = hh_fit_function_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_0_20 = math.sqrt((hh_fit_function_0_20.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_away_integral_0_20 = hh_fit_function_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_h_total_integral_error_0_20 = hh_von_fit_0_20.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_0_20 = hh_von_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_0_20 = hh_v2_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_0_20 = hh_v2_fit_0_20.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_20_halved = h_h_ue_integral_error_0_20/2
    h_h_near_integral_error_0_20 = math.sqrt((hh_von_fit_0_20.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_near_integral_0_20 = hh_von_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_0_20.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_0_20 = math.sqrt((hh_von_fit_0_20.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_away_integral_0_20 = hh_von_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_0_20.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_dphi_0_20_zeroed = h_h_dphi_0_20.Clone("h_h_dphi_0_20_zeroed")
    h_h_dphi_0_20_zeroed.Add(hh_v2_fit_0_20, -1)

elif USE_V2:

    DPHI_BINS = h_h_dphi_0_20.GetNbinsX()
    h_h_total_integral_0_20 = 0
    h_h_near_integral_0_20 = 0
    h_h_away_integral_0_20 = 0
    h_h_ue_integral_0_20 = hh_v2_fit_0_20.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_error_0_20 = 0
    h_h_near_integral_error_0_20 = 0
    h_h_away_integral_error_0_20 = 0
    h_h_ue_integral_error_0_20 = hh_v2_fit_0_20.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_20_halved = h_h_ue_integral_error_0_20/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_0_20 += h_h_dphi_0_20.GetBinContent(bin_num)
        h_h_total_integral_error_0_20 += (h_h_dphi_0_20.GetBinError(bin_num))**2
        part = h_h_dphi_0_20.GetBinContent(bin_num) - hh_v2_fit_0_20.Eval(h_h_dphi_0_20.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_0_20 += part
            h_h_near_integral_error_0_20 += (h_h_dphi_0_20.GetBinError(bin_num))**2
        else:
            h_h_away_integral_0_20 += part
            h_h_away_integral_error_0_20 += (h_h_dphi_0_20.GetBinError(bin_num))**2
    h_h_total_integral_error_0_20 = math.sqrt(h_h_total_integral_error_0_20)
    h_h_near_integral_error_0_20 = math.sqrt(h_h_near_integral_error_0_20)
    h_h_near_integral_error_0_20 = math.sqrt(h_h_near_integral_error_0_20**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_away_integral_error_0_20 = math.sqrt(h_h_away_integral_error_0_20)
    h_h_away_integral_error_0_20 = math.sqrt(h_h_away_integral_error_0_20**2 + h_h_ue_integral_error_0_20_halved**2)

else: 

    h_h_dphi_0_20_zeroed = h_h_dphi_0_20.Clone("h_h_dphi_0_20_zeroed")
    h_h_dphi_0_20_zeroed.Add(hh_ue_line_0_20, -1)

    DPHI_BINS = h_h_dphi_0_20.GetNbinsX()
    h_h_total_integral_0_20 = 0
    h_h_near_integral_0_20 = 0
    h_h_away_integral_0_20 = 0
    h_h_ue_integral_0_20 = hh_ue_avg_0_20*DPHI_BINS
    h_h_total_integral_error_0_20 = 0
    h_h_near_integral_error_0_20 = 0
    h_h_away_integral_error_0_20 = 0
    h_h_ue_integral_error_0_20 = hh_ue_avg_error_0_20*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_20_halved = h_h_ue_integral_error_0_20/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_0_20 += h_h_dphi_0_20.GetBinContent(bin_num)
        h_h_total_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
        part = h_h_dphi_0_20.GetBinContent(bin_num) - hh_ue_avg_0_20
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_0_20 += part
            h_h_near_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
        else:
            h_h_away_integral_0_20 += part
            h_h_away_integral_error_0_20 += h_h_dphi_0_20.GetBinError(bin_num)**2
    h_h_total_integral_error_0_20 = math.sqrt(h_h_total_integral_error_0_20)
    h_h_near_integral_error_0_20 = math.sqrt(h_h_near_integral_error_0_20)
    h_h_near_integral_error_0_20 = math.sqrt(h_h_near_integral_error_0_20**2 + h_h_ue_integral_error_0_20_halved**2)
    h_h_away_integral_error_0_20 = math.sqrt(h_h_away_integral_error_0_20)
    h_h_away_integral_error_0_20 = math.sqrt(h_h_away_integral_error_0_20**2 + h_h_ue_integral_error_0_20_halved**2)


output_file.cd()
h_lambda_2d_subtracted_0_20.Write()

h_lambda_dphi_subtracted_0_20.Write()
h_h_dphi_0_20.Write()



if USE_FIT:
    fit_function_0_20.Write()
    v2_fit_0_20.Write()
    hh_fit_function_0_20.Write()
    hh_ue_avg_fit_0_20.Write()
    h_lambda_dphi_subtracted_with_fit_0_20.Write()
    h_h_dphi_with_fit_0_20.Write()
    h_lambda_dphi_subtracted_0_20_zeroed.Write()
    h_h_dphi_0_20_zeroed.Write()
elif USE_V2:
    v2_fit_0_20.Write()
    hh_v2_fit_0_20.Write()
elif USE_VON:
    von_fit_0_20.Write()
    hh_von_fit_0_20.Write()
    v2_fit_0_20.Write()
    hh_v2_fit_0_20.Write()
    h_lambda_dphi_subtracted_with_fit_0_20.Write()
    h_h_dphi_with_fit_0_20.Write()
    h_lambda_dphi_subtracted_0_20_zeroed.Write()
    h_h_dphi_0_20_zeroed.Write()
else:
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

    h_lambda_dphi_subtracted_0_20_zeroed.Write()
    h_h_dphi_0_20_zeroed.Write()


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
############################################# 20-50 CENTRALITY ##############################################
############################################################################################################
############################################################################################################

if IS_NORMAL_PID:
    input_file_20_50 = rt.TFile("../online/v0/output/v0_normal_20_50.root")
elif IS_WIDE_PID:
    input_file_20_50 = rt.TFile("../online/v0/output/v0_wide_20_50.root")
elif IS_NARROW_PID:
    input_file_20_50 = rt.TFile("../online/v0/output/v0_narrow_20_50.root")
elif IS_TOF_PID:
    input_file_20_50 = rt.TFile("../online/v0/output/v0_tof_20_50.root")
else:
    print("ERROR: No PID cut specified")
    exit(1)

input_list_20_50 = input_file_20_50.Get("h-lambda")
input_file_20_50.Close()

trig_dist_20_50 = input_list_20_50.FindObject("fTriggerDistEff_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fTriggerDistEff")
lambda_dist_20_50 = input_list_20_50.FindObject("fTriggeredLambdaDist")


h_h_20_50 = input_list_20_50.FindObject("fDphiHHEff_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHHEff")
h_h_mixed_20_50 = input_list_20_50.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHHMixed")

h_lambda_20_50 = input_list_20_50.FindObject("fDphiHLambdaEff_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHLambdaEff")
h_lambda_mixed_20_50 = input_list_20_50.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_20_50.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_20_50.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_20_50.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_mixed_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_20_50.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_20_50.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_20_50.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


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

bin_1 = lambda_mass_dist_20_50.FindBin(1.102)
bin_2 = lambda_mass_dist_20_50.FindBin(1.132)
point_one = [1.102, lambda_mass_dist_20_50.GetBinContent(bin_1)]
point_two = [1.132, lambda_mass_dist_20_50.GetBinContent(bin_2)]
bg_starting_params_20_50 = get_straight_line(point_one, point_two)
lambda_mass_fit_20_50 = rt.TF1("lambda_mass_fit_20_50", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)", 1.102, 1.127)
lambda_mass_fit_20_50.SetNpx(1000)
lambda_mass_fit_20_50.SetParameter(0, 1.74e03)
lambda_mass_fit_20_50.SetParameter(1, 1.116)
lambda_mass_fit_20_50.SetParameter(2, 1.30576e-03 )
lambda_mass_fit_20_50.SetParameter(3, 1.804166e-03)
lambda_mass_fit_20_50.SetParameter(4, bg_starting_params_20_50[0])
lambda_mass_fit_20_50.SetParameter(5, bg_starting_params_20_50[1])
lambda_mass_dist_fit_20_50 = lambda_mass_dist_20_50.Clone("lambda_mass_dist_fit_20_50")


lambda_mass_bg_fit_20_50 = rt.TF1("bg_fit_20_50", "pol1", 1.096, 1.136)
lambda_mass_bg_fit_20_50.SetNpx(1000)
lambda_mass_bg_fit_20_50.SetParameter(0, bg_starting_params_20_50[0])
lambda_mass_bg_fit_20_50.SetParameter(1, bg_starting_params_20_50[1])
lambda_mass_dist_fit_20_50.Fit(lambda_mass_fit_20_50, "RS")

voigt_fit_20_50 = rt.TF1("voigt_fit_20_50", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.102, 1.127)
voigt_fit_20_50.SetNpx(1000)
voigt_fit_20_50.SetParameter(0, lambda_mass_fit_20_50.GetParameter(0))
voigt_fit_20_50.SetParameter(1, lambda_mass_fit_20_50.GetParameter(1))
voigt_fit_20_50.SetParameter(2, lambda_mass_fit_20_50.GetParameter(2))
voigt_fit_20_50.SetParameter(3, lambda_mass_fit_20_50.GetParameter(3))


residual_20_50 = lambda_mass_dist_20_50.Clone("residual_20_50")
residual_20_50.GetXaxis().SetRangeUser(1.094, 1.142)
residual_20_50.Add(lambda_mass_bg_fit_20_50, -1)


RSB_region_20_50 = rt.TBox(RSB_MIN, 0, RSB_MAX, lambda_mass_dist_20_50.GetMaximum()*1.055)
RSB_region_20_50.SetLineColor(rt.kRed)
RSB_region_20_50.SetFillColor(rt.kRed)
RSB_region_20_50.SetFillStyle(3003)


RSB_min_line_20_50 = rt.TLine(RSB_MIN, 0, RSB_MIN, lambda_mass_dist_20_50.GetMaximum()*1.05)
RSB_min_line_20_50.SetLineColor(rt.kRed)
RSB_min_line_20_50.SetLineWidth(2)
RSB_min_line_20_50.SetLineStyle(2)


RSB_max_line_20_50 = rt.TLine(RSB_MAX, 0, RSB_MAX, lambda_mass_dist_20_50.GetMaximum()*1.05)
RSB_max_line_20_50.SetLineColor(rt.kRed)
RSB_max_line_20_50.SetLineWidth(2)
RSB_max_line_20_50.SetLineStyle(2)

left_signal_bin_20_50 = lambda_mass_dist_20_50.FindBin(SIG_MIN)
right_signal_bin_20_50 = lambda_mass_dist_20_50.FindBin(SIG_MAX)

full_region_bin_1 = residual_20_50.FindBin(FULL_REGION_MIN)
full_region_bin_2 = residual_20_50.FindBin(FULL_REGION_MAX)

lambda_bg_20_50 = 0
lambda_total_20_50 = 0
for bin_num in range(left_signal_bin_20_50, right_signal_bin_20_50 + 1):
    lambda_total_20_50 += lambda_mass_dist_20_50.GetBinContent(bin_num)
    lambda_bg_20_50 += lambda_mass_bg_fit_20_50.Eval(lambda_mass_dist_20_50.GetBinCenter(bin_num))


lambda_signal_20_50 = lambda_total_20_50 - lambda_bg_20_50
lambda_signal_total_ratio_20_50 = lambda_signal_20_50/lambda_total_20_50

signal_scale_20_50 = residual_20_50.Integral(full_region_bin_1, full_region_bin_2)/residual_20_50.Integral(left_signal_bin_20_50, right_signal_bin_20_50)

output_file.cd()
RSB_region_20_50.Write("RSB_region_20_50")
RSB_min_line_20_50.Write("RSB_min_line_20_50")
RSB_max_line_20_50.Write("RSB_max_line_20_50")

lambda_mass_dist_20_50.Write()
lambda_mass_dist_fit_20_50.Write()
lambda_mass_fit_20_50.Write()
lambda_mass_bg_fit_20_50.Write()
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
h_lambda_2d_mixcor_rsb_20_50 = make_mixed_corrections(h_lambda_20_50, h_lambda_mixed_20_50, RSB_MIN, RSB_MAX)

h_h_2d_mixcor_20_50 = make_mixed_corrections(h_h_20_50, h_h_mixed_20_50, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_20_50.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_20_50.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_h_2d_mixcor_20_50.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)

h_lambda_2d_mixcor_sig_20_50.SetName("h_lambda_2d_mixcor_sig_20_50")
h_lambda_2d_mixcor_rsb_20_50.SetName("h_lambda_2d_mixcor_rsb_20_50")
h_h_2d_mixcor_20_50.SetName("h_h_2d_mixcor_20_50")




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_20_50.Scale(1.0/num_trigs_20_50)
h_lambda_2d_mixcor_rsb_20_50.Scale(1.0/num_trigs_20_50)
h_h_2d_mixcor_20_50.Scale(1.0/num_trigs_20_50)

output_file.cd()
h_lambda_2d_nomixcor_20_50.Write()
h_lambda_mixed_2d_20_50.Write()
h_h_2d_nomixcor_20_50.Write()
h_h_mixed_2d_20_50.Write()
h_lambda_2d_mixcor_sig_20_50.Write()
h_lambda_2d_mixcor_rsb_20_50.Write()
h_h_2d_mixcor_20_50.Write()

# SIDEBAND SUBTRACTION SECTION
if DO_SIDEBAND_SUBTRACTION:

    # Normalize side band to 1
    h_lambda_2d_mixcor_rsb_20_50.Scale(1/h_lambda_2d_mixcor_rsb_20_50.Integral())


    # using RSB for sideband subtraction
    h_lambda_2d_subtracted_20_50 = h_lambda_2d_mixcor_sig_20_50.Clone("h_lambda_2d_subtracted_20_50")
    bg_integral_20_50 = (1 - lambda_signal_total_ratio_20_50)*h_lambda_2d_subtracted_20_50.Integral()
    h_lambda_2d_subtracted_20_50.Add(h_lambda_2d_mixcor_rsb_20_50, -bg_integral_20_50)

    # save the RSB deltaphi distribution
    h_lambda_dphi_rsb_20_50 = h_lambda_2d_mixcor_rsb_20_50.ProjectionY("h_lambda_dphi_rsb_20_50")
    output_file.cd()
    h_lambda_dphi_rsb_20_50.Write()

else:
    h_lambda_2d_subtracted_20_50 = h_lambda_2d_mixcor_sig_20_50.Clone("h_lambda_2d_subtracted_20_50")


# Correct for signal scale
h_lambda_2d_subtracted_20_50.Scale(signal_scale_20_50)
# Correct for PID loss
h_lambda_2d_subtracted_20_50.Scale(PID_CORRECTION)
# Correct for branching ratio
h_lambda_2d_subtracted_20_50.Scale(1/BRANCHING_RATIO)

# Correct for two-track efficiency
h_lambda_2d_subtracted_20_50 = apply_two_track_correction(h_lambda_2d_subtracted_20_50, TWOTRACK_TEMPLATE)



# INTEGRAL AND RATIO SECTION
h_lambda_dphi_subtracted_20_50 = h_lambda_2d_subtracted_20_50.ProjectionY("h_lambda_dphi_subtracted_20_50")

# this is completely insane but i need to scale the PID stuff exactly prior to fitting
if IS_WIDE_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_20_50.Scale(0.9270985488699187)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_20_50.Scale(0.9275249876914017)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_20_50.Scale(0.9298246568065354)
elif IS_NARROW_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_20_50.Scale(1.1844748339350628)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_20_50.Scale(1.1625641331436003)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_20_50.Scale(1.1977363386098983)
elif IS_TOF_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_20_50.Scale(0.9911803993835522)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_20_50.Scale(1.7238151003786188)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_20_50.Scale(0.8182694751096258)


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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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

elif USE_FIT:

    # fitting to v2 component + four gaussians
    fit_string_20_50 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    fit_string_20_50 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    fit_function_20_50 = rt.TF1("fit_function_20_50", fit_string_20_50, -2, 6)
    fit_function_20_50.SetNpx(1000)
    if PT_MODE == 0:
        fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][0])
        fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][2])
        fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][0])
        fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][2])
        fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][0])
        fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][2])
        fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    fit_function_20_50.SetParLimits(3, 5e-4, 0.6)
    fit_function_20_50.SetParameter(3, 0.009)
    fit_function_20_50.FixParameter(4, 0)
    fit_function_20_50.SetParLimits(5, 0.1, 2)
    fit_function_20_50.SetParameter(5, 0.8)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    fit_function_20_50.SetParLimits(6, 0,  0.6)
    fit_function_20_50.SetParameter(6, 0)
    fit_function_20_50.FixParameter(7, 2*rt.TMath.Pi())
    fit_function_20_50.SetParLimits(8, 0.1, 2)
    fit_function_20_50.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    fit_function_20_50.SetParLimits(9, 5e-4, 0.6)
    fit_function_20_50.SetParameter(9, 0.003)
    fit_function_20_50.FixParameter(10, rt.TMath.Pi())
    fit_function_20_50.SetParLimits(11, 0.1, 2)
    fit_function_20_50.SetParameter(11, 1)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    fit_function_20_50.SetParLimits(12, 0, 0.6)
    fit_function_20_50.SetParameter(12, 0)
    fit_function_20_50.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    fit_function_20_50.SetParLimits(14, 0.1, 2)
    fit_function_20_50.SetParameter(14, 0.1)

    # setting parameters for constant background
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

    fit_function_20_50.SetParameter(12, ue_avg_20_50)


    h_lambda_dphi_subtracted_with_fit_20_50 = h_lambda_dphi_subtracted_20_50.Clone("h_lambda_dphi_subtracted_with_fit_20_50")
    fit_result_20_50 = h_lambda_dphi_subtracted_with_fit_20_50.Fit(fit_function_20_50, "RS")

    ue_avg_fit_20_50 = rt.TF1("ue_avg_fit_20_50", "pol0", -2, 6)
    ue_avg_fit_20_50.SetParameter(0, ue_avg_20_50)
    ue_avg_fit_20_50.SetParError(0, ue_avg_error_20_50)

elif USE_V2:
    ue_avg_20_50 = (h_lambda_dphi_subtracted_20_50.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(2)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(7)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(8)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(9)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(16))/6

    v2_fit_20_50 = rt.TF1("v2_fit_20_50", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    ue_avg_20_50 = (h_lambda_dphi_subtracted_20_50.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(2)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(7)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(8)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(9)
                   + h_lambda_dphi_subtracted_20_50.GetBinContent(16))/6

    v2_fit_20_50 = rt.TF1("v2_fit_20_50", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"

    von_fit_20_50 = rt.TF1("von_fit_20_50", von_fit_string, -2, 6)

    scale = 0.8
    if PT_MODE == 0:
        von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][0])
        von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][2])
        von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][3])
        # von_fit_20_50.SetParLimits(3, scale*0.015/3, scale*0.015*3)
        von_fit_20_50.SetParameter(3, scale*0.015)
        # von_fit_20_50.SetParLimits(4, (7.6/scale)/3, (7.6/scale)*3)
        von_fit_20_50.SetParameter(4, 7.6/scale)
        # von_fit_20_50.SetParLimits(5, scale*0.014/3, scale*0.014*3)
        von_fit_20_50.SetParameter(5, scale*0.014)
        # von_fit_20_50.SetParLimits(6, (2.2/scale)/3, (2.2/scale)*3)
        von_fit_20_50.SetParameter(6, 2.2/scale)
    elif PT_MODE == 1:
        von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][0])
        von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][2])
        von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][3])
        # von_fit_20_50.SetParLimits(3, scale*0.02/6, scale*0.02*6)
        von_fit_20_50.SetParameter(3, scale*0.02)
        # von_fit_20_50.SetParLimits(4, (4.6/scale)/3, (4.6/scale)*3)
        von_fit_20_50.SetParameter(4, 4.6/scale)
        # von_fit_20_50.SetParLimits(5, scale*0.057/6, scale*0.057*6)
        von_fit_20_50.SetParameter(5, scale*0.057)
        # von_fit_20_50.SetParLimits(6, (1/scale), 3*(1/scale))
        von_fit_20_50.SetParameter(6, 1/scale)
    elif PT_MODE == 2:
        von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][0])
        von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][2])
        von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][3])
        # von_fit_20_50.SetParLimits(3, scale*0.009/3, scale*0.009*3)
        von_fit_20_50.SetParameter(3, scale*0.009)
        # von_fit_20_50.SetParLimits(4, (8.88/scale)/3, (8.88/scale)*3)
        von_fit_20_50.SetParameter(4, 8.88/scale)
        # von_fit_20_50.SetParLimits(5, scale*0.0076/6, scale*0.0076*6)
        von_fit_20_50.SetParameter(5, scale*0.0076)
        # von_fit_20_50.SetParLimits(6, (2.55/scale)/3, (2.55/scale)*3)
        von_fit_20_50.SetParameter(6, 2.55/scale)

    h_lambda_dphi_subtracted_with_fit_20_50 = h_lambda_dphi_subtracted_20_50.Clone("h_lambda_dphi_subtracted_with_fit_20_50")

    if PT_MODE == 1: # low pt has away-side width too large for von mises fit, only works if we symmetrize the peaks
        fit_result_20_50 = h_lambda_dphi_subtracted_with_fit_20_50.Fit(von_fit_20_50, "R", "", -rt.TMath.Pi()/2, rt.TMath.Pi() - EPSILON)
    else:
        fit_result_20_50 = h_lambda_dphi_subtracted_with_fit_20_50.Fit(von_fit_20_50, "R")


    if PT_MODE == 0:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][0])
        v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][1])
        v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][2])
        v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_lambda_cent_20_50_trigger_4_8_assoc_25_4"][3])


else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_lambda_dphi_subtracted_20_50_zeroed = h_lambda_dphi_subtracted_20_50.Clone("h_lambda_dphi_subtracted_20_50_zeroed")
    h_lambda_dphi_subtracted_20_50_zeroed.Add(ue_avg_fit_20_50, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_lambda_total_integral_error_20_50 = fit_function_20_50.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_20_50 = fit_function_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_20_50 = ue_avg_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_20_50 = ue_avg_error_20_50*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_20_50_halved = h_lambda_ue_integral_error_20_50/2
    h_lambda_near_integral_error_20_50 = math.sqrt((fit_function_20_50.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_near_integral_20_50 = fit_function_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_20_50 = math.sqrt((fit_function_20_50.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_away_integral_20_50 = fit_function_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_lambda_total_integral_error_20_50 = von_fit_20_50.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_20_50 = von_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_20_50 = v2_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_20_50 = v2_fit_20_50.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_20_50_halved = h_lambda_ue_integral_error_20_50/2
    h_lambda_near_integral_error_20_50 = math.sqrt((von_fit_20_50.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_near_integral_20_50 = von_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_20_50 = math.sqrt((von_fit_20_50.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_away_integral_20_50 = von_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_dphi_subtracted_20_50_zeroed = h_lambda_dphi_subtracted_20_50.Clone("h_lambda_dphi_subtracted_20_50_zeroed")
    h_lambda_dphi_subtracted_20_50_zeroed.Add(v2_fit_20_50, -1)

elif USE_V2:

    DPHI_BINS = h_lambda_dphi_subtracted_20_50.GetNbinsX()
    h_lambda_total_integral_20_50 = 0
    h_lambda_near_integral_20_50 = 0
    h_lambda_away_integral_20_50 = 0
    h_lambda_ue_integral_20_50 = v2_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

    h_lambda_total_integral_error_20_50 = 0
    h_lambda_near_integral_error_20_50 = 0
    h_lambda_away_integral_error_20_50 = 0
    h_lambda_ue_integral_error_20_50 = v2_fit_20_50.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_20_50_halved = h_lambda_ue_integral_error_20_50/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_20_50 += h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num)
        h_lambda_total_integral_error_20_50 += (h_lambda_dphi_subtracted_20_50.GetBinError(bin_num))**2
        part = h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num) - v2_fit_20_50.Eval(h_lambda_dphi_subtracted_20_50.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_20_50 += part
            h_lambda_near_integral_error_20_50 += (h_lambda_dphi_subtracted_20_50.GetBinError(bin_num))**2
        else:
            h_lambda_away_integral_20_50 += part
            h_lambda_away_integral_error_20_50 += (h_lambda_dphi_subtracted_20_50.GetBinError(bin_num))**2
        
    h_lambda_total_integral_error_20_50 = math.sqrt(h_lambda_total_integral_error_20_50)
    h_lambda_near_integral_error_20_50 = math.sqrt(h_lambda_near_integral_error_20_50)
    h_lambda_near_integral_error_20_50 = math.sqrt(h_lambda_near_integral_error_20_50**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_away_integral_error_20_50 = math.sqrt(h_lambda_away_integral_error_20_50)
    h_lambda_away_integral_error_20_50 = math.sqrt(h_lambda_away_integral_error_20_50**2 + h_lambda_ue_integral_error_20_50_halved**2)

else: 

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
    h_lambda_ue_integral_error_20_50 = ue_avg_error_20_50*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_20_50_halved = h_lambda_ue_integral_error_20_50/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_20_50 += h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num)
        h_lambda_total_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2
        part = h_lambda_dphi_subtracted_20_50.GetBinContent(bin_num) - ue_avg_20_50
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_20_50 += part
            h_lambda_near_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2
        else:
            h_lambda_away_integral_20_50 += part
            h_lambda_away_integral_error_20_50 += h_lambda_dphi_subtracted_20_50.GetBinError(bin_num)**2

    h_lambda_total_integral_error_20_50 = math.sqrt(h_lambda_total_integral_error_20_50)
    h_lambda_near_integral_error_20_50 = math.sqrt(h_lambda_near_integral_error_20_50)
    h_lambda_near_integral_error_20_50 = math.sqrt(h_lambda_near_integral_error_20_50**2 + h_lambda_ue_integral_error_20_50_halved**2)
    h_lambda_away_integral_error_20_50 = math.sqrt(h_lambda_away_integral_error_20_50)
    h_lambda_away_integral_error_20_50 = math.sqrt(h_lambda_away_integral_error_20_50**2 + h_lambda_ue_integral_error_20_50_halved**2)


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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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
elif USE_FIT:

    # fitting to v2 component + four gaussians
    hh_fit_string_20_50 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_fit_string_20_50 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    hh_fit_function_20_50 = rt.TF1("hh_fit_function_20_50", hh_fit_string_20_50, -2, 6)
    hh_fit_function_20_50.SetNpx(1000)
    if PT_MODE == 0:
        hh_fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][0])
        hh_fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][2])
        hh_fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][0])
        hh_fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][2])
        hh_fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_fit_function_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][0])
        hh_fit_function_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][2])
        hh_fit_function_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    hh_fit_function_20_50.SetParLimits(3, 0.02, 1)
    hh_fit_function_20_50.SetParameter(3, 0.16)
    hh_fit_function_20_50.FixParameter(4, 0)
    hh_fit_function_20_50.SetParLimits(5, 0.1, 1)
    hh_fit_function_20_50.SetParameter(5, 0.3)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    hh_fit_function_20_50.SetParLimits(6, 0, 1)
    hh_fit_function_20_50.SetParameter(6, 0)
    hh_fit_function_20_50.FixParameter(7, 2*rt.TMath.Pi())
    hh_fit_function_20_50.SetParLimits(8, 0.1, 2)
    hh_fit_function_20_50.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    hh_fit_function_20_50.SetParLimits(9, 0.02, 1)
    hh_fit_function_20_50.SetParameter(9, 0.3)
    hh_fit_function_20_50.FixParameter(10, rt.TMath.Pi())
    hh_fit_function_20_50.SetParLimits(11, 0.3, 2)
    hh_fit_function_20_50.SetParameter(11, 0.7)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    hh_fit_function_20_50.SetParLimits(12, 0, 1)
    hh_fit_function_20_50.SetParameter(12, 0)
    hh_fit_function_20_50.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    hh_fit_function_20_50.SetParLimits(14, 0.1, 3)
    hh_fit_function_20_50.SetParameter(14, 0.1)


    # setting parameters for constant background
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

    hh_fit_function_20_50.SetParameter(12, hh_ue_avg_20_50)


    h_h_dphi_with_fit_20_50 = h_h_dphi_20_50.Clone("h_h_dphi_with_fit_20_50")
    h_h_dphi_with_fit_20_50.Fit(hh_fit_function_20_50, "R")

    hh_ue_avg_fit_20_50 = rt.TF1("hh_ue_avg_fit_20_50", "pol0", -2, 6)
    hh_ue_avg_fit_20_50.SetParameter(0, hh_ue_avg_20_50)
    hh_ue_avg_fit_20_50.SetParError(0, hh_ue_avg_error_20_50)

elif USE_V2:
    hh_ue_avg_20_50 = (h_h_dphi_20_50.GetBinContent(1) 
                   + h_h_dphi_20_50.GetBinContent(2)
                   + h_h_dphi_20_50.GetBinContent(7)
                   + h_h_dphi_20_50.GetBinContent(8)
                   + h_h_dphi_20_50.GetBinContent(9)
                   + h_h_dphi_20_50.GetBinContent(16))/6

    hh_v2_fit_20_50 = rt.TF1("hh_v2_fit_20_50", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    hh_ue_avg_20_50 = (h_h_dphi_20_50.GetBinContent(1) 
                   + h_h_dphi_20_50.GetBinContent(2)
                   + h_h_dphi_20_50.GetBinContent(7)
                   + h_h_dphi_20_50.GetBinContent(8)
                   + h_h_dphi_20_50.GetBinContent(9)
                   + h_h_dphi_20_50.GetBinContent(16))/6

    hh_v2_fit_20_50 = rt.TF1("hh_v2_fit_20_50", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    hh_von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    hh_von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
    hh_von_fit_20_50 = rt.TF1("hh_von_fit_20_50", hh_von_fit_string, -2, 6)

    if PT_MODE == 0:
        hh_von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][0])
        hh_von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][2])
        hh_von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][0])
        hh_von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][2])
        hh_von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_von_fit_20_50.FixParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][0])
        hh_von_fit_20_50.FixParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][2])
        hh_von_fit_20_50.FixParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][3])

    hh_von_fit_20_50.SetParLimits(3, 0, 1)
    hh_von_fit_20_50.SetParameter(3, 0.02)
    hh_von_fit_20_50.SetParLimits(4, 0, 100)
    hh_von_fit_20_50.SetParameter(4, 1)
    hh_von_fit_20_50.SetParLimits(5, 0, 1)
    hh_von_fit_20_50.SetParameter(5, 0.01)
    hh_von_fit_20_50.SetParLimits(6, 0, 100)
    hh_von_fit_20_50.SetParameter(6, 1)

    h_h_dphi_with_fit_20_50 = h_h_dphi_20_50.Clone("h_h_dphi_with_fit_20_50")
    hh_fit_result_20_50 = h_h_dphi_with_fit_20_50.Fit(hh_von_fit_20_50, "R")

    if PT_MODE == 0:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_20_50.SetParameter(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_20_50.SetParError(0, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_20_50.SetParameter(1, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_20_50.SetParameter(2, v2_fitpars_dict["h_h_cent_20_50_trigger_4_8_assoc_25_4"][3])

else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_h_dphi_20_50_zeroed = h_h_dphi_20_50.Clone("h_h_dphi_20_50_zeroed")
    h_h_dphi_20_50_zeroed.Add(hh_ue_avg_fit_20_50, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_h_total_integral_error_20_50 = hh_fit_function_20_50.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_20_50 = hh_fit_function_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_20_50 = hh_ue_avg_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_20_50 = hh_ue_avg_error_20_50*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_20_50_halved = h_h_ue_integral_error_20_50/2
    h_h_near_integral_error_20_50 = math.sqrt((hh_fit_function_20_50.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_near_integral_20_50 = hh_fit_function_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_20_50 = math.sqrt((hh_fit_function_20_50.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_away_integral_20_50 = hh_fit_function_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_h_total_integral_error_20_50 = hh_von_fit_20_50.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_20_50 = hh_von_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_20_50 = hh_v2_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_20_50 = hh_v2_fit_20_50.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_20_50_halved = h_h_ue_integral_error_20_50/2
    h_h_near_integral_error_20_50 = math.sqrt((hh_von_fit_20_50.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_near_integral_20_50 = hh_von_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_20_50.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_20_50 = math.sqrt((hh_von_fit_20_50.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_away_integral_20_50 = hh_von_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_20_50.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_dphi_20_50_zeroed = h_h_dphi_20_50.Clone("h_h_dphi_20_50_zeroed")
    h_h_dphi_20_50_zeroed.Add(hh_v2_fit_20_50, -1)

elif USE_V2:

    DPHI_BINS = h_h_dphi_20_50.GetNbinsX()
    h_h_total_integral_20_50 = 0
    h_h_near_integral_20_50 = 0
    h_h_away_integral_20_50 = 0
    h_h_ue_integral_20_50 = hh_v2_fit_20_50.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_error_20_50 = 0
    h_h_near_integral_error_20_50 = 0
    h_h_away_integral_error_20_50 = 0
    h_h_ue_integral_error_20_50 = hh_v2_fit_20_50.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_20_50_halved = h_h_ue_integral_error_20_50/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_20_50 += h_h_dphi_20_50.GetBinContent(bin_num)
        h_h_total_integral_error_20_50 += (h_h_dphi_20_50.GetBinError(bin_num))**2
        part = h_h_dphi_20_50.GetBinContent(bin_num) - hh_v2_fit_20_50.Eval(h_h_dphi_20_50.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_20_50 += part
            h_h_near_integral_error_20_50 += (h_h_dphi_20_50.GetBinError(bin_num))**2
        else:
            h_h_away_integral_20_50 += part
            h_h_away_integral_error_20_50 += (h_h_dphi_20_50.GetBinError(bin_num))**2
    h_h_total_integral_error_20_50 = math.sqrt(h_h_total_integral_error_20_50)
    h_h_near_integral_error_20_50 = math.sqrt(h_h_near_integral_error_20_50)
    h_h_near_integral_error_20_50 = math.sqrt(h_h_near_integral_error_20_50**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_away_integral_error_20_50 = math.sqrt(h_h_away_integral_error_20_50)
    h_h_away_integral_error_20_50 = math.sqrt(h_h_away_integral_error_20_50**2 + h_h_ue_integral_error_20_50_halved**2)

else: 

    h_h_dphi_20_50_zeroed = h_h_dphi_20_50.Clone("h_h_dphi_20_50_zeroed")
    h_h_dphi_20_50_zeroed.Add(hh_ue_line_20_50, -1)

    DPHI_BINS = h_h_dphi_20_50.GetNbinsX()
    h_h_total_integral_20_50 = 0
    h_h_near_integral_20_50 = 0
    h_h_away_integral_20_50 = 0
    h_h_ue_integral_20_50 = hh_ue_avg_20_50*DPHI_BINS
    h_h_total_integral_error_20_50 = 0
    h_h_near_integral_error_20_50 = 0
    h_h_away_integral_error_20_50 = 0
    h_h_ue_integral_error_20_50 = hh_ue_avg_error_20_50*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_20_50_halved = h_h_ue_integral_error_20_50/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_20_50 += h_h_dphi_20_50.GetBinContent(bin_num)
        h_h_total_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
        part = h_h_dphi_20_50.GetBinContent(bin_num) - hh_ue_avg_20_50
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_20_50 += part
            h_h_near_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
        else:
            h_h_away_integral_20_50 += part
            h_h_away_integral_error_20_50 += h_h_dphi_20_50.GetBinError(bin_num)**2
    h_h_total_integral_error_20_50 = math.sqrt(h_h_total_integral_error_20_50)
    h_h_near_integral_error_20_50 = math.sqrt(h_h_near_integral_error_20_50)
    h_h_near_integral_error_20_50 = math.sqrt(h_h_near_integral_error_20_50**2 + h_h_ue_integral_error_20_50_halved**2)
    h_h_away_integral_error_20_50 = math.sqrt(h_h_away_integral_error_20_50)
    h_h_away_integral_error_20_50 = math.sqrt(h_h_away_integral_error_20_50**2 + h_h_ue_integral_error_20_50_halved**2)


output_file.cd()
h_lambda_2d_subtracted_20_50.Write()

h_lambda_dphi_subtracted_20_50.Write()
h_h_dphi_20_50.Write()



if USE_FIT:
    fit_function_20_50.Write()
    ue_avg_fit_20_50.Write()
    hh_fit_function_20_50.Write()
    hh_ue_avg_fit_20_50.Write()
    h_lambda_dphi_subtracted_with_fit_20_50.Write()
    h_h_dphi_with_fit_20_50.Write()
    h_lambda_dphi_subtracted_20_50_zeroed.Write()
    h_h_dphi_20_50_zeroed.Write()
elif USE_V2:
    v2_fit_20_50.Write()
    hh_v2_fit_20_50.Write()
elif USE_VON:
    von_fit_20_50.Write()
    hh_von_fit_20_50.Write()
    v2_fit_20_50.Write()
    hh_v2_fit_20_50.Write()
    h_lambda_dphi_subtracted_with_fit_20_50.Write()
    h_h_dphi_with_fit_20_50.Write()
    h_lambda_dphi_subtracted_20_50_zeroed.Write()
    h_h_dphi_20_50_zeroed.Write()
else:
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

    h_lambda_dphi_subtracted_20_50_zeroed.Write()
    h_h_dphi_20_50_zeroed.Write()


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
############################################# 50-80 CENTRALITY ##############################################
############################################################################################################
############################################################################################################

if IS_NORMAL_PID:
    input_file_50_80 = rt.TFile("../online/v0/output/v0_normal_50_80.root")
elif IS_WIDE_PID:
    input_file_50_80 = rt.TFile("../online/v0/output/v0_wide_50_80.root")
elif IS_NARROW_PID:
    input_file_50_80 = rt.TFile("../online/v0/output/v0_narrow_50_80.root")
elif IS_TOF_PID:
    input_file_50_80 = rt.TFile("../online/v0/output/v0_tof_50_80.root")
else:
    print("ERROR: No PID cut specified")
    exit(1)

input_list_50_80 = input_file_50_80.Get("h-lambda")
input_file_50_80.Close()

trig_dist_50_80 = input_list_50_80.FindObject("fTriggerDistEff_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fTriggerDistEff")
lambda_dist_50_80 = input_list_50_80.FindObject("fTriggeredLambdaDist")


h_h_50_80 = input_list_50_80.FindObject("fDphiHHEff_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHHEff")
h_h_mixed_50_80 = input_list_50_80.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHHMixed")

h_lambda_50_80 = input_list_50_80.FindObject("fDphiHLambdaEff_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHLambdaEff")
h_lambda_mixed_50_80 = input_list_50_80.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_50_80.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_50_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_50_80.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_mixed_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_50_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_50_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_50_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


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

bin_1 = lambda_mass_dist_50_80.FindBin(1.102)
bin_2 = lambda_mass_dist_50_80.FindBin(1.132)
point_one = [1.102, lambda_mass_dist_50_80.GetBinContent(bin_1)]
point_two = [1.132, lambda_mass_dist_50_80.GetBinContent(bin_2)]
bg_starting_params_50_80 = get_straight_line(point_one, point_two)
lambda_mass_fit_50_80 = rt.TF1("lambda_mass_fit_50_80", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)", 1.102, 1.127)
lambda_mass_fit_50_80.SetNpx(1000)
lambda_mass_fit_50_80.SetParameter(0, 1.74e03)
lambda_mass_fit_50_80.SetParameter(1, 1.116)
lambda_mass_fit_50_80.SetParameter(2, 1.30576e-03 )
lambda_mass_fit_50_80.SetParameter(3, 1.804166e-03)
lambda_mass_fit_50_80.SetParameter(4, bg_starting_params_50_80[0])
lambda_mass_fit_50_80.SetParameter(5, bg_starting_params_50_80[1])
lambda_mass_dist_fit_50_80 = lambda_mass_dist_50_80.Clone("lambda_mass_dist_fit_50_80")


lambda_mass_bg_fit_50_80 = rt.TF1("bg_fit_50_80", "pol1", 1.096, 1.136)
lambda_mass_bg_fit_50_80.SetNpx(1000)
lambda_mass_bg_fit_50_80.SetParameter(0, bg_starting_params_50_80[0])
lambda_mass_bg_fit_50_80.SetParameter(1, bg_starting_params_50_80[1])
lambda_mass_dist_fit_50_80.Fit(lambda_mass_fit_50_80, "RS")

voigt_fit_50_80 = rt.TF1("voigt_fit_50_80", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.102, 1.127)
voigt_fit_50_80.SetNpx(1000)
voigt_fit_50_80.SetParameter(0, lambda_mass_fit_50_80.GetParameter(0))
voigt_fit_50_80.SetParameter(1, lambda_mass_fit_50_80.GetParameter(1))
voigt_fit_50_80.SetParameter(2, lambda_mass_fit_50_80.GetParameter(2))
voigt_fit_50_80.SetParameter(3, lambda_mass_fit_50_80.GetParameter(3))


residual_50_80 = lambda_mass_dist_50_80.Clone("residual_50_80")
residual_50_80.GetXaxis().SetRangeUser(1.094, 1.142)
residual_50_80.Add(lambda_mass_bg_fit_50_80, -1)


RSB_region_50_80 = rt.TBox(RSB_MIN, 0, RSB_MAX, lambda_mass_dist_50_80.GetMaximum()*1.055)
RSB_region_50_80.SetLineColor(rt.kRed)
RSB_region_50_80.SetFillColor(rt.kRed)
RSB_region_50_80.SetFillStyle(3003)


RSB_min_line_50_80 = rt.TLine(RSB_MIN, 0, RSB_MIN, lambda_mass_dist_50_80.GetMaximum()*1.05)
RSB_min_line_50_80.SetLineColor(rt.kRed)
RSB_min_line_50_80.SetLineWidth(2)
RSB_min_line_50_80.SetLineStyle(2)


RSB_max_line_50_80 = rt.TLine(RSB_MAX, 0, RSB_MAX, lambda_mass_dist_50_80.GetMaximum()*1.05)
RSB_max_line_50_80.SetLineColor(rt.kRed)
RSB_max_line_50_80.SetLineWidth(2)
RSB_max_line_50_80.SetLineStyle(2)

left_signal_bin_50_80 = lambda_mass_dist_50_80.FindBin(SIG_MIN)
right_signal_bin_50_80 = lambda_mass_dist_50_80.FindBin(SIG_MAX)

full_region_bin_1 = residual_50_80.FindBin(FULL_REGION_MIN)
full_region_bin_2 = residual_50_80.FindBin(FULL_REGION_MAX)

lambda_bg_50_80 = 0
lambda_total_50_80 = 0
for bin_num in range(left_signal_bin_50_80, right_signal_bin_50_80 + 1):
    lambda_total_50_80 += lambda_mass_dist_50_80.GetBinContent(bin_num)
    lambda_bg_50_80 += lambda_mass_bg_fit_50_80.Eval(lambda_mass_dist_50_80.GetBinCenter(bin_num))


lambda_signal_50_80 = lambda_total_50_80 - lambda_bg_50_80
lambda_signal_total_ratio_50_80 = lambda_signal_50_80/lambda_total_50_80

signal_scale_50_80 = residual_50_80.Integral(full_region_bin_1, full_region_bin_2)/residual_50_80.Integral(left_signal_bin_50_80, right_signal_bin_50_80)

output_file.cd()
RSB_region_50_80.Write("RSB_region_50_80")
RSB_min_line_50_80.Write("RSB_min_line_50_80")
RSB_max_line_50_80.Write("RSB_max_line_50_80")

lambda_mass_dist_50_80.Write()
lambda_mass_dist_fit_50_80.Write()
lambda_mass_fit_50_80.Write()
lambda_mass_bg_fit_50_80.Write()
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
h_lambda_2d_mixcor_rsb_50_80 = make_mixed_corrections(h_lambda_50_80, h_lambda_mixed_50_80, RSB_MIN, RSB_MAX)

h_h_2d_mixcor_50_80 = make_mixed_corrections(h_h_50_80, h_h_mixed_50_80, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_50_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_50_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_h_2d_mixcor_50_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)

h_lambda_2d_mixcor_sig_50_80.SetName("h_lambda_2d_mixcor_sig_50_80")
h_lambda_2d_mixcor_rsb_50_80.SetName("h_lambda_2d_mixcor_rsb_50_80")
h_h_2d_mixcor_50_80.SetName("h_h_2d_mixcor_50_80")




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_50_80.Scale(1.0/num_trigs_50_80)
h_lambda_2d_mixcor_rsb_50_80.Scale(1.0/num_trigs_50_80)
h_h_2d_mixcor_50_80.Scale(1.0/num_trigs_50_80)

output_file.cd()
h_lambda_2d_nomixcor_50_80.Write()
h_lambda_mixed_2d_50_80.Write()
h_h_2d_nomixcor_50_80.Write()
h_h_mixed_2d_50_80.Write()
h_lambda_2d_mixcor_sig_50_80.Write()
h_lambda_2d_mixcor_rsb_50_80.Write()
h_h_2d_mixcor_50_80.Write()

# SIDEBAND SUBTRACTION SECTION
if DO_SIDEBAND_SUBTRACTION:

    # Normalize side band to 1
    h_lambda_2d_mixcor_rsb_50_80.Scale(1/h_lambda_2d_mixcor_rsb_50_80.Integral())


    # using RSB for sideband subtraction
    h_lambda_2d_subtracted_50_80 = h_lambda_2d_mixcor_sig_50_80.Clone("h_lambda_2d_subtracted_50_80")
    bg_integral_50_80 = (1 - lambda_signal_total_ratio_50_80)*h_lambda_2d_subtracted_50_80.Integral()
    h_lambda_2d_subtracted_50_80.Add(h_lambda_2d_mixcor_rsb_50_80, -bg_integral_50_80)

    # save the RSB deltaphi distribution
    h_lambda_dphi_rsb_50_80 = h_lambda_2d_mixcor_rsb_50_80.ProjectionY("h_lambda_dphi_rsb_50_80")
    output_file.cd()
    h_lambda_dphi_rsb_50_80.Write()

else:
    h_lambda_2d_subtracted_50_80 = h_lambda_2d_mixcor_sig_50_80.Clone("h_lambda_2d_subtracted_50_80")


# Correct for signal scale
h_lambda_2d_subtracted_50_80.Scale(signal_scale_50_80)
# Correct for PID loss
h_lambda_2d_subtracted_50_80.Scale(PID_CORRECTION)
# Correct for branching ratio
h_lambda_2d_subtracted_50_80.Scale(1/BRANCHING_RATIO)

# Correct for two-track efficiency
h_lambda_2d_subtracted_50_80 = apply_two_track_correction(h_lambda_2d_subtracted_50_80, TWOTRACK_TEMPLATE)



# INTEGRAL AND RATIO SECTION
h_lambda_dphi_subtracted_50_80 = h_lambda_2d_subtracted_50_80.ProjectionY("h_lambda_dphi_subtracted_50_80")

# this is completely insane but i need to scale the PID stuff exactly prior to fitting
if IS_WIDE_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_50_80.Scale(0.9215763047724077)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_50_80.Scale(0.9263212321482822)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_50_80.Scale(0.919115638819075)
elif IS_NARROW_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_50_80.Scale(1.1569937866811322)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_50_80.Scale(1.1573270204439357)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_50_80.Scale(1.1450656317216028)
elif IS_TOF_PID:
    if PT_MODE == 0:
        h_lambda_dphi_subtracted_50_80.Scale(0.9938083447210315)
    elif PT_MODE == 1:
        h_lambda_dphi_subtracted_50_80.Scale(1.7314056184347324)
    elif PT_MODE == 2:
        h_lambda_dphi_subtracted_50_80.Scale(0.831187045177973)

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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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

elif USE_FIT:

    # fitting to v2 component + four gaussians
    fit_string_50_80 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    fit_string_50_80 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    fit_function_50_80 = rt.TF1("fit_function_50_80", fit_string_50_80, -2, 6)
    fit_function_50_80.SetNpx(1000)
    if PT_MODE == 0:
        fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][0])
        fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][2])
        fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][0])
        fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][2])
        fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][0])
        fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][2])
        fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    fit_function_50_80.SetParLimits(3, 5e-4, 0.6)
    fit_function_50_80.SetParameter(3, 0.009)
    fit_function_50_80.FixParameter(4, 0)
    fit_function_50_80.SetParLimits(5, 0.1, 2)
    fit_function_50_80.SetParameter(5, 0.8)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    fit_function_50_80.SetParLimits(6, 0,  0.6)
    fit_function_50_80.SetParameter(6, 0)
    fit_function_50_80.FixParameter(7, 2*rt.TMath.Pi())
    fit_function_50_80.SetParLimits(8, 0.1, 2)
    fit_function_50_80.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    fit_function_50_80.SetParLimits(9, 5e-4, 0.6)
    fit_function_50_80.SetParameter(9, 0.003)
    fit_function_50_80.FixParameter(10, rt.TMath.Pi())
    fit_function_50_80.SetParLimits(11, 0.1, 2)
    fit_function_50_80.SetParameter(11, 1)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    fit_function_50_80.SetParLimits(12, 0, 0.6)
    fit_function_50_80.SetParameter(12, 0)
    fit_function_50_80.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    fit_function_50_80.SetParLimits(14, 0.1, 2)
    fit_function_50_80.SetParameter(14, 0.1)

    # setting parameters for constant background
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

    fit_function_50_80.SetParameter(12, ue_avg_50_80)


    h_lambda_dphi_subtracted_with_fit_50_80 = h_lambda_dphi_subtracted_50_80.Clone("h_lambda_dphi_subtracted_with_fit_50_80")
    fit_result_50_80 = h_lambda_dphi_subtracted_with_fit_50_80.Fit(fit_function_50_80, "R")

    ue_avg_fit_50_80 = rt.TF1("ue_avg_fit_50_80", "pol0", -2, 6)
    ue_avg_fit_50_80.SetParameter(0, ue_avg_50_80)
    ue_avg_fit_50_80.SetParError(0, ue_avg_error_50_80)

elif USE_V2:
    ue_avg_50_80 = (h_lambda_dphi_subtracted_50_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(16))/6

    v2_fit_50_80 = rt.TF1("v2_fit_50_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    ue_avg_50_80 = (h_lambda_dphi_subtracted_50_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_50_80.GetBinContent(16))/6

    v2_fit_50_80 = rt.TF1("v2_fit_50_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"

    von_fit_50_80 = rt.TF1("von_fit_50_80", von_fit_string, -2, 6)

    scale = 0.6
    if PT_MODE == 0:
        von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][0])
        von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][2])
        von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][3])
        # von_fit_50_80.SetParLimits(3, scale*0.015/3, scale*0.015*3)
        von_fit_50_80.SetParameter(3, scale*0.015)
        # von_fit_50_80.SetParLimits(4, (7.6/scale)/3, (7.6/scale)*3)
        von_fit_50_80.SetParameter(4, 7.6/scale)
        # von_fit_50_80.SetParLimits(5, scale*0.014/3, scale*0.014*3)
        von_fit_50_80.SetParameter(5, scale*0.014)
        # von_fit_50_80.SetParLimits(6, 1, 10)
        von_fit_50_80.SetParameter(6, 5)
    elif PT_MODE == 1:
        von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][0])
        von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][2])
        von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][3])
        # von_fit_50_80.SetParLimits(3, scale*0.02/6, scale*0.02*6)
        von_fit_50_80.SetParameter(3, scale*0.02)
        # von_fit_50_80.SetParLimits(4, (4.6/scale)/3, (4.6/scale)*3)
        von_fit_50_80.SetParameter(4, 4.6/scale)
        # von_fit_50_80.SetParLimits(5, scale*0.057/6, 6*scale*0.057)
        von_fit_50_80.SetParameter(5, scale*0.057)
        # von_fit_50_80.SetParLimits(6, 1, 10)
        von_fit_50_80.SetParameter(6, 3)
    elif PT_MODE == 2:
        von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][0])
        von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][2])
        von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][3])
        # von_fit_50_80.SetParLimits(3, scale*0.009/3, scale*0.009*3)
        von_fit_50_80.SetParameter(3, scale*0.009)
        # von_fit_50_80.SetParLimits(4, (8.88/scale)/6, (8.88/scale)*6)
        von_fit_50_80.SetParameter(4, 8.88/scale)
        # von_fit_50_80.SetParLimits(5, scale*0.0076/6, scale*0.0076*6)
        von_fit_50_80.SetParameter(5, scale*0.0076)
        # von_fit_50_80.SetParLimits(6, 1, 9)
        von_fit_50_80.SetParameter(6, 3)

    h_lambda_dphi_subtracted_with_fit_50_80 = h_lambda_dphi_subtracted_50_80.Clone("h_lambda_dphi_subtracted_with_fit_50_80")
    if PT_MODE == 1: # low pt has away-side width too large for von mises fit, only works if we symmetrize the peaks
        fit_result_50_80 = h_lambda_dphi_subtracted_with_fit_50_80.Fit(von_fit_50_80, "R", "", -rt.TMath.Pi()/2, rt.TMath.Pi() - EPSILON)
    else:
        fit_result_50_80 = h_lambda_dphi_subtracted_with_fit_50_80.Fit(von_fit_50_80, "R")


    if PT_MODE == 0:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][0])
        v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][1])
        v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][2])
        v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_50_80_trigger_4_8_assoc_25_4"][3])


else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_lambda_dphi_subtracted_50_80_zeroed = h_lambda_dphi_subtracted_50_80.Clone("h_lambda_dphi_subtracted_50_80_zeroed")
    h_lambda_dphi_subtracted_50_80_zeroed.Add(ue_avg_fit_50_80, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_lambda_total_integral_error_50_80 = fit_function_50_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_50_80 = fit_function_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_50_80 = ue_avg_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_50_80 = ue_avg_error_50_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_50_80_halved = h_lambda_ue_integral_error_50_80/2
    h_lambda_near_integral_error_50_80 = math.sqrt((fit_function_50_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_near_integral_50_80 = fit_function_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_50_80 = math.sqrt((fit_function_50_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_away_integral_50_80 = fit_function_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_lambda_total_integral_error_50_80 = von_fit_50_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_50_80 = von_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_50_80 = v2_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_50_80 = v2_fit_50_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_50_80_halved = h_lambda_ue_integral_error_50_80/2
    h_lambda_near_integral_error_50_80 = math.sqrt((von_fit_50_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_near_integral_50_80 = von_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_50_80 = math.sqrt((von_fit_50_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_away_integral_50_80 = von_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_dphi_subtracted_50_80_zeroed = h_lambda_dphi_subtracted_50_80.Clone("h_lambda_dphi_subtracted_50_80_zeroed")
    h_lambda_dphi_subtracted_50_80_zeroed.Add(v2_fit_50_80, -1)

elif USE_V2:

    DPHI_BINS = h_lambda_dphi_subtracted_50_80.GetNbinsX()
    h_lambda_total_integral_50_80 = 0
    h_lambda_near_integral_50_80 = 0
    h_lambda_away_integral_50_80 = 0
    h_lambda_ue_integral_50_80 = v2_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

    h_lambda_total_integral_error_50_80 = 0
    h_lambda_near_integral_error_50_80 = 0
    h_lambda_away_integral_error_50_80 = 0
    h_lambda_ue_integral_error_50_80 = v2_fit_50_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_50_80_halved = h_lambda_ue_integral_error_50_80/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_50_80 += h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num)
        h_lambda_total_integral_error_50_80 += (h_lambda_dphi_subtracted_50_80.GetBinError(bin_num))**2
        part = h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num) - v2_fit_50_80.Eval(h_lambda_dphi_subtracted_50_80.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_50_80 += part
            h_lambda_near_integral_error_50_80 += (h_lambda_dphi_subtracted_50_80.GetBinError(bin_num))**2
        else:
            h_lambda_away_integral_50_80 += part
            h_lambda_away_integral_error_50_80 += (h_lambda_dphi_subtracted_50_80.GetBinError(bin_num))**2
        
    h_lambda_total_integral_error_50_80 = math.sqrt(h_lambda_total_integral_error_50_80)
    h_lambda_near_integral_error_50_80 = math.sqrt(h_lambda_near_integral_error_50_80)
    h_lambda_near_integral_error_50_80 = math.sqrt(h_lambda_near_integral_error_50_80**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_away_integral_error_50_80 = math.sqrt(h_lambda_away_integral_error_50_80)
    h_lambda_away_integral_error_50_80 = math.sqrt(h_lambda_away_integral_error_50_80**2 + h_lambda_ue_integral_error_50_80_halved**2)

else: 

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
    h_lambda_ue_integral_error_50_80 = ue_avg_error_50_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_50_80_halved = h_lambda_ue_integral_error_50_80/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_50_80 += h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num)
        h_lambda_total_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2
        part = h_lambda_dphi_subtracted_50_80.GetBinContent(bin_num) - ue_avg_50_80
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_50_80 += part
            h_lambda_near_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2
        else:
            h_lambda_away_integral_50_80 += part
            h_lambda_away_integral_error_50_80 += h_lambda_dphi_subtracted_50_80.GetBinError(bin_num)**2

    h_lambda_total_integral_error_50_80 = math.sqrt(h_lambda_total_integral_error_50_80)
    h_lambda_near_integral_error_50_80 = math.sqrt(h_lambda_near_integral_error_50_80)
    h_lambda_near_integral_error_50_80 = math.sqrt(h_lambda_near_integral_error_50_80**2 + h_lambda_ue_integral_error_50_80_halved**2)
    h_lambda_away_integral_error_50_80 = math.sqrt(h_lambda_away_integral_error_50_80)
    h_lambda_away_integral_error_50_80 = math.sqrt(h_lambda_away_integral_error_50_80**2 + h_lambda_ue_integral_error_50_80_halved**2)


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

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
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
elif USE_FIT:

    # fitting to four gaussians + a constant background
    hh_fit_function_50_80 = rt.TF1("hh_fit_function_50_80", "gaus(0) + gaus(3) + gaus(6) + gaus(9) + pol0(12)", -2, 6)
    hh_fit_function_50_80.SetNpx(1000)

    # fitting to v2 component + four gaussians
    hh_fit_string_50_80 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_fit_string_50_80 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    hh_fit_function_50_80 = rt.TF1("hh_fit_function_50_80", hh_fit_string_50_80, -2, 6)
    hh_fit_function_50_80.SetNpx(1000)
    if PT_MODE == 0:
        hh_fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][0])
        hh_fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][2])
        hh_fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][0])
        hh_fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][2])
        hh_fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_fit_function_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][0])
        hh_fit_function_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][2])
        hh_fit_function_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    hh_fit_function_50_80.SetParLimits(3, 0.02, 1)
    hh_fit_function_50_80.SetParameter(3, 0.16)
    hh_fit_function_50_80.FixParameter(4, 0)
    hh_fit_function_50_80.SetParLimits(5, 0.1, 1)
    hh_fit_function_50_80.SetParameter(5, 0.3)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    hh_fit_function_50_80.SetParLimits(6, 0, 1)
    hh_fit_function_50_80.SetParameter(6, 0)
    hh_fit_function_50_80.FixParameter(7, 2*rt.TMath.Pi())
    hh_fit_function_50_80.SetParLimits(8, 0.1, 2)
    hh_fit_function_50_80.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    hh_fit_function_50_80.SetParLimits(9, 0.02, 1)
    hh_fit_function_50_80.SetParameter(9, 0.3)
    hh_fit_function_50_80.FixParameter(10, rt.TMath.Pi())
    hh_fit_function_50_80.SetParLimits(11, 0.4, 2)
    hh_fit_function_50_80.SetParameter(11, 0.7)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    hh_fit_function_50_80.SetParLimits(12, 0, 1)
    hh_fit_function_50_80.SetParameter(12, 0)
    hh_fit_function_50_80.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    hh_fit_function_50_80.SetParLimits(14, 0.1, 3)
    hh_fit_function_50_80.SetParameter(14, 0.1)


    # setting parameters for constant background
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

    hh_fit_function_50_80.SetParameter(12, hh_ue_avg_50_80)


    h_h_dphi_with_fit_50_80 = h_h_dphi_50_80.Clone("h_h_dphi_with_fit_50_80")
    h_h_dphi_with_fit_50_80.Fit(hh_fit_function_50_80, "R")

    hh_ue_avg_fit_50_80 = rt.TF1("hh_ue_avg_fit_50_80", "pol0", -2, 6)
    hh_ue_avg_fit_50_80.SetParameter(0, hh_ue_avg_50_80)
    hh_ue_avg_fit_50_80.SetParError(0, hh_ue_avg_error_50_80)

elif USE_V2:
    hh_ue_avg_50_80 = (h_h_dphi_50_80.GetBinContent(1) 
                   + h_h_dphi_50_80.GetBinContent(2)
                   + h_h_dphi_50_80.GetBinContent(7)
                   + h_h_dphi_50_80.GetBinContent(8)
                   + h_h_dphi_50_80.GetBinContent(9)
                   + h_h_dphi_50_80.GetBinContent(16))/6

    hh_v2_fit_50_80 = rt.TF1("hh_v2_fit_50_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    hh_ue_avg_50_80 = (h_h_dphi_50_80.GetBinContent(1) 
                   + h_h_dphi_50_80.GetBinContent(2)
                   + h_h_dphi_50_80.GetBinContent(7)
                   + h_h_dphi_50_80.GetBinContent(8)
                   + h_h_dphi_50_80.GetBinContent(9)
                   + h_h_dphi_50_80.GetBinContent(16))/6

    hh_v2_fit_50_80 = rt.TF1("hh_v2_fit_50_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    hh_von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    hh_von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
    hh_von_fit_50_80 = rt.TF1("hh_von_fit_50_80", hh_von_fit_string, -2, 6)

    if PT_MODE == 0:
        hh_von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][0])
        hh_von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][2])
        hh_von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][0])
        hh_von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][2])
        hh_von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_von_fit_50_80.FixParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][0])
        hh_von_fit_50_80.FixParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][2])
        hh_von_fit_50_80.FixParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][3])

    hh_von_fit_50_80.SetParLimits(3, 0, 1)
    hh_von_fit_50_80.SetParameter(3, 0.02)
    hh_von_fit_50_80.SetParLimits(4, 0, 100)
    hh_von_fit_50_80.SetParameter(4, 1)
    hh_von_fit_50_80.SetParLimits(5, 0, 1)
    hh_von_fit_50_80.SetParameter(5, 0.01)
    hh_von_fit_50_80.SetParLimits(6, 0, 100)
    hh_von_fit_50_80.SetParameter(6, 1)

    h_h_dphi_with_fit_50_80 = h_h_dphi_50_80.Clone("h_h_dphi_with_fit_50_80")
    hh_fit_result_50_80 = h_h_dphi_with_fit_50_80.Fit(hh_von_fit_50_80, "R")

    if PT_MODE == 0:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_50_80.SetParameter(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_50_80.SetParError(0, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_50_80.SetParameter(1, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_50_80.SetParameter(2, v2_fitpars_dict["h_h_cent_50_80_trigger_4_8_assoc_25_4"][3])

else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_h_dphi_50_80_zeroed = h_h_dphi_50_80.Clone("h_h_dphi_50_80_zeroed")
    h_h_dphi_50_80_zeroed.Add(hh_ue_avg_fit_50_80, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_h_total_integral_error_50_80 = hh_fit_function_50_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_50_80 = hh_fit_function_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_50_80 = hh_ue_avg_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_50_80 = hh_ue_avg_error_50_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_50_80_halved = h_h_ue_integral_error_50_80/2
    h_h_near_integral_error_50_80 = math.sqrt((hh_fit_function_50_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_near_integral_50_80 = hh_fit_function_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_50_80 = math.sqrt((hh_fit_function_50_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_away_integral_50_80 = hh_fit_function_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_h_total_integral_error_50_80 = hh_von_fit_50_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_50_80 = hh_von_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_50_80 = hh_v2_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_50_80 = hh_v2_fit_50_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_50_80_halved = h_h_ue_integral_error_50_80/2
    h_h_near_integral_error_50_80 = math.sqrt((hh_von_fit_50_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_near_integral_50_80 = hh_von_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_50_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_50_80 = math.sqrt((hh_von_fit_50_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_away_integral_50_80 = hh_von_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_50_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_dphi_50_80_zeroed = h_h_dphi_50_80.Clone("h_h_dphi_50_80_zeroed")
    h_h_dphi_50_80_zeroed.Add(hh_v2_fit_50_80, -1)

elif USE_V2:

    DPHI_BINS = h_h_dphi_50_80.GetNbinsX()
    h_h_total_integral_50_80 = 0
    h_h_near_integral_50_80 = 0
    h_h_away_integral_50_80 = 0
    h_h_ue_integral_50_80 = hh_v2_fit_50_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_error_50_80 = 0
    h_h_near_integral_error_50_80 = 0
    h_h_away_integral_error_50_80 = 0
    h_h_ue_integral_error_50_80 = hh_v2_fit_50_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_50_80_halved = h_h_ue_integral_error_50_80/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_50_80 += h_h_dphi_50_80.GetBinContent(bin_num)
        h_h_total_integral_error_50_80 += (h_h_dphi_50_80.GetBinError(bin_num))**2
        part = h_h_dphi_50_80.GetBinContent(bin_num) - hh_v2_fit_50_80.Eval(h_h_dphi_50_80.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_50_80 += part
            h_h_near_integral_error_50_80 += (h_h_dphi_50_80.GetBinError(bin_num))**2
        else:
            h_h_away_integral_50_80 += part
            h_h_away_integral_error_50_80 += (h_h_dphi_50_80.GetBinError(bin_num))**2
    h_h_total_integral_error_50_80 = math.sqrt(h_h_total_integral_error_50_80)
    h_h_near_integral_error_50_80 = math.sqrt(h_h_near_integral_error_50_80)
    h_h_near_integral_error_50_80 = math.sqrt(h_h_near_integral_error_50_80**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_away_integral_error_50_80 = math.sqrt(h_h_away_integral_error_50_80)
    h_h_away_integral_error_50_80 = math.sqrt(h_h_away_integral_error_50_80**2 + h_h_ue_integral_error_50_80_halved**2)

else: 

    h_h_dphi_50_80_zeroed = h_h_dphi_50_80.Clone("h_h_dphi_50_80_zeroed")
    h_h_dphi_50_80_zeroed.Add(hh_ue_line_50_80, -1)

    DPHI_BINS = h_h_dphi_50_80.GetNbinsX()
    h_h_total_integral_50_80 = 0
    h_h_near_integral_50_80 = 0
    h_h_away_integral_50_80 = 0
    h_h_ue_integral_50_80 = hh_ue_avg_50_80*DPHI_BINS
    h_h_total_integral_error_50_80 = 0
    h_h_near_integral_error_50_80 = 0
    h_h_away_integral_error_50_80 = 0
    h_h_ue_integral_error_50_80 = hh_ue_avg_error_50_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_50_80_halved = h_h_ue_integral_error_50_80/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_50_80 += h_h_dphi_50_80.GetBinContent(bin_num)
        h_h_total_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
        part = h_h_dphi_50_80.GetBinContent(bin_num) - hh_ue_avg_50_80
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_50_80 += part
            h_h_near_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
        else:
            h_h_away_integral_50_80 += part
            h_h_away_integral_error_50_80 += h_h_dphi_50_80.GetBinError(bin_num)**2
    h_h_total_integral_error_50_80 = math.sqrt(h_h_total_integral_error_50_80)
    h_h_near_integral_error_50_80 = math.sqrt(h_h_near_integral_error_50_80)
    h_h_near_integral_error_50_80 = math.sqrt(h_h_near_integral_error_50_80**2 + h_h_ue_integral_error_50_80_halved**2)
    h_h_away_integral_error_50_80 = math.sqrt(h_h_away_integral_error_50_80)
    h_h_away_integral_error_50_80 = math.sqrt(h_h_away_integral_error_50_80**2 + h_h_ue_integral_error_50_80_halved**2)


output_file.cd()
h_lambda_2d_subtracted_50_80.Write()

h_lambda_dphi_subtracted_50_80.Write()
h_h_dphi_50_80.Write()



if USE_FIT:
    fit_function_50_80.Write()
    ue_avg_fit_50_80.Write()
    hh_fit_function_50_80.Write()
    hh_ue_avg_fit_50_80.Write()
    h_lambda_dphi_subtracted_with_fit_50_80.Write()
    h_h_dphi_with_fit_50_80.Write()
    h_lambda_dphi_subtracted_50_80_zeroed.Write()
    h_h_dphi_50_80_zeroed.Write()
elif USE_V2:
    v2_fit_50_80.Write()
    hh_v2_fit_50_80.Write()
elif USE_VON:
    von_fit_50_80.Write()
    hh_von_fit_50_80.Write()
    v2_fit_50_80.Write()
    hh_v2_fit_50_80.Write()
    h_lambda_dphi_subtracted_with_fit_50_80.Write()
    h_h_dphi_with_fit_50_80.Write()
    h_lambda_dphi_subtracted_50_80_zeroed.Write()
    h_h_dphi_50_80_zeroed.Write()
else:
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

    h_lambda_dphi_subtracted_50_80_zeroed.Write()
    h_h_dphi_50_80_zeroed.Write()


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
############################################# MB CENTRALITY ################################################
############################################################################################################
############################################################################################################

if IS_NORMAL_PID:
    input_file_0_80 = rt.TFile("../online/v0/output/v0_normal.root")
elif IS_WIDE_PID:
    input_file_0_80 = rt.TFile("../online/v0/output/v0_wide.root")
elif IS_NARROW_PID:
    input_file_0_80 = rt.TFile("../online/v0/output/v0_narrow.root")
elif IS_TOF_PID:
    input_file_0_80 = rt.TFile("../online/v0/output/v0_tof.root")
else:
    print("ERROR: No PID cut specified")
    exit(1)

input_list_0_80 = input_file_0_80.Get("h-lambda")
input_file_0_80.Close()

trig_dist_0_80 = input_list_0_80.FindObject("fTriggerDistEff_highestPt") if DO_HIGHEST_PT else input_list_0_80.FindObject("fTriggerDistEff")
lambda_dist_0_80 = input_list_0_80.FindObject("fTriggeredLambdaDist")


h_h_0_80 = input_list_0_80.FindObject("fDphiHHEff_highestPt") if DO_HIGHEST_PT else input_list_0_80.FindObject("fDphiHHEff")
h_h_mixed_0_80 = input_list_0_80.FindObject("fDphiHHMixed_highestPt") if DO_HIGHEST_PT else input_list_0_80.FindObject("fDphiHHMixed")

h_lambda_0_80 = input_list_0_80.FindObject("fDphiHLambdaEff_highestPt") if DO_HIGHEST_PT else input_list_0_80.FindObject("fDphiHLambdaEff")
h_lambda_mixed_0_80 = input_list_0_80.FindObject("fDphiHLambdaMixed_highestPt") if DO_HIGHEST_PT else input_list_0_80.FindObject("fDphiHLambdaMixed")


# Setting the trigger Pt (this is never changed again)
trig_dist_0_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_0_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_h_mixed_0_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_0_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)
h_lambda_mixed_0_80.GetAxis(0).SetRangeUser(TRIG_PT_LOW, TRIG_PT_HIGH)


# Setting the associated Pt (this is never changed again)
lambda_dist_0_80.GetAxis(0).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_0_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_h_mixed_0_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_0_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)
h_lambda_mixed_0_80.GetAxis(1).SetRangeUser(ASSOC_PT_LOW, ASSOC_PT_HIGH)


### TRIGGER SECTION ### 

trig_dist_0_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)
lambda_dist_0_80.GetAxis(2).SetRangeUser(ETA_MIN, ETA_MAX)


# Getting the single-particle trigger distributions
trig_pt_dist_0_80 = trig_dist_0_80.Projection(0).Clone("trig_pt_dist_0_80")
trig_phi_dist_0_80 = trig_dist_0_80.Projection(1).Clone("trig_phi_dist_0_80")
trig_eta_dist_0_80 = trig_dist_0_80.Projection(2).Clone("trig_eta_dist_0_80")
trig_2d_dist_0_80 = trig_dist_0_80.Projection(0, 3).Clone("trig_2d_dist_0_80")

# Get total number of triggers (used for per-trigger normalization)
num_trigs_0_80 = trig_2d_dist_0_80.Integral()

output_file.cd()
trig_pt_dist_0_80.Write()
trig_phi_dist_0_80.Write()
trig_eta_dist_0_80.Write()
trig_2d_dist_0_80.Write()

### SIGNAL ANALYSIS SECTION ###

lambda_mass_dist_0_80 = lambda_dist_0_80.Projection(3).Clone("lambda_mass_dist_0_80")

bin_1 = lambda_mass_dist_0_80.FindBin(1.102)
bin_2 = lambda_mass_dist_0_80.FindBin(1.132)
point_one = [1.102, lambda_mass_dist_0_80.GetBinContent(bin_1)]
point_two = [1.132, lambda_mass_dist_0_80.GetBinContent(bin_2)]
bg_starting_params_0_80 = get_straight_line(point_one, point_two)
lambda_mass_fit_0_80 = rt.TF1("lambda_mass_fit_0_80", "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)", 1.102, 1.127)
lambda_mass_fit_0_80.SetNpx(1000)
lambda_mass_fit_0_80.SetParameter(0, 1.74e03)
lambda_mass_fit_0_80.SetParameter(1, 1.116)
lambda_mass_fit_0_80.SetParameter(2, 1.30576e-03 )
lambda_mass_fit_0_80.SetParameter(3, 1.804166e-03)
lambda_mass_fit_0_80.SetParameter(4, bg_starting_params_0_80[0])
lambda_mass_fit_0_80.SetParameter(5, bg_starting_params_0_80[1])
lambda_mass_dist_fit_0_80 = lambda_mass_dist_0_80.Clone("lambda_mass_dist_fit_0_80")


lambda_mass_bg_fit_0_80 = rt.TF1("bg_fit_0_80", "pol1", 1.096, 1.136)
lambda_mass_bg_fit_0_80.SetNpx(1000)
lambda_mass_bg_fit_0_80.SetParameter(0, bg_starting_params_0_80[0])
lambda_mass_bg_fit_0_80.SetParameter(1, bg_starting_params_0_80[1])
lambda_mass_dist_fit_0_80.Fit(lambda_mass_fit_0_80, "RS")

voigt_fit_0_80 = rt.TF1("voigt_fit_0_80", "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.102, 1.127)
voigt_fit_0_80.SetNpx(1000)
voigt_fit_0_80.SetParameter(0, lambda_mass_fit_0_80.GetParameter(0))
voigt_fit_0_80.SetParameter(1, lambda_mass_fit_0_80.GetParameter(1))
voigt_fit_0_80.SetParameter(2, lambda_mass_fit_0_80.GetParameter(2))
voigt_fit_0_80.SetParameter(3, lambda_mass_fit_0_80.GetParameter(3))


residual_0_80 = lambda_mass_dist_0_80.Clone("residual_0_80")
residual_0_80.GetXaxis().SetRangeUser(1.094, 1.142)
residual_0_80.Add(lambda_mass_bg_fit_0_80, -1)


RSB_region_0_80 = rt.TBox(RSB_MIN, 0, RSB_MAX, lambda_mass_dist_0_80.GetMaximum()*1.055)
RSB_region_0_80.SetLineColor(rt.kRed)
RSB_region_0_80.SetFillColor(rt.kRed)
RSB_region_0_80.SetFillStyle(3003)


RSB_min_line_0_80 = rt.TLine(RSB_MIN, 0, RSB_MIN, lambda_mass_dist_0_80.GetMaximum()*1.05)
RSB_min_line_0_80.SetLineColor(rt.kRed)
RSB_min_line_0_80.SetLineWidth(2)
RSB_min_line_0_80.SetLineStyle(2)


RSB_max_line_0_80 = rt.TLine(RSB_MAX, 0, RSB_MAX, lambda_mass_dist_0_80.GetMaximum()*1.05)
RSB_max_line_0_80.SetLineColor(rt.kRed)
RSB_max_line_0_80.SetLineWidth(2)
RSB_max_line_0_80.SetLineStyle(2)

left_signal_bin_0_80 = lambda_mass_dist_0_80.FindBin(SIG_MIN)
right_signal_bin_0_80 = lambda_mass_dist_0_80.FindBin(SIG_MAX)

full_region_bin_1 = residual_0_80.FindBin(FULL_REGION_MIN)
full_region_bin_2 = residual_0_80.FindBin(FULL_REGION_MAX)

lambda_bg_0_80 = 0
lambda_total_0_80 = 0
for bin_num in range(left_signal_bin_0_80, right_signal_bin_0_80 + 1):
    lambda_total_0_80 += lambda_mass_dist_0_80.GetBinContent(bin_num)
    lambda_bg_0_80 += lambda_mass_bg_fit_0_80.Eval(lambda_mass_dist_0_80.GetBinCenter(bin_num))


lambda_signal_0_80 = lambda_total_0_80 - lambda_bg_0_80
lambda_signal_total_ratio_0_80 = lambda_signal_0_80/lambda_total_0_80

signal_scale_0_80 = residual_0_80.Integral(full_region_bin_1, full_region_bin_2)/residual_0_80.Integral(left_signal_bin_0_80, right_signal_bin_0_80)

output_file.cd()
RSB_region_0_80.Write("RSB_region_0_80")
RSB_min_line_0_80.Write("RSB_min_line_0_80")
RSB_max_line_0_80.Write("RSB_max_line_0_80")

lambda_mass_dist_0_80.Write()
lambda_mass_dist_fit_0_80.Write()
lambda_mass_fit_0_80.Write()
lambda_mass_bg_fit_0_80.Write()
residual_0_80.Write()
RSB_region_0_80.Write()
RSB_min_line_0_80.Write()
RSB_max_line_0_80.Write()

### MIXED EVENT CORRECTION SECTION ###

axes = arr.array('i', [2, 3, 4, 5])
h_lambda_0_80 = h_lambda_0_80.Projection(4, axes)
h_lambda_mixed_0_80 = h_lambda_mixed_0_80.Projection(4, axes)

h_h_0_80 = h_h_0_80.Projection(2, 3, 4)
h_h_mixed_0_80 = h_h_mixed_0_80.Projection(2, 3, 4)


# Setting up 2-d correlation plots before the mixed event correction
h_lambda_2d_nomixcor_0_80 = h_lambda_0_80.Projection(0, 1).Clone("h_lambda_2d_nomixcor_0_80")
h_lambda_mixed_2d_0_80 = h_lambda_mixed_0_80.Projection(0, 1).Clone("h_lambda_mixed_2d_0_80")

h_h_2d_nomixcor_0_80 = h_h_0_80.Project3D("xye").Clone("h_h_2d_nomixcor_0_80")
h_h_mixed_2d_0_80 = h_h_mixed_0_80.Project3D("xye").Clone("h_h_mixed_2d_0_80")

h_lambda_2d_mixcor_sig_0_80 = make_mixed_corrections(h_lambda_0_80, h_lambda_mixed_0_80, SIG_MIN, SIG_MAX)
h_lambda_2d_mixcor_rsb_0_80 = make_mixed_corrections(h_lambda_0_80, h_lambda_mixed_0_80, RSB_MIN, RSB_MAX)

h_h_2d_mixcor_0_80 = make_mixed_corrections(h_h_0_80, h_h_mixed_0_80, SIG_MIN, SIG_MAX, is_hh=True)


h_lambda_2d_mixcor_sig_0_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_lambda_2d_mixcor_rsb_0_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)
h_h_2d_mixcor_0_80.GetXaxis().SetRangeUser(DELTA_ETA_MIN, DELTA_ETA_MAX)

h_lambda_2d_mixcor_sig_0_80.SetName("h_lambda_2d_mixcor_sig_0_80")
h_lambda_2d_mixcor_rsb_0_80.SetName("h_lambda_2d_mixcor_rsb_0_80")
h_h_2d_mixcor_0_80.SetName("h_h_2d_mixcor_0_80")




# per-trigger normalization done here
h_lambda_2d_mixcor_sig_0_80.Scale(1.0/num_trigs_0_80)
h_lambda_2d_mixcor_rsb_0_80.Scale(1.0/num_trigs_0_80)
h_h_2d_mixcor_0_80.Scale(1.0/num_trigs_0_80)

output_file.cd()
h_lambda_2d_nomixcor_0_80.Write()
h_lambda_mixed_2d_0_80.Write()
h_h_2d_nomixcor_0_80.Write()
h_h_mixed_2d_0_80.Write()
h_lambda_2d_mixcor_sig_0_80.Write()
h_lambda_2d_mixcor_rsb_0_80.Write()
h_h_2d_mixcor_0_80.Write()

# SIDEBAND SUBTRACTION SECTION
if DO_SIDEBAND_SUBTRACTION:

    # Normalize side band to 1
    h_lambda_2d_mixcor_rsb_0_80.Scale(1/h_lambda_2d_mixcor_rsb_0_80.Integral())


    # using RSB for sideband subtraction
    h_lambda_2d_subtracted_0_80 = h_lambda_2d_mixcor_sig_0_80.Clone("h_lambda_2d_subtracted_0_80")
    bg_integral_0_80 = (1 - lambda_signal_total_ratio_0_80)*h_lambda_2d_subtracted_0_80.Integral()
    h_lambda_2d_subtracted_0_80.Add(h_lambda_2d_mixcor_rsb_0_80, -bg_integral_0_80)

    # save the RSB deltaphi distribution
    h_lambda_dphi_rsb_0_80 = h_lambda_2d_mixcor_rsb_0_80.ProjectionY("h_lambda_dphi_rsb_0_80")
    output_file.cd()
    h_lambda_dphi_rsb_0_80.Write()

else:
    h_lambda_2d_subtracted_0_80 = h_lambda_2d_mixcor_sig_0_80.Clone("h_lambda_2d_subtracted_0_80")


# Correct for signal scale
h_lambda_2d_subtracted_0_80.Scale(signal_scale_0_80)
# Correct for PID loss
h_lambda_2d_subtracted_0_80.Scale(PID_CORRECTION)
# Correct for branching ratio
h_lambda_2d_subtracted_0_80.Scale(1/BRANCHING_RATIO)

# Correct for two-track efficiency
h_lambda_2d_subtracted_0_80 = apply_two_track_correction(h_lambda_2d_subtracted_0_80, TWOTRACK_TEMPLATE)



# INTEGRAL AND RATIO SECTION
h_lambda_dphi_subtracted_0_80 = h_lambda_2d_subtracted_0_80.ProjectionY("h_lambda_dphi_subtracted_0_80")

if USE_AVG_4:

    ue_line_0_80 = rt.TF1("ue_line_0_80", "pol0", -2, 6)
    ue_upper_line_0_80 = rt.TF1("ue_upper_line_0_80", "pol0", -2, 6)
    ue_lower_line_0_80 = rt.TF1("ue_lower_line_0_80", "pol0", -2, 6)
    zero_line_0_80 = rt.TF1("zero_line_0_80", "pol0", -2, 6)
    zero_upper_line_0_80 = rt.TF1("zero_upper_line_0_80", "pol0", -2, 6)
    zero_lower_line_0_80 = rt.TF1("zero_lower_line_0_80", "pol0", -2, 6)

    ue_avg_0_80 = (h_lambda_dphi_subtracted_0_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(16))/4
    ue_avg_error_0_80 = (1/4)*(math.sqrt(h_lambda_dphi_subtracted_0_80.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_0_80.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(16)**2))


    ue_line_0_80.SetParameter(0, ue_avg_0_80)
    ue_upper_line_0_80.SetParameter(0, ue_avg_0_80 + ue_avg_error_0_80)
    ue_lower_line_0_80.SetParameter(0, ue_avg_0_80 - ue_avg_error_0_80)

    zero_line_0_80.SetParameter(0, 0)
    zero_upper_line_0_80.SetParameter(0, ue_avg_error_0_80)
    zero_lower_line_0_80.SetParameter(0, -ue_avg_error_0_80)

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
    ue_line_0_80 = rt.TF1("ue_line_0_80", "pol0", -2, 6)
    ue_upper_line_0_80 = rt.TF1("ue_upper_line_0_80", "pol0", -2, 6)
    ue_lower_line_0_80 = rt.TF1("ue_lower_line_0_80", "pol0", -2, 6)
    zero_line_0_80 = rt.TF1("zero_line_0_80", "pol0", -2, 6)
    zero_upper_line_0_80 = rt.TF1("zero_upper_line_0_80", "pol0", -2, 6)
    zero_lower_line_0_80 = rt.TF1("zero_lower_line_0_80", "pol0", -2, 6)
    ue_avg_0_80 = (h_lambda_dphi_subtracted_0_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(16))/6

    ue_avg_error_0_80 = (1/6)*(math.sqrt(h_lambda_dphi_subtracted_0_80.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_0_80.GetBinError(2)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(7)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(16)**2))


    ue_line_0_80.SetParameter(0, ue_avg_0_80)
    ue_upper_line_0_80.SetParameter(0, ue_avg_0_80 + ue_avg_error_0_80)
    ue_lower_line_0_80.SetParameter(0, ue_avg_0_80 - ue_avg_error_0_80)

    zero_line_0_80.SetParameter(0, 0)
    zero_upper_line_0_80.SetParameter(0, ue_avg_error_0_80)
    zero_lower_line_0_80.SetParameter(0, -ue_avg_error_0_80)

elif USE_ZYAM:
    ue_line_0_80 = rt.TF1("ue_line_0_80", "pol0", -2, 6)
    ue_upper_line_0_80 = rt.TF1("ue_upper_line_0_80", "pol0", -2, 6)
    ue_lower_line_0_80 = rt.TF1("ue_lower_line_0_80", "pol0", -2, 6)
    zero_line_0_80 = rt.TF1("zero_line_0_80", "pol0", -2, 6)
    zero_upper_line_0_80 = rt.TF1("zero_upper_line_0_80", "pol0", -2, 6)
    zero_lower_line_0_80 = rt.TF1("zero_lower_line_0_80", "pol0", -2, 6)
    min_bin = h_lambda_dphi_subtracted_0_80.GetMinimumBin()
    ue_avg_0_80 = h_lambda_dphi_subtracted_0_80.GetBinContent(min_bin)
    ue_avg_error_0_80 = h_lambda_dphi_subtracted_0_80.GetBinError(min_bin)


    ue_line_0_80.SetParameter(0, ue_avg_0_80)
    ue_upper_line_0_80.SetParameter(0, ue_avg_0_80 + ue_avg_error_0_80)
    ue_lower_line_0_80.SetParameter(0, ue_avg_0_80 - ue_avg_error_0_80)

    zero_line_0_80.SetParameter(0, 0)
    zero_upper_line_0_80.SetParameter(0, ue_avg_error_0_80)
    zero_lower_line_0_80.SetParameter(0, -ue_avg_error_0_80)

elif USE_FIT:

    # fitting to v2 component + four gaussians
    fit_string_0_80 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    fit_string_0_80 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    fit_function_0_80 = rt.TF1("fit_function_0_80", fit_string_0_80, -2, 6)
    fit_function_0_80.SetNpx(1000)
    if PT_MODE == 0:
        fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][0])
        fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][2])
        fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][0])
        fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][2])
        fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][0])
        fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][2])
        fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    fit_function_0_80.SetParLimits(3, 5e-4, 0.6)
    fit_function_0_80.SetParameter(3, 0.009)
    fit_function_0_80.FixParameter(4, 0)
    fit_function_0_80.SetParLimits(5, 0.1, 2)
    fit_function_0_80.SetParameter(5, 0.8)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    fit_function_0_80.SetParLimits(6, 0,  0.6)
    fit_function_0_80.SetParameter(6, 0)
    fit_function_0_80.FixParameter(7, 2*rt.TMath.Pi())
    fit_function_0_80.SetParLimits(8, 0.1, 2)
    fit_function_0_80.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    fit_function_0_80.SetParLimits(9, 5e-4, 0.6)
    fit_function_0_80.SetParameter(9, 0.003)
    fit_function_0_80.FixParameter(10, rt.TMath.Pi())
    fit_function_0_80.SetParLimits(11, 0.1, 2)
    fit_function_0_80.SetParameter(11, 1)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    fit_function_0_80.SetParLimits(12, 0, 0.6)
    fit_function_0_80.SetParameter(12, 0)
    fit_function_0_80.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    fit_function_0_80.SetParLimits(14, 0.1, 2)
    fit_function_0_80.SetParameter(14, 0.1)

    # setting parameters for constant background
    ue_avg_0_80 = (h_lambda_dphi_subtracted_0_80.GetBinContent(1) 
                + h_lambda_dphi_subtracted_0_80.GetBinContent(2)
                + h_lambda_dphi_subtracted_0_80.GetBinContent(7)
                + h_lambda_dphi_subtracted_0_80.GetBinContent(8)
                + h_lambda_dphi_subtracted_0_80.GetBinContent(9)
                + h_lambda_dphi_subtracted_0_80.GetBinContent(16))/6

    ue_avg_error_0_80 = (1/6)*(math.sqrt(h_lambda_dphi_subtracted_0_80.GetBinError(1)**2 
                   + h_lambda_dphi_subtracted_0_80.GetBinError(2)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(7)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(8)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(9)**2
                   + h_lambda_dphi_subtracted_0_80.GetBinError(16)**2))

    fit_function_0_80.SetParameter(12, ue_avg_0_80)


    h_lambda_dphi_subtracted_with_fit_0_80 = h_lambda_dphi_subtracted_0_80.Clone("h_lambda_dphi_subtracted_with_fit_0_80")
    fit_result_0_80 = h_lambda_dphi_subtracted_with_fit_0_80.Fit(fit_function_0_80, "R")


    ue_avg_fit_0_80 = rt.TF1("ue_avg_fit_0_80", "pol0", -2, 6)
    ue_avg_fit_0_80.SetParameter(0, ue_avg_0_80)
    ue_avg_fit_0_80.SetParError(0, ue_avg_error_0_80)

elif USE_V2:
    ue_avg_0_80 = (h_lambda_dphi_subtracted_0_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(16))/6

    v2_fit_0_80 = rt.TF1("v2_fit_0_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    ue_avg_0_80 = (h_lambda_dphi_subtracted_0_80.GetBinContent(1) 
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(2)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(7)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(8)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(9)
                   + h_lambda_dphi_subtracted_0_80.GetBinContent(16))/6

    v2_fit_0_80 = rt.TF1("v2_fit_0_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"

    von_fit_0_80 = rt.TF1("von_fit_0_80", von_fit_string, -2, 6)

    if PT_MODE == 0:
        von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][0])
        von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][2])
        von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][3])
        von_fit_0_80.SetParLimits(3, 0.015/3, 0.015*3)
        von_fit_0_80.SetParameter(3, 0.015)
        von_fit_0_80.SetParLimits(4, 7.6/3, 7.6*3)
        von_fit_0_80.SetParameter(4, 7.6)
        von_fit_0_80.SetParLimits(5, 0.014/3, 0.014*3)
        von_fit_0_80.SetParameter(5, 0.014)
        von_fit_0_80.SetParLimits(6, 2.2/3, 2.2*3)
        von_fit_0_80.SetParameter(6, 2.2)
    elif PT_MODE == 1:
        von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][0])
        von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][2])
        von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][3])
        von_fit_0_80.SetParLimits(3, 0.02/3, 0.02*3)
        von_fit_0_80.SetParameter(3, 0.02)
        von_fit_0_80.SetParLimits(4, 4.6/3, 4.6*3)
        von_fit_0_80.SetParameter(4, 4.6)
        von_fit_0_80.SetParLimits(5, 0.057/3, 0.057*3)
        von_fit_0_80.SetParameter(5, 0.057)
        von_fit_0_80.SetParLimits(6, 1/3, 3)
        von_fit_0_80.SetParameter(6, 1)
    elif PT_MODE == 2:
        von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][0])
        von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][2])
        von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][3])
        von_fit_0_80.SetParLimits(3, 0.009/3, 0.009*3)
        von_fit_0_80.SetParameter(3, 0.009)
        von_fit_0_80.SetParLimits(4, 8.88/3, 8.88*3)
        von_fit_0_80.SetParameter(4, 8.88)
        von_fit_0_80.SetParLimits(5, 0.0076/3, 0.0076*3)
        von_fit_0_80.SetParameter(5, 0.0076)
        von_fit_0_80.SetParLimits(6, 2.55/3, 2.55*3)
        von_fit_0_80.SetParameter(6, 2.55)

    h_lambda_dphi_subtracted_with_fit_0_80 = h_lambda_dphi_subtracted_0_80.Clone("h_lambda_dphi_subtracted_with_fit_0_80")
    if PT_MODE == 1: # low pt has away-side width too large for von mises fit, only works if we symmetrize the peaks
        fit_result_0_80 = h_lambda_dphi_subtracted_with_fit_0_80.Fit(von_fit_0_80, "R", "", 0, rt.TMath.Pi())
    else:
        fit_result_0_80 = h_lambda_dphi_subtracted_with_fit_0_80.Fit(von_fit_0_80, "R")


    if PT_MODE == 0:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][0])
        v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][1])
        v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][2])
        v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_lambda_cent_0_80_trigger_4_8_assoc_25_4"][3])


else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_lambda_dphi_subtracted_0_80_zeroed = h_lambda_dphi_subtracted_0_80.Clone("h_lambda_dphi_subtracted_0_80_zeroed")
    h_lambda_dphi_subtracted_0_80_zeroed.Add(ue_avg_fit_0_80, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_lambda_total_integral_error_0_80 = fit_function_0_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_0_80 = fit_function_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_0_80 = ue_avg_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_0_80 = ue_avg_error_0_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_80_halved = h_lambda_ue_integral_error_0_80/2
    h_lambda_near_integral_error_0_80 = math.sqrt((fit_function_0_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_near_integral_0_80 = fit_function_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_0_80 = math.sqrt((fit_function_0_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_away_integral_0_80 = fit_function_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - ue_avg_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_lambda_total_integral_error_0_80 = von_fit_0_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_total_integral_0_80 = von_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_0_80 = v2_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_ue_integral_error_0_80 = v2_fit_0_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_80_halved = h_lambda_ue_integral_error_0_80/2
    h_lambda_near_integral_error_0_80 = math.sqrt((von_fit_0_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_near_integral_0_80 = von_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_away_integral_error_0_80 = math.sqrt((von_fit_0_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_away_integral_0_80 = von_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - v2_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_lambda_dphi_subtracted_0_80_zeroed = h_lambda_dphi_subtracted_0_80.Clone("h_lambda_dphi_subtracted_0_80_zeroed")
    h_lambda_dphi_subtracted_0_80_zeroed.Add(v2_fit_0_80, -1)

elif USE_V2:

    DPHI_BINS = h_lambda_dphi_subtracted_0_80.GetNbinsX()
    h_lambda_total_integral_0_80 = 0
    h_lambda_near_integral_0_80 = 0
    h_lambda_away_integral_0_80 = 0
    h_lambda_ue_integral_0_80 = v2_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

    h_lambda_total_integral_error_0_80 = 0
    h_lambda_near_integral_error_0_80 = 0
    h_lambda_away_integral_error_0_80 = 0
    h_lambda_ue_integral_error_0_80 = v2_fit_0_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_80_halved = h_lambda_ue_integral_error_0_80/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_0_80 += h_lambda_dphi_subtracted_0_80.GetBinContent(bin_num)
        h_lambda_total_integral_error_0_80 += (h_lambda_dphi_subtracted_0_80.GetBinError(bin_num))**2
        part = h_lambda_dphi_subtracted_0_80.GetBinContent(bin_num) - v2_fit_0_80.Eval(h_lambda_dphi_subtracted_0_80.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_0_80 += part
            h_lambda_near_integral_error_0_80 += (h_lambda_dphi_subtracted_0_80.GetBinError(bin_num))**2
        else:
            h_lambda_away_integral_0_80 += part
            h_lambda_away_integral_error_0_80 += (h_lambda_dphi_subtracted_0_80.GetBinError(bin_num))**2
        
    h_lambda_total_integral_error_0_80 = math.sqrt(h_lambda_total_integral_error_0_80)
    h_lambda_near_integral_error_0_80 = math.sqrt(h_lambda_near_integral_error_0_80)
    h_lambda_near_integral_error_0_80 = math.sqrt(h_lambda_near_integral_error_0_80**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_away_integral_error_0_80 = math.sqrt(h_lambda_away_integral_error_0_80)
    h_lambda_away_integral_error_0_80 = math.sqrt(h_lambda_away_integral_error_0_80**2 + h_lambda_ue_integral_error_0_80_halved**2)

else: 

    h_lambda_dphi_subtracted_0_80_zeroed = h_lambda_dphi_subtracted_0_80.Clone("h_lambda_dphi_subtracted_0_80_zeroed")
    h_lambda_dphi_subtracted_0_80_zeroed.Add(ue_line_0_80, -1)

    DPHI_BINS = h_lambda_dphi_subtracted_0_80.GetNbinsX()

    h_lambda_total_integral_0_80 = 0
    h_lambda_near_integral_0_80 = 0
    h_lambda_away_integral_0_80 = 0
    h_lambda_ue_integral_0_80 = ue_avg_0_80*DPHI_BINS

    h_lambda_total_integral_error_0_80 = 0
    h_lambda_near_integral_error_0_80 = 0
    h_lambda_away_integral_error_0_80 = 0
    h_lambda_ue_integral_error_0_80 = ue_avg_error_0_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_lambda_ue_integral_error_0_80_halved = h_lambda_ue_integral_error_0_80/2

    for bin_num in range(1, DPHI_BINS + 1):
        h_lambda_total_integral_0_80 += h_lambda_dphi_subtracted_0_80.GetBinContent(bin_num)
        h_lambda_total_integral_error_0_80 += h_lambda_dphi_subtracted_0_80.GetBinError(bin_num)**2
        part = h_lambda_dphi_subtracted_0_80.GetBinContent(bin_num) - ue_avg_0_80
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_lambda_near_integral_0_80 += part
            h_lambda_near_integral_error_0_80 += h_lambda_dphi_subtracted_0_80.GetBinError(bin_num)**2
        else:
            h_lambda_away_integral_0_80 += part
            h_lambda_away_integral_error_0_80 += h_lambda_dphi_subtracted_0_80.GetBinError(bin_num)**2

    h_lambda_total_integral_error_0_80 = math.sqrt(h_lambda_total_integral_error_0_80)
    h_lambda_near_integral_error_0_80 = math.sqrt(h_lambda_near_integral_error_0_80)
    h_lambda_near_integral_error_0_80 = math.sqrt(h_lambda_near_integral_error_0_80**2 + h_lambda_ue_integral_error_0_80_halved**2)
    h_lambda_away_integral_error_0_80 = math.sqrt(h_lambda_away_integral_error_0_80)
    h_lambda_away_integral_error_0_80 = math.sqrt(h_lambda_away_integral_error_0_80**2 + h_lambda_ue_integral_error_0_80_halved**2)


h_h_dphi_0_80 = h_h_2d_mixcor_0_80.ProjectionY("h_h_dphi_0_80")

if USE_AVG_4:
    hh_ue_line_0_80 = rt.TF1("hh_ue_line_0_80", "pol0", -2, 6)
    hh_ue_upper_line_0_80 = rt.TF1("hh_ue_upper_line_0_80", "pol0", -2, 6)
    hh_ue_lower_line_0_80 = rt.TF1("hh_ue_lower_line_0_80", "pol0", -2, 6)
    hh_zero_line_0_80 = rt.TF1("hh_zero_line_0_80", "pol0", -2, 6)
    hh_zero_upper_line_0_80 = rt.TF1("hh_zero_upper_line_0_80", "pol0", -2, 6)
    hh_zero_lower_line_0_80 = rt.TF1("hh_zero_lower_line_0_80", "pol0", -2, 6)
    hh_ue_avg_0_80 = (h_h_dphi_0_80.GetBinContent(1) 
                   + h_h_dphi_0_80.GetBinContent(8)
                   + h_h_dphi_0_80.GetBinContent(9)
                   + h_h_dphi_0_80.GetBinContent(16))/4

    hh_ue_avg_error_0_80 = (1/4)*(math.sqrt(h_h_dphi_0_80.GetBinError(1)**2 
                   + h_h_dphi_0_80.GetBinError(8)**2
                   + h_h_dphi_0_80.GetBinError(9)**2
                   + h_h_dphi_0_80.GetBinError(16)**2))


    hh_ue_line_0_80.SetParameter(0, hh_ue_avg_0_80)
    hh_ue_upper_line_0_80.SetParameter(0, hh_ue_avg_0_80 + hh_ue_avg_error_0_80)
    hh_ue_lower_line_0_80.SetParameter(0, hh_ue_avg_0_80 - hh_ue_avg_error_0_80)

    hh_zero_line_0_80.SetParameter(0, 0)
    hh_zero_upper_line_0_80.SetParameter(0, hh_ue_avg_error_0_80)
    hh_zero_lower_line_0_80.SetParameter(0, -hh_ue_avg_error_0_80)

elif USE_AVG_6 or USE_AVG_6_NONNEGATIVE:
    hh_ue_line_0_80 = rt.TF1("hh_ue_line_0_80", "pol0", -2, 6)
    hh_ue_upper_line_0_80 = rt.TF1("hh_ue_upper_line_0_80", "pol0", -2, 6)
    hh_ue_lower_line_0_80 = rt.TF1("hh_ue_lower_line_0_80", "pol0", -2, 6)
    hh_zero_line_0_80 = rt.TF1("hh_zero_line_0_80", "pol0", -2, 6)
    hh_zero_upper_line_0_80 = rt.TF1("hh_zero_upper_line_0_80", "pol0", -2, 6)
    hh_zero_lower_line_0_80 = rt.TF1("hh_zero_lower_line_0_80", "pol0", -2, 6)
    hh_ue_avg_0_80 = (h_h_dphi_0_80.GetBinContent(1) 
                   + h_h_dphi_0_80.GetBinContent(2)
                   + h_h_dphi_0_80.GetBinContent(7)
                   + h_h_dphi_0_80.GetBinContent(8)
                   + h_h_dphi_0_80.GetBinContent(9)
                   + h_h_dphi_0_80.GetBinContent(16))/6

    hh_ue_avg_error_0_80 = (1/6)*(math.sqrt(h_h_dphi_0_80.GetBinError(1)**2 
                   + h_h_dphi_0_80.GetBinError(2)**2
                   + h_h_dphi_0_80.GetBinError(7)**2
                   + h_h_dphi_0_80.GetBinError(8)**2
                   + h_h_dphi_0_80.GetBinError(9)**2
                   + h_h_dphi_0_80.GetBinError(16)**2))

    hh_ue_line_0_80.SetParameter(0, hh_ue_avg_0_80)
    hh_ue_upper_line_0_80.SetParameter(0, hh_ue_avg_0_80 + hh_ue_avg_error_0_80)
    hh_ue_lower_line_0_80.SetParameter(0, hh_ue_avg_0_80 - hh_ue_avg_error_0_80)

    hh_zero_line_0_80.SetParameter(0, 0)
    hh_zero_upper_line_0_80.SetParameter(0, hh_ue_avg_error_0_80)
    hh_zero_lower_line_0_80.SetParameter(0, -hh_ue_avg_error_0_80)

elif USE_ZYAM:
    hh_ue_line_0_80 = rt.TF1("hh_ue_line_0_80", "pol0", -2, 6)
    hh_ue_upper_line_0_80 = rt.TF1("hh_ue_upper_line_0_80", "pol0", -2, 6)
    hh_ue_lower_line_0_80 = rt.TF1("hh_ue_lower_line_0_80", "pol0", -2, 6)
    hh_zero_line_0_80 = rt.TF1("hh_zero_line_0_80", "pol0", -2, 6)
    hh_zero_upper_line_0_80 = rt.TF1("hh_zero_upper_line_0_80", "pol0", -2, 6)
    hh_zero_lower_line_0_80 = rt.TF1("hh_zero_lower_line_0_80", "pol0", -2, 6)
    
    min_bin = h_h_dphi_0_80.GetMinimumBin()
    hh_ue_avg_0_80 = h_h_dphi_0_80.GetBinContent(min_bin)
    hh_ue_avg_error_0_80 = h_h_dphi_0_80.GetBinError(min_bin)

    hh_ue_line_0_80.SetParameter(0, hh_ue_avg_0_80)
    hh_ue_upper_line_0_80.SetParameter(0, hh_ue_avg_0_80 + hh_ue_avg_error_0_80)
    hh_ue_lower_line_0_80.SetParameter(0, hh_ue_avg_0_80 - hh_ue_avg_error_0_80)

    hh_zero_line_0_80.SetParameter(0, 0)
    hh_zero_upper_line_0_80.SetParameter(0, hh_ue_avg_error_0_80)
    hh_zero_lower_line_0_80.SetParameter(0, -hh_ue_avg_error_0_80)
elif USE_FIT:

    # fitting to v2 component + four gaussians
    hh_fit_string_0_80 = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_fit_string_0_80 += " + gaus(3) + gaus(6) + gaus(9) + gaus(12)"
    hh_fit_function_0_80 = rt.TF1("hh_fit_function_0_80", hh_fit_string_0_80, -2, 6)
    hh_fit_function_0_80.SetNpx(1000)
    if PT_MODE == 0:
        hh_fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][0])
        hh_fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][2])
        hh_fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][0])
        hh_fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][2])
        hh_fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_fit_function_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][0])
        hh_fit_function_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][2])
        hh_fit_function_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][3])

    # setting parameters for first gaussian (centered at 0)
    hh_fit_function_0_80.SetParLimits(3, 0.02, 1)
    hh_fit_function_0_80.SetParameter(3, 0.16)
    hh_fit_function_0_80.FixParameter(4, 0)
    hh_fit_function_0_80.SetParLimits(5, 0.1, 1)
    hh_fit_function_0_80.SetParameter(5, 0.3)
    # setting parameters for second gaussian (centered at 0 + 2pi)
    hh_fit_function_0_80.SetParLimits(6, 0, 1)
    hh_fit_function_0_80.SetParameter(6, 0)
    hh_fit_function_0_80.FixParameter(7, 2*rt.TMath.Pi())
    hh_fit_function_0_80.SetParLimits(8, 0.1, 2)
    hh_fit_function_0_80.SetParameter(8, 0.1)
    # setting parameters for third gaussian (centered at pi)
    hh_fit_function_0_80.SetParLimits(9, 0.02, 1)
    hh_fit_function_0_80.SetParameter(9, 0.3)
    hh_fit_function_0_80.FixParameter(10, rt.TMath.Pi())
    hh_fit_function_0_80.SetParLimits(11, 0.3, 2)
    hh_fit_function_0_80.SetParameter(11, 0.7)
    # setting parameters for fourth gaussian (centered at pi + 2pi)
    hh_fit_function_0_80.SetParLimits(12, 0, 1)
    hh_fit_function_0_80.SetParameter(12, 0)
    hh_fit_function_0_80.FixParameter(13, rt.TMath.Pi() - 2*rt.TMath.Pi())
    hh_fit_function_0_80.SetParLimits(14, 0.1, 3)
    hh_fit_function_0_80.SetParameter(14, 0.1)


    # setting parameters for constant background
    hh_ue_avg_0_80 = (h_h_dphi_0_80.GetBinContent(1) 
                + h_h_dphi_0_80.GetBinContent(2)
                + h_h_dphi_0_80.GetBinContent(7)
                + h_h_dphi_0_80.GetBinContent(8)
                + h_h_dphi_0_80.GetBinContent(9)
                + h_h_dphi_0_80.GetBinContent(16))/6

    hh_ue_avg_error_0_80 = (1/6)*(math.sqrt(h_h_dphi_0_80.GetBinError(1)**2 
                   + h_h_dphi_0_80.GetBinError(2)**2
                   + h_h_dphi_0_80.GetBinError(7)**2
                   + h_h_dphi_0_80.GetBinError(8)**2
                   + h_h_dphi_0_80.GetBinError(9)**2
                   + h_h_dphi_0_80.GetBinError(16)**2))

    hh_fit_function_0_80.SetParameter(12, hh_ue_avg_0_80)


    h_h_dphi_with_fit_0_80 = h_h_dphi_0_80.Clone("h_h_dphi_with_fit_0_80")
    h_h_dphi_with_fit_0_80.Fit(hh_fit_function_0_80, "R")

    hh_ue_avg_fit_0_80 = rt.TF1("hh_ue_avg_fit_0_80", "pol0", -2, 6)
    hh_ue_avg_fit_0_80.SetParameter(0, hh_ue_avg_0_80)
    hh_ue_avg_fit_0_80.SetParError(0, hh_ue_avg_error_0_80)

elif USE_V2:
    hh_ue_avg_0_80 = (h_h_dphi_0_80.GetBinContent(1) 
                   + h_h_dphi_0_80.GetBinContent(2)
                   + h_h_dphi_0_80.GetBinContent(7)
                   + h_h_dphi_0_80.GetBinContent(8)
                   + h_h_dphi_0_80.GetBinContent(9)
                   + h_h_dphi_0_80.GetBinContent(16))/6

    hh_v2_fit_0_80 = rt.TF1("hh_v2_fit_0_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)

    if PT_MODE == 0:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][3])


elif USE_VON:
    hh_ue_avg_0_80 = (h_h_dphi_0_80.GetBinContent(1) 
                   + h_h_dphi_0_80.GetBinContent(2)
                   + h_h_dphi_0_80.GetBinContent(7)
                   + h_h_dphi_0_80.GetBinContent(8)
                   + h_h_dphi_0_80.GetBinContent(9)
                   + h_h_dphi_0_80.GetBinContent(16))/6

    hh_v2_fit_0_80 = rt.TF1("hh_v2_fit_0_80", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -2, 6)
    hh_von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
    hh_von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
    hh_von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
    hh_von_fit_0_80 = rt.TF1("hh_von_fit_0_80", hh_von_fit_string, -2, 6)

    if PT_MODE == 0:
        hh_von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][0])
        hh_von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][2])
        hh_von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][0])
        hh_von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][2])
        hh_von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_von_fit_0_80.FixParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][0])
        hh_von_fit_0_80.FixParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][2])
        hh_von_fit_0_80.FixParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][3])

    hh_von_fit_0_80.SetParLimits(3, 0, 1)
    hh_von_fit_0_80.SetParameter(3, 0.02)
    hh_von_fit_0_80.SetParLimits(4, 0, 100)
    hh_von_fit_0_80.SetParameter(4, 1)
    hh_von_fit_0_80.SetParLimits(5, 0, 1)
    hh_von_fit_0_80.SetParameter(5, 0.01)
    hh_von_fit_0_80.SetParLimits(6, 0, 100)
    hh_von_fit_0_80.SetParameter(6, 1)

    h_h_dphi_with_fit_0_80 = h_h_dphi_0_80.Clone("h_h_dphi_with_fit_0_80")
    hh_fit_result_0_80 = h_h_dphi_with_fit_0_80.Fit(hh_von_fit_0_80, "R")

    if PT_MODE == 0:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_2_4"][3])
    elif PT_MODE == 1:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_15_25"][3])
    elif PT_MODE == 2:
        hh_v2_fit_0_80.SetParameter(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][0])
        hh_v2_fit_0_80.SetParError(0, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][1])
        hh_v2_fit_0_80.SetParameter(1, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][2])
        hh_v2_fit_0_80.SetParameter(2, v2_fitpars_dict["h_h_cent_0_80_trigger_4_8_assoc_25_4"][3])

else:
    raise NotImplementedError("UE line mode not supported")


if USE_FIT:
    h_h_dphi_0_80_zeroed = h_h_dphi_0_80.Clone("h_h_dphi_0_80_zeroed")
    h_h_dphi_0_80_zeroed.Add(hh_ue_avg_fit_0_80, -1)


    # we are only allowed to call IntegralError (and have it make sense) if the most recent fit corresponds to the TF1 being integrated

    h_h_total_integral_error_0_80 = hh_fit_function_0_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_0_80 = hh_fit_function_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_0_80 = hh_ue_avg_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_0_80 = hh_ue_avg_error_0_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_80_halved = h_h_ue_integral_error_0_80/2
    h_h_near_integral_error_0_80 = math.sqrt((hh_fit_function_0_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_near_integral_0_80 = hh_fit_function_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_0_80 = math.sqrt((hh_fit_function_0_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_away_integral_0_80 = hh_fit_function_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_ue_avg_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH

elif USE_VON:

    h_h_total_integral_error_0_80 = hh_von_fit_0_80.IntegralError(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_0_80 = hh_von_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_0_80 = hh_v2_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_ue_integral_error_0_80 = hh_v2_fit_0_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_80_halved = h_h_ue_integral_error_0_80/2
    h_h_near_integral_error_0_80 = math.sqrt((hh_von_fit_0_80.IntegralError(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_near_integral_0_80 = hh_von_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_0_80.Integral(-rt.TMath.Pi()/2, rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_away_integral_error_0_80 = math.sqrt((hh_von_fit_0_80.IntegralError(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH)**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_away_integral_0_80 = hh_von_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH - hh_v2_fit_0_80.Integral(rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_dphi_0_80_zeroed = h_h_dphi_0_80.Clone("h_h_dphi_0_80_zeroed")
    h_h_dphi_0_80_zeroed.Add(hh_v2_fit_0_80, -1)

elif USE_V2:

    DPHI_BINS = h_h_dphi_0_80.GetNbinsX()
    h_h_total_integral_0_80 = 0
    h_h_near_integral_0_80 = 0
    h_h_away_integral_0_80 = 0
    h_h_ue_integral_0_80 = hh_v2_fit_0_80.Integral(-rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)/BIN_WIDTH
    h_h_total_integral_error_0_80 = 0
    h_h_near_integral_error_0_80 = 0
    h_h_away_integral_error_0_80 = 0
    h_h_ue_integral_error_0_80 = hh_v2_fit_0_80.GetParError(0)*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_80_halved = h_h_ue_integral_error_0_80/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_0_80 += h_h_dphi_0_80.GetBinContent(bin_num)
        h_h_total_integral_error_0_80 += (h_h_dphi_0_80.GetBinError(bin_num))**2
        part = h_h_dphi_0_80.GetBinContent(bin_num) - hh_v2_fit_0_80.Eval(h_h_dphi_0_80.GetBinCenter(bin_num))
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_0_80 += part
            h_h_near_integral_error_0_80 += (h_h_dphi_0_80.GetBinError(bin_num))**2
        else:
            h_h_away_integral_0_80 += part
            h_h_away_integral_error_0_80 += (h_h_dphi_0_80.GetBinError(bin_num))**2
    h_h_total_integral_error_0_80 = math.sqrt(h_h_total_integral_error_0_80)
    h_h_near_integral_error_0_80 = math.sqrt(h_h_near_integral_error_0_80)
    h_h_near_integral_error_0_80 = math.sqrt(h_h_near_integral_error_0_80**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_away_integral_error_0_80 = math.sqrt(h_h_away_integral_error_0_80)
    h_h_away_integral_error_0_80 = math.sqrt(h_h_away_integral_error_0_80**2 + h_h_ue_integral_error_0_80_halved**2)

else: 

    h_h_dphi_0_80_zeroed = h_h_dphi_0_80.Clone("h_h_dphi_0_80_zeroed")
    h_h_dphi_0_80_zeroed.Add(hh_ue_line_0_80, -1)

    DPHI_BINS = h_h_dphi_0_80.GetNbinsX()
    h_h_total_integral_0_80 = 0
    h_h_near_integral_0_80 = 0
    h_h_away_integral_0_80 = 0
    h_h_ue_integral_0_80 = hh_ue_avg_0_80*DPHI_BINS
    h_h_total_integral_error_0_80 = 0
    h_h_near_integral_error_0_80 = 0
    h_h_away_integral_error_0_80 = 0
    h_h_ue_integral_error_0_80 = hh_ue_avg_error_0_80*(2*rt.TMath.Pi())/BIN_WIDTH
    h_h_ue_integral_error_0_80_halved = h_h_ue_integral_error_0_80/2
    for bin_num in range(1, DPHI_BINS + 1):
        h_h_total_integral_0_80 += h_h_dphi_0_80.GetBinContent(bin_num)
        h_h_total_integral_error_0_80 += h_h_dphi_0_80.GetBinError(bin_num)**2
        part = h_h_dphi_0_80.GetBinContent(bin_num) - hh_ue_avg_0_80
        if part < 0 and USE_AVG_6_NONNEGATIVE:
            part = 0
            continue
        if bin_num < 9:
            h_h_near_integral_0_80 += part
            h_h_near_integral_error_0_80 += h_h_dphi_0_80.GetBinError(bin_num)**2
        else:
            h_h_away_integral_0_80 += part
            h_h_away_integral_error_0_80 += h_h_dphi_0_80.GetBinError(bin_num)**2
    h_h_total_integral_error_0_80 = math.sqrt(h_h_total_integral_error_0_80)
    h_h_near_integral_error_0_80 = math.sqrt(h_h_near_integral_error_0_80)
    h_h_near_integral_error_0_80 = math.sqrt(h_h_near_integral_error_0_80**2 + h_h_ue_integral_error_0_80_halved**2)
    h_h_away_integral_error_0_80 = math.sqrt(h_h_away_integral_error_0_80)
    h_h_away_integral_error_0_80 = math.sqrt(h_h_away_integral_error_0_80**2 + h_h_ue_integral_error_0_80_halved**2)


output_file.cd()
h_lambda_2d_subtracted_0_80.Write()

h_lambda_dphi_subtracted_0_80.Write()
h_h_dphi_0_80.Write()



if USE_FIT:
    fit_function_0_80.Write()
    ue_avg_fit_0_80.Write()
    hh_fit_function_0_80.Write()
    hh_ue_avg_fit_0_80.Write()
    h_lambda_dphi_subtracted_with_fit_0_80.Write()
    h_h_dphi_with_fit_0_80.Write()
    h_lambda_dphi_subtracted_0_80_zeroed.Write()
    h_h_dphi_0_80_zeroed.Write()
elif USE_V2:
    v2_fit_0_80.Write()
    hh_v2_fit_0_80.Write()
elif USE_VON:
    von_fit_0_80.Write()
    hh_von_fit_0_80.Write()
    v2_fit_0_80.Write()
    hh_v2_fit_0_80.Write()
    h_lambda_dphi_subtracted_with_fit_0_80.Write()
    h_h_dphi_with_fit_0_80.Write()
    h_lambda_dphi_subtracted_0_80_zeroed.Write()
    h_h_dphi_0_80_zeroed.Write()
else:
    ue_upper_line_0_80.Write()
    ue_line_0_80.Write()
    ue_lower_line_0_80.Write()
    zero_upper_line_0_80.Write()
    zero_line_0_80.Write()
    zero_lower_line_0_80.Write()

    hh_ue_upper_line_0_80.Write()
    hh_ue_line_0_80.Write()
    hh_ue_lower_line_0_80.Write()
    hh_zero_upper_line_0_80.Write()
    hh_zero_line_0_80.Write()
    hh_zero_lower_line_0_80.Write()

    h_lambda_dphi_subtracted_0_80_zeroed.Write()
    h_h_dphi_0_80_zeroed.Write()


near_ratio_0_80 = h_lambda_near_integral_0_80/h_h_near_integral_0_80
away_ratio_0_80 = h_lambda_away_integral_0_80/h_h_away_integral_0_80
ue_ratio_0_80 = h_lambda_ue_integral_0_80/h_h_ue_integral_0_80
total_ratio_0_80 = h_lambda_total_integral_0_80/h_h_total_integral_0_80

near_ratio_error_0_80 = near_ratio_0_80*math.sqrt((h_lambda_near_integral_error_0_80/h_lambda_near_integral_0_80)**2
                                                 + (h_h_near_integral_error_0_80/h_h_near_integral_0_80)**2)
away_ratio_error_0_80 = away_ratio_0_80*math.sqrt((h_lambda_away_integral_error_0_80/h_lambda_away_integral_0_80)**2
                                                 + (h_h_away_integral_error_0_80/h_h_away_integral_0_80)**2)
ue_ratio_error_0_80 = ue_ratio_0_80*math.sqrt((h_lambda_ue_integral_error_0_80/h_lambda_ue_integral_0_80)**2
                                                 + (h_h_ue_integral_error_0_80/h_h_ue_integral_0_80)**2)
total_ratio_error_0_80 = total_ratio_0_80*math.sqrt((h_lambda_total_integral_error_0_80/h_lambda_total_integral_0_80)**2
                                                 + (h_h_total_integral_error_0_80/h_h_total_integral_0_80)**2)

############################################################################################################
############################################################################################################
############################################ FINAL RATIOS ##################################################
############################################################################################################
############################################################################################################

mult_list = arr.array('d', [35, 65, 90])
mult_error_list = arr.array('d', [15, 15, 10])

mult_list_minbias = arr.array('d', [40])
mult_error_list_minbias = arr.array('d', [40])


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

h_lambda_near_list_minbias = arr.array('d', [h_lambda_near_integral_0_80])
h_lambda_near_error_list_minbias = arr.array('d', [h_lambda_near_integral_error_0_80])
h_lambda_away_list_minbias = arr.array('d', [h_lambda_away_integral_0_80])
h_lambda_away_error_list_minbias = arr.array('d', [h_lambda_away_integral_error_0_80])
h_lambda_ue_list_minbias = arr.array('d', [h_lambda_ue_integral_0_80])
h_lambda_ue_error_list_minbias = arr.array('d', [h_lambda_ue_integral_error_0_80])
h_lambda_total_list_minbias = arr.array('d', [h_lambda_total_integral_0_80])
h_lambda_total_error_list_minbias = arr.array('d', [h_lambda_total_integral_error_0_80])

h_h_near_list_minbias = arr.array('d', [h_h_near_integral_0_80])
h_h_near_error_list_minbias = arr.array('d', [h_h_near_integral_error_0_80])
h_h_away_list_minbias = arr.array('d', [h_h_away_integral_0_80])
h_h_away_error_list_minbias = arr.array('d', [h_h_away_integral_error_0_80])
h_h_ue_list_minbias = arr.array('d', [h_h_ue_integral_0_80])
h_h_ue_error_list_minbias = arr.array('d', [h_h_ue_integral_error_0_80])
h_h_total_list_minbias = arr.array('d', [h_h_total_integral_0_80])
h_h_total_error_list_minbias = arr.array('d', [h_h_total_integral_error_0_80])


near_ratio_list = arr.array('d', [near_ratio_50_80, near_ratio_20_50, near_ratio_0_20])
near_ratio_error_list = arr.array('d', [near_ratio_error_0_20, near_ratio_error_20_50, near_ratio_error_50_80])
away_ratio_list = arr.array('d', [away_ratio_50_80, away_ratio_20_50, away_ratio_0_20])
away_ratio_error_list = arr.array('d', [away_ratio_error_50_80, away_ratio_error_20_50, away_ratio_error_0_20])
ue_ratio_list = arr.array('d', [ue_ratio_50_80, ue_ratio_20_50, ue_ratio_0_20])
ue_ratio_error_list = arr.array('d', [ue_ratio_error_50_80, ue_ratio_error_20_50, ue_ratio_error_0_20])
total_ratio_list = arr.array('d', [total_ratio_50_80, total_ratio_20_50, total_ratio_0_20])
total_ratio_error_list = arr.array('d', [total_ratio_error_50_80, total_ratio_error_20_50, total_ratio_error_0_20])

print(near_ratio_list, "near ratio list")
print(away_ratio_list, "away ratio list")

near_ratio_list_minbias = arr.array('d', [near_ratio_0_80])
near_ratio_error_list_minbias = arr.array('d', [near_ratio_error_0_80])
away_ratio_list_minbias = arr.array('d', [away_ratio_0_80])
away_ratio_error_list_minbias = arr.array('d', [away_ratio_error_0_80])
ue_ratio_list_minbias = arr.array('d', [ue_ratio_0_80])
ue_ratio_error_list_minbias = arr.array('d', [ue_ratio_error_0_80])
total_ratio_list_minbias = arr.array('d', [total_ratio_0_80])
total_ratio_error_list_minbias = arr.array('d', [total_ratio_error_0_80])


h_lambda_near_graph = rt.TGraphErrors(3, mult_list, h_lambda_near_list, mult_error_list, h_lambda_near_error_list)
h_lambda_away_graph = rt.TGraphErrors(3, mult_list, h_lambda_away_list, mult_error_list, h_lambda_away_error_list)
h_lambda_ue_graph = rt.TGraphErrors(3, mult_list, h_lambda_ue_list, mult_error_list, h_lambda_ue_error_list)
h_lambda_total_graph = rt.TGraphErrors(3, mult_list, h_lambda_total_list, mult_error_list, h_lambda_total_error_list)
h_h_near_graph = rt.TGraphErrors(3, mult_list, h_h_near_list, mult_error_list, h_h_near_error_list)
h_h_away_graph = rt.TGraphErrors(3, mult_list, h_h_away_list, mult_error_list, h_h_away_error_list)
h_h_ue_graph = rt.TGraphErrors(3, mult_list, h_h_ue_list, mult_error_list, h_h_ue_error_list)
h_h_total_graph = rt.TGraphErrors(3, mult_list, h_h_total_list, mult_error_list, h_h_total_error_list)

h_lambda_near_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_lambda_near_list_minbias, mult_error_list_minbias, h_lambda_near_error_list_minbias)
h_lambda_away_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_lambda_away_list_minbias, mult_error_list_minbias, h_lambda_away_error_list_minbias)
h_lambda_ue_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_lambda_ue_list_minbias, mult_error_list_minbias, h_lambda_ue_error_list_minbias)
h_lambda_total_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_lambda_total_list_minbias, mult_error_list_minbias, h_lambda_total_error_list_minbias)
h_h_near_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_h_near_list_minbias, mult_error_list_minbias, h_h_near_error_list_minbias)
h_h_away_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_h_away_list_minbias, mult_error_list_minbias, h_h_away_error_list_minbias)
h_h_ue_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_h_ue_list_minbias, mult_error_list_minbias, h_h_ue_error_list_minbias)
h_h_total_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, h_h_total_list_minbias, mult_error_list_minbias, h_h_total_error_list_minbias)


near_ratio_graph = rt.TGraphErrors(3, mult_list, near_ratio_list, mult_error_list, near_ratio_error_list)
away_ratio_graph = rt.TGraphErrors(3, mult_list, away_ratio_list, mult_error_list, away_ratio_error_list)
ue_ratio_graph = rt.TGraphErrors(3, mult_list, ue_ratio_list, mult_error_list, ue_ratio_error_list)
total_ratio_graph = rt.TGraphErrors(3, mult_list, total_ratio_list, mult_error_list, total_ratio_error_list)


near_ratio_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, near_ratio_list_minbias, mult_error_list_minbias, near_ratio_error_list_minbias)
away_ratio_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, away_ratio_list_minbias, mult_error_list_minbias, away_ratio_error_list_minbias)
ue_ratio_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, ue_ratio_list_minbias, mult_error_list_minbias, ue_ratio_error_list_minbias)
total_ratio_graph_minbias = rt.TGraphErrors(1, mult_list_minbias, total_ratio_list_minbias, mult_error_list_minbias, total_ratio_error_list_minbias)

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

h_lambda_near_graph_minbias.Write("h_lambda_near_graph_minbias")
h_lambda_away_graph_minbias.Write("h_lambda_away_graph_minbias")
h_lambda_ue_graph_minbias.Write("h_lambda_ue_graph_minbias")
h_lambda_total_graph_minbias.Write("h_lambda_total_graph_minbias")
h_h_near_graph_minbias.Write("h_h_near_graph_minbias")
h_h_away_graph_minbias.Write("h_h_away_graph_minbias")
h_h_ue_graph_minbias.Write("h_h_ue_graph_minbias")
h_h_total_graph_minbias.Write("h_h_total_graph_minbias")
near_ratio_graph_minbias.Write("near_ratio_graph_minbias")
away_ratio_graph_minbias.Write("away_ratio_graph_minbias")
ue_ratio_graph_minbias.Write("ue_ratio_graph_minbias")
total_ratio_graph_minbias.Write("total_ratio_graph_minbias")

output_file.Close()
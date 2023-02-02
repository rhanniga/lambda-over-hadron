import math

import array as arr
import ROOT as rt

def get_rms(ratios):
    rms = 0
    for ratio in ratios:
        rms += ratio**2
    rms = math.sqrt(rms/len(ratios))
    return rms

def ratio_error(X, Y, X_err, Y_err):
    return math.sqrt((X_err/X)**2 + (Y_err/Y)**2) * (X/Y)


LOW_PT = False
HIGH_PT = False
NORMAL_PT = True


CENTRAL_TECHNIQUE = "6 bin avg"
N_DPHI_BINS = 16

PRINT_YIELDS = False
PRINT_SYSTEMATICS = False

total_dphi_sys_0_20 = math.sqrt(0.01**2 + 0.001**2 + 0.011**2 + 0.052**2)
total_dphi_sys_20_50 = math.sqrt(0.009**2 + 0.003**2 + 0.018**2 + 0.052**2)
total_dphi_sys_50_80 = math.sqrt(0.019**2 + 0.008**2 + 0.037**2 + 0.052**2)

hh_total_dphi_sys_0_20 = 0.051
hh_total_dphi_sys_20_50 = 0.051
hh_total_dphi_sys_50_80 = 0.051


if LOW_PT:
    avg6_file = rt.TFile("output/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    avg4_file = rt.TFile("output/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    avg6nonneg_file = rt.TFile("output/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    fullfit_file = rt.TFile("output/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    v2_file = rt.TFile("output/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    zyam_file = rt.TFile("output/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    von_file = rt.TFile("output/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal.root")
    file_dict = {
        "6 bin avg": avg6_file,
        "4 bin avg": avg4_file,
        "6 bin avg nonneg": avg6nonneg_file,
        "full fit": fullfit_file,
        "v2": v2_file,
        "zyam": zyam_file,
        "von": von_file
    }

elif HIGH_PT:
    avg6_file = rt.TFile("output/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    avg4_file = rt.TFile("output/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    avg6nonneg_file = rt.TFile("output/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    fullfit_file = rt.TFile("output/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    v2_file = rt.TFile("output/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    zyam_file = rt.TFile("output/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    von_file = rt.TFile("output/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal.root")
    file_dict = {
        "6 bin avg": avg6_file,
        "4 bin avg": avg4_file,
        "6 bin avg nonneg": avg6nonneg_file,
        "full fit": fullfit_file,
        "v2": v2_file,
        "zyam": zyam_file,
        "von": von_file
    }

elif NORMAL_PT:
    avg6_file = rt.TFile("output/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    avg4_file = rt.TFile("output/v0_avg4_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    avg6nonneg_file = rt.TFile("output/v0_avg6nonneg_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    fullfit_file = rt.TFile("output/v0_fullfit_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    v2_file = rt.TFile("output/v0_v2_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    zyam_file = rt.TFile("output/v0_zyam_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    von_file = rt.TFile("output/v0_von_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
    file_dict = {
        "6 bin avg": avg6_file,
        "4 bin avg": avg4_file,
        "6 bin avg nonneg": avg6nonneg_file,
        "full fit": fullfit_file,
        "v2": v2_file,
        "zyam": zyam_file,
        "von": von_file
    }
else:
    print("PT RANGE NOT DEFINED, EXITING")
    quit()


h_lambda_near_dict = {}
h_lambda_away_dict = {}
h_lambda_ue_dict = {}
h_lambda_total_dict = {}

h_h_near_dict = {}
h_h_away_dict = {}
h_h_ue_dict = {}
h_h_total_dict = {}

for key in file_dict:
    h_lambda_near_dict[key] = file_dict[key].Get("h_lambda_near_graph")
    h_lambda_away_dict[key] = file_dict[key].Get("h_lambda_away_graph")
    h_lambda_ue_dict[key] = file_dict[key].Get("h_lambda_ue_graph")
    h_lambda_total_dict[key] = file_dict[key].Get("h_lambda_total_graph")

    h_h_near_dict[key] = file_dict[key].Get("h_h_near_graph")
    h_h_away_dict[key] = file_dict[key].Get("h_h_away_graph")
    h_h_ue_dict[key] = file_dict[key].Get("h_h_ue_graph")
    h_h_total_dict[key] = file_dict[key].Get("h_h_total_graph")


BIN_WIDTH = 2*rt.TMath.Pi()/N_DPHI_BINS

for key in file_dict:

    hl_near_graph = h_lambda_near_dict[key]
    hl_away_graph = h_lambda_away_dict[key]
    hl_ue_graph = h_lambda_ue_dict[key]
    hl_total_graph = h_lambda_total_dict[key]

    hh_near_graph = h_h_near_dict[key]
    hh_away_graph = h_h_away_dict[key]
    hh_ue_graph = h_h_ue_dict[key]
    hh_total_graph = h_h_total_dict[key]

    if key == "full fit" or key == "von":
        for i in range(3):

            hl_near_graph.SetPoint(i, hl_near_graph.GetX()[i], hl_near_graph.GetY()[i]/BIN_WIDTH)
            hl_away_graph.SetPoint(i, hl_away_graph.GetX()[i], hl_away_graph.GetY()[i]/BIN_WIDTH)
            hl_ue_graph.SetPoint(i, hl_ue_graph.GetX()[i], hl_ue_graph.GetY()[i]/BIN_WIDTH)
            hl_total_graph.SetPoint(i, hl_total_graph.GetX()[i], hl_total_graph.GetY()[i]/BIN_WIDTH)

            hh_near_graph.SetPoint(i, hh_near_graph.GetX()[i], hh_near_graph.GetY()[i]/BIN_WIDTH)
            hh_away_graph.SetPoint(i, hh_away_graph.GetX()[i], hh_away_graph.GetY()[i]/BIN_WIDTH)
            hh_ue_graph.SetPoint(i, hh_ue_graph.GetX()[i], hh_ue_graph.GetY()[i]/BIN_WIDTH)
            hh_total_graph.SetPoint(i, hh_total_graph.GetX()[i], hh_total_graph.GetY()[i]/BIN_WIDTH)
            
    if key == "v2":
        for i in range(3):
            hl_ue_graph.SetPoint(i, hl_ue_graph.GetX()[i], hl_ue_graph.GetY()[i]/BIN_WIDTH)
            hh_ue_graph.SetPoint(i, hh_ue_graph.GetX()[i], hh_ue_graph.GetY()[i]/BIN_WIDTH)



if PRINT_YIELDS:
    for key in file_dict:
        print(f"-----------{key}---------------")
        print("----h lambda----")

        if(key == CENTRAL_TECHNIQUE):
            print(f"0-20\% & {h_lambda_near_dict[key].GetY()[2]:.2e} $\pm$ {h_lambda_near_dict[key].GetEY()[2]:.2e} & {h_lambda_away_dict[key].GetY()[2]:.2e} $\pm$ {h_lambda_away_dict[key].GetEY()[2]:.2e} & {h_lambda_ue_dict[key].GetY()[2]:.2e} $\pm$ {h_lambda_ue_dict[key].GetEY()[2]:.2e} & {h_lambda_total_dict[key].GetY()[2]:.2e} $\pm$ {h_lambda_total_dict[key].GetEY()[2]:.2e} \\\\")
            print(f"20-50\% & {h_lambda_near_dict[key].GetY()[1]:.2e} $\pm$ {h_lambda_near_dict[key].GetEY()[1]:.2e} & {h_lambda_away_dict[key].GetY()[1]:.2e} $\pm$ {h_lambda_away_dict[key].GetEY()[1]:.2e} & {h_lambda_ue_dict[key].GetY()[1]:.2e} $\pm$ {h_lambda_ue_dict[key].GetEY()[1]:.2e} & {h_lambda_total_dict[key].GetY()[1]:.2e} $\pm$ {h_lambda_total_dict[key].GetEY()[1]:.2e} \\\\")
            print(f"50-80\% & {h_lambda_near_dict[key].GetY()[0]:.2e} $\pm$ {h_lambda_near_dict[key].GetEY()[0]:.2e} & {h_lambda_away_dict[key].GetY()[0]:.2e} $\pm$ {h_lambda_away_dict[key].GetEY()[0]:.2e} & {h_lambda_ue_dict[key].GetY()[0]:.2e} $\pm$ {h_lambda_ue_dict[key].GetEY()[0]:.2e} & {h_lambda_total_dict[key].GetY()[0]:.2e} $\pm$ {h_lambda_total_dict[key].GetEY()[0]:.2e} \\\\")
        else:
            print(f"0-20\% & {h_lambda_near_dict[key].GetY()[2]:.2e}  & {h_lambda_away_dict[key].GetY()[2]:.2e}  & {h_lambda_ue_dict[key].GetY()[2]:.2e} & {h_lambda_total_dict[key].GetY()[2]:.2e} \\\\")
            print(f"20-50\% & {h_lambda_near_dict[key].GetY()[1]:.2e} & {h_lambda_away_dict[key].GetY()[1]:.2e}  & {h_lambda_ue_dict[key].GetY()[1]:.2e} & {h_lambda_total_dict[key].GetY()[1]:.2e} \\\\")
            print(f"50-80\% & {h_lambda_near_dict[key].GetY()[0]:.2e} & {h_lambda_away_dict[key].GetY()[0]:.2e}  & {h_lambda_ue_dict[key].GetY()[0]:.2e} & {h_lambda_total_dict[key].GetY()[0]:.2e} \\\\")
            
        print("----h h----")

        if(key == CENTRAL_TECHNIQUE):
            print(f"0-20\% & {h_h_near_dict[key].GetY()[2]:.2e} $\pm$ {h_h_near_dict[key].GetEY()[2]:.2e} & {h_h_away_dict[key].GetY()[2]:.2e} $\pm$ {h_h_away_dict[key].GetEY()[2]:.2e} & {h_h_ue_dict[key].GetY()[2]:.2e} $\pm$ {h_h_ue_dict[key].GetEY()[2]:.2e} & {h_h_total_dict[key].GetY()[2]:.2e} $\pm$ {h_h_total_dict[key].GetEY()[2]:.2e} \\\\")
            print(f"20-50\% & {h_h_near_dict[key].GetY()[1]:.2e} $\pm$ {h_h_near_dict[key].GetEY()[1]:.2e} & {h_h_away_dict[key].GetY()[1]:.2e} $\pm$ {h_h_away_dict[key].GetEY()[1]:.2e} & {h_h_ue_dict[key].GetY()[1]:.2e} $\pm$ {h_h_ue_dict[key].GetEY()[1]:.2e} & {h_h_total_dict[key].GetY()[1]:.2e} $\pm$ {h_h_total_dict[key].GetEY()[1]:.2e} \\\\")
            print(f"50-80\% & {h_h_near_dict[key].GetY()[0]:.2e} $\pm$ {h_h_near_dict[key].GetEY()[0]:.2e} & {h_h_away_dict[key].GetY()[0]:.2e} $\pm$ {h_h_away_dict[key].GetEY()[0]:.2e} & {h_h_ue_dict[key].GetY()[0]:.2e} $\pm$ {h_h_ue_dict[key].GetEY()[0]:.2e} & {h_h_total_dict[key].GetY()[0]:.2e} $\pm$ {h_h_total_dict[key].GetEY()[0]:.2e} \\\\")
        else:
            print(f"0-20\% & {h_h_near_dict[key].GetY()[2]:.2e}  & {h_h_away_dict[key].GetY()[2]:.2e}  & {h_h_ue_dict[key].GetY()[2]:.2e} & {h_h_total_dict[key].GetY()[2]:.2e} \\\\")
            print(f"20-50\% & {h_h_near_dict[key].GetY()[1]:.2e} & {h_h_away_dict[key].GetY()[1]:.2e}  & {h_h_ue_dict[key].GetY()[1]:.2e} & {h_h_total_dict[key].GetY()[1]:.2e} \\\\")
            print(f"50-80\% & {h_h_near_dict[key].GetY()[0]:.2e} & {h_h_away_dict[key].GetY()[0]:.2e}  & {h_h_ue_dict[key].GetY()[0]:.2e} & {h_h_total_dict[key].GetY()[0]:.2e} \\\\")



h_lambda_near_ratios_0_20 = []
central_value = 0
for key in h_lambda_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_near_dict[key].GetY()[2]
    else:
        h_lambda_near_ratios_0_20.append((h_lambda_near_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_lambda_near_ratios_0_20)
h_lambda_near_systematic_0_20 = rms_0_20*central_value
h_lambda_near_ratios_20_50 = []
central_value = 0
for key in h_lambda_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_near_dict[key].GetY()[1]
    else:
        h_lambda_near_ratios_20_50.append((h_lambda_near_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_lambda_near_ratios_20_50)
h_lambda_near_systematic_20_50 = rms_20_50*central_value
h_lambda_near_ratios_50_80 = []
central_value = 0
for key in h_lambda_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_near_dict[key].GetY()[0]
    else:
        h_lambda_near_ratios_50_80.append((h_lambda_near_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_lambda_near_ratios_50_80)
near_yield_0_20_systematic = math.sqrt(rms_0_20**2 + total_dphi_sys_0_20**2)
near_yield_0_20_upper_systematic = math.sqrt((near_yield_0_20_systematic/2)**2 + 0.08**2)
near_yield_0_20_lower_systematic = (near_yield_0_20_systematic/2)
near_yield_20_50_systematic = math.sqrt(rms_20_50**2 + total_dphi_sys_20_50**2)
near_yield_20_50_upper_systematic = math.sqrt((near_yield_20_50_systematic/2)**2 + 0.08**2)
near_yield_20_50_lower_systematic = (near_yield_20_50_systematic/2)
near_yield_50_80_systematic = math.sqrt(rms_50_80**2 + total_dphi_sys_50_80**2)
near_yield_50_80_upper_systematic = math.sqrt((near_yield_50_80_systematic/2)**2 + 0.08**2)
near_yield_50_80_lower_systematic = (near_yield_50_80_systematic/2)


h_lambda_away_ratios_0_20 = []
central_value = 0
for key in h_lambda_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_away_dict[key].GetY()[2]
    else:
        h_lambda_away_ratios_0_20.append((h_lambda_away_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_lambda_away_ratios_0_20)
h_lambda_away_ratios_20_50 = []
central_value = 0
for key in h_lambda_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_away_dict[key].GetY()[1]
    else:
        h_lambda_away_ratios_20_50.append((h_lambda_away_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_lambda_away_ratios_20_50)
h_lambda_away_ratios_50_80 = []
central_value = 0
for key in h_lambda_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_away_dict[key].GetY()[0]
    else:
        h_lambda_away_ratios_50_80.append((h_lambda_away_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_lambda_away_ratios_50_80)


away_yield_0_20_systematic = math.sqrt(rms_0_20**2 + total_dphi_sys_0_20**2)
away_yield_20_50_systematic = math.sqrt(rms_20_50**2 + total_dphi_sys_20_50**2)
away_yield_50_80_systematic = math.sqrt(rms_50_80**2 + total_dphi_sys_50_80**2)

h_lambda_ue_ratios_0_20 = []
central_value = 0
for key in h_lambda_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_ue_dict[key].GetY()[2]
    else:
        h_lambda_ue_ratios_0_20.append((h_lambda_ue_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_lambda_ue_ratios_0_20)
h_lambda_ue_systematic_0_20 = rms_0_20*central_value
h_lambda_ue_ratios_20_50 = []
central_value = 0
for key in h_lambda_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_ue_dict[key].GetY()[1]
    else:
        h_lambda_ue_ratios_20_50.append((h_lambda_ue_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_lambda_ue_ratios_20_50)
h_lambda_ue_systematic_20_50 = rms_20_50*central_value
h_lambda_ue_ratios_50_80 = []
central_value = 0
for key in h_lambda_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_lambda_ue_dict[key].GetY()[0]
    else:
        h_lambda_ue_ratios_50_80.append((h_lambda_ue_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_lambda_ue_ratios_50_80)
h_lambda_ue_systematic_50_80 = rms_50_80*central_value

ue_yield_0_20_systematic = math.sqrt(rms_0_20**2 + total_dphi_sys_0_20**2)
ue_yield_20_50_systematic = math.sqrt(rms_20_50**2 + total_dphi_sys_20_50**2)
ue_yield_50_80_systematic = math.sqrt(rms_50_80**2 + total_dphi_sys_50_80**2)

h_lambda_total_ratios_0_20 = []
central_valtotal = 0
for key in h_lambda_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_valtotal = h_lambda_total_dict[key].GetY()[2]
    else:
        h_lambda_total_ratios_0_20.append((h_lambda_total_dict[key].GetY()[2]/central_valtotal) - 1)
rms_0_20 = get_rms(h_lambda_total_ratios_0_20)
h_lambda_total_systematic_0_20 = rms_0_20*central_valtotal
h_lambda_total_ratios_20_50 = []
central_valtotal = 0
for key in h_lambda_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_valtotal = h_lambda_total_dict[key].GetY()[1]
    else:
        h_lambda_total_ratios_20_50.append((h_lambda_total_dict[key].GetY()[1]/central_valtotal) - 1)
rms_20_50 = get_rms(h_lambda_total_ratios_20_50)
h_lambda_total_systematic_20_50 = rms_20_50*central_valtotal
h_lambda_total_ratios_50_80 = []
central_valtotal = 0
for key in h_lambda_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_valtotal = h_lambda_total_dict[key].GetY()[0]
    else:
        h_lambda_total_ratios_50_80.append((h_lambda_total_dict[key].GetY()[0]/central_valtotal) - 1)
rms_50_80 = get_rms(h_lambda_total_ratios_50_80)
h_lambda_total_systematic_50_80 = rms_50_80*central_valtotal

total_yield_0_20_systematic = math.sqrt(rms_0_20**2 + total_dphi_sys_0_20**2)
total_yield_20_50_systematic = math.sqrt(rms_20_50**2 + total_dphi_sys_20_50**2)
total_yield_50_80_systematic = math.sqrt(rms_50_80**2 + total_dphi_sys_50_80**2)

h_h_near_ratios_0_20 = []
central_value = 0
for key in h_h_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_near_dict[key].GetY()[2]
    else:
        h_h_near_ratios_0_20.append((h_h_near_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_h_near_ratios_0_20)
h_h_near_systematic_0_20 = rms_0_20*central_value
h_h_near_ratios_20_50 = []
central_value = 0
for key in h_h_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_near_dict[key].GetY()[1]
    else:
        h_h_near_ratios_20_50.append((h_h_near_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_h_near_ratios_20_50)
h_h_near_systematic_20_50 = rms_20_50*central_value
h_h_near_ratios_50_80 = []
central_value = 0
for key in h_h_near_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_near_dict[key].GetY()[0]
    else:
        h_h_near_ratios_50_80.append((h_h_near_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_h_near_ratios_50_80)
hh_near_yield_0_20_systematic = math.sqrt(rms_0_20**2 + hh_total_dphi_sys_0_20**2)
hh_near_yield_20_50_systematic = math.sqrt(rms_20_50**2 + hh_total_dphi_sys_20_50**2)
hh_near_yield_50_80_systematic = math.sqrt(rms_50_80**2 + hh_total_dphi_sys_50_80**2)

h_h_away_ratios_0_20 = []
central_value = 0
for key in h_h_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_away_dict[key].GetY()[2]
    else:
        h_h_away_ratios_0_20.append((h_h_away_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_h_away_ratios_0_20)
h_h_away_systematic_0_20 = rms_0_20*central_value
h_h_away_ratios_20_50 = []
central_value = 0
for key in h_h_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_away_dict[key].GetY()[1]
    else:
        h_h_away_ratios_20_50.append((h_h_away_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_h_away_ratios_20_50)
h_h_away_systematic_20_50 = rms_20_50*central_value
h_h_away_ratios_50_80 = []
central_value = 0
for key in h_h_away_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_away_dict[key].GetY()[0]
    else:
        h_h_away_ratios_50_80.append((h_h_away_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_h_away_ratios_50_80)
hh_away_yield_0_20_systematic = math.sqrt(rms_0_20**2 + hh_total_dphi_sys_0_20**2)
hh_away_yield_20_50_systematic = math.sqrt(rms_20_50**2 + hh_total_dphi_sys_20_50**2)
hh_away_yield_50_80_systematic = math.sqrt(rms_50_80**2 + hh_total_dphi_sys_50_80**2)

h_h_ue_ratios_0_20 = []
central_value = 0
for key in h_h_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_ue_dict[key].GetY()[2]
    else:
        h_h_ue_ratios_0_20.append((h_h_ue_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_h_ue_ratios_0_20)
h_h_ue_systematic_0_20 = rms_0_20*central_value
h_h_ue_ratios_20_50 = []
central_value = 0
for key in h_h_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_ue_dict[key].GetY()[1]
    else:
        h_h_ue_ratios_20_50.append((h_h_ue_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_h_ue_ratios_20_50)
h_h_ue_systematic_20_50 = rms_20_50*central_value
h_h_ue_ratios_50_80 = []
central_value = 0
for key in h_h_ue_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_ue_dict[key].GetY()[0]
    else:
        h_h_ue_ratios_50_80.append((h_h_ue_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_h_ue_ratios_50_80)
hh_ue_yield_0_20_systematic = math.sqrt(rms_0_20**2 + hh_total_dphi_sys_0_20**2)
hh_ue_yield_20_50_systematic = math.sqrt(rms_20_50**2 + hh_total_dphi_sys_20_50**2)
hh_ue_yield_50_80_systematic = math.sqrt(rms_50_80**2 + hh_total_dphi_sys_50_80**2)

h_h_total_ratios_0_20 = []
central_value = 0
for key in h_h_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_total_dict[key].GetY()[2]
    else:
        h_h_total_ratios_0_20.append((h_h_total_dict[key].GetY()[2]/central_value) - 1)
rms_0_20 = get_rms(h_h_total_ratios_0_20)
h_h_total_systematic_0_20 = rms_0_20*central_value
h_h_total_ratios_20_50 = []
central_value = 0
for key in h_h_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_total_dict[key].GetY()[1]
    else:
        h_h_total_ratios_20_50.append((h_h_total_dict[key].GetY()[1]/central_value) - 1)
rms_20_50 = get_rms(h_h_total_ratios_20_50)
h_h_total_systematic_20_50 = rms_20_50*central_value
h_h_total_ratios_50_80 = []
central_value = 0
for key in h_h_total_dict:
    if key == CENTRAL_TECHNIQUE:
        central_value = h_h_total_dict[key].GetY()[0]
    else:
        h_h_total_ratios_50_80.append((h_h_total_dict[key].GetY()[0]/central_value) - 1)
rms_50_80 = get_rms(h_h_total_ratios_50_80)
hh_total_yield_0_20_systematic = math.sqrt(rms_0_20**2 + hh_total_dphi_sys_0_20**2)
hh_total_yield_20_50_systematic = math.sqrt(rms_20_50**2 + hh_total_dphi_sys_20_50**2)
hh_total_yield_50_80_systematic = math.sqrt(rms_50_80**2 + hh_total_dphi_sys_50_80**2)

if PRINT_SYSTEMATICS:
    print("---------------h lambda yield systematics----------------")
    print(f"0-20\% & +{near_yield_0_20_upper_systematic:.2e}, -{near_yield_0_20_lower_systematic:.2e} & {away_yield_0_20_systematic:.2e}  & {ue_yield_0_20_systematic:.2e} & {total_yield_0_20_systematic:.2e} \\\\")
    print(f"20-50\% & +{near_yield_20_50_upper_systematic:.2e}, -{near_yield_20_50_lower_systematic:.2e}  & {away_yield_20_50_systematic:.2e}  & {ue_yield_20_50_systematic:.2e} & {total_yield_20_50_systematic:.2e} \\\\")
    print(f"50-80\% & +{near_yield_50_80_upper_systematic:.2e}, -{near_yield_50_80_lower_systematic:.2e}  & {away_yield_50_80_systematic:.2e}  & {ue_yield_50_80_systematic:.2e} & {total_yield_50_80_systematic:.2e} \\\\")

    print("---------------h h yield systematics----------------")
    print(f"0-20\% & {hh_near_yield_0_20_systematic:.2e}   & {hh_away_yield_0_20_systematic:.2e}  & {hh_ue_yield_0_20_systematic:.2e} & {hh_total_yield_0_20_systematic:.2e} \\\\")
    print(f"20-50\% & {hh_near_yield_20_50_systematic:.2e} & {hh_away_yield_20_50_systematic:.2e}  & {hh_ue_yield_20_50_systematic:.2e} & {hh_total_yield_20_50_systematic:.2e} \\\\")
    print(f"50-80\% & {hh_near_yield_50_80_systematic:.2e} & {hh_away_yield_50_80_systematic:.2e}  & {hh_ue_yield_50_80_systematic:.2e} & {hh_total_yield_50_80_systematic:.2e} \\\\")


mult_list = arr.array('d', [35, 65, 90])
mult_error_list_stat = arr.array('d', [15, 15, 10])
mult_error_list_sys = arr.array('d', [3, 3, 3])


near_graph = h_lambda_near_dict[CENTRAL_TECHNIQUE]
away_graph = h_lambda_away_dict[CENTRAL_TECHNIQUE]
ue_graph = h_lambda_ue_dict[CENTRAL_TECHNIQUE]
total_graph = h_lambda_total_dict[CENTRAL_TECHNIQUE]


hh_near_graph = h_h_near_dict[CENTRAL_TECHNIQUE]
hh_away_graph = h_h_away_dict[CENTRAL_TECHNIQUE]
hh_ue_graph = h_h_ue_dict[CENTRAL_TECHNIQUE]
hh_total_graph = h_h_total_dict[CENTRAL_TECHNIQUE]

near_yield_values = arr.array('d', [near_graph.GetY()[0], near_graph.GetY()[1], near_graph.GetY()[2]])
near_yield_errors = arr.array('d', [near_graph.GetEY()[0], near_graph.GetEY()[1], near_graph.GetEY()[2]])

away_yield_values = arr.array('d', [away_graph.GetY()[0], away_graph.GetY()[1], away_graph.GetY()[2]])
away_yield_errors = arr.array('d', [away_graph.GetEY()[0], away_graph.GetEY()[1], away_graph.GetEY()[2]])

ue_yield_values = arr.array('d', [ue_graph.GetY()[0], ue_graph.GetY()[1], ue_graph.GetY()[2]])
ue_yield_errors = arr.array('d', [ue_graph.GetEY()[0], ue_graph.GetEY()[1], ue_graph.GetEY()[2]])

total_yield_values = arr.array('d', [total_graph.GetY()[0], total_graph.GetY()[1], total_graph.GetY()[2]])
total_yield_errors = arr.array('d', [total_graph.GetEY()[0], total_graph.GetEY()[1], total_graph.GetEY()[2]])

hh_near_yield_values = arr.array('d', [hh_near_graph.GetY()[0], hh_near_graph.GetY()[1], hh_near_graph.GetY()[2]])
hh_near_yield_errors = arr.array('d', [hh_near_graph.GetEY()[0], hh_near_graph.GetEY()[1], hh_near_graph.GetEY()[2]])

hh_away_yield_values = arr.array('d', [hh_away_graph.GetY()[0], hh_away_graph.GetY()[1], hh_away_graph.GetY()[2]])
hh_away_yield_errors = arr.array('d', [hh_away_graph.GetEY()[0], hh_away_graph.GetEY()[1], hh_away_graph.GetEY()[2]])

hh_ue_yield_values = arr.array('d', [hh_ue_graph.GetY()[0], hh_ue_graph.GetY()[1], hh_ue_graph.GetY()[2]])
hh_ue_yield_errors = arr.array('d', [hh_ue_graph.GetEY()[0], hh_ue_graph.GetEY()[1], hh_ue_graph.GetEY()[2]])

hh_total_yield_values = arr.array('d', [hh_total_graph.GetY()[0], hh_total_graph.GetY()[1], hh_total_graph.GetY()[2]])
hh_total_yield_errors = arr.array('d', [hh_total_graph.GetEY()[0], hh_total_graph.GetEY()[1], hh_total_graph.GetEY()[2]])

near_yield_systematics_upper = arr.array('d', [near_yield_50_80_upper_systematic*near_graph.GetY()[0], near_yield_20_50_upper_systematic*near_graph.GetY()[1], near_yield_0_20_upper_systematic*near_graph.GetY()[2]])
near_yield_systematics_lower = arr.array('d', [near_yield_50_80_lower_systematic*near_graph.GetY()[0], near_yield_20_50_lower_systematic*near_graph.GetY()[1], near_yield_0_20_lower_systematic*near_graph.GetY()[2]])
away_yield_systematics = arr.array('d', [away_yield_50_80_systematic*away_graph.GetY()[0], away_yield_20_50_systematic*away_graph.GetY()[1], away_yield_0_20_systematic*away_graph.GetY()[2]])
ue_yield_systematics = arr.array('d', [ue_yield_50_80_systematic*ue_graph.GetY()[0], ue_yield_20_50_systematic*ue_graph.GetY()[1], ue_yield_0_20_systematic*ue_graph.GetY()[2]])
total_yield_systematics = arr.array('d', [total_yield_50_80_systematic*total_graph.GetY()[0], total_yield_20_50_systematic*total_graph.GetY()[1], total_yield_0_20_systematic*total_graph.GetY()[2]])

hh_near_yield_systematics = arr.array('d', [hh_near_yield_50_80_systematic*hh_near_graph.GetY()[0], hh_near_yield_20_50_systematic*hh_near_graph.GetY()[1], hh_near_yield_0_20_systematic*hh_near_graph.GetY()[2]])
hh_away_yield_systematics = arr.array('d', [hh_away_yield_50_80_systematic*hh_away_graph.GetY()[0], hh_away_yield_20_50_systematic*hh_away_graph.GetY()[1], hh_away_yield_0_20_systematic*hh_away_graph.GetY()[2]])
hh_ue_yield_systematics = arr.array('d', [hh_ue_yield_50_80_systematic*hh_ue_graph.GetY()[0], hh_ue_yield_20_50_systematic*hh_ue_graph.GetY()[1], hh_ue_yield_0_20_systematic*hh_ue_graph.GetY()[2]])
hh_total_yield_systematics = arr.array('d', [hh_total_yield_50_80_systematic*hh_total_graph.GetY()[0], hh_total_yield_20_50_systematic*hh_total_graph.GetY()[1], hh_total_yield_0_20_systematic*hh_total_graph.GetY()[2]])

near_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, near_yield_values, mult_error_list_sys, mult_error_list_sys, near_yield_systematics_lower, near_yield_systematics_upper)
away_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, away_yield_values, mult_error_list_sys, mult_error_list_sys, away_yield_systematics, away_yield_systematics)
ue_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, ue_yield_values, mult_error_list_sys, mult_error_list_sys, ue_yield_systematics, ue_yield_systematics)
total_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, total_yield_values, mult_error_list_sys, mult_error_list_sys, total_yield_systematics, total_yield_systematics)

hh_near_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, hh_near_yield_values, mult_error_list_sys, mult_error_list_sys, hh_near_yield_systematics, hh_near_yield_systematics)
hh_away_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, hh_away_yield_values, mult_error_list_sys, mult_error_list_sys, hh_away_yield_systematics, hh_away_yield_systematics)
hh_ue_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, hh_ue_yield_values, mult_error_list_sys, mult_error_list_sys, hh_ue_yield_systematics, hh_ue_yield_systematics)
hh_total_graph_final_syst = rt.TGraphAsymmErrors( 3, mult_list, hh_total_yield_values, mult_error_list_sys, mult_error_list_sys, hh_total_yield_systematics, hh_total_yield_systematics)

near_ratio_graph = file_dict[CENTRAL_TECHNIQUE].Get("near_ratio_graph")
away_ratio_graph = file_dict[CENTRAL_TECHNIQUE].Get("away_ratio_graph")
ue_ratio_graph = file_dict[CENTRAL_TECHNIQUE].Get("ue_ratio_graph")
total_ratio_graph = file_dict[CENTRAL_TECHNIQUE].Get("total_ratio_graph")

near_ratio_graph_values = arr.array('d', [near_ratio_graph.GetY()[i] for i in range(3)])
near_ratio_graph_sys = arr.array('d', [ratio_error(near_graph_final_syst.GetY()[i], hh_near_graph_final_syst.GetY()[i], near_graph_final_syst.GetEYhigh()[i], hh_near_graph_final_syst.GetEYhigh()[i]) for i in range(3)])
away_ratio_graph_values = arr.array('d', [away_ratio_graph.GetY()[i] for i in range(3)])
away_ratio_graph_sys = arr.array('d', [ratio_error(away_graph_final_syst.GetY()[i], hh_away_graph_final_syst.GetY()[i], away_graph_final_syst.GetEYhigh()[i], hh_away_graph_final_syst.GetEYhigh()[i]) for i in range(3)])
ue_ratio_graph_values = arr.array('d', [ue_ratio_graph.GetY()[i] for i in range(3)])
ue_ratio_graph_sys = arr.array('d', [ratio_error(ue_graph_final_syst.GetY()[i], hh_ue_graph_final_syst.GetY()[i], ue_graph_final_syst.GetEYhigh()[i], hh_ue_graph_final_syst.GetEYhigh()[i]) for i in range(3)])
total_ratio_graph_values = arr.array('d', [total_ratio_graph.GetY()[i] for i in range(3)])
total_ratio_graph_sys = arr.array('d', [ratio_error(total_graph_final_syst.GetY()[i], hh_total_graph_final_syst.GetY()[i], total_graph_final_syst.GetEYhigh()[i], hh_total_graph_final_syst.GetEYhigh()[i]) for i in range(3)])
near_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, near_ratio_graph_values, mult_error_list_sys, near_ratio_graph_sys)
away_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, away_ratio_graph_values, mult_error_list_sys, away_ratio_graph_sys)
ue_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, ue_ratio_graph_values, mult_error_list_sys, ue_ratio_graph_sys)
total_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, total_ratio_graph_values, mult_error_list_sys, total_ratio_graph_sys)

justin_infile = rt.TFile("output/fitsyst_fullpt.root")

ratiosNear = justin_infile.Get("ratiosNear")
ratiosAway = justin_infile.Get("ratiosAway")
ratiosUE = justin_infile.Get("ratiosBulk")
ratiosTotal = justin_infile.Get("ratiosTot")

ratiosNearSyst = justin_infile.Get("ratiosNearSyst")
ratiosAwaySyst = justin_infile.Get("ratiosAwaySyst")
ratiosUESyst = justin_infile.Get("ratiosBulkSyst")
ratiosTotalSyst = justin_infile.Get("ratiosTotSyst")

justin_near_ratio_graph_values = arr.array('d', [near_ratio_graph.GetY()[i]/ratiosNear.GetY()[i] for i in range(3)])
justin_away_ratio_graph_values = arr.array('d', [away_ratio_graph.GetY()[i]/ratiosAway.GetY()[i] for i in range(3)])
justin_ue_ratio_graph_values = arr.array('d', [ue_ratio_graph.GetY()[i]/ratiosUE.GetY()[i] for i in range(3)])
justin_total_ratio_graph_values = arr.array('d', [total_ratio_graph.GetY()[i]/ratiosTotal.GetY()[i] for i in range(3)])

justin_near_ratio_graph_stat_errors = arr.array('d', [ratio_error(near_ratio_graph.GetY()[i], ratiosNear.GetY()[i], near_ratio_graph.GetEY()[i], ratiosNear.GetEY()[i]) for i in range(3)])
justin_away_ratio_graph_stat_errors = arr.array('d', [ratio_error(away_ratio_graph.GetY()[i], ratiosAway.GetY()[i], away_ratio_graph.GetEY()[i], ratiosAway.GetEY()[i]) for i in range(3)])
justin_ue_ratio_graph_stat_errors = arr.array('d', [ratio_error(ue_ratio_graph.GetY()[i], ratiosUE.GetY()[i], ue_ratio_graph.GetEY()[i], ratiosUE.GetEY()[i]) for i in range(3)])
justin_total_ratio_graph_stat_errors = arr.array('d', [ratio_error(total_ratio_graph.GetY()[i], ratiosTotal.GetY()[i], total_ratio_graph.GetEY()[i], ratiosTotal.GetEY()[i]) for i in range(3)])

justin_near_ratio_graph_sys = arr.array('d', [ratio_error(near_ratio_graph.GetY()[i], ratiosNear.GetY()[i], near_ratio_graph_final_syst.GetEY()[i], ratiosNearSyst.GetEY()[i]) for i in range(3)])
justin_away_ratio_graph_sys = arr.array('d', [ratio_error(away_ratio_graph.GetY()[i], ratiosAway.GetY()[i], away_ratio_graph_final_syst.GetEY()[i], ratiosAwaySyst.GetEY()[i]) for i in range(3)])
justin_ue_ratio_graph_sys = arr.array('d', [ratio_error(ue_ratio_graph.GetY()[i], ratiosUE.GetY()[i], ue_ratio_graph_final_syst.GetEY()[i], ratiosUESyst.GetEY()[i]) for i in range(3)])
justin_total_ratio_graph_sys = arr.array('d', [ratio_error(total_ratio_graph.GetY()[i], ratiosTotal.GetY()[i], total_ratio_graph_final_syst.GetEY()[i], ratiosTotalSyst.GetEY()[i]) for i in range(3)])


justin_near_ratio_graph = rt.TGraphErrors(3, mult_list, justin_near_ratio_graph_values, mult_error_list_stat, justin_near_ratio_graph_stat_errors)
justin_away_ratio_graph = rt.TGraphErrors(3, mult_list, justin_away_ratio_graph_values, mult_error_list_stat, justin_away_ratio_graph_stat_errors)
justin_ue_ratio_graph = rt.TGraphErrors(3, mult_list, justin_ue_ratio_graph_values, mult_error_list_stat, justin_ue_ratio_graph_stat_errors)
justin_total_ratio_graph = rt.TGraphErrors(3, mult_list, justin_total_ratio_graph_values, mult_error_list_stat, justin_total_ratio_graph_stat_errors)

justin_near_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, justin_near_ratio_graph_values, mult_error_list_sys, justin_near_ratio_graph_sys)
justin_away_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, justin_away_ratio_graph_values, mult_error_list_sys, justin_away_ratio_graph_sys)
justin_ue_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, justin_ue_ratio_graph_values, mult_error_list_sys, justin_ue_ratio_graph_sys)
justin_total_ratio_graph_final_syst = rt.TGraphErrors(3, mult_list, justin_total_ratio_graph_values, mult_error_list_sys, justin_total_ratio_graph_sys)

outfile_string = "output/final_yield_ratio_syst"
if LOW_PT:
    outfile_string += "_lowpt"
elif HIGH_PT:
    outfile_string += "_highpt"

outfile_string += ".root"

out_file = rt.TFile(outfile_string, "RECREATE")
out_file.cd()

near_graph.Write("near_yield_graph")
away_graph.Write("away_yield_graph")
ue_graph.Write("ue_yield_graph")
total_graph.Write("total_yield_graph")
near_graph_final_syst.Write("near_yield_graph_final_syst")
away_graph_final_syst.Write("away_yield_graph_final_syst")
ue_graph_final_syst.Write("ue_yield_graph_final_syst")
total_graph_final_syst.Write("total_yield_graph_final_syst")

hh_near_graph.Write("hh_near_yield_graph")
hh_away_graph.Write("hh_away_yield_graph")
hh_ue_graph.Write("hh_ue_yield_graph")
hh_total_graph.Write("hh_total_yield_graph")
hh_near_graph_final_syst.Write("hh_near_yield_graph_final_syst")
hh_away_graph_final_syst.Write("hh_away_yield_graph_final_syst")
hh_ue_graph_final_syst.Write("hh_ue_yield_graph_final_syst")
hh_total_graph_final_syst.Write("hh_total_yield_graph_final_syst")



near_ratio_graph.Write("near_ratio_graph")
away_ratio_graph.Write("away_ratio_graph")
ue_ratio_graph.Write("ue_ratio_graph")
total_ratio_graph.Write("total_ratio_graph")
near_ratio_graph_final_syst.Write("near_ratio_graph_final_syst")
away_ratio_graph_final_syst.Write("away_ratio_graph_final_syst")
ue_ratio_graph_final_syst.Write("ue_ratio_graph_final_syst")
total_ratio_graph_final_syst.Write("total_ratio_graph_final_syst")



justin_near_ratio_graph.Write("lambda_phi_near_ratio_graph")
justin_away_ratio_graph.Write("lambda_phi_away_ratio_graph")
justin_ue_ratio_graph.Write("lambda_phi_ue_ratio_graph")
justin_total_ratio_graph.Write("lambda_phi_total_ratio_graph")
justin_near_ratio_graph_final_syst.Write("lambda_phi_near_ratio_graph_final_syst")
justin_away_ratio_graph_final_syst.Write("lambda_phi_away_ratio_graph_final_syst")
justin_ue_ratio_graph_final_syst.Write("lambda_phi_ue_ratio_graph_final_syst")
justin_total_ratio_graph_final_syst.Write("lambda_phi_total_ratio_graph_final_syst")

out_file.Close()
#! /usr/bin/env python
import math
import logging

import array as arr
import ROOT as rt

import strangehelper as sh

# define an epsilon to avoid bin edge effects
EPSILON = 0.0001

# define the bins that will *likely* never change (index 0 = name, index 1 = lower bound, index 2 = upper bound)
TRIGGER_PT_BINS = [["trig_4_8", 4.0, 8.0 - EPSILON]]

ASSOCIATED_PT_BINS = [["assoc_2_4", 2.0, 4.0 - EPSILON], 
                      ["assoc_15_25", 1.5, 2.5 - EPSILON],
                      ["assoc_25_4", 2.5, 4.0 - EPSILON]]

CENTRALITY_BINS = [["cent_0_20", 0.0, 20.0 - EPSILON],
                   ["cent_20_50", 20.0, 50.0 - EPSILON],
                   ["cent_50_80", 50.0, 80.0 - EPSILON],
                   ["cent_0_80", 0.0, 80.0 - EPSILON]]

DELTA_ETA_BINS = [["delta_eta_12", -1.2, 1.2 - EPSILON]]

# define bins that are varied for systematics (default is always bin 0)
SIGNAL_BINS = [["signal_default", 1.102, 1.130 - EPSILON],
               ["signal_narrow", 1.108, 1.124 - EPSILON],
               ["signal_narrower", 1.112, 1.120 - EPSILON],
               ["signal_wide", 1.100, 1.132 - EPSILON],
               ["signal_wider", 1.098, 1.134 - EPSILON]]

SIDEBAND_BINS = [["sideband_default",1.135, 1.15 - EPSILON],
                 ["sideband_narrow", 1.135, 1.145 - EPSILON],
                 ["sideband_wide", 1.135, 1.16 - EPSILON],
                 ["sideband_leftshift", 1.084, 1.096 - EPSILON],
                 ["sideband_rightshift", 1.14, 1.155 - EPSILON]]

FULL_REGION_BINS = [["full_region_default", 1.10, 1.132 - EPSILON]]
ETA_BINS = [["eta_default",-0.8, 0.8 - EPSILON]]
PID_FILE_NAMES = [["pid_default", "../online/output/v0_central.root"]]

# define the default intput fulename
DEFAULT_FILE_NAME = "../online/output/v0_central.root"

# define the output file name
OUTPUT_FILE_NAME = "offline_analysis_test.root"

def get_dphi_dists(input_list,
                   input_name,
                   centrality_bin,
                   trigger_pt_bin,
                   associated_pt_bin,
                   delta_eta_bin,
                   eta_bin,
                   signal_bin,
                   sideband_bin,
                   full_region_bin,
                   output_file = None,
                   do_dihadron = True):
    
    logging.info("Getting dphi distributions for...")
    logging.info("\tInput file: " + input_name)
    logging.info("\tCentrality: " + centrality_bin[0])
    logging.info("\tTrigger pt: " + trigger_pt_bin[0])
    logging.info("\tAssociated pt: " + associated_pt_bin[0])
    logging.info("\tDelta eta: " + delta_eta_bin[0])
    logging.info("\tEta: " + eta_bin[0])
    logging.info("\tSignal bin: " + signal_bin[0])
    logging.info("\tSideband bin: " + sideband_bin[0])
    logging.info("\tFull region bin: " + full_region_bin[0])
    logging.info("\tOutput file: " + output_file.GetName() if output_file else "None")

    # get unique string for cloning to avoid name conflicts
    unique_clone_string = input_name + "_" + \
                        centrality_bin[0] + "_" + \
                        trigger_pt_bin[0] + "_" + \
                        associated_pt_bin[0] + "_" + \
                        delta_eta_bin[0] + "_" + \
                        eta_bin[0] + "_" + \
                        signal_bin[0] + "_" + \
                        sideband_bin[0] + "_" + \
                        full_region_bin[0]


    # projecting trigger distribution within the specified centrality and trigger pt bins down to 1d (pt)
    trigger_dist = input_list.FindObject("fTriggerDistEff")
    trigger_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[1], trigger_pt_bin[2])
    trigger_dist.GetAxis(2).SetRangeUser(eta_bin[1], eta_bin[2])
    trigger_dist.GetAxis(4).SetRangeUser(centrality_bin[1], centrality_bin[2])
    trigger_pt_dist = trigger_dist.Projection(0).Clone("trigger_pt_dist_" + unique_clone_string)

    # get the number of triggers
    num_triggers = trigger_pt_dist.Integral()
    logging.info("\t\tNumber of triggers: " + str(num_triggers))

    # projecting lambda distribution within the specified centrality and associated pt bins down to 1d (mass)
    lambda_dist = input_list.FindObject("fTriggeredLambdaDist")
    lambda_dist.GetAxis(0).SetRangeUser(associated_pt_bin[1], associated_pt_bin[2])
    lambda_dist.GetAxis(2).SetRangeUser(eta_bin[1], eta_bin[2])
    lambda_dist.GetAxis(4).SetRangeUser(centrality_bin[1], centrality_bin[2])
    lambda_mass_dist = lambda_dist.Projection(3).Clone("lambda_mass_dist_" + unique_clone_string)


    # do lambda fitting to extract signal and background
    lambda_signal, lambda_background, signal_scale = sh.fit_lambda_mass(lambda_mass_dist,
                                                        signal_bin,
                                                        sideband_bin,
                                                        full_region_bin,
                                                        unique_clone_string,
                                                        output_file)

    logging.info("\t\tLambda signal: " + str(lambda_signal))
    logging.info("\t\tLambda background: " + str(lambda_background))
    logging.info("\t\tLambda signal scale: " + str(signal_scale))
    lambda_signal_over_total = lambda_signal / (lambda_signal + lambda_background)

    if do_dihadron:
        # projecting h-h distribution within the specified centrality, trigger pt, and associated pt bins down to 3d for mixed event correction
        h_h_dist = input_list.FindObject("fDphiHHEff")    
        h_h_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[1], trigger_pt_bin[2])
        h_h_dist.GetAxis(1).SetRangeUser(associated_pt_bin[1], associated_pt_bin[2])
        h_h_dist.GetAxis(5).SetRangeUser(centrality_bin[1], centrality_bin[2])
        h_h_dist = h_h_dist.Projection(2, 3, 4).Clone("h_h_dist_" + unique_clone_string)

        h_h_mixed_dist = input_list.FindObject("fDphiHHMixed")
        h_h_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[1], trigger_pt_bin[2])
        h_h_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin[1], associated_pt_bin[2])
        h_h_mixed_dist.GetAxis(5).SetRangeUser(centrality_bin[1], centrality_bin[2])
        h_h_mixed_dist = h_h_mixed_dist.Projection(2, 3, 4).Clone("h_h_mixed_dist_" + unique_clone_string)

        # do dihadron mixed event correction
        h_h_2d_mixcor = sh.make_h_h_mixed_corrections(h_h_dist, h_h_mixed_dist, unique_clone_string, output_file)

        # do per-trigger scaling
        h_h_2d_mixcor.Scale(1.0/num_triggers)

        # project down to dPhi in given dEta range
        h_h_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bin[1], delta_eta_bin[2])
        h_h_dphi = h_h_2d_mixcor.ProjectionY("h_h_dphi_" + unique_clone_string)
    else:
        # function returns h_h and h_lambda distributions
        h_h_dphi = None
    

    # project h-lambda distribution within the specified centrality, trigger pt, and associated pt bins down to 4d for mixed event correction
    h_lambda_dist_axes = arr.array('i', [2, 3, 4, 5])

    h_lambda_dist = input_list.FindObject("fDphiHLambdaEff")
    h_lambda_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[1], trigger_pt_bin[2])
    h_lambda_dist.GetAxis(1).SetRangeUser(associated_pt_bin[1], associated_pt_bin[2])
    h_lambda_dist.GetAxis(6).SetRangeUser(centrality_bin[1], centrality_bin[2])
    h_lambda_dist = h_lambda_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_dist_" + unique_clone_string)

    h_lambda_mixed_dist = input_list.FindObject("fDphiHLambdaMixed")
    h_lambda_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[1], trigger_pt_bin[2])
    h_lambda_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin[1], associated_pt_bin[2])
    h_lambda_mixed_dist.GetAxis(6).SetRangeUser(centrality_bin[1], centrality_bin[2])
    h_lambda_mixed_dist = h_lambda_mixed_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_mixed_dist_" + unique_clone_string) 

    # do h-lambda mixed event correction
    h_lambda_2d_mixcor_signal, h_lambda_2d_mixcor_sideband = sh.make_h_lambda_mixed_corrections(h_lambda_dist, h_lambda_mixed_dist,
                                                            signal_bin, sideband_bin, unique_clone_string, output_file)

    # do per-trigger scaling 
    h_lambda_2d_mixcor_signal.Scale(1.0/num_triggers)
    h_lambda_2d_mixcor_sideband.Scale(1.0/num_triggers)

    # set deta range before subtracting off h-lambda background
    h_lambda_2d_mixcor_signal.GetXaxis().SetRangeUser(delta_eta_bin[1], delta_eta_bin[2])
    h_lambda_2d_mixcor_sideband.GetXaxis().SetRangeUser(delta_eta_bin[1], delta_eta_bin[2])
    
    # subtract off h-lambda background using sideband
    h_lambda_2d_subtracted = sh.apply_sideband_subtraction(h_lambda_2d_mixcor_signal, h_lambda_2d_mixcor_sideband,
                                                            lambda_signal_over_total, unique_clone_string, output_file)
    
    # correct for tails in lambda mass distribution
    h_lambda_2d_subtracted.Scale(signal_scale)

    # correct for PID scaling
    two_sigma_normal = 0.9544
    three_sigma_normal = 0.9974
    pid_correction = (1.0 / two_sigma_normal) * (1.0 / three_sigma_normal)
    h_lambda_2d_subtracted.Scale(pid_correction)

    # correct for branching ratio
    lambda_br = 0.639
    h_lambda_2d_subtracted.Scale(1/lambda_br)

    # correct for two-track efficiency
    h_lambda_2d_subtracted = sh.apply_two_track_correction(h_lambda_2d_subtracted, associated_pt_bin[0])

    h_lambda_dphi = h_lambda_2d_subtracted.ProjectionY("h_lambda_dphi_" + unique_clone_string)

    return h_lambda_dphi, h_h_dphi


if __name__ == "__main__":

    # initialize the logging
    logging.basicConfig(filename='offline_analysis.log',
                        filemode='w',
                        format='%(message)s', 
                        level=logging.DEBUG)

    logging.info("Starting the offline analysis...")

    # initialize the output file
    logging.info("Initializing the output file under the name: " + OUTPUT_FILE_NAME)
    output_file = rt.TFile(OUTPUT_FILE_NAME, "RECREATE")

    # initialize the central input file
    central_file = rt.TFile(DEFAULT_FILE_NAME, "READ")
    central_list = central_file.Get("h-lambda")

    # run through the centrality, trigger pt, and associated pt bins for our central values
    for trigger_pt_bin in TRIGGER_PT_BINS:
        for associated_pt_bin in ASSOCIATED_PT_BINS:
            for centrality_bin in CENTRALITY_BINS:
                h_lambda_dphi_default, h_h_dphi_default = get_dphi_dists(central_list,
                                                        "default",
                                                         centrality_bin,
                                                         trigger_pt_bin,
                                                         associated_pt_bin, 
                                                         DELTA_ETA_BINS[0],
                                                         ETA_BINS[0],
                                                         SIGNAL_BINS[0],
                                                         SIDEBAND_BINS[0],
                                                         FULL_REGION_BINS[0],
                                                         output_file)
                
                h_lambda_signal_variations = {signal_bin[0] : get_dphi_dists(central_list, 
                                                                     "default",
                                                                     centrality_bin,
                                                                     trigger_pt_bin,
                                                                     associated_pt_bin,
                                                                     DELTA_ETA_BINS[0],
                                                                     ETA_BINS[0],
                                                                     signal_bin,
                                                                     SIDEBAND_BINS[0],
                                                                     FULL_REGION_BINS[0],
                                                                     output_file=None,
                                                                     do_dihadron=False)[0] for signal_bin in SIGNAL_BINS[1:]}

                h_lambda_sideband_variations = {sideband_bin[0] : get_dphi_dists(central_list,
                                                                        "default",
                                                                        centrality_bin,
                                                                        trigger_pt_bin,
                                                                        associated_pt_bin,
                                                                        DELTA_ETA_BINS[0],
                                                                        ETA_BINS[0],
                                                                        SIGNAL_BINS[0],
                                                                        sideband_bin,
                                                                        FULL_REGION_BINS[0],
                                                                        output_file=None,
                                                                        do_dihadron=False)[0] for sideband_bin in SIDEBAND_BINS[1:]} 

                h_lambda_dphi_signal_sys = sh.get_systematic_uncertainty_dphi(h_lambda_signal_variations, h_lambda_dphi_default, output_file)
                print(f"Lambda signal systematic uncertainty: {h_lambda_dphi_signal_sys} for assoc pt bin: {associated_pt_bin} and trigger pt bin: {trigger_pt_bin}")
                quit()
                h_lambda_dphi_sideband_sys = sh.get_systematic_uncertainty_dphi(h_lambda_sideband_variations, h_lambda_dphi_default, output_file)
                print(f"Lambda sideband systematic uncertainty: {h_lambda_dphi_sideband_sys} for assoc pt bin: {associated_pt_bin} and trigger pt bin: {trigger_pt_bin}")
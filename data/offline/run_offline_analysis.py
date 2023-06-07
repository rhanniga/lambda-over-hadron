#!/usr/bin/env python
import logging

from sys import argv

import array as arr
import ROOT as rt

import strangehelper as sh
from systematics_helper import AnalysisBin, AnalysisBins, DphiSystematicHelper, WidthSystematicHelper
from yield_extractor import YieldExtractor, FitType

# define an epsilon to avoid bin edge effects
EPSILON = 0.0001

# define the bins that will *likely* never change
TRIGGER_PT_BINS = AnalysisBins("trigger_pt_bins", 
                               [AnalysisBin("trigger_4_8", 4.0, 8.0 - EPSILON)])
ASSOCIATED_PT_BINS = AnalysisBins("associated_pt_bins",
                                [AnalysisBin("assoc_2_4", 2.0, 4.0 - EPSILON),
                                AnalysisBin("assoc_15_25", 1.5, 2.5 - EPSILON),
                                AnalysisBin("assoc_25_4", 2.5, 4.0 - EPSILON)])
CENTRALITY_BINS = AnalysisBins("centrality_bins",
                                [AnalysisBin("cent_0_20", 0.0, 20.0 - EPSILON),
                                AnalysisBin("cent_20_50", 20.0, 50.0 - EPSILON),
                                AnalysisBin("cent_50_80", 50.0, 80.0 - EPSILON),
                                AnalysisBin("cent_0_80", 0.0, 80.0 - EPSILON)])
DELTA_ETA_BINS = AnalysisBins("delta_eta_bins",
                                [AnalysisBin("delta_eta_12", -1.2, 1.2 - EPSILON)])

# define bins that are varied for systematics (default is always bin 0)
SIGNAL_BINS = AnalysisBins("signal_bins",
                            [AnalysisBin("signal_default", 1.102, 1.130 - EPSILON),
                            AnalysisBin("signal_narrow", 1.108, 1.124 - EPSILON),
                            AnalysisBin("signal_narrower", 1.112, 1.120 - EPSILON)])
                            # AnalysisBin("signal_wide", 1.100, 1.132 - EPSILON),
                            # AnalysisBin("signal_wider", 1.098, 1.134 - EPSILON)])

SIDEBAND_BINS = AnalysisBins("sideband_bins",
                            [AnalysisBin("sideband_default", 1.135, 1.15 - EPSILON),
                            # AnalysisBin("sideband_narrow", 1.135, 1.145 - EPSILON),
                            # AnalysisBin("sideband_wide", 1.135, 1.16 - EPSILON)])
                            AnalysisBin("sideband_leftshift", 1.084, 1.096 - EPSILON),
                            AnalysisBin("sideband_rightshift", 1.14, 1.155 - EPSILON)])

ETA_BINS = AnalysisBins("eta_bins",
                        [AnalysisBin("eta_default", -0.8, 0.8 - EPSILON)])

FULL_REGION_BINS = AnalysisBins("full_region_bins",
                                [AnalysisBin("full_region_default", 1.10, 1.132 - EPSILON)])


PID_BINS = AnalysisBins("pid_bins",
                        [AnalysisBin("pid_default", 1, 1, "../online/output/v0_central.root"),
                        AnalysisBin("pid_narrow", 0.6, 0.6, "../online/output/v0_pid_narrow.root"),
                        AnalysisBin("pid_wide", 1.4, 1.4, "../online/output/v0_pid_wide.root")])

ALL_VARIATIONS = [SIGNAL_BINS, SIDEBAND_BINS, ETA_BINS, FULL_REGION_BINS, PID_BINS]

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
    logging.info("\tInput name: " + input_name)
    logging.info("\tCentrality: " + centrality_bin.name)
    logging.info("\tTrigger pt: " + trigger_pt_bin.name)
    logging.info("\tAssociated pt: " + associated_pt_bin.name)

    logging.info("\tDelta eta: " + delta_eta_bin.name)
    logging.info("\tEta: " + eta_bin.name)
    logging.info("\tSignal bin: " + signal_bin.name)
    logging.info("\tSideband bin: " + sideband_bin.name)
    logging.info("\tFull region bin: " + full_region_bin.name)
    logging.info("\tOutput file: " + output_file.GetName() if output_file else "\tNone")

    # get unique string for cloning to avoid name conflicts
    unique_clone_string = input_name + "_" + \
                        centrality_bin.name + "_" + \
                        trigger_pt_bin.name + "_" + \
                        associated_pt_bin.name + "_" + \
                        delta_eta_bin.name + "_" + \
                        eta_bin.name + "_" + \
                        signal_bin.name + "_" + \
                        sideband_bin.name + "_" + \
                        full_region_bin.name

    # projecting trigger distribution within the specified centrality and trigger pt bins down to 1d (pt)
    trigger_dist = input_list.FindObject("fTriggerDistEff")
    trigger_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    trigger_dist.GetAxis(2).SetRangeUser(eta_bin.lower_bound, eta_bin.upper_bound)
    trigger_dist.GetAxis(4).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    trigger_pt_dist = trigger_dist.Projection(0).Clone("trigger_pt_dist_" + unique_clone_string)

    # get the number of triggers
    num_triggers = trigger_pt_dist.Integral()
    logging.info("\t\tNumber of triggers: " + str(num_triggers))

    # projecting lambda distribution within the specified centrality and associated pt bins down to 1d (mass)
    lambda_dist = input_list.FindObject("fTriggeredLambdaDist")
    lambda_dist.GetAxis(0).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    lambda_dist.GetAxis(2).SetRangeUser(eta_bin.lower_bound, eta_bin.upper_bound)
    lambda_dist.GetAxis(4).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
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
        h_h_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
        h_h_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
        h_h_dist.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
        h_h_dist = h_h_dist.Projection(2, 3, 4).Clone("h_h_dist_" + unique_clone_string)

        h_h_mixed_dist = input_list.FindObject("fDphiHHMixed")
        h_h_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
        h_h_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
        h_h_mixed_dist.GetAxis(5).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
        h_h_mixed_dist = h_h_mixed_dist.Projection(2, 3, 4).Clone("h_h_mixed_dist_" + unique_clone_string)

        # do dihadron mixed event correction
        h_h_2d_mixcor = sh.make_h_h_mixed_corrections(h_h_dist, h_h_mixed_dist, unique_clone_string, output_file)

        # do per-trigger scaling
        h_h_2d_mixcor.Scale(1.0/num_triggers)

        # project down to dPhi in given dEta range
        h_h_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
        h_h_dphi = h_h_2d_mixcor.ProjectionY("h_h_dphi_" + unique_clone_string)
    else:
        # function returns h_h and h_lambda distributions
        h_h_dphi = None
    

    # project h-lambda distribution within the specified centrality, trigger pt, and associated pt bins down to 4d for mixed event correction
    h_lambda_dist_axes = arr.array('i', [2, 3, 4, 5])

    h_lambda_dist = input_list.FindObject("fDphiHLambdaEff")
    h_lambda_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_lambda_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_lambda_dist.GetAxis(6).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_lambda_dist = h_lambda_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_dist_" + unique_clone_string)

    h_lambda_mixed_dist = input_list.FindObject("fDphiHLambdaMixed")
    h_lambda_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin.lower_bound, trigger_pt_bin.upper_bound)
    h_lambda_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin.lower_bound, associated_pt_bin.upper_bound)
    h_lambda_mixed_dist.GetAxis(6).SetRangeUser(centrality_bin.lower_bound, centrality_bin.upper_bound)
    h_lambda_mixed_dist = h_lambda_mixed_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_mixed_dist_" + unique_clone_string) 

    # do h-lambda mixed event correction
    h_lambda_2d_mixcor_signal, h_lambda_2d_mixcor_sideband = sh.make_h_lambda_mixed_corrections(h_lambda_dist, h_lambda_mixed_dist,
                                                            signal_bin, sideband_bin, unique_clone_string, output_file)

    # do per-trigger scaling 
    h_lambda_2d_mixcor_signal.Scale(1.0/num_triggers)
    h_lambda_2d_mixcor_sideband.Scale(1.0/num_triggers)

    # set deta range before subtracting off h-lambda background
    h_lambda_2d_mixcor_signal.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
    h_lambda_2d_mixcor_sideband.GetXaxis().SetRangeUser(delta_eta_bin.lower_bound, delta_eta_bin.upper_bound)
    
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
    h_lambda_2d_subtracted = sh.apply_two_track_correction(h_lambda_2d_subtracted, associated_pt_bin.name)

    h_lambda_dphi = h_lambda_2d_subtracted.ProjectionY("h_lambda_dphi_" + unique_clone_string)

    return h_lambda_dphi, h_h_dphi

if __name__ == "__main__":

    # define the output file name
    OUTPUT_FILE_NAME = argv[1] + ".root" if len(argv) > 1 else "offline_analysis.root"
    LOG_FILE_NAME = argv[1] + ".log" if len(argv) > 1 else "offline_analysis.log"

    # initialize the logging
    logging.basicConfig(filename=LOG_FILE_NAME,
                        filemode='w',
                        format='%(message)s', 
                        level=logging.DEBUG)

    logging.info("Starting the offline analysis...")

    # initialize the output file
    logging.info("Initializing the output file under the name: " + OUTPUT_FILE_NAME)
    output_file = rt.TFile(OUTPUT_FILE_NAME, "RECREATE")

    # get the default list for the central values
    default_list = PID_BINS.analysis_bins[0].object_list

    # run through the centrality, trigger pt, and associated pt bins for our central values
    for trigger_pt_bin in TRIGGER_PT_BINS.analysis_bins:
        for associated_pt_bin in ASSOCIATED_PT_BINS.analysis_bins:
            for centrality_bin in CENTRALITY_BINS.analysis_bins:

                # first run the default to get things that can be plotted
                h_lambda_dphi_default, h_h_dphi_default = get_dphi_dists(default_list,
                                                    "default",
                                                    centrality_bin,
                                                    trigger_pt_bin,
                                                    associated_pt_bin, 
                                                    DELTA_ETA_BINS.analysis_bins[0],
                                                    ETA_BINS.analysis_bins[0],
                                                    SIGNAL_BINS.analysis_bins[0],
                                                    SIDEBAND_BINS.analysis_bins[0],
                                                    FULL_REGION_BINS.analysis_bins[0],
                                                    output_file)
            
                default_dphi_list = rt.TList()
                default_dphi_list.Add(h_lambda_dphi_default)
                default_dphi_list.Add(h_h_dphi_default)

                h_lambda_dphi_yield_extractor = YieldExtractor(h_lambda_dphi_default, centrality_bin.name, trigger_pt_bin.name, associated_pt_bin.name)
                h_lambda_dphi_yield_extractor.extract_yield(fit_type=FitType.AVG_SIX)
                h_lambda_dphi_yield_extractor.save_fits("h_lambda_dphi_yield_extractor_" + centrality_bin.name + "_" + trigger_pt_bin.name + "_" + associated_pt_bin.name, output_file)
                print(h_lambda_dphi_yield_extractor.yields[FitType.AVG_SIX])
                # handle yield systematics with yield extractor (can be done for each variation as well)
                # handle dphi systematics using dphi dists 
                # handle width systematics with ??? (TODO)
                


                output_file.WriteObject(default_dphi_list, "default_dphi_list_" + centrality_bin.name + "_" + trigger_pt_bin.name + "_" + associated_pt_bin.name)

                h_lambda_variations = {}

                for index, VARIATIONS in enumerate(ALL_VARIATIONS):
                    h_lambda_variations[VARIATIONS.name] = []
                    for variation in VARIATIONS.analysis_bins[1:]:
                        if VARIATIONS.name == "signal_bins":
                            h_lambda_variations[VARIATIONS.name].append(get_dphi_dists(default_list,
                                                             VARIATIONS.name,
                                                             centrality_bin,
                                                             trigger_pt_bin,
                                                             associated_pt_bin,
                                                             DELTA_ETA_BINS.analysis_bins[0],
                                                             ETA_BINS.analysis_bins[0],
                                                             variation,
                                                             SIDEBAND_BINS.analysis_bins[0],
                                                             FULL_REGION_BINS.analysis_bins[0],
                                                             None,
                                                             False)[0])
                        elif VARIATIONS.name == "sideband_bins":
                            h_lambda_variations[VARIATIONS.name].append(get_dphi_dists(default_list,
                                                             VARIATIONS.name,
                                                             centrality_bin,
                                                             trigger_pt_bin,
                                                             associated_pt_bin,
                                                             DELTA_ETA_BINS.analysis_bins[0],
                                                             ETA_BINS.analysis_bins[0],
                                                             SIGNAL_BINS.analysis_bins[0],
                                                             variation,
                                                             FULL_REGION_BINS.analysis_bins[0],
                                                             None,
                                                             False)[0])
                        elif VARIATIONS.name == "full_region_bins":
                            h_lambda_variations[VARIATIONS.name].append(get_dphi_dists(default_list,
                                                             VARIATIONS.name,
                                                             centrality_bin,
                                                             trigger_pt_bin,
                                                             associated_pt_bin,
                                                             DELTA_ETA_BINS.analysis_bins[0],
                                                             ETA_BINS.analysis_bins[0],
                                                             SIGNAL_BINS.analysis_bins[0],
                                                             SIDEBAND_BINS.analysis_bins[0],
                                                             variation,
                                                             None,
                                                             False)[0])
                        elif VARIATIONS.name == "pid_bins":
                            tmp_dist = get_dphi_dists(variation.object_list,
                                                             variation.name,
                                                             centrality_bin,
                                                             trigger_pt_bin,
                                                             associated_pt_bin,
                                                             DELTA_ETA_BINS.analysis_bins[0],
                                                             ETA_BINS.analysis_bins[0],
                                                             SIGNAL_BINS.analysis_bins[0],
                                                             SIDEBAND_BINS.analysis_bins[0],
                                                             FULL_REGION_BINS.analysis_bins[0],
                                                             None,
                                                             False)[0]


                            h_lambda_variations[VARIATIONS.name].append(tmp_dist)

                h_lambda_dphi_systematics = DphiSystematicHelper(h_lambda_dphi_default, h_lambda_variations)                     
                h_lambda_dphi_systematics.calculate_systematics()
                print(h_lambda_dphi_systematics.contributions)
                h_lambda_width_systematics = WidthSystematicHelper(h_lambda_dphi_default, h_lambda_variations, 
                                                                   centrality_bin.name, trigger_pt_bin.name, associated_pt_bin.name)
                h_lambda_width_systematics.extract_all_widths()
                h_lambda_width_systematics.calculate_systematics()
                print(h_lambda_width_systematics.contributions)
                quit()
    
        
    # close the output file
    output_file.Close()

    # close all files corresponding to variations that had separate files
    for variation in PID_BINS.analysis_bins:
        variation.input_file.Close()
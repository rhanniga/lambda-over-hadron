#! /usr/bin/env python
import math
import logging

import array as arr
import ROOT as rt

import strangehelper as sh

# initialize the logging
logging.basicConfig(filename='offline_analysis.log',
                    filemode='w',
                    format='%(message)s', 
                    level=logging.DEBUG)

logging.info("Starting the offline analysis...")

# initialize the output file
logging.info("Initializing the output file under the name offline_analysis.root...")
output_file = rt.TFile("offline_analysis.root", "RECREATE")

# epsilon for upper bound of bins (to avoid double counting, root defaults to upper bound exclusive)
epsilon = 0.00001

# branching ratio for lambda -> p + pi
lambda_br = 0.639

# default PID scaling
two_sigma_normal = 0.9544
three_sigma_normal = 0.9974
pid_correction = 1/(two_sigma_normal * three_sigma_normal)

# define the bins that will *likely* never change
trigger_pt_bins = {"trig_4_8" : [4.0, 8.0 - epsilon]}

associated_pt_bins = { "assoc_2_4" : [2.0, 4.0 - epsilon], 
                    "assoc_15_25" : [1.5, 2.5 - epsilon],
                    "assoc_25_4" : [2.5, 4.0 - epsilon]}

centrality_bins = {"cent_0_20" : [0.0, 20.0 - epsilon],
                "cent_20_50" : [20.0, 50.0 - epsilon],
                "cent_50_80" : [50.0, 80.0 - epsilon],
                "cent_0_80" : [0.0, 80.0 - epsilon]}

delta_eta_bins = {"delta_eta_12" : [-1.2, 1.2 - epsilon]}


# define the bins that will change for systematics
signal_bins = {"signal_central_value" : [1.102, 1.130 - epsilon]}
sideband_bins = {"sideband_central_value" : [1.135, 1.15 - epsilon]}
full_region_bins = {"full_region_central_value" : [1.10, 1.132 - epsilon]}
eta_bins = {"eta_central_value" : [-0.8, 0.8 - epsilon]}

central_file_name = "../online/output/v0_central.root"
central_file = rt.TFile(central_file_name)

# use mult selection list to get number of events (for now)
mult_selection_list = central_file.Get("MultSelection")
mult_selection_list = mult_selection_list.Get("cListMultSelection")
event_counter_dist = mult_selection_list.FindObject("fHistEventCounter")
num_events = event_counter_dist.GetBinContent(1)

logging.info("Number of events in central file: " + str(num_events))


central_list = central_file.Get("h-lambda")

for centrality_name, centrality_bin in centrality_bins.items():
    logging.info("Processing centrality bin: " + centrality_name)
    for trigger_pt_name, trigger_pt_bin in trigger_pt_bins.items():
        logging.info("\tProcessing trigger bin: " + trigger_pt_name)

        single_dist_axes = arr.array('i', [0, 1, 2, 3])

        # projecting trigger distribution within the specified centrality and trigger pt bins down to 1d (pt)
        trigger_dist = central_list.FindObject("fTriggerDistEff")
        trigger_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[0], trigger_pt_bin[1])
        trigger_dist.GetAxis(2).SetRangeUser(eta_bins["eta_central_value"][0], eta_bins["eta_central_value"][1])
        trigger_dist.GetAxis(4).SetRangeUser(centrality_bin[0], centrality_bin[1])
        trigger_pt_dist = trigger_dist.Projection(0).Clone("trigger_pt_dist_" + centrality_name + "_" +  trigger_pt_name)

        # get the number of triggers
        num_triggers = trigger_pt_dist.Integral()
        logging.info("\tNumber of triggers: " + str(num_triggers))
        logging.info("\tTriggers per min-bias event: " + str(num_triggers / num_events))

        for associated_pt_name, associated_pt_bin in associated_pt_bins.items():

            logging.info("\t\tProcessing associated bin: " + associated_pt_name)

            # projecting lambda distribution within the specified centrality and associated pt bins down to 1d (mass)
            lambda_dist = central_list.FindObject("fTriggeredLambdaDist")
            lambda_dist.GetAxis(0).SetRangeUser(associated_pt_bin[0], associated_pt_bin[1])
            lambda_dist.GetAxis(2).SetRangeUser(eta_bins["eta_central_value"][0], eta_bins["eta_central_value"][1])
            lambda_dist.GetAxis(4).SetRangeUser(centrality_bin[0], centrality_bin[1])
            lambda_mass_dist = lambda_dist.Projection(3).Clone("lambda_mass_dist_" + centrality_name + "_" +  associated_pt_name)

            # do lambda fitting to extract signal and background
            lambda_signal, lambda_background, signal_scale = sh.fit_lambda_mass(lambda_mass_dist,
                                                                centrality_name, 
                                                                associated_pt_name,
                                                                signal_bins["signal_central_value"],
                                                                full_region_bins["full_region_central_value"], 
                                                                sideband_bins["sideband_central_value"],
                                                                output_file)
            logging.info("\t\t\tLambda signal: " + str(lambda_signal))
            logging.info("\t\t\tLambda background: " + str(lambda_background))
            logging.info("\t\t\tLambda signal scale: " + str(signal_scale))
            logging.info("\t\t\tLambdas per min-bias event: " + str(lambda_signal / num_events))
            lambda_signal_over_total = lambda_signal / (lambda_signal + lambda_background)

            # projecting h-h distribution within the specified centrality, trigger pt, and associated pt bins down to 3d for mixed event correction
            h_h_dist = central_list.FindObject("fDphiHHEff")    
            h_h_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[0], trigger_pt_bin[1])
            h_h_dist.GetAxis(1).SetRangeUser(associated_pt_bin[0], associated_pt_bin[1])
            h_h_dist.GetAxis(5).SetRangeUser(centrality_bin[0], centrality_bin[1])
            h_h_dist = h_h_dist.Projection(2, 3, 4).Clone("h_h_dist_" + centrality_name + "_" +  trigger_pt_name + "_" + associated_pt_name)

            h_h_mixed_dist = central_list.FindObject("fDphiHHMixed")
            h_h_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[0], trigger_pt_bin[1])
            h_h_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin[0], associated_pt_bin[1])
            h_h_mixed_dist.GetAxis(5).SetRangeUser(centrality_bin[0], centrality_bin[1])
            h_h_mixed_dist = h_h_mixed_dist.Projection(2, 3, 4).Clone("h_h_mixed_dist_" + centrality_name + "_" +  trigger_pt_name + "_" + associated_pt_name)

            # do dihadron mixed event correction
            h_h_2d_mixcor = sh.make_h_h_mixed_corrections(h_h_dist, h_h_mixed_dist, 
                                                          centrality_name, trigger_pt_name, associated_pt_name, 
                                                          output_file)
            
            # project h-lambda distribution within the specified centrality, trigger pt, and associated pt bins down to 4d for mixed event correction
            h_lambda_dist_axes = arr.array('i', [2, 3, 4, 5])

            h_lambda_dist = central_list.FindObject("fDphiHLambdaEff")
            h_lambda_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[0], trigger_pt_bin[1])
            h_lambda_dist.GetAxis(1).SetRangeUser(associated_pt_bin[0], associated_pt_bin[1])
            h_lambda_dist.GetAxis(6).SetRangeUser(centrality_bin[0], centrality_bin[1])
            h_lambda_dist = h_lambda_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_dist_" + 
                                                                                  centrality_name + "_" +  
                                                                                  trigger_pt_name + "_" + 
                                                                                  associated_pt_name)

            h_lambda_mixed_dist = central_list.FindObject("fDphiHLambdaMixed")
            h_lambda_mixed_dist.GetAxis(0).SetRangeUser(trigger_pt_bin[0], trigger_pt_bin[1])
            h_lambda_mixed_dist.GetAxis(1).SetRangeUser(associated_pt_bin[0], associated_pt_bin[1])
            h_lambda_mixed_dist.GetAxis(6).SetRangeUser(centrality_bin[0], centrality_bin[1])
            h_lambda_mixed_dist = h_lambda_mixed_dist.Projection(4, h_lambda_dist_axes).Clone("h_lambda_mixed_dist_" + 
                                                                                              centrality_name + "_" +  
                                                                                              trigger_pt_name + "_" + 
                                                                                              associated_pt_name)

            # do h-lambda mixed event correction
            h_lambda_2d_mixcor_signal, h_lambda_2d_mixcor_sideband = sh.make_h_lambda_mixed_corrections(h_lambda_dist, h_lambda_mixed_dist,
                                                                    centrality_name, trigger_pt_name, associated_pt_name,
                                                                    "signal_central_value",  "sideband_central_value",
                                                                    signal_bins["signal_central_value"], sideband_bins["sideband_central_value"],
                                                                    output_file)

            # per-trigger scaling for h-lambda and h-h distributions
            h_h_2d_mixcor.Scale(1.0/num_triggers)
            h_lambda_2d_mixcor_signal.Scale(1.0/num_triggers)
            h_lambda_2d_mixcor_sideband.Scale(1.0/num_triggers)

            # project h-h and h-lambda distributions down to 1d for fitting and extracting yields
            h_h_2d_mixcor.GetXaxis().SetRangeUser(delta_eta_bins["delta_eta_12"][0], delta_eta_bins["delta_eta_12"][1])
            h_lambda_2d_mixcor_signal.GetXaxis().SetRangeUser(delta_eta_bins["delta_eta_12"][0], delta_eta_bins["delta_eta_12"][1])
            h_lambda_2d_mixcor_sideband.GetXaxis().SetRangeUser(delta_eta_bins["delta_eta_12"][0], delta_eta_bins["delta_eta_12"][1])
            
            # subtract off h-lambda background using sideband
            h_lambda_2d_subtracted = sh.apply_sideband_subtraction(h_lambda_2d_mixcor_signal, h_lambda_2d_mixcor_sideband,
                                                                   lambda_signal_over_total,
                                                                   "h_lambda_2d_subtracted_" + centrality_name + "_" + 
                                                                   trigger_pt_name + "_" + associated_pt_name)
            
            # correct for tails in lambda mass distribution
            h_lambda_2d_subtracted.Scale(signal_scale)

            # correct for PID scaling
            h_lambda_2d_subtracted.Scale(pid_correction)

            # correct for branching ratio
            h_lambda_2d_subtracted.Scale(1/lambda_br)

            # correct for two-track efficiency
            h_lambda_2d_subtracted = sh.apply_two_track_correction(h_lambda_2d_subtracted, associated_pt_name)



            h_h_dphi = h_h_2d_mixcor.ProjectionY("h_h_dphi_" + centrality_name + "_" + 
                                                 trigger_pt_name + "_" + associated_pt_name + 
                                                 "_delta_eta_12")

            h_lambda_dphi = h_lambda_2d_subtracted.ProjectionY("h_lambda_dphi_" + centrality_name + "_" + 
                                                               trigger_pt_name + "_" + associated_pt_name +
                                                               "_delta_eta_12")

output_file.Close()
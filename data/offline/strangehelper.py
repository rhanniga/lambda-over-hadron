import logging
import math

import ROOT as rt
import fit_helper as fh


# a function to calculate parabola that passes through all three input points
def get_parabola(point_one, point_two, point_three):
    x1, x2, x3 = point_one[0], point_two[0], point_three[0]
    y1, y2, y3 = point_one[1], point_two[1], point_three[1]
    
    denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
    
    A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    B = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
    C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
    
    return C, B, A


# a function to calculate the straight line that passes through all two input points
def get_straight_line(point_one, point_two):
    x1, x2 = point_one[0], point_two[0]
    y1, y2 = point_one[1], point_two[1]
    
    m = (y2 - y1) / (x2 - x1)
    c = y1 - m * x1
    
    return c, m

# a function that fits the input lambda mass distribution with a voigt + pol1 background
def fit_lambda_mass(lambda_mass_dist,
                    signal_bin,
                    sideband_bin,
                    full_region_bin,
                    unique_clone_string, 
                    outfile):

    logging.info("\t\tFitting lambda mass dist with voigt + pol1...")

    bg_fit_left_bin = lambda_mass_dist.FindBin(1.102)
    bg_fit_right_bin = lambda_mass_dist.FindBin(1.132)
    bg_fit_left_point = [1.102, lambda_mass_dist.GetBinContent(bg_fit_left_bin)]
    bg_fit_right_point = [1.132, lambda_mass_dist.GetBinContent(bg_fit_right_bin)]
    bg_starting_params = get_straight_line(bg_fit_left_point, bg_fit_right_point)

    lambda_mass_fit = rt.TF1("lambda_mass_fit_" + unique_clone_string,
                             "[0]*TMath::Voigt(x - [1], [2], [3], 4) + pol1(4)", 1.102, 1.127)

    lambda_mass_fit.SetNpx(1000)
    lambda_mass_fit.SetParameter(0, 1.74e03)
    lambda_mass_fit.SetParameter(1, 1.116)
    lambda_mass_fit.SetParameter(2, 1.30576e-03 )
    lambda_mass_fit.SetParameter(3, 1.804166e-03)
    lambda_mass_fit.SetParameter(4, bg_starting_params[0])
    lambda_mass_fit.SetParameter(5, bg_starting_params[1])

    lambda_mass_dist_fit = lambda_mass_dist.Clone("lambda_mass_dist_fit_" + unique_clone_string)
    lambda_mass_fit_result = lambda_mass_dist_fit.Fit(lambda_mass_fit, "RS")
    lambda_mass_fit_status = lambda_mass_fit_result.Status()
    if lambda_mass_fit_status != 0:
        logging.error("\t\tERROR: Lambda mass fit failed! Status: " + str(lambda_mass_fit_status))
    else:
        logging.info("\t\tLambda mass fit successful!")
        logging.info("\t\t\tChi2/NDF: " + str(lambda_mass_fit.GetChisquare()/lambda_mass_fit.GetNDF()))
        logging.info("\t\t\tMean: " + str(lambda_mass_fit.GetParameter(1)) + " +/- " + str(lambda_mass_fit.GetParError(1)))
        logging.info("\t\t\tSigma: " + str(lambda_mass_fit.GetParameter(2)) + " +/- " + str(lambda_mass_fit.GetParError(2)))
        logging.info("\t\t\tGamma: " + str(lambda_mass_fit.GetParameter(3)) + " +/- " + str(lambda_mass_fit.GetParError(3)))

    lambda_mass_bg_fit = rt.TF1("lambda_mass_bg_fit_" + unique_clone_string, "pol1", 1.096, 1.136)
    lambda_mass_bg_fit.SetNpx(1000)
    lambda_mass_bg_fit.SetParameter(0, bg_starting_params[0])
    lambda_mass_bg_fit.SetParameter(1, bg_starting_params[1])

    voigt_fit = rt.TF1("lambda_mass_voigt_fit_" + unique_clone_string,
                       "[0]*TMath::Voigt(x - [1], [2], [3], 4)", 1.102, 1.127)
    voigt_fit.SetNpx(1000)
    voigt_fit.SetParameter(0, lambda_mass_fit.GetParameter(0))
    voigt_fit.SetParameter(1, lambda_mass_fit.GetParameter(1))
    voigt_fit.SetParameter(2, lambda_mass_fit.GetParameter(2))
    voigt_fit.SetParameter(3, lambda_mass_fit.GetParameter(3))


    residual = lambda_mass_dist.Clone("lambda_mass_residual_" + unique_clone_string)
    residual.GetXaxis().SetRangeUser(1.094, 1.142)
    residual.Add(lambda_mass_bg_fit, -1)

    RSB_region = rt.TBox(sideband_bin.lower_bound, 0, sideband_bin.upper_bound, lambda_mass_dist.GetMaximum()*1.055)
    RSB_region.SetLineColor(rt.kRed)
    RSB_region.SetFillColor(rt.kRed)
    RSB_region.SetFillStyle(3003)

    RSB_min_line = rt.TLine(sideband_bin.lower_bound, 0, sideband_bin.upper_bound, lambda_mass_dist.GetMaximum()*1.05)
    RSB_min_line.SetLineColor(rt.kRed)
    RSB_min_line.SetLineWidth(2)
    RSB_min_line.SetLineStyle(2)

    RSB_max_line = rt.TLine(sideband_bin.lower_bound, 0, sideband_bin.upper_bound, lambda_mass_dist.GetMaximum()*1.05)
    RSB_max_line.SetLineColor(rt.kRed)
    RSB_max_line.SetLineWidth(2)
    RSB_max_line.SetLineStyle(2)

    left_signal_bin = lambda_mass_dist.FindBin(signal_bin.lower_bound)
    right_signal_bin = lambda_mass_dist.FindBin(signal_bin.upper_bound)

    left_full_region_bin = residual.FindBin(full_region_bin.lower_bound)
    right_full_region_bin = residual.FindBin(full_region_bin.upper_bound)

    lambda_bg = 0
    lambda_total = 0
    for bin_num in range(left_signal_bin, right_signal_bin + 1):
        lambda_total += lambda_mass_dist.GetBinContent(bin_num)
        lambda_bg += lambda_mass_bg_fit.Eval(lambda_mass_dist.GetBinCenter(bin_num))

    lambda_signal = lambda_total - lambda_bg
    signal_scale = residual.Integral(left_full_region_bin, right_full_region_bin)/residual.Integral(left_signal_bin, right_signal_bin)

    lambda_fit_info_string = unique_clone_string + "\n"
    lambda_fit_info_string += "lambda signal: " + str(lambda_signal) + "\n"
    lambda_fit_info_string += "lambda bg: " + str(lambda_bg) + "\n"
    lambda_fit_info_string += "lambda total: " + str(lambda_total) + "\n"
    lambda_fit_info_string += "signal scale: " + str(signal_scale) + "\n"

    lambda_fit_info_string = rt.TObjString(lambda_fit_info_string)
    
    if outfile:
        output_list = rt.TList()
        output_list.Add(lambda_mass_dist)
        output_list.Add(lambda_mass_dist_fit)
        output_list.Add(lambda_mass_bg_fit)
        output_list.Add(voigt_fit)
        output_list.Add(residual)
        output_list.Add(RSB_region)
        output_list.Add(RSB_min_line)
        output_list.Add(RSB_max_line)
        output_list.Add(lambda_fit_info_string)


        outfile.cd()
        outfile.WriteObject(output_list, "lambda_mass_fit_" + unique_clone_string)

    return lambda_signal, lambda_bg, signal_scale


# a function to perform the dihadron acceptance corrections
def make_h_h_mixed_corrections(same3d, mixed3d, unique_clone_string, outfile):

    logging.info("\t\t\tDoing dihadron mixed event corrections...")

    num_zbins = same3d.GetNbinsZ()
    same2d_uncorrected = same3d.Project3D("xye").Clone("h_h_2d_nomixcor_" + unique_clone_string)
    mixed2d_uncorrected = mixed3d.Project3D("xye").Clone("h_h_2d_mixed_" + unique_clone_string)

    for zbin in range(num_zbins):

        same3d.GetZaxis().SetRange(zbin+1, zbin+1)
        same2d = same3d.Project3D("xye")
        same2d.SetName(f"same2dproj_zbin_{zbin}")

        mixed3d.GetZaxis().SetRange(zbin+1, zbin+1)
        mixed2d = mixed3d.Project3D("xye")
        mixed2d.SetName(f"mix2dproj_zbin_{zbin}")

        scale = 0.5*(mixed2d.Integral(mixed2d.GetXaxis().FindBin(-0.01),    #xmin
                                    mixed2d.GetXaxis().FindBin(0.01),     #xmax 
                                    mixed2d.GetYaxis().FindBin(0.0),      #ymin
                                    mixed2d.GetYaxis().FindBin(0.0)))     #ymax
        same2d.Divide(mixed2d)
        same2d.Scale(scale)
        
        if zbin == 0:
            same2d_total = same2d.Clone("h_h_2d_mixcor_" + unique_clone_string)
        else:
            same2d_total.Add(same2d)

    if outfile:
        output_list = rt.TList()
        output_list.Add(same2d_uncorrected)
        output_list.Add(mixed2d_uncorrected)
        output_list.Add(same2d_total)

        outfile.cd()
        outfile.WriteObject(output_list, "h_h_mixed_corrections_" + unique_clone_string)
    
    logging.info("\t\t\tDone with dihadron mixed event corrections.")

    return same2d_total


# a function to perform h-lambda mixed event acceptance corrections
def make_h_lambda_mixed_corrections(same, mixed, signal_bin, sideband_bin, unique_clone_string, outfile):

    logging.info("\t\t\tDoing h-lambda mixed event corrections...")

    same.GetAxis(2).SetRangeUser(signal_bin.lower_bound, signal_bin.upper_bound)
    mixed.GetAxis(2).SetRangeUser(signal_bin.lower_bound, signal_bin.upper_bound)
    same3d_signal = same.Projection(0, 1, 3).Clone("same3d_signal")
    mixed3d_signal = mixed.Projection(0, 1, 3).Clone("mixed3d_signal")

    same.GetAxis(2).SetRangeUser(sideband_bin.lower_bound, sideband_bin.upper_bound)
    mixed.GetAxis(2).SetRangeUser(sideband_bin.lower_bound, sideband_bin.upper_bound)
    same3d_sideband = same.Projection(0, 1, 3).Clone("same3d_sideband")
    mixed3d_sideband = mixed.Projection(0, 1, 3).Clone("mixed3d_sideband")

    num_zbins = same3d_signal.GetNbinsZ()

    same2d_signal_uncorrected = same3d_signal.Project3D("xye").Clone("h_lambda_2d_nomixcor_" + unique_clone_string)
    mixed2d_signal_uncorrected = mixed3d_signal.Project3D("xye").Clone("h_lambda_2d_mixed_" + unique_clone_string)

    same2d_sideband_uncorrected = same3d_sideband.Project3D("xye").Clone("h_lambda_2d_nomixcor_" + unique_clone_string)
    mixed2d_sideband_uncorrected = mixed3d_sideband.Project3D("xye").Clone("h_lambda_2d_mixed_" + unique_clone_string)
    
    for zbin in range(num_zbins):

        same3d_signal.GetZaxis().SetRange(zbin+1, zbin+1)
        same2d_signal = same3d_signal.Project3D("xye")
        same2d_signal.SetName(f"same2dproj_signal_zbin_{zbin}")
        mixed3d_signal.GetZaxis().SetRange(zbin+1, zbin+1)
        mixed2d_signal = mixed3d_signal.Project3D("xye")
        mixed2d_signal.SetName(f"mix2dproj_signal_zbin_{zbin}")
        scale_signal = 0.5*(mixed2d_signal.Integral(mixed2d_signal.GetXaxis().FindBin(-0.01),    #xmin
                                    mixed2d_signal.GetXaxis().FindBin(0.01),     #xmax 
                                    mixed2d_signal.GetYaxis().FindBin(0.0),      #ymin
                                    mixed2d_signal.GetYaxis().FindBin(0.0)))     #ymax
        same2d_signal.Divide(mixed2d_signal)
        same2d_signal.Scale(scale_signal)

        same3d_sideband.GetZaxis().SetRange(zbin+1, zbin+1)
        same2d_sideband = same3d_sideband.Project3D("xye")
        same2d_sideband.SetName(f"same2dproj_sideband_zbin_{zbin}")
        mixed3d_sideband.GetZaxis().SetRange(zbin+1, zbin+1)
        mixed2d_sideband = mixed3d_sideband.Project3D("xye")
        mixed2d_sideband.SetName(f"mix2dproj_sideband_zbin_{zbin}")
        scale_sideband = 0.5*(mixed2d_sideband.Integral(mixed2d_sideband.GetXaxis().FindBin(-0.01),    #xmin
                                    mixed2d_sideband.GetXaxis().FindBin(0.01),     #xmax 
                                    mixed2d_sideband.GetYaxis().FindBin(0.0),      #ymin
                                    mixed2d_sideband.GetYaxis().FindBin(0.0)))     #ymax
        same2d_sideband.Divide(mixed2d_sideband)
        same2d_sideband.Scale(scale_sideband)
        
        if zbin == 0:
            same2d_signal_total = same2d_signal.Clone("h_lambda_2d_mixcor_" + unique_clone_string)
            same2d_sideband_total = same2d_sideband.Clone("h_lambda_2d_mixcor_"  + unique_clone_string)

        else:
            same2d_signal_total.Add(same2d_signal)
            same2d_sideband_total.Add(same2d_sideband)
    
    if outfile:
        output_list = rt.TList()
        output_list.Add(same2d_signal_uncorrected)
        output_list.Add(mixed2d_signal_uncorrected)
        output_list.Add(same2d_signal_total)
        output_list.Add(same2d_sideband_uncorrected)
        output_list.Add(mixed2d_sideband_uncorrected)
        output_list.Add(same2d_sideband_total)

        outfile.cd()
        outfile.WriteObject(output_list, "h_lambda_mixed_corrections_" + unique_clone_string)
    logging.info("\t\t\tDone with h-lambda mixed event corrections.")

    return same2d_signal_total, same2d_sideband_total

# a function that uses sideband distribution to subtract background from signal distribution
def apply_sideband_subtraction(signal_dist, sideband_dist, ratio, unique_clone_string, outfile):

    logging.info("\t\t\tApplying sideband subtraction...")

    subtracted_dist = signal_dist.Clone("h_lambda_2d_subtracted_" + unique_clone_string)
    sideband_dist.Scale(1/sideband_dist.Integral())
    bg_integral = (1 - ratio)*subtracted_dist.Integral()
    subtracted_dist.Add(sideband_dist, -bg_integral)

    logging.info("\t\t\tDone with sideband subtraction.")

    if outfile:
        output_list = rt.TList()
        output_list.Add(signal_dist)
        output_list.Add(sideband_dist)
        output_list.Add(subtracted_dist)

        outfile.cd()
        outfile.WriteObject(output_list, "h_lambda_sideband_subtraction_" + unique_clone_string)

    return subtracted_dist

# a function that applies two-track efficiency correction using MC generated templates
def apply_two_track_correction(uncorrected_dist, associated_pt_name):

    logging.info("\t\t\tApplying two-track efficiency correction...")

    # make copy of uncorrected dist before entering input file context
    corrected_dist = uncorrected_dist.Clone("h_lambda_2d_corrected_" + associated_pt_name)

    # get the correction template
    if associated_pt_name == "assoc_2_4":
        infile = rt.TFile("templates/twotrack_template.root", "READ")
    elif associated_pt_name == "assoc_15_25":
        infile = rt.TFile("templates/twotrack_template_lowpt.root", "READ")
    elif associated_pt_name == "assoc_25_4":
        infile = rt.TFile("templates/twotrack_template_highpt.root", "READ")
    else:
        print(f"ERROR: invalid associated pt bin for template correction: {associated_pt_name}")
        return
    correction_template = infile.Get("twotrack_template")

    deta_min_bin = correction_template.GetXaxis().FindBin(-1.4)
    deta_max_bin = correction_template.GetXaxis().FindBin(1.4 - 0.001)

    # loop over all bins in the input distribution
    for xbin in range(deta_min_bin, deta_max_bin+1):
        for ybin in range(uncorrected_dist.GetNbinsY()):
            # get the bin content of the input distribution
            bin_content = uncorrected_dist.GetBinContent(xbin, ybin+1)
            bin_error = uncorrected_dist.GetBinError(xbin, ybin+1)
            
            # get the correction factor from the correction template
            correction_factor = 1/correction_template.GetBinContent(xbin, ybin+1)
            
            # apply the correction factor to the bin content
            corrected_bin_content = bin_content * correction_factor
            corrected_bin_error = bin_error * correction_factor
            
            # set the bin content of the output distribution
            corrected_dist.SetBinContent(xbin, ybin+1, corrected_bin_content)
            corrected_dist.SetBinError(xbin, ybin+1, corrected_bin_error)
    
    infile.Close()

    logging.info("\t\t\tDone with two-track efficiency correction.")

    return corrected_dist

def get_systematic_uncertainty_dphi(variations, default, outfile):
    rms = 0
    n = 0
    for name, variation in variations.items():
        ratio = variation/default
        for i in range(1, ratio.GetNbinsX() + 1):
            rms += (ratio.GetBinContent(i) - 1)**2
            n += 1
    return math.sqrt(rms/n)

def fit_and_extract_yields(dphi_dist, fit_type, starting_params):

    if fit_type == fh.FitType.AVG_SIX:
        fit_dist = fh.fit_avg_six(dphi_dist, starting_params)
    elif fit_type == fh.FitType.AVG_FOUR:
        fit_dist = fh.fit_avg_four(dphi_dist, starting_params)
    elif fit_type == fh.FitType.ZYAM:
        fit_dist = fh.fit_zyam(dphi_dist, starting_params)
    elif fit_type == fh.FitType.V2:
        fit_dist = fh.fit_v2(dphi_dist, starting_params)
    elif fit_type == fh.FitType.DOUBLE_GAUS:
        fit_dist = fh.fit_double_gaus(dphi_dist, starting_params)
    elif fit_type == fh.FitType.VON_MISES:
        fit_dist = fh.fit_von_mises(dphi_dist, starting_params)
    else:
        logging.error("Invalid fit type: " + fit_type)
        return None

    return fit_dist

def do_all_fits(dphi_dist, trigger_pt_bin, associated_pt_bin, centrality_bin, is_dihadron=False):

    

    return True
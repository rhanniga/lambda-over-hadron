import math
import ROOT as rt

from enum import Enum

class FitType(Enum):
    AVG_SIX = 1
    AVG_FOUR = 2
    # V2 = 3
    # ZYAM = 4
    # DOUBLE_GAUS = 5
    # VON_MISES = 6

class YieldExtractor:

    def __init__(self, dphi_dist, 
                        centrality_string=None, 
                        trigger_pt_string=None, 
                        associated_pt_string=None, 
                        is_dihadron=False):

        self.dphi_dist = dphi_dist
        self.centrality_string = centrality_string
        self.trigger_pt_string = trigger_pt_string
        self.associated_pt_string = associated_pt_string
        self.is_dihadron = is_dihadron

        self.fits = rt.TList()
        self.yields = {}
    
    def extract_all_yields(self, list_name=None, output_file=None):

        for fit_type in FitType:
            self.extract_yield(fit_type)
        
        # self.save_fits(list_name, output_file)
    
    def extract_yield(self, fit_type):

        if fit_type == FitType.AVG_SIX:
            self.extract_yield_avg(n_bins=6, accept_negative_contributions=True)
        elif fit_type == FitType.AVG_FOUR:
            self.extract_yield_avg(n_bins=4, accept_negative_contributions=True)
        
    
    def extract_yield_avg(self, n_bins, accept_negative_contributions):

        # setup the fitting
        if n_bins == 4:
            avg_bins = [1, 8, 9, 16]
            fit_type = FitType.AVG_FOUR
            self.yields[fit_type] = {}
        elif n_bins == 6:
            avg_bins = [1, 2, 7, 8, 9, 16]
            fit_type = FitType.AVG_SIX
            self.yields[fit_type] = {}
        else:
            raise ValueError("n_bins must be 4 or 6 for now")

        ue_fit_name = f"avg_{n_bins}"
        ue_fit = rt.TF1(ue_fit_name, "[0]", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
        
        ue_avg = sum([self.dphi_dist.GetBinContent(bin) for bin in avg_bins])/n_bins
        ue_avg_err = math.sqrt(sum([self.dphi_dist.GetBinError(bin)**2 for bin in avg_bins]))/n_bins

        ue_fit.SetParameter(0, ue_avg)
        ue_fit.SetParError(0, ue_avg_err)

        self.fits.Add(ue_fit)

        # extract the yields from the fit
        ns_yield = 0
        ns_yield_err = 0

        as_yield = 0
        as_yield_err = 0

        ue_yield = ue_avg*self.dphi_dist.GetNbinsX()
        ue_yield_err = ue_avg_err*self.dphi_dist.GetNbinsX()

        total_yield = 0
        total_yield_err = 0

        for bin_num in range(1, self.dphi_dist.GetNbinsX() + 1):
            total_yield += self.dphi_dist.GetBinContent(bin_num)
            total_yield_err += self.dphi_dist.GetBinError(bin_num)**2

            part = self.dphi_dist.GetBinContent(bin_num) - ue_fit.Eval(self.dphi_dist.GetBinCenter(bin_num))
            if part < 0 and not accept_negative_contributions:
                part = 0
                continue
            if bin_num < 9:
                ns_yield += part
                ns_yield_err += self.dphi_dist.GetBinError(bin_num)**2
            else:
                as_yield += part
                as_yield_err += self.dphi_dist.GetBinError(bin_num)**2
        
        total_yield_err = math.sqrt(total_yield_err)
        ns_yield_err = math.sqrt(ns_yield_err)
        as_yield_err = math.sqrt(as_yield_err)

        ns_yield_err = math.sqrt(ns_yield_err**2 + (ue_yield_err/2)**2)
        as_yield_err = math.sqrt(as_yield_err**2 + (ue_yield_err/2)**2)
        
        self.yields[fit_type]["ns"] = (ns_yield, ns_yield_err)
        self.yields[fit_type]["as"] = (as_yield, as_yield_err)
        self.yields[fit_type]["ue"] = (ue_yield, ue_yield_err)
        self.yields[fit_type]["total"] = (total_yield, total_yield_err)
    

    def save_fits(self, list_name, output_file):
        output_file.WriteObject(self.fits, list_name)

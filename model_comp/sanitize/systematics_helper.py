import math
import ROOT as rt

class AnalysisBin:
    def __init__(self, name, lower_bound, upper_bound, file_name=None):
        self.name = name
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        if file_name:
            self.input_file = rt.TFile.Open(file_name)
        else:
            self.input_file = None
            self.object_list = None

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

    def close_file(self):
        if self.input_file:
            self.input_file.Close()
        
    
class AnalysisBins:
    def __init__(self, name, analysis_bins):
        self.name = name
        self.analysis_bins =  analysis_bins

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

class DphiSystematicHelper:
    def __init__(self, default_dist, variations):
        self.default_dist = default_dist
        self.variations = variations
        self.contributions = {}
        self.total_systematic = None
    
    def add_static_contribution(self, variation_name, contribution):
        self.contributions[variation_name] = contribution

    def get_contribution(self, variation_name):
        return self.contributions[variation_name]

    def get_rms(self, variation_name):

        if not self.variations[variation_name]:
            return 0

        rms = 0
        n = 0
        for variation in self.variations[variation_name]:
            ratio = variation/self.default_dist
            for i in range(1, ratio.GetNbinsX() + 1):
                rms += (ratio.GetBinContent(i) - 1)**2
                n += 1
        rms = math.sqrt(rms/n)
        self.contributions[variation_name] = rms
        return rms
    
    def get_total_systematic(self):
        total = 0
        for variation_name, contribution in self.contributions.items():
            total += contribution**2
        total = math.sqrt(total)
        self.total_systematic = total
        return total    

    def calculate_systematics(self):
        for variation_name in self.variations.keys():
            self.get_rms(variation_name)
        self.get_total_systematic()
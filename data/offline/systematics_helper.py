import math
import ROOT as rt

V2_FITPARS_FILEPATH = "v2_fitpars.txt"
# funny and inefficienct dict comprehension to practice my python prowess
# the exact ordering is: key: pedestal, pedestal error, trigger v2 (fixed to weighted avg), associated v2 (fixed to weighted avg)
V2_FITPARS_DICT = {
    line[0] : [float(line[1]), float(line[2]), float(line[3]), float(line[4])] for line in [line.split() for line in open(V2_FITPARS_FILEPATH).readlines()]
}

class AnalysisBin:
    def __init__(self, name, lower_bound, upper_bound, file_name=None):
        self.name = name
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        if file_name:
            self.input_file = rt.TFile.Open(file_name)
            self.object_list = self.input_file.Get("h-lambda")
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
    


class WidthSystematicHelper:
    def __init__(self, default_dist, variations, cent_name, trig_pt_name, assoc_pt_name, is_dihadron=False):
        self.default_dist = default_dist
        self.variations = variations
        self.cent_name = cent_name
        self.trig_pt_name = trig_pt_name
        self.assoc_pt_name = assoc_pt_name
        self.is_dihadron = is_dihadron

        self.default_width = self.extract_widths(default_dist)

        self.width_variations = {}
        self.contributions = {}
        self.total_systematic = {"ns_width": 0, "as_width": 0}

    def get_total_systematic(self):
        ns_total, as_total = 0, 0
        for variation_name, contribution in self.contributions.items():
            ns_total += contribution["ns_width"]**2
            as_total += contribution["as_width"]**2
        ns_total = math.sqrt(ns_total)
        as_total = math.sqrt(as_total)
        self.total_systematic["ns_width"] = ns_total
        self.total_systematic["as_width"] = as_total
        return ns_total, as_total    

    def get_rms(self, variation_name):
        self.contributions[variation_name] = {}
        if not self.width_variations[variation_name]:
            return 0, 0
        ns_rms, as_rms = 0, 0
        n = 0
        for variation in self.width_variations[variation_name]:
            ns_ratio = variation["ns_width"][0]/self.default_width["ns_width"][0]
            as_ratio = variation["as_width"][0]/self.default_width["as_width"][0]
            ns_rms += (ns_ratio - 1)**2
            as_rms += (as_ratio - 1)**2
            n += 1
        ns_rms = math.sqrt(ns_rms/n)
        as_rms = math.sqrt(as_rms/n)
        self.contributions[variation_name]["ns_width"] = ns_rms
        self.contributions[variation_name]["as_width"] = as_rms
        return ns_rms, as_rms

    def calculate_systematics(self):
        for variation_name in self.width_variations:
            self.get_rms(variation_name)
        self.get_total_systematic()

    def extract_all_widths(self):
        for variation_name, variation in self.variations.items():
            self.width_variations[variation_name] = []
            for dist in variation:
                self.width_variations[variation_name].append(self.extract_widths(dist))
        print(self.width_variations)

    def get_v2_parameters(self):
        if self.is_dihadron:
            key = "h_h_"
        else:
            key = "h_lambda_"
        key += self.cent_name + "_" + self.trig_pt_name + "_" + self.assoc_pt_name
        return V2_FITPARS_DICT[key]

    def get_width_from_kappa(self, kappa):
        return rt.TMath.Sqrt(-2*rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa)))

    def get_width_error_from_kappa(self, kappa, kappa_error):
        deriv = (rt.TMath.BesselI0(kappa)**2 + rt.TMath.BesselI0(kappa)*rt.TMath.BesselI(2, kappa) - 2*rt.TMath.BesselI1(kappa)**2)/(2*rt.TMath.Sqrt(2)*rt.TMath.BesselI0(kappa)*rt.TMath.BesselI1(kappa)*rt.TMath.Sqrt(-rt.TMath.Log(rt.TMath.BesselI1(kappa)/rt.TMath.BesselI0(kappa))))
        return deriv * kappa_error

    def extract_widths(self, dist, use_gaus=False):

        v2_parameters = self.get_v2_parameters()

        if use_gaus:
            fit_function_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
            fit_function_string += " + gaus(3) + gaus(6)"
            fit_function = rt.TF1("fit_function", fit_function_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

            # fixing v2 parameters (obtained from separate macro)
            fit_function.SetParameter(0, v2_parameters[0])
            fit_function.FixParameter(1, v2_parameters[2])
            fit_function.FixParameter(2, v2_parameters[3])
            # setting parameters for near-side gaussian (centered at 0)
            fit_function.SetParLimits(3, 0.0, 1)
            fit_function.SetParameter(3, 0.1)
            fit_function.FixParameter(4, 0)
            fit_function.SetParLimits(5, 0.1, 2)
            fit_function.SetParameter(5, 0.5)
            # setting parameters for near-side mirror gaussian (centered at 0 + 2pi)
            fit_function.SetParLimits(6, 0.0, 1)
            fit_function.SetParameter(6, 0.1)
            fit_function.FixParameter(7, rt.TMath.Pi())
            fit_function.SetParLimits(8, 0.1, 10)
            fit_function.SetParameter(8, 0.5)

        else:
            fit_function_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
            fit_function_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
            fit_function_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"
            fit_function = rt.TF1("fit_function", fit_function_string, -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)

            # fixing v2 parameters (obtained from separate macro)
            fit_function.SetParameter(0, v2_parameters[0])
            fit_function.FixParameter(1, v2_parameters[2])
            fit_function.FixParameter(2, v2_parameters[3])
            # setting parameters for near-side von-mises
            fit_function.SetParLimits(3, 0, 1)
            fit_function.SetParameter(3, 0.02)
            fit_function.SetParLimits(4, 0, 100)
            fit_function.SetParameter(4, 1)
            # setting parameters for away-side von-mises
            fit_function.SetParLimits(5, 0, 1)
            fit_function.SetParameter(5, 0.01)
            fit_function.SetParLimits(6, 0, 100)
            fit_function.SetParameter(6, 1)
        
        dist.Fit(fit_function, "BR")
        v2_comp = rt.TF1("v2_comp", "[0]*(1 + 2*([1]*[2]*cos(2*x)))", -rt.TMath.Pi()/2, 3*rt.TMath.Pi()/2)
        v2_comp.SetParameter(0, fit_function.GetParameter(0))
        v2_comp.SetParameter(1, fit_function.GetParameter(1))
        v2_comp.SetParameter(2, fit_function.GetParameter(2))

        dist.GetYaxis().SetRangeUser(0.8*dist.GetMinimum(), 1.2*dist.GetMaximum())
        test_c = rt.TCanvas("test_c", "test_c", 800, 600)
        dist.Draw()
        v2_comp.SetLineColor(rt.kBlue + 2)
        v2_comp.Draw("same")
        test_c.SaveAs("test.png")

        return_dict = {}
        if use_gaus:
            return_dict["ns_width"] = [fit_function.GetParameter(5), fit_function.GetParError(5)]
            return_dict["as_width"] = [fit_function.GetParameter(8), fit_function.GetParError(8)]
        else:
            return_dict["ns_width"] = [self.get_width_from_kappa(fit_function.GetParameter(4)), 
                                       self.get_width_error_from_kappa(fit_function.GetParameter(4), fit_function.GetParError(4))]
            return_dict["as_width"] = [self.get_width_from_kappa(fit_function.GetParameter(6)),
                                        self.get_width_error_from_kappa(fit_function.GetParameter(6), fit_function.GetParError(6))]
        
        return return_dict
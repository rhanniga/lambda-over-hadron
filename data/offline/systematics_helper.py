from ROOT import TFile

class AnalysisBin:
    def __init__(self, name, lower_bound, upper_bound, file_name=None):
        self.name = name
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        if file_name:
            self.input_file = TFile.Open(file_name)
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

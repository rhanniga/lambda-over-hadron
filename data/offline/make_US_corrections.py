import sys
from enum import IntEnum

import ROOT as rt

LSB_RANGE = []
SIGNAL_RANGE = []
RSB_RANGE = []

class Methods(IntEnum):
    FIT = 0
    LS = 1
    ROTATED = 2

def get_scale_factor(file_name, method=Methods.FIT, side="RIGHT"):
    i = rt.TFile(file_name)
    l = i.Get("h-lambda") 

    
    lsb = l.FindObject("")
    rsb = 
    signal =

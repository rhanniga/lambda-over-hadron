from ROOT import TFile, TCanvas
import array as arr
import configparser

####################################### CONFIGURATION START #######################################
config = configparser.ConfigParser()
config.read("mix_corrections.ini")
input_file = TFile(config['Head']['InputFile'])
object_list = input_file.Get(config['Head']['ListName'])
particle_name = config['Head']['Particle']
trig_pt_low = float(config['Pt']['TrigPtLow'])
trig_pt_high = float(config['Pt']['TrigPtHigh'])
assoc_pt_low = float(config['Pt']['AssocPtLow']) 
assoc_pt_high = float(config['Pt']['AssocPtHigh']) 
cent_low = float(config['Cent']['CentLow'])
cent_high = float(config['Cent']['CentHigh'])
peak_low = float(config['MassRegions']['PeakLow'])
peak_high = float(config['MassRegions']['PeakHigh'])
rsb_low = float(config['MassRegions']['RightSidebandLow'])
rsb_high = float(config['MassRegions']['RightSidebandHigh'])
lsb_low = float(config['MassRegions']['LeftSidebandLow'])
lsb_high = float(config['MassRegions']['LeftSidebandHigh'])
####################################### CONFIGURATION END #######################################

def reset_ranges(hist, dim):
    #can we just talk about how stupid it is that TH1,2,3D's don't have a GetAxis attr
    if dim == 1:
        hist.GetXaxis().SetRangeUser(0, 0)

    elif dim == 2:
        hist.GetXaxis().SetRangeUser(0, 0)
        hist.GetYaxis().SetRangeUser(0, 0)

    elif dim == 3:
        hist.GetXaxis().SetRangeUser(0, 0)
        hist.GetYaxis().SetRangeUser(0, 0)
        hist.GetZaxis().SetRangeUser(0, 0)

    else:
        for axis in range(dim):
            hist.GetAxis(axis).SetRangeUser(0, 0)

def make_mixed_corrections_full(same, mixed):
    reset_ranges(same, 4)
    reset_ranges(mixed, 4)

    same.GetAxis(2).SetRange(1, same.GetAxis(2).GetNbins())
    same.Sumw2()
    mixed.GetAxis(2).SetRange(1, mixed.GetAxis(2).GetNbins())
    mixed.Sumw2()

    for zbin in range(10):

        same.GetAxis(3).SetRange(zbin+1, zbin+1)
        same3D = same.Projection(0, 1, 2)
        same3D.SetName(f"same3dproj_zbin_{zbin}")

        mixed.GetAxis(3).SetRange(zbin+1, zbin+1)
        mixed3D = mixed.Projection(0, 1, 2)
        mixed3D.SetName(f"mix3dproj_zbin_{zbin}")

        #scaling by average of bins adjacent to 0
        scale = 0.5*(mixed3D.Integral(mixed3D.GetXaxis().FindBin(-0.01),    #xmin
                                      mixed3D.GetXaxis().FindBin(0.01),     #xmax 
                                      mixed3D.GetYaxis().FindBin(0.0),      #ymin
                                      mixed3D.GetYaxis().FindBin(0.0),      #ymax
                                      1,                                    #zmin
                                      mixed3D.GetNbinsZ()))                 #zmax

        same3D.Divide(mixed3D)
        same3D.Scale(scale)

        if zbin == 0:
            same3D_total = same3D.Clone("3dproj_total")
        else:
            same3D_total.Add(same3D)
    
    return same3D_total

def make_mixed_corrections(same, mixed, lowmass=1.11, highmass=1.12, is_hh=False):

    if not is_hh:
        reset_ranges(same, 4)
        reset_ranges(mixed, 4)

        same.GetAxis(2).SetRangeUser(lowmass, highmass)
        mixed.GetAxis(2).SetRangeUser(lowmass, highmass)
        same3D = same.Projection(0, 1, 3)
        same3D.Sumw2()
        mixed3D = mixed.Projection(0, 1, 3)
        mixed3D.Sumw2()
    else:
        reset_ranges(same, 3)
        reset_ranges(mixed, 3)

        same3D = same
        same3D.Sumw2()
        mixed3D = mixed
        mixed3D.Sumw2()

    for zbin in range(10):

        same3D.GetZaxis().SetRange(zbin+1, zbin+1)
        same2D = same3D.Project3D("xye")
        same2D.SetName(f"same2dproj_zbin_{zbin}")

        mixed3D.GetZaxis().SetRange(zbin+1, zbin+1)
        mixed2D = mixed3D.Project3D("xye")
        mixed2D.SetName(f"mix2dproj_zbin_{zbin}")

        #scaling by average of bins adjacent to 0
        scale = 0.5*(mixed2D.Integral(mixed2D.GetXaxis().FindBin(-0.01),    #xmin
                                      mixed2D.GetXaxis().FindBin(0.01),     #xmax 
                                      mixed2D.GetYaxis().FindBin(0.0),      #ymin
                                      mixed2D.GetYaxis().FindBin(0.0)))     #ymax

        same2D.Divide(mixed2D)
        same2D.Scale(scale)
        if zbin == 0:
            same2D_total = same2D.Clone("2dproj_total")
        else:
            same2D_total.Add(same2D)

    return same2D_total

def set_sparse_pt_ranges():

    #setting pt range on trigger distribution
    object_list.FindObject("fTriggerDist").GetAxis(0).SetRangeUser(trig_pt_low, trig_pt_high)

    #setting trigger pt range on relevant dphi distributions
    object_list.FindObject("fDphiHLambda").GetAxis(0).SetRangeUser(trig_pt_low, trig_pt_high)
    object_list.FindObject("fDphiHLambdaMixed").GetAxis(0).SetRangeUser(trig_pt_low, trig_pt_high)
    object_list.FindObject("fDphiHLambdaLS").GetAxis(0).SetRangeUser(trig_pt_low, trig_pt_high)
    object_list.FindObject("fDphiHLambdaLSMixed").GetAxis(0).SetRangeUser(trig_pt_low, trig_pt_high)

    #setting associated pt range on relevant dphi distributions
    object_list.FindObject("fDphiHLambda").GetAxis(1).SetRangeUser(assoc_pt_low, assoc_pt_high)
    object_list.FindObject("fDphiHLambdaMixed").GetAxis(1).SetRangeUser(assoc_pt_low, assoc_pt_high)
    object_list.FindObject("fDphiHLambdaLS").GetAxis(1).SetRangeUser(assoc_pt_low, assoc_pt_high)
    object_list.FindObject("fDphiHLambdaLSMixed").GetAxis(1).SetRangeUser(assoc_pt_low, assoc_pt_high)

def setup():

    set_sparse_pt_ranges()

    hist_dict = {}

    hist_dict["trig_dist"] = object_list.FindObject("fTriggerDist").Projection(0, 3)

    #need to pass an actual array to projection (not a list)
    axes = arr.array('i', [2, 3, 4, 5])
    hist_dict["dphi_h_strange"] = object_list.FindObject("fDphiHLambda").Projection(4, axes)
    hist_dict["dphi_h_strange_mixed"] = object_list.FindObject("fDphiHLambdaMixed").Projection(4, axes)
    hist_dict["dphi_h_strange_LS"] = object_list.FindObject("fDphiHLambdaLS").Projection(4, axes)
    hist_dict["dphi_h_strange_LS_mixed"] = object_list.FindObject("fDphiHLambdaLSMixed").Projection(4, axes)
    hist_dict["dphi_h_h"] = object_list.FindObject("fDphiHH").Projection(2, 3, 4)
    hist_dict["dphi_h_h"].Sumw2()
    hist_dict["dphi_h_h_mixed"] = object_list.FindObject("fDphiHHMixed").Projection(2, 3, 4)
    hist_dict["dphi_h_h_mixed"].Sumw2()

    return hist_dict

def generate_output_string():
    #didn't want to have this in one line
    out = "trig_"
    out += str(int(trig_pt_low)) + "_"
    out += str(int(trig_pt_high)) + "_"
    out += "assoc_"
    out += str(int(assoc_pt_low)) + "_"
    out += str(int(assoc_pt_high)) + "_"
    out += "cent_"
    out += str(int(cent_low)) + "_"
    out += str(int(cent_high)) + "_"
    out += "mixcorr_"
    out += "h" + particle_name + ".root"
    return out
    
def output_to_file(output_file_string, output_list):
    output_file = TFile(output_file_string, "RECREATE")
    output_file.cd()
    for thing in output_list:
        thing.Write()
    output_file.Close()

hist_dict = setup()

#making the corrections, appending corrected stuff to dictionary
hist_dict["hLambda3Dfull"] = make_mixed_corrections_full(hist_dict["dphi_h_strange"], hist_dict["dphi_h_strange_mixed"])
hist_dict["hLambda2Dfull_test"] = make_mixed_corrections(hist_dict["dphi_h_strange"], hist_dict["dphi_h_strange_mixed"], 1.06, 1.16-0.001)
hist_dict["hLambda2Dpeak"] = make_mixed_corrections(hist_dict["dphi_h_strange"], hist_dict["dphi_h_strange_mixed"], peak_low, peak_high)
hist_dict["hLambdaLS2Dpeak"] = make_mixed_corrections(hist_dict["dphi_h_strange_LS"], hist_dict["dphi_h_strange_LS_mixed"], peak_low, peak_high)
hist_dict["hLambda2DLside"] = make_mixed_corrections(hist_dict["dphi_h_strange"], hist_dict["dphi_h_strange_mixed"], lsb_low, lsb_high)
hist_dict["hLambdaLS2DLside"] = make_mixed_corrections(hist_dict["dphi_h_strange_LS"], hist_dict["dphi_h_strange_LS_mixed"], lsb_low, lsb_high)
hist_dict["hLambda2DRside"] = make_mixed_corrections(hist_dict["dphi_h_strange"], hist_dict["dphi_h_strange_mixed"], rsb_low, rsb_high)
hist_dict["hLambdaLS2DRside"] = make_mixed_corrections(hist_dict["dphi_h_strange_LS"], hist_dict["dphi_h_strange_LS_mixed"], rsb_low, rsb_high)
hist_dict["hh2D"] = make_mixed_corrections(hist_dict["dphi_h_h"], hist_dict["dphi_h_h_mixed"], is_hh=True) 

#list containing stuff to be written to output file
output_list = []

for hist in hist_dict:
    hist_dict[hist].SetName(hist)
    output_list.append(hist_dict[hist])

output_string = generate_output_string()
output_to_file(output_string, output_list)
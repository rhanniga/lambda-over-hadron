# simple script to project the multiplicity bins down into separate root files
import array as arr
import ROOT as rt

from sys import argv

epsilon = 0.0001
mult_low = 20
mult_high = 50 - epsilon

mult_string = "_20_50.root"

file_location = str(argv[1])
outfile_location = str(argv[2])
in_file = rt.TFile(file_location, "READ")

in_list = in_file.Get("h-lambda")

trig_dist = in_list.FindObject("fTriggerDistEff")
trig_dist_highest_pt = in_list.FindObject("fTriggerDistEff_highestPt")

lambda_dist = in_list.FindObject("fTriggeredLambdaDist")

h_h = in_list.FindObject("fDphiHHEff")
h_h_mixed = in_list.FindObject("fDphiHHMixed")
h_lambda = in_list.FindObject("fDphiHLambdaEff")
h_lambda_mixed = in_list.FindObject("fDphiHLambdaMixed")

h_h_highest_pt = in_list.FindObject("fDphiHHEff_highestPt")
h_h_mixed_highest_pt = in_list.FindObject("fDphiHHMixed_highestPt")
h_lambda_highest_pt = in_list.FindObject("fDphiHLambdaEff_highestPt")
h_lambda_mixed_highest_pt = in_list.FindObject("fDphiHLambdaMixed_highestPt")

hl_in_dists = [
    h_lambda,
    h_lambda_highest_pt,
    h_lambda_mixed,
    h_lambda_mixed_highest_pt
]

hh_in_dists = [
    h_h,
    h_h_highest_pt,
    h_h_mixed,
    h_h_mixed_highest_pt
]

out_list = rt.TList()
out_list.SetName("h-lambda")


for dist in [trig_dist, trig_dist_highest_pt, lambda_dist]:
    name = dist.GetName()
    dist.GetAxis(4).SetRangeUser(mult_low, mult_high)
    axes = arr.array('i', [0, 1, 2, 3])
    proj_dist = dist.Projection(4, axes)
    proj_dist.SetName(name)
    out_list.Add(proj_dist)

for dist in hl_in_dists:
    name = dist.GetName()
    dist.GetAxis(6).SetRangeUser(mult_low, mult_high)
    axes = arr.array('i', [0, 1, 2, 3, 4, 5])
    proj_dist = dist.Projection(6, axes)
    proj_dist.SetName(name)
    out_list.Add(proj_dist)

for dist in hh_in_dists:
    name = dist.GetName()
    dist.GetAxis(5).SetRangeUser(mult_low, mult_high)
    axes = arr.array('i', [0, 1, 2, 3, 4])
    proj_dist = dist.Projection(5, axes)
    proj_dist.SetName(name)
    out_list.Add(proj_dist)

out_file =  rt.TFile(outfile_location, "RECREATE")
out_file.WriteObject(out_list, "h-lambda")
out_file.Close()
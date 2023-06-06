import ROOT as rt
rt.gStyle.SetOptStat(0)

# Get the percentiles for the mult distribution ()
def get_mult_percentile_location(mult_dist, percentile):

    current_sum = 0
    total_sum = mult_dist.Integral(2, mult_dist.GetNbinsX())

    for b in range(mult_dist.GetNbinsX(), 0, -1):
        current_sum += mult_dist.GetBinContent(b)
        if current_sum/total_sum >= percentile:
            print("Percentile: ", percentile, " Location: ", mult_dist.GetBinCenter(b), " Bin: ", b)
            return b
    print("Given percentile not found")
    print("Percentile: ", percentile)
    return -1

dpmjet_infile = rt.TFile("output/dpmjet_17f2b_fast.root")
dpmjet_inlist = dpmjet_infile.Get("h-lambda")
dpmjet_mult_dist = dpmjet_inlist.FindObject("fMultDist")

cutoffs = [get_mult_percentile_location(dpmjet_mult_dist, 0.2),
           get_mult_percentile_location(dpmjet_mult_dist, 0.5),
           get_mult_percentile_location(dpmjet_mult_dist, 0.8),
           get_mult_percentile_location(dpmjet_mult_dist, 1.0),]

can = rt.TCanvas("c", "c", 800, 600)
can.SetLogy()

colored_hists = []
colors = [rt.kRed + 2, rt.kOrange + 2, rt.kBlue + 2, rt.kGray]
colors = colors[::-1]
test = rt.THStack("test", "")
for i, c in enumerate(cutoffs[::-1]):
    print(i, c)
    dpmjet_mult_dist.GetXaxis().SetRange(c, 250)
    dpmjet_mult_dist.SetTitle(";N_{ch} in V0A;N_{events}")
    dpmjet_mult_dist_clone = dpmjet_mult_dist.Clone(f"dpmjet_mult_dist_{c}")
    dpmjet_mult_dist_clone.SetFillColor(colors[i])
    test.Add(dpmjet_mult_dist_clone)

test.Draw("NOSTACK")
test.GetXaxis().SetLimits(2, 250)
can.SaveAs("dpmjet_mult_dist.png")

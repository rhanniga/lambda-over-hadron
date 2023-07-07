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

phsd_infile = rt.TFile("output/phsd_out_50mil_betterdists.root")
phsd_mult_dist = phsd_infile.Get("fMultDist")
dpmjet_mult_dist = phsd_infile.Get("fMultDist")

dpmjet_cutoffs = [get_mult_percentile_location(dpmjet_mult_dist, 0.2),
           get_mult_percentile_location(dpmjet_mult_dist, 0.5),
           get_mult_percentile_location(dpmjet_mult_dist, 0.8),
           get_mult_percentile_location(dpmjet_mult_dist, 1.0),]

phsd_cutoffs = [get_mult_percentile_location(phsd_mult_dist, 0.2),
           get_mult_percentile_location(phsd_mult_dist, 0.5),
           get_mult_percentile_location(phsd_mult_dist, 0.8),
           get_mult_percentile_location(phsd_mult_dist, 1.0),]

can = rt.TCanvas("c", "c", 800, 600)
can.SetLogy()

colored_hists = []
colors = [rt.kRed + 2, rt.kOrange + 2, rt.kBlue + 2, rt.kGray]
colors = colors[::-1]

# for i, c in enumerate(phsd_cutoffs[::-1]):
#     print("--------phsd--------")
#     print(i, c)

dpmjet_mult_dist.SetTitle(";N_{ch} in V0A;N_{events}")
dpmjet_mult_dist.GetXaxis().SetRange(2, 250)
dpmjet_mult_dist.Draw()
test_dist1 = rt.TH1D("test_dist1", "test_dist", 100, 0, 100)
test_dist2 = rt.TH1D("test_dist2", "test_dist", 100, 0, 100)
test_dist3 = rt.TH1D("test_dist3", "test_dist", 100, 0, 100)
test_dist4 = rt.TH1D("test_dist4", "test_dist", 100, 0, 100)
leg = rt.TLegend(0.6, 0.6, 0.9, 0.9)
for i, c in enumerate(dpmjet_cutoffs[::-1]):

    dpmjet_mult_dist_clone = dpmjet_mult_dist.Clone(f"dpmjet_mult_dist_{c}")
    dpmjet_mult_dist_clone.GetXaxis().SetRange(c, 250)
    dpmjet_mult_dist_clone.SetFillColor(colors[i])
    dpmjet_mult_dist_clone.DrawCopy("same")

    if i == 0:
        test_dist1.SetFillColor(colors[i])
        leg.AddEntry(test_dist1, "80-100%", "f")
    elif i == 1:
        test_dist2.SetFillColor(colors[i])
        leg.AddEntry(test_dist2, "50-80%", "f")
    elif i == 2:
        test_dist3.SetFillColor(colors[i])
        leg.AddEntry(test_dist3, "20-50%", "f")
    elif i == 3:
        test_dist4.SetFillColor(colors[i])
        leg.AddEntry(test_dist4, "0-20%", "f")

leg.Draw()
can.SaveAs("dpmjet_mult_dist.png")



import ROOT as rt
rt.gStyle.SetOptStat(0)

dpmjet_infile = rt.TFile("output/AnalysisResults.root")
dpmjet_inlist = dpmjet_infile.Get("h-lambda")

mult_dist = dpmjet_inlist.FindObject("fMultDist")
num_events = mult_dist.Integral(2, mult_dist.GetNbinsX())/20

cent_0_20_events = num_events*0.2
cent_20_50_events = num_events*0.3
cent_50_80_events = num_events*0.3

hadron_dist = dpmjet_inlist.FindObject("fTriggeredTriggerDist_MC")
test_dist = dpmjet_inlist.FindObject("fTriggerDist_MC")


hadron_dist.GetAxis(0).SetRangeUser(0.15, 100)

for n_events, cent_bin in zip([cent_0_20_events, cent_20_50_events, cent_50_80_events], [4, 3, 2]):
    hadron_dist.GetAxis(4).SetRange(cent_bin, cent_bin)
    hadron_eta_dist = hadron_dist.Projection(2)
    c = rt.TCanvas("c", "c", 800, 600)
    c.cd()
    hadron_eta_dist.GetXaxis().SetRangeUser(-1, 1-0.0001)
    hadron_eta_dist.Scale(1/n_events, "width")

    left_bin = hadron_eta_dist.FindBin(-0.5)
    right_bin = hadron_eta_dist.FindBin(0.5-0.0001)
    integral = 0
    nbins = 0
    for b in range(left_bin, right_bin+1):
        integral += hadron_eta_dist.GetBinContent(b)
        nbins += 1
    avg = integral/nbins
    
    print("Centrality bin: {}".format(cent_bin))
    print("dN/dEta: {}".format(avg))


import ROOT as rt
rt.gStyle.SetOptStat(0)

c = rt.TCanvas("c", "", 800, 600)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.10)
c.SetTopMargin(0.05)

BRANCHING_RATIO = 0.639

data_infile = rt.TFile("output/data_infile.root")
phsd_infile = rt.TFile("output/phsd_out_50mil_betterdists.root")
epos_infile = rt.TFile("output/epos_17f2a_fast.root")
dpmjet_infile = rt.TFile("output/dpmjet_17f2b_fast.root")

phsd_mult_dist = phsd_infile.Get("fMultDist")
phsd_nevents = phsd_mult_dist.Integral()
epos_mult_dist = epos_infile.Get("h-lambda").FindObject("fMultDist")
epos_nevents = epos_mult_dist.Integral()
dpmjet_mult_dist = dpmjet_infile.Get("h-lambda").FindObject("fMultDist")
dpmjet_nevents = dpmjet_mult_dist.Integral()

data_nevents = 248341207

print(f"Num events phsd: {phsd_nevents}")
print(f"Num events epos: {epos_nevents}")
print(f"Num events dpmjet: {dpmjet_nevents}")

phsd_trigger_dist = phsd_infile.Get("fTriggerDist_MC")
epos_trigger_dist = epos_infile.Get("h-lambda").FindObject("fTriggerDist_MC")
dpmjet_trigger_dist = dpmjet_infile.Get("h-lambda").FindObject("fTriggerDist_MC")

phsd_lambda_dist = phsd_infile.Get("fLambdaDist_MC")
epos_lambda_dist = epos_infile.Get("h-lambda").FindObject("fLambdaDist_MC")
dpmjet_lambda_dist = dpmjet_infile.Get("h-lambda").FindObject("fLambdaDist_MC")

phsd_lambda_dist_no_eta_cuts = phsd_infile.Get("fLambdaDist_MC_no_eta_cut")
epos_lambda_dist_no_eta_cuts = epos_infile.Get("h-lambda").FindObject("fLambdaDist_MC_no_eta_cut")
dpmjet_lambda_dist_no_eta_cuts = dpmjet_infile.Get("h-lambda").FindObject("fLambdaDist_MC_no_eta_cut")

triggered_phsd_lambda_dist = phsd_infile.Get("fTriggeredLambdaDist_MC")
triggered_epos_lambda_dist = epos_infile.Get("h-lambda").FindObject("fTriggeredLambdaDist_MC")
triggered_dpmjet_lambda_dist = dpmjet_infile.Get("h-lambda").FindObject("fTriggeredLambdaDist_MC")


phsd_mult_dist.GetXaxis().SetRange(1, 250)
phsd_mult_dist.SetLineColor(rt.kRed + 2)
phsd_mult_dist.SetLineWidth(2)
phsd_mult_dist.SetTitle("")
phsd_mult_dist.GetXaxis().SetTitle("N_{ch} in V0A acceptance")
phsd_mult_dist.GetYaxis().SetTitle("N_{events}")
phsd_mult_dist.Scale(epos_nevents/phsd_nevents)

epos_mult_dist.GetXaxis().SetRange(1, 250)
epos_mult_dist.SetLineColor(rt.kGreen + 2)
epos_mult_dist.SetLineWidth(2)
epos_mult_dist.SetTitle("")
epos_mult_dist.GetXaxis().SetTitle("N_{ch} in V0A acceptance")
epos_mult_dist.GetYaxis().SetTitle("N_{events}")
epos_mult_dist.Scale(epos_nevents/epos_nevents)

dpmjet_mult_dist.GetXaxis().SetRange(1, 250)
dpmjet_mult_dist.SetLineColor(rt.kBlue + 2)
dpmjet_mult_dist.SetLineWidth(2)
dpmjet_mult_dist.SetTitle("")
dpmjet_mult_dist.GetXaxis().SetTitle("N_{ch} in V0A acceptance")
dpmjet_mult_dist.GetYaxis().SetTitle("N_{events}")
dpmjet_mult_dist.Scale(epos_nevents/dpmjet_nevents)


c.SetLogy()

phsd_mult_dist.Draw()
epos_mult_dist.Draw("same")
dpmjet_mult_dist.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(phsd_mult_dist, "PHSD", "l")
leg.AddEntry(epos_mult_dist, "EPOS", "l")
leg.AddEntry(dpmjet_mult_dist, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/mult_dist.pdf")


data_trigger_pt = data_infile.Get("trig_pt_dist_0_80")
phsd_trigger_pt = phsd_trigger_dist.Projection(0).Clone("phsd_trigger_pt")
epos_trigger_pt = epos_trigger_dist.Projection(0).Clone("epos_trigger_pt")
dpmjet_trigger_pt = dpmjet_trigger_dist.Projection(0).Clone("dpmjet_trigger_pt")

data_trigger_pt.SetLineColor(rt.kBlack)
data_trigger_pt.SetLineWidth(2)
data_trigger_pt.SetTitle("")
data_trigger_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
data_trigger_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
data_trigger_pt.Scale(1/data_nevents)

phsd_trigger_pt.SetLineColor(rt.kRed + 2)
phsd_trigger_pt.SetLineWidth(2)
phsd_trigger_pt.SetTitle("")
phsd_trigger_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
phsd_trigger_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
phsd_trigger_pt.Scale(1/phsd_nevents)

epos_trigger_pt.SetLineColor(rt.kGreen + 2)
epos_trigger_pt.SetLineWidth(2)
epos_trigger_pt.SetTitle("")
epos_trigger_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
epos_trigger_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
epos_trigger_pt.Scale(1/epos_nevents)

dpmjet_trigger_pt.SetLineColor(rt.kBlue + 2)
dpmjet_trigger_pt.SetLineWidth(2)
dpmjet_trigger_pt.SetTitle("")
dpmjet_trigger_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
dpmjet_trigger_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
dpmjet_trigger_pt.Scale(1/dpmjet_nevents)

phsd_trigger_pt.Draw()
data_trigger_pt.Draw("same")
epos_trigger_pt.Draw("same")
dpmjet_trigger_pt.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(data_trigger_pt, "Data", "l")
leg.AddEntry(phsd_trigger_pt, "PHSD", "l")
leg.AddEntry(epos_trigger_pt, "EPOS", "l")
leg.AddEntry(dpmjet_trigger_pt, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/pt_dist.pdf")

phsd_trigger_pt.GetXaxis().SetRangeUser(4, 7.999)
epos_trigger_pt.GetXaxis().SetRangeUser(4, 7.999)
dpmjet_trigger_pt.GetXaxis().SetRangeUser(4, 7.999)

c.SetLogy(0)

epos_trigger_pt.Draw()
data_trigger_pt.Draw("same")
phsd_trigger_pt.Draw("same")
dpmjet_trigger_pt.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(data_trigger_pt, "Data", "l")
leg.AddEntry(phsd_trigger_pt, "PHSD", "l")
leg.AddEntry(epos_trigger_pt, "EPOS", "l")
leg.AddEntry(dpmjet_trigger_pt, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/pt_dist_zoom.pdf")

phsd_ntriggers = phsd_trigger_pt.Integral()
epos_ntriggers = epos_trigger_pt.Integral()
dpmjet_ntriggers = dpmjet_trigger_pt.Integral()

print("PHSD trig/event: ", phsd_ntriggers)
print("EPOS trig/event: ", epos_ntriggers)
print("DPMJET trig/event: ", dpmjet_ntriggers)

phsd_lambda_pt = phsd_lambda_dist.Projection(0).Clone("phsd_lambda_pt")
epos_lambda_pt = epos_lambda_dist.Projection(0).Clone("epos_lambda_pt")
dpmjet_lambda_pt = dpmjet_lambda_dist.Projection(0).Clone("dpmjet_lambda_pt")

phsd_lambda_pt.SetLineColor(rt.kRed + 2)
phsd_lambda_pt.SetLineWidth(2)
phsd_lambda_pt.SetTitle("")
phsd_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
phsd_lambda_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
phsd_lambda_pt.Scale(1/phsd_nevents)

epos_lambda_pt.SetLineColor(rt.kGreen + 2)
epos_lambda_pt.SetLineWidth(2)
epos_lambda_pt.SetTitle("")
epos_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
epos_lambda_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
epos_lambda_pt.Scale(1/(epos_nevents*BRANCHING_RATIO))

dpmjet_lambda_pt.SetLineColor(rt.kBlue + 2)
dpmjet_lambda_pt.SetLineWidth(2)
dpmjet_lambda_pt.SetTitle("")
dpmjet_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
dpmjet_lambda_pt.GetYaxis().SetTitle("1/N_{events} dN/dp_{T} (GeV/c)^{-1}")
dpmjet_lambda_pt.Scale(1/(dpmjet_nevents*BRANCHING_RATIO))

epos_lambda_pt.Draw()
phsd_lambda_pt.Draw("same")
dpmjet_lambda_pt.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(phsd_lambda_pt, "PHSD", "l")
leg.AddEntry(epos_lambda_pt, "EPOS", "l")
leg.AddEntry(dpmjet_lambda_pt, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/lambda_pt_dist.pdf")

phsd_lambda_pt.GetXaxis().SetRangeUser(1.5, 3.999)
epos_lambda_pt.GetXaxis().SetRangeUser(1.5, 3.999)
dpmjet_lambda_pt.GetXaxis().SetRangeUser(1.5, 3.999)

c.SetLogy(0)

epos_lambda_pt.Draw()
phsd_lambda_pt.Draw("same")
dpmjet_lambda_pt.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(phsd_lambda_pt, "PHSD", "l")
leg.AddEntry(epos_lambda_pt, "EPOS", "l")
leg.AddEntry(dpmjet_lambda_pt, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/lambda_pt_dist_zoom.pdf")

phsd_nlambdas = phsd_lambda_pt.Integral()
epos_nlambdas = epos_lambda_pt.Integral()
dpmjet_nlambdas = dpmjet_lambda_pt.Integral()

print("PHSD lambda/event: ", phsd_nlambdas)
print("EPOS lambda/event: ", epos_nlambdas)
print("DPMJET lambda/event: ", dpmjet_nlambdas)

phsd_lambda_eta = phsd_lambda_dist_no_eta_cuts.Projection(2).Clone("phsd_lambda_eta")
epos_lambda_eta = epos_lambda_dist_no_eta_cuts.Projection(2).Clone("epos_lambda_eta")
dpmjet_lambda_eta = dpmjet_lambda_dist_no_eta_cuts.Projection(2).Clone("dpmjet_lambda_eta")

phsd_lambda_eta.SetLineColor(rt.kRed + 2)
phsd_lambda_eta.SetLineWidth(2)
phsd_lambda_eta.SetTitle("")
phsd_lambda_eta.GetXaxis().SetTitle("#eta")
phsd_lambda_eta.GetYaxis().SetTitle("1/N_{events} dN/d#eta")
phsd_lambda_eta.Scale(1/phsd_nevents)

epos_lambda_eta.SetLineColor(rt.kGreen + 2)
epos_lambda_eta.SetLineWidth(2)
epos_lambda_eta.SetTitle("")
epos_lambda_eta.GetXaxis().SetTitle("#eta")
epos_lambda_eta.GetYaxis().SetTitle("1/N_{events} dN/d#eta")
epos_lambda_eta.Scale(1/(epos_nevents*BRANCHING_RATIO))

dpmjet_lambda_eta.SetLineColor(rt.kBlue + 2)
dpmjet_lambda_eta.SetLineWidth(2)
dpmjet_lambda_eta.SetTitle("")
dpmjet_lambda_eta.GetXaxis().SetTitle("#eta")
dpmjet_lambda_eta.GetYaxis().SetTitle("1/N_{events} dN/d#eta")
dpmjet_lambda_eta.Scale(1/(dpmjet_nevents*BRANCHING_RATIO))

epos_lambda_eta.Draw()
phsd_lambda_eta.Draw("same")
dpmjet_lambda_eta.Draw("same")

leg = rt.TLegend(0.75, 0.75, 0.9, 0.9)
leg.AddEntry(phsd_lambda_eta, "PHSD", "l")
leg.AddEntry(epos_lambda_eta, "EPOS", "l")
leg.AddEntry(dpmjet_lambda_eta, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/lambda_eta_dist.pdf")


triggered_phsd_lambda_pt = triggered_phsd_lambda_dist.Projection(0).Clone("triggered_phsd_lambda_pt")
triggered_epos_lambda_pt = triggered_epos_lambda_dist.Projection(0).Clone("triggered_epos_lambda_pt")
triggered_dpmjet_lambda_pt = triggered_dpmjet_lambda_dist.Projection(0).Clone("triggered_dpmjet_lambda_pt")

phsd_ntriggers *= phsd_nevents
epos_ntriggers *= epos_nevents
dpmjet_ntriggers *= dpmjet_nevents

triggered_phsd_lambda_pt.SetLineColor(rt.kRed + 2)
triggered_phsd_lambda_pt.SetLineWidth(2)
triggered_phsd_lambda_pt.SetTitle("")
triggered_phsd_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
triggered_phsd_lambda_pt.GetYaxis().SetTitle("1/N_{trig} dN/dp_{T} (GeV/c)^{-1}")
triggered_phsd_lambda_pt.Scale(1/phsd_ntriggers)

triggered_epos_lambda_pt.SetLineColor(rt.kGreen + 2)
triggered_epos_lambda_pt.SetLineWidth(2)
triggered_epos_lambda_pt.SetTitle("")
triggered_epos_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
triggered_epos_lambda_pt.GetYaxis().SetTitle("1/N_{trig} dN/dp_{T} (GeV/c)^{-1}")
triggered_epos_lambda_pt.Scale(1/(epos_ntriggers*BRANCHING_RATIO))

triggered_dpmjet_lambda_pt.SetLineColor(rt.kBlue + 2)
triggered_dpmjet_lambda_pt.SetLineWidth(2)
triggered_dpmjet_lambda_pt.SetTitle("")
triggered_dpmjet_lambda_pt.GetXaxis().SetTitle("p_{T} (GeV/c)")
triggered_dpmjet_lambda_pt.GetYaxis().SetTitle("1/N_{trig} dN/dp_{T} (GeV/c)^{-1}")
triggered_dpmjet_lambda_pt.Scale(1/(dpmjet_ntriggers*BRANCHING_RATIO))

triggered_epos_lambda_pt.Draw()
triggered_phsd_lambda_pt.Draw("same")
triggered_dpmjet_lambda_pt.Draw("same")

leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(triggered_phsd_lambda_pt, "PHSD", "l")
leg.AddEntry(triggered_epos_lambda_pt, "EPOS", "l")
leg.AddEntry(triggered_dpmjet_lambda_pt, "DPMJET", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw("same")

c.SaveAs("figures/triggered_lambda_pt_dist.pdf")

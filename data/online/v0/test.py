import ROOT as rt


c = rt.TCanvas("c", "c", 800, 600)

infile = rt.TFile("AnalysisResults.root")

inlist = infile.Get("h-lambda")

indist = inlist.FindObject("fEventSelection")

indist.GetYaxis().SetRangeUser(0.0, 19.99)
indist = indist.ProjectionX()

print(indist.GetBinContent(3)/indist.GetBinContent(1))
indist.Draw()
c.SaveAs("test.png")

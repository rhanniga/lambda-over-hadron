import ROOT as rt


c = rt.TCanvas("c", "c", 800, 600)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.10)

colors = [rt.kRed - 3, rt.kBlue - 3, rt.kGreen - 3, rt.kMagenta - 3, rt.kCyan - 3, rt.kOrange - 3, rt.kYellow - 3, rt.kGray + 2, rt.kBlack]

outfile = rt.TFile("h_phi_fits.root", "RECREATE")

# for filename in ["fitsyst_hightest.root", "fitsyst_lowtest.root", "fitsyst_2_4.root"]:
for filename in ["fitsyst_2_4.root"]:

    if "high" in filename:
        pt_string = "_highpt"
    elif "low" in filename:
        pt_string = "_lowpt"
    else:
        pt_string = ""
    

    infile = rt.TFile(filename, "READ")

    color_index = 0
    for cent_bin in ["_0_20", "_20_50", "_50_80"]:

        h_phi_dphi = infile.Get("hPhidphi_Eff" + cent_bin)

        h_phi_dphi.SetTitle("")
        h_phi_dphi.SetLineColor(colors[color_index])
        h_phi_dphi.SetMarkerColor(colors[color_index])
        
        h_phi_dphi.GetYaxis().SetRangeUser(0.8*h_phi_dphi.GetMinimum(), 1.3*h_phi_dphi.GetMaximum())
        h_phi_dphi.Draw()
        c.SaveAs("figures/h_phi_dphi" + cent_bin + pt_string + ".pdf")

        von_fit_string = "[0]*(1 + 2*([1]*[2]*cos(2*x)))"
        von_fit_string += " + [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))"
        von_fit_string += " + [5]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- TMath::Pi()))"

        von_fit = rt.TF1("von_fit", von_fit_string, -2, 6)
        von_fit.SetNpx(1000)

        if "_0_20" in cent_bin:
            von_fit.FixParameter(1, 0.15)
            von_fit.FixParameter(2, 0.15)

            von_fit.SetParameter(0, 0.006)
            von_fit.SetParLimits(0, 0.005, 0.008)
            von_fit.SetParLimits(3, 0, 0.01)
            von_fit.SetParameter(3, 0.02)
            von_fit.SetParLimits(4, 0, 100)
            von_fit.SetParameter(4, 1)
            von_fit.SetParLimits(5, 0, 0.01)
            von_fit.SetParameter(5, 0.01)
            von_fit.SetParLimits(6, 0, 100)
            von_fit.SetParameter(6, 1)
        elif "_20_50" in cent_bin:
            von_fit.FixParameter(1, 0.15*0.85)
            von_fit.FixParameter(2, 0.15*0.85)

            von_fit.SetParameter(0, 0.003)
            von_fit.SetParLimits(0, 0.003, 0.00368)
            von_fit.SetParLimits(3, 0.0001, 0.0008)
            von_fit.SetParameter(3, 0.0002)
            von_fit.SetParLimits(4, 0.1, 50)
            von_fit.SetParameter(4, 1)
            von_fit.SetParLimits(5, 0.0001, 0.0015)
            von_fit.SetParameter(5, 0.0001)
            von_fit.SetParLimits(6, 0.1, 50)
            von_fit.SetParameter(6, 1)

        elif "_50_80" in cent_bin:
            von_fit.FixParameter(1, 0.15*0.50)
            von_fit.FixParameter(2, 0.15*0.50)

            von_fit.SetParameter(0, 0.002)
            von_fit.SetParLimits(0, 0.0015, 0.0025)
            von_fit.SetParLimits(3, 0.0001, 0.001)
            von_fit.SetParameter(3, 0.0002)
            von_fit.SetParLimits(4, 0.1, 15)
            von_fit.SetParameter(4, 10)
            von_fit.SetParLimits(5, 0.0001, 0.001)
            von_fit.SetParameter(5, 0.0002)
            von_fit.SetParLimits(6, 0.1, 40)
            von_fit.SetParameter(6, 10)


        outfile.cd()
        h_phi_dphi.Fit("von_fit", "R")

        c.cd()


        von_fit.SetLineColor(colors[color_index])
        von_fit.SetLineWidth(2)
        von_fit.Draw("same")

        von_fit.Write("h_phi_von_fit" + cent_bin + pt_string)

        c.SaveAs("figures/h_phi_dphi_fit" + cent_bin + pt_string + ".pdf")
        c.SaveAs("figures/h_phi_dphi_fit" + cent_bin + pt_string + ".png")

        color_index += 1

outfile.Close()
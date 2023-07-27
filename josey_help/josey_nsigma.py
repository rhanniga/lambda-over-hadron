import ROOT as rt

c = rt.TCanvas("c", "", 800, 600)

bethe_bloch_nsigma_electron = rt.TH1D("bethe_bloch_nsigma_electron", ";p_{T} (GeV/c);n#sigma_{e}", 25, 0.0, 10.0)
bethe_bloch_nsigma_pion = rt.TH1D("bethe_bloch_nsigma_pion", ";p_{T} (GeV/c);n#sigma_{e}", 25, 0.0, 10.0)
bethe_bloch_nsigma_kaon = rt.TH1D("bethe_bloch_nsigma_kaon", ";p_{T} (GeV/c);n#sigma_{e}", 25, 0.0, 10.0)
bethe_bloch_nsigma_proton = rt.TH1D("bethe_bloch_nsigma_proton", ";p_{T} (GeV/c);n#sigma_{e}", 25, 0.0, 10.0)

def main(filepath="AnalysisResults.root"):

    infile = rt.TFile(filepath)

    inlist = infile.Get("nsigma-electron")

    big_dist = inlist.FindObject("fBigDist")

    for momentum_bin in range(1, big_dist.GetAxis(0).GetNbins() + 1, 4):

        big_dist.GetAxis(0).SetRange(momentum_bin, momentum_bin + 4)

        zero_nsigma_electron_bin = big_dist.GetAxis(2).FindBin(0.0)
        big_dist.GetAxis(2).SetRange(zero_nsigma_electron_bin, zero_nsigma_electron_bin)

        tpc_signal = big_dist.Projection(1)
        tpc_average_elec_0 = tpc_signal.GetMean(1)

        one_nsigma_electron_bin = big_dist.GetAxis(2).FindBin(1.0)
        big_dist.GetAxis(2).SetRange(one_nsigma_electron_bin, one_nsigma_electron_bin)

        tpc_signal = big_dist.Projection(1)
        tpc_average_elec_1 = tpc_signal.GetMean()

        if not (tpc_average_elec_0 > 0 and tpc_average_elec_1 > 0):
            continue

        width = tpc_average_elec_1 - tpc_average_elec_0

        big_dist.GetAxis(2).SetRange(0, -1)

        zero_nsigma_pion_bin = big_dist.GetAxis(3).FindBin(0.0)
        big_dist.GetAxis(3).SetRange(zero_nsigma_pion_bin, zero_nsigma_pion_bin)
        tpc_signal = big_dist.Projection(1)
        tpc_average_pion_0 = tpc_signal.GetMean()
        nsigma_pion = (tpc_average_pion_0 - tpc_average_elec_0) / width
        bethe_bloch_nsigma_pion.SetBinContent(int(momentum_bin/4) + 1, nsigma_pion)
        big_dist.GetAxis(3).SetRange(0, -1)

        zero_nsigma_kaon_bin = big_dist.GetAxis(4).FindBin(0.0)
        big_dist.GetAxis(4).SetRange(zero_nsigma_kaon_bin, zero_nsigma_kaon_bin)
        tpc_signal = big_dist.Projection(1)
        tpc_average_kaon_0 = tpc_signal.GetMean()
        nsigma_kaon = (tpc_average_kaon_0 - tpc_average_elec_0) / width
        bethe_bloch_nsigma_kaon.SetBinContent(int(momentum_bin/4)+1, nsigma_kaon)
        big_dist.GetAxis(4).SetRange(0, -1)

        zero_nsigma_proton_bin = big_dist.GetAxis(5).FindBin(0.0)
        big_dist.GetAxis(5).SetRange(zero_nsigma_proton_bin, zero_nsigma_proton_bin)
        tpc_signal = big_dist.Projection(1)
        tpc_average_proton_0 = tpc_signal.GetMean()
        nsigma_proton = (tpc_average_proton_0 - tpc_average_elec_0) / width
        bethe_bloch_nsigma_proton.SetBinContent(int(momentum_bin/4) + 1, nsigma_proton)
        big_dist.GetAxis(5).SetRange(0, -1)
        

    big_dist.GetAxis(0).SetRange(0,-1)
    nsigma_electron = big_dist.Projection(2, 0)

    bethe_bloch_nsigma_pion.SetLineColor(rt.kRed)
    bethe_bloch_nsigma_kaon.SetLineColor(rt.kMagenta)
    bethe_bloch_nsigma_proton.SetLineColor(rt.kGreen)

    bethe_bloch_nsigma_pion.SetLineWidth(2)
    bethe_bloch_nsigma_kaon.SetLineWidth(2)
    bethe_bloch_nsigma_proton.SetLineWidth(2)

    nsigma_electron.GetXaxis().SetRangeUser(2, 6)
    nsigma_electron.Draw("COLZ")
    bethe_bloch_nsigma_pion.Draw("l same")
    bethe_bloch_nsigma_kaon.Draw("l same")
    bethe_bloch_nsigma_proton.Draw("l same")


    c.SaveAs("shitty_plot.png")

    # # get electron, pion, kaon, proton nsigma 
    # nsigma_electron_bethe = inlist.FindObject("fNSigmaElectronBethe")
    # nsigma_pion_bethe = inlist.FindObject("fNSigmaPionBethe")
    # nsigma_kaon_bethe = inlist.FindObject("fNSigmaKaonBethe")
    # nsigma_proton_bethe = inlist.FindObject("fNSigmaProtonBethe")



    nsigma_electron = big_dist.Projection(2)

    # nsigma_pion_bethe.GetXaxis().SetRangeUser(2.5, 3.0 - 0.0001)
    # nsigma_electron_bethe_1d = nsigma_pion_bethe.ProjectionY()
    # nsigma_electron_bethe_1d.SetLineColor(rt.kRed)

    # nsigma_electron.Draw()
    # nsigma_electron_bethe_1d.Draw("same")
    # c.SaveAs("nsigma_electron.png")


if __name__ == "__main__":
    main()
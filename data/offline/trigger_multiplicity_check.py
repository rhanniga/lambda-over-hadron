import ROOT as rt
rt.gStyle.SetOptStat(0)

from sys import argv


def main(input_file_path, trig_pt_min, trig_pt_max):

    input_file = rt.TFile(input_file_path, 'READ')

    trig_pt_dist_0_20 = input_file.Get('trig_pt_dist_0_20')
    trig_pt_dist_20_50 = input_file.Get('trig_pt_dist_20_50')
    trig_pt_dist_50_80 = input_file.Get('trig_pt_dist_50_80')


    trig_pt_dict = {
        '0-20%': trig_pt_dist_0_20,
        '20-50%': trig_pt_dist_20_50,
        '50-80%': trig_pt_dist_50_80
    }

    colors = [rt.kGreen + 2, rt.kBlue + 2, rt.kRed + 2]

    c = rt.TCanvas('c', '', 800, 600)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.10)
    c.SetTopMargin(0.05)
    c.SetLogy()

    legend = rt.TLegend(0.6, 0.65, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    for key, trig_pt_dist in trig_pt_dict.items():

        trig_pt_dist.GetXaxis().SetRangeUser(trig_pt_min, trig_pt_max)
        trig_pt_dist.Scale(1.0 / trig_pt_dist.Integral())

        mean = trig_pt_dist.GetMean()
        mean_error = trig_pt_dist.GetMeanError()

        trig_pt_dist.SetTitle('')
        trig_pt_dist.GetXaxis().SetTitle('p_{T}^{trig} [GeV/c]')
        trig_pt_dist.GetXaxis().SetTitleOffset(1.2)
        trig_pt_dist.GetYaxis().SetTitle('Self-normalized counts')
        trig_pt_dist.SetLineColor(colors.pop())
        trig_pt_dist.SetMarkerColor(trig_pt_dist.GetLineColor())
        trig_pt_dist.SetMarkerStyle(20)
        trig_pt_dist.SetMarkerSize(0.5)

        if key == '0-20%':
            trig_pt_dist.Draw('e1')
        else:
            trig_pt_dist.Draw('e1 same')

        legend.AddEntry(trig_pt_dist, key + f", <p_{{T}}> = {mean:.2f} #pm {mean_error:.6f}", 'lp')

    legend.Draw()

    c.SaveAs('trig_pt_dist_{}_{}.pdf'.format(trig_pt_min, trig_pt_max))


if __name__ == '__main__':
    main(str(argv[1]), float(argv[2]), float(argv[3]))
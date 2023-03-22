import math

import ROOT as rt
import array as arr

from strangehelper import get_parabola

rt.gStyle.SetOptStat(0)


def get_ratio_error(numerator, denominator, numerator_error, denominator_error):
    return math.sqrt((numerator_error/denominator)**2 + (numerator*denominator_error/denominator**2)**2)


c = rt.TCanvas("c", "c", 800, 600)
c.SetRightMargin(0.05)
c.SetLeftMargin(0.15)
c.SetBottomMargin(0.17)
c.SetTopMargin(0.07)
c.SetTicks(1, 1)


lambda_phi_ratio = False
use_uncor_syst = True


low_pt_slopes = []
high_pt_slopes = []

colors = [rt.kBlue + 2, rt.kMagenta + 2, rt.kGreen + 3]

for PT_MODE in [0, 1, 2]:
    for use_new_x_axis in [True, False]:
        for lambda_phi_ratio in [True, False]:
            if PT_MODE == 0:
                graph_infile = rt.TFile("output/systematics/final_yield_ratio_syst_newnch.root")
            elif PT_MODE == 1:
                graph_infile = rt.TFile("output/systematics/final_yield_ratio_syst_newnch_lowpt.root")
            elif PT_MODE == 2:
                graph_infile = rt.TFile("output/systematics/final_yield_ratio_syst_newnch_highpt.root")

            if use_new_x_axis:

                near_graph = graph_infile.Get("near_yield_graph_new_x_axis")
                away_graph = graph_infile.Get("away_yield_graph_new_x_axis")
                ue_graph = graph_infile.Get("ue_yield_graph_new_x_axis")
                total_graph = graph_infile.Get("total_yield_graph_new_x_axis")

                hh_near_graph = graph_infile.Get("hh_near_yield_graph_new_x_axis")
                hh_away_graph = graph_infile.Get("hh_away_yield_graph_new_x_axis")
                hh_ue_graph = graph_infile.Get("hh_ue_yield_graph_new_x_axis")
                hh_total_graph = graph_infile.Get("hh_total_yield_graph_new_x_axis")

                near_graph_final_nch_dep_syst = graph_infile.Get("near_yield_graph_final_nch_dep_syst_new_x_axis")
                away_graph_final_nch_dep_syst = graph_infile.Get("away_yield_graph_final_nch_dep_syst_new_x_axis")
                ue_graph_final_nch_dep_syst = graph_infile.Get("ue_yield_graph_final_nch_dep_syst_new_x_axis")
                total_graph_final_nch_dep_syst = graph_infile.Get("total_yield_graph_final_nch_dep_syst_new_x_axis")
                hh_near_graph_final_nch_dep_syst = graph_infile.Get("hh_near_yield_graph_final_nch_dep_syst_new_x_axis")
                hh_away_graph_final_nch_dep_syst = graph_infile.Get("hh_away_yield_graph_final_nch_dep_syst_new_x_axis")
                hh_ue_graph_final_nch_dep_syst = graph_infile.Get("hh_ue_yield_graph_final_nch_dep_syst_new_x_axis")
                hh_total_graph_final_nch_dep_syst = graph_infile.Get("hh_total_yield_graph_final_nch_dep_syst_new_x_axis")

                near_graph_final_syst = graph_infile.Get("near_yield_graph_final_syst_new_x_axis")
                away_graph_final_syst = graph_infile.Get("away_yield_graph_final_syst_new_x_axis")
                ue_graph_final_syst = graph_infile.Get("ue_yield_graph_final_syst_new_x_axis")
                total_graph_final_syst = graph_infile.Get("total_yield_graph_final_syst_new_x_axis")
                hh_near_graph_final_syst = graph_infile.Get("hh_near_yield_graph_final_syst_new_x_axis")
                hh_away_graph_final_syst = graph_infile.Get("hh_away_yield_graph_final_syst_new_x_axis")
                hh_ue_graph_final_syst = graph_infile.Get("hh_ue_yield_graph_final_syst_new_x_axis")
                hh_total_graph_final_syst = graph_infile.Get("hh_total_yield_graph_final_syst_new_x_axis")


                if lambda_phi_ratio:
                    near_ratio_graph = graph_infile.Get("lambda_phi_near_ratio_graph_new_x_axis")
                    away_ratio_graph = graph_infile.Get("lambda_phi_away_ratio_graph_new_x_axis")
                    ue_ratio_graph = graph_infile.Get("lambda_phi_ue_ratio_graph_new_x_axis")
                    total_ratio_graph = graph_infile.Get("lambda_phi_total_ratio_graph_new_x_axis")

                    near_ratio_graph_final_syst = graph_infile.Get("lambda_phi_near_ratio_graph_final_syst_new_x_axis")
                    away_ratio_graph_final_syst = graph_infile.Get("lambda_phi_away_ratio_graph_final_syst_new_x_axis")
                    ue_ratio_graph_final_syst = graph_infile.Get("lambda_phi_ue_ratio_graph_final_syst_new_x_axis")
                    total_ratio_graph_final_syst = graph_infile.Get("lambda_phi_total_ratio_graph_final_syst_new_x_axis")

                    near_ratio_graph_final_nch_dep_syst = graph_infile.Get("lambda_phi_near_ratio_graph_final_syst_new_x_axis")
                    away_ratio_graph_final_nch_dep_syst = graph_infile.Get("lambda_phi_away_ratio_graph_final_syst_new_x_axis")
                    ue_ratio_graph_final_nch_dep_syst = graph_infile.Get("lambda_phi_ue_ratio_graph_final_syst_new_x_axis")
                    total_ratio_graph_final_nch_dep_syst = graph_infile.Get("lambda_phi_total_ratio_graph_final_syst_new_x_axis")


                else:

                    near_ratio_graph = graph_infile.Get("near_ratio_graph_new_x_axis")
                    away_ratio_graph = graph_infile.Get("away_ratio_graph_new_x_axis")
                    ue_ratio_graph = graph_infile.Get("ue_ratio_graph_new_x_axis")
                    total_ratio_graph = graph_infile.Get("total_ratio_graph_new_x_axis")

                    near_ratio_graph_final_syst = graph_infile.Get("near_ratio_graph_final_nch_dep_syst_new_x_axis")
                    away_ratio_graph_final_syst = graph_infile.Get("away_ratio_graph_final_nch_dep_syst_new_x_axis")
                    ue_ratio_graph_final_syst = graph_infile.Get("ue_ratio_graph_final_nch_dep_syst_new_x_axis")
                    total_ratio_graph_final_syst = graph_infile.Get("total_ratio_graph_final_nch_dep_syst_new_x_axis")

                    near_ratio_graph_final_nch_dep_syst = graph_infile.Get("near_ratio_graph_final_nch_dep_syst_new_x_axis")
                    away_ratio_graph_final_nch_dep_syst = graph_infile.Get("away_ratio_graph_final_nch_dep_syst_new_x_axis")
                    ue_ratio_graph_final_nch_dep_syst = graph_infile.Get("ue_ratio_graph_final_nch_dep_syst_new_x_axis")
                    total_ratio_graph_final_nch_dep_syst = graph_infile.Get("total_ratio_graph_final_nch_dep_syst_new_x_axis")

                    near_ratio_graph_final_syst = graph_infile.Get("near_ratio_graph_final_syst_new_x_axis")
                    away_ratio_graph_final_syst = graph_infile.Get("away_ratio_graph_final_syst_new_x_axis")
                    ue_ratio_graph_final_syst = graph_infile.Get("ue_ratio_graph_final_syst_new_x_axis")
                    total_ratio_graph_final_syst = graph_infile.Get("total_ratio_graph_final_syst_new_x_axis")

            else:
                near_graph = graph_infile.Get("near_yield_graph")
                away_graph = graph_infile.Get("away_yield_graph")
                ue_graph = graph_infile.Get("ue_yield_graph")
                total_graph = graph_infile.Get("total_yield_graph")

                hh_near_graph = graph_infile.Get("hh_near_yield_graph")
                hh_away_graph = graph_infile.Get("hh_away_yield_graph")
                hh_ue_graph = graph_infile.Get("hh_ue_yield_graph")
                hh_total_graph = graph_infile.Get("hh_total_yield_graph")

                near_graph_final_syst = graph_infile.Get("near_yield_graph_final_syst")
                away_graph_final_syst = graph_infile.Get("away_yield_graph_final_syst")
                ue_graph_final_syst = graph_infile.Get("ue_yield_graph_final_syst")
                total_graph_final_syst = graph_infile.Get("total_yield_graph_final_syst")

                hh_near_graph_final_syst = graph_infile.Get("hh_near_yield_graph_final_syst")
                hh_away_graph_final_syst = graph_infile.Get("hh_away_yield_graph_final_syst")
                hh_ue_graph_final_syst = graph_infile.Get("hh_ue_yield_graph_final_syst")
                hh_total_graph_final_syst = graph_infile.Get("hh_total_yield_graph_final_syst")

                if lambda_phi_ratio:
                    near_ratio_graph = graph_infile.Get("lambda_phi_near_ratio_graph")
                    away_ratio_graph = graph_infile.Get("lambda_phi_away_ratio_graph")
                    ue_ratio_graph = graph_infile.Get("lambda_phi_ue_ratio_graph")
                    total_ratio_graph = graph_infile.Get("lambda_phi_total_ratio_graph")
                    near_ratio_graph_final_syst = graph_infile.Get("lambda_phi_near_ratio_graph_final_syst")
                    away_ratio_graph_final_syst = graph_infile.Get("lambda_phi_away_ratio_graph_final_syst")
                    ue_ratio_graph_final_syst = graph_infile.Get("lambda_phi_ue_ratio_graph_final_syst")
                    total_ratio_graph_final_syst = graph_infile.Get("lambda_phi_total_ratio_graph_final_syst")

                else:
                    near_ratio_graph = graph_infile.Get("near_ratio_graph")
                    away_ratio_graph = graph_infile.Get("away_ratio_graph")
                    ue_ratio_graph = graph_infile.Get("ue_ratio_graph")
                    total_ratio_graph = graph_infile.Get("total_ratio_graph")
                    near_ratio_graph_final_syst = graph_infile.Get("near_ratio_graph_final_syst")
                    away_ratio_graph_final_syst = graph_infile.Get("away_ratio_graph_final_syst")
                    ue_ratio_graph_final_syst = graph_infile.Get("ue_ratio_graph_final_syst")
                    total_ratio_graph_final_syst = graph_infile.Get("total_ratio_graph_final_syst")
                    near_ue_ratio_graph = graph_infile.Get("near_ue_ratio_graph")
                    away_ue_ratio_graph = graph_infile.Get("away_ue_ratio_graph")
                    near_ue_ratio_graph_final_syst = graph_infile.Get("near_ue_ratio_graph_final_syst")
                    away_ue_ratio_graph_final_syst = graph_infile.Get("away_ue_ratio_graph_final_syst")


            near_ratio_graph.SetMarkerStyle(20)
            near_ratio_graph.SetMarkerSize(1.5)
            near_ratio_graph.SetMarkerColor(rt.kRed+1)
            near_ratio_graph.SetLineColor(rt.kRed+2)
            near_ratio_graph.SetLineWidth(2)

            near_ratio_graph_final_syst.SetMarkerStyle(20)
            near_ratio_graph_final_syst.SetMarkerSize(0)
            near_ratio_graph_final_syst.SetMarkerColor(rt.kRed+1)
            near_ratio_graph_final_syst.SetLineColor(rt.kRed+2)
            near_ratio_graph_final_syst.SetLineWidth(2)
            near_ratio_graph_final_syst.SetFillStyle(0)


            away_ratio_graph.SetMarkerStyle(21)
            away_ratio_graph.SetMarkerSize(1.5)
            away_ratio_graph.SetMarkerColor(rt.kBlue+1)
            away_ratio_graph.SetLineColor(rt.kBlue+2)
            away_ratio_graph.SetLineWidth(2)

            away_ratio_graph_final_syst.SetMarkerStyle(21)
            away_ratio_graph_final_syst.SetMarkerSize(0)
            away_ratio_graph_final_syst.SetMarkerColor(rt.kBlue+1)
            away_ratio_graph_final_syst.SetLineColor(rt.kBlue+2)
            away_ratio_graph_final_syst.SetLineWidth(2)
            away_ratio_graph_final_syst.SetFillStyle(0)

            ue_ratio_graph.SetMarkerStyle(34)
            ue_ratio_graph.SetMarkerSize(2)
            ue_ratio_graph.SetMarkerColor(rt.kGreen+2)
            ue_ratio_graph.SetLineColor(rt.kGreen+3)
            ue_ratio_graph.SetLineWidth(2)

            ue_ratio_graph_final_syst.SetMarkerSize(0)
            ue_ratio_graph_final_syst.SetMarkerColor(rt.kGreen+2)
            ue_ratio_graph_final_syst.SetLineColor(rt.kGreen+3)
            ue_ratio_graph_final_syst.SetLineWidth(2)
            ue_ratio_graph_final_syst.SetFillStyle(0)

            total_ratio_graph.SetMarkerStyle(47)
            total_ratio_graph.SetMarkerSize(2)
            total_ratio_graph.SetMarkerColor(rt.kMagenta+2)
            total_ratio_graph.SetLineColor(rt.kMagenta+3)
            total_ratio_graph.SetLineWidth(2)
            total_ratio_graph.SetFillColor(rt.kMagenta+1)
            total_ratio_graph.SetFillStyle(3144)

            total_ratio_graph_final_syst.SetMarkerSize(0)
            total_ratio_graph_final_syst.SetMarkerColor(rt.kMagenta+2)
            total_ratio_graph_final_syst.SetLineColor(rt.kMagenta+3)
            total_ratio_graph_final_syst.SetLineWidth(2)
            total_ratio_graph_final_syst.SetFillColor(rt.kMagenta+1)
            total_ratio_graph_final_syst.SetFillStyle(0)


            if use_new_x_axis:
                near_ratio_graph_final_nch_dep_syst.SetMarkerStyle(20)
                near_ratio_graph_final_nch_dep_syst.SetMarkerSize(0)
                near_ratio_graph_final_nch_dep_syst.SetMarkerColor(rt.kRed+1)
                near_ratio_graph_final_nch_dep_syst.SetLineColor(rt.kRed+2)
                near_ratio_graph_final_nch_dep_syst.SetLineWidth(0)
                near_ratio_graph_final_nch_dep_syst.SetFillColor(rt.kRed+1)
                near_ratio_graph_final_nch_dep_syst.SetFillStyle(3144)

                away_ratio_graph_final_nch_dep_syst.SetMarkerStyle(21)
                away_ratio_graph_final_nch_dep_syst.SetMarkerSize(0)
                away_ratio_graph_final_nch_dep_syst.SetMarkerColor(rt.kBlue+1)
                away_ratio_graph_final_nch_dep_syst.SetLineColor(rt.kBlue+2)
                away_ratio_graph_final_nch_dep_syst.SetLineWidth(0)
                away_ratio_graph_final_nch_dep_syst.SetFillColor(rt.kBlue+1)
                away_ratio_graph_final_nch_dep_syst.SetFillStyle(3144)

                ue_ratio_graph_final_nch_dep_syst.SetMarkerSize(0)
                ue_ratio_graph_final_nch_dep_syst.SetMarkerColor(rt.kGreen+2)
                ue_ratio_graph_final_nch_dep_syst.SetLineColor(rt.kGreen+3)
                ue_ratio_graph_final_nch_dep_syst.SetLineWidth(0)
                ue_ratio_graph_final_nch_dep_syst.SetFillColor(rt.kGreen+1)
                ue_ratio_graph_final_nch_dep_syst.SetFillStyle(3144)

                total_ratio_graph_final_nch_dep_syst.SetMarkerSize(0)
                total_ratio_graph_final_nch_dep_syst.SetMarkerColor(rt.kMagenta+2)
                total_ratio_graph_final_nch_dep_syst.SetLineColor(rt.kMagenta+3)
                total_ratio_graph_final_nch_dep_syst.SetLineWidth(0)
                total_ratio_graph_final_nch_dep_syst.SetFillColor(rt.kMagenta+1)
                total_ratio_graph_final_nch_dep_syst.SetFillStyle(3144)


            ratios_legend = rt.TLegend(0.58, 0.65, 0.93, 0.88)
            ratios_legend.SetMargin(0.35)
            ratios_legend.AddEntry(ue_ratio_graph, "Underlying Event", "pl")
            ratios_legend.AddEntry(away_ratio_graph, "Away-side (Jet)", "pl")
            ratios_legend.AddEntry(near_ratio_graph, "Near-side (Jet)", "pl")
            ratios_legend.AddEntry(total_ratio_graph, "Total (Jet + UE)", "pl")
            ratios_legend.SetLineWidth(0)

            if use_new_x_axis:
                plotting_hist = rt.TH1D("plotting_hist", "", 65, 0, 65)
                plotting_hist.SetMarkerStyle(20)
                plotting_hist.SetMarkerSize(1)
                plotting_hist.SetMarkerColor(rt.kRed+1)
                plotting_hist.SetLineColor(rt.kRed+2)
                plotting_hist.SetLineWidth(2)
                plotting_hist.GetXaxis().SetTitle("#LT#frac{d#it{N}_{ch}}{d#it{#eta}}#GT_{|#it{#eta}_{lab}| < 0.5}")
                plotting_hist.GetXaxis().SetTitleSize(0.05)
                plotting_hist.GetXaxis().SetLabelSize(0.045)
                plotting_hist.GetYaxis().SetLabelSize(0.045)
                plotting_hist.GetXaxis().SetTitleOffset(1.35)
                plotting_hist.GetXaxis().SetRangeUser(5, 55)
                if lambda_phi_ratio:
                    plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minus#phi}#right)")
                else:
                    plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minush}#right)")

                plotting_hist.GetYaxis().SetTitleSize(0.05)
                plotting_hist.GetYaxis().SetTitleOffset(1.3)
                plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
                plotting_hist.Draw("PE")
            else:
                mult_bin_widths = arr.array('d', [0.0, 20.0, 50.0, 80.0, 100.0])
                plotting_hist = rt.TH1D("plotting_hist", "", 4, mult_bin_widths)
                plotting_hist.SetMarkerStyle(20)
                plotting_hist.SetMarkerSize(1)
                plotting_hist.SetMarkerColor(rt.kRed+1)
                plotting_hist.SetLineColor(rt.kRed+2)
                plotting_hist.SetLineWidth(2)
                plotting_hist.GetXaxis().SetTitle("Multiplicity (%)")
                plotting_hist.GetXaxis().SetTitleSize(0.05)
                plotting_hist.GetXaxis().SetLabelSize(0.045)
                plotting_hist.GetYaxis().SetLabelSize(0.045)
                plotting_hist.GetXaxis().SetTitleOffset(1.2)
                plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
                plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
                plotting_hist.SetStats(0)

                if lambda_phi_ratio:
                    plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minus#phi}#right)")
                else:
                    plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minush}#right)")

                plotting_hist.GetYaxis().SetTitleSize(0.05)
                plotting_hist.GetYaxis().SetTitleOffset(1.3)
                plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
                plotting_hist.SetStats(0)
                plotting_hist.Draw("PE")
            
            label_x_start = 0.19
            label_y_start = 0.88
            label_text_space = 0.059
            alice_data_label = rt.TLatex()
            alice_data_label.SetNDC()
            alice_data_label.SetTextSize(0.045)
            alice_data_label.SetTextAlign(13)
            alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
            alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}")
            alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{4.0 <   #it{p}_{T,trig}    < 8.0 GeV/#it{c}}")
            if PT_MODE == 0:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{2.0 <   #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")
            elif PT_MODE == 1:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{1.5 <   #it{p}_{T,assoc} < 2.5 GeV/#it{c}}")
            elif PT_MODE == 2:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{2.5 <   #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")
            alice_data_label.DrawLatex(label_x_start, label_y_start - 4*label_text_space, "#bf{|#Delta#it{#eta}| < 1.2}")

            if lambda_phi_ratio:
                plotting_hist.GetYaxis().SetRangeUser(0.0, 30)
            else:
                plotting_hist.GetYaxis().SetRangeUser(0.0, 0.40)

            plotting_hist.SetStats(0)

            near_ratio_graph.Draw("PE SAME")
            near_ratio_graph_final_syst.Draw("E2 SAME")
            away_ratio_graph_final_syst.Draw("E2 SAME")
            away_ratio_graph.Draw("PE SAME")
            ue_ratio_graph_final_syst.Draw("E2 SAME")
            ue_ratio_graph.Draw("PE SAME")
            total_ratio_graph_final_syst.Draw("E2 SAME")
            total_ratio_graph.Draw("PE SAME")
            ratios_legend.Draw()
            c.Draw()

            if lambda_phi_ratio:
                if use_new_x_axis:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_lowpt.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_highpt.eps")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis.eps")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_lowpt.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_highpt.eps")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot.eps")
            else:
                if use_new_x_axis:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_new_x_axis_lowpt.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_new_x_axis_highpt.eps")
                    else:
                        c.SaveAs("figures/ratio_plot_new_x_axis.eps")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_lowpt.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_highpt.eps")
                    else:
                        c.SaveAs("figures/ratio_plot.eps")

            if lambda_phi_ratio:
                if use_new_x_axis:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_lowpt.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_highpt.pdf")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis.pdf")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_lowpt.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_highpt.pdf")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot.pdf")
            else:
                if use_new_x_axis:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_new_x_axis_lowpt.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_new_x_axis_highpt.pdf")
                    else:
                        c.SaveAs("figures/ratio_plot_new_x_axis.pdf")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_lowpt.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_highpt.pdf")
                    else:
                        c.SaveAs("figures/ratio_plot.pdf")



            if use_new_x_axis:
                near_fit = rt.TF1("near_fit", "pol1", 10, 55)
                away_fit = rt.TF1("away_fit", "pol1", 10, 55)
                ue_fit = rt.TF1("ue_fit", "pol1", 10, 55)
                total_fit = rt.TF1("total_fit", "pol1", 10, 55)

                near_fit.SetLineColor(rt.kRed+2)
                near_fit.SetLineWidth(2)
                near_fit.SetLineStyle(2)

                away_fit.SetLineColor(rt.kBlue+2)
                away_fit.SetLineWidth(2)
                away_fit.SetLineStyle(2)

                ue_fit.SetLineColor(rt.kGreen+2)
                ue_fit.SetLineWidth(2)
                ue_fit.SetLineStyle(2)

                total_fit.SetLineColor(rt.kMagenta+2)
                total_fit.SetLineWidth(2)
                total_fit.SetLineStyle(2)

                near_ratio_graph_final_nch_dep_syst.Fit(near_fit, "R")
                away_ratio_graph_final_nch_dep_syst.Fit(away_fit, "R")
                ue_ratio_graph_final_nch_dep_syst.Fit(ue_fit, "R")
                total_ratio_graph_final_nch_dep_syst.Fit(total_fit, "R")

                legend_addition = near_ratio_graph_final_nch_dep_syst.Clone("legend_addition1")
                legend_addition.SetFillColor(rt.kGray + 2)
                ratios_legend.AddEntry(legend_addition, "Mult. uncor. syst", "F")

                near_fit.Draw("SAME")
                away_fit.Draw("SAME")
                ue_fit.Draw("SAME")
                total_fit.Draw("SAME")

                near_ratio_graph_final_nch_dep_syst.Draw("E2 SAME")
                away_ratio_graph_final_nch_dep_syst.Draw("E2 SAME")
                ue_ratio_graph_final_nch_dep_syst.Draw("E2 SAME")
                total_ratio_graph_final_nch_dep_syst.Draw("E2 SAME")

                c.Draw()

                if lambda_phi_ratio:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_lowpt_with_fits.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_highpt_with_fits.eps")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_with_fits.eps")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_new_x_axis_lowpt_with_fits.eps")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_new_x_axis_highpt_with_fits.eps")
                    else:
                        c.SaveAs("figures/ratio_plot_new_x_axis_with_fits.eps")

                if lambda_phi_ratio:
                    if PT_MODE == 1:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_lowpt_with_fits.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_highpt_with_fits.pdf")
                    else:
                        c.SaveAs("figures/lambda_phi_ratio_plot_new_x_axis_with_fits.pdf")
                else:
                    if PT_MODE == 1:
                        c.SaveAs("figures/ratio_plot_new_x_axis_lowpt_with_fits.pdf")
                    elif PT_MODE == 2:
                        c.SaveAs("figures/ratio_plot_new_x_axis_highpt_with_fits.pdf")
                    else:
                        c.SaveAs("figures/ratio_plot_new_x_axis_with_fits.pdf")

                near_fit_slope = near_fit.GetParameter(1)
                near_fit_slope_err = near_fit.GetParError(1)

                away_fit_slope = away_fit.GetParameter(1)
                away_fit_slope_err = away_fit.GetParError(1)

                ue_fit_slope = ue_fit.GetParameter(1)
                ue_fit_slope_err = ue_fit.GetParError(1)

                total_fit_slope = total_fit.GetParameter(1)
                total_fit_slope_err = total_fit.GetParError(1)

                if use_new_x_axis and not lambda_phi_ratio:
                    if PT_MODE == 1:
                        low_pt_slopes.append([near_fit_slope, near_fit_slope_err])
                        low_pt_slopes.append([away_fit_slope, away_fit_slope_err])
                        low_pt_slopes.append([ue_fit_slope, ue_fit_slope_err])
                        low_pt_slopes.append([total_fit_slope, total_fit_slope_err])
                    elif PT_MODE == 2:
                        high_pt_slopes.append([near_fit_slope, near_fit_slope_err])
                        high_pt_slopes.append([away_fit_slope, away_fit_slope_err])
                        high_pt_slopes.append([ue_fit_slope, ue_fit_slope_err])
                        high_pt_slopes.append([total_fit_slope, total_fit_slope_err])

                ratio_slope_hist = rt.TH1F("ratio_slope_hist", "", 4, 0, 4)

                ratio_slope_hist.GetXaxis().SetBinLabel(1, "Near-side (jet)")
                ratio_slope_hist.GetXaxis().SetBinLabel(2, "Away-side (jet)")
                ratio_slope_hist.GetXaxis().SetBinLabel(3, "Underlying event")
                ratio_slope_hist.GetXaxis().SetBinLabel(4, "Total")

                ratio_slope_hist.SetBinContent(1, near_fit_slope)
                ratio_slope_hist.SetBinError(1, near_fit_slope_err)

                ratio_slope_hist.SetBinContent(2, away_fit_slope)
                ratio_slope_hist.SetBinError(2, away_fit_slope_err)

                ratio_slope_hist.SetBinContent(3, ue_fit_slope)
                ratio_slope_hist.SetBinError(3, ue_fit_slope_err)

                ratio_slope_hist.SetBinContent(4, total_fit_slope)
                ratio_slope_hist.SetBinError(4, total_fit_slope_err)

                ratio_slope_hist.SetMarkerStyle(20)
                ratio_slope_hist.SetMarkerSize(1.5)
                ratio_slope_hist.SetMarkerColor(colors[PT_MODE])
                ratio_slope_hist.SetLineColor(colors[PT_MODE])
                ratio_slope_hist.SetLineWidth(2)

                ratio_slope_hist.GetYaxis().SetTitle("Slope")
                ratio_slope_hist.GetYaxis().SetTitleOffset(1.2)
                ratio_slope_hist.GetYaxis().SetTitleSize(0.05)
                ratio_slope_hist.GetYaxis().SetLabelSize(0.05)
                ratio_slope_hist.GetYaxis().SetMaxDigits(3)

                ratio_slope_hist.GetXaxis().SetTitle("Kinematic region")
                ratio_slope_hist.GetXaxis().SetTitleOffset(1.2)

                ratio_slope_legend = rt.TLegend(0.2, 0.7, 0.5, 0.9)
                ratio_slope_legend.SetBorderSize(0)
                if PT_MODE == 0:
                    ratio_slope_legend.AddEntry(ratio_slope_hist, "Slopes for 2.0 < p_{T, assoc}^{#Lambda, h} < 4.0 GeV/c", "LEP")
                elif PT_MODE == 1:
                    ratio_slope_legend.AddEntry(ratio_slope_hist, "Slopes for 1.5 < p_{T, assoc}^{#Lambda, h} < 2.5 GeV/c", "LEP")
                elif PT_MODE == 2:
                    ratio_slope_legend.AddEntry(ratio_slope_hist, "Slopes for 2.5 < p_{T, assoc}^{#Lambda, h} < 4.0 GeV/c", "LEP")
                
                ratio_slope_hist.GetYaxis().SetRangeUser(-0.001, 0.005)

                ratio_slope_hist.Draw()
                ratio_slope_legend.Draw("SAME")

                if PT_MODE == 1:    
                    c.SaveAs("figures/ratio_slope_plot_lowpt.eps")
                elif PT_MODE == 2:
                    c.SaveAs("figures/ratio_slope_plot_highpt.eps")
                else:
                    c.SaveAs("figures/ratio_slope_plot.eps")
            
            





            hh_near_graph_final_syst.SetMarkerSize(0)
            hh_near_graph_final_syst.SetLineColor(rt.kRed + 3)
            hh_near_graph_final_syst.SetFillColor(rt.kRed + 3)
            hh_near_graph_final_syst.SetFillStyle(0)
            hh_near_graph_final_syst.SetLineWidth(2)

            hh_near_graph.SetMarkerStyle(20)
            hh_near_graph.SetMarkerSize(1)
            hh_near_graph.SetMarkerColor(rt.kRed + 3)
            hh_near_graph.SetLineColor(rt.kRed + 3)
            hh_near_graph.SetFillColor(rt.kRed + 3)
            hh_near_graph.SetFillStyle(3144)
            hh_near_graph.SetLineWidth(2)

            near_graph_final_syst.SetMarkerSize(0)
            near_graph_final_syst.SetFillColor(rt.kPink-1)
            near_graph_final_syst.SetLineColor(rt.kPink-1)
            near_graph_final_syst.SetFillStyle(0)
            near_graph_final_syst.SetLineWidth(2)

            near_graph.SetMarkerStyle(21)
            near_graph.SetMarkerSize(1)
            near_graph.SetMarkerColor(rt.kPink-1)
            near_graph.SetLineColor(rt.kPink-1)
            near_graph.SetFillColor(rt.kPink-1)
            near_graph.SetFillStyle(3144)
            near_graph.SetLineWidth(2)

            hh_away_graph_final_syst.SetMarkerSize(0)
            hh_away_graph_final_syst.SetLineWidth(2)
            hh_away_graph_final_syst.SetLineColor(rt.kViolet - 6)
            hh_away_graph_final_syst.SetFillColor(rt.kViolet - 6)
            hh_away_graph_final_syst.SetFillStyle(0)

            hh_away_graph.SetMarkerStyle(20)
            hh_away_graph.SetMarkerSize(1)
            hh_away_graph.SetMarkerColor(rt.kViolet - 6)
            hh_away_graph.SetLineColor(rt.kViolet - 6)
            hh_away_graph.SetLineWidth(2)
            hh_away_graph.SetFillColor(rt.kViolet - 6)
            hh_away_graph.SetFillStyle(3144)

            away_graph_final_syst.SetMarkerSize(0)
            away_graph_final_syst.SetLineWidth(2)
            away_graph_final_syst.SetLineColor(rt.kBlue-2)
            away_graph_final_syst.SetFillColor(rt.kBlue - 2)
            away_graph_final_syst.SetFillStyle(0)

            away_graph.SetMarkerStyle(21)
            away_graph.SetMarkerSize(1)
            away_graph.SetMarkerColor(rt.kBlue-2)
            away_graph.SetLineColor(rt.kBlue-2)
            away_graph.SetLineWidth(2)
            away_graph.SetFillColor(rt.kBlue - 2)
            away_graph.SetFillStyle(3144)

            if use_new_x_axis:
                hh_near_graph_final_nch_dep_syst.SetMarkerSize(0)
                hh_near_graph_final_nch_dep_syst.SetLineWidth(0)
                hh_near_graph_final_nch_dep_syst.SetLineColor(rt.kRed + 3)
                hh_near_graph_final_nch_dep_syst.SetFillColor(rt.kRed + 3)
                hh_near_graph_final_nch_dep_syst.SetFillStyle(3144)

                hh_away_graph_final_nch_dep_syst.SetMarkerSize(0)
                hh_away_graph_final_nch_dep_syst.SetLineWidth(0)
                hh_away_graph_final_nch_dep_syst.SetLineColor(rt.kViolet - 6)
                hh_away_graph_final_nch_dep_syst.SetFillColor(rt.kViolet - 6)
                hh_away_graph_final_nch_dep_syst.SetFillStyle(3144)

                near_graph_final_nch_dep_syst.SetMarkerSize(0)
                near_graph_final_nch_dep_syst.SetLineWidth(0)
                near_graph_final_nch_dep_syst.SetLineColor(rt.kPink-1)
                near_graph_final_nch_dep_syst.SetFillColor(rt.kPink-1)
                near_graph_final_nch_dep_syst.SetFillStyle(3144)

                away_graph_final_nch_dep_syst.SetMarkerSize(0)
                away_graph_final_nch_dep_syst.SetLineWidth(0)
                away_graph_final_nch_dep_syst.SetLineColor(rt.kBlue-2)
                away_graph_final_nch_dep_syst.SetFillColor(rt.kBlue - 2)
                away_graph_final_nch_dep_syst.SetFillStyle(3144)


            ratios_legend = rt.TLegend(0.57, 0.68, 0.79, 0.90)
            ratios_legend.SetMargin(0.3)
            ratios_legend.AddEntry(near_graph, "h#minus#Lambda near-side yield (x20)", "pl")
            ratios_legend.AddEntry(away_graph, "h#minus#Lambda away-side yield (x20)" , "pl")
            ratios_legend.AddEntry(hh_near_graph, "h#minush near-side yield", "pl")
            ratios_legend.AddEntry(hh_away_graph, "h#minush away-side yield" , "pl")


            ratios_legend.SetLineWidth(0)
            ratios_legend.SetFillStyle(0)
            ratios_legend.SetTextSize(0.037)


            if use_new_x_axis:

                plotting_hist = rt.TH1D("plotting_hist", "", 65, 0, 65)
                plotting_hist.SetMarkerStyle(20)
                plotting_hist.SetMarkerSize(1)
                plotting_hist.SetMarkerColor(rt.kRed+1)
                plotting_hist.SetLineColor(rt.kRed+2)
                plotting_hist.SetLineWidth(2)
                plotting_hist.GetXaxis().SetTitle("#LT#frac{d#it{N}_{ch}}{d#it{#eta}}#GT_{|#it{#eta}_{lab}| < 0.5}")
                plotting_hist.GetXaxis().SetTitleSize(0.05)
                plotting_hist.GetXaxis().SetLabelSize(0.045)
                plotting_hist.GetYaxis().SetLabelSize(0.045)
                plotting_hist.GetXaxis().SetTitleOffset(1.35)
                plotting_hist.GetXaxis().SetRangeUser(5, 55)

                if lambda_phi_ratio:
                    plotting_hist.GetYaxis().SetTitle("Per-trigger pairwise yield")
                else:
                    plotting_hist.GetYaxis().SetTitle("Per-trigger pairwise yield")

                plotting_hist.GetYaxis().SetTitleSize(0.05)
                plotting_hist.GetYaxis().SetTitleOffset(1.1)
                plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)

                plotting_hist.Draw("PE")
            else:

                mult_bin_widths = arr.array('d', [0.0, 20.0, 50.0, 80.0, 100.0])
                plotting_hist = rt.TH1D("plotting_hist", "", 4, mult_bin_widths)
                plotting_hist.SetMarkerStyle(20)
                plotting_hist.SetMarkerSize(1)
                plotting_hist.SetMarkerColor(rt.kRed+1)
                plotting_hist.SetLineColor(rt.kRed+2)
                plotting_hist.SetLineWidth(2)
                plotting_hist.GetXaxis().SetTitle("Multiplicity (%)")
                plotting_hist.GetXaxis().SetTitleSize(0.05)
                plotting_hist.GetXaxis().SetLabelSize(0.045)
                plotting_hist.GetYaxis().SetLabelSize(0.045)
                plotting_hist.GetXaxis().SetTitleOffset(1.2)
                plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
                plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
                plotting_hist.SetStats(0)

                if lambda_phi_ratio:
                    plotting_hist.GetYaxis().SetTitle("Per-trigger pairwise yield")
                else:
                    plotting_hist.GetYaxis().SetTitle("Per-trigger pairwise yield")
                plotting_hist.GetYaxis().SetTitleSize(0.05)
                plotting_hist.GetYaxis().SetTitleOffset(1.1)
                plotting_hist.SetStats(0)
                plotting_hist.Draw("PE")


            label_x_start = 0.19
            label_y_start = 0.88
            label_text_space = 0.059
            alice_data_label = rt.TLatex()
            alice_data_label.SetNDC()
            alice_data_label.SetTextSize(0.045)
            alice_data_label.SetTextAlign(13)
            alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
            alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}")
            alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{4.0 <   #it{p}_{T,trig}    < 8.0 GeV/#it{c}}")
            if PT_MODE == 0:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{2.0 <   #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")
            elif PT_MODE == 1:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{1.5 <   #it{p}_{T,assoc} < 2.5 GeV/#it{c}}")
            elif PT_MODE == 2:
                alice_data_label.DrawLatex(label_x_start, label_y_start - 3*label_text_space, "#bf{2.5 <   #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")
            alice_data_label.DrawLatex(label_x_start, label_y_start - 4*label_text_space, "#bf{|#Delta#it{#eta}| < 1.2}")

            # scale the h-lambda graph points
            SCALE_FACTOR = 20
            for i in range(3):

                    near_graph_final_syst.SetPoint(i, near_graph_final_syst.GetX()[i], near_graph_final_syst.GetY()[i] * SCALE_FACTOR)
                    near_graph_final_nch_dep_syst.SetPoint(i, near_graph_final_nch_dep_syst.GetX()[i], near_graph_final_nch_dep_syst.GetY()[i] * SCALE_FACTOR)
                    near_graph.SetPoint(i, near_graph.GetX()[i], near_graph.GetY()[i] * SCALE_FACTOR)
                    near_graph_final_syst.SetPointError(i, near_graph_final_syst.GetEXlow()[i], near_graph_final_syst.GetEXhigh()[i], near_graph_final_syst.GetEYlow()[i] * SCALE_FACTOR, near_graph_final_syst.GetEYhigh()[i] * SCALE_FACTOR)
                    near_graph_final_nch_dep_syst.SetPointError(i, near_graph_final_nch_dep_syst.GetEXlow()[i], near_graph_final_nch_dep_syst.GetEXhigh()[i], near_graph_final_nch_dep_syst.GetEYlow()[i] * SCALE_FACTOR, near_graph_final_nch_dep_syst.GetEYhigh()[i] * SCALE_FACTOR)
                    near_graph.SetPointError(i, near_graph.GetErrorX(i), near_graph.GetErrorY(i) * SCALE_FACTOR)

                    away_graph_final_syst.SetPoint(i, away_graph_final_syst.GetX()[i], away_graph_final_syst.GetY()[i] * SCALE_FACTOR)
                    away_graph_final_nch_dep_syst.SetPoint(i, away_graph_final_nch_dep_syst.GetX()[i], away_graph_final_nch_dep_syst.GetY()[i] * SCALE_FACTOR)
                    away_graph.SetPoint(i, away_graph.GetX()[i], away_graph.GetY()[i] * SCALE_FACTOR)
                    away_graph_final_syst.SetPointError(i, away_graph_final_syst.GetEXlow()[i], away_graph_final_syst.GetEXhigh()[i], away_graph_final_syst.GetEYlow()[i] * SCALE_FACTOR, away_graph_final_syst.GetEYhigh()[i] * SCALE_FACTOR)
                    away_graph_final_nch_dep_syst.SetPointError(i, away_graph_final_nch_dep_syst.GetEXlow()[i], away_graph_final_nch_dep_syst.GetEXhigh()[i], away_graph_final_nch_dep_syst.GetEYlow()[i] * SCALE_FACTOR, away_graph_final_nch_dep_syst.GetEYhigh()[i] * SCALE_FACTOR)
                    away_graph.SetPointError(i, away_graph.GetErrorX(i), away_graph.GetErrorY(i) * SCALE_FACTOR)


            hh_near_graph_final_syst.Draw("E2")
            hh_near_graph.Draw("PE SAME")
            hh_away_graph_final_syst.Draw("E2")
            hh_away_graph.Draw("PE SAME")


            near_graph_final_syst.Draw("E2 SAME")
            near_graph.Draw("PE SAME")
            away_graph_final_syst.Draw("E2 SAME")
            away_graph.Draw("PE SAME")
            ratios_legend.Draw()
            c.Draw()

            if use_new_x_axis:
                    if PT_MODE == 2:
                            c.SaveAs("figures/pairwise_plot_new_x_axis_highpt.eps")
                    elif PT_MODE == 1:
                            c.SaveAs("figures/pairwise_plot_new_x_axis_lowpt.eps")
                    else:
                            c.SaveAs("figures/pairwise_plot_new_x_axis.eps")
            else:
                    if PT_MODE == 2:
                            c.SaveAs("figures/pairwise_plot_highpt.eps")
                    elif PT_MODE == 1:
                            c.SaveAs("figures/pairwise_plot_lowpt.eps")
                    else:
                            c.SaveAs("figures/pairwise_plot.eps")

            if use_new_x_axis:
                    if PT_MODE == 2:
                            c.SaveAs("figures/pairwise_plot_new_x_axis_highpt.pdf")
                    elif PT_MODE == 1:
                            c.SaveAs("figures/pairwise_plot_new_x_axis_lowpt.pdf")
                    else:
                            c.SaveAs("figures/pairwise_plot_new_x_axis.pdf")
            else:
                    if PT_MODE == 2:
                            c.SaveAs("figures/pairwise_plot_highpt.pdf")
                    elif PT_MODE == 1:
                            c.SaveAs("figures/pairwise_plot_lowpt.pdf")
                    else:
                            c.SaveAs("figures/pairwise_plot.pdf")

            if use_new_x_axis:
                near_fit = rt.TF1("near_fit", "pol1", 10, 55)
                away_fit = rt.TF1("away_fit", "pol1", 10, 55)

                hh_near_fit = rt.TF1("hh_near_fit", "pol1", 10, 55)
                hh_away_fit = rt.TF1("hh_away_fit", "pol1", 10, 55)

                near_fit.SetLineColor(rt.kPink - 1)
                near_fit.SetLineWidth(2)
                near_fit.SetLineStyle(2)
                away_fit.SetLineColor(rt.kBlue-2)
                away_fit.SetLineWidth(2)
                away_fit.SetLineStyle(2)

                hh_near_fit.SetLineColor(rt.kRed + 3)
                hh_near_fit.SetLineWidth(2)
                hh_near_fit.SetLineStyle(2)
                hh_away_fit.SetLineColor(rt.kViolet - 6)
                hh_away_fit.SetLineWidth(2)
                hh_away_fit.SetLineStyle(2)



                near_graph_final_nch_dep_syst.Fit(near_fit, "R")
                away_graph_final_nch_dep_syst.Fit(away_fit, "R")
                hh_near_graph_final_nch_dep_syst.Fit(hh_near_fit, "R")
                hh_away_graph_final_nch_dep_syst.Fit(hh_away_fit, "R")

                legend_addition = near_graph_final_nch_dep_syst.Clone("legend_addition2")
                legend_addition.SetFillColor(rt.kGray + 2)
                legend_addition.SetFillStyle(3144)
                ratios_legend.AddEntry(legend_addition, "Mult. uncor. syst", "F")

                near_fit.Draw("SAME")
                away_fit.Draw("SAME")
                hh_near_fit.Draw("SAME")
                hh_away_fit.Draw("SAME")

                near_graph_final_nch_dep_syst.Draw("E2 SAME")
                away_graph_final_nch_dep_syst.Draw("E2 SAME")
                hh_near_graph_final_nch_dep_syst.Draw("E2 SAME")
                hh_away_graph_final_nch_dep_syst.Draw("E2 SAME")

                c.Draw()
                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_lowpt_with_fits.eps")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_highpt_with_fits.eps")
                else:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_with_fits.eps")

                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_lowpt_with_fits.pdf")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_highpt_with_fits.pdf")
                else:
                    c.SaveAs("figures/pairwise_plot_new_x_axis_with_fits.pdf")

                near_fit_slope = near_fit.GetParameter(1)/SCALE_FACTOR
                near_fit_slope_err = near_fit.GetParError(1)/SCALE_FACTOR
                away_fit_slope = away_fit.GetParameter(1)/SCALE_FACTOR
                away_fit_slope_err = away_fit.GetParError(1)/SCALE_FACTOR
                hh_near_fit_slope = hh_near_fit.GetParameter(1)
                hh_near_fit_slope_err = hh_near_fit.GetParError(1)
                hh_away_fit_slope = hh_away_fit.GetParameter(1)
                hh_away_fit_slope_err = hh_away_fit.GetParError(1)
                h_lambda_yield_slope_hist = rt.TH1F("h_lambda_yield_slope_hist", "", 2, 0, 2)
                h_h_yield_slope_hist = rt.TH1F("h_h_yield_slope_hist", "", 2, 0, 2)
                h_lambda_yield_slope_hist.SetBinContent(1, near_fit_slope)
                h_lambda_yield_slope_hist.SetBinError(1, near_fit_slope_err)
                h_lambda_yield_slope_hist.SetBinContent(2, away_fit_slope)
                h_lambda_yield_slope_hist.SetBinError(2, away_fit_slope_err)
                h_h_yield_slope_hist.SetBinContent(1, hh_near_fit_slope)
                h_h_yield_slope_hist.SetBinError(1, hh_near_fit_slope_err)
                h_h_yield_slope_hist.SetBinContent(2, hh_away_fit_slope)
                h_h_yield_slope_hist.SetBinError(2, hh_away_fit_slope_err)
                h_lambda_yield_slope_hist.SetMarkerStyle(20)
                h_lambda_yield_slope_hist.SetMarkerSize(1.5)
                h_lambda_yield_slope_hist.SetMarkerColor(colors[PT_MODE])
                h_lambda_yield_slope_hist.SetLineColor(colors[PT_MODE])
                h_lambda_yield_slope_hist.SetLineWidth(2)
                h_h_yield_slope_hist.SetMarkerStyle(21)
                h_h_yield_slope_hist.SetMarkerSize(1.5)
                h_h_yield_slope_hist.SetMarkerColor(colors[PT_MODE] - 6)
                h_h_yield_slope_hist.SetLineColor(colors[PT_MODE] - 6)
                h_h_yield_slope_hist.SetLineWidth(2)
                h_lambda_yield_slope_hist.GetXaxis().SetBinLabel(1, "Near-side (jet)")
                h_lambda_yield_slope_hist.GetXaxis().SetBinLabel(2, "Away-side (jet)")
                h_lambda_yield_slope_hist.GetXaxis().SetTitle("Kinematic region")
                h_lambda_yield_slope_hist.GetYaxis().SetTitle("Slope")
                h_lambda_yield_slope_hist.GetYaxis().SetTitleOffset(1.2)
                h_lambda_yield_slope_hist.GetYaxis().SetMaxDigits(3)
                h_lambda_yield_slope_hist.GetYaxis().SetRangeUser(-0.001, 0.005)
                h_lambda_yield_slope_hist.GetXaxis().SetLabelSize(0.05)
                h_lambda_yield_slope_hist.GetXaxis().SetTitleOffset(1.2)
                leg = rt.TLegend(0.2, 0.7, 0.5, 0.9)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                if PT_MODE == 1:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda yield slopes (low p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h yield slopes (low p_{T} bin)", "PL")
                elif PT_MODE == 2:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda yield slopes (high p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h yield slopes (high p_{T} bin)", "PL")
                elif PT_MODE == 0:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda yield slopes (central p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h yield slopes (central p_{T} bin)", "PL")
                h_lambda_yield_slope_hist.Draw()
                h_h_yield_slope_hist.Draw("SAME")
                leg.Draw("SAME")
                c.Draw()

                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_slopes_lowpt.eps")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_slopes_highpt.eps")
                else:
                    c.SaveAs("figures/pairwise_plot_slopes.eps")

                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_slopes_lowpt.pdf")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_slopes_highpt.pdf")
                else:
                    c.SaveAs("figures/pairwise_plot_slopes.pdf")



                print(f"near slope: {near_fit_slope:.2e} $\\pm$ {near_fit_slope_err:.2e}")
                print(f"away slope: {away_fit_slope:.2e} $\\pm$ {away_fit_slope_err:.2e}")

                print(f"hh near slope: {hh_near_fit_slope:.2e} $\\pm$ {hh_near_fit_slope_err:.2e}")
                print(f"hh away slope: {hh_away_fit_slope:.2e} $\\pm$ {hh_away_fit_slope_err:.2e}")



                # For the yield plots, the better metric is percent increase

                near_fit_slope = ((near_graph_final_nch_dep_syst.GetY()[2]/near_graph_final_nch_dep_syst.GetY()[0]))
                near_fit_slope_err = get_ratio_error(near_graph_final_nch_dep_syst.GetY()[2], near_graph_final_nch_dep_syst.GetY()[0], near_graph_final_nch_dep_syst.GetEYhigh()[2], near_graph_final_nch_dep_syst.GetEYhigh()[0])
                away_fit_slope = ((away_graph_final_nch_dep_syst.GetY()[2]/away_graph_final_nch_dep_syst.GetY()[0]))
                away_fit_slope_err = get_ratio_error(away_graph_final_nch_dep_syst.GetY()[2], away_graph_final_nch_dep_syst.GetY()[0], away_graph_final_nch_dep_syst.GetEYhigh()[2], away_graph_final_nch_dep_syst.GetEYhigh()[0])

                hh_near_fit_slope = ((hh_near_graph_final_nch_dep_syst.GetY()[2]/hh_near_graph_final_nch_dep_syst.GetY()[0]))
                hh_near_fit_slope_err = get_ratio_error(hh_near_graph_final_nch_dep_syst.GetY()[2], hh_near_graph_final_nch_dep_syst.GetY()[0], hh_near_graph_final_nch_dep_syst.GetEYhigh()[2], hh_near_graph_final_nch_dep_syst.GetEYhigh()[0])
                hh_away_fit_slope = ((hh_away_graph_final_nch_dep_syst.GetY()[2]/hh_away_graph_final_nch_dep_syst.GetY()[0]))
                hh_away_fit_slope_err = get_ratio_error(hh_away_graph_final_nch_dep_syst.GetY()[2], hh_away_graph_final_nch_dep_syst.GetY()[0], hh_away_graph_final_nch_dep_syst.GetEYhigh()[2], hh_away_graph_final_nch_dep_syst.GetEYhigh()[0])

                h_lambda_yield_slope_hist = rt.TH1F("h_lambda_yield_slope_hist", "", 2, 0, 2)
                h_h_yield_slope_hist = rt.TH1F("h_h_yield_slope_hist", "", 2, 0, 2)
                h_lambda_yield_slope_hist.SetBinContent(1, near_fit_slope)
                h_lambda_yield_slope_hist.SetBinError(1, near_fit_slope_err)
                h_lambda_yield_slope_hist.SetBinContent(2, away_fit_slope)
                h_lambda_yield_slope_hist.SetBinError(2, away_fit_slope_err)
                h_h_yield_slope_hist.SetBinContent(1, hh_near_fit_slope)
                h_h_yield_slope_hist.SetBinError(1, hh_near_fit_slope_err)
                h_h_yield_slope_hist.SetBinContent(2, hh_away_fit_slope)
                h_h_yield_slope_hist.SetBinError(2, hh_away_fit_slope_err)
                h_lambda_yield_slope_hist.SetMarkerStyle(20)
                h_lambda_yield_slope_hist.SetMarkerSize(1.5)
                h_lambda_yield_slope_hist.SetMarkerColor(colors[PT_MODE])
                h_lambda_yield_slope_hist.SetLineColor(colors[PT_MODE])
                h_lambda_yield_slope_hist.SetLineWidth(2)
                h_h_yield_slope_hist.SetMarkerStyle(21)
                h_h_yield_slope_hist.SetMarkerSize(1.5)
                h_h_yield_slope_hist.SetMarkerColor(colors[PT_MODE] - 6)
                h_h_yield_slope_hist.SetLineColor(colors[PT_MODE] - 6)
                h_h_yield_slope_hist.SetLineWidth(2)
                h_lambda_yield_slope_hist.GetXaxis().SetBinLabel(1, "Near-side (jet)")
                h_lambda_yield_slope_hist.GetXaxis().SetBinLabel(2, "Away-side (jet)")
                h_lambda_yield_slope_hist.GetXaxis().SetTitle("Kinematic region")
                h_lambda_yield_slope_hist.GetYaxis().SetTitle("Yield_{0-20%}/Yield_{50-80%}")
                h_lambda_yield_slope_hist.GetYaxis().SetTitleOffset(1.2)
                h_lambda_yield_slope_hist.GetYaxis().SetMaxDigits(3)
                h_lambda_yield_slope_hist.GetYaxis().SetRangeUser(0.8, 3)
                h_lambda_yield_slope_hist.GetXaxis().SetLabelSize(0.05)
                h_lambda_yield_slope_hist.GetXaxis().SetTitleOffset(1.2)
                leg = rt.TLegend(0.2, 0.7, 0.5, 0.9)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                if PT_MODE == 1:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda (low p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h (low p_{T} bin)", "PL")
                elif PT_MODE == 2:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda (high p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h (high p_{T} bin)", "PL")
                elif PT_MODE == 0:
                    leg.AddEntry(h_lambda_yield_slope_hist, "h-#Lambda (central p_{T} bin)", "PL")
                    leg.AddEntry(h_h_yield_slope_hist, "h-h (central p_{T} bin)", "PL")
                h_lambda_yield_slope_hist.Draw()
                h_h_yield_slope_hist.Draw("SAME")
                leg.Draw("SAME")
                c.Draw()

                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_increase_lowpt.eps")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_increase_highpt.eps")
                else:
                    c.SaveAs("figures/pairwise_plot_increase.eps")

                if PT_MODE == 1:
                    c.SaveAs("figures/pairwise_plot_increase_lowpt.pdf")
                elif PT_MODE == 2:
                    c.SaveAs("figures/pairwise_plot_increase_highpt.pdf")
                else:
                    c.SaveAs("figures/pairwise_plot_increase.pdf")




slope_plotting_hist = rt.TH1F("slope_plotting_hist", "", 7, 1.0, 4.5)

slope_plotting_hist.SetLineWidth(0)
slope_plotting_hist.SetMarkerSize(0)
slope_plotting_hist.GetXaxis().SetTitle("#it{p}_{T}^{assoc.} (GeV/#it{c})")
# slope_plotting_hist.GetYaxis().SetTitle("Yield Ratio #left(#frac{h#minus#Lambda}{h#minush}#right) vs. #LT#frac{d#it{N}_{ch}}{d#it{#eta}}#GT_{|#it{#eta}_{lab}| < 0.5} slope")
slope_plotting_hist.GetYaxis().SetTitle("Yield ratio slope")
slope_plotting_hist.GetXaxis().SetTitleSize(0.05)
slope_plotting_hist.GetXaxis().SetLabelSize(0.045)
slope_plotting_hist.GetYaxis().SetLabelSize(0.045)
slope_plotting_hist.GetXaxis().SetTitleOffset(1.2)
slope_plotting_hist.GetXaxis().SetRangeUser(0.0, 100.0)
slope_plotting_hist.GetYaxis().SetRangeUser(0.0, 1.5)
slope_plotting_hist.SetStats(0)
slope_plotting_hist.GetYaxis().SetTitleSize(0.05)
slope_plotting_hist.GetYaxis().SetTitleOffset(0.9)
slope_plotting_hist.GetYaxis().SetRangeUser(0.0, 6e-3)
slope_plotting_hist.GetYaxis().SetMaxDigits(3)
slope_plotting_hist.GetXaxis().SetRangeUser(1.01, 4.49)
slope_plotting_hist.SetStats(0)
slope_plotting_hist.Draw()


pt_bins = arr.array('d', [1.5, 2.5, 4.0])
near_slope_hist = rt.TH1F("near_slope_hist", "", 2, pt_bins)
away_slope_hist = rt.TH1F("away_slope_hist", "", 2, pt_bins)
ue_slope_hist = rt.TH1F("ue_slope_hist", "", 2, pt_bins)
total_slope_hist = rt.TH1F("total_slope_hist", "", 2, pt_bins)

near_slope_hist.SetMarkerStyle(20)
near_slope_hist.SetMarkerSize(1.5)
near_slope_hist.SetMarkerColor(rt.kRed+1)
near_slope_hist.SetLineColor(rt.kRed+2)
near_slope_hist.SetLineWidth(2)
near_slope_hist.SetFillColor(rt.kRed+1)
near_slope_hist.SetFillStyle(3144)

away_slope_hist.SetMarkerStyle(21)
away_slope_hist.SetMarkerSize(1.5)
away_slope_hist.SetMarkerColor(rt.kBlue+1)
away_slope_hist.SetLineColor(rt.kBlue+2)
away_slope_hist.SetLineWidth(2)
away_slope_hist.SetFillColor(rt.kBlue+1)
away_slope_hist.SetFillStyle(3144)

ue_slope_hist.SetMarkerStyle(34)
ue_slope_hist.SetMarkerSize(1.5)
ue_slope_hist.SetMarkerColor(rt.kGreen+2)
ue_slope_hist.SetLineColor(rt.kGreen+3)
ue_slope_hist.SetLineWidth(2)
ue_slope_hist.SetFillColor(rt.kGreen+1)
ue_slope_hist.SetFillStyle(3144)

total_slope_hist.SetMarkerStyle(47)
total_slope_hist.SetMarkerSize(1.5)
total_slope_hist.SetMarkerColor(rt.kMagenta+2)
total_slope_hist.SetLineColor(rt.kMagenta+3)
total_slope_hist.SetLineWidth(2)
total_slope_hist.SetFillColor(rt.kMagenta+1)
total_slope_hist.SetFillStyle(3144)


near_slope_hist.SetBinContent(1, low_pt_slopes[0][0])
near_slope_hist.SetBinError(1, low_pt_slopes[0][1])
print(low_pt_slopes)
near_slope_hist.SetBinContent(2, high_pt_slopes[0][0])
near_slope_hist.SetBinError(2, high_pt_slopes[0][1])

away_slope_hist.SetBinContent(1, low_pt_slopes[1][0])
away_slope_hist.SetBinError(1, low_pt_slopes[1][1])
away_slope_hist.SetBinContent(2, high_pt_slopes[1][0])
away_slope_hist.SetBinError(2, high_pt_slopes[1][1])

ue_slope_hist.SetBinContent(1, low_pt_slopes[2][0])
ue_slope_hist.SetBinError(1, low_pt_slopes[2][1])
ue_slope_hist.SetBinContent(2, high_pt_slopes[2][0])
ue_slope_hist.SetBinError(2, high_pt_slopes[2][1])

total_slope_hist.SetBinContent(1, low_pt_slopes[3][0])
total_slope_hist.SetBinError(1, low_pt_slopes[3][1])
total_slope_hist.SetBinContent(2, high_pt_slopes[3][0])
total_slope_hist.SetBinError(2, high_pt_slopes[3][1])


slopes_legend = rt.TLegend(0.58, 0.65, 0.93, 0.88)
slopes_legend.SetMargin(0.35)
slopes_legend.AddEntry(ue_slope_hist, "Underlying Event", "pl")
slopes_legend.AddEntry(away_slope_hist, "Away-side (Jet)", "pl")
slopes_legend.AddEntry(near_slope_hist, "Near-side (Jet)", "pl")
slopes_legend.AddEntry(total_slope_hist, "Total (Jet + UE)", "pl")
slopes_legend.SetLineWidth(0)

label_x_start = 0.19
label_y_start = 0.88
label_text_space = 0.059
alice_data_label = rt.TLatex()
alice_data_label.SetNDC()
alice_data_label.SetTextSize(0.045)
alice_data_label.SetTextAlign(13)
alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Preliminary")
alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}")
alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{4.0 <   #it{p}_{T}^{trig.}    < 8.0 GeV/#it{c}}")
alice_data_label.DrawLatex(label_x_start, label_y_start - 3.2*label_text_space, "#bf{|#Delta#it{#eta}| < 1.2}")

slopes_legend.Draw("SAME")
total_slope_hist.Draw("SAME E1")
near_slope_hist.Draw("SAME E1")
away_slope_hist.Draw("SAME E1")
ue_slope_hist.Draw("SAME E1")

c.Draw()
c.SaveAs("figures/consolidated_slope_plot.pdf")
c.SaveAs("figures/consolidated_slope_plot.eps")
import ROOT as rt
rt.gStyle.SetOptStat(0)


c = rt.TCanvas("c", "c", 800, 600)

for PT_MODE in [0, 1, 2]:
    for MIX_COR_MODE in [0, 1, 2]:
        for FANCY_LABEL in [True, False]:
            if PT_MODE == 0:
                infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_20_40_delta_eta_12_normal.root")
            elif PT_MODE == 1:
                infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_15_25_delta_eta_12_normal_lowpt.root")
            elif PT_MODE == 2:
                infile = rt.TFile("output/central_value/v0_avg6_sideband_subtraction_rsb_1135_115_sig_1102_113_trig_40_80_assoc_25_40_delta_eta_12_normal_highpt.root")

            mult_strings = ["_0_20", "_20_50", "_50_80"]

            for mult_string in mult_strings:

                if MIX_COR_MODE == 0:
                    h_lambda_2d_dist = infile.Get("h_lambda_2d_mixcor_sig" + mult_string)
                    h_h_2d_dist = infile.Get("h_h_2d_mixcor" + mult_string)
                elif MIX_COR_MODE == 1:
                    h_lambda_2d_dist = infile.Get("h_lambda_2d_nomixcor" + mult_string)
                    h_h_2d_dist = infile.Get("h_h_2d_nomixcor" + mult_string)
                elif MIX_COR_MODE == 2:
                    h_lambda_2d_dist = infile.Get("h_lambda_mixed_2d" + mult_string)
                    h_h_2d_dist = infile.Get("h_h_mixed_2d" + mult_string)

                h_lambda_2d_dist.SetTitle("")
                h_lambda_2d_dist.GetZaxis().SetMaxDigits(3)
                if MIX_COR_MODE == 0:
                    h_lambda_2d_dist.GetZaxis().SetTitle("#frac{1}{N_{trig}}(h-#Lambda pairs per (#Delta#eta, #Delta#varphi) bin)")
                else:
                    h_lambda_2d_dist.GetZaxis().SetTitle("h-#Lambda pairs per (#Delta#eta, #Delta#varphi) bin")
                h_lambda_2d_dist.GetZaxis().SetTitleOffset(1.4)
                h_lambda_2d_dist.GetZaxis().SetLabelOffset(0.01)
                h_lambda_2d_dist.GetYaxis().SetTitle("#Delta#varphi_{h-#Lambda}")
                h_lambda_2d_dist.GetXaxis().SetTitle("#Delta#eta_{h-#Lambda}")
                h_lambda_2d_dist.GetXaxis().SetTitleSize(0.05)
                h_lambda_2d_dist.GetXaxis().SetTitleOffset(1.3)
                h_lambda_2d_dist.GetYaxis().SetTitleSize(0.05)
                h_lambda_2d_dist.GetYaxis().SetTitleOffset(1.3)

                h_h_2d_dist.SetTitle("")
                h_h_2d_dist.GetZaxis().SetMaxDigits(3)
                if MIX_COR_MODE == 0:
                    h_h_2d_dist.GetZaxis().SetTitle("#frac{1}{N_{trig}}(h-h pairs per (#Delta#eta, #Delta#varphi) bin)")
                else:
                    h_h_2d_dist.GetZaxis().SetTitle("h-h pairs per (#Delta#eta, #Delta#varphi) bin")
                h_h_2d_dist.GetZaxis().SetTitleOffset(1.3)
                h_h_2d_dist.GetZaxis().SetLabelOffset(0.01)
                h_h_2d_dist.GetYaxis().SetTitle("#Delta#varphi_{h-h}")
                h_h_2d_dist.GetXaxis().SetTitle("#Delta#eta_{h-h}")
                h_h_2d_dist.GetXaxis().SetTitleSize(0.05)
                h_h_2d_dist.GetXaxis().SetTitleOffset(1.3)
                h_h_2d_dist.GetYaxis().SetTitleSize(0.05)
                h_h_2d_dist.GetYaxis().SetTitleOffset(1.3)


                c.SetTheta(18)
                c.SetPhi(50)
                c.SetRightMargin(0.1)
                c.SetTopMargin(0.15)
                c.SetLeftMargin(0.13)

                h_lambda_2d_dist.Draw("SURF1")

                if FANCY_LABEL:
                    label_y_start = 0.96 
                    label_x_start = 0.68
                    pt_label_x_start = 0.07
                    label_text_space = 0.05
                    alice_data_label = rt.TLatex()
                    alice_data_label.SetNDC()
                    alice_data_label.SetTextSize(0.04)
                    alice_data_label.SetTextAlign(13)
                    alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Performance")
                    alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p-Pb #sqrt{s_{NN}} = 5.02 TeV}")
                    if mult_string == "_0_20":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{0-20% V0A centrality}")
                    elif mult_string == "_20_50":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{20-50% V0A centrality}")
                    elif mult_string == "_50_80":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{50-80% V0A centrality}")
                    alice_data_label.DrawLatex(pt_label_x_start, label_y_start, "#bf{4.0 < #it{p}_{T, trig}^{h} < 8.0 GeV/c}")
                    if PT_MODE == 0:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{2.0 < #it{p}_{T, assoc}^{#Lambda} < 4.0 GeV/c}")
                    elif PT_MODE == 1:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{1.5 < #it{p}_{T, assoc}^{#Lambda} < 2.5 GeV/c}")
                    elif PT_MODE == 2:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{2.5 < #it{p}_{T, assoc}^{#Lambda} < 4.0 GeV/c}")

                if FANCY_LABEL:
                    if PT_MODE == 0:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor_fancy_label" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor_fancy_label" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed_fancy_label" + mult_string + ".pdf")
                    if PT_MODE == 1:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor_fancy_label" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor_fancy_label" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed_fancy_label" + mult_string + "_lowpt.pdf")
                    if PT_MODE == 2:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor_fancy_label" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor_fancy_label" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed_fancy_label" + mult_string + "_highpt.pdf")
                else:
                    if PT_MODE == 0:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed" + mult_string + ".pdf")
                    if PT_MODE == 1:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed" + mult_string + "_lowpt.pdf")
                    if PT_MODE == 2:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_lambda_2d_mixcor" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_lambda_2d_nomixcor" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_lambda_2d_mixed" + mult_string + "_highpt.pdf")

                h_h_2d_dist.Draw("SURF1")
                if FANCY_LABEL:
                    label_y_start = 0.96 
                    label_x_start = 0.68
                    pt_label_x_start = 0.07
                    label_text_space = 0.05
                    alice_data_label = rt.TLatex()
                    alice_data_label.SetNDC()
                    alice_data_label.SetTextSize(0.04)
                    alice_data_label.SetTextAlign(13)
                    alice_data_label.DrawLatex(label_x_start, label_y_start, "ALICE Performance")
                    alice_data_label.DrawLatex(label_x_start, label_y_start - label_text_space, "#bf{p-Pb #sqrt{s_{NN}} = 5.02 TeV}")
                    if mult_string == "_0_20":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{0-20% V0A centrality}")
                    elif mult_string == "_20_50":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{20-50% V0A centrality}")
                    elif mult_string == "_50_80":
                        alice_data_label.DrawLatex(label_x_start, label_y_start - 2*label_text_space, "#bf{50-80% V0A centrality}")
                    alice_data_label.DrawLatex(pt_label_x_start, label_y_start, "#bf{4.0 < #it{p}_{T, trig}^{h} < 8.0 GeV/c}")
                    if PT_MODE == 0:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{2.0 < #it{p}_{T, assoc}^{h} < 4.0 GeV/c}")
                    elif PT_MODE == 1:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{1.5 < #it{p}_{T, assoc}^{h} < 2.5 GeV/c}")
                    elif PT_MODE == 2:
                        alice_data_label.DrawLatex(pt_label_x_start, label_y_start - 1*label_text_space - 0.01, "#bf{2.5 < #it{p}_{T, assoc}^{h} < 4.0 GeV/c}")

                if FANCY_LABEL:
                    if PT_MODE == 0:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor_fancy_label" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor_fancy_label" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed_fancy_label" + mult_string + ".pdf")
                    if PT_MODE == 1:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor_fancy_label" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor_fancy_label" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed_fancy_label" + mult_string + "_lowpt.pdf")
                    if PT_MODE == 2:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor_fancy_label" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor_fancy_label" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed_fancy_label" + mult_string + "_highpt.pdf")
                else:
                    if PT_MODE == 0:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor" + mult_string + ".pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed" + mult_string + ".pdf")
                    if PT_MODE == 1:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor" + mult_string + "_lowpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed" + mult_string + "_lowpt.pdf")
                    if PT_MODE == 2:
                        if MIX_COR_MODE == 0:
                            c.SaveAs("figures/h_h_2d_mixcor" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 1:
                            c.SaveAs("figures/h_h_2d_nomixcor" + mult_string + "_highpt.pdf")
                        elif MIX_COR_MODE == 2:
                            c.SaveAs("figures/h_h_2d_mixed" + mult_string + "_highpt.pdf")
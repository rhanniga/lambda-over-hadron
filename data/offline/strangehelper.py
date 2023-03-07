# a function to calculate parabola that passes through all three input points
def get_parabola(point_one, point_two, point_three):
    x1, x2, x3 = point_one[0], point_two[0], point_three[0]
    y1, y2, y3 = point_one[1], point_two[1], point_three[1]
    
    denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
    
    A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    B = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
    C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
    
    return C, B, A

# a function to calculate the straight line that passes through all two input points
def get_straight_line(point_one, point_two):
    x1, x2 = point_one[0], point_two[0]
    y1, y2 = point_one[1], point_two[1]
    
    m = (y2 - y1) / (x2 - x1)
    c = y1 - m * x1
    
    return c, m

# a function to apply the two-track correction template to the input th2d
def apply_two_track_correction(uncorrected_dist, correction_template):

    # make a copy of the input distribution
    corrected_dist = uncorrected_dist.Clone()

    deta_min_bin = correction_template.GetXaxis().FindBin(-1.4)
    deta_max_bin = correction_template.GetXaxis().FindBin(1.4 - 0.001)

    # loop over all bins in the input distribution
    for xbin in range(deta_min_bin, deta_max_bin+1):
        for ybin in range(uncorrected_dist.GetNbinsY()):
            # get the bin content of the input distribution
            bin_content = uncorrected_dist.GetBinContent(xbin, ybin+1)
            bin_error = uncorrected_dist.GetBinError(xbin, ybin+1)
            
            # get the correction factor from the correction template
            correction_factor = 1/correction_template.GetBinContent(xbin, ybin+1)
            
            # apply the correction factor to the bin content
            corrected_bin_content = bin_content * correction_factor
            corrected_bin_error = bin_error * correction_factor
            
            # set the bin content of the output distribution
            corrected_dist.SetBinContent(xbin, ybin+1, corrected_bin_content)
            corrected_dist.SetBinError(xbin, ybin+1, corrected_bin_error)

    return corrected_dist

# a function to perform mixed event acceptance correction
def make_mixed_corrections(same, mixed, mass_low=1.11, mass_high=1.12, is_hh=False):
    if is_hh:
        same3d = same
        same3d.Sumw2()
        mixed3d = mixed
        mixed3d.Sumw2()
    else:
        same.GetAxis(2).SetRangeUser(mass_low, mass_high)
        mixed.GetAxis(2).SetRangeUser(mass_low, mass_high)
        same3d = same.Projection(0, 1, 3)
        same3d.Sumw2()
        mixed3d = mixed.Projection(0, 1, 3)
        mixed3d.Sumw2()
        
    for zbin in range(10):
        same3d.GetZaxis().SetRange(zbin+1, zbin+1)
        same2d = same3d.Project3D("xye")
        same2d.SetName(f"same2dproj_zbin_{zbin}")

        mixed3d.GetZaxis().SetRange(zbin+1, zbin+1)
        mixed2d = mixed3d.Project3D("xye")
        mixed2d.SetName(f"mix2dproj_zbin_{zbin}")

        #scaling by average of bins adjacent to 0
        scale = 0.5*(mixed2d.Integral(mixed2d.GetXaxis().FindBin(-0.01),    #xmin
                                    mixed2d.GetXaxis().FindBin(0.01),     #xmax 
                                    mixed2d.GetYaxis().FindBin(0.0),      #ymin
                                    mixed2d.GetYaxis().FindBin(0.0)))     #ymax
        same2d.Divide(mixed2d)
        same2d.Scale(scale)
        
        if zbin == 0:
            same2d_total = same2d.Clone("2dproj_total")
            mixed2d_total = mixed2d.Clone("2dproj_total_mixed")
        else:
            same2d_total.Add(same2d)
            mixed2d_total.Add(mixed2d)

    return same2d_total
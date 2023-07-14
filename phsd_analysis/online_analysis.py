
import ROOT as rt
import math

from multiprocessing import Process as worker

from sys import argv
from array import array

# constants for various things
COL_ENERGY = 4000 + (82/208)*4000
COL_MOMENTUM = 4000 - (82/208)*4000
COL_LVECTOR = rt.TLorentzVector(0, 0, COL_MOMENTUM, COL_ENERGY)

# particle class, mostly for keeping track of indices within event stack
class Particle:
    def __init__(self, pdg, pt, eta, phi, index):
        self.pdg = pdg
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.index = index

# main class for analysis
class PHSDAnalyzer:

    def __init__(self):
        self.initialize_output_objects()

    def initialize_output_objects(self):

        self.output_list = rt.TList()
        self.output_list.SetOwner(True)

        # mult dist (in V0A acceptance) for determining percentiles
        self.mult_dist = rt.TH1D("fMultDist", "Charged pat. dist in V0A acceptance", 1000, 0, 1000)
        self.output_list.Add(self.mult_dist)

        # single particle distributions, axes are (pt, phi, eta, zvtx = 0 since PHSD has no PV info, mult bin)
        single_dist_bins = array("i", [100, 16, 100, 1, 4])
        single_dist_mins = array("d", [0.0, 0, -13, 0, 0])
        single_dist_maxes = array("d", [10.0, 6.28, 13, 1, 4])

        self.trigger_dist = rt.THnSparseD("fTriggerDist_MC", "Trigger distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.trigger_dist.Sumw2()
        self.output_list.Add(self.trigger_dist)

        self.hadron_dist = rt.THnSparseD("fAssociatedDist_MC", "Assoc. hadron distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.hadron_dist.Sumw2()
        self.output_list.Add(self.hadron_dist)

        self.lambda_dist = rt.THnSparseD("fLambdaDist_MC", "Lambda distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_dist.Sumw2()
        self.output_list.Add(self.lambda_dist)

        self.lambda_from_sigma_dist = rt.THnSparseD("fLambdaFromSigmaDist_MC", "Lambda distribution (from sigmas)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_sigma_dist.Sumw2()
        self.output_list.Add(self.lambda_from_sigma_dist)

        self.lambda_from_omega_dist = rt.THnSparseD("fLambdaFromOmegaDist_MC", "Lambda distribution (from omegas)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_omega_dist.Sumw2()
        self.output_list.Add(self.lambda_from_omega_dist)

        self.lambda_from_xi_dist = rt.THnSparseD("fLambdaFromXiDist_MC", "Lambda distribution (from xis)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_xi_dist.Sumw2()
        self.output_list.Add(self.lambda_from_xi_dist)

        self.phi_dist = rt.THnSparseD("fPhiDist_MC", "Phi distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.phi_dist.Sumw2()
        self.output_list.Add(self.phi_dist)

        self.trigger_dist_no_eta_cut = rt.THnSparseD("fTriggerDist_MC_no_eta_cut", "Trigger distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.trigger_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.trigger_dist_no_eta_cut)

        self.hadron_dist_no_eta_cut = rt.THnSparseD("fAssociatedDist_MC_no_eta_cut", "Assoc. hadron distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.hadron_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.hadron_dist_no_eta_cut)

        self.lambda_dist_no_eta_cut = rt.THnSparseD("fLambdaDist_MC_no_eta_cut", "Lambda distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.lambda_dist_no_eta_cut)

        self.lambda_from_sigma_dist_no_eta_cut = rt.THnSparseD("fLambdaFromSigmaDist_MC_no_eta_cut", "Lambda distribution (from sigmas,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_sigma_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.lambda_from_sigma_dist_no_eta_cut)

        self.lambda_from_omega_dist_no_eta_cut = rt.THnSparseD("fLambdaFromOmegaDist_MC_no_eta_cut", "Lambda distribution (from omegas,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_omega_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.lambda_from_omega_dist_no_eta_cut)

        self.lambda_from_xi_dist_no_eta_cut = rt.THnSparseD("fLambdaFromXiDist_MC_no_eta_cut", "Lambda distribution (from xis,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.lambda_from_xi_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.lambda_from_xi_dist_no_eta_cut)

        self.phi_dist_no_eta_cut = rt.THnSparseD("fPhiDist_MC_no_eta_cut", "Phi distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.phi_dist_no_eta_cut.Sumw2()
        self.output_list.Add(self.phi_dist_no_eta_cut)

        self.triggered_trigger_dist = rt.THnSparseD("fTriggeredTriggerDist_MC", "Trigger distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.triggered_trigger_dist.Sumw2()
        self.output_list.Add(self.triggered_trigger_dist)

        self.triggered_associated_dist = rt.THnSparseD("fTriggeredAssociatedDist_MC", "Associated distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.triggered_associated_dist.Sumw2()
        self.output_list.Add(self.triggered_associated_dist)

        self.triggered_lambda_dist = rt.THnSparseD("fTriggeredLambdaDist_MC", "Lambda distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.triggered_lambda_dist.Sumw2()
        self.output_list.Add(self.triggered_lambda_dist)

        self.triggered_phi_dist = rt.THnSparseD("fTriggeredPhiDist_MC", "Phi distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
        self.triggered_phi_dist.Sumw2()
        self.output_list.Add(self.triggered_phi_dist)

        # correlation distributions axes are (trigger pt, associated pt, delta phi, delta eta, zvtx = 0, mult)
        correlation_dist_bins = array("i", [18, 16, 16, 20, 1, 4])
        correlation_dist_mins = array("d", [1, 0, -rt.TMath.Pi()/2, -2, 0, 0])
        correlation_dist_maxes = array("d", [10, 4, 3*rt.TMath.Pi()/2, 2, 1, 4])

        self.h_h_dist = rt.THnSparseD("fDphiHH_MC", "H-hadron dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes)
        self.h_h_dist.Sumw2()
        self.output_list.Add(self.h_h_dist)

        self.h_lambda_dist = rt.THnSparseD("fDphiHLambda_MC", "H-lambda dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes)
        self.h_lambda_dist.Sumw2()
        self.output_list.Add(self.h_lambda_dist)

        self.h_primary_lambda_dist = rt.THnSparseD("fDphiHPrimaryLambda_MC", "H-primarylambda dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes)
        self.h_primary_lambda_dist.Sumw2()
        self.output_list.Add(self.h_primary_lambda_dist)

        self.h_phi_dist = rt.THnSparseD("fDphiHPhi_MC", "H-phi dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes) 
        self.h_phi_dist.Sumw2()
        self.output_list.Add(self.h_phi_dist)

    # check PDG code to see if particle is charged hadron (or electron/muon, assuming production has been enabled)
    def is_charged_hadron(self, pdg):
        if rt.TMath.Abs(pdg) == 211 or rt.TMath.Abs(pdg) == 321 or rt.TMath.Abs(pdg) == 2212 or rt.TMath.Abs(pdg) == 11 or rt.TMath.Abs(pdg) == 13:
            return True
        else:
            return False

    # get multiplicity bin from pre-calculated percentiles
    def get_mult_bin(self, num_tracks_v0a):
        if num_tracks_v0a > 51.5: return 3 # 0-20%
        elif 28.5 < num_tracks_v0a < 51.5: return 2 # 20-50%
        elif 12.5 < num_tracks_v0a < 28.5: return 1 # 50-80%
        elif num_tracks_v0a < 12.5: return 0 # 80-100%

    def fill_single_dist(self, particle_list, dist, mult_bin):
        for particle in particle_list:
            single_dist_array = array("d", [particle.pt, particle.phi, particle.eta, 0.5, mult_bin])
            dist.Fill(single_dist_array)
        
    # fill correlation distribution for given trigger and associated particle lists
    def fill_correlation_dist(self, trigger_list, associated_list, dist, mult_bin):
        for trigger in trigger_list:
            for associated in associated_list:

                # do not correlate particle with itself
                if trigger.index == associated.index:
                    continue

                # keep dphi in range [-pi/2, 3pi/2]
                delta_phi = trigger.phi - associated.phi
                if delta_phi < -rt.TMath.Pi()/2:
                    delta_phi += 2*rt.TMath.Pi()
                elif delta_phi > 3*rt.TMath.Pi()/2:
                    delta_phi -= 2*rt.TMath.Pi()

                delta_eta = trigger.eta - associated.eta

                correlation_dist_array = array("d", [trigger.pt, associated.pt, delta_phi, delta_eta, 0.5, mult_bin])
                dist.Fill(correlation_dist_array)



    def process_event(self, event):

        trigger_list = []

        hadron_list = []
        lambda_list = []
        lambda_from_sigma_list = []
        lambda_from_xi_list = []
        lambda_from_omega_list = []
        phi_list = []

        hadron_list_no_eta_cut = []
        lambda_list_no_eta_cut = []
        lambda_from_sigma_list_no_eta_cut = []
        lambda_from_xi_list_no_eta_cut = []
        lambda_from_omega_list_no_eta_cut = []
        phi_list_no_eta_cut = []

        # counting number of charged hadrons in V0A acceptance
        num_tracks_in_v0a_acceptance = 0
            
        # bool to keep track of whether or not event has trigger (> 4 GeV charged)
        is_triggered_event = False

        # loop over tracks in event
        for index, track in enumerate(event):

            track_info = track.split()

            track_pdg = int(track_info[0])

            track_lorentz_vector = rt.TLorentzVector(float(track_info[2]), float(track_info[3]), -float(track_info[4]), float(track_info[5]))
            track_lorentz_vector.Boost(-COL_LVECTOR.BoostVector())

            track_pt = track_lorentz_vector.Pt()
            track_eta = track_lorentz_vector.Eta()
            track_phi = track_lorentz_vector.Phi()

            particle = Particle(track_pdg, track_pt, track_eta, track_phi, index)

            if particle.pt <= 0.15:
                continue


            if rt.TMath.Abs(particle.pdg) == 3122: # lambdas
                lambda_list_no_eta_cut.append(particle)
                if rt.TMath.Abs(particle.eta) < 0.8:
                    lambda_list.append(particle)

            elif rt.TMath.Abs(particle.pdg) == 3212: # sigmas
                decay_event = rt.TGenPhaseSpace()
                decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.0]))
                decay_event.Generate()
                lambda_lorentz_vector = decay_event.GetDecay(0)

                lambda_pt = lambda_lorentz_vector.Pt()
                lambda_eta = lambda_lorentz_vector.Eta()
                lambda_phi = lambda_lorentz_vector.Phi()

                lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)

                lambda_from_sigma_list_no_eta_cut.append(lambda_particle)
                if rt.TMath.Abs(lambda_particle.eta) < 0.8:
                    lambda_from_sigma_list.append(lambda_particle)

            elif rt.TMath.Abs(particle.pdg) == 3334: # omegas
                decay_event = rt.TGenPhaseSpace()
                rand_number = rt.gRandom.Uniform(0, 1)
                if rand_number <= 0.678:
                    decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.493677]))
                    decay_event.Generate()
                    lambda_lorentz_vector = decay_event.GetDecay(0)

                    lambda_pt = lambda_lorentz_vector.Pt()
                    lambda_eta = lambda_lorentz_vector.Eta()
                    lambda_phi = lambda_lorentz_vector.Phi()

                    lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)

                    lambda_from_omega_list_no_eta_cut.append(lambda_particle)
                    if rt.TMath.Abs(lambda_particle.eta) < 0.8:
                        lambda_from_omega_list.append(lambda_particle)

            elif rt.TMath.Abs(particle.pdg) == 3322: # neutral xis
                decay_event = rt.TGenPhaseSpace()
                decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.134976]))
                decay_event.Generate()
                lambda_lorentz_vector = decay_event.GetDecay(0)

                lambda_pt = lambda_lorentz_vector.Pt()
                lambda_eta = lambda_lorentz_vector.Eta()
                lambda_phi = lambda_lorentz_vector.Phi()

                lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)

                lambda_from_xi_list_no_eta_cut.append(lambda_particle)
                if rt.TMath.Abs(lambda_particle.eta) < 0.8:
                    lambda_from_xi_list.append(lambda_particle)
            
            elif rt.TMath.Abs(particle.pdg) == 3312: # charged xis
                decay_event = rt.TGenPhaseSpace()
                decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.139570]))
                decay_event.Generate()
                lambda_lorentz_vector = decay_event.GetDecay(0)

                lambda_pt = lambda_lorentz_vector.Pt()
                lambda_eta = lambda_lorentz_vector.Eta()
                lambda_phi = lambda_lorentz_vector.Phi()

                lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)

                lambda_from_xi_list_no_eta_cut.append(lambda_particle)
                if rt.TMath.Abs(lambda_particle.eta) < 0.8:
                    lambda_from_xi_list.append(lambda_particle)

            elif rt.TMath.Abs(particle.pdg) == 333: # phis
                phi_list_no_eta_cut.append(particle)
                if rt.TMath.Abs(particle.eta) < 0.8:
                    phi_list.append(particle)

            elif self.is_charged_hadron(particle.pdg):
                hadron_list_no_eta_cut.append(particle)
                if rt.TMath.Abs(particle.eta) < 0.8:
                    trigger_list.append(particle)
                    hadron_list.append(particle)

                    if particle.pt > 4:
                        is_triggered_event = True

                if particle.eta > 2.8 and particle.eta < 5.1:
                    num_tracks_in_v0a_acceptance += 1

        self.mult_dist.Fill(num_tracks_in_v0a_acceptance)

        mult_bin = self.get_mult_bin(num_tracks_in_v0a_acceptance)

        lambda_list_total = lambda_list + lambda_from_sigma_list + lambda_from_omega_list + lambda_from_xi_list
        lambda_list_total_no_eta_cut = lambda_list_no_eta_cut + lambda_from_sigma_list_no_eta_cut + lambda_from_omega_list_no_eta_cut + lambda_from_xi_list_no_eta_cut

        self.fill_single_dist(trigger_list, self.trigger_dist, mult_bin)
        self.fill_single_dist(hadron_list, self.hadron_dist, mult_bin)
        self.fill_single_dist(lambda_list, self.lambda_dist, mult_bin)
        self.fill_single_dist(lambda_from_sigma_list, self.lambda_from_sigma_dist, mult_bin)
        self.fill_single_dist(lambda_from_omega_list, self.lambda_from_omega_dist, mult_bin)
        self.fill_single_dist(lambda_from_xi_list, self.lambda_from_xi_dist, mult_bin)

        self.fill_single_dist(hadron_list_no_eta_cut, self.hadron_dist_no_eta_cut, mult_bin)
        self.fill_single_dist(lambda_list_no_eta_cut, self.lambda_dist_no_eta_cut, mult_bin)
        self.fill_single_dist(lambda_from_sigma_list_no_eta_cut, self.lambda_from_sigma_dist_no_eta_cut, mult_bin)
        self.fill_single_dist(lambda_from_omega_list_no_eta_cut, self.lambda_from_omega_dist_no_eta_cut, mult_bin)
        self.fill_single_dist(lambda_from_xi_list_no_eta_cut, self.lambda_from_xi_dist_no_eta_cut, mult_bin)


        self.fill_correlation_dist(trigger_list, lambda_list_no_eta_cut, self.h_primary_lambda_dist, mult_bin)
        self.fill_correlation_dist(trigger_list, lambda_list_total_no_eta_cut, self.h_lambda_dist, mult_bin)
        self.fill_correlation_dist(trigger_list, hadron_list_no_eta_cut, self.h_h_dist, mult_bin)
        self.fill_correlation_dist(trigger_list, phi_list_no_eta_cut, self.h_phi_dist, mult_bin)

    def parse_and_process(self, file_name):

        # we have to read line by line since the file is so large
        with open(file_name, "r") as f:
            event = []
            line_num = 0
            num_tracks = -999
            event_num = 0
            for line in f:
                if line_num == 0 or line_num == num_tracks + 2:
                    num_tracks = int(line.split()[0])
                    line_num = 0
                    self.process_event(event)
                    if event_num % 1000 == 0:
                        print("Processed event " + str(event_num))
                    event_num += 1
                    event = []
                elif line_num > 1:
                    # ignore elastic tracks
                    if line.split()[6] == "-1":
                        line_num += 1
                        continue
                    event.append(line)
                line_num += 1

    def write_out(self, starting_file_name):

        outfile_name = starting_file_name
        outfile_name = outfile_name[6:]
        outfile_name = "root_out" + outfile_name
        outfile_name = outfile_name[:-3]
        outfile_name = outfile_name + "root"

        outfile = rt.TFile(outfile_name, "RECREATE")

        outfile.WriteObject(self.output_list, "h-lambda")

        outfile.Close()


if __name__ == "__main__":

    # initialize analyzer
    phsd_analyzer = PHSDAnalyzer()

    # parse and process the input file
    phsd_analyzer.parse_and_process(file_name=argv[1])

    # write to output file
    phsd_analyzer.write_out(starting_file_name=argv[1])

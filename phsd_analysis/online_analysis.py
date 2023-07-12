# this is going to be ugly, but don't let it bother you

import ROOT as rt
import math

from multiprocessing import Process as worker

from sys import argv
from array import array


# particle class, mostly for keeping track of indices within event stack
class Particle:
    def __init__(self, pdg, pt, eta, phi, index):
        self.pdg = pdg
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.index = index

# fill correlation distribution for given trigger and associated particle lists
def fill_correlation_dist(trigger_list, associated_list, dist, mult_bin):
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

# return [pt, eta, phi] tuple from input momentum list (px, py, pz)
def get_kinematics(p):
    pt = rt.TMath.Sqrt(p[0]**2 + p[1]**2)
    # eta is shifted by 0.465 due to collision asymmetry
    eta = rt.TMath.ATanH(p[2]/(rt.TMath.Sqrt(p[0]**2 + p[1]**2 + p[2]**2))) + 0.465
    # eta is reversed from LHC convention
    eta *= -1
    phi = rt.TMath.ATan2(p[1], p[0])
    if phi < 0:
        phi += 2*rt.TMath.Pi()
    elif phi > 2*rt.TMath.Pi():
        phi -= 2*rt.TMath.Pi()

    return [pt, eta, phi]

# check PDG code to see if particle is charged hadron (or electron/muon, assuming production has been enabled)
def is_charged_hadron(pdg):
    if rt.TMath.Abs(pdg) == 211 or rt.TMath.Abs(pdg) == 321 or rt.TMath.Abs(pdg) == 2212 or rt.TMath.Abs(pdg) == 11 or rt.TMath.Abs(pdg) == 13:
        return True
    else:
        return False

# mult dist (in V0A acceptance) for determining percentiles
mult_dist = rt.TH1D("fMultDist", "Charged pat. dist in V0A acceptance", 1000, 0, 1000)

# single particle distributions, axes are (pt, phi, eta, zvtx = 0 since PHSD has no PV info)
single_dist_bins = array("i", [100, 16, 100, 1, 4])
single_dist_mins = array("d", [0.0, 0, -13, 0, 0])
single_dist_maxes = array("d", [10.0, 6.28, 13, 1, 4])

trigger_dist = rt.THnSparseD("fTriggerDist_MC", "Trigger distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
trigger_dist.Sumw2()
associated_dist = rt.THnSparseD("fAssociatedDist_MC", "Associated distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
associated_dist.Sumw2()
lambda_dist = rt.THnSparseD("fLambdaDist_MC", "Lambda distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
lambda_dist.Sumw2()
phi_dist = rt.THnSparseD("fPhiDist_MC", "Phi distribution", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
phi_dist.Sumw2()

trigger_dist_no_eta_cut = rt.THnSparseD("fTriggerDist_MC_no_eta_cut", "Trigger distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
trigger_dist_no_eta_cut.Sumw2()
associated_dist_no_eta_cut = rt.THnSparseD("fAssociatedDist_MC_no_eta_cut", "Associated distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
associated_dist_no_eta_cut.Sumw2()
lambda_dist_no_eta_cut = rt.THnSparseD("fLambdaDist_MC_no_eta_cut", "Lambda distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
lambda_dist_no_eta_cut.Sumw2()
lambda_from_sigma_dist_no_eta_cut = rt.THnSparseD("fLambdaFromSigmaDist_MC_no_eta_cut", "Lambda distribution (from sigmas,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
lambda_from_sigma_dist_no_eta_cut.Sumw2()
lambda_from_omega_dist_no_eta_cut = rt.THnSparseD("fLambdaFromOmegaDist_MC_no_eta_cut", "Lambda distribution (from omegas,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
lambda_from_omega_dist_no_eta_cut.Sumw2()
lambda_from_xi_dist_no_eta_cut = rt.THnSparseD("fLambdaFromXiDist_MC_no_eta_cut", "Lambda distribution (from xis,no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
lambda_from_xi_dist_no_eta_cut.Sumw2()
phi_dist_no_eta_cut = rt.THnSparseD("fPhiDist_MC_no_eta_cut", "Phi distribution (no eta cut)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
phi_dist_no_eta_cut.Sumw2()
triggered_trigger_dist = rt.THnSparseD("fTriggeredTriggerDist_MC", "Trigger distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
triggered_trigger_dist.Sumw2()
triggered_associated_dist = rt.THnSparseD("fTriggeredAssociatedDist_MC", "Associated distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
triggered_associated_dist.Sumw2()
triggered_lambda_dist = rt.THnSparseD("fTriggeredLambdaDist_MC", "Lambda distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
triggered_lambda_dist.Sumw2()
triggered_phi_dist = rt.THnSparseD("fTriggeredPhiDist_MC", "Phi distribution (triggered event)", 5, single_dist_bins, single_dist_mins, single_dist_maxes)
triggered_phi_dist.Sumw2()

# correlation distributions axes are (trigger pt, associated pt, delta phi, delta eta, zvtx = 0, mult)
correlation_dist_bins = array("i", [18, 16, 16, 20, 1, 4])
correlation_dist_mins = array("d", [1, 0, -rt.TMath.Pi()/2, -2, 0, 0])
correlation_dist_maxes = array("d", [10, 4, 3*rt.TMath.Pi()/2, 2, 1, 4])

h_h_dist = rt.THnSparseD("fDphiHH_MC", "H-hadron dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes)
h_h_dist.Sumw2()
h_lambda_dist = rt.THnSparseD("fDphiHLambda_MC", "H-lambda dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes)
h_lambda_dist.Sumw2()
h_phi_dist = rt.THnSparseD("fDphiHPhi_MC", "H-phi dist", 6, correlation_dist_bins, correlation_dist_mins, correlation_dist_maxes) 
h_phi_dist.Sumw2()

def get_mult_bin(num_tracks_v0a):
    if num_tracks_v0a > 51.5: return 3
    elif 28.5 < num_tracks_v0a < 51.5: return 2
    elif 12.5 < num_tracks_v0a < 28.5: return 1
    elif num_tracks_v0a < 12.5: return 0

def process_event(event):

    trigger_list = []

    associated_list_no_eta_cut = []
    lambda_list_no_eta_cut = []
    phi_list_no_eta_cut = []


    collision_energy = 4000 + (82/208)*4000
    collision_momentum = 4000 - (82/208)*4000

    collision_lorentz_vector = rt.TLorentzVector(0, 0, collision_momentum, collision_energy)

    # counting number of charged hadrons in V0A acceptance
    num_tracks_in_v0a_acceptance = 0

    # loop over tracks in event to get number of charged hadrons in V0A acceptance
    for index, track in enumerate(event):

        track_info = track.split()
        track_pdg = int(track_info[0])
        track_fourmomentum = [float(track_info[2]), float(track_info[3]), -float(track_info[4]), float(track_info[5])]
        track_lorentz_vector = rt.TLorentzVector(track_fourmomentum[0], track_fourmomentum[1], track_fourmomentum[2], track_fourmomentum[3])
        track_lorentz_vector.Boost(-collision_lorentz_vector.BoostVector())
        track_pt = track_lorentz_vector.Pt()
        track_eta = track_lorentz_vector.Eta()
        track_phi = track_lorentz_vector.Phi()
        particle = Particle(track_pdg, track_pt, track_eta, track_phi, index)

        if particle.pt <= 0.15:
            continue
        if is_charged_hadron(particle.pdg):
            if particle.eta > 2.8 and particle.eta < 5.1:
                num_tracks_in_v0a_acceptance += 1


    mult_dist.Fill(num_tracks_in_v0a_acceptance)

    mult_bin = get_mult_bin(num_tracks_in_v0a_acceptance)
        
    # bool to keep track of whether or not event has trigger (> 4 GeV charged)
    is_triggered_event = False

    # loop over tracks in event
    for index, track in enumerate(event):
        track_info = track.split()
        track_pdg = int(track_info[0])
        track_fourmomentum = [float(track_info[2]), float(track_info[3]), -float(track_info[4]), float(track_info[5])]
        track_lorentz_vector = rt.TLorentzVector(track_fourmomentum[0], track_fourmomentum[1], track_fourmomentum[2], track_fourmomentum[3])
        track_lorentz_vector.Boost(-collision_lorentz_vector.BoostVector())
        track_pt = track_lorentz_vector.Pt()
        track_eta = track_lorentz_vector.Eta()
        track_phi = track_lorentz_vector.Phi()
        particle = Particle(track_pdg, track_pt, track_eta, track_phi, index)

        # track_info = track.split()
        # track_pdg = int(track_info[0])
        # track_momentum = [float(track_info[2]), float(track_info[3]), float(track_info[4])]
        # track_kinematics = get_kinematics(track_momentum)
        # track_pt = track_kinematics[0]
        # track_eta = track_kinematics[1]
        # track_phi = track_kinematics[2]
        # particle = Particle(track_pdg, track_pt, track_eta, track_phi, index)

        if particle.pt <= 0.15:
            continue

        single_dist_array = array("d", [particle.pt, particle.phi, particle.eta, 0.5, mult_bin])

        if rt.TMath.Abs(particle.pdg) == 3122: # lambdas
            lambda_list_no_eta_cut.append(particle)
            lambda_dist_no_eta_cut.Fill(single_dist_array)
            if rt.TMath.Abs(particle.eta) < 0.8:
                lambda_dist.Fill(single_dist_array)

        elif rt.TMath.Abs(particle.pdg) == 3212: # sigmas
            decay_event = rt.TGenPhaseSpace()
            decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.0]))
            decay_event.Generate()
            lambda_lorentz_vector = decay_event.GetDecay(0)

            # lambda_momentum = [lambda_lorentz_vector.Px(), lambda_lorentz_vector.Py(), lambda_lorentz_vector.Pz()]
            # lambda_kinematics = get_kinematics(lambda_momentum)

            lambda_pt = lambda_lorentz_vector.Pt()
            lambda_eta = lambda_lorentz_vector.Eta()
            lambda_phi = lambda_lorentz_vector.Phi()

            lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)
            lambda_list_no_eta_cut.append(lambda_particle)
            lambda_from_sigma_dist_no_eta_cut.Fill(array("d", [lambda_particle.pt, lambda_particle.phi, lambda_particle.eta, 0.5]))

        elif rt.TMath.Abs(particle.pdg) == 3334: # omegas
            decay_event = rt.TGenPhaseSpace()
            rand_number = rt.gRandom.Uniform(0, 1)
            if rand_number <= 0.678:
                decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.493677]))
                decay_event.Generate()
                lambda_lorentz_vector = decay_event.GetDecay(0)

                # lambda_momentum = [lambda_lorentz_vector.Px(), lambda_lorentz_vector.Py(), lambda_lorentz_vector.Pz()]
                # lambda_kinematics = get_kinematics(lambda_momentum)

                lambda_pt = lambda_lorentz_vector.Pt()
                lambda_eta = lambda_lorentz_vector.Eta()
                lambda_phi = lambda_lorentz_vector.Phi()

                lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)
                lambda_list_no_eta_cut.append(lambda_particle)
                lambda_from_omega_dist_no_eta_cut.Fill(array("d", [lambda_particle.pt, lambda_particle.phi, lambda_particle.eta, 0.5]))

        elif rt.TMath.Abs(particle.pdg) == 3322: # neutral xis
            decay_event = rt.TGenPhaseSpace()
            decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.134976]))
            decay_event.Generate()
            lambda_lorentz_vector = decay_event.GetDecay(0)

            # lambda_momentum = [lambda_lorentz_vector.Px(), lambda_lorentz_vector.Py(), lambda_lorentz_vector.Pz()]
            # lambda_kinematics = get_kinematics(lambda_momentum)
            # lambda_pt = lambda_kinematics[0]
            # lambda_eta = lambda_kinematics[1]
            # lambda_phi = lambda_kinematics[2]

            lambda_pt = lambda_lorentz_vector.Pt()
            lambda_eta = lambda_lorentz_vector.Eta()
            lambda_phi = lambda_lorentz_vector.Phi()

            lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)
            lambda_list_no_eta_cut.append(lambda_particle)
            lambda_from_xi_dist_no_eta_cut.Fill(array("d", [lambda_particle.pt, lambda_particle.phi, lambda_particle.eta, 0.5]))
        
        elif rt.TMath.Abs(particle.pdg) == 3312: # charged xis
            decay_event = rt.TGenPhaseSpace()
            decay_event.SetDecay(track_lorentz_vector, 2, array("d", [1.115683, 0.139570]))
            decay_event.Generate()
            lambda_lorentz_vector = decay_event.GetDecay(0)

            # lambda_momentum = [lambda_lorentz_vector.Px(), lambda_lorentz_vector.Py(), lambda_lorentz_vector.Pz()]
            # lambda_kinematics = get_kinematics(lambda_momentum)
            # lambda_pt = lambda_kinematics[0]
            # lambda_eta = lambda_kinematics[1]
            # lambda_phi = lambda_kinematics[2]

            lambda_pt = lambda_lorentz_vector.Pt()
            lambda_eta = lambda_lorentz_vector.Eta()
            lambda_phi = lambda_lorentz_vector.Phi()

            lambda_particle = Particle(3122, lambda_pt, lambda_eta, lambda_phi, index)
            lambda_list_no_eta_cut.append(lambda_particle)
            lambda_from_xi_dist_no_eta_cut.Fill(array("d", [lambda_particle.pt, lambda_particle.phi, lambda_particle.eta, 0.5]))

        elif rt.TMath.Abs(particle.pdg) == 333: # phis
            phi_list_no_eta_cut.append(particle)
            phi_dist_no_eta_cut.Fill(single_dist_array)
            if rt.TMath.Abs(particle.eta) < 0.8:
                phi_dist.Fill(single_dist_array)

        elif is_charged_hadron(particle.pdg):
            associated_list_no_eta_cut.append(particle)
            associated_dist_no_eta_cut.Fill(single_dist_array)
            trigger_dist_no_eta_cut.Fill(single_dist_array)
            if rt.TMath.Abs(particle.eta) < 0.8:
                trigger_list.append(particle)
                trigger_dist.Fill(single_dist_array)
                associated_dist.Fill(single_dist_array)
                if particle.pt > 4:
                    is_triggered_event = True



    fill_correlation_dist(trigger_list, lambda_list_no_eta_cut, h_lambda_dist, mult_bin)
    fill_correlation_dist(trigger_list, associated_list_no_eta_cut, h_h_dist, mult_bin)
    fill_correlation_dist(trigger_list, phi_list_no_eta_cut, h_phi_dist, mult_bin)

    # too lazy to re-write code so we do this instead
    if is_triggered_event:
        # loop over tracks in event
        for index, track in enumerate(event):
            track_info = track.split()
            track_pdg = int(track_info[0])
            track_fourmomentum = [float(track_info[2]), float(track_info[3]), -float(track_info[4]), float(track_info[5])]
            track_lorentz_vector = rt.TLorentzVector(track_fourmomentum[0], track_fourmomentum[1], track_fourmomentum[2], track_fourmomentum[3])
            track_lorentz_vector.Boost(-collision_lorentz_vector.BoostVector())
            track_pt = track_lorentz_vector.Pt()
            track_eta = track_lorentz_vector.Eta()
            track_phi = track_lorentz_vector.Phi()
            particle = Particle(track_pdg, track_pt, track_eta, track_phi, index)
            if particle.pt <= 0.15:
                continue
            single_dist_array = array("d", [particle.pt, particle.phi, particle.eta, 0.5, mult_bin])
            if rt.TMath.Abs(particle.pdg) == 3122:
                if rt.TMath.Abs(particle.eta) < 0.8:
                    triggered_lambda_dist.Fill(single_dist_array)
            elif rt.TMath.Abs(particle.pdg) == 333:
                if rt.TMath.Abs(particle.eta) < 0.8:
                    triggered_phi_dist.Fill(single_dist_array)
            elif is_charged_hadron(particle.pdg):
                if rt.TMath.Abs(particle.eta) < 0.8:
                    triggered_trigger_dist.Fill(single_dist_array)
                    triggered_associated_dist.Fill(single_dist_array)



# we have to read line by line since the file is so large
with open(argv[1], "r") as f:
    event = []
    line_num = 0
    num_tracks = -999
    event_num = 0
    for line in f:
        if line_num == 0 or line_num == num_tracks + 2:
            num_tracks = int(line.split()[0])
            line_num = 0
            process_event(event)
            if event_num % 1000 == 0:
                print("Processed event " + str(event_num))
            event_num += 1

            event = []
        elif line_num > 1:
            # if line.split()[6] == "-1":
            #     line_num += 1
            #     continue
            event.append(line)
        line_num += 1

outfile_name = argv[1]
outfile_name = outfile_name[6:]
outfile_name = "root_out" + outfile_name
outfile_name = outfile_name[:-3]
outfile_name = outfile_name + "root"
outfile = rt.TFile(outfile_name, "RECREATE")
mult_dist.Write()

trigger_dist.Write()
associated_dist.Write()
lambda_dist.Write()
phi_dist.Write()

trigger_dist_no_eta_cut.Write()
associated_dist_no_eta_cut.Write()
lambda_dist_no_eta_cut.Write()
lambda_from_sigma_dist_no_eta_cut.Write()
lambda_from_omega_dist_no_eta_cut.Write()
lambda_from_xi_dist_no_eta_cut.Write()
phi_dist_no_eta_cut.Write()

triggered_trigger_dist.Write()
triggered_associated_dist.Write()
triggered_lambda_dist.Write()
triggered_phi_dist.Write()

h_h_dist.Write()
h_lambda_dist.Write()
h_phi_dist.Write()

outfile.Close()

/************************************************************************* 
 * Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. * 
 *                                                                        * 
 * Author: Ryan Hannigan
 * Contributors are mentioned in the code where appropriate.              * 
 *                                                                        * 
 * Permission to use, copy, modify and distribute this software and its   * 
 * documentation strictly for non-commercial purposes is hereby granted   * 
 * without fee, provided that the above copyright notice appears in all   * 
 * copies and that both the copyright notice and this permission notice   * 
 * appear in the supporting documentation. The authors make no claims     * 
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  * 
 **************************************************************************/

#include <iostream>
#include <set>

#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCFParticle.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliAODv0.h"


#include "AliAnalysisTaskLambdaHadronCloseTrackEff.h"

ClassImp(AliAnalysisTaskLambdaHadronCloseTrackEff);

AliAnalysisTaskLambdaHadronCloseTrackEff::AliAnalysisTaskLambdaHadronCloseTrackEff() :
    AliAnalysisTaskSE(),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fTriggerBit(0.0),
    fCorPoolMgr(0x0),
    fMCCorPoolMgr(0x0)
{
}

AliAnalysisTaskLambdaHadronCloseTrackEff::AliAnalysisTaskLambdaHadronCloseTrackEff(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fMCArray(0x0),
    fOutputList(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fCentEstimator(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fTriggerBit(0.0),
    fCorPoolMgr(0x0),
    fMCCorPoolMgr(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskLambdaHadronCloseTrackEff::~AliAnalysisTaskLambdaHadronCloseTrackEff()
{
    if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {fMultLow, fMultHigh};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    fMCCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fMCCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    //Correlation axes are: q (relative momenta), p (avg momenta), dphi, zvtx
    int hp_cor_bins[5] = {100, 100, 32, 32, 10};
    double hp_cor_mins[5] = {0, 0, -1.0*TMath::Pi()/2.0, -2, -10};
    double hp_cor_maxes[5] = {50, 50, 3.0*TMath::Pi()/2.0, 2, 10};

    fRecoHProton = new THnSparseF("fRecoHProton", "fRecoHProton", 5, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRecoHProton->Sumw2();
    fRecoHProton->GetAxis(0)->SetTitle("avg_xy_dist");
    fRecoHProton->GetAxis(1)->SetTitle("avg_z_dist");
    fRecoHProton->GetAxis(2)->SetTitle("dphi");
    fRecoHProton->GetAxis(3)->SetTitle("deta");
    fRecoHProton->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRecoHProton);

    fRecoHProtonMixed = new THnSparseF("fRecoHProtonMixed", "fRecoHProtonMixed", 5, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRecoHProtonMixed->Sumw2();
    fRecoHProtonMixed->GetAxis(0)->SetTitle("avg_xy_dist");
    fRecoHProtonMixed->GetAxis(1)->SetTitle("avg_z_dist");
    fRecoHProtonMixed->GetAxis(2)->SetTitle("dphi");
    fRecoHProtonMixed->GetAxis(3)->SetTitle("deta");
    fRecoHProtonMixed->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRecoHProtonMixed);

    fRealHProton = new THnSparseF("fRealHProton", "fRealHProton", 5, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRealHProton->Sumw2();
    fRealHProton->GetAxis(0)->SetTitle("avg_xy_dist");
    fRealHProton->GetAxis(1)->SetTitle("avg_z_dist");
    fRealHProton->GetAxis(2)->SetTitle("dphi");
    fRealHProton->GetAxis(3)->SetTitle("deta");
    fRealHProton->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRealHProton);

    fRealHProtonMixed = new THnSparseF("fRealHProtonMixed", "fRealHProtonMixed", 5, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fRealHProtonMixed->Sumw2();
    fRealHProtonMixed->GetAxis(0)->SetTitle("avg_xy_dist");
    fRealHProtonMixed->GetAxis(1)->SetTitle("avg_z_dist");
    fRealHProtonMixed->GetAxis(2)->SetTitle("dphi");
    fRealHProtonMixed->GetAxis(3)->SetTitle("deta");
    fRealHProtonMixed->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fRealHProtonMixed);

    fHProtonNotFound = new THnSparseF("fHProtonNotFound", "fHProtonNotFound", 5, hp_cor_bins, hp_cor_mins, hp_cor_maxes);
    fHProtonNotFound->Sumw2();
    fHProtonNotFound->GetAxis(0)->SetTitle("avg_xy_dist");
    fHProtonNotFound->GetAxis(1)->SetTitle("avg_z_dist");
    fHProtonNotFound->GetAxis(2)->SetTitle("dphi");
    fHProtonNotFound->GetAxis(3)->SetTitle("deta");
    fHProtonNotFound->GetAxis(4)->SetTitle("zvtx");
    fOutputList->Add(fHProtonNotFound);





    fMinDistanceAll = new TH1D("fMinDistanceAll", "fMinDistanceAll", 500, 0, 50);
    fOutputList->Add(fMinDistanceAll);

    fMinDistanceNotFound = new TH1D("fMinDistanceNotFound", "fMinDistanceNotFound", 500, 0, 50);
    fOutputList->Add(fMinDistanceNotFound);

    PostData(1, fOutputList);

}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::GetXYZatR(AliMixingParticle *track, float radius, double *xyz) {

    Double_t alpha=0.0;

    double bz = (double)fAOD->GetMagneticField();

    double position[3];

    position[0] = track->Xv();
    position[1] = track->Yv();
    position[2] = track->Zv();

    Double_t radPos2 = position[0]*position[0]+position[1]*position[1];  
    Double_t radMax  = 45.;

    if (radPos2 < radMax*radMax) { 
        alpha = track->Phi(); 
    } else { 
        Float_t phiPos = TMath::Pi()+TMath::ATan2(-position[1], -position[0]);
        alpha = TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
    }

    TVector3 ver(position[0],position[1],position[2]);
    TVector3 mom(track->Px(),track->Py(),track->Pz());

    ver.RotateZ(-alpha);
    mom.RotateZ(-alpha);

    Double_t fx = ver.X();
    Double_t fy = ver.Y();
    Double_t fz = ver.Z();
    Double_t sn = TMath::Sin(mom.Phi());
    Double_t tgl = mom.Pz()/mom.Pt();
    Double_t crv = TMath::Sign(1/mom.Pt(),(Double_t)track->Charge())*bz*kB2C;

    if ( (TMath::Abs(bz))<kAlmost0Field ) crv=0.;
    double tR = 1./crv;
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double phi0 = TMath::ATan2(y0,x0);
    if (tR<0) phi0 += TMath::Pi();
    if      (phi0 > TMath::Pi()) phi0 -= 2.*TMath::Pi();
    else if (phi0 <-TMath::Pi()) phi0 += 2.*TMath::Pi();
    double cs0 = TMath::Cos(phi0);
    double sn0 = TMath::Sin(phi0);
    double r0 = x0*cs0 + y0*sn0 - tR;
    double r2R = 1.+r0/tR;

    if (r2R<kAlmost0) return kFALSE;
    double xr2R = radius/tR;
    double r2Ri = 1./r2R;
    double cosT = 0.5*(r2R + (1-xr2R*xr2R)*r2Ri);
    if ( TMath::Abs(cosT)>kAlmost1 ) {
        return kFALSE;
    }
    double t = TMath::ACos(cosT);
    if (tR<0) t = -t;
    double xyzi[3];
    xyzi[0] = x0 - tR*TMath::Cos(t+phi0);
    xyzi[1] = y0 - tR*TMath::Sin(t+phi0);
    if (xyz) {
        double t0 = TMath::ATan2(cs,-sn) - phi0;
        double z0 = fz - t0*tR*tgl;    
        xyzi[2] = z0 + tR*t*tgl;
    }
    else xyzi[2] = 0;

    Double_t cs_a=TMath::Cos(alpha), sn_a=TMath::Sin(alpha), x_a=xyzi[0];
    xyzi[0]=x_a*cs_a - xyzi[1]*sn_a; xyzi[1]=x_a*sn_a + xyzi[1]*cs_a;

    if (xyz) {
        xyz[0] = xyzi[0];
        xyz[1] = xyzi[1];
        xyz[2] = xyzi[2];
    }

    return true;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::GetXYZatR(AliAODMCParticle *track, float radius, double *xyz) {

    Double_t alpha=0.0;

    double bz = (double)fAOD->GetMagneticField();

    double position[3];

    position[0] = track->Xv();
    position[1] = track->Yv();
    position[2] = track->Zv();

    Double_t radPos2 = position[0]*position[0]+position[1]*position[1];  
    Double_t radMax  = 45.;

    if (radPos2 < radMax*radMax) { 
        alpha = track->Phi(); 
    } else { 
        Float_t phiPos = TMath::Pi()+TMath::ATan2(-position[1], -position[0]);
        alpha = TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
    }

    TVector3 ver(position[0],position[1],position[2]);
    TVector3 mom(track->Px(),track->Py(),track->Pz());

    ver.RotateZ(-alpha);
    mom.RotateZ(-alpha);

    Double_t fx = ver.X();
    Double_t fy = ver.Y();
    Double_t fz = ver.Z();
    Double_t sn = TMath::Sin(mom.Phi());
    Double_t tgl = mom.Pz()/mom.Pt();
    Double_t crv = TMath::Sign(1/mom.Pt(),(Double_t)track->Charge())*bz*kB2C;

    if ( (TMath::Abs(bz))<kAlmost0Field ) crv=0.;
    double tR = 1./crv;
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double phi0 = TMath::ATan2(y0,x0);
    if (tR<0) phi0 += TMath::Pi();
    if      (phi0 > TMath::Pi()) phi0 -= 2.*TMath::Pi();
    else if (phi0 <-TMath::Pi()) phi0 += 2.*TMath::Pi();
    double cs0 = TMath::Cos(phi0);
    double sn0 = TMath::Sin(phi0);
    double r0 = x0*cs0 + y0*sn0 - tR;
    double r2R = 1.+r0/tR;

    if (r2R<kAlmost0) {
        std::cout << "r2R is too small" << std::endl;
        return kFALSE;
    }

    double xr2R = radius/tR;
    double r2Ri = 1./r2R;

    double cosT = 0.5*(r2R + (1-xr2R*xr2R)*r2Ri);
    if ( TMath::Abs(cosT)>kAlmost1 ) {
        std::cout << "cosT is too big" << std::endl;
        return kFALSE;
    }

    double t = TMath::ACos(cosT);
    if (tR<0) t = -t;

    double xyzi[3];
    xyzi[0] = x0 - tR*TMath::Cos(t+phi0);
    xyzi[1] = y0 - tR*TMath::Sin(t+phi0);
    if (xyz) {
        double t0 = TMath::ATan2(cs,-sn) - phi0;
        double z0 = fz - t0*tR*tgl;    
        xyzi[2] = z0 + tR*t*tgl;
    }
    else xyzi[2] = 0;

    Double_t cs_a=TMath::Cos(alpha), sn_a=TMath::Sin(alpha), x_a=xyzi[0];
    xyzi[0]=x_a*cs_a - xyzi[1]*sn_a; xyzi[1]=x_a*sn_a + xyzi[1]*cs_a;

    if (xyz) {
        xyz[0] = xyzi[0];
        xyz[1] = xyzi[1];
        xyz[2] = xyzi[2];
    }

    return true;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::GetXYZatR(AliAODTrack *track, float radius, double *xyz) {

    Double_t alpha=0.0;

    double bz = (double)fAOD->GetMagneticField();

    double position[3];
    track->GetPosition(position);

    Double_t radPos2 = position[0]*position[0]+position[1]*position[1];  
    Double_t radMax  = 45.;

    if (radPos2 < radMax*radMax) { 
        alpha = track->Phi(); 
    } else { 
        Float_t phiPos = TMath::Pi()+TMath::ATan2(-position[1], -position[0]);
        alpha = TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
    }

    TVector3 ver(position[0],position[1],position[2]);
    TVector3 mom(track->Px(),track->Py(),track->Pz());

    ver.RotateZ(-alpha);
    mom.RotateZ(-alpha);

    Double_t fx = ver.X();
    Double_t fy = ver.Y();
    Double_t fz = ver.Z();
    Double_t sn = TMath::Sin(mom.Phi());
    Double_t tgl = mom.Pz()/mom.Pt();
    Double_t crv = TMath::Sign(1/mom.Pt(),(Double_t)track->Charge())*bz*kB2C;

    if ( (TMath::Abs(bz))<kAlmost0Field ) crv=0.;
    double tR = 1./crv;
    double cs = TMath::Sqrt((1-sn)*(1+sn));
    double x0 = fx - sn*tR;
    double y0 = fy + cs*tR;
    double phi0 = TMath::ATan2(y0,x0);
    if (tR<0) phi0 += TMath::Pi();
    if      (phi0 > TMath::Pi()) phi0 -= 2.*TMath::Pi();
    else if (phi0 <-TMath::Pi()) phi0 += 2.*TMath::Pi();
    double cs0 = TMath::Cos(phi0);
    double sn0 = TMath::Sin(phi0);
    double r0 = x0*cs0 + y0*sn0 - tR;
    double r2R = 1.+r0/tR;

    if (r2R<kAlmost0) return kFALSE;
    double xr2R = radius/tR;
    double r2Ri = 1./r2R;
    double cosT = 0.5*(r2R + (1-xr2R*xr2R)*r2Ri);
    if ( TMath::Abs(cosT)>kAlmost1 ) {
        return kFALSE;
    }
    double t = TMath::ACos(cosT);
    if (tR<0) t = -t;
    double xyzi[3];
    xyzi[0] = x0 - tR*TMath::Cos(t+phi0);
    xyzi[1] = y0 - tR*TMath::Sin(t+phi0);
    if (xyz) {
        double t0 = TMath::ATan2(cs,-sn) - phi0;
        double z0 = fz - t0*tR*tgl;    
        xyzi[2] = z0 + tR*t*tgl;
    }
    else xyzi[2] = 0;

    Double_t cs_a=TMath::Cos(alpha), sn_a=TMath::Sin(alpha), x_a=xyzi[0];
    xyzi[0]=x_a*cs_a - xyzi[1]*sn_a; xyzi[1]=x_a*sn_a + xyzi[1]*cs_a;

    if (xyz) {
        xyz[0] = xyzi[0];
        xyz[1] = xyzi[1];
        xyz[2] = xyzi[2];
    }
    
    return true;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeSameHProtonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> proton_list, THnSparse* corr_dist, double zVtx)
{
    double corr_point[5];

    double TPC_INNER = 85;
    double TPC_MID = 167.5;
    double TPC_OUTER = 250;


    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        // for now we put the pt cuts here, maybe at axes later
        if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;

        for(int i = 0; i < (int)proton_list.size(); i++) {
            auto proton = proton_list[i];
            // for now we put the pt cuts here, maybe at axes later
            if(proton->Pt() < 2 || proton->Pt() > 4) continue;

            // make sure current trigger isn't current proton
            if(trigger->GetLabel() == proton->GetLabel()) continue;

            // make sure that the charge is opposite
            if(trigger->Charge() == proton->Charge()) continue;

            double trig_pos[3];
            double proton_pos[3];

            bool goodTrigDist = GetXYZatR(trigger, TPC_INNER, trig_pos);
            bool goodProtonDist = GetXYZatR(proton, TPC_INNER, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_inner = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_inner = TMath::Abs(trig_pos[2] - proton_pos[2]);

            goodTrigDist = GetXYZatR(trigger, TPC_MID, trig_pos);
            goodProtonDist = GetXYZatR(proton, TPC_MID, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_mid = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_mid = TMath::Abs(trig_pos[2] - proton_pos[2]);

            goodTrigDist = GetXYZatR(trigger, TPC_OUTER, trig_pos);
            goodProtonDist = GetXYZatR(proton, TPC_OUTER, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_outer = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_outer = TMath::Abs(trig_pos[2] - proton_pos[2]);


            corr_point[0] = (dR_inner + dR_mid + dR_outer) / 3.;; 
            corr_point[1] = (dZ_inner + dZ_mid + dZ_outer) / 3.;; 

            corr_point[2] = trigger->Phi() - proton->Phi();

            if(corr_point[2] < -TMath::Pi()/2.0) {
                corr_point[2] += 2.0*TMath::Pi();
            }
            else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                corr_point[2] -= 2.0*TMath::Pi();
            }

            corr_point[3] = trigger->Eta() - proton->Eta();

            corr_point[4] = zVtx;

            corr_dist->Fill(corr_point);

        }
    }
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeSameMCHProtonCorrelations(std::vector<AliAODMCParticle*> trigger_list, std::vector<AliAODMCParticle*> proton_list, THnSparse* corr_dist, double zVtx)
{
    double corr_point[5];

    double TPC_INNER = 85;
    double TPC_MID = 167.5;
    double TPC_OUTER = 250;

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        // for now we put the pt cuts here, maybe at axes later
        if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;

        for(int i = 0; i < (int)proton_list.size(); i++) {
            auto proton = proton_list[i];
            // for now we put the pt cuts here, maybe at axes later
            if(proton->Pt() < 2 || proton->Pt() > 4) continue;
            // make sure current trigger isn't current proton
            if(trigger->GetLabel() == proton->GetLabel()) continue;
            // make sure that the charge is opposite
            if(trigger->Charge() == proton->Charge()) continue;

            double trig_pos[3];
            double proton_pos[3];

            bool goodTrigDist = GetXYZatR(trigger, TPC_INNER, trig_pos);
            bool goodProtonDist = GetXYZatR(proton, TPC_INNER, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_inner = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_inner = TMath::Abs(trig_pos[2] - proton_pos[2]);

            goodTrigDist = GetXYZatR(trigger, TPC_MID, trig_pos);
            goodProtonDist = GetXYZatR(proton, TPC_MID, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_mid = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_mid = TMath::Abs(trig_pos[2] - proton_pos[2]);

            goodTrigDist = GetXYZatR(trigger, TPC_OUTER, trig_pos);
            goodProtonDist = GetXYZatR(proton, TPC_OUTER, proton_pos);

            if(!goodTrigDist || !goodProtonDist) continue;

            double dR_outer = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
            double dZ_outer = TMath::Abs(trig_pos[2] - proton_pos[2]);


            corr_point[0] = (dR_inner + dR_mid + dR_outer) / 3.;; 
            corr_point[1] = (dZ_inner + dZ_mid + dZ_outer) / 3.;; 

            corr_point[2] = trigger->Phi() - proton->Phi();

            if(corr_point[2] < -TMath::Pi()/2.0) {
                corr_point[2] += 2.0*TMath::Pi();
            }
            else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                corr_point[2] -= 2.0*TMath::Pi();
            }

            corr_point[3] = trigger->Eta() - proton->Eta();

            corr_point[4] = zVtx;

            corr_dist->Fill(corr_point);

        }
    }
}



void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeMixedHProtonCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> proton_list , THnSparse* corr_dist, double zVtx)
{
    
    double corr_point[5];

    double TPC_INNER = 85;
    double TPC_MID = 167.5;
    double TPC_OUTER = 250;

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {

            AliMixingParticle *trigger = (AliMixingParticle*) tracks->At(i);
            if(!trigger) continue;
            // for now we put the pt cuts here, maybe at axes later
            if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;
            for(int j = 0; j < (int)proton_list.size(); j++) {
                auto proton = proton_list[j];
                // for now we put the pt cuts here, maybe at axes later
                if(proton->Pt() < 2 || proton->Pt() > 4) continue;
                // make sure that the charge is opposite
                if(trigger->Charge() == proton->Charge()) continue;

                double trig_pos[3];
                double proton_pos[3];

                bool goodTrigDist = GetXYZatR(trigger, TPC_INNER, trig_pos);
                bool goodProtonDist = GetXYZatR(proton, TPC_INNER, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_inner = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_inner = TMath::Abs(trig_pos[2] - proton_pos[2]);

                goodTrigDist = GetXYZatR(trigger, TPC_MID, trig_pos);
                goodProtonDist = GetXYZatR(proton, TPC_MID, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_mid = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_mid = TMath::Abs(trig_pos[2] - proton_pos[2]);

                goodTrigDist = GetXYZatR(trigger, TPC_OUTER, trig_pos);
                goodProtonDist = GetXYZatR(proton, TPC_OUTER, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_outer = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_outer = TMath::Abs(trig_pos[2] - proton_pos[2]);


                corr_point[0] = (dR_inner + dR_mid + dR_outer) / 3.;; 
                corr_point[1] = (dZ_inner + dZ_mid + dZ_outer) / 3.;; 

                corr_point[2] = trigger->Phi() - proton->Phi();

                if(corr_point[2] < -TMath::Pi()/2.0) {
                    corr_point[2] += 2.0*TMath::Pi();
                }
                else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                    corr_point[2] -= 2.0*TMath::Pi();
                }

                corr_point[3] = trigger->Eta() - proton->Eta();

                corr_point[4] = zVtx;

                corr_dist->Fill(corr_point);


            }
        }
    }
}
void AliAnalysisTaskLambdaHadronCloseTrackEff::MakeMixedMCHProtonCorrelations(AliEventPool* fPool, std::vector<AliAODMCParticle*> proton_list , THnSparse* corr_dist, double zVtx)
{
    
    double corr_point[5];

    double TPC_INNER = 85;
    double TPC_MID = 167.5;
    double TPC_OUTER = 250;

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliMixingParticle *trigger = (AliMixingParticle*) tracks->At(i);
            if(!trigger) continue;
            // for now we put the pt cuts here, maybe at axes later
            if(trigger->Pt() < 4 || trigger->Pt() > 8) continue;
            for(int j = 0; j < (int)proton_list.size(); j++) {
                auto proton = proton_list[j];
                // for now we put the pt cuts here, maybe at axes later
                if(proton->Pt() < 2 || proton->Pt() > 4) continue;
                // make sure that the charge is opposite
                if(trigger->Charge() == proton->Charge()) continue;

                double trig_pos[3];
                double proton_pos[3];

                bool goodTrigDist = GetXYZatR(trigger, TPC_INNER, trig_pos);
                bool goodProtonDist = GetXYZatR(proton, TPC_INNER, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_inner = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_inner = TMath::Abs(trig_pos[2] - proton_pos[2]);

                goodTrigDist = GetXYZatR(trigger, TPC_MID, trig_pos);
                goodProtonDist = GetXYZatR(proton, TPC_MID, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_mid = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_mid = TMath::Abs(trig_pos[2] - proton_pos[2]);

                goodTrigDist = GetXYZatR(trigger, TPC_OUTER, trig_pos);
                goodProtonDist = GetXYZatR(proton, TPC_OUTER, proton_pos);

                if(!goodTrigDist || !goodProtonDist) continue;

                double dR_outer = TMath::Sqrt(TMath::Power(trig_pos[0] - proton_pos[0], 2) + TMath::Power(trig_pos[1] - proton_pos[1], 2));
                double dZ_outer = TMath::Abs(trig_pos[2] - proton_pos[2]);


                corr_point[0] = (dR_inner + dR_mid + dR_outer) / 3.;; 
                corr_point[1] = (dZ_inner + dZ_mid + dZ_outer) / 3.;; 

                corr_point[2] = trigger->Phi() - proton->Phi();

                if(corr_point[2] < -TMath::Pi()/2.0) {
                    corr_point[2] += 2.0*TMath::Pi();
                }
                else if(corr_point[2] > 3.0*TMath::Pi()/2.0) {
                    corr_point[2] -= 2.0*TMath::Pi();
                }

                corr_point[3] = trigger->Eta() - proton->Eta();

                corr_point[4] = zVtx;

                corr_dist->Fill(corr_point);

            }
        }
    }
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassDaughterCuts(AliAODTrack *track){

    if(track->GetID() < 0) return false;

    bool pass = true;

    pass = pass && (TMath::Abs(track->Eta()) < 0.8);
    pass = pass && (track->Pt() > 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 80);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);
    // pass = pass && track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA);

    int label = track->GetLabel();


    if(label < 0) return false;

    AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(label);

    // we are only concerned with secondary protons
    if(mcpart->IsPhysicalPrimary()) return false;

    return pass;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassTriggerCuts(AliAODTrack *track, bool checkMC){

    if(!(TMath::Abs(track->Eta()) < 0.8)) return false;
    if(!(track->Pt() > 0.15)) return false;
    if(!track->TestBit(fTriggerBit)) return false;

    if(checkMC) {
        int label = track->GetLabel();
        if(label < 0) return false;

        AliAODMCParticle* mcpart = (AliAODMCParticle*)fMCArray->At(label);

        int pdg = mcpart->GetPdgCode();
        if(!IsMCChargedHadron(pdg)) return false;
        if(!mcpart->IsPhysicalPrimary()) return false;
    }

    return true;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::IsMCChargedHadron(int pdg_code) {
    // checks to see if pdg code matches pion, proton, kaon, electron or muon
    if((TMath::Abs(pdg_code) == 321)
        || (TMath::Abs(pdg_code) == 211)
        || (TMath::Abs(pdg_code) == 2212)
        || (TMath::Abs(pdg_code) == 11)
        || (TMath::Abs(pdg_code) == 13)) return true;

    return false;
}

bool AliAnalysisTaskLambdaHadronCloseTrackEff::PassMCTriggerCuts(AliAODMCParticle *mc_particle){


    if(!IsMCChargedHadron(mc_particle->PdgCode())) return false;
    if(!mc_particle->IsPhysicalPrimary()) return false; // for now trigger is physical primary, could change
    if(!(TMath::Abs(mc_particle->Eta()) < 0.8)) return false;
    if(!(mc_particle->Pt() > 0.15)) return false;

    return true;
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::UserExec(Option_t*)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
    }

    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
        AliError("Array of MC particles not found");
        return;
    }

    fpidResponse = fInputHandler->GetPIDResponse();

    //Event cuts
    TString cent_estimator = fCentEstimator;
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < fMultLow || multPercentile > fMultHigh) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();

    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> proton_list;
    
    std::set<int> trigger_index_set;
    std::set<int> proton_index_set;

    std::vector<AliAODMCParticle*> proton_list_mc;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    //MC trigger list used for event mixing
    TObjArray* fMixedMCTrackObjArray = new TObjArray;
    fMixedMCTrackObjArray->SetOwner(kTRUE);

    // RECO SAME EVENT SECTION

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {

    
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;
        
        //Filter for trigger particles
        if(PassTriggerCuts(track, true)) {
            trigger_list.push_back(track);
            trigger_index_set.insert(track->GetLabel());
            double position[3];
            track->GetPosition(position);
            int mc_label = track->GetLabel();

            auto *mixing_track = new AliMixingParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), position);
            fMixedTrackObjArray->Add(mixing_track);
            }

        if(PassDaughterCuts(track)) {
            int mc_label = track->GetLabel();
            if(mc_label >= 0) {
                auto mc_part = (AliAODMCParticle*)fMCArray->At(mc_label);
                if(TMath::Abs(mc_part->GetPdgCode()) == 2212) {
                    // guaranteed secondary protons that pass loose daughter cuts
                    proton_list.push_back(track);
                    proton_list_mc.push_back(mc_part);
                    proton_index_set.insert(track->GetLabel());
                }
            }
        }
    }

    MakeSameHProtonCorrelations(trigger_list, proton_list, fRecoHProton, primZ);

    // MC SAME EVENT SECTION
    std::vector<AliAODMCParticle*> real_trigger_list;
    std::vector<AliAODMCParticle*> real_proton_list;

    std::set<int> real_proton_index_set;
    std::set<int> real_trigger_index_set;

    for(int mc_index = 0; mc_index < fMCArray->GetEntriesFast(); mc_index++) {

        AliAODMCParticle* mc_particle = (AliAODMCParticle*)fMCArray->At(mc_index);

        if(PassMCTriggerCuts(mc_particle)) {
            real_trigger_list.push_back(mc_particle);
            real_trigger_index_set.insert(mc_index);
            double position[3] = {mc_particle->Xv(), mc_particle->Yv(), mc_particle->Zv()};
            auto mixing_mc_particle = new AliMixingParticle(mc_particle->Pt(), mc_particle->Eta(), mc_particle->Phi(), mc_particle->Charge(), position);
            fMixedMCTrackObjArray->Add(mixing_mc_particle);
        }
        if(TMath::Abs(mc_particle->Eta()) < 0.8 && mc_particle->Pt() > 0.15) {
            // again selecting real protons, but not physical primaries
            if(TMath::Abs(mc_particle->GetPdgCode()) == 2212 && !mc_particle->IsPhysicalPrimary())  {
                real_proton_index_set.insert(mc_index);
                real_proton_list.push_back(mc_particle);
            }
        }
    }

    MakeSameMCHProtonCorrelations(real_trigger_list, real_proton_list, fRealHProton, primZ);


    if(proton_list.size() > 0 ) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fCorPool->IsReady()) {
                MakeMixedHProtonCorrelations(fCorPool, proton_list, fRecoHProtonMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }


    if(real_proton_list.size() > 0) {
        AliEventPool *fMCCorPool = fMCCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fMCCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }
        else {
            if(fMCCorPool->IsReady()) {
                MakeMixedMCHProtonCorrelations(fMCCorPool, real_proton_list, fRealHProtonMixed, primZ);
            }
            if(fMixedMCTrackObjArray->GetEntries() > 0) {
                fMCCorPool->UpdatePool(fMixedMCTrackObjArray);
            }
        }
    }



    // find intersection of trigger_index_set and real_trigger_index_set
    std::set<int> trigger_intersection;
    std::set_intersection(trigger_index_set.begin(), trigger_index_set.end(), real_trigger_index_set.begin(), real_trigger_index_set.end(), std::inserter(trigger_intersection, trigger_intersection.begin()));

    // subract intersection from real_trigger_index_set
    std::set<int> triggers_not_found;
    std::set_difference(real_trigger_index_set.begin(), real_trigger_index_set.end(), trigger_intersection.begin(), trigger_intersection.end(), std::inserter(triggers_not_found, triggers_not_found.begin()));

    // find intersection of proton_index_set and real_proton_index_set
    std::set<int> proton_intersection;
    std::set_intersection(proton_index_set.begin(), proton_index_set.end(), real_proton_index_set.begin(), real_proton_index_set.end(), std::inserter(proton_intersection, proton_intersection.begin()));

    // subract intersection from real_proton_index_set
    std::set<int> protons_not_found;
    std::set_difference(real_proton_index_set.begin(), real_proton_index_set.end(), proton_intersection.begin(), proton_intersection.end(), std::inserter(protons_not_found, protons_not_found.begin()));

    // correlate the not_found triggers with not_found protons
    std::vector<AliAODMCParticle*> not_found_trigger_list;
    std::vector<AliAODMCParticle*> not_found_proton_list;

    for(auto trigger_index : triggers_not_found) {
        not_found_trigger_list.push_back((AliAODMCParticle*)fMCArray->At(trigger_index));
    }
    for(auto proton_index : protons_not_found) {
        not_found_proton_list.push_back((AliAODMCParticle*)fMCArray->At(proton_index));
    }

    MakeSameMCHProtonCorrelations(not_found_trigger_list, not_found_proton_list, fHProtonNotFound, primZ);


    PostData(1, fOutputList);
}

void AliAnalysisTaskLambdaHadronCloseTrackEff::Terminate(Option_t *option)
{
}

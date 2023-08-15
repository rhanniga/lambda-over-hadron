#include <cmath>
#include <TMath.h>

#include <AliCDBEntry.h>
#include <AliTPCParamSR.h>

#include <iostream>

void plotBetheBloch2();
void test(){
    cout<< "HELLOOOOO"<< endl;
  
   // void plotBetheBloch2();
 plotBetheBloch2();
}


void plotBetheBloch2(){
    cout<<"MADE IT HERE!"<< endl;
    TFile *file1 = new TFile("AnalysisResults072423.root","READ");
    TList *hadronElec1 = (TList*)file1->Get("h-lambda"); //list with 0-10% centrality results
    TH2D *dedx = (TH2D*)hadronElec1->FindObject("dedx");

   // plotting
    int NUM_BINS = 100;
    double x[NUM_BINS];
    
    for(int i = 0; i < NUM_BINS; i++) {
        double increment=10.0/NUM_BINS;
        x[i] =(i+0.5)*increment;
    }
    //define the four curves
    double ydEdxProton[NUM_BINS],ydEdxPion[NUM_BINS],ydEdxKaon[NUM_BINS],ydEdxElec[NUM_BINS];
    
  
    for(int i = 0; i < NUM_BINS; i++) {
        
        double momentum = x[i];
       
        double BetheBlochFunction( double momentum, double mass);
       
        std::cout << "Momentum: " << momentum << std::endl; 
        //masses in MeV
        double massElec=0.51099895e-3;
        double massProton=938.27208816e-3;
        double massPion= 134.9768e-3;
        double massKaon=493.677e-3;
        
        ydEdxElec[i]=  BetheBlochFunction(momentum, massElec);//returns the functions value at that momentum and mass
        ydEdxProton[i]= BetheBlochFunction(momentum, massProton);//returns the functions value at that momentum and mass
        ydEdxPion[i]=  BetheBlochFunction(momentum, massPion);//returns the functions value at that momentum and mass
        ydEdxKaon[i]=  BetheBlochFunction(momentum, massKaon);//returns the functions value at that momentum and mass
      
    }
    auto grdEdxProton = new TGraph (NUM_BINS, x, ydEdxProton);
    auto grdEdxPion = new TGraph (NUM_BINS, x, ydEdxPion);
    auto grdEdxKaon = new TGraph (NUM_BINS, x, ydEdxKaon);
    auto grdEdxElec= new TGraph (NUM_BINS, x, ydEdxElec);
    
    //Draw the Canvas
    TCanvas *canvasBetheBloch = new TCanvas("canvasBetheBloch","Bethe Bloch ",60,60,900,500);
    canvasBetheBloch->cd()->SetLogz();
    canvasBetheBloch->cd()->SetLogx();
    dedx->GetXaxis()->SetTitle("p (GeV/c)");
    dedx->GetXaxis()->SetRangeUser(0,10);
    dedx->GetYaxis()->SetTitle("dE/dx in TPC");
    dedx->Draw("colz");
    grdEdxProton->Draw("SAME");
    grdEdxProton->SetLineColor(1);
    grdEdxPion->Draw("SAME");
    grdEdxPion->SetLineColor(1);
    grdEdxKaon->Draw("SAME");
    grdEdxKaon->SetLineColor(1);
    grdEdxElec->Draw("SAME");
    grdEdxElec->SetLineColor(1);
    
    //LEGEND
    void ProcessLegend(TLegend *leg, Int_t font);
    
    auto legenddEdx = new TLegend(.6,0.7,.85,0.9);
    legenddEdx->AddEntry((TObject*)0, "ALICE     0-100%", "");
    legenddEdx->AddEntry((TObject*)0, "p-Pb, #sqrt{s_{NN}}= 5.02 TeV", "");
    
    Int_t font=42;
    ProcessLegend(legenddEdx, font);
   
    auto legendElec = new TLegend(0.0001,0.5,.85,0.9);
    legendElec->AddEntry((TObject*)0, "e", "");
    ProcessLegend(legendElec, font);
    
    auto legendPion= new TLegend(0.0001,0.1,.85,0.9);
    legendPion->AddEntry((TObject*)0, "#it{#pi}", "");
    ProcessLegend(legendPion, font);
    
    auto legendKaon= new TLegend(0.05,0.85,.85,0.9);
    legendKaon->AddEntry((TObject*)0, "K", "");
    ProcessLegend(legendKaon, font);
    
    auto legendProton= new TLegend(0.25,0.85,.85,0.9);
    legendProton->AddEntry((TObject*)0, "p", "");
    ProcessLegend(legendProton, font);
    
    legenddEdx->Draw("SAME");
    legendElec->Draw("SAME");
    legendPion->Draw("SAME");
    legendKaon->Draw("SAME");
    

    
    
}
double BetheBlochFunction(double momentum, double mass){
    TFile *file = new TFile("Run0_999999999_v0_s3.root","READ");

    AliCDBEntry *entry = (AliCDBEntry*)file->Get("AliCDBEntry");

    AliTPCParamSR *tpcParam = (AliTPCParamSR*)entry->GetObject();

    TVectorD *betheBloch = tpcParam->GetBetheBlochParamAlice();

    for(int i = 0; i < betheBloch->GetNoElements(); i++) {
   //  std::cout << i << " " << betheBloch->operator[](i) << std::endl;
    }
    
    double p1=betheBloch->operator[](0);
    double p2=betheBloch->operator[](1);
    double p3=betheBloch->operator[](2);
    double p4=betheBloch->operator[](3);
    double p5=betheBloch->operator[](4);
    
    
    
   double Energy= TMath::Sqrt(momentum*momentum+mass*mass);
    double gamma= Energy/mass;
    double beta = TMath::Sqrt(1 - (1/(gamma*gamma)));

    
    double function=(p1/TMath::Power(beta,p4))*(p2-TMath::Power(beta,p4)-TMath::Log(p3+(1/TMath::Power(beta*gamma,p5))));
   // cout<< "Value of function: "<< function<< endl;
    return 55*function;
}


void ProcessLegend(TLegend *leg, Int_t font)
    {
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->SetTextFont(font);
        leg->SetFillStyle(0);
        leg->Draw();
   }

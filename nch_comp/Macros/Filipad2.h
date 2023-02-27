#include <iostream>
//bk
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
//endbk
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLatex.h"

class Filipad2 {

    public:

        // ---- Costructor ------------
        Filipad2(int inID=1, float inRelSize=1.1, float inR = 0.4, int inXOffset = 100, int inYOffset=100, float inAspect=0.7, int ichop=5, int ichopPt=3){
            aspektCanvas = inAspect;
            sizeCanvas   = 300*inRelSize;

            MarginLeft   = 0.15;
            MarginBottom = 0.08; 
            MarginRight  = 0.03;
            MarginTop    = 0.02;
            ID           = inID+ichopPt*100;
            mcpad        = Form("c%d",inID);
            ratio = inR;

            //sdxCanvas    = inXOffset;
            int inID0 = inID-1;
            //sdxCanvas    = inXOffset*(inID0%ichop)+10;
            //sdyCanvas    = inYOffset*(inID0-inID0%ichop)/ichop+10;
            //sdyCanvas    = inYOffset*(inID0-inID0%ichop)/ichop+10;
	    sdxCanvas = inXOffset*((inID-ichopPt)/10.)+10;
	    sdyCanvas = inYOffset*(ichopPt)+10;
	
            //cout <<"sdxCanvas= "<<  sdxCanvas <<" sdyCanvas="<< sdyCanvas <<endl; 
            space        = 0;
        }

        // ---- Destructor ------------
        ~Filipad2(){
            cout<<"Destructor"<<endl;
            if(C)      delete C;
            if(toppad) delete toppad;
        }


        TPad* GetPad(int padID){ return (TPad*) toppad->cd(padID);}//coordinates 0,0 = upper left 

        //-------------------------------------------------------------------------------------
        void Draw(){ 

            char name[200];
            // cout<<"Draw"<<endl;
            C = new TCanvas(mcpad, mcpad, sdxCanvas, sdyCanvas, sizeCanvas*aspektCanvas, sizeCanvas);//the main canvas
						C -> SetName(Form("%s",mcpad.Data()));
            C->SetFillStyle(4000); C->SetFillColor(10);
            gStyle->SetOptStat(0);    gStyle->SetOptTitle(0);
            C->SetTopMargin(0.); C->SetBottomMargin(0.);//our pads will have no margin
            C->SetLeftMargin(0.);C->SetRightMargin(0.);
            //C->cd();
            C->Draw();

            toppad = gPad;
            toppad->Clear();

            TPad   *pp      = NULL;
            pp = new TPad(Form("UPad%d",ID),Form("UPad%d",ID), 0, 0 + ratio + space  , 1, 1, 0); //create pad
        
            pp->SetNumber(1);   //assign a number to it. Possible to access it via :  toppad->cd(ih);  
            pp->SetTopMargin(MarginTop/(1-ratio)); pp->SetBottomMargin(0.0015);//our pads will have no margin
            pp->SetLeftMargin(MarginLeft);pp->SetRightMargin(MarginRight);
            pp->Draw();

            pp = new TPad(Form("LPad%d",ID),Form("LPad%d",ID),0, 0, 1, ratio ,0); //create pad
            pp->SetNumber(2);   //assign a number to it. Possible to access it via :  toppad->cd(ih); 
            pp->SetTopMargin(0.0015); pp->SetBottomMargin(MarginBottom/ratio);//our pads will have no margin
            pp->SetLeftMargin(MarginLeft);pp->SetRightMargin(MarginRight);
            pp->Draw();
        }


        //-------------------------------------------------------------------------------------
        void Hset(TH1* hid, TString xtit="", TString ytit="",
                double titoffx = 2.5, double titoffy = 1.5,
                double titsizex = 20, double titsizey = 20,
                double labeloffx = 0.01, double labeloffy = 0.001,
                double labelsizex = 16, double labelsizey = 16,
                int divx = 505, int divy=505){

            hid->GetXaxis()->CenterTitle(1);
            hid->GetYaxis()->CenterTitle(1);

            hid->GetXaxis()->SetTitleOffset(titoffx);
            hid->GetYaxis()->SetTitleOffset(titoffy);

            hid->GetXaxis()->SetTitleFont(43);
            hid->GetYaxis()->SetTitleFont(43);
            hid->GetXaxis()->SetTitleSize(titsizex);
            hid->GetYaxis()->SetTitleSize(titsizey);

            hid->GetXaxis()->SetLabelOffset(labeloffx);
            hid->GetYaxis()->SetLabelOffset(labeloffy);

            hid->GetXaxis()->SetLabelFont(43);
            hid->GetYaxis()->SetLabelFont(43);
            hid->GetXaxis()->SetLabelSize(labelsizex);
            hid->GetYaxis()->SetLabelSize(labelsizey);

            hid->GetXaxis()->SetNdivisions(divx);
            hid->GetYaxis()->SetNdivisions(divy);

            hid->GetXaxis()->SetTitle(xtit);
            hid->GetYaxis()->SetTitle(ytit);
        }
        

        void SetMarginLeft(float x){ MarginLeft = x;}
        void SetMarginRight(float x){ MarginRight = x;}
        void SetMarginTop(float x){ MarginTop = x;}
        void SetMarginBottom(float x){ MarginBottom = x;}



        //  M A I N     C A N V A S 
        int   ID;
        float aspektCanvas;      //size and positioning of the main canvas   
        int   sizeCanvas;
        int   sdxCanvas; 
        int   sdyCanvas;

        float MarginLeft; //Margins around the latice of nx times ny pads 
        float MarginBottom;
        float MarginRight;
        float MarginTop;

        float space;
        float ratio;
        TString mcpad;
        TCanvas *C;
        TVirtualPad *toppad;
};



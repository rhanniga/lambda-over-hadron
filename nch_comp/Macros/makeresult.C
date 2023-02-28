#include "BSHelper.cxx"
#include "dndetaCommon.h"
#include "Filipad2.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
//#include <boost/progress.hpp>

const float inf = 1e20;
//TString dataname="MulLHC15nHighMult";
//TString dataname="MulLHC16l";


//TString dataname = "MulLHC17p";
//TString p6name = "MulLHC17pMCPYTHIA8";
//TString p8name = "MulLHC17pMCPYTHIA8";

TString dataname = "data";
TString p6name = "MC";
TString p8name = "MC";

/*
TString dataname = "PtCutLHC16lHM";
TString p6name = "PtCutLHC16lHMMCPYTHIA8";
TString p8name = "PtCutLHC16lHMMCPYTHIA8";
*/
/*
TString dataname = "PtCutLHC16lHMTest";
TString p6name = "PtCutLHC16lHMMCPYTHIA8Test";
TString p8name = "PtCutLHC16lHMMCPYTHIA8Test";
*/
/*
TString dataname = "PtCutLHC16l";
TString p6name = "PtCutLHC16lMCPYTHIA8";
TString p8name = "PtCutLHC16lMCPYTHIA8";
*/
/*
TString dataname = "PtCutLHC16l_Full_LessBins";
TString p6name = "PtCutLHC16lMCPYTHIA8_Full_LessBins";
TString p8name = "PtCutLHC16lMCPYTHIA8_Full_LessBins";
*/
//TString dataname = "PtCutLHC16l";
//TString p6name = "PtCutLHC16lMCEPOS";
//TString p8name = "PtCutLHC16lMCEPOS";

//TString dataname = "MulLHC17p";
//TString p6name = "MulLHC17pMCPYTHIA8";
//TString p8name = "MulLHC17pMCPYTHIA8";
//TString p8name = "MulLHC16lMCEPOS";

//TString dataname = "MulLHC17p";
//TString p6name   = "MulLHC17pMCPYTHIA8";
//TString p8name   = "MulLHC17pMCPYTHIA8";
//TString p8name = "MulLHC16lMCEPOSHighMult";
//TString p6name= "MulLHC16lMCPYTHIA8HighMult";
//TString p8name = "MulLHC16lMCPYTHIA8HighMult";

//TString dataname = "MulLHC15f";
//TString p6name = "MulLHC15fMCPYTHIA8";
//TString p8name = "MulLHC15fMCPYTHIA8";

//TString dataname = "MulLHC15n";
//TString p6name = "MulLHC15nMCPYTHIA8";
//TString p8name = "MulLHC15nMCPYTHIA8";

//TString dataname = "MulLHC10d";
//TString p6name = "MulLHC10dMCPYTHIA6";
//TString p8name = "MulLHC10dMCPYTHIA6";

enum
{
	kECbegin = 1,
	kDATA = 1,
	kINEL,
	kINELg0,
	kINEL300,
	kINEL100,
	kINEL200,
	kECend
};
enum
{
  	kTrigbegin = 1,
	kHMMBAND = 1,
	kHMMBAND300,
	kMBAND015,
	kMBAND300,
	kMBAND100,
	kMBAND200,
	kTrigend
};
enum
{
	kParTypebegin = 1,
	kParDATA = 1,
	kMotherStrange,
	kBkg,
	kPion,
	kKaon,
	kProtonBK,
	kOPar,
	kParTypeend
};
enum
{
	kNoStrVar = 1,
	kStrVarUp,
	kStrVarDw
};


    enum
    {
        kNoTrkCutVar = 1,
        kHybrid = 1 ,
        kITSTPC2011,
        kITSTPC2011dcazdw, 
        kITSTPC2011dcazup,
        kITSTPC2011dcardw, 
        kITSTPC2011dcarup,
        kITSTPC2011nclutpcdw, 
        kITSTPC2011nclutpcup,
        kITSTPC2011chitpcdw, 
        kITSTPC2011chitpcup,
        kITSTPC2011globalconsdw, 
        kITSTPC2011globalconsup,
        kTrackCutend
    };

enum
{
	kCentClassBinBegin = 1,
	kV0M = 1,
	kV0A,
	kV0C,
	kSPDMult,
	kRefMult05,
	kRefMult08,
	kCentClassBinEnd
};

Double1D vtxbin = {-10, 10};
Double1D etabin = {-0.8, 0.8};

//--> For comparison with twiki values (MB)
//Int_t eventclass = kINEL;
//Int_t triggtype = kMBAND015;

//--> For comparison with twiki values (HM)
//Int_t eventclass = kINEL;
//Int_t triggtype = kHMMBAND;

//--> To caompute dndeta for my analysis (MB)
Int_t eventclass = kINEL300;
Int_t triggtype = kMBAND300;

//--> To caompute dndeta for my analysis (MB)
//Int_t eventclass = kINEL300;
//Int_t triggtype = kHMMBAND300;

Int_t centtype = kV0M;
Double1D cent = {0,20}; 
Double1D mccent =  cent;
//Double1D  mccent={1,5};
//Double1D mccent =  { (cent[0]==0 && cent[1] == 0.01 ) ? 0.001 : cent[0], cent[1]<=1 ? 1 : cent[1]};

Bool_t DoSystematic = true;
Bool_t SaveHistogramsAndResults = true;
Bool_t SaveResultsOnly = true;
Bool_t MakeSymm = false;

//Double1D  cent={0,0.001,0.01,0.1,0.5,1,5,10,15,20,30,40,50,70,100};
// 7 TeV 0-0.01 = 0 - 0.5, 0.01 - 0.1 = 0 - 0.5, 50-70 = 40-50
Int_t iptvar = kITSTPC2011; // no ptunseen variation
Int1D colors = {kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2, kPink + 10, kRed, kBlue, kGreen + 3, kMagenta + 2, kPink + 10, kBlack, kRed, kBlue};
Int1D markers = {1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 24, 24, 24, 24, 24};
Int1D lstyle = {1, 1, 1, 1, 1, 9, 9, 9, 9, 9};
Bool_t DoDiffTune = false;

int n = 0;
const char *mstring = "d#it{N}_{ch}/d#it{#eta}";
const char *avmstring = "#LT d#it{N}_{ch}/d#it{#eta} #GT";
const char *avmstringincl = "#LT d#it{N}_{ch}/d#it{#eta} #GT / #LT d#it{N}_{ch}/d#it{#eta} #GT_{Inclusive}";
const char *mstringincl = "(d#it{N}_{ch}/d#it{#eta})/(d#it{N}_{ch}/d#it{#eta})_{Inclusive}";

Double_t rkaon = 1;
Double_t rproton = 1;
Double_t ropar = 1;
Int_t phirange = -1;

// When wants to save hists in ../Figures and final results in ../FinalResults
// During macro test, set to false!

Int_t modelcounter = 0;
Int1D modelcolors = {2, 4, 5, 6, 7};

TGraphAsymmErrors *trackcutsystematic = nullptr;
TGraphAsymmErrors *diffractionratiosystematic = nullptr;
TGraphAsymmErrors *diffractionsystematic = nullptr;
TGraphAsymmErrors *particleratiosystematic = nullptr;
TGraphAsymmErrors *stranegeparticlesystematic = nullptr;
TGraphAsymmErrors *vertexsystematic = nullptr;
TGraphAsymmErrors *centralityestimatorsystematic = nullptr;
TGraphAsymmErrors *phicutsystematic = nullptr;
TGraphAsymmErrors *triggerbiassystematic = nullptr;
TGraphAsymmErrors *symsystematic = nullptr;

void AddPoint(TGraphErrors *graph, Float_t x, Float_t y, Float_t xe, Float_t ye)
{
	graph->SetPoint(graph->GetN(), x, y);
	graph->SetPointError(graph->GetN() - 1, xe, ye);
}

void setpad2(TVirtualPad *pad)
{
	pad->SetTopMargin(0.02);
	pad->SetLeftMargin(0.13);
	pad->SetRightMargin(0.02);
	pad->SetBottomMargin(0.15);
	pad->SetName(Form("c%d", ++n));
}

void ReverseXAxis(TGraph *g)
{
	// Remove the current axis
	g->GetXaxis()->SetLabelOffset(999);
	g->GetXaxis()->SetTickLength(0);

	// Redraw the new axis
	gPad->Update();
	TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
								 gPad->GetUymin(),
								 gPad->GetUxmin(),
								 gPad->GetUymin(),
								 g->GetXaxis()->GetXmin(),
								 g->GetXaxis()->GetXmax(),
								 510, "-SDH");
	newaxis->SetLabelOffset(-0.03);
	newaxis->Draw();
}

void ReverseXGraph(TGraph *g)
{
	// Create a new graph
	Int_t n = g->GetN();
	Double_t *x = g->GetX();
	Double_t *y = g->GetY();
	Double_t xr[100];
	Double_t dx = g->GetXaxis()->GetXmin() + g->GetXaxis()->GetXmax();
	for (Int_t i = 0; i < n; i++)
	{
		xr[i] = -x[i] + dx;
	}

	TGraph *gr = new TGraph(n, xr, y);
	gr->SetMarkerStyle(20);
	gr->SetLineColor(kRed);
	gr->SetMarkerColor(kRed);
	gr->Draw("PL");
}

void setpad(TVirtualPad *pad)
{
	pad->SetTopMargin(0.02);
	pad->SetLeftMargin(0.13);
	pad->SetRightMargin(0.2);
	pad->SetBottomMargin(0.15);
	pad->SetName(Form("c%d", ++n));
}
void hset(TH1 &hid, TString xtit = "", TString ytit = "",
		  double titoffx = 0.9, double titoffy = 1.2,
		  double titsizex = 0.06, double titsizey = 0.06,
		  double labeloffx = 0.01, double labeloffy = 0.001,
		  double labelsizex = 0.05, double labelsizey = 0.05,
		  int divx = 510, int divy = 510)
{
	hid.SetStats(0);

	hid.GetXaxis()->CenterTitle(1);
	hid.GetYaxis()->CenterTitle(1);

	hid.GetXaxis()->SetTitleOffset(titoffx);
	hid.GetYaxis()->SetTitleOffset(titoffy);

	hid.GetXaxis()->SetTitleSize(titsizex);
	hid.GetYaxis()->SetTitleSize(titsizey);

	hid.GetXaxis()->SetLabelOffset(labeloffx);
	hid.GetYaxis()->SetLabelOffset(labeloffy);

	hid.GetXaxis()->SetLabelSize(labelsizex);
	hid.GetYaxis()->SetLabelSize(labelsizey);

	hid.GetXaxis()->SetNdivisions(divx);
	hid.GetYaxis()->SetNdivisions(divy);

	hid.GetXaxis()->SetTitle(xtit);
	hid.GetYaxis()->SetTitle(ytit);
}
TH1D *BaseCorrection(TString mcname, double SDR, double DDR, bool dopariclevariation, Int_t strangevar);
TH1D *MakeModelResults(TString mcname, Bool_t makesymmetric);
TGraphAsymmErrors *CalcRatio(TGraphAsymmErrors *num, TGraphAsymmErrors *den);
TGraphAsymmErrors *CCCC(TGraphAsymmErrors *num, TGraphAsymmErrors *den, double temp);
TGraphAsymmErrors *CalcRatio(TGraphAsymmErrors *num, TH1 *den);

TGraphAsymmErrors *DiffractionSystematic()
{
	double SDR = 0.20;
	double DDR = 0.12;

	const int nvar = 2;
	TH1D *diff[nvar];
	TH1D *ratios[nvar];

	vector<TString> mcnames = {"MulLHC15nMBMCPYTHIA6DiffTune" , "MulLHC15nMBMCPYTHIA6"};
	vector<TString> legendnames = {"Diffraction tune", "Default"};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	TString p6defaultname = p6name;
	for (int i = 0; i < nvar; i++)
	{
		cout << "\n Current status : diffraction ratio variation ongoing" << endl;
		p6name = mcnames[i];
		diff[i] = BaseCorrection(p6name.Data(), SDR, DDR, false, kNoStrVar);
		hset(*diff[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
		diff[i]->SetLineColor(colors[i]);
		diff[i]->SetLineStyle(lstyle[i]);
		diff[i]->Draw("histsame");
		leg->AddEntry(diff[i],legendnames[i], "l");
		ratios[i] = (TH1D *)diff[i]->Clone();
		ratios[i]->Divide(diff[0]);
	}
	leg->Draw();
	p6name = p6defaultname;

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.92);
		ratios[i]->SetMaximum(1.06);
		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *temp = new TH1D("temp", "", 40, 0.8, 1.2);
	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(ratios[1]);
	TF1 *f = new TF1("f", "[0]", -2, 2);
	Error->Fit("f", "QN");

	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		y = f->Eval(x);
		Error->SetPointEYhigh(k, y > 1 ? y - 1 : 0);
		Error->SetPointEYlow(k, y<1 ? 1-y : 0);
	}

	delete temp;
	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_diffractionsys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}


TGraphAsymmErrors *DiffractionRatioSystematic()
{
	double SDR = 0.20;
	double DDR = 0.12;
	double sdscale[] = {SDR, SDR * 1.3, SDR * 1.15, SDR * 1, SDR * 0.85,
						SDR * 0.7, SDR * 0.85, SDR * 1, SDR * 1.15};
	double ddscale[] = {DDR, DDR * 1, DDR * 1.15, DDR * 1.3, DDR * 1.15,
						DDR * 1, DDR * 0.85, DDR * 0.7, DDR * 0.85};

	const int nvar = (sizeof(sdscale) / sizeof(sdscale[0]));
	TH1D *diff[nvar];
	TH1D *ratios[nvar];

	TString mcname;
	mcname = p6name.Data();

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);
	for (int i = 0; i < nvar; i++)
	{
		cout << "\n Current status : diffraction ratio variation ongoing" << endl;
		diff[i] = BaseCorrection(mcname.Data(), sdscale[i], ddscale[i], false, kNoStrVar);
		hset(*diff[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
		diff[i]->SetLineColor(colors[i]);
		diff[i]->SetLineStyle(lstyle[i]);
		diff[i]->Draw("histsame");
		leg->AddEntry(diff[i], Form("SD : %3.2f DD: %3.2f", sdscale[i], ddscale[i]), "l");
		ratios[i] = (TH1D *)diff[i]->Clone();
		ratios[i]->Divide(diff[0]);
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.9);
		ratios[i]->SetMaximum(1.1);
		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *temp = new TH1D("temp", "", 40, 0.8, 1.2);
	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		temp->Reset();
		for (int j = 0; j < nvar; j++)
		{
			temp->Fill(ratios[j]->GetBinContent(i));
		}
		error->SetBinError(i, temp->GetRMS());
	}

	delete temp;
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_diffratiosys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}

TGraphAsymmErrors *TriggerBiasSystematic()
{
	TString defaultdata = dataname;
	//Int_t defaulttriggerclass = triggtype;
	//Int1D triggerclasses = {defaulttriggerclass, kMBORg0, defaulttriggerclass};
	vector<TString> mcnames;
	vector<TString> datanames;
	vector<TString> tits;

	if (dataname.Contains("15n"))
	{
		mcnames = {p6name.Data(), p6name.Data()};
		datanames = {dataname.Data(), "MulLHC15nHighMu"};
		tits = {"Total", "#mu ~ 0.05"};
	}
	else if (dataname.Contains("16l"))
	{
		mcnames = {p6name.Data(), p6name.Data()};
		datanames = {dataname.Data(), "MulLHC16lHighMultHighMu"};
		tits = {"Total", "#mu ~ 0.02"};
	}
	else if (dataname.Contains("10d"))
	{
		mcnames = {p6name.Data(), p6name.Data()};
		datanames = {dataname.Data(), "MulLHC10dHighMu"};
		tits = {"Total", "#mu ~ 0.1"};
	}
	else if (dataname.Contains("17p"))
	{
		mcnames = {p6name.Data(), p6name.Data()};
		datanames = {dataname.Data(), "MulLHC17pHighMu"};
		tits = {"Total", "#mu ~ 0.02"};
	}

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (datanames.size());
	TH1D *hists[nvar];
	TH1D *ratios[nvar];

	cout << "\nCurrent status : trigger variation ongoing" << endl;
	//boost::progress_display show_progress( nvar );
	for (int i = 0; i < nvar; i++)
	{
		dataname = datanames[i];
		TString mcname = mcnames[i];
		hists[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, kNoStrVar);
		hset(*hists[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
		hists[i]->SetMaximum(hists[i]->GetBinContent(hists[i]->GetXaxis()->FindBin(0.05)) * 1.5);
		hists[i]->SetMinimum(hists[i]->GetBinContent(hists[i]->GetXaxis()->FindBin(0.05)) * 0.9);
		hists[i]->SetLineColor(colors[i]);
		hists[i]->SetLineStyle(lstyle[i]);
		hists[i]->Draw("histsame][");
		leg->AddEntry(hists[i], tits[i].Data(), "l");
		ratios[i] = (TH1D *)hists[i]->Clone();
		ratios[i]->Divide(hists[0]);
		//++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.85);
		ratios[i]->SetMaximum(1.15);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame][");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	TH1D *errorh = (TH1D *)error->Clone();
	TH1D *errorl = (TH1D *)error->Clone();
	Double_t MAX = 0, MIN = 0;
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double maxe = 0, mine = 0;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
		if (bincenter<-0.5) continue;
		if (bincenter>0.5) continue;

		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			double vale = ratios[j]->GetBinError(i);

			if (val > max)
			{
				max = val;
				maxe = vale;
			}
			if (val < min)
			{
				min = val;
				mine = vale;
			}
			if (abs(bincenter) < 1.5)
			{
				if (val > MAX)
				{
					MAX = val;
				}
				if (val < MIN)
				{
					MIN = val;
				}
			}
		}

		//if (abs(max-1)<abs(1-min)) max = (1-min)+1;
		max = (1 - min) + 1;
		maxe = mine;
		errorh->SetBinContent(i, max);
		errorh->SetBinError(i, maxe);
		errorl->SetBinContent(i, 2 - max);
		errorl->SetBinError(i, maxe);
	}

	if (abs(MAX - 1) < abs(1 - MIN))
		MAX = (1 - MIN) + 1;
	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		//Error->SetPointEYhigh(k,0.04);
		//Error->SetPointEYlow(k,0.04);
	}

	TF1 *fh = new TF1("fh", "[0]", -2, 2);
	TF1 *fl = new TF1("fl", "[0]", -2, 2);
	errorh->Fit("fh", "QN");
	errorl->Fit("fl", "QN");
	fh->SetLineWidth(3);
	fl->SetLineWidth(3);
	fh->Draw("same");
	fl->Draw("same");

	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		//Error->SetPointEYhigh(k,fh->Eval(x)-1);
		//Error->SetPointEYlow(k,1-fl->Eval(x));
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	cout << "Trigger class and data name set to defaults" << endl;
	dataname = defaultdata;
	//	triggtype = defaulttriggerclass;

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_triggerbias.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));
	return Error;
}




void makeCharacteristics()
{

	TString mcname = p6name.Data();
	auto clistmc = LoaddndetaResultList(mcname.Data(), "output");
	auto clistdata = LoaddndetaResultList(dataname.Data(), "output");
	auto hnmc = (TH1D *)clistmc->FindObject("hEventNumbers");
	Double_t nmc = hnmc->GetBinContent(5);
	auto hndata = (TH1D *)clistdata->FindObject("hEventNumbers");
	Double_t ndata = hndata->GetBinContent(5);

	cout << nmc << endl;
	cout << ndata << endl;

	auto hz = (TH1D *)clistdata->FindObject("zdata");
	auto hzmc = (TH1D *)clistmc->FindObject("zdata");
	auto hztrue = (TH1D *)clistmc->FindObject("ztruth");

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(1);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	hz->Scale(1. / ndata, "width");
	hz->SetLineColor(1);
	hz->SetMarkerColor(1);
	hz->SetMarkerStyle(7);
	hzmc->Scale(1. / nmc, "width");
	hzmc->SetLineColor(2);
	hzmc->SetMarkerColor(2);
	hzmc->SetMarkerStyle(24);
	hztrue->Scale(1. / nmc, "width");
	hztrue->SetLineColor(4);
	hztrue->SetMarkerColor(4);
	hztrue->SetMarkerStyle(8);
	hz->SetMaximum(2);
	hz->SetMinimum(1e-7);
	hset(*hz, "#it{z}_{vtx}", "Normalised entries", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
	hz->Draw("PZsame");
	hzmc->Draw("PZsame");
	hztrue->Draw("PZsame");

	TLegend *leg = new TLegend(0.150718, 0.0817391, 0.827751, 0.332174, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);
	leg->AddEntry(hz, "DATA");
	leg->AddEntry(hzmc, "MC-REC");
	leg->AddEntry(hztrue, "MC-TRUTH");
	leg->Draw();

	auto zratio = (TH1D *)hz->Clone();
	zratio->Divide(hzmc);

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetGridy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	zratio->SetMinimum(0.1);
	zratio->SetMaximum(10);
	hset(*zratio, "#it{z}_{vtx} (cm)", "Ratio (DATA/MC-REC)", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	zratio->Draw("PZsame");

	Filipad2 *canvas2 = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas2->Draw();
	TPad *p2 = canvas2->GetPad(1); //upper pad
	p2->SetRightMargin(0.1);
	p2->SetTickx();
	p2->SetGridy(0);
	p2->SetGridx(1);
	p2->SetGridy(1);
	p2->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto phieta = (TH2D *)clistdata->FindObject("hPhiEta");
	phieta->Scale(1. / phieta->Integral(), "width");
	hset(*phieta, "#it{#phi}", "#it{#eta}, DATA", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	phieta->Draw("colz");

	p2 = canvas2->GetPad(2); //lower pad
	p2->SetTickx();
	p2->SetGridy(1);
	p2->SetGridx(1);
	p2->SetGridy(0);
	p2->cd();
	p2->SetRightMargin(0.1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto phieta2 = (TH2D *)clistmc->FindObject("hPhiEta");
	phieta2->Scale(1. / phieta2->Integral(), "width");
	hset(*phieta2, "#it{#phi}", "#it{#eta}, MC-REC", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	phieta2->Draw("colz");

	Filipad2 *canvas3 = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas3->Draw();
	TPad *p3 = canvas3->GetPad(1); //upper pad
	p3->SetTickx();
	p3->SetGridy(0);
	p3->SetGridx(1);
	p3->SetGridy(1);
	p3->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	auto phi = phieta->ProjectionX("_data", 0, -1, "");
	auto phi2 = phieta2->ProjectionX("_mc", 0, -1, "");
	hset(*phi, "#it{#phi}", "Norm entries", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	phi->Draw("histsame");
	phi->SetLineWidth(2);
	phi->SetLineColor(1);

	phi2->SetLineColor(2);
	phi2->SetMarkerStyle(7);
	phi2->SetMarkerColor(2);
	phi2->Draw("pzsame");

	leg = new TLegend(0.167464, 0.0226087, 0.772727, 0.186087, NULL, "brNDC");
	leg->SetTextSize(0.045);
	leg->SetBorderSize(0);
	leg->AddEntry(phi, "DATA");
	leg->AddEntry(phi2, "MC-REC");
	leg->Draw();

	auto phiratio = (TH1D *)phi->Clone();
	phiratio->Divide(phi2);

	p3 = canvas3->GetPad(2); //lower pad
	p3->SetTickx();
	p3->SetGridy(1);
	p3->SetGridx(1);
	p3->SetLogy(0);
	p3->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	hset(*phiratio, "#it{#phi}", "Ratio (DATA/MC-REC)", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	phiratio->Draw("histzsame");

	new TCanvas;
	setpad2(gPad);
	auto cent = (TH1D *)clistdata->FindObject("hcentV0Mzcut20");
	hset(*cent, "Multiplicity percentile (%)", "Entries, DATA", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	cent->SetLineWidth(2);
	cent->SetLineColor(1);
	cent->SetMinimum(0);
	cent->Draw("histsame");
	/*
	canvas3 = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1 );
	canvas3->Draw();
 	p3 = canvas3->GetPad(1); //upper pad
	p3->SetTickx(); p3->SetGridy(0); p3->SetGridx(1); p3->SetGridy(1); p3->cd();
	gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
	
	auto cent = (TH1D*) clistdata -> FindObject("hcentV0Mzcut10");
	auto cent2 = (TH1D*) clistmc -> FindObject("hcentV0Mzcut10");
	hset(*cent, "centrality", "Entries, DATA" ,0.9,0.9, 0.08,0.08, 0.01,0.01, 0.07,0.07, 515,505);
	cent -> SetLineWidth(2);
	cent -> SetLineColor(1);
	cent -> SetMinimum(0);
	cent -> Draw("histsame");




	p3 = canvas3->GetPad(2); //lower pad
	p3->SetTickx(); p3->SetGridy(1); p3->SetGridx(1); p3->SetLogy(0); p3->cd();
	gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
	hset(*cent2, "centrality (V0M)", "Entries, MC-REC" ,0.9,0.9, 0.08,0.08, 0.01,0.01, 0.07,0.07, 515,505);

	cent2 -> SetLineWidth(2);
	cent2 -> SetLineColor(1);
	cent2 -> Draw("histsame");
*/

	TFile::Open("characteristics.root", "recreate");
	hz->SetName("data_z");
	hzmc->SetName("mcrec_z");
	hztrue->SetName("mctrue_Z");
	hz->Write();
	hzmc->Write();
	hztrue->Write();

	phieta->SetName("data_phieta");
	phieta2->SetName("mcrec_phieta");
	phi->SetName("data_phi");
	phi2->SetName("mcrec_phi");
	cent->SetName("data_cent_V0M");
	//cent2 -> SetName("mcrec_cent_V0M");

	phieta->Write();
	phieta2->Write();
	phi->Write();
	phi2->Write();
	cent->Write();
	//cent2 -> Write();
}
TGraphAsymmErrors *VertexSystematic()
{

	TString mcname;
	mcname = p6name.Data();

	Double1D defaultvtxbin = vtxbin;
	Double1D vtxvarl = {-10, -7, -4, -7, -4, -4, -10 , -7, -10,-15};
	Double1D vtxvarr = {10,   7,  4,  4,  7, 10,   4 , 10, 7,  15};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (vtxvarl.size());
	TH1D *vertices[nvar];
	TH1D *ratios[nvar];

	cout << "Current status : vertex range variation ongoing" << endl;
	//boost::progress_display show_progress( nvar );
	for (int i = 0; i < vtxvarl.size(); i++)
	{
		vtxbin[0] = vtxvarl[i];
		vtxbin[1] = vtxvarr[i];
		//cout<<"Makeing a result for "<<Form("%3.0f cm < z <%3.0f cm",vtxvarl[i],vtxvarr[i])<<endl;
		vertices[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, kNoStrVar);
		hset(*vertices[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);

		vertices[i]->SetLineColor(colors[i]);
		vertices[i]->SetLineStyle(lstyle[i]);
		vertices[i]->Draw("histsame][");
		leg->AddEntry(vertices[i], Form("%3.0f cm < #it{z}_{vtx} < %3.0f cm", vtxvarl[i], vtxvarr[i]), "l");
		ratios[i] = (TH1D *)vertices[i]->Clone();
		ratios[i]->Divide(vertices[0]);
		//++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.95);
		ratios[i]->SetMaximum(1.05);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame][");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	TH1D *errorh = (TH1D *)error->Clone();
	TH1D *errorl = (TH1D *)error->Clone();
	double globalmax = 0;
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double maxe = 0, mine = 0;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
 		if (bincenter<-0.5) continue;
		if (bincenter>0.5) continue;
		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			double vale = ratios[j]->GetBinError(i);
			if (val > max)
			{
				max = val;
				maxe = vale;
			}
			if (val < min)
			{
				min = val;
				mine = vale;
			}
		}
		errorh->SetBinContent(i, max);
		errorh->SetBinError(i, maxe);
		errorl->SetBinContent(i, min);
		errorl->SetBinError(i, mine);
		double_t localmax =  0;
		if (fabs(bincenter)<0.8) localmax = TMath::Max(max - 1, 1 - min);
		if (fabs(bincenter) < 0.8 && globalmax < localmax)
			globalmax = localmax;
	}

	TF1 *fh = new TF1("fh", "[0]", -2, 2);
	TF1 *fl = new TF1("fl", "[0]", -2, 2);
	errorh->Fit("fh", "QN");
	errorl->Fit("fl", "QN");
	fh->SetLineWidth(3);
	fl->SetLineWidth(3);
	fh -> Draw("same");
	fl -> Draw("same");

	//Chiara
	double meanChiarah=0;
	double meanChiaral=0;
	int counter=0;
	for (int k = 0; k < errorh->GetNbinsX(); k++)
	{
	  double x, yh, yl;
	  x=errorh->GetBinCenter(k);
	  if (x< -0.5) continue;
	  if (x> 0.5) continue;
	  yh = errorh->GetBinContent(k);
	  yl = errorl->GetBinContent(k);
	  meanChiarah += yh;
	  meanChiaral += yl;
	  counter+=1;
	}
	meanChiarah = meanChiarah/counter;
	meanChiaral = meanChiaral/counter;

	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		//Error->SetPointEYhigh(k, fh->Eval(x) - 1);
		//Error->SetPointEYlow(k, 1 - fl->Eval(x));
		//Error->SetPointEYhigh(k, globalmax);
		//Error->SetPointEYlow(k, globalmax);
		Error->SetPointEYhigh(k, meanChiarah -1); 
		Error->SetPointEYlow(k, 1-meanChiaral); 

	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	vtxbin = defaultvtxbin;
	cout << "vtx bins are set to default ones " << Form("%3.0f cm < %s < %3.0f cm", vtxbin[0], "#it{z}_{vtx}", vtxbin[1]) << endl;

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_vertexsys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}

/*

TGraphAsymmErrors* VertexSystematic(){

	
	TString mcname;
	if (DoDiffTune) mcname = p6difftunename.Data();
	else mcname = p6name.Data();

	Double1D defaultvtxbin = vtxbin;
	//Double1D vtxvarl = {-10,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9};
	//Double1D vtxvarr = { 10, -9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1 );
	canvas->Draw();
	TPad* p = canvas->GetPad(1); //upper pad
	p->SetTickx(); p->SetGridy(0); p->SetGridx(0); p->SetLogy(0); p->cd();
	gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299,0.457,0.925,0.9, NULL,"brNDC");
	leg -> SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (vtxvarl.size());
	TH1D * vertices[nvar];
	TH1D * ratios[nvar];

	cout<<"Current status : vertex range variation ongoing"<<endl; 
	//boost::progress_display show_progress( nvar );
	for (int i=0; i<vtxvarl.size(); i++){
		vtxbin[0] = vtxvarl[i];
		vtxbin[1] = vtxvarr[i];
		//cout<<"Makeing a result for "<<Form("%3.0f cm < z <%3.0f cm",vtxvarl[i],vtxvarr[i])<<endl;
		vertices[i] = BaseCorrection(mcname.Data(),0.2,0.12, false, kNoStrVar);
		hset(*vertices[i], "#it{#eta}", mstring ,0.9,0.9, 0.08,0.08, 0.01,0.01, 0.07,0.07, 510,510);

		vertices[i] -> SetLineColor(colors[i]);
		vertices[i] -> SetLineStyle(lstyle[i]);
		vertices[i] -> Draw("histsame][");
		leg->AddEntry(vertices[i],Form("%3.0f cm < #it{z}_{vtx} < %3.0f cm",vtxvarl[i],vtxvarr[i]),"l");
		ratios[i] = (TH1D*) vertices[i] -> Clone();
		ratios[i] -> Divide(vertices[0]);
	//	++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx(); p->SetGridy(1); p->SetGridx(1); p->SetLogy(0); p->cd();
	gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
	for (int i=0; i<nvar; i++){
	    ratios[i] -> SetMinimum(0.9);
	    ratios[i] -> SetMaximum(1.1);
	    
	    hset(*ratios[i], "#it{#eta}", "Ratio" ,0.9,0.9, 0.08,0.08, 0.01,0.01, 0.07,0.07, 515,505);
	    ratios[i] -> Draw("histsame][");
	}


	TH1D* temp = new TH1D("temp","",40,0.5,1.5);
	TH1D *error = (TH1D*) ratios[0] ->Clone();
	error->SetFillColorAlpha(1, 0.35);
	for (int i=1; i<=ratios[0]->GetNbinsX(); i++){
		temp -> Reset();
		for (int j=0; j<nvar; j++){
			temp->Fill(ratios[j]->GetBinContent(i));
		}
		error -> SetBinContent(i,temp->GetRMS());
	}
	delete temp;

	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);


	TH1D* errorh = (TH1D*) error -> Clone();
	TH1D* errorl = (TH1D*) error -> Clone();
	for (int i=1; i<=error->GetNbinsX(); i++){
		double bincenter = error -> GetXaxis() -> GetBinCenter(i);
		errorh -> SetBinContent(i,error->GetBinContent(i)+1);
		errorl -> SetBinContent(i,1-error->GetBinContent(i));
	}

	TF1* fh = new TF1("fh","[0]",-2,2) ;
	TF1* fl = new TF1("fl","[0]",-2,2) ;
	errorh -> Fit("fh","QN");
	errorl -> Fit("fl","QN");
	fh -> SetLineWidth(3);
	fl -> SetLineWidth(3);
	//fh -> Draw("same");
	//fl -> Draw("same");
	
	for (int k=0; k<Error->GetN();k++){
		double x, y;
		Error->GetPoint(k,x,y);
		Error->SetPoint(k,x,1);
		Error->SetPointEYhigh(k,fh->Eval(x)-1);
		Error->SetPointEYlow(k,1-fl->Eval(x));
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	vtxbin = defaultvtxbin;
	cout<<"vtx bins are set to default ones "<<Form("%3.0f cm < %s < %3.0f cm",vtxbin[0],"#it{z}_{vtx}",vtxbin[1])<<endl;
	
	TString title = "INELg0";
	if (eventclass == kNSD) title = "NSD";
	if (eventclass == kINEL) title = "INEL";
	if (SaveHistogramsAndResults) gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_vertexsys.pdf",dataname.Data(),title.Data(),centtype,cent[0],cent[1]));

	return Error;

}*/

TGraphAsymmErrors *PhiCutSystematic()
{

	TString mcname;
	mcname = p6name.Data();

	Double1D phivar = {-1, 1, 2, 3}; // 0 -  2pi/3, 2pi/3 - 4pi/3, 4pi/3 - 2pi
	vector<TString> tits = {" 0 - 2 #pi", "0 - 2/3 #pi", "2/3 #pi - 4/3 #pi", "4/3 #pi - 2 #pi"};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (phivar.size());
	TH1D *phiranges[nvar];
	TH1D *ratios[nvar];

	cout << "\nAzimuthal acceptance variation ongoing " << endl;
	//boost::progress_display show_progress( nvar );
	for (int i = 0; i < phivar.size(); i++)
	{
		phirange = phivar.at(i);
		phiranges[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, kNoStrVar);
		if (i > 0)
			phiranges[i]->Scale(3.);
		hset(*phiranges[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);

		phiranges[i]->SetLineColor(colors[i]);
		phiranges[i]->SetLineStyle(lstyle[i]);
		phiranges[i]->Draw("histsame][");
		leg->AddEntry(phiranges[i], tits.at(i), "l");
		ratios[i] = (TH1D *)phiranges[i]->Clone();
		ratios[i]->Divide(phiranges[0]);
		//++show_progress;
	}
	leg->Draw();
	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.85);
		ratios[i]->SetMaximum(1.15);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame][");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	TH1D *errorh = (TH1D *)error->Clone();
	TH1D *errorl = (TH1D *)error->Clone();
	double globalmax = 0;
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double maxe = 0, mine = 0;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
		if (bincenter<-0.5) continue;
		if (bincenter>0.5) continue;

		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			double vale = ratios[j]->GetBinError(i);
			if (val > max)
			{
				max = val;
				maxe = vale;
			}
			if (val < min)
			{
				min = val;
				mine = vale;
			}
		}
		double_t localmax =  0;
		if (fabs(bincenter)<0.8) localmax = TMath::Max(max - 1, 1 - min);
		if (fabs(bincenter) < 0.8 && globalmax < localmax)
			globalmax = localmax;
		errorh->SetBinContent(i, max);
		errorh->SetBinError(i, maxe);
		errorl->SetBinContent(i, min);
		errorl->SetBinError(i, mine);
	}

	TF1 *fh = new TF1("fh", "[0]", -0.5, 0.5);
	TF1 *fl = new TF1("fl", "[0]", -0.5, 0.5);
	errorh->Fit("fh", "QN");
	fh->SetLineWidth(3);
	fl->SetLineWidth(3);
	fh->Draw("same");
	fl->Draw("same");

	//Chiara
	double meanChiarah=0;
	double meanChiaral=0;
	int counter=0;
	for (int k = 0; k < errorh->GetNbinsX(); k++)
	{
	  double x, yh, yl;
	  x=errorh->GetBinCenter(k);
	  if (x< -0.5) continue;
	  if (x> 0.5) continue;
	  yh = errorh->GetBinContent(k);
	  yl = errorl->GetBinContent(k);
	  meanChiarah += yh;
	  meanChiaral += yl;
	  counter+=1;
	}
	meanChiarah = meanChiarah/counter;
	meanChiaral = meanChiaral/counter;

	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		Double_t mean = ((errorh->Interpolate(x)-1)+(1-errorl->Interpolate(x)))/2.;
		//Error->SetPointEYhigh(k, fh->Eval(x) - 1);
		//Error->SetPointEYlow(k, fl->Eval(x) - 1);
		Error->SetPointEYhigh(k, meanChiarah -1); 
		Error->SetPointEYlow(k, 1-meanChiaral); 
		//		Error->SetPointEYhigh(k, globalmax);
		//Error->SetPointEYlow(k, globalmax);
		//Error->SetPointEYhigh(k, mean);
	//	Error->SetPointEYlow(k, mean);
	}
	//smoothing
	for (int k = 1; k < Error->GetN()-1; k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);

		Double_t yme = Error->GetErrorYhigh(k - 1);
		Double_t ymp = Error->GetErrorYhigh(k + 1);
		Double_t mean = (yme + ymp) / 2.;
		Error->SetPointEYhigh(k, mean);
		Error->SetPointEYlow(k, mean);
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	phirange = -1;
	cout << "phi cut is set to thedefault one " << endl;
	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_phicut.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}

TGraphAsymmErrors *TrackCutSystematic()
{

	TString mcname;
	mcname = p6name.Data();

	std::vector<Int_t> trackcutvar = {
		kITSTPC2011,
		kHybrid,
		kITSTPC2011dcazdw,
		kITSTPC2011dcazup,
		kITSTPC2011dcardw,
		kITSTPC2011dcarup,
		kITSTPC2011nclutpcdw,
		kITSTPC2011nclutpcup,
		kITSTPC2011chitpcdw,
		kITSTPC2011chitpcup,
		kITSTPC2011globalconsdw,
		kITSTPC2011globalconsup}; //default, up 100% down 50%

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.18, 0.457, 0.81, 0.89, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (trackcutvar.size());
	TH1D *h[nvar];
	TH1D *ratios[nvar];

	cout << "\nCurrent status : Trackcut variation ongoing" << endl;
	//boost::progress_display show_progress( nvar );
	Int_t originaltrackcut = iptvar;
	std::vector<TString> trklegend = {"ITSTPC2011","Hybrid","dcazd","dcazu","dcard","dcaru","ncld","nclu","chitpcd","chitpcup","glocondw","gloconup"};

	for (int i = 0; i < nvar; i++)
	{
		iptvar = trackcutvar[i];
		h[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, kNoStrVar);
		hset(*h[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);

		h[i]->SetLineColor(colors[i]);
		h[i]->SetLineStyle(lstyle[i]);
		h[i]->Draw("histsame");
		leg->AddEntry(h[i], trklegend.at(i).Data(), "l");
		ratios[i] = (TH1D *)h[i]->Clone();
		ratios[i]->Divide(h[0]);
		//++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.95);
		ratios[i]->SetMaximum(1.05);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	TH1D *errorh = (TH1D *)error->Clone();
	TH1D *errorl = (TH1D *)error->Clone();

	double globalmax = 0;
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double maxe = 0, mine = 0;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
		//		if (bincenter<-0.5) continue;
		//		if (bincenter>0.5) continue;

		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			double vale = ratios[j]->GetBinError(i);
			if (val>1.1 || val<0.9)
				continue;
			if (val > max && val < 1.1 && val>1)
			{
				max = val;
				maxe = vale;
			}
			if (val < min && val>0.9 && val<1.)
			{
				min = val;
				mine = vale;
			}
			//max = TMath::Max(max - 1, 1 - min);
			//max = 1 + max;
			//min = 1 - max;
		}
		errorh->SetBinContent(i, max);
		errorh->SetBinError(i, maxe);
		errorl->SetBinContent(i, min);
		errorl->SetBinError(i, mine);
		double_t localmax =  0;
		if (fabs(bincenter) < 0.8)
			localmax = TMath::Max(max - 1, 1 - min);
		if (fabs(bincenter) < 0.8 && globalmax < localmax)
			globalmax = localmax;
	}
	TF1 *fh = new TF1("fh", "[0]", -0.5, 0.5); //fit in eta range I am interested in
	TF1 *fl = new TF1("fl", "[0]", -0.5, 0.5);//fit in eta range I am interested in
	errorh->Fit("fh", "QN");
	errorl->Fit("fl", "QN");
	fh->SetLineWidth(3);
	fl->SetLineWidth(3);
	fh->Draw("same");
	fl->Draw("same");

	double meanChiarah=0;
	double meanChiaral=0;
	int counter=0;
	for (int k = 0; k < errorh->GetNbinsX(); k++)
	{
	  double x, yh, yl;
	  x=errorh->GetBinCenter(k);
	  if (x< -0.5) continue;
	  if (x> 0.5) continue;
	  yh = errorh->GetBinContent(k);
	  yl = errorl->GetBinContent(k);
	  meanChiarah += yh;
	  meanChiaral += yl;
	  counter+=1;
	}
	meanChiarah = meanChiarah/counter;
	meanChiaral = meanChiaral/counter;

	for (int k = 0; k < Error->GetN(); k++)
	{
		double x, y;
		Error->GetPoint(k, x, y);
		Error->SetPoint(k, x, 1);
		Double_t mean = ((errorh->Interpolate(x)-1)+(1-errorl->Interpolate(x)))/2.;
		//		Error->SetPointEYhigh(k, fh->Eval(x) - 1); 
		//		Error->SetPointEYlow(k, -fl->Eval(x) + 1); 
		Error->SetPointEYhigh(k, meanChiarah -1); 
		Error->SetPointEYlow(k, 1-meanChiaral); 
		//Error->SetPointEYhigh(k, (1-fl->Eval(x) ));
		//Error->SetPointEYlow(k, (1-fl->Eval(x) ));
		//Error->SetPointEYhigh(k, globalmax);
		//Error->SetPointEYlow(k, globalmax);
		//Error->SetPointEYhigh(k, fh -> Eval(x)-1  );
		//Error->SetPointEYlow(k, fh -> Eval(x)-1 );
		//Error->SetPointEYhigh(k, mean);
	//	Error->SetPointEYlow(k, mean);
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");
	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_trackcutsys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	cout << "Pt undetectable particle variation bin is set to default ones " << endl;
	iptvar = originaltrackcut;
	return Error;
}

TGraphAsymmErrors *StrangeParticleSystematic()
{

	TString mcname;
	mcname = p6name.Data();

	Int1D var = {kNoStrVar, kStrVarUp, kStrVarDw}; //default, up 100% down 50%
	vector<TString> varc = {"SEF #times 1.0", "SEF #times 1.3", "SEF #times 0.7"};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (var.size());
	TH1D *h[nvar];
	TH1D *ratios[nvar];

	cout << "\nCurrent status : Strange particle variation ongoing" << endl;
	//boost::progress_display show_progress( nvar );
	for (int i = 0; i < nvar; i++)
	{
		h[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, var[i]);
		hset(*h[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);

		h[i]->SetLineColor(colors[i]);
		h[i]->SetLineStyle(lstyle[i]);
		h[i]->Draw("histsame");
		leg->AddEntry(h[i], varc[i].Data(), "l");
		ratios[i] = (TH1D *)h[i]->Clone();
		ratios[i]->Divide(h[0]);
		//++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.994);
		ratios[i]->SetMaximum(1.006);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			if (val > max)
				max = val;
			if (val < min)
				min = val;
		}
		for (int k = 0; k < Error->GetN(); k++)
		{
			double x, y;
			Error->GetPoint(k, x, y);
			Error->SetPoint(k, x, 1);
			if (abs(x - bincenter) < 0.01)
			{
				Error->SetPointEYhigh(k, max - 1);
				Error->SetPointEYlow(k, 1 - min);
			}
		}
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_strangeparsys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}

TGraphAsymmErrors *CentralityEstimatorSystematic()
{
	TString mcname;
	mcname = p6name.Data();

	Int1D var = {kV0M, kV0A, kV0C};
	vector<TString> varc = {"V0M", "V0A", "V0C"};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (var.size());
	TH1D *h[nvar];
	TH1D *ratios[nvar];

	Int_t defaultvalue = centtype;
	for (int i = 0; i < nvar; i++)
	{
		cout << "\n Current status : Centralrity estimator variation ongoing" << endl;
		centtype = var[i];

		h[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, false, kNoStrVar);
		hset(*h[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
		h[i]->SetLineColor(colors[i]);
		h[i]->SetLineStyle(lstyle[i]);
		h[i]->Draw("histsame");
		leg->AddEntry(h[i], varc[i].Data(), "l");
		ratios[i] = (TH1D *)h[i]->Clone();
		ratios[i]->Divide(h[0]);
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.8);
		ratios[i]->SetMaximum(1.2);

		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *temp = new TH1D("temp", "", 40, 0.5, 1.5);
	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		temp->Reset();
		for (int j = 0; j < nvar; j++)
		{
			temp->Fill(ratios[j]->GetBinContent(i));
		}
		error->SetBinError(i, temp->GetRMS());
	}
	delete temp;

	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_centestimatorsys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	cout << "Estimator set to the default value" << endl;
	centtype = defaultvalue;

	return Error;
}

TGraphAsymmErrors *ParticleRatioSystematic()
{

	TString mcname;
	mcname = p6name.Data();

	Double1D kratios = {1., 0.7, 1.3, 1., 1., 1., 1., 1.3, 0.7, 1.3, 0.7};
	Double1D pratios = {1., 1., 1., 0.7, 1.3, 1., 1., 1.3, 0.7, 0.7, 1.3};
	Double1D oratios = {1., 1., 1., 1., 1., 0.7, 1.3, 1., 1., 1., 1.};

	Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.18, 0.45, 0.81, 0.89, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	const int nvar = (kratios.size());
	TH1D *h[nvar];
	TH1D *ratios[nvar];

	cout << "\nCurrent status : particle ratio variation ongoing" << endl;
	//boost::progress_display show_progress( nvar );
	for (int i = 0; i < nvar; i++)
	{
		rkaon = kratios[i];
		rproton = pratios[i];
		ropar = oratios[i];
		h[i] = BaseCorrection(mcname.Data(), 0.2, 0.12, true, kNoStrVar);
		hset(*h[i], "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
		h[i]->SetLineColor(colors[i]);
		h[i]->SetLineStyle(lstyle[i]);
		h[i]->Draw("histsame");
		leg->AddEntry(h[i], Form("kaon #times %2.1f, proton #times %2.1f, other #times %2.1f", rkaon, rproton, ropar), "l");
		ratios[i] = (TH1D *)h[i]->Clone();
		ratios[i]->Divide(h[0]);
		//++show_progress;
	}
	leg->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	for (int i = 0; i < nvar; i++)
	{
		ratios[i]->SetMinimum(0.95);
		ratios[i]->SetMaximum(1.05);
		hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
		ratios[i]->Draw("histsame");
	}

	TH1D *error = (TH1D *)ratios[0]->Clone();
	error->SetFillColorAlpha(1, 0.35);
	TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
	for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
	{
		double max = 1, min = 1;
		double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
		for (int j = 0; j < nvar; j++)
		{
			double val = ratios[j]->GetBinContent(i);
			if (val > max)
				max = val;
			if (val < min)
				min = val;
		}
		for (int k = 0; k < Error->GetN(); k++)
		{
			double x, y;
			Error->GetPoint(k, x, y);
			Error->SetPoint(k, x, 1);
			if (abs(x - bincenter) < 0.01)
			{
				Error->SetPointEYhigh(k, max - 1);
				Error->SetPointEYlow(k, 1 - min);
			}
		}
	}

	Error->SetFillColorAlpha(1, 0.35);
	Error->Draw("2same");

	cout << "Particle ratio variation bins are set to default ones " << endl;
	rkaon = 1.;
	rproton = 1.;
	ropar = 1.;

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%s_cent_%d_%f_%f_parratiosys.pdf", dataname.Data(), title.Data(), centtype, cent[0], cent[1]));

	return Error;
}

//void makeresult(Double1D centMin=0, Double1D centMax=100)
void makeresult(Bool_t DoPhiSyst=1, Bool_t DoVertexSyst = 1, Bool_t DoTrackCutSyst=1, Bool_t DoParticleRatioSyst=1, Bool_t DoStrParticleSyst=1)
{

        //Double1D cent[2] = {centMin, centMax};
	TH1D *corrbyp6 = BaseCorrection(p6name.Data(), 0.2, 0.12, false, kNoStrVar);
	TH1D *corrbyp8 = BaseCorrection(p8name.Data(), 0.2, 0.12, false, kNoStrVar);

	auto avg = (TH1D *)corrbyp6->Clone("avg");
	avg->Add(corrbyp8);
	avg->Scale(0.5);
	avg->SetTitle("");

	TH1D *finalresulthist = avg;

	Filipad2 *can = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
	can->Draw();
	TPad *p = can->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TString title = "INELg0";

	//TLegend  *collisionstitle = new TLegend(0.233254,0.824348,0.766746,0.893913, "pp, #sqrt{#it{s}} = 13 TeV", "brNDC");
	TLegend  *collisionstitle = new TLegend(0.233254,0.824348,0.766746,0.893913, "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "brNDC");
	collisionstitle->SetTextSize(0.0869565);
	collisionstitle->SetBorderSize(0);

	TLegend *ectitle = new TLegend(0.511962,0.713043,0.942584,0.810435, title.Data(), "brNDC");
	ectitle->SetTextSize(0.08);
	ectitle->SetBorderSize(0);


	TLegend *leg = new TLegend(0.217703,0.41913,0.851675,0.70087, "", "brNDC");
	leg->SetTextSize(0.0556522);
	leg->SetBorderSize(0);

	avg->SetMaximum(avg->GetMaximum() * 1.6);
	avg->SetMinimum(avg->GetMinimum() * 0.8);
	hset(*avg, "#it{#eta}", mstring, 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);

	finalresulthist->GetXaxis()->SetRangeUser(-2.5, 2.5);
	finalresulthist->Draw();

	double dndetap6 = corrbyp6->Integral(corrbyp6->GetXaxis()->FindBin(-0.49), corrbyp6->GetXaxis()->FindBin(0.49), "width");
	double dndetap8 = corrbyp8->Integral(corrbyp8->GetXaxis()->FindBin(-0.49), corrbyp8->GetXaxis()->FindBin(0.49), "width");
	cout << "pythia dndeta : " << dndetap6 << " epos dndeta : " << dndetap8 << " difference : " << dndetap6 / dndetap8 - 1 << endl;
	double dndetaavg = avg->Integral(avg->GetXaxis()->FindBin(-0.49), avg->GetXaxis()->FindBin(0.49), "width");
	TH1D *h = (dndetap6 > dndetap8) ? corrbyp6 : corrbyp8;

	TH1D *rh = (TH1D *)h->Clone();
	rh->Divide(avg);

	TH1D *ratios[2];

	if (MakeSymm)
	{
		TH1D *avgtemp = (TH1D *)avg->Clone();
		TH1D *mirrored = (TH1D *)avg->Clone();
		for (auto i = 0; i < avg->GetNbinsX(); i++)
		{
			double ocent = avg->GetXaxis()->GetBinCenter(i);
			double orig = avg->GetBinContent(i);
			double mirror = avgtemp->GetBinContent(avgtemp->GetXaxis()->FindBin(-1 * ocent));
			mirrored->SetBinContent(i, mirror);
			avg->SetBinContent(i, (orig + mirror) / 2.);
		}

		Filipad2 *canvas = new Filipad2(++n, 2, 0.5, 100, 50, 0.7, 1, 1);
		canvas->Draw();
		TPad *p = canvas->GetPad(1); //upper pad
		p->SetTickx();
		p->SetGridy(0);
		p->SetGridx(0);
		p->SetLogy(0);
		p->cd();
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		TLegend *leg = new TLegend(0.248804,0.38087,0.868421,0.902609, NULL, "brNDC");
		leg->SetTextSize(0.062);
		leg->SetBorderSize(0);
		avg->Draw("pzsame");
		mirrored->Draw("pzsame");
		avgtemp->Draw("pzsame");

		ratios[0] = (TH1D *)avgtemp->Clone();
		ratios[1] = (TH1D *)mirrored->Clone();
		ratios[0]->Divide(avg);
		ratios[1]->Divide(avg);

		TH1D *error = (TH1D *)ratios[0]->Clone();
		error->SetFillColorAlpha(1, 0.35);
		TGraphAsymmErrors *Error = new TGraphAsymmErrors(error);
		TH1D *errorh = (TH1D *)error->Clone();
		TH1D *errorl = (TH1D *)error->Clone();
		for (int i = 1; i <= ratios[0]->GetNbinsX(); i++)
		{
			double max = 1, min = 1;
			double maxe = 0, mine = 0;
			double bincenter = ratios[0]->GetXaxis()->GetBinCenter(i);
 			if (bincenter<-0.5) continue;
			if (bincenter>0.5) continue;

			for (int j = 0; j < 2; j++)
			{
				double val = ratios[j]->GetBinContent(i);
				double vale = ratios[j]->GetBinError(i);
				if (val > max)
				{
					max = val;
					maxe = vale;
				}
				if (val < min)
				{
					min = val;
					mine = vale;
				}
			}
			errorh->SetBinContent(i, max);
			errorh->SetBinError(i, maxe);
			errorl->SetBinContent(i, min);
			errorl->SetBinError(i, mine);
		}

		TF1 *fh = new TF1("fh", "[0]+[1]*x*x", -2, 2);
		TF1 *fl = new TF1("fl", "[0]+[1]*x*x", -2, 2);
		errorh->Fit("fh", "QN");
		errorl->Fit("fl", "QN");
		fh->SetLineWidth(3);
		fl->SetLineWidth(3);

		for (int k = 0; k < Error->GetN(); k++)
		{
			double x, y;
			Error->GetPoint(k, x, y);
			Error->SetPoint(k, x, 1);
			Error->SetPointEYhigh(k, fh->Eval(x) - 1);
			Error->SetPointEYlow(k, 1 - fl->Eval(x));
		}

		p = canvas->GetPad(2); //lower pad
		p->SetTickx();
		p->SetGridy(1);
		p->SetGridx(1);
		p->SetLogy(0);
		p->cd();
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		for (int i = 0; i < 2; i++)
		{
			ratios[i]->SetMinimum(0.85);
			ratios[i]->SetMaximum(1.15);

			hset(*ratios[i], "#it{#eta}", "Ratio", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
			ratios[i]->Draw("histsame][");
		}

		Error->SetFillColorAlpha(1, 0.35);
		Error->Draw("2same");

		symsystematic = Error;
	}

	double dndetaavgdt = avg->Integral(avg->GetXaxis()->FindBin(-0.49), avg->GetXaxis()->FindBin(0.49), "width");
	double dndetaavgdt10 = avg->Integral(avg->GetXaxis()->FindBin(-0.79), avg->GetXaxis()->FindBin(0.79), "width")/1.6;
	//double dndetaavgdt = avg->Integral(avg->GetXaxis()->FindBin(-0.99),avg->GetXaxis()->FindBin(0.99),"width");

	TGraphAsymmErrors *AVG = new TGraphAsymmErrors(avg);
	TGraphAsymmErrors *sysh = new TGraphAsymmErrors(rh);

	for (int k = 0; k < AVG->GetN(); k++)
	{
		double ex, ey;
		sysh->GetPoint(k, ex, ey);
		double x, y;
		AVG->GetPoint(k, x, y);
		if (x == -1.9)
		{
			AVG->SetPoint(k - 1, -2.1, y);
			AVG->SetPointEYhigh(k - 1, (ey - 1) * y);
			AVG->SetPointEYlow(k - 1, (ey - 1) * y);
		}
		AVG->SetPointEYhigh(k, (ey - 1) * y);
		AVG->SetPointEYlow(k, (ey - 1) * y);
	}
	AVG->SetFillColorAlpha(1, 0.35);

	Int_t keq19 = 0;

	TGraphAsymmErrors *finalresult = AVG;
	TGraphAsymmErrors *finalresultcorr = (TGraphAsymmErrors *)finalresult->Clone();
	TGraphAsymmErrors *finalresultuncorr = (TGraphAsymmErrors *)finalresult->Clone();

	Double_t errorh = 0;
	Double_t errorhTrackCut = 0;
	Double_t errorhVertex = 0;
	Double_t errorhPhi = 0;
	Double_t errorhParticleRatio = 0;
	Double_t errorhStrParticle = 0;
	Double_t errorl = 0;
	Double_t errorlTrackCut = 0;
	Double_t errorlVertex = 0;
	Double_t errorlPhi = 0;
	Double_t errorlParticleRatio = 0;
	Double_t errorlStrParticle = 0;
	Double_t errorh10 = 0;
	Double_t errorl10 = 0;
	auto *erroronly = (TGraphAsymmErrors *)finalresult->Clone();
	if (DoSystematic)
	{
        if (DoTrackCutSyst) trackcutsystematic = TrackCutSystematic();
		if (DoParticleRatioSyst) particleratiosystematic = ParticleRatioSystematic();
		if (DoStrParticleSyst) stranegeparticlesystematic = StrangeParticleSystematic();
		if (DoVertexSyst) vertexsystematic = VertexSystematic();
		if (DoPhiSyst) phicutsystematic = PhiCutSystematic();
		if (DoDiffTune){
			diffractionsystematic = DiffractionSystematic();
			diffractionratiosystematic = DiffractionRatioSystematic();
		}

		//triggerbiassystematic = TriggerBiasSystematic();
		//centralityestimatorsystematic = CentralityEstimatorSystematic();

		p = can->GetPad(1); //lower pad
		for (int k = 0; k < finalresult->GetN(); k++)
		{
			double eyhuncorr, eyluncorr;
			double x, y;
			finalresult->GetPoint(k, x, y);
			eyhuncorr = finalresult->GetErrorYhigh(k);
			eyluncorr = finalresult->GetErrorYlow(k);
			eyhuncorr /= y;
			eyluncorr /= y;
			double ediffh = 0;
			double ediffl = 0;
			double ediffh2 = 0;
			double ediffl2 = 0;
			double ematching = 0;

			Double_t egendependence = 0;
			Double_t material = 0.002;

			if (dataname.Contains("16l")){
			  if ( triggtype == kMBAND015 )
			    egendependence = 0.009;
			  /*
			  else if (triggtype == kMBAND050)
			    egendependence = 0.013;
			  */
			  else if (triggtype == kMBAND100)
			    egendependence = 0.019;
			  else if (triggtype == kMBAND200)
			    egendependence = 0.017;
			}
			if (dataname.Contains("17p"))
			  {
			    if ( triggtype == kMBAND015 )
			      egendependence = 0.01;
			    /*
			    else if (triggtype == kMBAND050)
			      egendependence = 0.012;
			    */
			    else if (triggtype == kMBAND100)
			      egendependence = 0.02;
			    else if (triggtype == kMBAND200)
			      egendependence = 0.02;
			  }

			if (DoDiffTune)
			{
				ediffh = diffractionratiosystematic->GetErrorYhigh(k);
				ediffl = diffractionratiosystematic->GetErrorYlow(k);
				ediffh2 = diffractionsystematic->GetErrorYhigh(k);
				ediffl2 = diffractionsystematic->GetErrorYlow(k);
			}
			double esymh = 0, esyml = 0;
			if (MakeSymm)
			{
				esymh = symsystematic->GetErrorYhigh(k);
				esyml = symsystematic->GetErrorYlow(k);
			}

			double eptunseenh = 0;
			double eptunseenl = 0;
			if (DoTrackCutSyst){
			  eptunseenh = trackcutsystematic->GetErrorYhigh(k);
			  eptunseenl = trackcutsystematic->GetErrorYlow(k);
			}
			double eparratioh =0;
			double eparratiol =0;
			if (DoParticleRatioSyst){
			  eparratioh = particleratiosystematic->GetErrorYhigh(k);
			  eparratiol = particleratiosystematic->GetErrorYlow(k);
			}
			double estrparh =0;
			double estrparl =0;
			if (DoStrParticleSyst){
			  estrparh = stranegeparticlesystematic->GetErrorYhigh(k);
			  estrparl = stranegeparticlesystematic->GetErrorYlow(k);
			}
			double evertexh = 0;
			double evertexl = 0;
			if (DoVertexSyst){
			  evertexh = vertexsystematic->GetErrorYhigh(k);
			  evertexl = vertexsystematic->GetErrorYlow(k);
			}
			//double eestimatorh = centralityestimatorsystematic->GetErrorYhigh(k);
			//double eestimatorl = centralityestimatorsystematic->GetErrorYlow(k);
			double ephicuth=0;
			double ephicutl=0;
			if (DoPhiSyst){
			  ephicuth = phicutsystematic->GetErrorYhigh(k);
			  ephicutl = phicutsystematic->GetErrorYlow(k);
			}
			//double etriggerh = triggerbiassystematic->GetErrorYhigh(k);
			//double etriggerl = triggerbiassystematic->GetErrorYlow(k);
			eyhuncorr = sqrt(egendependence*egendependence  + evertexh * evertexh + ephicuth * ephicuth);
			eyluncorr = sqrt(egendependence*egendependence  + evertexl * evertexl + ephicutl * ephicutl);
			double eyhcorr = sqrt( material*material + eptunseenh * eptunseenh+ eparratioh * eparratioh + estrparh * estrparh);
			double eylcorr = sqrt( material*material  + eptunseenl * eptunseenl+ eparratiol * eparratiol + estrparl * estrparl);
			double eyhall = sqrt(eyhcorr * eyhcorr + eyhuncorr * eyhuncorr);
			double eylall = sqrt(eylcorr * eylcorr + eyluncorr * eyluncorr);
			if (x>-0.5 && x < 0.5){
			  /*
			  cout << "\nCoordinates: " << x << " - " << y << endl;
			  cout << std::setprecision(4);
			  cout << "\n\e[35meyhall " << eyhall << endl;
			  cout << "It's composed of " << endl;
			  cout << "Material " << material << endl;
			  cout << "Track cut systematic " << eptunseenh << endl;
			  cout << "Particle Ratio systematic " << eparratioh << endl;
			  cout << "Strange particle systematic " << estrparh << endl;
			  cout << "Phi systematic " << ephicuth << endl;
			  cout << "Vertex " << evertexh << endl;
			  cout << "Energy dependence " << egendependence << "\e[39m " << endl;
			  */
			}
			erroronly->SetPoint(k, x, 1);
			erroronly->SetPointEYhigh(k, eyhall);
			erroronly->SetPointEYlow(k, eylall);
			finalresult->SetPointEYhigh(k, eyhall * y);
			finalresult->SetPointEYlow(k, eylall * y);
			finalresultcorr->SetPointEYhigh(k, eyhcorr * y);
			finalresultcorr->SetPointEYlow(k, eylcorr * y);
			finalresultuncorr->SetPointEYhigh(k, eyhuncorr * y);
			finalresultuncorr->SetPointEYlow(k, eyluncorr * y);
			if (x > -0.5 && x < 0.5)
			{
				errorh += finalresult->GetErrorYhigh(k) * 0.1;
				errorl += finalresult->GetErrorYlow(k) * 0.1;
				if (DoTrackCutSyst){
				  errorhTrackCut += trackcutsystematic->GetErrorYhigh(k) * 0.1;
				  errorlTrackCut += trackcutsystematic->GetErrorYlow(k) * 0.1;
				}
				if (DoVertexSyst){
				  errorhVertex += vertexsystematic->GetErrorYhigh(k) * 0.1;
				  errorlVertex += vertexsystematic->GetErrorYlow(k) * 0.1;
				}
				if (DoPhiSyst){
				  errorhPhi += phicutsystematic->GetErrorYhigh(k) * 0.1;
				  errorlPhi += phicutsystematic->GetErrorYlow(k) * 0.1;
				}
				if (DoStrParticleSyst){
				  errorhStrParticle = stranegeparticlesystematic->GetErrorYhigh(k) *0.1;
				  errorlStrParticle = stranegeparticlesystematic->GetErrorYlow(k)*0.1;
				}
				if (DoParticleRatioSyst){
				  errorhParticleRatio = particleratiosystematic->GetErrorYhigh(k) *0.1;
				  errorlParticleRatio = particleratiosystematic->GetErrorYlow(k)*0.1;
				}
			}
			if (x > -0.8 && x < 0.8)
			{
				errorh10 += finalresult->GetErrorYhigh(k) * 0.1;
				errorl10 += finalresult->GetErrorYlow(k) * 0.1;
			}
		}
	}

	finalresult->SetFillColorAlpha(1, 0.35);
	finalresult->Draw("2psame");

	//errorh /=1.2;
	//errorl /=1.2;
	auto p6 = MakeModelResults(p6name, false);
	auto p8 = MakeModelResults(p8name, false);

	p = can->GetPad(1); //lower pad
	p6->SetMarkerColor(2);
	p6->SetLineColor(2);
	p8->SetMarkerColor(2);
	p8->SetLineColor(2);

	p6->Draw("histsame");
	if (dataname.Contains("16l")) p8->Draw("histsame");
	double dndetaat0p6 = p6->Integral(p6->GetXaxis()->FindBin(-0.49), p6->GetXaxis()->FindBin(0.49), "width");
	double dndetaat0p8 = p8->Integral(p8->GetXaxis()->FindBin(-0.49), p8->GetXaxis()->FindBin(0.49), "width");
	double dndetaat0p62 = p6->Integral(p6->GetXaxis()->FindBin(-0.79), p6->GetXaxis()->FindBin(0.79), "width")/1.6;
	double dndetaat0p82 = p8->Integral(p8->GetXaxis()->FindBin(-0.79), p8->GetXaxis()->FindBin(0.79), "width")/1.6;
	//leg->AddEntry(finalresult, Form("ALICE, |#eta|<0.5, %4.2f^{+%2.2f}_{-%2.2f}", dndetaavgdt, errorh, errorl), "lp3");
	//leg->AddEntry(finalresult, Form("ALICE, %4.2f#pm%2.2f (|#it{#eta}|<0.8)", dndetaavgdt10, errorh10/1.6, errorl10/1.6), "lp3");
	leg->AddEntry(finalresult, Form("ALICE"), "lp3");
	leg->AddEntry(p6, Form("EPOS-LHC"), "lp3");
	cout << "dNdeta in 0.5 : " << dndetaavgdt << " +" << errorh << " -" << errorl << endl;
	cout << "dNdeta in 1.0 : " << dndetaavgdt10  << " " << errorh10 / 2 << " " << errorl10 / 2 << endl;
	cout << "Relative errors: " << endl;
	cout << "dNdeta in 0.5 (track cut syst) : " << dndetaavgdt << " " << errorhTrackCut << " " << errorlTrackCut << endl;
	cout << "dNdeta in 0.5 (vertex syst) : " << dndetaavgdt << " " << errorhVertex << " " << errorlVertex << endl;
	cout << "dNdeta in 0.5 (phi syst) : " << dndetaavgdt << " " << errorhPhi << " " << errorlPhi << endl;
	cout << "dNdeta in 0.5 (particle ratio syst) : " << dndetaavgdt << " " << errorhParticleRatio << " " << errorlParticleRatio << endl;
	cout << "dNdeta in 0.5 (strange particle syst) : " << dndetaavgdt << " " << errorhStrParticle << " " << errorlStrParticle << endl;

	//p6->Scale(dndetaavgdt10 / dndetaat0p62 / 1.18);
	//p8->Scale(dndetaavgdt10 / dndetaat0p62 / 1.18);
	//p6->Scale(dndetaavgdt10 / dndetaat0p62 );
	//p8->Scale(dndetaavgdt10 / dndetaat0p62 );
	double epos25 = p6->Integral(p6->GetXaxis()->FindBin(-2.499), p6->GetXaxis()->FindBin(2.499), "width");
	cout << "EPOS result: " << epos25 << endl;

	if (dataname.Contains("16l"))
	{
		//leg->AddEntry(p6, Form("PYTHIA8 Monash 2013 |#eta|<0.5, %4.2f", dndetaat0p6), "l");
		//leg->AddEntry(p8, Form("EPOS LHC |#eta|<0.5, %4.2f", dndetaat0p8), "l");
		//leg->AddEntry(p6, Form("PYTHIA8 Monash 2013, %4.2f (|#it{#eta}|<0.8)", dndetaat0p62), "l");
		//leg->AddEntry(p8, Form("EPOS LHC, %4.2f (|#it{#eta}|<0.8)", dndetaat0p82), "l");
		//leg->AddEntry(p6, Form("PYTHIA8 Monash 2013"), "l");
		//leg->AddEntry(p8, Form("EPOS LHC"), "l");
	}
	else if (dataname.Contains("17p"))
	{
		//leg->AddEntry(p6, Form("PYTHIA8-Monash 2013 |#eta|<0.5, %4.2f", dndetaat0p8), "l");
		//leg->AddEntry(p6, Form("PYTHIA8-Monash 2013 |#eta|<0.8, %4.2f", dndetaat0p82), "l");
		leg->AddEntry(p6, Form("PYTHIA8 Monash 2013, %4.2f (|#it{#eta}|<0.8)", dndetaat0p62), "l");
	}
	

	cout << "Final result summy table ==========================" << endl;
	cout << "Muldiplicity bin : " << Form("%4.1f - %4.1f  %s", cent[0], cent[1], "%") << endl;
	cout << "Vertex bin       : " << Form("%2.0f - %2.0f cm", vtxbin[0], vtxbin[1]) << endl;
	cout << "\n\n"
		 << endl;
	cout << "    eta          data         PYTHIA6         PYTHIA8 " << endl;
	for (int k = 0; k < finalresult->GetN(); k++) //this is a loop over eta values
	{
		double x, y;
		finalresult->GetPoint(k, x, y);
		double eyh = finalresult->GetErrorYhigh(k);
		double eyl = finalresult->GetErrorYlow(k);
		double p6value = p6->GetBinContent(p6->GetXaxis()->FindBin(x));
		double p8value = p8->GetBinContent(p8->GetXaxis()->FindBin(x));
		if(y == 0) continue;
		cout << Form("    %2.1f   ", x) << Form("   %4.2f + %.2f (%4.2f) - %.2f (%4.2f)  ", y, eyh, eyh / y * 100, eyl, eyl / y * 100) << Form("   %4.2f   %4.2f   ", p6value, p8value) << endl;
	}

	TGraphAsymmErrors *atlasdata = nullptr;
	leg->Draw();
	TLegend *ectitle2 = new TLegend(0.511962, 0.613043, 0.942584, 0.71, Form ("Cent: %d-%d", int (cent[0]), int(cent[1])), "brNDC");
	ectitle2->SetTextSize(0.08);
	ectitle2->SetBorderSize(0);
	ectitle2->Draw();
	collisionstitle->Draw();
	ectitle->Draw();

	p = can->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	auto rp6 = (TH1D *)p6->Clone();
	rp6->Divide(finalresulthist);
	auto rp8 = (TH1D *)p8->Clone();
	rp8->Divide(finalresulthist);

	hset(*rp6, "#it{#eta}", "Ratio (/ALICE)", 0.9, 0.9, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 515, 505);
	rp6->SetMinimum(0.92);
	rp6->SetMaximum(1.06);
	auto eventgeneratordependence = (TH1D *)finalresulthist->Clone();
	eventgeneratordependence->Divide(corrbyp6);

	//rp6->GetXaxis()->SetRangeUser(etabin[0], etabin[1]);
	rp6->GetXaxis()->SetRangeUser(-2.5, 2.5);
	rp6->SetMaximum(1.3);
	rp6->SetMinimum(0.3);
	rp6->Draw("histsame");
	if (dataname.Contains("16l")) rp8->Draw("histsame");
	//eventgeneratordependence->Draw("histsame");
	erroronly->Draw("2same");

	TGraphAsymmErrors *atlasoveralice = nullptr;
	if (SaveResultsOnly)
	{
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%d.pdf", dataname.Data(),eventclass));
		TFile *f = new TFile(Form("../FinalResults/%s_%d_cent_%d_%d.root", dataname.Data(), eventclass,int(cent[0]),int(cent[1])), "recreate");
		finalresult->SetName("dataall");
		finalresultcorr->SetName("datacorr");
		finalresultuncorr->SetName("datuncorr");
		if (dataname.Contains("16l"))
		{
			p6->SetName("p8");
			p8->SetName("epos");
		}
		else
		{
			p6->SetName("p8");
			p8->SetName("p8");
		}
		finalresult->Write();
		finalresultcorr->Write();
		finalresultuncorr->Write();
		p6->Write();
		p8->Write();
	}
}

void makeresult(double centl, double centr, int ctype, int trigtype, TString dname, TString mname)
{
	cent[0] = centl;
	cent[1] = centr;
	centtype = ctype;
	if (centtype == 1)
		mccent = {(cent[0] == 0 && cent[1] == 0.01) ? 0.001 : cent[0], cent[1] <= 1 ? 1 : cent[1]};
	else
		mccent = {(cent[0] == 0 && cent[1] == 1) ? 0.1 : cent[0], cent[1] <= 1 ? 1 : cent[1]};
	triggtype = trigtype;
	dataname = dname;
	p6name = mname;
	p8name = mname;
	makeresult();
}

void allcentmakeresult()
{

	Double1D allcent = {0, 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	//Double1D  allcent={40,50,70,100};
	//Double1D  allcent={0,1,5,10,15,20,30,40,50,70,100};
	//Double1D allcent={0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
	for (auto i = 0; i < allcent.size() - 1; i++)
	{
		cent.at(0) = allcent.at(i);
		cent.at(1) = allcent.at(i + 1);
		mccent = cent;
		makeresult();
	}
}

TH1D *MakeModelResults(TString mcname, Bool_t makesymmetric)
{
	auto clist = LoaddndetaResultList(mcname.Data(), "output");
	Int_t centselection = (cent[0] == 0 && cent[1] == 100) ? -1 : 1;
	auto Hkinezvtx = BSTHnSparseHelper::Load("hkinezvtx", clist);
	Hkinezvtx.SetBin("Cent", mccent);
	auto hz = Hkinezvtx.GetTH1(Form("hmodel%s-", dataname.Data()), 2, {eventclass, centselection, -1, centtype});
	auto normfac = hz->Integral(hz->GetXaxis()->FindBin(vtxbin[0]), hz->GetXaxis()->FindBin(vtxbin[1]));

	auto Hkinedndeta = BSTHnSparseHelper::Load("hkinedndeta", clist);
	Hkinedndeta.SetBin("ParticleType", {kPion, kOPar});
	Hkinedndeta.SetBin("Z", vtxbin);
	Hkinedndeta.SetBin("Cent", mccent);
	auto fr = Hkinedndeta.GetTH1(Form("htp%s-", dataname.Data()), 5, {eventclass, centselection, 1, 1, 1, -1, centtype, -1});

	fr->SetName(Form("Model%s_dndeta", mcname.Data()));
	fr->Scale(1. / normfac, "width");
	fr->SetLineColor(modelcolors[modelcounter]);
	fr->SetMarkerColor(modelcolors[modelcounter]);
	fr->SetLineWidth(3);

	modelcounter++;

	if (makesymmetric)
	{
		auto hdndetaReflection = (TH1D *)fr->Clone();
		hdndetaReflection->Reset();
		for (auto i = 1; i <= hdndetaReflection->GetNbinsX(); i++)
		{
			double bincent = hdndetaReflection->GetXaxis()->GetBinCenter(i);
			int bin = fr->GetXaxis()->FindBin(-1 * bincent);
			hdndetaReflection->SetBinContent(i, fr->GetBinContent(bin));
			hdndetaReflection->SetBinError(i, fr->GetBinError(bin));
		}
		fr->Add(hdndetaReflection);
		fr->Scale(0.5);
	}

	return fr;
}

/*

TH1D* MakeModelResults(TString mcname, Bool_t makesymmetric){
	auto clist = LoaddndetaResultList( mcname.Data(), "output");	
	Int_t centselection = (cent[0]==0 && cent[1] == 100) ? -1: 1;
	auto HcentGeo = BSTHnSparseHelper::Load( "hcentGeo", clist );
	auto hcentGeo = HcentGeo.GetTH1(Form("hmodel%s-",dataname.Data()),1,{centtype, -1 });
	auto normfac = hcentGeo->Integral(hcentGeo->GetXaxis()->FindBin(cent[0]+0.0001),hcentGeo->GetXaxis()->FindBin(cent[1]-0.0001));
	//new TCanvas;
	hcentGeo -> Scale(1.,"width");
	hcentGeo -> SetMinimum(0);
	//hcentGeo->Draw();
	
	auto Hkinedndeta = BSTHnSparseHelper::Load( "hkinedNdEtaGeoCent", clist );
  //CreateTHnSparse( "hkinedNdEtaGeoCent","",3,{binCent,binEta,binCentClass},"s");
	//Hkinedndeta.SetBin("ParticleType",{kPion,kOPar});
	Hkinedndeta.SetBin("Cent",cent);
	auto fr = Hkinedndeta.GetTH1(Form("htp%s-",dataname.Data()),1,{1, -1,1});	
	fr -> Scale(1./normfac,"width");
	

	fr->SetName(Form("Model%s_dndeta",mcname.Data()));
	fr->SetLineColor(modelcolors[modelcounter]);
	fr->SetMarkerColor(modelcolors[modelcounter]);
	fr->SetLineWidth(3);
	

	if (makesymmetric) {
		auto hdndetaReflection = (TH1D*) fr -> Clone();
		hdndetaReflection -> Reset();
		for (auto i = 1; i <= hdndetaReflection -> GetNbinsX() ; i++){
			double bincent  = hdndetaReflection -> GetXaxis() -> GetBinCenter(i);
			int bin =  fr -> GetXaxis()-> FindBin(-1*bincent);
			hdndetaReflection -> SetBinContent (i,fr->GetBinContent(bin));
			hdndetaReflection -> SetBinError (i,fr->GetBinError(bin));
		}
		fr -> Add(hdndetaReflection);
		fr -> Scale(0.5);
	}
	
	modelcounter++;

	return fr;
}
*/
void TriggerEfficiency()
{

	auto clistmc = LoaddndetaResultList(p6name.Data(), "output");
	auto clistdata = LoaddndetaResultList(dataname.Data(), "output");
	//Int_t centselection = (cent[0]==0 && cent[1] == 100) ? -1: 1;
	Int_t centselection = (cent[0] == 0 && cent[1] == 100) ? 1 : 1;

	// Load z vtx distributions
	auto Hkinezvtx = BSTHnSparseHelper::Load("hkinezvtx", clistmc);
	auto Hmczvtx = BSTHnSparseHelper::Load("hreczvtx", clistmc);
	auto Hdatazvtx = BSTHnSparseHelper::Load("hreczvtx", clistdata);

	//{ 0EventClass, 1Cent 2z_vtx }
	Double1D centbin = {0, 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	centbin = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	//centbin = {0,5,10,15,20,30,40,50,70,100};
	Hkinezvtx.SetBin("Cent", centbin);
	Hmczvtx.SetBin("Cent", centbin);

	double newtotal = 0;

	for (auto i : Hkinezvtx.BinRange("Cent"))
	{
		auto hkinez = Hkinezvtx.GetTH1(Form("h%s-", dataname.Data()), 2, {eventclass, i, -1, centtype});
		auto hmcz = Hmczvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {eventclass, triggtype, i, -1, centtype});

		double bw = (centbin[i] - centbin[i - 1]);
		cout << "binwidth : " << bw << endl;
		double te = hmcz->Integral(hmcz->GetXaxis()->FindBin(-9.99), hmcz->GetXaxis()->FindBin(9.99)) / hkinez->Integral(hkinez->GetXaxis()->FindBin(-9.99), hkinez->GetXaxis()->FindBin(9.99));
		if (i == 1)
			te = 1;
		cout << "Trigger Efficiency = " << te << endl;
		newtotal += bw / te;
	}

	cout << "newtotal " << newtotal << endl;
	double accumbinwidth = 0;

	for (auto i : Hkinezvtx.BinRange("Cent"))
	{
		auto hkinez = Hkinezvtx.GetTH1(Form("h%s-", dataname.Data()), 2, {eventclass, i, -1, centtype});
		auto hmcz = Hmczvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {eventclass, triggtype, i, -1, centtype});

		double bw = (centbin[i] - centbin[i - 1]);
		cout << "bin : " << centbin[i - 1] << "--" << centbin[i] << endl;
		double te = hmcz->Integral(hmcz->GetXaxis()->FindBin(-9.99), hmcz->GetXaxis()->FindBin(9.99)) / hkinez->Integral(hkinez->GetXaxis()->FindBin(-9.99), hkinez->GetXaxis()->FindBin(9.99));
		if (i == 1)
			te = 1;
		double newbinwidth = bw / te / newtotal * 100;
		cout << accumbinwidth << endl;
		accumbinwidth += newbinwidth;
		cout << accumbinwidth << endl;
		//cout<<"newbinwidth "<<newbinwidth*100<<endl;
	}

	Hkinezvtx.SetBin("Cent", {0, 100});
	Hmczvtx.SetBin("Cent", {0, 100});
	auto hkinez = Hkinezvtx.GetTH1(Form("h%s-", dataname.Data()), 2, {eventclass, -1, -1, centtype});
	auto hmcz = Hmczvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {eventclass, triggtype, -1, -1, centtype});
	cout << " 0-100 % trigger Efficiency : " << hmcz->Integral(hmcz->GetXaxis()->FindBin(-9.99), hmcz->GetXaxis()->FindBin(9.99)) / hkinez->Integral(hkinez->GetXaxis()->FindBin(-9.99), hkinez->GetXaxis()->FindBin(9.99)) << endl;
}

TH1D *BaseCorrection(TString mcname, double SDR, double DDR, bool dopariclevariation, Int_t strangevar)
{

	auto clistmc = LoaddndetaResultList(mcname.Data(), "output");
	auto clistdata = LoaddndetaResultList(dataname.Data(), "output");
	//Int_t centselection = (cent[0]==0 && cent[1] == 100) ? -1: 1;
	Int_t centselection = (cent[0] == 0 && cent[1] == 100) ? -1 : 1;

	// Load z vtx distributions
	auto Hkinezvtx = BSTHnSparseHelper::Load("hkinezvtx", clistmc);
	auto Hmczvtx = BSTHnSparseHelper::Load("hreczvtx", clistmc);
	auto Hdatazvtx = BSTHnSparseHelper::Load("hreczvtx", clistdata);

	//{ 0EventClass, 1Cent 2z_vtx }
	Hkinezvtx.SetBin("Cent", mccent);
	//--------------------
	auto hkinez = Hkinezvtx.GetTH1(Form("h%s-", dataname.Data()), 2, {eventclass, centselection, -1, centtype});
	//0EventClass, 1TriggClass  2Cent 3z_vtx }
	Hmczvtx.SetBin("Cent", cent);
	auto hmcz = Hmczvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {eventclass, triggtype, centselection, -1,centtype});

	//---------------------

	//auto hkinez = Hkinezvtx.GetTH1(Form("h%s-", dataname.Data()), 2, {eventclass, centselection, -1, centtype});
	//Hmczvtx.SetBin("Cent", mccent);
	//auto hmcz = Hmczvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {eventclass, triggtype, centselection, -1, centtype});

	cout << "Trigger Efficiency : " << hmcz->Integral(hmcz->GetXaxis()->FindBin(-9.99), hmcz->GetXaxis()->FindBin(9.99)) / hkinez->Integral(hkinez->GetXaxis()->FindBin(-9.99), hkinez->GetXaxis()->FindBin(9.99)) << endl;

	Hdatazvtx.SetBin("Cent", cent);
	auto hdataz = Hdatazvtx.GetTH1(Form("h%s-", dataname.Data()), 3, {kDATA, triggtype, centselection, -1, centtype});

	cout << "# of event : " << hdataz->Integral(hdataz->GetXaxis()->FindBin(-9.99), hdataz->GetXaxis()->FindBin(9.99)) << endl;

	auto hcorrfevt = (TH1D *)hkinez->Clone();
	hcorrfevt->Divide(hmcz);
	hdataz->Multiply(hcorrfevt);

	//Calculation of strange particle enhancement/suppression in MC
	auto Hmcstrange = BSTHnSparseHelper::Load("hv0eta", clistmc);
	Hmcstrange.SetBin("Cent", mccent);
	// 0TriggClass, 1binCent 2binV0Type 3binEta
	auto hmcstrange = Hmcstrange.GetTH1("hstrmc", 3, {triggtype, centselection, -1, -1, centtype});
	auto Hecmc = BSTHnSparseHelper::Load("heventcount", clistmc);
	auto hecmc = Hecmc.GetTH1("hecmc", 1, {triggtype, -1, centtype});
	double necmc = 0;
	if (centselection == -1)
		necmc = hecmc->Integral(0, -1);
	else
		necmc = hecmc->Integral(hecmc->GetXaxis()->FindBin(cent[0]), hecmc->GetXaxis()->FindBin(cent[1]));
	hmcstrange->Scale(1. / necmc, "width");

	auto Hdatastrange = BSTHnSparseHelper::Load("hv0eta", clistdata);
	// 0TriggClass, 1binCent 2binV0Type 3binEta
	Hdatastrange.SetBin("Cent", cent);

	auto hdatastrange = Hdatastrange.GetTH1("hstrdata", 3, {triggtype, centselection, -1, -1, centtype});
	auto Hecdata = BSTHnSparseHelper::Load("heventcount", clistdata);
	auto hecdata = Hecdata.GetTH1("hecdata", 1, {triggtype, -1, centtype});
	double necdata = 0;
	if (centselection == -1)
		necdata = hecdata->Integral(0, -1);
	else
		necdata = hecdata->Integral(hecdata->GetXaxis()->FindBin(cent[0]), hecdata->GetXaxis()->FindBin(cent[1]));
	hdatastrange->Scale(1. / necdata, "width");

	hdatastrange->Divide(hmcstrange);
	TF1 *constant = new TF1("constant", "[0]", -0.8, 0.8);
	hdatastrange->Fit("constant", "QN", "same", -0.8, 0.8);

	if (0)
	{
		new TCanvas;
		setpad(gPad);
		gPad->SetTitle("SCF");
		hdatastrange->SetMinimum(0);
		hdatastrange->SetMaximum(3);
		hset(*hdatastrange, "#it{#eta}", "SSF", 1.0, 0.9, 0.07, 0.07, 0.01, 0.001, 0.06, 0.06, 510, 510);
		hdatastrange->SetTitle("");
		hdatastrange->Draw();
		constant->Draw("same");
		auto tex = new TLatex(-0.4, 2.4, "MC and ALICE DATA pp, #sqrt{s} = 5.02 TeV");
		tex->SetTextSize(0.04);
		tex->Draw();
		tex = new TLatex(0.1, 0.9, Form("Fitted SSF = %4.3f #pm %4.3f", constant->GetParameter(0), constant->GetParError(0)));
		tex->SetTextSize(0.04);
		tex->Draw();
	}
	//end of Calculation of strange particle enhancement/suppression in MC

	// Load Kine (eta,z)
	auto Hkinedndeta = BSTHnSparseHelper::Load("hkinedndeta", clistmc);
	Hkinedndeta.SetBin("Cent", mccent);
	// order : 0evt_class 1centrality 2zvtx_mcgen 3PID 4ptvar 5eta_mcgen
	TH2D *hkinedndetapion = nullptr;
	hkinedndetapion = Hkinedndeta.GetTH2(Form("htp1%s-", dataname.Data()), 5, 2, {eventclass, centselection, 1, kPion, 1, -1, centtype, phirange});
	auto hkinedndetakaon = Hkinedndeta.GetTH2(Form("htp2%s-", dataname.Data()), 5, 2, {eventclass, centselection, 1, kKaon, 1, -1, centtype, phirange});
	auto hkinedndetaproton = Hkinedndeta.GetTH2(Form("htp3%s-", dataname.Data()), 5, 2, {eventclass, centselection, 1, kProtonBK, 1, -1, centtype, phirange});
	auto hkinedndetaoparticle = Hkinedndeta.GetTH2(Form("htp4%s-", dataname.Data()), 5, 2, {eventclass, centselection, 1, kOPar, 1, -1, centtype, phirange});
	hkinedndetakaon->Scale(rkaon);
	hkinedndetaproton->Scale(rproton);
	hkinedndetaoparticle->Scale(ropar);
	hkinedndetapion->Add(hkinedndetakaon);
	hkinedndetapion->Add(hkinedndetaproton);
	hkinedndetapion->Add(hkinedndetaoparticle);

	Hkinedndeta.SetBin("ParticleType", {kPion - 0.5, kOPar + 0.5});
	Hkinedndeta.SetBin("Z", vtxbin);
	Hkinedndeta.SetBin("Cent", mccent);
	auto hkinedndeta = Hkinedndeta.GetTH2(Form("ht%s-", dataname.Data()), 5, 2, {eventclass, centselection, 1, 1, 1, -1, centtype, phirange});
	cout << "bkkim test : " << hkinedndeta->Integral() << endl;
	if (dopariclevariation)
		hkinedndeta = hkinedndetapion;

			// Load MC_REC (eta,z)
	auto Hmcdndeta = BSTHnSparseHelper::Load("hrecdndeta", clistmc);
	Hmcdndeta.SetBin("Cent", mccent);
	auto Hmcdndetastran = BSTHnSparseHelper::Load("hrecdndeta", clistmc);
	auto Hmcdndetabkg = BSTHnSparseHelper::Load("hrecdndeta", clistmc);
	TH2D *hmcdndetapion = nullptr;
	hmcdndetapion = Hmcdndeta.GetTH2(Form("hr1%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kPion, iptvar, -1, centtype, phirange});
	auto hmcdndetakaon = Hmcdndeta.GetTH2(Form("hr2%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kKaon, iptvar, -1, centtype, phirange});
	auto hmcdndetaproton = Hmcdndeta.GetTH2(Form("hr3%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kProtonBK, iptvar, -1, centtype, phirange});
	auto hmcdndetaoparticle = Hmcdndeta.GetTH2(Form("hr4%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kOPar, iptvar, -1, centtype, phirange});
	hmcdndetakaon->Scale(rkaon);
	hmcdndetaproton->Scale(rproton);
	hmcdndetaoparticle->Scale(ropar);
	hmcdndetapion->Add(hmcdndetakaon);
	hmcdndetapion->Add(hmcdndetaproton);
	hmcdndetapion->Add(hmcdndetaoparticle);

	Hmcdndeta.SetBin("ParticleType", {kPion - 0.5, kOPar + 0.5});
	Hmcdndeta.SetBin("Z", vtxbin);
	// order : 0evt_class 1trigger_class 2centrality 3zvtx_mcrecD 4PID 5ptvar 6eta_mcrec
	auto hmcdndeta = Hmcdndeta.GetTH2(Form("hr%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, 1, iptvar, -1, centtype, phirange});


	if (dopariclevariation)
		hmcdndeta = hmcdndetapion;


	Hmcdndetastran.SetBin("Z", vtxbin);
	Hmcdndetastran.SetBin("Cent", mccent);
	auto hmcdndetastran = Hmcdndetastran.GetTH2(Form("hrstr%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kMotherStrange, iptvar, -1, centtype, phirange});

	Hmcdndetabkg.SetBin("Z", vtxbin);
	Hmcdndetabkg.SetBin("Cent", mccent);
	auto hmcdndetabkg = Hmcdndetabkg.GetTH2(Form("hrbkg%s-", dataname.Data()), 6, 3, {eventclass, triggtype, centselection, 1, kBkg, iptvar, -1, centtype, phirange});

	//if (strangevar == kNoStrVar) hmcdndeta -> Add(hmcdndetastran,constant->GetParameter(0));
	if (strangevar == kNoStrVar && dataname.Contains("LHC15f"))
		hmcdndeta->Add(hmcdndetastran, 1);
	else if (strangevar == kNoStrVar)
		hmcdndeta->Add(hmcdndetastran, 1);
	else if (strangevar == kStrVarUp)
		hmcdndeta->Add(hmcdndetastran, 1.3);
	else if (strangevar == kStrVarDw)
		hmcdndeta->Add(hmcdndetastran, 0.7);
	//else if (strangevar == kNoStrVar) hmcdndeta -> Add(hmcdndetastran,constant->GetParameter(0));
	//else if (strangevar == kStrVarUp) hmcdndeta -> Add(hmcdndetastran,constant->GetParameter(0) * 1.1) ;
	//else if (strangevar == kStrVarDw) hmcdndeta -> Add(hmcdndetastran,constant->GetParameter(0) * 0.9) ;
	hmcdndeta->Add(hmcdndetabkg, 1.);

	// Load DATA (eta,z)
	auto Hdatadndeta = BSTHnSparseHelper::Load("hrecdndeta", clistdata);
	// order : 0evt_class 1trigger_class 2centrality 3zvtx_data 4PID 5ptvar 6eta_data
	Hdatadndeta.SetBin("Z", vtxbin);
	Hdatadndeta.SetBin("Cent", cent);
	auto hdatadndeta = Hdatadndeta.GetTH2(Form("hd%s-", dataname.Data()), 6, 3, {kDATA, triggtype, centselection, 1, kParDATA, iptvar, -1, centtype, phirange});
	cout << "It has " << hdatadndeta->GetEntries() << " entries " <<  endl;
	auto hcorrf = (TH2D *)hkinedndeta->Clone();
	//WRONG auto hcorrf = (TH2D *)hmcdndeta->Clone();
	hcorrf->Divide(hmcdndeta);
	cout << "hmcdndeta (=hrecdndeta in MC) has " << hmcdndeta->GetEntries() << " entries " <<  endl;
	cout << "hCorrf (= from hkinedndeta in MC) has " <<  hcorrf->GetEntries() << " entries " <<  endl;
	cout <<hcorrf->GetNbinsX() << " " << hcorrf->GetNbinsY() << endl;
	for (auto x = 1; x <= hcorrf->GetNbinsX(); x++)
		for (auto y = 1; y <= hcorrf->GetNbinsY(); y++)
		{
			if (hcorrf->GetBinContent(x, y) > 4)
			{
			  // cout << "bin with content > 4 " << x << " - " << y << " content: " <<  hcorrf->GetBinContent(x, y) << endl;
			  //REMOVE hcorrf->SetBinContent(x, y, 0);
			  hcorrf->SetBinContent(x, y, 0);
			  //REMOVEhcorrf->SetBinError(x, y, 0);
			  hcorrf->SetBinError(x, y, 0);
			}

			if (hcorrf->GetBinContent(x, y) > 4)
			{
				///hcorrf->SetBinContent(x, y, 4);
				//hcorrf->SetBinError(x, y, hcorrf->GetBinError(x,y)*3);
			}
		}

	//new TCanvas;
	//hcorrf -> Draw("colz");

	hdatadndeta->Multiply(hcorrf);

	auto fr = hdatadndeta->ProjectionX("_x", 0, -1, "");
	fr->Reset();
	auto letabin = hdatadndeta->GetXaxis()->FindBin(etabin[0]);
	auto retabin = hdatadndeta->GetXaxis()->FindBin(etabin[1]);
	//	cout << "letabin " << letabin << " retabin " << retabin<< endl;
	for (auto ieta = letabin; ieta <= retabin; ieta++)
	{
		auto h = hdatadndeta->ProjectionY(Form("_y%s-", dataname.Data()), ieta, ieta, "e");
		
		for (Int_t b= 1; b<= h->GetNbinsX(); b++){ 
		  //		  cout << h->GetBinContent(b) << " in bin with center " << h->GetBinCenter(b) << endl; 
		}
		
		auto lbincut = hdatadndeta->GetYaxis()->FindBin(vtxbin[0]);
		auto rbincut = hdatadndeta->GetYaxis()->FindBin(vtxbin[1]);
		auto lbin = TMath::Max(h->FindFirstBinAbove(), lbincut);
		auto rbin = TMath::Min(h->FindLastBinAbove(), rbincut);
		if (lbin > rbin) {
		  //cout << lbin << " " << rbin << endl;
		  continue;
		}
		Double_t totale = 0;
		auto total = h->IntegralAndError(lbin, rbin, totale);
		auto norm = hdataz->Integral(lbin, rbin);
		total = total / norm / hdatadndeta->GetXaxis()->GetBinWidth(ieta);
		totale = totale / norm / hdatadndeta->GetXaxis()->GetBinWidth(ieta);
		fr->SetBinContent(ieta, total);
		fr->SetBinError(ieta, totale);
		//cout << "fr -> GetBinContent " << fr->GetBinContent(ieta) << " (eta = " <<  fr->GetBinCenter(ieta) << ")" <<endl;
	}

	fr->SetLineWidth(2);
	fr->GetXaxis()->SetRangeUser(etabin[0], etabin[1]);
	fr->SetStats(false);
	fr->SetMarkerStyle(28);
	fr->SetMarkerColor(1);
	fr->SetLineColor(1);
	fr->SetMinimum(fr->GetMinimum() / 2.);
	fr->SetMaximum(fr->GetMaximum() * 2.);
	//cout<<"dndeta at eta=0"<< fr->Integral(fr->GetXaxis()->FindBin(-0.49),fr->GetXaxis()->FindBin(0.49),"width")<<endl;
	//new TCanvas;
	//fr->Draw("e");
	fr->SetName(Form("%s%s", dataname.Data(), mcname.Data()));

	//delete objects
	delete clistmc;
	delete clistdata;

	return fr;
}

void showlowbinmcresults()
{

	TString rname = "MulLHC16l";
	TString energy = "13000";
	vector<TString> centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100"};
	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};
	vector<Double_t> binwidths = {1, 4, 5, 5, 10, 10, 10, 10, 20, 30};

	TFile *fepos = TFile::Open(Form("../QuickGenEPOS/AnalysisResults_EPOSLHC_%s.root", energy.Data()), "open");
	TH1D *hn = (TH1D *)gROOT->FindObject("hNevents");
	Double_t norm = hn->GetBinContent(1);

	TH1D *epos = (TH1D *)gROOT->FindObject(Form("etaDist10"));
	epos->Scale(1. / (norm / 10.), "width");
	TH1D *eposmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	eposmb->Scale(1. / norm, "width");
	epos->Divide(eposmb);

	TFile *fcron = TFile::Open(Form("../QuickGenCRON/AnalysisResults_cron_%s.root", energy.Data()), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);
	TH1D *cron = (TH1D *)gROOT->FindObject(Form("etaDist10"));
	cron->Scale(1. / (norm / 10.), "width");
	TH1D *cronmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	cronmb->Scale(1. / norm, "width");
	cron->Divide(cronmb);

	TFile *fcroff = TFile::Open(Form("../QuickGenCROFF/AnalysisResults_croff_%s.root", energy.Data()), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);
	TH1D *croff = (TH1D *)gROOT->FindObject(Form("etaDist10"));
	croff->Scale(1. / (norm / 10.), "width");
	TH1D *croffmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	croffmb->Scale(1. / norm, "width");
	croff->Divide(croffmb);

	Filipad2 *canvas = new Filipad2(++n, 2, 0.4, 100, 50, 1, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	epos->GetXaxis()->SetRangeUser(-1.5, 1.5);
	epos->SetMaximum(1.35);
	epos->SetMinimum(0.6);
	epos->SetTitle(0);
	epos->SetStats(0);
	epos->SetLineColor(4);
	cron->SetLineColor(2);
	croff->SetLineColor(kGreen + 3);
	epos->SetLineWidth(2);
	cron->SetLineWidth(2);
	croff->SetLineWidth(2);
	hset(*epos, "#it{#eta}", mstringincl, 0.8, 0.85, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
	epos->Draw("histsame");
	cron->Draw("histsame");
	croff->Draw("histsame");

	TFile *fdata = nullptr;
	TFile *fdatamb = nullptr;
	TGraphAsymmErrors *data = nullptr;
	TGraphAsymmErrors *datamb = nullptr;

	if (!energy.Contains("SPD"))
	{
		fdata = TFile::Open(Form("../FinalResults/%s_INELg0_cent_1_40.000000_50.000000.root", rname.Data()));
		data = (TGraphAsymmErrors *)gROOT->FindObject("data");
		fdatamb = TFile::Open(Form("../FinalResults/%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
		datamb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	}
	else
	{
		fdata = TFile::Open(Form("../FinalResults/%s_INELg0_cent_4_40.000000_50.000000.root", rname.Data()));
		data = (TGraphAsymmErrors *)gROOT->FindObject("data");
		fdatamb = TFile::Open(Form("../FinalResults/%s_INELg0_cent_4_0.000000_100.000000.root", rname.Data()));
		datamb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	}
	auto r = CalcRatio(data, datamb);
	r->Draw("pe1same");

	TString energy2;
	if (energy.Contains("7000"))
		energy2 = "7";
	else if (energy.Contains("13000"))
		energy2 = "13";
	else if (energy.Contains("5020"))
		energy2 = "5.02";

	TLegend *legendtitle = new TLegend(0.187799, 0.839721, 0.549043, 0.912892, Form("ALICE, pp"));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.08);
	legendtitle->Draw();

	legendtitle = new TLegend(0.188995, 0.745645, 0.550239, 0.804878, Form("#sqrt{s} = %s TeV", energy2.Data()));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.0627178);
	legendtitle->Draw();

	TString fowardorcentral;
	if (energy.Contains("SPD"))
		fowardorcentral = "Central Mult. Class : 40-50%";
	else
		fowardorcentral = "Forward Mult. Class : 40-50%";

	TLegend *legendper = new TLegend(0.565217, 0.674797, 0.819398, 0.915796, fowardorcentral.Data(), "brNDC");
	legendper->SetTextSize(0.0522648);
	legendper->SetBorderSize(0);
	legendper->AddEntry(r, "data", "lp");
	legendper->AddEntry(cron, "PYTHIA8 Monash", "l");
	legendper->AddEntry(croff, "PYTHIA8 Monash no CR", "l");
	legendper->AddEntry(epos, "EPOS LHC", "l");
	legendper->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TH1D *dummy = (TH1D *)epos->Clone();
	dummy->Reset();
	dummy->SetMaximum(1.65);
	dummy->SetMinimum(1);
	hset(*dummy, "#it{#eta}", "Data / model", 0.7, 0.62, 0.13, 0.12, 0.01, 0.01, 0.1, 0.1, 510, 505);

	dummy->Draw();
	TGraphAsymmErrors *repos = CalcRatio(r, epos);
	repos->Draw("e1same");

	TGraphAsymmErrors *rcron = CalcRatio(r, cron);
	rcron->Draw("e1same");

	TGraphAsymmErrors *rcroff = CalcRatio(r, croff);
	rcroff->Draw("e1same");

	canvas->GetPad(0)->SaveAs(Form("model_%s_lowbin.pdf", energy.Data()));
}

void showmcresults()
{

	TString rname = "MulLHC16l";
	TString energy = "13000_SPD";
	vector<TString> centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100"};
	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};
	vector<Double_t> binwidths = {1, 4, 5, 5, 10, 10, 10, 10, 20, 30};

	TFile *fepos = TFile::Open(Form("../QuickGenEPOS/AnalysisResults_EPOSLHC_%s.root", energy.Data()), "open");
	TH1D *hn = (TH1D *)gROOT->FindObject("hNevents");
	Double_t norm = hn->GetBinContent(1);

	TH1D *epos = (TH1D *)gROOT->FindObject(Form("etaDist0"));
	TH1D *h1 = (TH1D *)gROOT->FindObject(Form("etaDist1"));
	TH1D *h2 = (TH1D *)gROOT->FindObject(Form("etaDist2"));
	TH1D *h3 = (TH1D *)gROOT->FindObject(Form("etaDist3"));
	epos->Add(h1, 1);
	epos->Add(h2, 1);
	epos->Add(h3, 1);
	epos->Scale(1. / (norm / 100.), "width");
	TH1D *eposmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	eposmb->Scale(1. / norm, "width");
	epos->Divide(eposmb);

	TFile *fcron = TFile::Open(Form("../QuickGenCRON/AnalysisResults_cron_%s.root", energy.Data()), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);
	TH1D *cron = (TH1D *)gROOT->FindObject(Form("etaDist0"));
	h1 = (TH1D *)gROOT->FindObject(Form("etaDist1"));
	h2 = (TH1D *)gROOT->FindObject(Form("etaDist2"));
	h3 = (TH1D *)gROOT->FindObject(Form("etaDist3"));
	cron->Add(h1, 1);
	cron->Add(h2, 1);
	cron->Add(h3, 1);
	cron->Scale(1. / (norm / 100.), "width");
	TH1D *cronmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	cronmb->Scale(1. / norm, "width");
	cron->Divide(cronmb);

	TFile *fcroff = TFile::Open(Form("../QuickGenCROFF/AnalysisResults_croff_%s.root", energy.Data()), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);
	TH1D *croff = (TH1D *)gROOT->FindObject(Form("etaDist0"));
	h1 = (TH1D *)gROOT->FindObject(Form("etaDist1"));
	h2 = (TH1D *)gROOT->FindObject(Form("etaDist2"));
	h3 = (TH1D *)gROOT->FindObject(Form("etaDist3"));
	croff->Add(h1, 1);
	croff->Add(h2, 1);
	croff->Add(h3, 1);
	croff->Scale(1. / (norm / 100.), "width");
	TH1D *croffmb = (TH1D *)gROOT->FindObject(Form("etaDist"));
	croffmb->Scale(1. / norm, "width");
	croff->Divide(croffmb);

	Filipad2 *canvas = new Filipad2(++n, 2, 0.4, 100, 50, 1, 1, 1);
	canvas->Draw();
	TPad *p = canvas->GetPad(1); //upper pad
	p->SetTickx();
	p->SetGridy(0);
	p->SetGridx(0);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *leg = new TLegend(0.299, 0.457, 0.925, 0.9, NULL, "brNDC");
	leg->SetTextSize(0.062);
	leg->SetBorderSize(0);

	epos->GetXaxis()->SetRangeUser(-1.5, 1.5);
	epos->SetMaximum(5.5);
	epos->SetMinimum(3.1);

	if (energy.Contains("SPD"))
	{
		epos->SetMaximum(6.5);
		epos->SetMinimum(3.7);
	}
	epos->SetTitle(0);
	epos->SetStats(0);
	epos->SetLineColor(4);
	cron->SetLineColor(2);
	croff->SetLineColor(kGreen + 3);
	epos->SetLineWidth(2);
	cron->SetLineWidth(2);
	croff->SetLineWidth(2);
	hset(*epos, "#it{#eta}", mstringincl, 0.8, 0.85, 0.08, 0.08, 0.01, 0.01, 0.07, 0.07, 510, 510);
	epos->Draw("histsame");
	cron->Draw("histsame");
	croff->Draw("histsame");

	TFile *fdata = nullptr;
	TFile *fdatamb = nullptr;
	TGraphAsymmErrors *data = nullptr;
	TGraphAsymmErrors *datamb = nullptr;

	if (!energy.Contains("SPD"))
	{
		fdata = TFile::Open(Form("../FinalResults/%s_INELg0_cent_1_0.000000_1.000000.root", rname.Data()));
		data = (TGraphAsymmErrors *)gROOT->FindObject("data");
		fdatamb = TFile::Open(Form("../FinalResults/%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
		datamb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	}
	else
	{
		fdata = TFile::Open(Form("../FinalResults/%s_INELg0_cent_4_0.000000_1.000000.root", rname.Data()));
		data = (TGraphAsymmErrors *)gROOT->FindObject("data");
		fdatamb = TFile::Open(Form("../FinalResults/%s_INELg0_cent_4_0.000000_100.000000.root", rname.Data()));
		datamb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	}
	auto r = CalcRatio(data, datamb);
	r->Draw("pe1same");

	TString energy2;
	if (energy.Contains("7000"))
		energy2 = "7";
	else if (energy.Contains("13000"))
		energy2 = "13";
	else if (energy.Contains("5020"))
		energy2 = "5.02";

	TLegend *legendtitle = new TLegend(0.187799, 0.839721, 0.549043, 0.912892, Form("ALICE, pp"));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.08);
	legendtitle->Draw();

	legendtitle = new TLegend(0.188995, 0.745645, 0.550239, 0.804878, Form("#sqrt{s} = %s TeV", energy2.Data()));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.0627178);
	legendtitle->Draw();

	TString fowardorcentral;
	if (energy.Contains("SPD"))
		fowardorcentral = "Central Mult. Class : 0-1%";
	else
		fowardorcentral = "Forward Mult. Class : 0-1%";

	TLegend *legendper = new TLegend(0.565217, 0.674797, 0.819398, 0.915796, fowardorcentral.Data(), "brNDC");
	legendper->SetTextSize(0.0522648);
	legendper->SetBorderSize(0);
	legendper->AddEntry(r, "data", "lp");
	legendper->AddEntry(cron, "PYTHIA8 Monash", "l");
	legendper->AddEntry(croff, "PYTHIA8 Monash no CR", "l");
	legendper->AddEntry(epos, "EPOS LHC", "l");
	legendper->Draw();

	p = canvas->GetPad(2); //lower pad
	p->SetTickx();
	p->SetGridy(1);
	p->SetGridx(1);
	p->SetLogy(0);
	p->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TH1D *dummy = (TH1D *)epos->Clone();
	dummy->Reset();
	dummy->SetMaximum(1.45);
	dummy->SetMinimum(0.55);
	hset(*dummy, "#it{#eta}", "Data / model", 0.7, 0.62, 0.13, 0.12, 0.01, 0.01, 0.1, 0.1, 510, 505);

	dummy->Draw();
	TGraphAsymmErrors *repos = CalcRatio(r, epos);
	repos->Draw("e1same");

	TGraphAsymmErrors *rcron = CalcRatio(r, cron);
	rcron->Draw("e1same");

	TGraphAsymmErrors *rcroff = CalcRatio(r, croff);
	rcroff->Draw("e1same");

	canvas->GetPad(0)->SaveAs(Form("model_%s.pdf", energy.Data()));
}

void show5results(int option = 1)
{

	new TCanvas;
	setpad(gPad);

	TString rname = "MulLHC17p";

	vector<TString> filenames;

	if (option == 1)
		filenames = {
			Form("%s_INELg0_cent_1_0.000000_0.010000.root", rname.Data()), Form("%s_INELg0_cent_1_0.010000_0.050000.root", rname.Data()), Form("%s_INELg0_cent_1_0.050000_0.100000.root", rname.Data()), Form("%s_INELg0_cent_1_0.100000_0.500000.root", rname.Data()), Form("%s_INELg0_cent_1_0.500000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_1_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_1_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_1_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_1_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_1_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_1_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_1_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_1_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_1_70.000000_100.000000.root", rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_1.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_5.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_100.000000.root",rname.Data())
		};

	else if (option == 4)
		filenames = {
			Form("%s_INELg0_cent_4_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_4_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_4_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_4_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_4_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_4_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_4_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_4_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_4_70.000000_100.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_100.000000.root", rname.Data())};

	vector<TString> centtemp;
	if (option == 1)
		centtemp = {"0", "0.01", "0.05", "0.1", "0.5", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "1", "5", "100"};
	else if (option == 4)
		centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "5", "100"};

	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

	gStyle->SetTextFont(42);
	TH1D *dummy = new TH1D("dummy", "", 36, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstring, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(35);
	dummy->Draw();
	int nn = 0;
	TLatex latex;
	latex.SetTextSize(0.025);
	latex.SetTextFont(42);
	latex.SetTextAlign(21);
	latex.SetTextSize(0.04);
	//latex.DrawLatex(1.3,56,"Multiplicity bin");
	//latex.DrawLatex(1.3,53,"(SPD %)");
	TLegend *legendper = new TLegend(0.815186, 0.141053, 0.982808, 0.922105, "Multiplicity (%)", "brNDC");
	legendper->SetTextSize(0.04);
	legendper->SetBorderSize(0);
	legendper->SetTextAlign(22);
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("datacorr");
		h->SetFillColorAlpha(colortemp[nn], 0.5);
		h->SetMarkerStyle(7);
		h->SetMarkerColor(colortemp[nn]);
		h->GetXaxis()->SetRangeUser(-1.7, 1.7);
		h->Draw("p5same");
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		legendper->AddEntry(h, Form(" %s - %s ", centtemp.at(nn).Data(), centtemp.at(nn + 1).Data()), "pe3f");
		//cout<<"cent : "<<centtemp[nn]<<" - "<<centtemp[nn+1]<<" "<< " "<<dndeta<<"+"<<eh<<"-"<<el<<endl;
		//cout<<"cent : "<<centtemp[nn]<<" - "<<centtemp[nn+1]<<" "<< " "<<dndeta10/2<<"+"<<eh10/2<<"-"<<el10/2<< " in |eta|<1 " <<endl;
		printf(" %s - %s  %.2f\\pm%.2f  \n", centtemp[nn].Data(), centtemp[nn + 1].Data(), dndeta, (eh + el) / 2.);
		nn++;
	}

	nn = 0;
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		//printf("AddPoint(cent[%d], 5020, %4.2f,0,0,%.2f,%.2f, rcent[%d], 5.49, 0.08, 0.06); // %s - %s  \n",nn,dndeta,eh,el,nn, centtemp[nn].Data(),centtemp[nn+1].Data());
		printf(" %s - %s  %.2f\\pm%.2f  ", centtemp[nn].Data(), centtemp[nn + 1].Data(), dndeta, (eh + el) / 2.);
		nn++;
	}

	legendper->Draw();
	latex.SetTextSize(0.056);
	latex.SetTextColor(kBlack);
	latex.SetTextAlign(22);
	latex.DrawLatex(0, 56, "pp, #sqrt{#it{s}} = 5.02 TeV");
	latex.DrawLatex(0.72, 53, "ALICE");
	TLegend *lenergy = new TLegend(0.3, 0.87, 0.7, 0.95, "pp, #sqrt{#it{s}} = 5.02 TeV, ALICE", "brNDC");
	lenergy->SetTextSize(0.06);
	lenergy->SetBorderSize(0);
	lenergy->SetTextAlign(22);
	lenergy->Draw();

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/%s_%d_%s_allcent_5TeV.pdf", rname.Data(), option, title.Data()));

	new TCanvas;
	setpad(gPad);

	dummy = new TH1D("dummy", "", 30, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstringincl, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(6.05);
	dummy->Draw();

	nn = 0;
	TString mbfilename(Form("%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
	TFile *fmb = TFile::Open(Form("../FinalResults/%s", mbfilename.Data()));
	auto mb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		auto r = CalcRatio(h, mb);
		r->SetFillColorAlpha(colortemp[nn], 0.5);
		r->SetMarkerStyle(7);
		r->SetMarkerColor(colortemp[nn]);
		r->GetXaxis()->SetRangeUser(-1.5, 1.5);
		r->Draw("p5same");
		nn++;
	}

	legendper->Draw();
	lenergy->Draw();

	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC17p_%d_%s_rallcent_5TeV.pdf", centtype, title.Data()));
}

void show7results(int option = 1)
{

	new TCanvas;
	setpad(gPad);

	TString rname = "MulLHC10d";

	vector<TString> filenames;

	if (option == 1)
		filenames = {
			Form("%s_INELg0_cent_1_0.000000_0.010000.root", rname.Data()), Form("%s_INELg0_cent_1_0.010000_0.100000.root", rname.Data()), Form("%s_INELg0_cent_1_0.100000_0.500000.root", rname.Data()), Form("%s_INELg0_cent_1_0.500000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_1_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_1_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_1_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_1_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_1_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_1_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_1_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_1_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_1_70.000000_100.000000.root", rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_1.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_5.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_100.000000.root",rname.Data())
		};

	else if (option == 4)
		filenames = {
			Form("%s_INELg0_cent_4_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_4_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_4_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_4_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_4_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_4_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_4_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_4_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_4_70.000000_100.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_100.000000.root", rname.Data())};

	vector<TString> centtemp;
	if (option == 1)
		centtemp = {"0", "0.01", "0.1", "0.5", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "1", "5", "100"};
	else if (option == 4)
		centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "5", "100"};
	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

	TH1D *dummy = new TH1D("dummy", "", 37, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstring, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(42);
	dummy->Draw();
	int nn = 0;
	TLatex latex;
	latex.SetTextSize(0.025);
	latex.SetTextFont(42);
	latex.SetTextAlign(21);
	latex.SetTextSize(0.04);
	//latex.DrawLatex(1.3,56,"Multiplicity bin");
	//latex.DrawLatex(1.3,53,"(SPD %)");
	TLegend *legendper = new TLegend(0.815186, 0.141053, 0.982808, 0.922105, "Multiplicity (%)", "brNDC");
	legendper->SetTextSize(0.04);
	legendper->SetBorderSize(0);
	legendper->SetTextAlign(22);
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		auto hcorr = (TGraphAsymmErrors *)gROOT->FindObject("datacorr");
		auto hall = (TGraphAsymmErrors *)gROOT->FindObject("dataall");
		h->SetFillColorAlpha(colortemp[nn], 0.5);
		hall->SetFillColorAlpha(colortemp[nn], 0.5);
		//hcorr->SetFillColorAlpha(colortemp[nn], 0.5);
		h->SetMarkerStyle(7);
		h->SetMarkerColor(colortemp[nn]);
		hall->SetMarkerStyle(7);
		hall->SetMarkerColor(colortemp[nn]);
		hcorr->SetMarkerColor(colortemp[nn]);
		hcorr->SetMarkerStyle(7);
		hcorr->SetFillStyle(4050);
		hcorr->SetFillColor(0);
		h->GetXaxis()->SetRangeUser(-1.5, 1.5);
		hall->GetXaxis()->SetRangeUser(-1.5, 1.5);
		//hcorr->Draw("5same");
		hall->Draw("p5same");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
		}

		legendper->AddEntry(h, Form(" %s - %s ", centtemp.at(nn).Data(), centtemp.at(nn + 1).Data()), "pe3f");
		//latex.SetTextColor(colortemp[nn]);
		//latex.DrawLatex(1.3,foundy,Form(" %s - %s",centtemp.at(nn).Data(),centtemp.at(nn+1).Data()));
		cout << "cent : " << centtemp[nn] << " - " << centtemp[nn + 1] << " "
			 << " " << dndeta << " + " << eh << " - " << el << endl;
		nn++;
	}
	nn = 0;
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		//printf("AddPoint(cent[%d], 7000, %4.2f,0,0,%.2f,%.2f); // %s - %s  \n",nn,dndeta,eh,el,centtemp[nn].Data(),centtemp[nn+1].Data());
		//printf("AddPoint(cent[%d], 7000, %4.2f,0,0,%.2f,%.2f, rcent[%d], 5.91, 0.07, 0.05); // %s - %s  \n",nn,dndeta,eh,el,nn, centtemp[nn].Data(),centtemp[nn+1].Data());
		printf(" %s - %s  %.2f\\pm%.2f  \n", centtemp[nn].Data(), centtemp[nn + 1].Data(), dndeta, (eh + el) / 2.);
		nn++;
	}

	legendper->Draw();
	latex.SetTextSize(0.056);
	latex.SetTextColor(kBlack);
	latex.SetTextAlign(22);
	latex.DrawLatex(0, 56, "pp, #sqrt{#it{s}} = 7 TeV");
	latex.DrawLatex(0.72, 53, "ALICE");
	TLegend *lenergy = new TLegend(0.3, 0.87, 0.7, 0.95, "pp, #sqrt{#it{s}} = 7 TeV, ALICE", "brNDC");
	lenergy->SetTextSize(0.06);
	lenergy->SetBorderSize(0);
	lenergy->SetTextAlign(22);
	lenergy->Draw();

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC10d_%d_%s_allcent_7TeV.pdf", 1, title.Data()));

	new TCanvas;
	setpad(gPad);

	dummy = new TH1D("dummy", "", 30, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstringincl, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(6.5);
	dummy->Draw();

	nn = 0;
	TString mbfilename(Form("%s_INELg0_cent_1_0.000000_100.000000.root", "MulLHC10d"));
	//TString mbfilename(Form("%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
	TFile *fmb = TFile::Open(Form("../FinalResults/%s", mbfilename.Data()));
	auto mb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		auto r = CalcRatio(h, mb);
		r->SetFillColorAlpha(colortemp[nn], 0.5);
		r->SetMarkerStyle(7);
		r->SetMarkerColor(colortemp[nn]);
		r->GetXaxis()->SetRangeUser(-1.5, 1.5);
		r->Draw("p5same");
		nn++;
	}

	legendper->Draw();
	lenergy->Draw();

	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC10d_%d_%s_rallcent_7TeV.pdf", centtype, title.Data()));
}
void show8results(int option = 1)
{

	new TCanvas;
	setpad(gPad);

	TString rname = "MulLHC12h";

	vector<TString> filenames;

	if (option == 1)
		filenames = {
			Form("%s_INELg0_cent_1_0.000000_0.010000.root", rname.Data()), Form("%s_INELg0_cent_1_0.010000_0.100000.root", rname.Data()), Form("%s_INELg0_cent_1_0.100000_0.500000.root", rname.Data()), Form("%s_INELg0_cent_1_0.500000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_1_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_1_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_1_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_1_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_1_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_1_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_1_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_1_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_1_70.000000_100.000000.root", rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_1.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_5.000000.root",rname.Data())
			//,Form("%s_INELg0_cent_1_0.000000_100.000000.root",rname.Data())
		};

	else if (option == 4)
		filenames = {
			Form("%s_INELg0_cent_4_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_4_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_4_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_4_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_4_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_4_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_4_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_4_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_4_70.000000_100.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_0.000000_100.000000.root", rname.Data())};

	vector<TString> centtemp;
	if (option == 1)
		centtemp = {"0", "0.01", "0.1", "0.5", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "1", "5", "100"};
	else if (option == 4)
		centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100", "5", "100"};
	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

	TH1D *dummy = new TH1D("dummy", "", 37, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstring, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(42);
	dummy->Draw();
	int nn = 0;
	TLatex latex;
	latex.SetTextSize(0.025);
	latex.SetTextFont(42);
	latex.SetTextAlign(21);
	latex.SetTextSize(0.04);
	//latex.DrawLatex(1.3,56,"Multiplicity bin");
	//latex.DrawLatex(1.3,53,"(SPD %)");
	TLegend *legendper = new TLegend(0.815186, 0.141053, 0.982808, 0.922105, "Multiplicity (%)", "brNDC");
	legendper->SetTextSize(0.04);
	legendper->SetBorderSize(0);
	legendper->SetTextAlign(22);
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		h->SetFillColorAlpha(colortemp[nn], 0.5);
		h->SetMarkerStyle(7);
		h->SetMarkerColor(colortemp[nn]);
		h->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h->Draw("p5same");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
		}

		legendper->AddEntry(h, Form(" %s - %s ", centtemp.at(nn).Data(), centtemp.at(nn + 1).Data()), "pe3f");
		//latex.SetTextColor(colortemp[nn]);
		//latex.DrawLatex(1.3,foundy,Form(" %s - %s",centtemp.at(nn).Data(),centtemp.at(nn+1).Data()));
		cout << "cent : " << centtemp[nn] << " - " << centtemp[nn + 1] << " "
			 << " " << dndeta << " + " << eh << " - " << el << endl;
		nn++;
	}
	nn = 0;
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		//printf("AddPoint(cent[%d], 7000, %4.2f,0,0,%.2f,%.2f); // %s - %s  \n",nn,dndeta,eh,el,centtemp[nn].Data(),centtemp[nn+1].Data());
		//printf("AddPoint(cent[%d], 7000, %4.2f,0,0,%.2f,%.2f, rcent[%d], 5.91, 0.07, 0.05); // %s - %s  \n",nn,dndeta,eh,el,nn, centtemp[nn].Data(),centtemp[nn+1].Data());
		printf(" %s - %s  %.2f\\pm%.2f  \n", centtemp[nn].Data(), centtemp[nn + 1].Data(), dndeta, (eh + el) / 2.);
		nn++;
	}

	legendper->Draw();
	latex.SetTextSize(0.056);
	latex.SetTextColor(kBlack);
	latex.SetTextAlign(22);
	latex.DrawLatex(0, 56, "pp, #sqrt{#it{s}} = 7 TeV");
	latex.DrawLatex(0.72, 53, "ALICE");
	TLegend *lenergy = new TLegend(0.3, 0.87, 0.7, 0.95, "pp, #sqrt{#it{s}} = 8 TeV, ALICE", "brNDC");
	lenergy->SetTextSize(0.06);
	lenergy->SetBorderSize(0);
	lenergy->SetTextAlign(22);
	lenergy->Draw();

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC10d_%d_%s_allcent_7TeV.pdf", 1, title.Data()));

	new TCanvas;
	setpad(gPad);

	dummy = new TH1D("dummy", "", 30, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstringincl, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(6.5);
	dummy->Draw();

	nn = 0;
	TString mbfilename(Form("%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
	TFile *fmb = TFile::Open(Form("../FinalResults/%s", mbfilename.Data()));
	auto mb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		auto r = CalcRatio(h, mb);
		r->SetFillColorAlpha(colortemp[nn], 0.5);
		r->SetMarkerStyle(7);
		r->SetMarkerColor(colortemp[nn]);
		r->GetXaxis()->SetRangeUser(-1.5, 1.5);
		r->Draw("p5same");
		nn++;
	}

	legendper->Draw();
	lenergy->Draw();

	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC10d_%d_%s_rallcent_7TeV.pdf", centtype, title.Data()));
}

void show7results()
{

	new TCanvas;
	setpad(gPad);

	vector<TString> filenames = {
		"MulLHC10d_INELg0_cent_4_0.000000_1.000000.root", "MulLHC10d_INELg0_cent_4_1.000000_5.000000.root", "MulLHC10d_INELg0_cent_4_5.000000_10.000000.root", "MulLHC10d_INELg0_cent_4_10.000000_15.000000.root", "MulLHC10d_INELg0_cent_4_15.000000_20.000000.root", "MulLHC10d_INELg0_cent_4_20.000000_30.000000.root", "MulLHC10d_INELg0_cent_4_30.000000_40.000000.root", "MulLHC10d_INELg0_cent_4_40.000000_50.000000.root", "MulLHC10d_INELg0_cent_4_50.000000_70.000000.root", "MulLHC10d_INELg0_cent_4_70.000000_100.000000.root"};

	//vector<TString> centtemp={"0","0.01","0.1","0.5","1","5","10","15","20","30","40","50","70","100"};
	vector<TString> centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100"};
	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

	TH1D *dummy = new TH1D("dummy", "", 37, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstring, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(37);
	dummy->Draw();
	int nn = 0;
	TLatex latex;
	latex.SetTextSize(0.025);
	latex.SetTextFont(42);
	latex.SetTextAlign(21);
	latex.SetTextSize(0.04);
	//latex.DrawLatex(1.3,56,"Multiplicity bin");
	//latex.DrawLatex(1.3,53,"(SPD %)");
	TLegend *legendper = new TLegend(0.815186, 0.141053, 0.982808, 0.922105, "Multiplicity (%)", "brNDC");
	legendper->SetTextSize(0.04);
	legendper->SetBorderSize(0);
	legendper->SetTextAlign(22);
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		h->SetFillColorAlpha(colortemp[nn], 0.5);
		h->SetMarkerStyle(7);
		h->SetMarkerColor(colortemp[nn]);
		h->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h->Draw("p5same");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		legendper->AddEntry(h, Form(" %s - %s ", centtemp.at(nn).Data(), centtemp.at(nn + 1).Data()), "pe3f");
		//latex.SetTextColor(colortemp[nn]);
		//latex.DrawLatex(1.3,foundy,Form(" %s - %s",centtemp.at(nn).Data(),centtemp.at(nn+1).Data()));
		cout << "cent : " << centtemp[nn] << " - " << centtemp[nn + 1] << " "
			 << " " << dndeta << " + " << eh << " - " << el << endl;
		cout << "cent : " << centtemp[nn] << " - " << centtemp[nn + 1] << " "
			 << " " << dndeta10 / 2 << "+" << eh10 / 2 << "-" << el10 / 2 << " in |eta|<1 " << endl;
		nn++;
	}

	nn = 0;
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtempl[nn].Data(),centtemph[nn].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		printf("AddPoint(cent[%d], 7000, %4.2f,0,0,%.2f,%.2f); // %s - %s  \n", nn, dndeta, eh, el, centtemp[nn].Data(), centtemp[nn + 1].Data());
		nn++;
	}

	legendper->Draw();
	latex.SetTextSize(0.056);
	latex.SetTextColor(kBlack);
	latex.SetTextAlign(22);
	latex.DrawLatex(0, 56, "pp, #sqrt{#it{s}} = 7 TeV");
	latex.DrawLatex(0.72, 53, "ALICE");
	TLegend *lenergy = new TLegend(0.3, 0.87, 0.7, 0.95, "pp, #sqrt{#it{s}} = 7 TeV, ALICE", "brNDC");
	lenergy->SetTextSize(0.06);
	lenergy->SetBorderSize(0);
	lenergy->SetTextAlign(22);
	lenergy->Draw();

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC10d_%d_%s_allcent_7TeV.pdf", 1, title.Data()));
}

void show13results(int option = 1)
{

	new TCanvas;
	setpad(gPad);
	TString rname = "MulLHC16l";
	vector<TString> filenames;

	if (option == 1)
		filenames = {
			//Form("%s_INELg0_cent_1_0.000000_0.010000.root",rname.Data())
			Form("%s_INELg0_cent_1_0.000000_0.010000.root", rname.Data()), Form("%s_INELg0_cent_1_0.010000_0.100000.root", rname.Data()), Form("%s_INELg0_cent_1_0.100000_0.500000.root", rname.Data()), Form("%s_INELg0_cent_1_0.500000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_1_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_1_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_1_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_1_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_1_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_1_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_1_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_1_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_1_70.000000_100.000000.root", rname.Data())};

	else if (option == 2)
		filenames = {
			Form("%s_INELg0_cent_2_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_2_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_2_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_2_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_2_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_2_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_2_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_2_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_2_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_2_70.000000_100.000000.root", rname.Data())};

	else if (option == 4)
		filenames = {
			Form("%s_INELg0_cent_4_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_4_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_4_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_4_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_4_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_4_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_4_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_4_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_4_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_4_70.000000_100.000000.root", rname.Data())};
	else if (option == 5)
		filenames = {
			Form("%s_INELg0_cent_5_0.000000_1.000000.root", rname.Data()), Form("%s_INELg0_cent_5_1.000000_5.000000.root", rname.Data()), Form("%s_INELg0_cent_5_5.000000_10.000000.root", rname.Data()), Form("%s_INELg0_cent_5_10.000000_15.000000.root", rname.Data()), Form("%s_INELg0_cent_5_15.000000_20.000000.root", rname.Data()), Form("%s_INELg0_cent_5_20.000000_30.000000.root", rname.Data()), Form("%s_INELg0_cent_5_30.000000_40.000000.root", rname.Data()), Form("%s_INELg0_cent_5_40.000000_50.000000.root", rname.Data()), Form("%s_INELg0_cent_5_50.000000_70.000000.root", rname.Data()), Form("%s_INELg0_cent_5_70.000000_100.000000.root", rname.Data())};

	vector<TString> centtemp;
	vector<TString> centv;

	if (option == 1)
		centtemp = {"0", "0.01", "0.1", "0.5", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100"};
	else
		centtemp = {"0", "1", "5", "10", "15", "20", "30", "40", "50", "70", "100"};

	vector<int> colortemp = {kRed, kOrange, kGreen, kGreen + 3, kBlue, kBlue + 2, kViolet, kCyan + 2, kAzure, kMagenta, kMagenta + 3, kYellow + 3, kYellow + 2, kSpring - 7};

	TH1D *dummy = new TH1D("dummy", "", 50, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstring, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	//dummy->SetMaximum(50);
	dummy->SetMaximum(50);
	dummy->Draw();
	int nn = 0;
	TLatex latex;
	latex.SetTextSize(0.025);
	latex.SetTextFont(42);
	latex.SetTextAlign(21);
	latex.SetTextSize(0.04);
	TLegend *legendper = new TLegend(0.815186, 0.141053, 0.982808, 0.922105, "Multiplicity (%)", "brNDC");
	legendper->SetTextSize(0.04);
	legendper->SetBorderSize(0);
	legendper->SetTextAlign(22);

	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		h->SetFillColorAlpha(colortemp[nn], 0.5);
		h->SetMarkerStyle(7);
		h->SetMarkerColor(colortemp[nn]);
		h->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h->Draw("p5same");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		latex.SetTextColor(colortemp[nn]);
		legendper->AddEntry(h, Form(" %s - %s ", centtemp.at(nn).Data(), centtemp.at(nn + 1).Data()), "pe3f");
		//cout<<"cent : "<<centtemp[nn]<<" - "<<centtemp[nn+1]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//cout<<"cent : "<<centtemp[nn]<<" - "<<centtemp[nn+1]<<" "<< " "<<dndeta10/2<<"+"<<eh10/2<<"-"<<el10/2<< " in |eta|<1 " <<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		printf(" %s - %s  %.2f\\pm%.2f  \n", centtemp[nn].Data(), centtemp[nn + 1].Data(), dndeta, (eh + el) / 2.);
		nn++;
	}

	nn = 0;
	TGraphErrors *dataerror = new TGraphErrors;
	Double1D datadndeta;
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		double eh10 = 0, el10 = 0, dndeta10 = 0;
		for (int k = 0; k < h->GetN(); k++)
		{
			double x, y;
			h->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += h->GetErrorYhigh(k) * 0.1;
				el += h->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
			if (x > -1 && x < 1)
			{
				eh10 += h->GetErrorYhigh(k) * 0.1;
				el10 += h->GetErrorYlow(k) * 0.1;
				dndeta10 += y * 0.1;
			}
		}

		datadndeta.push_back(dndeta);
		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		//printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		//printf("AddPoint(cent[%d], 13000, %4.2f,0,0,%.2f,%.2f); // %s - %s  \n",nn,dndeta,eh,el,centtemp[nn].Data(),centtemp[nn+1].Data());
		printf("AddPoint(cent[%d], 13000, %4.2f,0,0,%.2f,%.2f, rcent[%d], 6.93, 0.10, 0.08); // %s - %s  \n", nn, dndeta, eh, el, nn, centtemp[nn].Data(), centtemp[nn + 1].Data());
		//AddPoint (dataerror, dndeta, dndeta, 0, (eh+el/2.));

		nn++;
	}

	legendper->Draw();
	latex.SetTextSize(0.056);
	latex.SetTextColor(kBlack);
	latex.SetTextAlign(22);
	latex.DrawLatex(0, 56, "pp, #sqrt{#it{s}} = 13 TeV");
	//latex.DrawLatex(0.72,53,"ALICE");
	TLegend *lenergy = new TLegend(0.3, 0.87, 0.7, 0.95, "pp, #sqrt{#it{s}} = 13TeV, ALICE", "brNDC");
	lenergy->SetTextSize(0.06);
	lenergy->SetBorderSize(0);
	lenergy->SetTextAlign(22);
	lenergy->Draw();

	TString title = "INELg0";
	if (eventclass == kINEL)
		title = "INEL";
	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC16l_%d_%s_allcent_13TeV.pdf", 1, title.Data()));

	new TCanvas;
	setpad(gPad);

	dummy = new TH1D("dummy", "", 30, -1.5, 1.5);
	hset(*dummy, "#it{#eta}", mstringincl, 0.9, 0.9, 0.07, 0.07, 0.01, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(7.1);
	dummy->Draw();

	nn = 0;
	TString mbfilename(Form("%s_INELg0_cent_1_0.000000_100.000000.root", rname.Data()));
	TFile *fmb = TFile::Open(Form("../FinalResults/%s", mbfilename.Data()));
	auto mb = (TGraphAsymmErrors *)gROOT->FindObject("data");
	for (auto name : filenames)
	{
		TFile *f = TFile::Open(Form("../FinalResults/%s", name.Data()));
		auto h = (TGraphAsymmErrors *)gROOT->FindObject("data");
		auto r = CalcRatio(h, mb);

		Double_t foundy = 0;
		double eh = 0, el = 0, dndeta = 0;
		for (int k = 0; k < r->GetN(); k++)
		{
			double x, y;
			r->GetPoint(k, x, y);
			if (abs(x - 0.05) < 0.01)
				foundy = y;
			if (x > -0.5 && x < 0.5)
			{
				eh += r->GetErrorYhigh(k) * 0.1;
				el += r->GetErrorYlow(k) * 0.1;
				dndeta += y * 0.1;
			}
		}

		double lbin = centtemp[nn].Atof();
		double rbin = centtemp[nn + 1].Atof();
		double mbin = (lbin + rbin / 2);
		AddPoint(dataerror, datadndeta.at(nn), dndeta, 0, (eh + el / 2.));

		r->SetFillColorAlpha(colortemp[nn], 0.5);
		r->SetMarkerStyle(7);
		r->SetMarkerColor(colortemp[nn]);
		r->GetXaxis()->SetRangeUser(-1.5, 1.5);
		r->Draw("p5same");
		nn++;
	}

	legendper->Draw();
	lenergy->Draw();

	if (SaveHistogramsAndResults)
		gPad->GetCanvas()->SaveAs(Form("../Figures/MulLHC16l_%d_%s_rallcent_13TeV.pdf", centtype, title.Data()));

	new TCanvas;
	setpad2(gPad);
	gPad->SetLogx(0);

	dummy = new TH1D("dummy", "", 105000, 0, 110);
	//const char *fractional = "#sigma/#sigma_{MB_{AND>0}}, #sigma/#sigma_{INEL_{>0}} (%)";
	const char *fractional = "#sigma/#sigma_{MB_{AND>0}}, #sigma/#sigma_{INEL_{>0}} (%)";
	hset(*dummy, avmstring, avmstringincl, 0.9, 0.9, 0.07, 0.07, -0.002, 0.01, 0.06, 0.06, 510, 510);
	dummy->SetMinimum(0);
	dummy->SetMaximum(9.2);
	//if (option==4) dummy->SetMaximum(6.5);
	dummy->GetXaxis()->SetRangeUser(0.003, 69);
	//if (option ==4) dummy->GetXaxis()->SetRangeUser(0.3,110);
	dummy->Draw();

	dataerror->SetMarkerStyle(20);
	dataerror->Draw("p1same");
	dataerror->SetMinimum(0);
	dataerror->SetMaximum(65);
	//ReverseXAxis(dataerror);
	//ReverseXGraph(dataerror);

	TLine *line = new TLine(0, 0, 5.5, 5.5);
	line->SetLineStyle(2);
	line->SetLineColor(2);
	line->SetLineWidth(2);
	line->SetLineColor(kBlack);

	//line -> Draw();

	TFile *fepos = TFile::Open(Form("../QuickGenEPOS/AnalysisResults_EPOSLHC_13000.root"), "open");
	if (option == 4)
		fepos = TFile::Open(Form("../QuickGenEPOS/AnalysisResults_EPOSLHC_13000_SPD.root"), "open");
	TH1D *hn = (TH1D *)gROOT->FindObject("hNevents");
	Double_t norm = hn->GetBinContent(1);
	vector<Double_t> binwidths = {0.01, 0.09, 0.4, 0.5, 4, 5, 5, 5, 10, 10, 10, 20, 30};
	if (option == 4)
		binwidths = {1, 4, 5, 5, 5, 10, 10, 10, 20, 30};

	Int1D binorder = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
	if (option == 4)
		binorder = {0, 4, 5, 6, 7, 8, 9, 10, 11, 12};

	TH1D *hmb = (TH1D *)gROOT->FindObject("etaDist");
	hmb->Scale(1. / norm, "width");
	double dndetamb = hmb->Integral(hmb->GetXaxis()->FindBin(-0.49), hmb->GetXaxis()->FindBin(0.49), "width");

	TGraphErrors *epose = new TGraphErrors;
	n = 0;
	for (auto i : binorder)
	{
		auto h = (TH1D *)gROOT->FindObject(Form("etaDist%d", i));
		if (option == 4 && i == 0)
		{
			auto h1 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 1));
			auto h2 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 2));
			auto h3 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 3));
			h->Add(h1);
			h->Add(h2);
			h->Add(h3);
			h->Scale(1. / (norm / 100.), "width");
		}
		else
			h->Scale(1. / (norm / 100. * binwidths.at(n)), "width");
		double dndeta = h->Integral(h->GetXaxis()->FindBin(-0.49), h->GetXaxis()->FindBin(0.49), "width");
		//AddPoint(epose,datadndeta.at(n),dndeta/dndetamb,0,0);
		AddPoint(epose, dndeta, dndeta / dndetamb, 0, 0);
		n++;
	}
	epose->SetLineColor(4);
	epose->SetLineWidth(5);
	epose->Draw("lsame");

	TFile *fcron = TFile::Open(Form("../QuickGenCRON/AnalysisResults_cron_13000.root"), "open");
	if (option == 4)
		fcron = TFile::Open(Form("../QuickGenCRON/AnalysisResults_cron_13000_SPD.root"), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);

	hmb = (TH1D *)gROOT->FindObject("etaDist");
	hmb->Scale(1. / norm, "width");
	dndetamb = hmb->Integral(hmb->GetXaxis()->FindBin(-0.49), hmb->GetXaxis()->FindBin(0.49), "width");

	TGraphErrors *crone = new TGraphErrors;
	n = 0;
	for (auto i : binorder)
	{
		auto h = (TH1D *)gROOT->FindObject(Form("etaDist%d", i));
		if (option == 4 && i == 0)
		{
			auto h1 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 1));
			auto h2 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 2));
			auto h3 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 3));
			h->Add(h1);
			h->Add(h2);
			h->Add(h3);
			h->Scale(1. / (norm / 100.), "width");
		}
		else
			h->Scale(1. / (norm / 100. * binwidths.at(n)), "width");

		double dndeta = h->Integral(h->GetXaxis()->FindBin(-0.49), h->GetXaxis()->FindBin(0.49), "width");
		//AddPoint(crone,datadndeta.at(n),dndeta/dndetamb,0,0);
		AddPoint(crone, dndeta, dndeta / dndetamb, 0, 0);
		n++;
	}
	crone->SetLineColor(2);
	crone->SetLineWidth(7);
	crone->SetLineStyle(7);
	crone->Draw("lsame");
	dataerror->Draw("p1same");

	TFile *fcroff = TFile::Open(Form("../QuickGenCROFF/AnalysisResults_croff_13000.root"), "open");
	if (option == 4)
		fcroff = TFile::Open(Form("../QuickGenCROFF/AnalysisResults_croff_13000_SPD.root"), "open");
	hn = (TH1D *)gROOT->FindObject("hNevents");
	norm = hn->GetBinContent(1);

	hmb = (TH1D *)gROOT->FindObject("etaDist");
	hmb->Scale(1. / norm, "width");
	dndetamb = hmb->Integral(hmb->GetXaxis()->FindBin(-0.49), hmb->GetXaxis()->FindBin(0.49), "width");

	TGraphErrors *croffe = new TGraphErrors;
	n = 0;
	for (auto i : binorder)
	{
		auto h = (TH1D *)gROOT->FindObject(Form("etaDist%d", i));
		if (option == 4 && i == 0)
		{
			auto h1 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 1));
			auto h2 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 2));
			auto h3 = (TH1D *)gROOT->FindObject(Form("etaDist%d", 3));
			h->Add(h1);
			h->Add(h2);
			h->Add(h3);
			h->Scale(1. / (norm / 100.), "width");
		}
		else
			h->Scale(1. / (norm / 100. * binwidths.at(n)), "width");
		double dndeta = h->Integral(h->GetXaxis()->FindBin(-0.49), h->GetXaxis()->FindBin(0.49), "width");
		//AddPoint(croffe,datadndeta.at(n),dndeta/dndetamb,0,0);
		AddPoint(croffe, dndeta, dndeta / dndetamb, 0, 0);
		n++;
	}
	croffe->SetLineColor(kGreen + 3);
	croffe->SetLineWidth(4);
	croffe->SetLineStyle(3);
	croffe->Draw("lsame");
	/*

	TFile * frope = TFile::Open(Form("../QuickGenRope/AnalysisResults_rope_13000.root"),"open");
	if (option==4) fcroff = TFile::Open(Form("../QuickGenCROFF/AnalysisResults_croff_13000_SPD.root"),"open");
	hn = (TH1D*) gROOT -> FindObject("hMultiPythia");
	norm = hn -> GetEntries();


	hmb = (TH1D*) gROOT->FindObject("etaDist"); 
	hmb -> Scale(1./norm,"width");
	dndetamb = hmb -> Integral(hmb->GetXaxis()->FindBin(-0.49),hmb->GetXaxis()->FindBin(0.49),"width");	

	TGraphErrors* ropee = new TGraphErrors;
	n =0;
	for (auto i : binorder){
		auto h =  (TH1D*) gROOT->FindObject(Form("etaDist%d",i));
		if ( option == 4 && i ==0 ){
			auto h1 = (TH1D*) gROOT->FindObject(Form("etaDist%d",1));
			auto h2 = (TH1D*) gROOT->FindObject(Form("etaDist%d",2));
			auto h3 = (TH1D*) gROOT->FindObject(Form("etaDist%d",3));
			h -> Add(h1);
			h -> Add(h2);
			h -> Add(h3);
			h -> Scale(1./(norm/100.), "width");

		}
		else h -> Scale(1./(norm/100.*binwidths.at(n)), "width");
		double dndeta = h -> Integral(h->GetXaxis()->FindBin(-0.49),h->GetXaxis()->FindBin(0.49),"width");	
		AddPoint(ropee,dndeta,dndeta/dndetamb,0,0);
		n++;
	}
	ropee -> SetLineColor(kRed+3);
	ropee -> SetLineWidth(2);*/
	//ropee -> Draw("lsame");

	TLegend *legendtitle = new TLegend(0.673324, 0.873418, 0.925473, 0.947257, Form("ALICE, pp"));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.07);
	legendtitle->Draw();

	legendtitle = new TLegend(0.678367, 0.793249, 0.90043, 0.835443, Form("#sqrt{s} = 13 TeV"));
	legendtitle->SetFillColor(0);
	legendtitle->SetBorderSize(0);
	legendtitle->SetTextSize(0.06);
	legendtitle->Draw();

	TString fowardorcentral;
	if (option == 4)
		fowardorcentral = "Central Mult. Class";
	else
		fowardorcentral = "Forward Mult. Class";

	legendper = new TLegend(0.179083, 0.675789, 0.43553, 0.949474, fowardorcentral.Data(), "brNDC");
	legendper->SetTextSize(0.0421053);
	legendper->SetBorderSize(0);
	legendper->AddEntry(dataerror, "Data", "lp");
	legendper->AddEntry(crone, "PYTHIA8 Monash", "l");
	legendper->AddEntry(croffe, "PYTHIA8 Monash no CR", "l");
	legendper->AddEntry(epose, "EPOS LHC", "l");
	//legendper -> AddEntry (ropee, "PYTHIA8 Rope", "l");
	legendper->Draw();

	/*
	if (option == 1)
		filenames = {
				 Form("%s_INELg0_cent_1_0.000000_1.000000.root",rname.Data())
				,Form("%s_INELg0_cent_1_0.000000_5.000000.root",rname.Data())
				,Form("%s_INELg0_cent_1_0.000000_100.000000.root",rname.Data())};
	else if (option ==4){
		filenames = {
				 Form("%s_INELg0_cent_1_0.000000_5.000000.root",rname.Data())
				,Form("%s_INELg0_cent_1_0.000000_100.000000.root",rname.Data())};
	
	}


	vector<TString> centtempl;
	vector<TString> centtemph;

	if (option ==1) {
		centtempl =  {"0","0","0"};
		centtemph =  {"1","5","100"};
	}
	else if (option == 4) {
		centtempl =  {"0","0"};
		centtemph =  {"5","100"};
	}

	nn =0;
	for (auto name : filenames){
		TFile* f = TFile::Open(Form("../FinalResults/%s",name.Data()));
		auto h = (TGraphAsymmErrors*) gROOT->FindObject("data");
		Double_t foundy = 0;
		double eh =0,el =0,dndeta=0;
		double eh10 =0,el10 =0,dndeta10=0;
		for (int k=0; k<h->GetN();k++){
			double x, y;
			h->GetPoint(k,x,y);
			if (abs(x-0.05)<0.01) foundy = y;
			if (x>-0.5 && x<0.5) {
				eh += h -> GetErrorYhigh(k)*0.1;
				el += h -> GetErrorYlow(k)*0.1;
				dndeta += y*0.1;
			}
			if (x>-1 && x<1) {
				eh10 += h -> GetErrorYhigh(k)*0.1;
				el10 += h -> GetErrorYlow(k)*0.1;
				dndeta10 += y*0.1;
			}
		}

		//cout<<"cent : "<<centtempl[nn]<<" - "<<centtemph[nn]<<" "<< " "<<dndeta<<" + "<<eh<<" - "<<el<<endl;
		printf("| %s - %s | | | %.2f | +%.2f -%.2f | | %.2f | +%.2f -%.2f | | | | \n ",centtemp[nn].Data(),centtemp[nn+1].Data(),dndeta,eh,el,dndeta10/2,eh10/2,el10/2);
		nn++;
	}
*/
}

TGraphAsymmErrors *CalcRatio(TGraphAsymmErrors *num, TGraphAsymmErrors *den)
{
	TGraphAsymmErrors *ratio = (TGraphAsymmErrors *)num->Clone();
	ratio->SetName(Form("%s_%d", num->GetName(), modelcounter++));
	ratio->SetTitle(Form("%s_%d", num->GetName(), modelcounter++));
	ratio->SetLineColor(num->GetLineColor());
	ratio->SetLineStyle(num->GetLineStyle());
	ratio->SetLineWidth(num->GetLineWidth());
	ratio->SetMarkerStyle(num->GetMarkerStyle());
	ratio->SetFillColorAlpha(num->GetFillColor(), 0.5);
	//ratio->SetFillStyle(0);
	//ratio->SetFillColor(0);
	//ratio->SetMarkerSize(0);

	for (Int_t i = 0; i < num->GetN(); i++)
	{
		double x, y;
		double xden, yden;
		num->GetPoint(i, x, y);
		den->GetPoint(i, xden, yden);
		if (yden != 0)
			ratio->SetPoint(i, x, y / yden);

		double numeh = num->GetErrorYhigh(i) / y;
		double numel = num->GetErrorYlow(i) / y;
		double nume = (numeh + numel) / 2.;

		double deneh = den->GetErrorYhigh(i) / yden;
		double denel = den->GetErrorYlow(i) / yden;
		double dene = (deneh + denel) / 2.;

		double esum = sqrt(nume * nume + dene * dene);

		ratio->SetPointEYhigh(i, esum * y / yden);
		ratio->SetPointEYlow(i, esum * y / yden);
	}
	return ratio;
}

TGraphAsymmErrors *CalcRatio(TGraphAsymmErrors *num, TH1 *den)
{
	TGraphAsymmErrors *ratio = (TGraphAsymmErrors *)num->Clone();
	ratio->SetName(Form("%s_%d", den->GetName(), modelcounter++));
	ratio->SetTitle(Form("%s_%d", den->GetName(), modelcounter++));
	ratio->SetLineColor(den->GetLineColor());
	ratio->SetLineStyle(den->GetLineStyle());
	ratio->SetLineWidth(den->GetLineWidth());
	ratio->SetMarkerStyle(den->GetMarkerStyle());
	//ratio->SetFillColorAlpha(den->GetFillColor(),0.5);
	//ratio->SetFillStyle(0);
	//ratio->SetFillColor(0);
	//ratio->SetMarkerSize(0);

	for (Int_t i = 0; i < num->GetN(); i++)
	{
		double x, y;
		double xden, yden;
		num->GetPoint(i, x, y);
		int ibin = den->GetXaxis()->FindBin(x);
		yden = den->GetBinContent(ibin);
		if (yden != 0)
			ratio->SetPoint(i, x, y / yden);

		double numeh = num->GetErrorYhigh(i) / y;
		double numel = num->GetErrorYlow(i) / y;
		double nume = (numeh + numel) / 2.;

		double dene = den->GetBinError(ibin);

		double esum = sqrt(nume * nume + dene * dene);

		ratio->SetPointEYhigh(i, esum * y / yden);
		ratio->SetPointEYlow(i, esum * y / yden);
	}
	return ratio;
}
TGraphAsymmErrors *CCCC(TGraphAsymmErrors *num, TGraphAsymmErrors *den, double temp){

	TGraphAsymmErrors *ratio = (TGraphAsymmErrors *)num->Clone();
	ratio->SetName(Form("%s_%d", num->GetName(), modelcounter++));
	ratio->SetTitle(Form("%s_%d", num->GetName(), modelcounter++));
	ratio->SetLineColor(num->GetLineColor());
	ratio->SetLineStyle(num->GetLineStyle());
	ratio->SetLineWidth(num->GetLineWidth());
	ratio->SetMarkerStyle(num->GetMarkerStyle());
	ratio->SetFillColorAlpha(num->GetFillColor(), 0.5);
	//ratio->SetFillStyle(0);
	//ratio->SetFillColor(0);
	//ratio->SetMarkerSize(0);

	for (Int_t i = 0; i < num->GetN(); i++)
	{
		double x, y;
		double xden, yden;
		num->GetPoint(i, x, y);
		yden = den->Eval(x);
		if (yden != 0)
			ratio->SetPoint(i, x, y / yden);

		double numeh = num->GetErrorYhigh(i) / y;
		double numel = num->GetErrorYlow(i) / y;


		ratio->SetPointEYhigh(i, numeh * y / yden);
		ratio->SetPointEYlow(i, numel * y / yden);
	}
	return ratio;
}

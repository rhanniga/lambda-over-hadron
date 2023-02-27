#ifndef BSHELPER_H
#define BSHELPER_H
#include <algorithm>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TCanvas.h>
#include <vector>
#include <map>
#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TGraph.h>
//---------------------------------------
//  typedef, const, global
//---------------------------------------

typedef std::vector<int>      Int1D;
typedef std::vector<Int1D>    Int2D;
typedef std::vector<Double_t> Double1D;
typedef std::vector<Double1D> Double2D;

const int min_int = std::numeric_limits<int>::min()+1;


//---------------------------------------
//  Utils
//---------------------------------------

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}
void __DEBUG_BS(int i, const char * fmt,... );

//---------------------------------------
//   Manage Root
//---------------------------------------

TFile * LoadRoot( TString filename );

//---------------------------------------
//  RANGER 
//---------------------------------------
class BSRanger {
public:
  typedef int value_type;
  struct iterator {
    iterator(size_t counter) : counter(counter) {}
    iterator operator++(){ return iterator(++counter); }
    bool operator!=(iterator o) { return counter != o.counter; }
    value_type operator*(){ return value_type{counter}; }
  private:
    value_type counter;
  };
  BSRanger(int begin, int end){ SetRange(begin,end); }
  void SetRange( int begin, int end ){if( begin>end )ErrorExit("being>end"); fBegin=begin;fEnd=end;}
  iterator begin(){ return iterator(fBegin); }
  iterator end()  { return iterator(fEnd+1);   }
private:
  int fBegin;
  int fEnd;
};

template< class T >
BSRanger range( T & t, int begin=0, int end=-1 ){ return BSRanger( begin, end<0?t.size()-1:end  ); }
BSRanger range( int end ){ return BSRanger( 0, end-1 ); }
BSRanger range( int begin, int end ){ return BSRanger( begin, end ); }

BSRanger bin_range( int begin, int end ){
  if( begin < 1 ) ErrorExit("begin for bin_range must larger than 0 ");
  return range( begin, end );
}
BSRanger bin_range( int end ){ return bin_range( 1, end ); }
template<class T>
BSRanger bin_range( const std::vector<T> &t, int begin=1, int end=min_int ){ return bin_range( begin, (end<0?std::min(int(t.size()-1),abs(end)):end));  }

/*
//---------------------------------------
//  RANGER - pair
//---------------------------------------
template< class Container, class Item >
class BSRangerPair {
public:
  struct item { size_t I; Item & V; };
  typedef item value_type;
  struct iterator {
    iterator(size_t _counter, Container& _con) : counter(_counter),fCont(_con) {}
    iterator operator++(){return iterator(++counter,fCont); }
    bool operator!=(iterator o) { return counter != o.counter; }
    value_type operator*(){ return value_type{counter,fCont[counter]}; }
  private:
    size_t counter;
    Container & fCont;
  };
  BSRangerPair( Container & cont,int begin,int end ):fCont(cont){ SetRange(begin, end ); };
  void SetRange( int begin, int end ){ if( begin>end ){ exit(1); }fBegin=begin;fEnd=end;}
  iterator begin(){ return iterator(fBegin, fCont); }
  iterator end()  { return iterator(fEnd, fCont);   }
private:
  int fBegin;
  int fEnd;
  Container & fCont;
};

template< class T >
BSRangerPair<vector<T>,T> range_pair( vector<T> & t, int begin=0, int end=-1 ){
  if( end >= int(t.size()) ) ErrorExit("\"end\" must smaller than size of array");
  return BSRangerPair<vector<T>,T>( t, begin, end<0?t.size():end );
}
*/

//====================================
// HistManager
//====================================

class BSHistHelper {
  // TODO
};

class BSTHnSparseHelper {
public:
  BSTHnSparseHelper(THnSparse*h):fH(h){LoadBinFromHist();}
  BSTHnSparseHelper(TObject*o);
  THnSparse * operator->(){ return fH; }
  THnSparse * Data(){ return fH; }
  THnSparse * Val(){ return fH; }
  THnSparse * Hist(){ return fH; }

  //================
  //  AXIS
  //================
  int  GetAxisID(TString name );
  TAxis *GetAxis(int i);
  TAxis *GetAxis(TString name){ return GetAxis(GetAxisID(name)); } 
  int  GetNdim(){ return fH->GetNdimensions(); };
  Int_t GetNbins(int i){ return GetAxis(i)->GetNbins();}
  Int_t GetNbins(TString n){ return GetNbins(GetAxisID(n));}
  Int_t GetNbinsU(int i){ return GetBin(i).size()-1;}
  Int_t GetNbinsU(TString n){ return GetNbinsU(GetAxisID(n));}
  void PrintAxis(Option_t * opt="");
  TH1D * GetTH1( TString name, Int_t xDim, Int1D bin, Option_t*opt="");
  TH1D * GetTH1( TString name, Int_t xDim, Int2D bin, Option_t*opt="");
  TH1D * GetTH1( TString name, Int1D bin, Option_t*opt=""){ return GetTH1( name, GetNdim()-1, bin, opt ); }
  TH2D * GetTH2( TString name, Int_t xDim,Int_t yDim, Int1D bin, Option_t*opt="");
  TH2D * GetTH2( TString name, Int_t xDim,Int_t yDim, Int2D bin, Option_t*opt="");

  BSRanger BinRange( int i ){ return bin_range( GetBin(i) ); }
  BSRanger BinRange( TString n ){ return bin_range( GetBin(n) ); }
  void ClearBin(){ LoadBinFromHist(); }
  void ClearBin(int i){ LoadBinFromHist(i); }
  const Double1D& GetBin( int i){ return fCBins[i] ; }
  const Double1D& GetBin( TString n ){ return GetBin(GetAxisID(n)); }
  void SetBin( Double2D bins );
  void SetBin( int i, Double1D bins );
  void SetBin( TString name, Double1D bins ){ SetBin( GetAxisID(name), bins); }
  static BSTHnSparseHelper Load( TString name, TObject * clist=nullptr );
  static BSTHnSparseHelper Load( TString name, TString fname, TString lname="" );
private:
  void LoadBinFromHist();
  void LoadBinFromHist( int iaxis );
  THnSparse * fH;
  Double2D    fCBins;
  Double2D    fCBinsB;
};

//__________________________________________________________
BSTHnSparseHelper::BSTHnSparseHelper(TObject*o){
  if(!o) ErrorExit("No Object");
  fH = dynamic_cast<THnSparse*>(o);
  if(!fH) ErrorExit(Form("%s is not THnSparse but %s", o->GetName(),o->ClassName()) );
  LoadBinFromHist();
}

//__________________________________________________________
TAxis * BSTHnSparseHelper::GetAxis(int i){
  if( i < 0 || i>=GetNdim() )
    ErrorExit(Form("Wrong Axis Index %d of %s", GetNdim(), fH->GetName() ));
  return fH->GetAxis(i); 
}
//__________________________________________________________
int BSTHnSparseHelper::GetAxisID( TString name ){
  for( int i=0;i<fH->GetNdimensions();i++ )
    if( name == fH->GetAxis(i)->GetName() )
      return i;
  fH->Print();PrintAxis();
  ErrorExit("No Axis "+name);
  return -1;
}

//__________________________________________________________
void BSTHnSparseHelper::PrintAxis(Option_t * opt){
  TString opts = opt;
  for( int i=0;i<fH->GetNdimensions();i++ ){
    TAxis * ax = fH->GetAxis(i);
    cout<<Form("%2d: %10s",i, fH->GetAxis(i)->GetName() );
    if( opts == "all" ){
      cout<<Form(" (nbin/min/max) = %4d %6.2f %8.2f ::",ax->GetNbins(),ax->GetXmin(), ax->GetXmax());
      for( int j=1;j<=ax->GetNbins()&&j<10;j++ ){
        const char * label = ax->GetBinLabel(j);
        if( !TString(label).IsNull() )
          cout<<Form(" %s", label );
        else
          cout<<Form(" %.2f", ax->GetBinLowEdge(j) );
      }
    }
    cout<<endl;
  }
}

//__________________________________________________________
TH1D * BSTHnSparseHelper::GetTH1( TString name, Int_t xDim, vector<int> bin, Option_t*opt){
  Int2D newbin;
  for( auto b : bin ) newbin.push_back({b,b});
  return GetTH1( name, xDim, newbin, opt );
}
//__________________________________________________________
TH1D * BSTHnSparseHelper::GetTH1( TString name, Int_t xDim, Int2D bin, Option_t*opt){
  //== Histogram Name and Title
  if( name.EndsWith("-") ) name+=Form("%sP%02d",fH->GetName(),xDim);
  TString title = fH->GetTitle();
  for( UInt_t i=0;i<bin.size();i++ ){
    Int1D b = bin[i];
    Int1D nb = { 0, 0 };
    int maxnbin = fCBinsB[i].size()-2;
    auto ax = GetAxis(i);
    if( b[0]<-1 || b[1]<-1 || b[0] > maxnbin || b[1] > maxnbin || b[0] > b[1] ) ErrorExit(Form("Wrong bin : %i : %s : %d %d",i,ax->GetName(),b[0],b[1])); 
    if( int(i) == xDim || b.size()==0 || ( b[0]<=0 && b[1]<=0 ) ){
    }else{
      if( b[0] > 0 ) nb[0] = fCBinsB[i][b[0]];
      if( b[1] > 0 ) nb[1] = fCBinsB[i][b[1]+1]-1;
      //if( b0 > b1 ) b1 = b1;
      name+=Form("%s%03d%03d",GetAxis(i)->GetName(), b[0],b[1]);
      TString label[2];
      label[0] = ax->GetBinLabel(nb[0]);
      label[1] = ax->GetBinLabel(nb[1]);
      if( label[0]=="" ) label[0] = Form("%.6f",ax->GetBinLowEdge( nb[0] )); else label[0] = Form("%d",nb[0]);// TODO Temp
      if( label[1]=="" ) label[1] = Form("%.6f",ax->GetBinUpEdge( nb[1] )); else label[1] = Form("%d",nb[1]); // TODO Temp
      TPMERegexp("\\.?0+$").Substitute(label[0],"");
      TPMERegexp("\\.?0+$").Substitute(label[1],"");
      title+=Form(" %s:%s", ax->GetName(), label[0].Data());
      if( label[0]!=label[1] ) title+=Form("-%s",label[1].Data());
    }
    fH->GetAxis(i)->SetRange(nb[0], nb[1]);
  }
  auto h = fH->Projection( xDim, opt );
  h->SetNameTitle(name,title);
  return h;
}
//__________________________________________________________
TH2D * BSTHnSparseHelper::GetTH2( TString name, Int_t xDim, Int_t yDim, vector<int> bin, Option_t*opt){
  Int2D newbin;
  for( auto b : bin ) newbin.push_back({b,b});
  return GetTH2( name, xDim,yDim, newbin, opt );
}
TH2D * BSTHnSparseHelper::GetTH2( TString name, Int_t xDim,Int_t yDim, Int2D bin, Option_t*opt){
  //== Histogram Name and Title
  if( name.EndsWith("-") ) name+=Form("%sP%02d",fH->GetName(),xDim);
  TString title = fH->GetTitle();
  for( UInt_t i=0;i<bin.size();i++ ){
    Int1D b = bin[i];
    Int1D nb = { 0, 0 };
    int maxnbin = fCBinsB[i].size()-2;
    auto ax = GetAxis(i);
    if( b[0]<-1 || b[1]<-1 || b[0] > maxnbin || b[1] > maxnbin || b[0] > b[1] ) ErrorExit(Form("Wrong bin : %i : %s : %d %d",i,ax->GetName(),b[0],b[1])); 
    if( int(i) == xDim || int(i) == yDim || b.size()==0 || ( b[0]<=0 && b[1]<=0 ) ){
    }else{
      if( b[0] > 0 ) nb[0] = fCBinsB[i][b[0]];
      if( b[1] > 0 ) nb[1] = fCBinsB[i][b[1]+1]-1;
      //if( b0 > b1 ) b1 = b1;
      name+=Form("%s%03d%03d",GetAxis(i)->GetName(), b[0],b[1]);
      TString label[2];
      label[0] = ax->GetBinLabel(nb[0]);
      label[1] = ax->GetBinLabel(nb[1]);
      if( label[0]=="" ) label[0] = Form("%.6f",ax->GetBinLowEdge( nb[0] )); else label[0] = Form("%d",nb[0]);// TODO Temp
      if( label[1]=="" ) label[1] = Form("%.6f",ax->GetBinUpEdge( nb[1] )); else label[1] = Form("%d",nb[1]); // TODO Temp
      TPMERegexp("\\.?0+$").Substitute(label[0],"");
      TPMERegexp("\\.?0+$").Substitute(label[1],"");
      title+=Form(" %s:%s", ax->GetName(), label[0].Data());
      if( label[0]!=label[1] ) title+=Form("-%s",label[1].Data());
    }
    fH->GetAxis(i)->SetRange(nb[0], nb[1]);
  }
  auto h = fH->Projection( yDim,xDim, opt );
  h->SetNameTitle(name,title);
  return h;
}



//__________________________________________________________
void BSTHnSparseHelper::SetBin( int iaxis, Double1D bins ){
  auto ax = GetAxis(iaxis);
  if( bins.size()==0 ) LoadBinFromHist ( iaxis );
  else{
    fCBins[iaxis] = bins;
    fCBinsB[iaxis] = { 0 };
    for( UInt_t ib=0;ib<bins.size();ib++ ){
      int bin = ax->FindBin( bins[ib] );
      //TODO if( ib>0 && bin <= fCBinsB[i][ib-1] ) ErrorExit("Wrong bin");
      fCBinsB[iaxis].push_back( ax->FindBin( bins[ib] ) );
    }
  }
}
//__________________________________________________________
void BSTHnSparseHelper::SetBin( Double2D bins ){
  for( UInt_t i=0;i<bins.size();i++ )
    SetBin( i, bins[i] );
}

//__________________________________________________________
BSTHnSparseHelper BSTHnSparseHelper::Load( TString name, TObject * clist){
  if( !clist ) clist = gDirectory;
  auto o = clist->FindObject(name);
  if( !o ) ErrorExit( "No THnSparse "+name );
  return BSTHnSparseHelper(o);
}
//__________________________________________________________
BSTHnSparseHelper BSTHnSparseHelper::Load( TString name, TString fname, TString lname){
  TFile* f = gROOT->GetFile(fname);
  if(!f) f = TFile::Open(fname);
  if(!f) ErrorExit( "No File : "+fname );
  TObject * cl = f;
  if( ! lname.IsNull() )
    cl = f->Get(lname);
  return Load(name, cl);
}

//__________________________________________________________
void BSTHnSparseHelper::LoadBinFromHist( int iaxis){
  auto ax = GetAxis(iaxis);
  fCBinsB[iaxis] = {0};
  for( int ib=1;ib<=ax->GetNbins()+1;ib++ ){
    fCBins[iaxis].push_back( ax->GetBinLowEdge(ib) );
    fCBinsB[iaxis].push_back( ib );
  }
}
//__________________________________________________________
void BSTHnSparseHelper::LoadBinFromHist(){
  fCBins.resize( GetNdim() );
  fCBinsB.resize( GetNdim() );
  for( int i=0;i<GetNdim();i++ )
    LoadBinFromHist(i);
}

//==========================
//  DRAWING
//==========================
class BSPlotManager {
public:
  BSPlotManager(TString name):fName(name){;}
  BSPlotManager* AddPadGroup(TString name);
  TString GetTitle(){ return fTitle; }
  TString GetTitle(int i){ return fTitles[i]; }
  TString GetCurrentTitle(){ return fTitles[fPadCursor]; }
  void    SetTitle(TString title){ fTitle=title; }
  TPad * AddPadInGrid(TString name, TString title,Option_t *opt, int nx, int ny, double wx=-1, double wy=-1){
    if(name.EndsWith("-")) name += Form("%02d%02d",nx,ny);
    return AddPad( name, title, Form("%sgrid;",opt),-1, -1, wx, wy, nx, ny );
  }
  TPad * AddPadInGridNext(TString name, TString title,Option_t *opt="", double wx=-1, double wy=-1){
    fGridCursorX++;
    if( fGridCursorX >= fGridXn ){
      fGridCursorX=0;
      fGridCursorY++;
    }
    return AddPadInGrid( name, title,opt, fGridCursorX,fGridCursorY, wx, wy );
  }
  TH2F* SetBaseHist(TString title, TString option="", TString xtitle="", TString ytitle="", double x0=0,double x1=0,double y0=0,double y1=0){
    auto hb = GetOrCreateBaseHist(title,option);
    hb->GetXaxis()->SetTitle(xtitle.Data());
    hb->GetYaxis()->SetTitle(ytitle.Data());
    if( x1-x0 > 1e-9 ) hb->GetXaxis()->Set(10, x0, x1 );
    if( y1-y0 > 1e-9 ) hb->GetYaxis()->Set(10, y0, y1 );
    hb->Draw();
    gPad->Update();
    return hb;
  }
  TH2F* GetOrCreateBaseHist(TString title, TString option=""){
    auto hb = fBaseHist[fPadCursor];
    if( !hb ) {
      if( title.IsNull() ) title = GetCurrentTitle();
      hb = fBaseHist[fPadCursor] = new TH2F(Form("hb%s_%d",fName.Data(),fPadCursor),title,10,0,10,10,0,10);
      hb->SetStats(kFALSE);
      hb->Draw();
    }
    return hb;
  }
  TH2F* GetOrCreateBaseHist(TString title, TString option, double x0, double x1, double y0, double y1){
    auto hb = SetRange( x0, x1, y0, y1 );
    hb->SetTitle(title);
    gPad->Update();
    return hb;
  }
  TH2F* GetBaseHist(){ return fBaseHist[fPadCursor]; }
  TH2F* SetRange( double x0, double x1, double y0,double y1){
    auto hb = GetOrCreateBaseHist("","");
    hb->GetXaxis()->Set(10, x0, x1 );
    hb->GetYaxis()->Set(10, y0, y1 );
    gPad->Update();
    return hb;
  }

  TH2F* ReRange(TString optymax,double ymax,TString optymin="",double ymin=1){
    //FIXME Graph?
    double _ymax,_ymin;
    int i = 0;
    for( auto hist : fHists[fPadCursor] ){
      auto cmax = ((TH1*)hist)->GetMaximum();
      auto cmin = ((TH1*)hist)->GetMinimum();
      if( i ==0 ){
        _ymax = cmax;_ymin=cmin;
      }else{
        if( _ymax < cmax ) _ymax=cmax;
        if( _ymin > cmin ) _ymin=cmin;
      }
      i++;
    }
    auto hb = GetOrCreateBaseHist("");
    if( optymax=="*" ) _ymax*=ymax;
    else if( optymax=="+" ) _ymax+=ymax;
    else _ymax = hb->GetYaxis()->GetXmax();
    if( optymin=="*" ) _ymin*=ymin;
    else if( optymin=="+" ) _ymin+=ymin;
    else _ymin = hb->GetYaxis()->GetXmin();
    hb->GetYaxis()->Set(10, _ymin, _ymax );
    gPad->Update();
    return hb;
  }

  array<double,4> GetObjectRange(TObject*o){
    if( o->InheritsFrom(TH1::Class()) ){
      auto h = dynamic_cast<TH1*>(o);
      return { h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),h->GetMinimum(),h->GetMaximum() };
    }else if(  o->InheritsFrom(TGraph::Class())){
      auto g = dynamic_cast<TGraph*>(o);
      return { g->GetXaxis()->GetXmin(), g->GetXaxis()->GetXmax(),g->GetYaxis()->GetXmin(),g->GetYaxis()->GetXmax() };
    }
    return { 0,10,0,10};

  }
  TH2F * SetRange( array<double,4> mm ){
    return SetRange( mm[0],mm[1],mm[2],mm[3] );
  }

  TH2F *SetRange( TObject *o ){
    return SetRange( GetObjectRange(o) );
  }
  TH2F *SetMain( TObject *h, Option_t* opt="", TString title="" ){
    auto hb = GetOrCreateBaseHist(title);
    TString opts = opt;
    if( opts.Contains("title") ) {
      if( title.IsNull() )
        hb->SetTitle( h->GetTitle() );
      else
        hb->SetTitle( title );
    }
    if( opts.Contains("stat") && h->InheritsFrom(TH1::Class())){
      auto hh = dynamic_cast<TH1*>(h);
      double  stat[10]; 
      hh->GetStats(stat);
      hb->PutStats(stat);
      hb->SetStats(kTRUE);
    }
    if( opts.Contains("range")) SetRange(h);
    gPad->Update();
    return hb;
  }
  TLegend* SetLegend(double x1, double y1, double x2, double y2, TString header="", Option_t*opt="brNDC"){
    auto l = fLegends[fPadCursor];
    l->SetHeader( header.Data());
    l->SetX1NDC(x1);l->SetX2NDC(x2);l->SetY1NDC(y1);l->SetY2NDC(y2);
    l->SetOption(opt);
    l->SetEntryOption("C");
    l->SetFillStyle(0); l->SetBorderSize(0); l->SetTextSize(0.04);
    return l;
  }
  TLegend* AddLegendEntry( TObject *h, TString title="", Option_t* opt="p"){
    auto l = fLegends[fPadCursor];
    if( title.IsNull() ) title=h->GetTitle();
    l->AddEntry(h,title,opt);
    return l;
  }
  TLegend* DrawLegend(Option_t *opt="brNDC"){
    auto l = fLegends[fPadCursor];
    l->Draw(opt);
    return l;
  }

  TH1* DrawHist( TH1 *h, Option_t* opt="", Option_t *gopt="",int color=0, int marker=20, TString legend="",Option_t* legendopt="p",TString xtitle="",TString ytitle=""){
    return dynamic_cast<TH1*>( DrawHist((TObject*)h,opt,gopt,color,marker,legend,legendopt,xtitle,ytitle));
  }
  TObject* DrawHist( TObject *h, Option_t* opt="", Option_t *gopt="",int color=0, int marker=20, TString legend="",Option_t* legendopt="p",TString xtitle="x",TString ytitle="y"){
    auto hb = GetBaseHist();
    auto hists = fHists[fPadCursor];
    TString gopts = gopt;gopts+=";";
    if( !hb ) hb=SetRange(h);
    if( !xtitle.IsNull() ){
      hb->GetXaxis()->SetTitle(xtitle.Data());
      hb->GetXaxis()->SetTitleOffset(0.9);
      hb->GetXaxis()->SetTitleSize(0.06);
      hb->GetXaxis()->SetLabelSize(0.05);
    }
    if( !ytitle.IsNull() ){
      hb->GetYaxis()->SetTitle(ytitle.Data());
      hb->GetYaxis()->SetTitleOffset(0.9);
      hb->GetYaxis()->SetTitleSize(0.06);
      hb->GetYaxis()->SetLabelSize(0.05);
    }
    if( !gopts.IsNull() ) SetMain(h,gopt);
    if(legend.IsNull()) legend = h->GetTitle();
    if(gopts.Contains("legend")) {
      AddLegendEntry(h,legend,legendopt);
      DrawLegend();
    }
    if(gopts.Contains("clone"))
      h=(TH1*)h->DrawClone(TString(opt)+"same");
    else
      h->Draw(TString(opt)+"same");
    if(gopts.Contains("*color;"))
      color = fColors[hists.size()];
    if(gopts.Contains("*marker;"))
      marker = fMarkers[hists.size()];
    if( h->InheritsFrom(TAttLine::Class()) ){
      auto hh = dynamic_cast<TAttLine*>(h);
      hh->SetLineColor(color);
    }
    if( h->InheritsFrom(TAttMarker::Class()) ){
      auto hh = dynamic_cast<TAttMarker*>(h);
      hh->SetMarkerColor(color);
      hh->SetMarkerStyle(marker);
    }
    gPad->Update();
    fHists[fPadCursor].push_back(h);
    return h;
  }


  TPad * AddPad( TString name, TString title, Option_t * opt, double x, double y, double wx, double wy, int nx=-1, int ny=-1 ){
    double xR, yR, wxR, wyR;
    TString opts = opt;
    //=== GRID
    if( opts.Contains("grid;") ){
      x = (fPadXw/fGridXn*nx)+fPadX;
      y = (fPadYw/fGridYn*ny)+fPadY;
      if( wx < 0 ) wx = fGridXw;
      if( wy < 0 ) wy = fGridYw;
      fGridCursorX = nx;
      fGridCursorY = ny;
    }
    //===  Real Size
    if( x > 1 ) xR = x;
    else xR = x*fScreenXw+fScreenX;
    if( y > 1 ) yR = y;
    else yR = y*fScreenYw+fScreenY;
    if( wx > 1 ) wxR = wx;
    else wxR = wx*fScreenXw;
    if( wy > 1 ) wyR = wy;
    else wyR = wy*fScreenYw;

    //=== Square Ratio : wy will be ignored
    if( opts.Contains("square;")){
      wyR=wxR;
    }
    auto pad = new TCanvas( name, "", xR, yR, wxR, wyR );
    fCanvas.push_back( pad );
    fBaseHist.push_back(nullptr);
    fLegends.push_back( new TLegend );
    fHists.push_back(vector<TObject*>());
    fPadCursor = fCanvas.size()-1;
    fTitles.push_back(title);
    pad->SetRightMargin(0.03);
    pad->SetTopMargin(0.02);
    pad->SetLeftMargin(0.12);
    pad->SetBottomMargin(0.12);
    return pad;
  }

  void Draw(){
    for( UInt_t i=0;i<fPads.size();i++ ){
      auto p = fPads[i];
      cout<<Form("%f %f %f %f", p->GetX1(), p->GetY1(), p->GetX2(), p->GetY2())<<endl;
      return;
      auto c = new TCanvas(Form("c%s",p->GetName()), "", p->GetX1()*fScreenXw+fScreenX, p->GetY1()*fScreenYw+fScreenY, (p->GetX2()-p->GetX1())*fScreenXw, (p->GetY2()-p->GetY1())*fScreenYw);
      c->Draw();
    }

  }

  TCanvas * DrawPad(TString name);
  TCanvas * GetPad(int i);
  TCanvas * GetCnavas(TString name);
  void SetGrid(int nx, int ny, double xw=0, double yw=0 ){
    fGridXn = nx;
    fGridYn = ny;
    fGridXw = (xw<1e-4)?1./nx:xw;
    fGridYw = (yw<1e-4)?1./ny:yw;
  }

  void SetColors(Int1D colors){ fColors=colors; };
  void SetMarker(Int1D markers){ fMarkers=markers; }

private:
  TString          fName;
  TString          fTitle;
  vector<TCanvas*> fCanvas;
  vector<TPad*>    fPads;
  vector<TH2F*>    fBaseHist;
  vector<TLegend*>    fLegends;
  vector<vector<TObject*>> fHists;
  vector<TString>  fTitles;
  //vector<vector<TPad*> fGrid;
  vector<BSPlotManager*> fSubManager;
  int    fGridXn=4;
  int    fGridYn=3;
  double fGridXw=1./fGridXn;
  double fGridYw=1./fGridYn;
  double fScreenX=50;
  double fScreenY=50;
  double fScreenXw=1400;
  double fScreenYw=900;
  double fPadX=0;
  double fPadY=0;
  double fPadXw=1;
  double fPadYw=1;
  int    fGridCursorX=-1;
  int    fGridCursorY=0;
  int    fPadCursor=-1;

  // Style Template;
  Int1D fColors={kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kBlack,kRed,kBlue,kGreen+3,kMagenta+2};
  Int1D fMarkers={24,25,26,27,28,30,20,21,22,33,34,24,25,26,27,28,30,20,21,22,33,34,24,25,26,27,28,30,20,21,22,33,34,24,25,26,27,28,30,20,21,22,33,34};
};

//====================================
//  Utils
//====================================
void __DEBUG_BS(int i, const char * fmt,... ){
  char temp[1024];
  va_list marker;
  va_start(marker, fmt);
  vsprintf(temp, fmt, marker);
  va_end(marker);

  std::cout<< Form("DEBUG %d : %s", i, temp)<<std::endl;
}

//====================================
//   Manage Root
//====================================

//__________________________________________________________
TFile * LoadRoot( TString filename ){
  auto f = gROOT->GetFile(filename);
  if(!f) f = TFile::Open(filename);
  if( !f ) ErrorExit("No File "+filename);
  return f;
}

std::array<double,4> GetMeanRMS(TH1* h){
  double r = h->GetRMS();
  double re = h->GetRMSError();
  double m = h->GetMean();
  double me = h->GetMeanError();
  double RMSReal = TMath::Sqrt(r*r+m*m);
  double RMSRealError =  TMath::Sqrt( r*r*re*re + m*m*me*me )/RMSReal;
  return { m, me, RMSReal, RMSRealError };
}
#endif

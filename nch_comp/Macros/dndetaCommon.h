#include "BSHelper.cxx"
//
TFile * LoaddndetaResults( TString name, TString runnum){
  if( !runnum.IsNull() ) runnum = "_"+runnum;
  auto fname = "../Results/"+name+"/AnalysisResults_"+name+runnum+".root";
  return  LoadRoot( fname ) ;
}
//__________________________________________________________
TObject*  LoaddndetaResultList( TFile *fh, TString clistname){
  auto l = fh->Get( clistname );
  if(!l) ErrorExit("No list "+clistname);
  return l;
}

//__________________________________________________________
TObject*  LoaddndetaResultList( TString fname, TString clistname, TString runnum=""){
  auto f = LoaddndetaResults( fname, runnum );
  return LoaddndetaResultList( f, clistname );
}

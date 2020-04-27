#include "chisquare_distr.C"

void test_chisquare_distr(int ndf = 8, int nevts = 1e5) { 
  
  int chi2max = ndf*3;
  //gROOT->LoadMacro("chisquare_distr.C");
  TF1 *fc = new TF1("fc",chisquare_distr,0,chi2max,2);
  fc->SetParameters(ndf,1);

  fc->Draw();

  TH1D *h = new TH1D("h","",int(chi2max*10),0,chi2max);
  h->FillRandom("fc",nevts);

  h->Draw();

  h->Fit("fc","R","",0.1,chi2max);
  gStyle->SetOptFit(1111);
}

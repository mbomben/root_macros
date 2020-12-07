#include "RooBinomial.h"
#include "RooRealVar.h"
#include "TROOT.h"
#include "RooPlot.h"



void BinomialDraw (int kval = 1, int Nval = 3, double pval = 0.5) {

  gROOT->ProcessLineSync(".x RooBinomial.cxx+");

  RooRealVar N("N","Trails",Nval); // number of trials

  RooRealVar p("p","Success prob.",pval,0.,1.); // probability

  RooRealVar k("k","Success",kval,0,Nval+1); // number of successes 

  k.setBins(Nval+1);

  RooBinomial pdf("pdf","Binomial PDF",k,N,p);

  // let's normalize it
  RooAbsReal *ipdfk = pdf.createIntegral(k);

  RooPlot* kframe = k.frame();

  pdf.plotOn(kframe);

  kframe->Draw();

}


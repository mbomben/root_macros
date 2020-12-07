#include "RooBinomial.h"
#include "RooRealVar.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "RooDataSet.h"



void BinomialFit  (int kval = 1, int Nval = 3, double pval = 0.5, int nevents = 1000)  {
  

  gROOT->ProcessLineSync(".x RooBinomial.cxx+");

  RooRealVar N("N","Trails",Nval); // number of trials

  RooRealVar p("p","Success prob.",pval,0.,1.); // probability

  RooRealVar k("k","Success",kval,0,Nval+1); // number of successes 

  k.setBins(Nval+1);

  RooBinomial pdf("pdf","Binomial PDF",k,N,p);

  RooPlot* kframe = k.frame();

  RooDataSet* dataset = pdf.generate(k,nevents);

  dataset->plotOn(kframe);

  pdf.fitTo(*dataset);

  pdf.plotOn(kframe);

  kframe->Draw();

}


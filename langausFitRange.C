#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TROOT.h"
#include <iostream>
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
using namespace RooFit;

void langausFitRange( const char* fileName, int HitDiscCnfg = 2, double xmax=-1 ) {
  gSystem->Load("libRooFit.so");

  TCanvas *c1 = new TCanvas("c1","",1024,768);
  c1->cd();


  TFile *file = new TFile(fileName,"READ");
  file->cd();
  gDirectory->cd("ClusteringAPIX/detector_20");
  TH1D* clusterSignal_d20 = (TH1D*)  gDirectory->Get("clusterSignal_d20");
  //clusterSignal_d20->Draw("e");

  // determining tot range
  int nbins = clusterSignal_d20->GetNbinsX();

  double xmin = clusterSignal_d20->GetBinLowEdge(1);
  if ( xmax < 0 ) { 
    double xmax = clusterSignal_d20->GetBinLowEdge(nbins+1);
  }
  double mean = clusterSignal_d20->GetMean();

  xmin += HitDiscCnfg;
  std::cout << "mean = " << mean << "\n";
  std::cout << "xmin = " << xmin << "\n";
  std::cout << "xmax = " << xmax << "\n";

  // the tot variable
  RooRealVar tot("tot","#tot [ns]",mean,xmin,xmax);
  // the tot ranges
  tot.setRange("no16Low",2.5,15.5);
  tot.setRange("no16High",16.5,30.5);
  char fit_range[256];
  sprintf(fit_range,"%s,%s","no16Low","no16High");
  std::cout << "fit_range = " << fit_range << "\n";
  // the data set
  RooDataHist* totEDH = new RooDataHist("totEDH","totEDH",tot,RooFit::Import(*clusterSignal_d20));
  // reduced data set
  //RooDataSet* data_reduce = (RooDataSet*)totEDH->reduce(RooFit::CutRange("no16Low"));
  //data_reduce->append(*((RooDataSet*)(totEDH->reduce(RooFit::CutRange("no16High")))));

  // Construct landau(t,landau_mean,landau_width) ;
  RooRealVar landau_mean("landau_mean","mean landau",5.,xmin,xmax) ;
  RooRealVar landau_width("landau_width","sigma landau",1,0.1,20) ;
  RooLandau landau("lx","lx",tot,landau_mean,landau_width) ;

  // Construct gauss(t,mg,gaussian_sigma)
  RooRealVar mg("mg","mg",0) ;
  RooRealVar gaussian_sigma("gaussian_sigma","gaussian_sigma",0.5,0.01,2) ;
  RooGaussian gauss("gauss","gauss",tot,mg,gaussian_sigma) ;


  // C o n s t r u c t   c o n v o l u t i o n   p d f 
  // ---------------------------------------

  // Set #bins to be used for FFT sampling to 10000
  tot.setBins(10000,"cache") ;


  // Construct landau (x) gauss
  RooFFTConvPdf lxg("lxg","landau (X) gauss",tot,landau,gauss) ;

  // Fit p.d.f. to data
  //RooFitResult *fitr = lxg.fitTo(*data_reduce,Range(fit_range)) ;
  RooFitResult *fitr = lxg.fitTo(*totEDH,RooFit::Range(fit_range)) ;
  // Plot data, landau pdf, landau (X) gauss pdf
  RooPlot* frame = tot.frame(RooFit::Title("landau (x) gauss convolution"),RooFit::Bins(nbins)) ;
  totEDH->plotOn(frame) ;
  //data_reduce->plotOn(frame) ;
  lxg.plotOn(frame) ;
  lxg.paramOn(frame);
  // Draw frame on canvas
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;

  // chi2
//  RooChi2Var chi2_lowstat("chi2_lowstat","chi2",lxg,*totEDH) ;
//  cout << chi2_lowstat.getVal() << endl ;
  
  TF1 * f = lxg.asTF( RooArgList(tot) );
  double mpv = f->GetMaximumX();
  std::cout << "MPV = " << mpv << "\n";

}

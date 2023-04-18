#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
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
#include "RooFitResult.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "TPaveText.h"
#include "fstream"
#include "TLegend.h"
void charge_LangausFit(const char* fileName, const char* histoName, double qmin, double qmax, int rebin=-1 )
{
  using namespace RooFit;

  gSystem->Load("libRooFit.so");

  TCanvas *c1 = new TCanvas("c1","",1024,768);
  c1->cd();

  TFile *file = new TFile(fileName,"READ");
  if ( file->IsZombie() ) { 
    std::cout << fileName << " file was not found/open.\nExiting...\n";
    exit( -1 );
  }

  TH1F* clusterSignal = (TH1F*)  gDirectory->Get(histoName);
  if ( clusterSignal == 0x0 ) {
    std::cout << histoName << " histo was not found/open.\nExiting...\n";
    exit( -2 );
  }

  if ( rebin > 0 ) {
    clusterSignal->Rebin(rebin);
  }

  //clusterSignal->Draw("e");
  int nentries = clusterSignal->GetEntries();

clusterSignal->Scale(1.0/nentries);
clusterSignal->SetStats(111);  
// determining q range
  int nbins = clusterSignal->GetNbinsX();
  //double xmax = qmax*1.1;
  double xmax = (qmax+clusterSignal->GetBinLowEdge(nbins+1))/2.0;
  double xmin = clusterSignal->GetBinLowEdge(1);
  if ( xmax < 0 ) { 
    xmax = clusterSignal->GetBinLowEdge(nbins+1);
  }
  double mean = clusterSignal->GetMean();

    //  xmin += HitDiscCnfg;


  // the cluster charge variable
  RooRealVar q("q","Cluster Charge [ke-]",mean,xmin,xmax);
  q.setRange("Range",qmin,qmax);
  char fit_range[256];
     sprintf(fit_range,"%s","Range");
  // the data set
  RooDataHist* qEDH = new RooDataHist("qEDH","qEDH",q,RooFit::Import(*clusterSignal));
  // reduced data set
  //RooDataSet* data_reduce = (RooDataSet*)qEDH->reduce(RooFit::CutRange("no16Low"));
  //data_reduce->append(*((RooDataSet*)(qEDH->reduce(RooFit::CutRange("no16High")))));

  // Construct landau(t,landau_mpv,landau_width) ;
  RooRealVar landau_mpv("landau_mpv","landau mpv",mean*0.8,qmin,qmax) ;
  RooRealVar landau_width("landau_width","landau width",1,0.01,15) ;
  //RooRealVar landau_width("landau_width","landau width",3) ;
  RooLandau landau("lx","lx",q,landau_mpv,landau_width) ;

  // Construct gauss(t,mg,gaussian_sigma)
  RooRealVar mg("mg","mg",0) ;
  RooRealVar gaussian_sigma("gaussian_sigma","gaussian_sigma",2,1,9) ;
  RooGaussian gauss("gauss","gauss",q,mg,gaussian_sigma) ;


  // C o n s t r u c t   c o n v o l u t i o n   p d f 
  // ---------------------------------------

  // Set #bins to be used for FFT sampling to 10000
  q.setBins(10000,"cache") ;


  // Construct landau (x) gauss
  RooFFTConvPdf lxg("lxg","landau (X) gauss",q,landau,gauss) ;

  // Fit p.d.f. to data
  RooFitResult *fitr = lxg.fitTo(*qEDH,RooFit::Range(fit_range),Save(),SumW2Error(kTRUE)) ;
  // Plot data, landau pdf, landau (X) gauss pdf
  RooPlot* frame = q.frame(RooFit::Title(" "),RooFit::Bins(nbins)) ;
  qEDH->plotOn(frame) ;
  //data_reduce->plotOn(frame) ;
  lxg.plotOn(frame) ;
  lxg.paramOn(frame,Layout(0.4, 0.7, 0.8),Format("NEU"));
  // Draw frame on canvas
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;
  gPad->SetTicks(1,1);
  // chi2
  //RooChi2Var chi2_lowstat("chi2_lowstat","chi2",lxg,*qEDH, RooFit::DataError(RooAbsData::RooAbsData::SumW2)) ;
  //cout << "\n\n\n\n\n\n\nChi2 = " << chi2_lowstat.getVal() << "\n\n\n\n\n\n\n" ;



  TString baseOutputFileName=fileName;

  std::cout << "baseOutputFileName = " << baseOutputFileName.Data() << "\n";

  baseOutputFileName.ReplaceAll(".root","");
  std::cout << "baseOutputFileName = " << baseOutputFileName.Data() << "\n";

  TString outputFileName;
  outputFileName.Form("%s_%s.txt",baseOutputFileName.Data(),histoName);
  std::cout << "outputFileName = " << outputFileName.Data() << "\n";

  
  TF1 * f = lxg.asTF( RooArgList(q) );
  double mpv = f->GetMaximumX();
  ofstream mpv_file;
  mpv_file.open(outputFileName.Data());
  mpv_file << landau_mpv.getVal() << "\t" << landau_mpv.getError() << "\n";
 
  mpv_file.close();




TPaveText* results_40V = new TPaveText(.4,.5,.7,.6,"NDC");
results_40V->SetBorderSize(0);
  TString mustring_40V;
//TString mustring_40V2;
 
 mustring_40V.Form("Convolution MPV = %1.2f ",mpv);
//mustring_40V2.Form("MPV2 irradiated %1.2f ToT",mpv2);
  
  results_40V->AddText(mustring_40V.Data());
 //results_40V->AddText(mustring_40V2.Data());
  
  results_40V->SetFillColor(0);
  results_40V->Draw("same");
  
  outputFileName.ReplaceAll(".txt",".pdf");
  c1->SaveAs(outputFileName.Data());


  //landau_mpv.Print();
}

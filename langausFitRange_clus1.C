#include "TVectorD.h"
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
void langausFitRange_clus1()
{
  using namespace RooFit;

  gSystem->Load("libRooFit.so");

  TCanvas *c1 = new TCanvas("c1","",1024,768);
  c1->cd();

TFile *file = new TFile("run000821-clustering-histo.root","READ");

 TH1F* clusterSignal_d20 = (TH1F*)  gDirectory->Get("ClusteringAPIX/detector_21/clusterSignal_d21");

  //clusterSignal_d20->Draw("e");
  int nentries = clusterSignal_d20->GetEntries();

clusterSignal_d20->Scale(1.0/nentries);
clusterSignal_d20->SetStats(111);  
// determining tot range
  int nbins = clusterSignal_d20->GetNbinsX();
  double xmax = 30;
  double xmin = clusterSignal_d20->GetBinLowEdge(1);
  if ( xmax < 0 ) { 
    double xmax = clusterSignal_d20->GetBinLowEdge(nbins+1);
  }
  double mean = clusterSignal_d20->GetMean();

    //  xmin += HitDiscCnfg;


  // the tot variable
  RooRealVar tot("tot","ToT [25ns]",mean,xmin,xmax);
  // the tot ranges
  tot.setRange("no16Low",1.5,15.5);
     tot.setRange("no16High",16.5,30.5);
  char fit_range[256];
     sprintf(fit_range,"%s,%s","no16Low","no16High");
  //sprintf(fit_range,"%s","no16Low");
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
  RooFitResult *fitr = lxg.fitTo(*totEDH,RooFit::Range(fit_range),Save()) ;
  // Plot data, landau pdf, landau (X) gauss pdf
  RooPlot* frame = tot.frame(RooFit::Title(" "),RooFit::Bins(nbins)) ;
  totEDH->plotOn(frame) ;
  //data_reduce->plotOn(frame) ;
  lxg.plotOn(frame) ;
  //lxg.paramOn(frame);
  // Draw frame on canvas
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;
  gPad->SetTicks(1,1);
  // chi2
//  RooChi2Var chi2_lowstat("chi2_lowstat","chi2",lxg,*totEDH) ;
//  cout << chi2_lowstat.getVal() << endl ;
  
  TF1 * f = lxg.asTF( RooArgList(tot) );
  double mpv = f->GetMaximumX();
ofstream mpv_file;
mpv_file.open("mpv_file.txt");
mpv_file << "MPV = " << mpv << "\n";



TPaveText* results_40V = new TPaveText(.6,.25,.8,.45,"NDC");
results_40V->SetBorderSize(0);
  TString mustring_40V;
//TString mustring_40V2;
 
 mustring_40V.Form("ToT MPV = %1.2f ",mpv);
//mustring_40V2.Form("MPV2 irradiated %1.2f ToT",mpv2);
  
  results_40V->AddText(mustring_40V.Data());
 //results_40V->AddText(mustring_40V2.Data());
  
  results_40V->SetFillColor(0);
  //results_40V->Draw("same");

 c1->SaveAs("output.root");
 c1->SaveAs("output.png");

 







  // Extract covariance and correlation matrix as TMatrixDSym
  const TMatrixDSym& cor = fitr->correlationMatrix() ;
  const TMatrixDSym& cov = fitr->covarianceMatrix() ;
  // Print correlation, covariance matrix
  std::cout << "correlation matrix" << "\n";
  cor.Print() ;
  std::cout << "covariance matrix" << "\n";
  cov.Print() ;

  // eigenvalues
  TMatrixDSymEigen ecov(cov);
  TVectorD evals_cov(ecov.GetEigenValues());
  evals_cov.Print();  
  // eigenvectors
  TMatrixD evecs_cov(ecov.GetEigenVectors());
  evecs_cov.Print();
  // extracting eigenvectors from evecs_cov (1 column per vector)
  TVectorD v1(3);
  v1[0] = evecs_cov[0][0]; 
  v1[1] = evecs_cov[1][0]; 
  v1[2] = evecs_cov[2][0]; 
  TVectorD v2(3);
  v2[0] = evecs_cov[0][1]; 
  v2[1] = evecs_cov[1][1]; 
  v2[2] = evecs_cov[2][1]; 
  TVectorD v3(3);
  v3[0] = evecs_cov[0][2]; 
  v3[1] = evecs_cov[1][2]; 
  v3[2] = evecs_cov[2][2]; 

  // a vector containing the values of the parameters determined by the fit
  TVectorD sol(3);
  sol[0] = landau_mean.getVal();
  sol[1] = landau_width.getVal();
  sol[2] = gaussian_sigma.getVal();

  // the square root of the eigenvalues
  double sigma1 = TMath::Sqrt(evals_cov[0]);
  double sigma2 = TMath::Sqrt(evals_cov[1]);
  double sigma3 = TMath::Sqrt(evals_cov[2]);

  //
  //  the sol vector is displaced by +/- each of the sigmas along the corresponding 
  //  eigenvector
  //  the resulting values are used to set the pdf parameters, which is then drawn
  //  and the new mpv value is obtained
  //

  // sol + sigma1*v1
  TVectorD lmp = sol+sigma1*v1;
  landau_mean.setVal(lmp[0]);
  landau_width.setVal(lmp[1]);
  gaussian_sigma.setVal(lmp[2]);
  lxg.plotOn(frame,LineColor(kRed)) ;
  RooFFTConvPdf lxg_plus_sigma1 = lxg;
  lxg_plus_sigma1.SetName("lxg_plus_sigma1");
  frame->Draw(); 
  TF1 * flmp = lxg.asTF( RooArgList(tot) );
  double mpv_lmp = flmp->GetMaximumX();
  flmp = lxg.asTF( RooArgList(tot) );
  mpv_lmp = flmp->GetMaximumX();
  // sol - sigma1*v1
  TVectorD lmm = sol-sigma1*v1;
  landau_mean.setVal(lmm[0]);
  landau_width.setVal(lmm[1]);
  gaussian_sigma.setVal(lmm[2]);
  lxg.plotOn(frame,LineColor(kBlack)) ;
  frame->Draw(); 
  TF1 * flmm = lxg.asTF( RooArgList(tot) );
  double mpv_lmm = flmm->GetMaximumX();
  flmm = lxg.asTF( RooArgList(tot) );
  mpv_lmm = flmm->GetMaximumX();
  // sol + sigma2*v2
  TVectorD lwp = sol+sigma2*v2;
  landau_mean.setVal(lwp[0]);
  landau_width.setVal(lwp[1]);
  gaussian_sigma.setVal(lwp[2]);
 lxg.plotOn(frame,LineColor(kRed),LineStyle(kDashed)) ;
 frame->Draw(); 
  TF1 * flwp = lxg.asTF( RooArgList(tot) );
  double mpv_lwp = flwp->GetMaximumX();
  flwp = lxg.asTF( RooArgList(tot) );
  mpv_lwp = flwp->GetMaximumX();
  // sol - sigma2*v2
  TVectorD lwm = sol-sigma2*v2;
  landau_mean.setVal(lwm[0]);
  landau_width.setVal(lwm[1]);
  gaussian_sigma.setVal(lwm[2]);
  lxg.plotOn(frame,LineColor(kBlack),LineStyle(kDashed)) ;
  frame->Draw(); 
  TF1 * flwm = lxg.asTF( RooArgList(tot) );
  double mpv_lwm = flwm->GetMaximumX();
  flwm = lxg.asTF( RooArgList(tot) );
  mpv_lwm = flwm->GetMaximumX();
  // sol + sigma3*v3
  TVectorD gwp = sol+sigma3*v3;
  landau_mean.setVal(gwp[0]);
  landau_width.setVal(gwp[1]);
  gaussian_sigma.setVal(gwp[2]);
 lxg.plotOn(frame,LineColor(kRed),LineStyle(kDotted)) ;
 frame->Draw(); 
  TF1 * fgwp = lxg.asTF( RooArgList(tot) );
  double mpv_gwp = fgwp->GetMaximumX();
  fgwp = lxg.asTF( RooArgList(tot) );
  mpv_gwp = fgwp->GetMaximumX();
  // sol - sigma3*v3
  TVectorD gwm = sol-sigma3*v3;
  landau_mean.setVal(gwm[0]);
  landau_width.setVal(gwm[1]);
  gaussian_sigma.setVal(gwm[2]);
  lxg.plotOn(frame,LineColor(kBlack),LineStyle(kDotted)) ;
  frame->Draw(); 
  TF1 * fgwm = lxg.asTF( RooArgList(tot) );
  double mpv_gwm = fgwm->GetMaximumX();
  fgwm = lxg.asTF( RooArgList(tot) );
  mpv_gwm = fgwm->GetMaximumX();

  double shifts[6] = {mpv_lmp-mpv,mpv_lmm-mpv,mpv_lwp-mpv,mpv_lwm-mpv,mpv_gwp-mpv,mpv_gwm-mpv};
  for (int i = 0; i<6; i++ ) {
    std::cout << "shift[" << i << "] = "  << shifts[i] << "\n";
  }
  std::cout << "MPV        = " << mpv << "\n";
  std::cout << "MPV+sigma1 = " << mpv_lmp << "\n";
  std::cout << "MPV-sigma1 = " << mpv_lmm << "\n";
  std::cout << "MPV+sigma2 = " << mpv_lwp << "\n";
  std::cout << "MPV-sigma2 = " << mpv_lwm << "\n";
  std::cout << "MPV+sigma3 = " << mpv_gwp << "\n";
  std::cout << "MPV-sigma3 = " << mpv_gwm << "\n";

  double sym_errors[3] = {(fabs(shifts[0])+fabs(shifts[1]))/2.0,(fabs(shifts[2])+fabs(shifts[3]))/2.0,(fabs(shifts[4])+fabs(shifts[5]))/2.0}; 
  for (int i = 0; i<3; i++ ) {
    std::cout << "Symmetric error #" << i << " = " << sym_errors[i] << "\n";
  }
  double total_error = TMath::Sqrt(sym_errors[0]*sym_errors[0]+sym_errors[1]*sym_errors[1]+sym_errors[2]*sym_errors[2]);
  std::cout << "MPV = (" << mpv << " +/- " << total_error << ") ToT\n";
  mpv_file << "MPVerror = " << total_error << "\n";
  mpv_file.close();

  
  TLegend *leg = new TLegend(.6,.5,.75,.85);
  clusterSignal_d20->SetMarkerStyle(20);
  clusterSignal_d20->SetMarkerSize(1.2);
  leg->AddEntry(clusterSignal_d20,"Data","pe");
  
  f->SetLineColor(kBlue);
  leg->AddEntry(f,"Fit","l");

  flmp->SetLineColor(kRed);
  leg->AddEntry(flmp,"Fit+#sigma_{1}","l");
  
  flmm->SetLineColor(kBlack);
  leg->AddEntry(flmm,"Fit-#sigma_{1}","l");
  
  flwp->SetLineColor(kRed);
  flwp->SetLineStyle(kDashed);
  leg->AddEntry(flwp,"Fit+#sigma_{2}","l");

  flwm->SetLineColor(kBlack);
  flwm->SetLineStyle(kDashed);
  leg->AddEntry(flwm,"Fit-#sigma_{2}","l");
  
  fgwp->SetLineColor(kRed);
  fgwp->SetLineStyle(kDotted);
  leg->AddEntry(fgwp,"Fit+#sigma_{3}","l");

  fgwm->SetLineColor(kBlack);
  fgwm->SetLineStyle(kDotted);
  leg->AddEntry(fgwm,"Fit-#sigma_{3}","l");


  leg->Draw();
  mustring_40V.Form("ToT MPV = (%1.1f +/- %1.1f) ",mpv,total_error);
  results_40V->DeleteText();
  results_40V->AddText(mustring_40V.Data());
  results_40V->Draw("same");

}

#include <iostream>
#include "TGraph.h"
#include "TMath.h" 
#include "TF1.h"
#include "TCanvas.h"

const double vsat = 2e-2; // cm/ns

double fitfunc(double *x, double *par) {
  double phi = x[0]; // cm^-2
  double beta = par[0]; // ns^-1*cm^2
  double w = par[1]; // mu
  double tau = 1./beta/phi; // ns
  double lambda = vsat*tau*1e4; // mu

  double expo = w/lambda;
  double cce = 1/expo*(1-TMath::Exp(-expo));

  return cce;
}

void CCE_vs_fluence(const char* fileName) {
  TGraph* g = new TGraph(fileName);
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->cd();
  TF1* f = new TF1("f",fitfunc,1e12,1e16,2);
  f->SetParameters(1e-16,200);
  f->SetParLimits(0,1e-19,1e-13);
  f->FixParameter(1,200);
  g->SetMarkerSize(2);
  g->SetMarkerStyle(24);
  g->Draw("AP");
  g->Fit("f");
  double b = f->GetParameter(0);
  double eb = f->GetParError(0);
  std::cout << "beta = ( " << b << " +/- " << eb << " ) cm^2 ns^-1\n";
}


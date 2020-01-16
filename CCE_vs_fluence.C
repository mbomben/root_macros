#include <iostream>
#include "TGraph.h"
#include "TMath.h" 
#include "TF1.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TPaveText.h"

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

void CCE_vs_fluence(const char* fileName, double w=200) {
  TGraph* g = new TGraph(fileName);
  g->SetTitle(";#Phi [n_{eq}/cm^{2}];CCE");
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->cd();
  TF1* f = new TF1("f",fitfunc,1e12,1e16,2);
  f->SetParameters(1e-16,300);
  //f->SetParLimits(0,1e-19,1e-13);
  f->FixParameter(1,w);
  g->SetMarkerSize(2);
  g->SetMarkerStyle(24);
  g->Draw("AP");
  g->Fit("f");
  double b = f->GetParameter(0);
  double eb = f->GetParError(0);
  std::cout << "beta = ( " << b << " +/- " << eb << " ) cm^2 ns^-1\n";
  TGaxis::SetMaxDigits(3);
  c1->SetGrid(1);
  c1->SetTicks(1,1);


  TPaveText *pt = new TPaveText(3e15,0.8,5e15,0.9);
  TString text;
  text = TString::Format("#beta = (%2.1e +/- %2.1e) cm^{2}/ns",b,eb);
  pt->AddText(text.Data());
  pt->Draw();
  TString saveFile(fileName);
  saveFile.ReplaceAll(".txt",".png");
  c1->SaveAs(saveFile.Data());
  saveFile = fileName;
  saveFile.ReplaceAll(".txt",".pdf");
  c1->SaveAs(saveFile.Data());


}


#include <iostream>
#include "TGraph.h"
#include "TMath.h" 
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <fstream>
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"

TCanvas* vdepl_vs_fluence( const char* fileName ) {

  std::ifstream file;
  file.open(fileName);
  if (!file.good()) { std::cout << "Error opening file\n"; exit(1); }

  double phi,vdepl,evdepl;
  const int nmax = 20;
  double *X = new double[nmax];
  double *Y = new double[nmax];
  double *EX = new double[nmax];
  double *EY = new double[nmax];
  int i = 0;
  while(file >> phi >> vdepl >> evdepl) {
    X[i]=phi;
    Y[i]=vdepl;
    EX[i]=0;
    EY[i]=evdepl;
    i++;
    if( i>= nmax ) break;
  }
  file.close();
  
  
  TGraphErrors* g = new TGraphErrors(i,X,Y,EX,EY);
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->cd();
  TF1 *f = new TF1("f","pol1",0,10);
  f->SetLineColor(kRed);
  g->SetTitle(";#Phi [1#times10^{15} n_{eq}/cm^{2}];V_{depl} [V]");
  g->SetName("vdepl_vs_fluence");
  g->SetMarkerStyle(24);
  g->SetMarkerSize(1);
  g->Draw("APE");
  g->Fit("f");
  //g->Draw("AE");
  gPad->Modified(); gPad->Update();  
  std::cout << "here\n"; 
  //g->Draw("AE");

  c1->SetGrid(1);
  c1->SetTicks(1,1);
  double slope  = f->GetParameter(1);
  double eslope = f->GetParError(1);

  TPaveText *pt = new TPaveText(1,500,3,600);
  TString text;
  text = TString::Format("#frac{dV_{depl}}{d#Phi} = (%2.2f +/- %2.2f) #frac{V}{1#times10^{15} n_{eq}/cm^{2}}",slope,eslope);
  pt->AddText(text.Data());
  pt->Draw();
  TString saveFile(fileName);
  saveFile.ReplaceAll(".txt",".png");
  c1->SaveAs(saveFile.Data());
  saveFile = fileName;
  saveFile.ReplaceAll(".txt",".pdf");
  c1->SaveAs(saveFile.Data());

  return c1;

}

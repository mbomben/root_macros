#include <fstream>
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TF1.h"


void CV(const char* fileName, const char* type, double A, double v1, double v2, double v3, double v4, double &vdepl, double &evdepl) {
  TString simsstring("SIMU");
  TString datastring("DATA");
  int NMAX = 500;
  if ( ! (simsstring.EqualTo(type) || datastring.EqualTo(type) ) ) {
    std::cerr << "type must be either SIMU or DATA, not -> " << type << "\n";
    exit(2);
  }
  double *C = new double[NMAX];
  double *V = new double[NMAX];
  double *logC = new double[NMAX];
  double *logV = new double[NMAX];
  double *C2 = new double[NMAX];
  double c,v;
  int i = 0;
  std::ifstream file;
  file.open(fileName);
  if ( !file.good() ) {
    std::cerr << "Problems opening file " << fileName << "\n";
    std::cerr << "rdstate =             " << file.rdstate() << "\n";
    exit(3);
  }

  while(1) {
    file >> v >> c;
    if (file.eof()) break;
    if (! (c>0) ) continue; 
    c=fabs(c);
    v=fabs(v);
    C[i]=c;
    V[i]=v;
    logC[i]=TMath::Log10(c);
    logV[i]=TMath::Log10(v);
    C2[i]=1./c/c;
    i++;
    if ( i > NMAX ) {
      std::cerr << "Too many lines: " << i << "\n";
      std::cerr << "Maximum is :    " << NMAX << "\n";
      exit(4);
    }

  } 

  TGraph* grCV = new TGraph(i,V,C);
  grCV->SetName("grCV");
  grCV->SetTitle("");
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->cd();
  grCV->SetMarkerColor(kBlue);
  grCV->SetMarkerStyle(24);
  grCV->SetMarkerSize(1.2);
  grCV->SetLineColor(kBlue);
  grCV->Draw("AP");
  grCV->GetYaxis()->SetTitle("C [F]");
  grCV->GetXaxis()->SetTitle("V_{bias} [V]");

  TGraph* grlogClogV = new TGraph(i,logV,logC);
  grlogClogV->SetName("grlogClogV");
  grlogClogV->SetTitle("");
  grlogClogV->SetMarkerColor(kBlack);
  grlogClogV->SetMarkerStyle(24);
  grlogClogV->SetMarkerSize(1.2);
  grlogClogV->SetLineColor(kBlack);
  grlogClogV->Draw("AP");
  grlogClogV->GetYaxis()->SetTitle("log_{10}(C/F)");
  grlogClogV->GetXaxis()->SetTitle("log_{10}(V_{bias}/V)");

  double logv1 = TMath::Log10(v1);
  double logv2 = TMath::Log10(v2);
  double logv3 = TMath::Log10(v3);
  double logv4 = TMath::Log10(v4);

  TF1 *loglin1 = new TF1("loglin1","[0]+[1]*x",logv1*0.8,logv2*1.2);
  TF1 *loglin2 = new TF1("loglin2","[0]+[1]*x",logv3*0.8,logv4*1.2);
  loglin1->SetLineColor(kRed);
  loglin2->SetLineColor(kBlue);
  grlogClogV->Fit("loglin1","R","",logv1,logv2);
  grlogClogV->Fit("loglin2","R+","",logv3,logv4);

  double logq1 = loglin1->GetParameter(0);
  double logq2 = loglin2->GetParameter(0);
  double logm1 = loglin1->GetParameter(1);
  double logm2 = loglin2->GetParameter(1);

  double elogq1 = loglin1->GetParError(0);
  double elogq2 = loglin2->GetParError(0);
  double elogm1 = loglin1->GetParError(1);
  double elogm2 = loglin2->GetParError(1);

  assert((logm1-logm2)!=0);
  double logvdep = (logq1-logq2)/(logm2-logm1);
  vdepl = TMath::Power(10.,logvdep);
  std::cout << "vdepl = " << vdepl << " V \n";
  double Deltam = sqrt(TMath::Power(elogm1,2.)+TMath::Power(elogm2,2.));
  double Deltaq = sqrt(TMath::Power(elogq1,2.)+TMath::Power(elogq2,2.));

  evdepl = sqrt(TMath::Power(Deltaq,2.)/TMath::Power((logm2-logm1),2.) + TMath::Power(Deltam,2.)*TMath::Power((logq1-logq2),2.)/TMath::Power((logm2-logm1),2.));

  evdepl = TMath::Power(10.,evdepl);
  std::cout << "evdepl = " << evdepl << " V \n";


  file.close();
} 

#include <fstream>
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TF1.h"

const double eR = 11.9;     // Silicon relative dielectric constant
const double e0 = 8.85e-14; // F/cm
const double q0 = 1.6e-19;  // C

void CV(const char* fileName, const char* type, double A, double v1, double v2, double v3, double v4, double &vdepl, double &evdepl, double &neff, double &eneff, double &w, double &ew) {

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
  double Deltam = sqrt(TMath::Power(elogm1,2.)+TMath::Power(elogm2,2.));
  double Deltaq = sqrt(TMath::Power(elogq1,2.)+TMath::Power(elogq2,2.));

  evdepl = sqrt(TMath::Power(Deltaq,2.)/TMath::Power((logm2-logm1),2.) + TMath::Power(Deltam,2.)*TMath::Power((logq1-logq2),2.)/TMath::Power((logm2-logm1),2.));

  evdepl = TMath::Power(10.,evdepl);

  TGraph *grC2V = new TGraph(i,V,C2);
  grC2V->SetName("grC2V");
  grC2V->SetTitle("");
  grC2V->SetMarkerColor(kBlue);
  grC2V->SetMarkerStyle(24);
  grC2V->SetMarkerSize(1.2);
  grC2V->SetLineColor(kBlue);
  grC2V->Draw("AP");
  grC2V->GetYaxis()->SetTitle("C^{-2} [F^{-2}]");
  grC2V->GetXaxis()->SetTitle("V_{bias} [V]");

  grC2V->Fit("pol1","R","",v1,v2);
  TF1 *lin = (TF1*)gROOT->GetFunction("pol1");
  lin->SetLineColor(kRed);

  double C2der = lin->GetParameter(1);
  double eC2der = lin->GetParError(1);

  neff = 2./A/A/q0/eR/e0/C2der;
  eneff = 2./A/A/q0/eR/e0/C2der/C2der*eC2der;

  w = TMath::Power(2.*eR*e0*vdepl/q0/neff,0.5);

  double ewA = 2*eR*e0/q0;
  ewA = TMath::Power(ewA,0.5);

  double ewNeff = 0.5*ewA*TMath::Power(vdepl,0.5)*TMath::Power(neff,-1.5)*eneff;
  double ewV = 0.5*ewA*TMath::Power(vdepl,-0.5)*TMath::Power(neff,-0.5)*evdepl;

  ew = TMath::Sqrt(ewNeff*ewNeff+ewV*ewV);

  std::cout << "v1 = " << v1 << " V\n"; 
  std::cout << "v2 = " << v2 << " V\n"; 
  std::cout << "v3 = " << v3 << " V\n"; 
  std::cout << "v4 = " << v4 << " V\n"; 

  std::cout << "vdepl = " << vdepl << " V \n";
  std::cout << "evdepl = " << evdepl << " V \n";

  std::cout << "neff = " << neff << " 1./cm^3\n";
  std::cout << "eneff = " << eneff << " 1./cm^3\n";
  std::cout << "w = " << w*1e+4 << " um\n";
  std::cout << "ew = " << ew*1e+4 << " um\n";

  file.close();
} 

#include <fstream>
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLine.h"
#include "TObjString.h"

const double eR = 11.9;     // Silicon relative dielectric constant
const double e0 = 8.85e-14; // F/cm
const double q0 = 1.6e-19;  // C

double rescale( double Iorig, double T) {
  const double offset = 273.15; // K
  double tref = 20.0; // ºC
  double Tref = tref + offset;

  const double kboltz = 8.617343e-5; // eV K−1
  const double egap = 1.21; //eV

  double scale = pow((Tref/T),2) * exp( (-egap/(2.*kboltz)) * (1./Tref -1./T));



  double scaledCurrent = scale*Iorig;

  return scaledCurrent;
}

       

void IV(const char* fileName, const char* type, double A, double w, double v1, double T, double fluence, double &alpha, double dV=50, double I0=0) {

  assert(fluence>0);
  TString simsstring("SIMU");
  TString datastring("DATA");
  int NMAX = 1000;
  if ( ! (simsstring.EqualTo(type) || datastring.EqualTo(type) ) ) {
    std::cerr << "type must be either SIMU or DATA, not -> " << type << "\n";
    exit(2);
  }
  double *I = new double[NMAX];
  double *V = new double[NMAX];
  double i,v;
  int j = 0;
  std::ifstream file;
  file.open(fileName);
  if ( !file.good() ) {
    std::cerr << "Problems opening file " << fileName << "\n";
    std::cerr << "rdstate =             " << file.rdstate() << "\n";
    exit(3);
  }
  while(1) {
    file >> v >> i;
    if (file.eof()) break;
    i=fabs(i);
    v=fabs(v);
    I[j]=i;
    V[j]=v;
    j++;
    if ( j > NMAX ) {
      std::cerr << "Too many lines: " << j << "\n";
      std::cerr << "Maximum is :    " << NMAX << "\n";
      exit(4);
    }

  } 
  TGraph* grIV = new TGraph(j,V,I);
  grIV->SetName("grIV");
  grIV->SetTitle("");
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->cd();
  grIV->SetMarkerColor(kBlue);
  grIV->SetMarkerStyle(24);
  grIV->SetMarkerSize(1.2);
  grIV->SetLineColor(kBlue);
  grIV->Draw("AP");
  grIV->GetYaxis()->SetTitle("I [A]");
  grIV->GetXaxis()->SetTitle("V_{bias} [V]");
  c1->SetTicks(1,1);

  double Ivdepl = grIV->Eval(v1+dV);
  Ivdepl -= I0;
  double Iscaled = rescale(Ivdepl,T);
  double currDens = Iscaled/A/w;
  alpha = currDens/fluence;      
  std::cout << "alpha = " << alpha << " A/cm\n";
 
  TLine *l = new TLine(v1+dV,grIV->GetYaxis()->GetXmin(),v1+dV,grIV->GetYaxis()->GetXmax());
  l->SetLineColor(kRed);
  l->SetLineWidth(2);
  l->SetLineStyle(kDashed);
  l->Draw();

  gPad->SetLogy(1);
 
  TString saveFile(fileName);
  TString globalName =  ((TObjString*)saveFile.Tokenize(".")->At(0))->GetString();
  saveFile.Form("%s.png",globalName.Data());
  c1->SaveAs(saveFile.Data());
  saveFile.Form("%s.pdf",globalName.Data());
  c1->SaveAs(saveFile.Data());
  
}

void IVanalysis(const char* fileName, const char* type, double A, double w, double v1, double eV, double T, double fluence, double &alpha, double dV=50, double I0=0) {
  
  double v;
  double alphas[3];
  v = v1 + eV;
  IV(fileName, type, A, w, v, T, fluence, alpha);
  alphas[0] = alpha;
  v = v1 - eV;
  IV(fileName, type, A, w, v, T, fluence, alpha);
  alphas[1] = alpha;
  v = v1;
  IV(fileName, type, A, w, v, T, fluence, alpha);
  alphas[2] = alpha;


  alpha  = TMath::Mean(3,alphas);
  double ealpha = TMath::RMS(3,alphas);

  std::cout << "alpha = " << alpha << " A/cm\n";
  std::cout << "ealpha = " << ealpha << " A/cm\n";
  std::cout << "alpha ealpha = " << alpha << " " << ealpha << " A/cm\n";
}


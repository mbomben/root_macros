using namespace RooFit;
using namespace RooStats;

/*double fitFunc( double *x, double *par ) {
  double x1= x[0];
  double x2= x[1];

  double mu1= par[0];
  double mu2= par[1];
  double  s1= par[2];
  double  s2= par[3];
  double rho= par[4];

  TMatrixD xV(2,1);
  xV[0][0] = x1;
  xV[1][0] = x2;

  TMatrixD V(2,2); // covariance matrix
  V[0][0] = s1*s1;
  V[1][1] = s2*s2;
  V[0][1] = s1*s2*rho;
  V[1][0] = V[0][1];

  TMatrixD muV(2,1);
  muV[0][0] = mu1;
  muV[1][0] = mu2;

  TMatrixD dV = xV-muV;

  TMatrixD dVT(2,1);
  dVT.Transpose(dV);

  TMatrixD IV = V;
  IV = IV.Invert();

  double Vdeterminant = V.Determinant();

  TMatrixD prod1  = IV * dV;
  TMatrixD prod2  = dVT * prod1;
  double   arg    = prod2[0][0];
           arg    = -1./2.*arg;

  double exponent = TMath::Exp(arg);

  double value    = exponent / 2. / TMath::Pi() / TMath::Sqrt(Vdeterminant);

  return value;

}*/

void MultiVariate( double mu1v = 0, double mu2v = 1, double s1v = 1, double s2v = 0.5, double rhov = 0.2, const char* drawOpt = "colz") {


  RooArgList xVec;
  RooArgList muVec;
  int dim = 2;

  // observables
  RooRealVar *x1 = new RooRealVar("x1","x1",0,-5,5);
  RooRealVar *x2 = new RooRealVar("x2","x2",0,-5,5);
  xVec.add(*x1);
  xVec.add(*x2);

  // means
  RooRealVar *mu1 = new RooRealVar("mu1","mu1",mu1v);
  RooRealVar *mu2 = new RooRealVar("mu2","mu2",mu2v);
  muVec.add(*mu1);
  muVec.add(*mu2);

  // covariance matrix
  TMatrixDSym cov(dim);
  cov[0][0] = s1v*s1v;
  cov[1][1] = s2v*s2v;
  cov[1][0] = s1v*s2v*rhov;
  cov[0][1] = cov[1][0];

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvg("mvg", "mvg", xVec, muVec, cov);

  // let's create data
  //RooDataSet *data2 = mvg.generate(RooArgSet(*x1, *x2), 100000);
  // let's get a 2D plot of them via TH2
  TCanvas *c1 = new TCanvas("c1","");
  c1->cd();
  TH2 *mvgplot = (TH2 *)mvg.createHistogram("x1,x2", 100, 100);
  mvgplot->SetTitle("Multivariate Gaussian");
  mvgplot->Draw(drawOpt);
  mvgplot->SetContour(999);
  c1->SetRightMargin(0.15);
  c1->SetTicks();
  mvgplot->GetZaxis()->SetTitle("Prob.");
  mvgplot->GetZaxis()->SetTitleOffset(1.5);

  gStyle->SetOptStat(0);

  c1->SaveAs("MultiVariate.pdf");

  *x1 = 4.;
  RooPlot *x2frame = x2->frame();
  mvg.plotOn(x2frame);
  x2frame->Draw();

  mvg.plotOn(x2frame,Project(*x1),LineStyle(kDashed),LineColor(kRed));
  x2frame->Draw();
  c1->SaveAs("sliceMultiVariate.pdf");
}

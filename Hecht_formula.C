
const double vsate = 1.08e7;  // cm/s
const double vsath = 1.08e7;  // cm/s

double func( double *x, double *par ) {
  double fluence = x[0];

  double betae = par[0];
  double w     = par[1];

  w *= 1e-4; // um -> cm conversion
  double de = vsate/fluence/betae; // average mean free path for electrons

  double lin = (de)/w; // linear term

  double quade = (de/w)*(de/w); // quadratic factor for electrons

  double expe = (1.0-TMath::Exp(-w/de)); // exponential factor for electrons

  double cce = 2.0*(lin - quade*expe);

  return cce;
}


double Hecht_formula( const char* fileName, double w = 130, double betae = 3.5e-7, double flmin=1, double flmax=3e16 ) {

  TGraphErrors *gre = new TGraphErrors(fileName);
  gre->SetMarkerStyle(24);
  gre->SetMarkerSize(1.4);
  gre->SetMarkerColor(3);

  gre->SetTitle(";#Phi [n_{eq}/cm^{2}];CCE");
  gre->Draw("APE");
  
  TF1* f = new TF1("f",func,1,2e16,2);
  f->FixParameter(1,w);
  f->SetParameter(0,betae);

  gre->Fit("f","R","",flmin,flmax);

  betae = f->GetParameter(0);

  return betae;
}


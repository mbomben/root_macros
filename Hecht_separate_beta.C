
const double vsate = 1.08e7;  // cm/s
const double vsath = 1.08e7;  // cm/s

double func( double *x, double *par ) {
  double fluence = x[0];

  double betae = par[0];
  double betah = par[1];
  double w     = par[2];

  w *= 1e-4; // um -> cm conversion
  double de = vsate/fluence/betae; // average mean free path for electrons
  double dh = vsath/fluence/betah; // average mean free path for holes

  double line = (de)/w; // linear term
  double linh = (dh)/w; // linear term

  double quade = (de/w)*(de/w); // quadratic factor for electrons
  double quadh = (dh/w)*(dh/w); // quadratic factor for electrons

  double expe = (1.0-TMath::Exp(-w/de)); // exponential factor for electrons
  double exph = (1.0-TMath::Exp(-w/dh)); // exponential factor for holes

  double cce = (line+linh) - quade*expe - quadh*exph;

  return cce;
}


double Hecht_formula( const char* fileName, double w = 130, double betae = 3.5e-7, double betah = 6.5e-7,  double flmin=1, double flmax=3e16 ) {

  TGraphErrors *gre = new TGraphErrors(fileName);
  gre->SetMarkerStyle(24);
  gre->SetMarkerSize(1.4);
  gre->SetMarkerColor(4);

  gre->SetTitle(";#Phi [n_{eq}/cm^{2}];CCE");
  gre->Draw("APE");
  gre->GetYaxis()->SetRangeUser(0,1.1);
  
  TF1* f = new TF1("f",func,1,2e16,3);
  f->SetParameter(0,betae);
  f->SetParameter(1,betah);
  f->FixParameter(2,w);

  gre->Fit("f","R","",flmin,flmax);

  betae = f->GetParameter(0);
  betah = f->GetParameter(1);

  double ebetae = f->GetParError(0);
  double ebetah = f->GetParError(1);

  std::cout << "betae = " << betae << " +/- " << ebetae << "\n";
  std::cout << "betah = " << betah << " +/- " << ebetah << "\n";

  gPad->SetTicks(1,1);
  
  TString sbetae;
  sbetae.Form("#beta_{e} = (%2.1e #pm %2.0e) cm^{2}/ns",betae/1e9,ebetae/1e9);
  TString sbetah;
  sbetah.Form("#beta_{h} = (%2.1e #pm %2.0e) cm^{2}/ns",betah/1e9,ebetah/1e9);
  std::cout << sbetae.Data() << "\n";
  std::cout << sbetah.Data() << "\n";

  TPaveText *tt = new TPaveText(.4,.6,.85,.85,"NDC");
  tt->AddText(sbetae.Data());
  tt->AddText(sbetah.Data());
  tt->Draw();

  TString baseName(fileName);
  TString pngName = baseName;
  TString pdfName = baseName;
  pngName.ReplaceAll(".txt",".png");
  gPad->SaveAs(pngName.Data());
  pdfName.ReplaceAll(".txt",".pdf");
  gPad->SaveAs(pdfName.Data());

  return f->GetChisquare()/f->GetNDF();
}


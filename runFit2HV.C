void runFit2HV(const char* fileName, float &rms, float prec = 0.02) 
{

  gROOT->Reset();
  gROOT->LoadMacro("Fit2HV.C+");

  TGraphErrors *gre = new TGraphErrors(fileName);

  Fit2HV(gre,prec,rms);

  gre->Draw("A*E");
  f1->Draw("LSAME");
  f2->Draw("LSAME");
}
  

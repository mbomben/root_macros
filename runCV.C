{
  gROOT->Reset();
  gSystem->CompileMacro("CV.C");
  gSystem->Load("CV_C.so");
  double v,e;
  CV("/Users/bomben/work/RD50/Natascha/traptunnel_newPetaseccaLow_CV/ninp_fei4_3pixels_fluence\=2e+15_bias\=2000_RD50_CV_pruned.dat","SIMU",1,10**1.6,10**2.1,10**2.1,10**2.7,v,e);
}

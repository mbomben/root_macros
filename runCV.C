{
  gROOT->Reset();
  gSystem->CompileMacro("CV.C");
  gSystem->Load("CV_C.so");
  double v,ev,neff,eneff,w,ew;
  CV("/Users/bomben/work/RD50/Natascha/traptunnel_newPetaseccaLow_CV/ninp_fei4_3pixels_fluence\=2e+15_bias\=2000_RD50_CV_pruned.dat","SIMU",5e-7,10**1.6,10**2.1,10**2.1,10**2.7,v,ev,neff,eneff,w,ew);
}

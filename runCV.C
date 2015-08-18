{
  gROOT->Reset();
  gSystem->CompileMacro("CV.C");
  gSystem->Load("CV_C.so");
  double v,ev,neff,eneff,w,ew;
  CV("ninp_fei4_3pixels_fluence=2e+15_bias=2000_RD50_CV_pruned.dat","SIMU",5e-7,40,110,110,280,v,ev,neff,eneff,w,ew);
}

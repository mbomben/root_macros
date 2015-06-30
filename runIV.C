{
  gROOT->Reset();
  gSystem->CompileMacro("IV.C");
  gSystem->Load("IV_C.so");
  double alpha;
  IV("/Users/bomben/work/RD50/Natascha/newPetaseccaLow/ninp_fei4_3pixels_fluence\=2e+15_bias\=2000_RD50_PX2_Current_pruned.dat","SIMU",50e-8,200e-4,130,253.15,2e15,alpha);
}

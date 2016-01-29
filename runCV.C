{
  gROOT->Reset();
  gROOT->LoadMacro("CV.C+");
  //gSystem->CompileMacro("CV.C");
  //gSystem->Load("CV_C.so");
  double v,ev,neff,eneff,w,ew;
  //CV("ninp_fei4_3pixels_fluence=2e+15_bias=2000_RD50_CV_pruned.dat","SIMU",5e-7,70,120,316,1000,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/temperature_studies/T20/small_steps_fl0_plateau/ninp_fei4_3pixels_fluence=0_bias=500_RD50_CV_pruned.dat","SIMU",5e-7,34, 60,100, 400,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/temperature_studies/T20/no_mult/ninp_fei4_3pixels_fluence=0_bias=500_RD50_CV_pruned.dat","SIMU",5e-7,4, 25,60, 100,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w50/small_steps_fl0_plateau/ninp_fei4_3pixels_fluence=0_bias=1000_RD50_CV_pruned.dat","SIMU",1.5e-6,33, 36,60, 100,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w50/small_steps_Petasecca_fl1e15_plateau/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,30, 40,60, 100,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/gain_studies/w200/small_steps_fl0_plateau/ninp_fei4_3pixels_fluence=0_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,30, 55,71, 100,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,20, 45,81, 200,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_traptunnel_fl1e15/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,20, 45,81, 200,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_interface_fl1e15/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,1, 10,81, 150,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_interface_traptunnel_fl1e15/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 40,81, 150,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 100,150, 200,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.08/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,150, 200,300, 500,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.09/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,110, 170,300, 500,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.10/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,50, 120,300, 500,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 100,300, 500,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.12/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 80,300, 500,v,ev,neff,eneff,w,ew,true);
  //CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.13/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 60,300, 500,v,ev,neff,eneff,w,ew,true);
  CV("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Petasecca_fl1e15_T20_Eg1.14/ninp_fei4_3pixels_fluence=1e+15_bias=500_RD50_CV_pruned.dat","SIMU",1.5e-6,10, 50,300, 500,v,ev,neff,eneff,w,ew,true);
}

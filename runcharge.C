{
  gROOT->Reset();
  gROOT->LoadMacro("charge.C+");
  double q = charge("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/no_mult_Perugia_fl0_Synopsys_Eg108/MIP-100V/ninp_fei4_3pixels_fluence=0_RD50-SEU_PAD_Current_pruned.dat",2e-11);

  std::cout << q << "\n";
}

{
  gROOT->Reset();
  gROOT->LoadMacro("vdepl_vs_fluence.C+");
  vdepl_vs_fluence("vdepl_Perugia2015.txt");
}
    

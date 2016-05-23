{
    gROOT->Reset();
    gROOT->LoadMacro("CCE_vs_fluence.C+");
    CCE_vs_fluence("/Users/mbomben/work/RD50/LGAD/thickness_studies/w200/cce_no_mult_Perugia_Synopsys_Eg108_900V.txt");
}

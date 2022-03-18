void makeFakeToTHistos() {
  TH1D* h1 = new TH1D("h1",";ToT;Entries",21,-0.5,20.5);
  TH1D* h2 = new TH1D("h2",";ToT;Entries",21,-0.5,20.5);

  TF1* f1 = new TF1("f1","[0]*TMath::Gaus(x,[1],[2],1)",0,20);
  TF1* f2 = new TF1("f2","[0]*TMath::Gaus(x,[1],[2],1)",0,20);

  f1->SetParameters(1,5,3);
  f2->SetParameters(1,4,4);

  h1->FillRandom("f1",19000);
  h2->FillRandom("f2",17000);

  h1->SetLineColor(kBlue);
  h1->Draw();
  h2->SetLineColor(kRed);
  h2->Draw("SAME");

  TFile* fFake = new TFile("fFake.root","RECREATE");
  fFake->cd();
  h1->Write();
  h2->Write();

  fFake->Close();
}


void ratioplot_ITkPixelRDOAnalysis() {
  string filename_f1 = "fFake.root";
  //string filename_f1 = "/Users/nakkalil/cernbox/athena_clones/21.9_MyAthenaWorkingFolder/run/PixelRDOAnalysis.root";
  //string filename_f2 = "/Users/nakkalil/cernbox/athena_clones/Athena08022022/run/1000_Evnts_ITkPixelRDOAnalysis.root";

  TFile* f1 = TFile::Open(filename_f1.c_str());
  //TFile* f2 = TFile::Open(filename_f2.c_str());

  string hist1 = "h1"; // h_belowThresh_ec,h_ToT,h_brlToT,h_ecToT,h_belowThresh_brl
  string hist2 = "h2";

  TH1F* h1 = static_cast<TH1F*>(f1->Get(hist1.c_str()));
  TH1F* h2 = static_cast<TH1F*>(f1->Get(hist2.c_str()));

   auto c1 = new TCanvas();
  
  //gStyle->SetOptStat(0);

  //  h1->SetTitle("");
  h1->SetLineColor(kBlue);
  // h1->GetXaxis()->SetRangeUser(0,20);
  // h1->GetYaxis()->SetRangeUser(0,3000);
  h1->Draw();
 

  h2->SetLineColor(kRed);
  //  h2->GetXaxis()->SetRangeUser(0,20);
  // h2->GetYaxis()->SetRangeUser(0,13000);
  h2->Draw("sames");

  TLegend* leg1 = new TLegend(0.6,0.7,0.9,0.9);
  leg1->AddEntry(h1,"R21.9");
  leg1->AddEntry(h2,"R22");
  leg1->Draw();

  auto rp = new TRatioPlot(h1,h2);
  c1->SetTicks(0,1);
  rp->SetH2DrawOpt("hist");
  rp->Draw();
  c1->Update();
 
  rp->GetUpperPad()->cd();
  h2->Draw("sames");
  //h2->Draw("][sames");
  rp->GetUpperPad()->Update();
  TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  ps1->SetX1NDC(0.65); ps1->SetX2NDC(0.85);
  ps1->SetTextColor(kBlue);
  ps1->SetLineColor(kBlue);
  TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);
  ps2->SetY1NDC(0.25); ps2->SetY2NDC(0.45);
  ps2->SetTextColor(kRed);
  ps2->SetLineColor(kRed);
}

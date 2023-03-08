void evaluate_EZ_asymmetry(const char* fileName="new_vs_old_E_FieldZ_asymmetry.txt") {
	TTree *t  = new TTree("t","");
	t->ReadFile("new_vs_old_E_FieldZ_asymmetry.txt","x:y:z:asymm");
	TCanvas *c1 = new TCanvas("c1","",1200,800);
	c1->cd();
	int nevts = t->GetEntries();
	t->Draw("asymm");
	gPad->SetLogy(1);
	gPad->SaveAs("new_vs_old_E_FieldZ_asymmetry.pdf");
	gPad->SaveAs("new_vs_old_E_FieldZ_asymmetry.png");
	// selecting events w/ asymmetry > 1%
	int nAbove1pc = t->Draw("asymm","asymm>0.01");
	gPad->SaveAs("new_vs_old_E_FieldZ_asymmetry_gt_1pc.pdf");
	gPad->SaveAs("new_vs_old_E_FieldZ_asymmetry_gt_1pc.png");
	double frac = (nAbove1pc*1.0)/(nevts*1.0);

	std::cout << "Fration of points w/ asymmetry > 1%: " << frac << "\n"; 
	TFile *outputFile = new TFile("new_vs_old_E_FieldZ_asymmetry.root","RECREATE");
	outputFile->cd();
	t->Write();
	outputFile->Close();
}



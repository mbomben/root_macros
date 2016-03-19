{
   //
   // This macro displays the Tree data structures
   //Author: Rene Brun

   gROOT->Reset();
  // create a TTree   
  TTree *tree = new TTree("tree","Electric field for Atlas pixels");

  Double_t x,y,z,E;

  // create a branch with energy 
  tree->Branch("x",&x);
  tree->Branch("y",&y);
  tree->Branch("z",&z);
  tree->Branch("E",&E);

  std::ifstream file;
  file.open("absEField-3D-map-fei4-200um-fl0-80V-noedge.txt");
  while( 1 ) {
    file >> x >> y >> z >> E;
    if ( file.eof() ) break;
    tree->Fill();
  } 

  TFile *f = new TFile("absEField-3D-map-fl6e15-600V-noedge.root","RECREATE");
  f->cd();
  tree->Write();
  f->Close();
}


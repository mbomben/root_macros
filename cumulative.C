void cumulative(TGraph* &gr) {
  int n = gr->GetN();
  double *X = new double[n];
  double *intY = new double[n];

  X = gr->GetX();

  for (int i = 1; i<n; i++ ) {
    intY[i] = gr->Integral(0,i);
    std::cout << intY[i] << "\n";
  }
  intY[n-1]= gr->Integral(0,-1);
  TGraph *gIntegral = new TGraph(n,X,intY);
  gIntegral->Draw("AC*");
}

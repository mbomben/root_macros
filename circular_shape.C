void circular_shape(double xc, double yc, double R, double theta0, double theta1, int N) {

  if ( N <=1 ) {
    std::cout << "N must be greater than 1\n";
    exit(2);
  }

  theta0 *= TMath::DegToRad();
  theta1 *= TMath::DegToRad();
  double delta_theta = theta1-theta0;
  std::cout << delta_theta/TMath::DegToRad() << "\n";

  double dtheta = delta_theta/double(N-1.0);

  int Npts = N;
  double *X = new double[Npts];
  double *Y = new double[Npts];

  double theta = 0;
  for (int i = 0; i < Npts; i++ ) {
    theta = theta0+(double)(i)*dtheta;
    X[i] = xc+R*TMath::Cos(theta);
    Y[i] = yc+R*TMath::Sin(theta);
  }

  TCanvas *c1 = new TCanvas("c1","",1000,1000);
  c1->cd();

  TGraph *gr = new TGraph(Npts,X,Y);
  gr->SetName("");
  gr->SetTitle(";X [#mum];Y [#mum]");
  gr->SetLineColor(kRed);
  gr->SetMarkerSize(1.0);
  gr->SetMarkerStyle(24);
  gr->Draw("ALP");
  gPad->SetGrid();
  gPad->SetTicks(1,1);

  std::cout << "\n\tpolygon=\"";
  std::cout << xc   << "," << yc   << " ";
  for (int i = 0; i < Npts-1; i++ ) {
    std::cout << X[i] << "," << Y[i] << " ";
  }
  std::cout << X[Npts-1] << "," << Y[Npts-1];
  std::cout << "\"\n\n";
}


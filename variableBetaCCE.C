
const double vsate = 1.06e7;  // cm/s
const double vsath = 1.06e7;  // cm/s


double variableBetaCCE( double phi = 1e14, double w = 200, double betae = 3.5e-7, double betah = 6.5e-7 ) {

  w *= 1e-4; // um -> cm conversion
  double de = vsate/phi/betae; // average mean free path for electrons
  double dh = vsath/phi/betah; // average mean free path for holes

  double lin = (de+dh)/w; // linear term

  double quade = (de/w)*(de/w); // quadratic factor for electrons
  double quadh = (dh/w)*(dh/w); // quadratic factor for holes

  double expe = (1.0-TMath::Exp(-w/de)); // exponential factor for electrons
  double exph = (1.0-TMath::Exp(-w/dh)); // exponential factor for holes

  double cce = lin - quade*expe - quadh*exph;
  
  std::cout << "For Phi = " << phi << "neq/cm2 in a "  << w*1e4<< " um thick detector the expected cce is " << cce << "\n";
  return cce;
}


double CV_correction( double f, double t, const char* tUnits = "C" ) {
  if ( std::strcmp( tUnits, "K" ) != 0 && std::strcmp( tUnits, "C" ) ) {
    usage();
    exit(1);
  }
  double coeff = 1.;

  double T = t;
  if ( std::strcmp( tUnits, "C" ) == 0 ) T += 273.15; // ºC to K

  double tFactor = 1.;
  double A = 0.0086;
  double delta = -0.114;
  double Ea = 0.635;


  tFactor = tempFactor(A,Ea,T);

  std::cout << "tFactor = " << tFactor << "\n";

  double fFactor = 1;

  fFactor = freqFactor(delta,f);

  std::cout << "fFactor = " << fFactor << "\n";

  coeff *= tFactor;
  coeff *= fFactor;
  std::cout << "coeff = " << coeff << "\n";

  double delta_error = 0.024;
  double A_error = 0.0016;
  double Ea_error = 0.050;

  double T_error = 1;
  double f_error = 50;

  double correction=systematics( A, Ea, T, delta, f, delta_error, A_error, Ea_error, T_error, f_error );
  std::cout << "Systematics = " << correction << "\n";
  return coeff;
}

void usage() {
  std::cout << "Usage: CV_correction( <frequence (Hz)>, <temperature>, [temperature units (C is default, K the other possible value)]\n";
}

double tempFactor(double A, double Ea, double T) {
  return (1+A*TMath::Exp(Ea/0.345))/(1+A*TMath::Exp((Ea/0.02354)*(((T-273.15)/T))));
}

double freqFactor( double delta, double f ) {
  return (1+delta)/(1+delta*TMath::Log10(f/1e3));
}

double systematics( double A, double Ea, double T, double delta, double f, double delta_error, double A_error, double Ea_error, double T_error, double f_error) {

  // central value
  double coeff = freqFactor(delta,f)*tempFactor(A,Ea,T);
  // syst. due to A
  double Aplus = A+A_error;
  double Aplus_coeff = freqFactor(delta,f)*tempFactor(Aplus,Ea,T);
  double Aminus = A-A_error;
  double Aminus_coeff = freqFactor(delta,f)*tempFactor(Aminus,Ea,T);

  double coeff_A_error = TMath::Max( TMath::Abs(Aplus_coeff-coeff), TMath::Abs(Aminus_coeff-coeff)  );
  // syst. due to Ea
  double Eaplus = Ea+Ea_error;
  double Eaplus_coeff = freqFactor(delta,f)*tempFactor(A,Eaplus,T);
  double Eaminus = Ea-Ea_error;
  double Eaminus_coeff = freqFactor(delta,f)*tempFactor(A,Eaminus,T);

  double coeff_Ea_error = TMath::Max( TMath::Abs(Eaplus_coeff-coeff), TMath::Abs(Eaminus_coeff-coeff)  );
  // syst. due to T
  double Tplus = T+T_error;
  double Tplus_coeff = freqFactor(delta,f)*tempFactor(A,Ea,Tplus);
  double Tminus = T-T_error;
  double Tminus_coeff = freqFactor(delta,f)*tempFactor(A,Ea,Tminus);

  double coeff_T_error = TMath::Max( TMath::Abs(Tplus_coeff-coeff), TMath::Abs(Tminus_coeff-coeff)  );
  // syst. due to f
  double fplus = f+f_error;
  double fplus_coeff = freqFactor(delta,fplus)*tempFactor(A,Ea,T);
  double fminus = f-f_error;
  double fminus_coeff = freqFactor(delta,fminus)*tempFactor(A,Ea,T);

  double coeff_f_error = TMath::Max( TMath::Abs(fplus_coeff-coeff), TMath::Abs(fminus_coeff-coeff)  );
  // syst. due to delta
  double deltaplus = delta+delta_error;
  double deltaplus_coeff = freqFactor(deltaplus,f)*tempFactor(A,Ea,T);
  double deltaminus = delta-delta_error;
  double deltaminus_coeff = freqFactor(deltaminus,f)*tempFactor(A,Ea,T);

  double coeff_delta_error = TMath::Max( TMath::Abs(deltaplus_coeff-coeff), TMath::Abs(deltaminus_coeff-coeff)  );

  double correction = TMath::Sqrt( coeff_delta_error*coeff_delta_error + coeff_f_error*coeff_f_error + coeff_T_error*coeff_T_error + 
                                   coeff_Ea_error*coeff_Ea_error + coeff_A_error*coeff_A_error );

  return correction;
}

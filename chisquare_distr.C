#include <iostream>
#include "TMath.h"


double chisquare_distr( double *x, double *par ) {


  double chi2 = x[0];

  if ( chi2 < 0 ) { 
    std::cout << "chi2 must be non-negative.\nExiting...\n";
    return -1;
  }

  double ndf = par[0];

  if ( ndf < 0 ) { 
    std::cout << "ndf must be non-negative.\nExiting...\n";
    return -2;
  }

  double nevts = par[1];

  if ( !(nevts > 0) ) { 
    std::cout << "number of events must be positive.\nExiting...\n";
    return -3;
  }

  double val = 0.0;

  val += TMath::Power(2.0,-ndf/2.0);
  val /= TMath::Gamma(ndf/2.0);
  val *= TMath::Power(chi2,(ndf-2.0)/2.0);
  val *= TMath::Exp(-chi2/2);
  val *= nevts;

  return val;

}

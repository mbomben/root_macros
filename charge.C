#include "TGraph.h"
#include <iostream>

double charge(const char* fileName, double t0 = 0.0, double t1 = -1.0) {
  if ( t0 < 0.0 ) {
    std::cout << "t0 must be >0!\nExiting...\n"; 
    exit(1);
  }
  if ( t1 < t0 && t1 > 0 ) {
    std::cout << "t1 must be > t0!\nExiting...\n"; 
    exit(2);
  }
  TGraph* g = new TGraph(fileName);
  int N = g->GetN();
  if (! N ) { 
    std::cout << "No entries!\nExiting...\n"; 
    exit(3);
  }
  double* X = new double[N];
  X = g->GetX();

  int bin0 = 0; 
  int bin1 = -1;

  if ( t0>0 ) {
    for ( int i = 0; i<N; i++ ) {
      if ( X[i] > t0 ) { bin0 = i-1; break; }
    }
  }

  //std::cout << "bin0 = " << bin0 << "\n";
  if ( t1>0 ) {
    for ( int i = 0; i<N; i++ ) {
      if ( X[i] > t1 ) bin1 = i-1;
    }
  }
  //std::cout << "bin1 = " << bin1 << "\n";
  double Q;
  Q = g->Integral(bin0);
  if ( bin1 >0 && bin0 >0) Q = g->Integral(bin0,bin1);
  if ( bin0 >0) Q = g->Integral(bin0,N-1);
  return Q;
}

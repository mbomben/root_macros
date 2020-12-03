 /*****************************************************************************
   * Project: RooFit                                                           *
   *                                                                           *
   * Binomial PDF                                                              *
   * author: Marco Bomben <marco.bomben@cern.ch>                               *
   * Inspired by RooPoisson by Kyle Cranmer <cranmer@cern.ch>                  * 
   *                                                                           *
   *****************************************************************************/
 
 /** \class RooBinomial
     \ingroup Roofit
 
 Binomial pdf
 **/
 
 #include "RooBinomial.h"
 
 #include "RooAbsReal.h"
 #include "RooAbsCategory.h"
 
 #include "RooRandom.h"
 #include "RooMath.h"
 #include "TMath.h"
 #include "Math/ProbFuncMathCore.h"
 #include "Math/PdfFuncMathCore.h"
 
 #include "BatchHelpers.h"
 #include "RooVDTHeaders.h"
 
 #include <limits>
 #include <cmath>
 
 using namespace std;
 
 ClassImp(RooBinomial);
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Constructor
 
 RooBinomial::RooBinomial(const char *name, const char *title,
              RooAbsReal& _k,
              RooAbsReal& _N,
              RooAbsReal& _p,
              Bool_t noRounding) :
   RooAbsPdf(name,title),
   k("k","k",this,_k),
   N("N","N",this,_N),
   p("p","p",this,_p),
   _noRounding(noRounding),
   _protectNegative(false)
 {
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Copy constructor
 
  RooBinomial::RooBinomial(const RooBinomial& other, const char* name) :
    RooAbsPdf(other,name),
    k("k",this,other.k),
    N("N",this,other.N),
    p("p",this,other.p),
    _noRounding(other._noRounding),
    _protectNegative(other._protectNegative)
 {
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Implementation in terms of the TMath::Binomial() function.
 
 Double_t RooBinomial::evaluate() const
 {
   Int_t nk = _noRounding ? k : floor(k);
   if(_protectNegative && p<0)
     return 1e-3;
   return ROOT::Math::beta_pdf(p, nk+1, N-nk+1);
 }
 
 
 
 namespace {
 
 template<class Tk, class Tp, class TN>
 void compute(const size_t n, double* __restrict output, Tk k, Tp p, TN N,
     const bool protectNegative, const bool noRounding) {
 
   for (size_t i = 0; i < n; ++i) { //CHECK_VECTORISE
     const double k_i = noRounding ? k[i] : floor(k[i]);
     const double N_i = noRounding ? N[i] : floor(N[i]);
     const double Nk_i = noRounding ? (N[i]-k[i]) : floor(N[i]-k[i]);
     // Do here all the calculations for the binomial coefficient
     // The std::lgamma yields different values than in the scalar implementation.
     // Need to check which one is more accurate.
 //    output[i] = std::lgamma(x_i + 1.);
     output[i] = TMath::LnGamma(N_i + 1.) - TMath::LnGamma(k_i + 1.) - TMath::LnGamma(Nk_i + 1.);
   }
 
 
   for (size_t i = 0; i < n; ++i) { //CHECK_VECTORISE
     const double k_i = noRounding ? k[i] : floor(k[i]);
     const double N_i = noRounding ? N[i] : floor(N[i]);
     const double Nk_i = noRounding ? (N[i]-k[i]) : floor(N[i]-k[i]);
     const double logp = _rf_fast_log(p[i]);
     const double log1minusp = _rf_fast_log(1.-p[i]);
     const double logBinomial = k_i * logp + Nk_i * log1minusp + output[i];
     output[i] = _rf_fast_exp(logBinomial);
 
     // Cosmetics
     if (k_i < 0.)
       output[i] = 0.;
     if (protectNegative && p[i] < 0.) // CHECK IF NEEDED
       output[i] = 1.E-3;
   }
 }
 
 }
 
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Compute Binomial values in batches.
 RooSpan<double> RooBinomial::evaluateBatch(std::size_t begin, std::size_t batchSize) const {
   using namespace BatchHelpers;
   auto kData = k.getValBatch(begin, batchSize);
   auto pData = p.getValBatch(begin, batchSize);
   auto NData = N.getValBatch(begin, batchSize);
   const bool batchk = !kData.empty();
   const bool batchp = !pData.empty();
   const bool batchN = !NData.empty();
 
   if (!batchp && !batchk && !batchN) {
     return {};
   }
   batchSize = findSize({ kData, pData, NData });
   auto output = _batchData.makeWritableBatchUnInit(begin, batchSize);
 
   if (batchk && !batchp && !batchN) {
     compute(batchSize, output.data(), kData, BracketAdapter<double>(p), BracketAdapter<double>(N),   _protectNegative, _noRounding);
   }
   else if (!batchk && batchp && !batchN ) {
     compute(batchSize, output.data(), BracketAdapter<double>(k), pData, BracketAdapter<double>(N), _protectNegative, _noRounding);
   }
   else if (!batchk && !batchp && batchN ) {
     compute(batchSize, output.data(), BracketAdapter<double>(k), BracketAdapter<double>(p), NData, _protectNegative, _noRounding);
   }
   else if (batchk && batchp && batchN ) {
     compute(batchSize, output.data(), kData, pData, NData, _protectNegative, _noRounding);
   }
   return output;
 }
 
 
 ////////////////////////////////////////////////////////////////////////////////
 
 Int_t RooBinomial::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
 {
   if (matchArgs(allVars,analVars,k)) return 1 ;
   if (matchArgs(allVars, analVars, p)) return 2;
   // We do not want to vary the total number of trials N
   return 0 ;
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 
 Double_t RooBinomial::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   R__ASSERT(code == 1 || code == 2) ;
 
   if(_protectNegative && p<0)
     return exp(-2*p); // make it fall quickly
 
   if (code == 1) {
     // Implement integral over k as summation. Add special handling in case
     // range boundaries are not on integer values of k
     const double kmin = std::max(0., k.min(rangeName));
     const double kmax = k.max(rangeName);
 
     if (kmax < 0. || kmax < kmin) {
       return 0.;
     }
     //if (!k.hasMax() || RooNumber::isInfinite(xmax)) { // k should be always less than N
     //  //Integrating the full Binomial distribution here 
     //  return 1.;
     //}
 
     // The range as integers. ikmin is included, ikmax outside.
     const unsigned int ikmin = kmin;
     const unsigned int ikmax = std::min(kmax + 1.,
                                         (double)std::numeric_limits<unsigned int>::max());
     
     // Sum from 0 to just before the bin outside of the range.
     if (ikmin == 0) {
       return ROOT::Math::binomial_cdf(ikmax - 1, p, N);
     }
     else {
       return ROOT::Math::binomial_cdf(ikmax - 1, p, N) - ROOT::Math::binomial_cdf(ikmin - 1, p, N);
     }
 
 
   } else if(code == 2) {
 
     // the integral with respect to the probability p is the integral of an incomplete beta function
     Double_t p_min = p.min(rangeName);
     Double_t p_max = p.max(rangeName);
 
     Double_t ik;
     if(_noRounding) ik = k + 1;
     else ik = Int_t(TMath::Floor(k)) + 1.0; // negative ik does not need protection (beta returns 0.0)
 
     return ROOT::Math::inc_beta(p_max, ik+1, N-ik+1) - ROOT::Math::inc_beta(p_max, ik+1, N-ik+1);
   }
 
   return 0;
 
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Advertise internal generator in k
 
 Int_t RooBinomial::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
 {
   if (matchArgs(directVars,generateVars,k)) return 1 ;
   return 0 ;
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Implement internal generator using TRandom::Binomial
 
 void RooBinomial::generateEvent(Int_t code)
 {
   R__ASSERT(code==1) ;
   Double_t kgen ;
   while(1) {
     kgen = RooRandom::randomGenerator()->Binomial(N,p);
     if (kgen<=k.max() && kgen>=k.min()) {
       k = kgen ;
       break;
     }
   }
   return;
 }

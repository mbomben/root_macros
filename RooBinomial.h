 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Binomial PDF                                                              *
  * author: Marco Bomben <marco.bomben@cern.ch>                               *
  * Inspired by RooPoisson by Kyle Cranmer <cranmer@cern.ch>                  * 
  *                                                                           * 
  *****************************************************************************/ 

#ifndef ROOBINOMIAL
#define ROOBINOMIAL

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooTrace.h"
 
class RooBinomial : public RooAbsPdf {
public:
  RooBinomial() { _noRounding = kFALSE ;   } ;
  RooBinomial(const char *name, const char *title, RooAbsReal& _k, RooAbsReal& _N, RooAbsReal& _p,  Bool_t noRounding=kFALSE);
  RooBinomial(const RooBinomial& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const override { return new RooBinomial(*this,newname); }
  inline virtual ~RooBinomial() {  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const override;
  void generateEvent(Int_t code) override;
  
  void setNoRounding(bool flag = kTRUE) {_noRounding = flag;}
  void protectNegativeMean(bool flag = kTRUE) {_protectNegative = flag;}  // is this needed?

protected:

  RooRealProxy k ;
  RooRealProxy N ;
  RooRealProxy p ;
  Bool_t  _noRounding ;
  Bool_t  _protectNegative ;
  
  Double_t evaluate() const override;
  RooSpan<double> evaluateBatch(std::size_t begin, std::size_t batchSize) const override;

  ClassDefOverride(RooBinomial,3) // A Binomial PDF
};
 
#endif

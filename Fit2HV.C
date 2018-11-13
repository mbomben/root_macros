#include <vector>
#include "TF1.h" 
#include "TGraph.h"
#include "TMath.h"
#include "TGraphErrors.h"
using namespace std;

TF1 *f1, *f2, *f0;

double finter(double *x, double *par) 
{
    return TMath::Abs(f1->EvalPar(x,par) - f2->EvalPar(x,par));
}

double fint(double start,double end) {                                
       TF1 *fint = new TF1("fint",finter,start,end,0);
       double xint = fint->GetMinimumX();
       return xint; 
}

Float_t Fit2HV(TGraphErrors * gr, Float_t tolerance, Float_t & rms){

  if ( fabs(tolerance-0.02) > 1e-6 ) {
    std::cout << "Recommended value for tolerance is 0.02\n";
  }
  vector<float> flat;
  
  Float_t xTest[gr->GetN()];
  Float_t yTest[gr->GetN()];
  
  for(Int_t ip = gr->GetN()-1; ip >= 0; ip--){
    Double_t x,y;
    gr->GetPoint(ip,x,y);
    
    if(ip == gr->GetN()-1 || ip == gr->GetN()-2){
      flat.push_back(ip);
      xTest[flat.size()-1] = x;
      yTest[flat.size()-1] = y;
      
    }else{
      
      TGraph * grTest = new TGraph(flat.size(), xTest, yTest);
      TF1 *Pol1 = new TF1("Pol1","pol1",0.,1000.);
      grTest->Fit(Pol1,"Q");

      Int_t n = flat.size()-1;
      
      if(TMath::Abs((y-Pol1->Eval(x))/y) < tolerance){
	flat.push_back(ip);
	xTest[flat.size()-1] = x;
	yTest[flat.size()-1] = y;
      }else{
	break;
      }
    }
    
  }
  
  Int_t nFlat = flat.size();
  Int_t nRise = gr->GetN() - nFlat;
  
  Float_t xRise[nRise];
  Float_t yRise[nRise];
  Float_t xRiseErr[nRise];
  Float_t yRiseErr[nRise];
  Float_t xFlat[nFlat];
  Float_t yFlat[nFlat];
  Float_t xFlatErr[nFlat];
  Float_t yFlatErr[nFlat];  

  for(Int_t ip = 0; ip < flat.size(); ip++){
    Double_t x,y;
    gr->GetPoint(flat.at(ip),x,y);
    //    cout << "Add Flat " << x << " " << y << endl;
    xFlat[ip] = x;
    yFlat[ip] = y;
    xFlatErr[ip] = gr->GetErrorX(flat.at(ip));
    yFlatErr[ip] = gr->GetErrorY(flat.at(ip));
  }
  
  for(Int_t ip = 0; ip < gr->GetN(); ip++){
    Double_t x,y;
    gr->GetPoint(ip,x,y);
    if(ip < nRise){
      //      cout << "Add Rise " << x << " " << y << endl;
      xRise[ip] = x;
      yRise[ip] = y;
      xRiseErr[ip] = gr->GetErrorX(ip);
      yRiseErr[ip] = gr->GetErrorY(ip);      
    }
  }

  //force (0,0)
  nRise++;
  xRise[nRise-1] = 0.;
  yRise[nRise-1] = 0.;
  xRiseErr[nRise-1] = 1.0;
  yRiseErr[nRise-1] = 0.10;
  
  TGraphErrors* grRise = new TGraphErrors(nRise, xRise, yRise, xRiseErr, yRiseErr);
  TGraphErrors* grFlat = new TGraphErrors(nFlat, xFlat, yFlat, xFlatErr, yFlatErr);

  f1 = new TF1("f1","[0]+[1]*sqrt(x)",0.,1000.);
  f2 = new TF1("f2","[0]+[1]*x",0.,1000.);
  
  grRise->Fit(f1,"Q");
  grFlat->Fit(f2,"Q");

  Float_t f1p0 = f1->GetParameter(0);
  Float_t f1p1 = f1->GetParameter(1);
  Float_t f2p0 = f2->GetParameter(0);
  Float_t f2p1 = f2->GetParameter(1);  

  Float_t f1p0Err = f1->GetParError(0);
  Float_t f1p1Err = f1->GetParError(1);
  Float_t f2p0Err = f2->GetParError(0);
  Float_t f2p1Err = f2->GetParError(1);  
  
  Float_t hv = fint(15.,750.);

  vector<float> hvv;
  
  for(Int_t i=0; i< 8; i++){
    if(i==0){
      f1->SetParameter(0,f1p0+f1p0Err);
      hvv.push_back(fint(5,750));
      f1->SetParameter(0,f1p0);
    }else if(i==1){
      f1->SetParameter(0,f1p0-f1p0Err);
      hvv.push_back(fint(5,750));
      f1->SetParameter(0,f1p0);
    }else if(i==2){
      f1->SetParameter(1,f1p1-f1p1Err);
      hvv.push_back(fint(5,750));
      f1->SetParameter(1,f1p1);
    }else if(i==3){
      f1->SetParameter(1,f1p1+f1p1Err);
      hvv.push_back(fint(5,750));
      f1->SetParameter(1,f1p1);            
    }else if(i==4){
      f2->SetParameter(0,f2p0-f2p0Err);
      hvv.push_back(fint(5,750));
      f2->SetParameter(0,f2p0);
    }else if(i==5){
      f2->SetParameter(0,f2p0+f2p0Err);
      hvv.push_back(fint(5,750));
      f2->SetParameter(0,f2p0);
    }else if(i==6){
      f2->SetParameter(1,f2p1-f2p1Err);
      hvv.push_back(fint(5,750));
      f2->SetParameter(1,f2p1);
    }else if(i==7){
      f2->SetParameter(1,f2p1+f2p1Err);
      hvv.push_back(fint(5,750));
      f2->SetParameter(1,f2p1);
    }
  }

  rms = 0.;
  for(Int_t i=0; i< 8; i++){
    rms += pow((hv-hvv.at(i)),2);
  }
  rms=sqrt(rms/8.);

  return hv;
}

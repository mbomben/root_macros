void convolution_box_gauss()

{

  // S e t u p   c o m p o n e n t   p d f s 

  // ---------------------------------------

 

  // Construct observable

  RooRealVar t("t","t",-100,100) ;

 

  // Construct landau(t,ml,sl) ;

  RooRealVar ml("ml","mean landau",5.,-20,20) ;

  RooRealVar sl("sl","sigma landau",1,0.1,10) ;

  RooLandau landau("lx","lx",t,ml,sl) ;

 

  // Construct gauss(t,mg,sg)

  RooRealVar mg("mg","mg",0) ;

  RooRealVar sg("sg","sg",200,0.1,10) ;

  RooGaussian gauss("gauss","gauss",t,mg,sg) ;

 

  // Construct Stepfuntion

  RooRealVar leftEdge("leftEdge","left edge",-50.,-100,100) ;

  RooRealVar rightEdge("rightEdge","right edge",50.,-100,100) ;

  RooGenericPdf stepFuncPDF("stepFuncPDF", "stepFuncPDF", "(@0 >= @1) && (@0 < @2)", RooArgList(t, leftEdge, rightEdge));

 

  // Construct gaussleft(t,mg,sg)

  RooRealVar mgl("mgl","mgl",-50, -20, 20) ;

  RooRealVar sgl("sgl","sgl",10, 0.1, 10) ;

  RooGaussian gaussl("gaussl","gauss",t,mgl,sgl) ;

 

  // Construct gaussright(t,mg,sg)

  RooRealVar mgr("mgr","mgr",50.,-20,20) ;

  RooRealVar sgr("sgr","sgr",10,0.1,10) ;

  RooGaussian gaussr("gaussr","gauss",t,mgr,sgr) ;

 

  // C o n s t r u c t   c o n v o l u t i o n   p d f 

  // ---------------------------------------

 

  // Set #bins to be used for FFT sampling to 10000

  t.setBins(10000,"cache") ; 

 

  // Construct step (x) gauss

  RooFFTConvPdf bxg("bxg","box (X) gauss",t,stepFuncPDF,gaussr) ;

 

  RooFFTConvPdf bxglr("bxglr","box (X) gausss",t,bxg,gaussl) ;

 

  RooDataSet* data = bxglr.generate(t,10000) ;

 

  // Visualize data set

//  RooPlot* frame = t.frame(Title("landau (x) gauss convolution")) ;

//  data->plotOn(frame) ;

//  frame->Draw() ;

 

  // Fit gxlx to daaa

  bxglr.fitTo(*data) ;

  // Plot data, landau pdf, landau (X) gauss pdf

  RooPlot* frame = t.frame(RooFit::Title("box (x) gauss convolution")) ;

  data->plotOn(frame) ;

  bxglr.plotOn(frame) ;

  stepFuncPDF.plotOn(frame,RooFit::LineStyle(kDashed)) ;

 

  // Draw frame on canvas

  new TCanvas("box_gauss_convolution","box_gauss_convolution",600,600) ;

  gPad->SetLeftMargin(0.15) ;

  frame->GetYaxis()->SetTitleOffset(1.4) ;

  frame->Draw() ;

 

}

using namespace RooStats;
using namespace RooFit;

void SimpleHypoTestInv( const char* infile =  "CountingModel.root", 
                        const char* workspaceName = "w",
                        const char* modelConfigName = "ModelConfig",
                        const char* dataName = "data" )
{
  /////////////////////////////////////////////////////////////
  // First part is just to access the workspace file 
  ////////////////////////////////////////////////////////////

  // open input file 
  TFile *file = TFile::Open(infile);
  if (!file) return;

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);


  // get the modelConfig out of the file
  RooStats::ModelConfig* mc = (RooStats::ModelConfig*) w->obj(modelConfigName);

  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);

  ModelConfig*  sbModel = (RooStats::ModelConfig*) w->obj(modelConfigName);
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
  bModel->SetName(TString(sbModel->GetName())+TString("_with_poi_0"));      
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );

  FrequentistCalculator  fc(*data, *bModel, *sbModel);
  //fc.SetToys(1000,500);  
  fc.SetToys(1000,500);   

  // asymptotic calculator
  AsymptoticCalculator  ac(*data, *bModel, *sbModel);
  ac.SetOneSided(true);  // for one-side tests (limits)
  //  ac->SetQTilde(true);
  AsymptoticCalculator::SetPrintLevel(-1);


  // create hypotest inverter 
  // passing the desired calculator 
  HypoTestInverter calc(ac);    // for asymptotic 
  //HypoTestInverter calc(fc);  // for frequentist

  // set confidence level (e.g. 95% upper limits)
  calc.SetConfidenceLevel(0.95);
  
  // for CLS
  bool useCLs = true;
  calc.UseCLs(useCLs);
  calc.SetVerbose(false);

  // configure ToyMC Samler (needed only for frequentit calculator)
  ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
   
  // profile likelihood test statistics 
  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  // for CLs (bounded intervals) use one-sided profile likelihood
  if (useCLs) profll.SetOneSided(true);

  // ratio of profile likelihood - need to pass snapshot for the alt 
  // RatioOfProfiledLikelihoodsTestStat ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   
  // set the test statistic to use 
  toymcs->SetTestStatistic(&profll);

  // if the pdf is not extended (e.g. in the Poisson model) 
  // we need to set the number of events
  if (!sbModel->GetPdf()->canBeExtended())
     toymcs->SetNEventsPerToy(1);


  int npoints = 40;  // number of points to scan
  // min and max (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  //poimin = 0; poimax=10;

  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;   calc.SetFixedScan(npoints,poimin,poimax);      HypoTestInverterResult * r = calc.GetInterval();   double upperLimit = r->UpperLimit();

  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
      // compute expected limit
   std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
   std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
      // plot now the result of the scan    
   HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",r);
   // plot in a new canvas with style
   TCanvas * c1 = new TCanvas("HypoTestInverter Scan");    c1->SetLogy(false);

  plot->Draw("CLb 2CL");  // plot also CLb and CLs+b 
  //plot->Draw("OBS");  // plot only observed p-value


  // plot also in a new canvas the test statistics distributions 
  
  // plot test statistics distributions for the two hypothesis
  // when distribution is generated (case of FrequentistCalculators)
  const int n = r->ArraySize();
  if (n> 0 &&  r->GetResult(0)->GetNullDistribution() ) { 
     TCanvas * c2 = new TCanvas("Test Statistic Distributions","",2);
     if (n > 1) {
        int ny = TMath::CeilNint( sqrt(n) );
        int nx = TMath::CeilNint(double(n)/ny);
        c2->Divide( nx,ny);
     }
     for (int i=0; i<n; i++) {         if (n > 1) c2->cd(i+1);
        SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
     }
  }


}

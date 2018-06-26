using namespace RooFit;
using namespace RooStats;

void CountingModel(  int nobs = 2,           // number of observed events
                     double b = 2.234,           // number of background events
                     double sigmab = 0.669 )   // relative uncertainty in b
{
   RooWorkspace w("w");
   
// make Poisson model * Gaussian constraint
   w.factory("sum:nexp(s[1,0,15],b[2.234,0,12])");
// Poisson of (n | s+b)
   w.factory("Poisson:pdf(nobs[0,15],nexp)");
   w.factory("Gaussian:constraint(b0[0,12],b,sigmab[1])");
   w.factory("PROD:model(pdf,constraint)");


   w.var("b0")->setVal(b);
   w.var("b0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmab")->setVal(sigmab*b);  
   

   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("s"));
   mc.SetObservables(*w.var("nobs"));
   mc.SetNuisanceParameters(*w.var("b"));

   // these are needed for the hypothesis tests
   mc.SetSnapshot(*w.var("s"));
   mc.SetGlobalObservables(*w.var("b0"));

   mc.Print();
   // import model in the workspace 
   w.import(mc);

   // make data set with the namber of observed events
   RooDataSet data("data","", *w.var("nobs"));
   w.var("nobs")->setVal(nobs);
   data.add(*w.var("nobs") );
   // import data set in workspace and save it in a file
   w.import(data);

   w.Print();

   TString fileName = "CountingModel.root"; 

   // write workspace in the file (recreate file if already existing)
   w.writeToFile(fileName, true);

}

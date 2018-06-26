#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <TCanvas.h>

//header to run Aplanarity
#include "Sphericity.cxx"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

std::map<TString, TH1F*> h1f_map;

bool passSR(Long64_t *cutflow,
	    //input variable
	    vector<double>& Jet20PT, vector<double>& Jet20ETA, vector<double>& Jet20PHI,
	    double mET, double dPhiMin20, double dPhiMin50,
	    double EtmissOverHT, double Aplanarity, double MeffIncl, 
	    //selections
	    TString SRname, int NjetsCut,
            double *jetPTCuts, double jetEtaCuts,
            double dPhiMin20Cut, double dPhiMin50Cut,
            double EtmissOverMeffCut, double EtmissOverHTCut,
            double AplanarityCut, double MeffInclCut);


TFile *fglobe=new TFile("Aplanar.root","RECREATE");
//----------------SETUP HISTOGRAM----------------------                                                                                     
TH1D *aplanarity4j1000 =  new TH1D("aplanarity4j1000","Aplanarity",24,0,0.5);
TH1D *aplanarity4j2600 =  new TH1D("aplanarity4j2600","Aplanarity",24,0,0.5);
TH1D *aplanarity4j3000 =  new TH1D("aplanarity4j3000","Aplanarity",24,0,0.5);

int jets_13analysis_ver2_Aplanarity(){
  //----------------READ INPUT FILE----------------------                              
  TFile *fin=new TFile("/phys/groups/tev/scratch4/users/paklim/\
LHT_Data_Backup/save_scan_jets_met_13TeV/\
f1000/jmetk2/tag_1_delphes_events.root","READONLY");
  TFile *fout=new TFile("signal3.root","RECREATE");
  //----------------SETUP HISTOGRAM----------------------                              
  TH1D *info = new TH1D("info","info",14,0,14);
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",50,0,2000);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",50,0,500);
  TH1D *h_m_average = new TH1D("h_m_average","m_average",200,0,200);
  TH1D *h_MET = new TH1D("h_MET","h_MET",50,0,2000);
  TH1D *h_dRmax = new TH1D("h_dRmax","dRmax",24,0,4);
  //----------------END OF HISTOGRAM---------------------                              
  TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  cout << nevent << endl;
  t->SetMakeClass(1);
  //----------------SETUP BRANCHES----------------------- 
  //MET                                 
  Int_t nMET = -1;
  Float_t MET[50];
  Float_t MET_PHI[50];
  t->SetBranchAddress("MissingET", &nMET);
  TBranch *met = t->GetBranch("MissingET.MET");
  met->SetAddress(&MET);
  TBranch *met_phi = t->GetBranch("MissingET.Phi");
  met_phi->SetAddress(&MET_PHI);
  //Jet               
  Int_t nJET = -1;
  Float_t JET_PT[100];
  Float_t JET_ETA[100];
  Float_t JET_PHI[100];
  Float_t JET_M[100];
  UInt_t JET_BTAG[100];
  //TLorentzVector all_vector[40];
  t->SetBranchAddress("Jet", &nJET);
  TBranch *jet_pt = t->GetBranch("Jet.PT");
  jet_pt->SetAddress(&JET_PT);
  TBranch *jet_eta = t->GetBranch("Jet.Eta");
  jet_eta->SetAddress(&JET_ETA);
  TBranch *jet_phi = t->GetBranch("Jet.Phi");
  jet_phi->SetAddress(&JET_PHI);
  TBranch *jet_m = t->GetBranch("Jet.Mass");
  jet_m->SetAddress(&JET_M);
  TBranch *jet_btag = t->GetBranch("Jet.BTag");
  jet_btag->SetAddress(&JET_BTAG);
  //Electron        
  Int_t nElec = -1;
  Float_t ELEC_PT[50];
  Float_t ELEC_ETA[50];
  Float_t ELEC_PHI[50];
  Float_t ELEC_SUMPT[50];
  Float_t ELEC_HOE[50];
  t->SetBranchAddress("Electron", &nElec);
  TBranch *e_pt = t->GetBranch("Electron.PT");
  e_pt->SetAddress(&ELEC_PT);
  TBranch *e_eta = t->GetBranch("Electron.Eta");
  e_eta->SetAddress(&ELEC_ETA);
  TBranch *e_phi = t->GetBranch("Electron.Phi");
  e_phi->SetAddress(&ELEC_PHI);
  TBranch *e_sumpt = t->GetBranch("Electron.SumPtCharged");
  e_sumpt->SetAddress(&ELEC_SUMPT);
  TBranch *e_hoe = t->GetBranch("Electron.EhadOverEem");
  e_hoe->SetAddress(&ELEC_HOE);
  //Muon                   
  Int_t nMUON = -1;
  Float_t MUON_PT[50];
  Float_t MUON_ETA[50];
  Float_t MUON_PHI[50];
  Float_t MUON_SUMPT[50];
  t->SetBranchAddress("Muon", &nMUON);
  TBranch *muon_pt = t->GetBranch("Muon.PT");
  muon_pt->SetAddress(&MUON_PT);
  TBranch *muon_eta = t->GetBranch("Muon.Eta");
  muon_eta->SetAddress(&MUON_ETA);
  TBranch *muon_phi = t->GetBranch("Muon.Phi");
  muon_phi->SetAddress(&MUON_PHI);
  TBranch *muon_sumpt = t->GetBranch("Muon.SumPtCharged");
  muon_sumpt->SetAddress(&MUON_SUMPT);

  /////////////////////////////////   
  // Start Cut Analysis          //       
  /////////////////////////////////
  TString signalRegions[11]={"3j-1300","4j-1000","4j-1400","4j-1800","4j-2200","4j-2600","4j-3000",
			     "5j-1700","5j-1600","5j-2000","5j-2600"};
  TString cutName[10]={"pTA", "pTB", "pTC", "eta", "dPhiMin",
		       "dPhiMin(MET,j50)","EtMiss/Meff", "EtMiss/HT","Aplanarity", "MeffIncl"};    

  Long64_t cut_flow_SR[24][15];
  for(int i=0; i< 24; i++)
    for(int j=0; j< 15; j++)
      cut_flow_SR[i][j]=0;

  //SR Cuts Defitions
  const int Nchannel=11;
  int NJetsCut[Nchannel];
  double jetPTCuts[Nchannel][10];
  double jetEtaCuts[Nchannel];
  double dPhiMin20Cut[Nchannel], dPhiMin50Cut[Nchannel];
  double EtmissOverHTCut[Nchannel], EtmissOverMeffCut[Nchannel];
  double AplanarityCut[Nchannel], MeffInclCut[Nchannel];

  for(int ch=0;ch<11;ch++){
    for(int j=0;j<10;j++){
      jetPTCuts[ch][j]=-1;
      jetEtaCuts[ch]=-1;
      dPhiMin20Cut[ch]=-1;
      dPhiMin50Cut[ch]=-1;
      EtmissOverHTCut[ch]=-1;
      EtmissOverMeffCut[ch]=-1;
      AplanarityCut[ch]=-1;
      MeffInclCut[ch]=-1;
    }
  }
  //2j-1200      
  //ch0-6 
 
  int ch=0;
  signalRegions[ch]="3j-1300";
  NJetsCut[ch]=3;
  jetPTCuts [ch][0]=700;   jetPTCuts[ch][1]=50;   jetPTCuts[ch][2]=50;  
  jetEtaCuts[ch]=2.8;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.2;
  EtmissOverHTCut[ch]=16; EtmissOverMeffCut[ch]=-1;
  AplanarityCut[ch]=-1; MeffInclCut[ch]=1300;  

  ch=1;
  signalRegions[ch]="4j-1000";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;    jetPTCuts[ch][3]=100;
  jetEtaCuts[ch]=1.2;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.3;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=1000;
  
  ch=2;
  signalRegions[ch]="4j-1400";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;    jetPTCuts[ch][3]=100;
  jetEtaCuts[ch]=2.0;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.25;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=1400;

  ch=3;
  signalRegions[ch]="4j-1800";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;    jetPTCuts[ch][3]=100;
  jetEtaCuts[ch]=2.0;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.25;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=1800;

  ch=4;
  signalRegions[ch]="4j-2200";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;    jetPTCuts[ch][3]=100;
  jetEtaCuts[ch]=2.0;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.25;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=2200;

  ch=5;
  signalRegions[ch]="4j-2600";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;    jetPTCuts[ch][3]=150;
  jetEtaCuts[ch]=2.0;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.2;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=2600;

  ch=6;
  signalRegions[ch]="4j-3000";
  NJetsCut[ch]=4;
  jetPTCuts [ch][0]=200;   jetPTCuts[ch][3]=150;
  jetEtaCuts[ch]=2.0;
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.2;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.2;
  AplanarityCut[ch]=0.04; MeffInclCut[ch]=3000;

  ch=7;
  signalRegions[ch]="5j-1700";
  NJetsCut[ch]=5;
  jetPTCuts [ch][0]=700;   jetPTCuts[ch][3]=50;     jetPTCuts[ch][4]=50;
  jetEtaCuts[ch]=2.8;  //2.8 instead of 5
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.2;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.3;
  AplanarityCut[ch]=-1; MeffInclCut[ch]=1700;

  ch=8;
  signalRegions[ch]="5j-1600";
  NJetsCut[ch]=5;
  jetPTCuts [ch][0]=200;   jetPTCuts[ch][5]=50;    
  jetEtaCuts[ch]=2.8;  //2.8 instead of 5                                                       
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.2;
  EtmissOverHTCut[ch]=-1; EtmissOverMeffCut[ch]=0.15;
  AplanarityCut[ch]=0.08; MeffInclCut[ch]=1600;

  ch=9;
  signalRegions[ch]="5j-2000";
  NJetsCut[ch]=5;
  jetPTCuts [ch][0]=200;   jetPTCuts[ch][5]=50;
  jetEtaCuts[ch]=2.8;  //2.8 instead of 5                                                       
  dPhiMin20Cut[ch]=0.4; dPhiMin50Cut[ch]=0.2;
  EtmissOverHTCut[ch]=15; EtmissOverMeffCut[ch]=-1;
  AplanarityCut[ch]=-1; MeffInclCut[ch]=2000;

  ch=10;
  signalRegions[ch]="5j-2600";
  NJetsCut[ch]=5;
  jetPTCuts [ch][0]=200;   jetPTCuts[ch][5]=50;
  jetEtaCuts[ch]=2.8;  //2.8 instead of 5                                                       
  dPhiMin20Cut[ch]=0.8; dPhiMin50Cut[ch]=0.4;
  EtmissOverHTCut[ch]=18; EtmissOverMeffCut[ch]=-1;
  AplanarityCut[ch]=-1; MeffInclCut[ch]=2600;



  Long64_t cut_flow[3]; std::fill_n(cut_flow, 3, 0);
  Int_t nevents = nevent;  //default set to nevent       
  for(Long64_t ievent = 0; ievent < nevents; ievent++){
    cut_flow[0]++;
    met->GetEntry(ievent);
    met_phi->GetEntry(ievent);
    jet_pt->GetEntry(ievent);
    jet_phi->GetEntry(ievent);
    jet_eta->GetEntry(ievent);
    jet_m->GetEntry(ievent);
    jet_btag->GetEntry(ievent);
    e_pt->GetEntry(ievent);
    e_eta->GetEntry(ievent);
    e_phi->GetEntry(ievent);
    e_sumpt->GetEntry(ievent);
    e_hoe->GetEntry(ievent);
    muon_pt->GetEntry(ievent);
    muon_eta->GetEntry(ievent);
    muon_phi->GetEntry(ievent);
    muon_sumpt->GetEntry(ievent);

    /////BaseCut/////    
    // letpon veto                                                               
    Int_t num_Lep = 0;
    for(int j=0; j<nMUON; j++){
      if(j>=10) break;
      if(MUON_PT[j] > 10 && fabs(MUON_ETA[j]) < 2.4) num_Lep++;
    }
    for(int j=0; j<nElec; j++){
      if(j>=10) break;
      if(ELEC_PT[j] > 20 && fabs(ELEC_ETA[j]) < 2.47) num_Lep++;
    }
    if(num_Lep>0) continue;
    cut_flow[1]++;

    // MET
    if(nMET<1) continue;
    cut_flow[2]++;    
    
    // select baseline jets
    vector<double> Jet20PT, Jet20ETA, Jet20PHI;
    vector<double> Jet50PT, Jet50ETA, Jet50PHI;
    for(int j=0; j<nJET; j++){
      if(j>=40) break;
      if(JET_PT[j]<=20 || fabs(JET_ETA[j])>=2.8) continue;
      Jet20PT.push_back (JET_PT[j]);
      Jet20ETA.push_back(JET_ETA[j]);
      Jet20PHI.push_back(JET_PHI[j]);
      
      if(JET_PT[j]<=50 || fabs(JET_ETA[j]>=2.8))continue;
      Jet50PT.push_back(JET_PT[j]);
      Jet50ETA.push_back(JET_ETA[j]);
      Jet50PHI.push_back(JET_PHI[j]);
    }
   
    //calculate global variables common to all signal regions
    // mindPhiMet20
    // For the 2-jet SR and Meff-5j-2600, up to three leading jets (if present)
    double dPhiMin20 = TMath::TwoPi();
    for(unsigned int j=0; j< Jet20PHI.size(); j++){
      if(j>=3) break;
      float dPhi = fmod( fabs(MET_PHI[0] - Jet20PHI.at(j)), TMath::TwoPi() );
      if(dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
      if(dPhi <  dPhiMin20) dPhiMin20 = dPhi;
    }

    // minPhiMet50
    // SR with at least 4,5,6, or more than three jets are present in 2-jet or 3-jet SRs
    double dPhiMin50 = TMath::TwoPi();
    for(unsigned int j=0; j< Jet50PHI.size(); j++){
      float dPhi = fmod( fabs(MET_PHI[0] - Jet50PHI.at(j)), TMath::TwoPi() );
      if(dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
      if(dPhi <  dPhiMin50) dPhiMin50 = dPhi;
    }

    // HT
    double HT=0;
    for(unsigned int j=0; j< Jet20PT.size(); j++){
      HT += Jet20PT.at(j);
    }
    double EtmissOverHT = MET[0]/sqrt(HT);

    // inclusive Meff     
    double MeffIncl = MET[0];
    for(unsigned int j=0; j<Jet50PT.size(); j++){
      MeffIncl += Jet50PT.at(j);
    }

    ///// *Aplanarity*/////
    double Sp=0;
    double ST=0;
    double Ap=0;

    vector<TLorentzVector> v_tlv;
    //prepare vector<TLorentzVector> of jets to use  
    for(size_t ijet=0; ijet< nJET; ijet++)  {
      if ( JET_PT[ijet] < 50. || fabs(JET_ETA[ijet])>2.8) continue;
      TLorentzVector jet;
      jet.SetPtEtaPhiM(JET_PT[ijet],
		       JET_ETA[ijet],
		       JET_PHI[ijet],
		       JET_M[ijet]);
      v_tlv.push_back(jet);
    }
    
    int njet = v_tlv.size();
    if(njet>2){
      Sphericity sp; //construct     
      sp.SetTLV(v_tlv, njet);
      sp.GetSphericity(Sp, ST, Ap); 
    }
    double Aplanarity = Ap;
 
    //cout << "Aplanarity " << Aplanarity << "event " << ievent << endl;
   
    
    
    ///// Signal Region Counts /////
    for(int k=0;k<Nchannel;k++){
      passSR(&cut_flow_SR[k][0],
	     //input variable
	     Jet20PT, Jet20ETA, Jet20PHI, 
	     MET[0], dPhiMin20, dPhiMin50,
	     EtmissOverHT, Aplanarity, MeffIncl,
	     //selections
	     signalRegions[k], NJetsCut[k],
	     &jetPTCuts[k][0], jetEtaCuts[k],
	     dPhiMin20Cut[k], dPhiMin50Cut[k],
	     EtmissOverMeffCut[k], EtmissOverHTCut[k],
	     AplanarityCut[k], MeffInclCut[k]);
    }
    
    
  } //end of events                                                                  
  
  ///// Output Cut Info  /////                                                       
  cout << "Output cut info ..." << endl;
  cout << " " << endl;
  ofstream fileout("cutinfo_local.dat");
  fileout << "total events: "       << " " << cut_flow[0] << endl;
  fileout << "lepton veto: "        << " " << cut_flow[1] << endl;
  for(int i=0; i< Nchannel; i++){
    for(int j=0; j< 10; j++){
      fileout << signalRegions[i] << ": " << cutName[j] << ": " << cut_flow_SR[i][j] << endl;
    }
  }
  fileout.close();
  fileout.close();

  // output histogram          
  int ipic = 0;
  int nhist = 11;
  TCanvas *c[nhist];
  TFile* fhistout=new TFile("hist_output.root","recreate");
  for( const auto& sm_pair : h1f_map ){
    TH1F* h1f = (TH1F*) sm_pair.second;
    h1f->Write();

    ipic +=1;
    c[ipic] = new TCanvas(Form("c%d",ipic));
    h1f->Draw();
    h1f->GetYaxis()->SetTitle("Aplanarity");
    c[ipic]->SetLogy();
    c[ipic]->SaveAs();
    
  }
  fhistout->Close();
 

  return 1;
}


void  fillHistogram(TString hname, double nBin, double min, double max, 
		    double value, double weight){
  
  auto hist =   h1f_map.find(hname);
  TH1F* h1f = 0;
  if( hist != h1f_map.end()) h1f = (TH1F*) hist->second;
  else{
    h1f = new TH1F(hname.Data(),"",nBin, min, max);
    h1f_map[hname]=h1f;
  }
  h1f->Fill(value, weight);
 
}


///// Setup CutFucntion passSR /////
bool passSR(Long64_t *cutflow,
            //input variable                  
            vector<double>& Jet20PT, vector<double>& Jet20ETA, vector<double>& Jet20PHI,
            double mET, double dPhiMin20, double dPhiMin50,
            double EtmissOverHT, double Aplanarity, double MeffIncl,
	    //selections               
            TString SRname, int NjetsCut,
            double *jetPTCuts, double jetEtaCuts,
            double dPhiMin20Cut, double dPhiMin50Cut,
            double EtmissOverMeffCut, double EtmissOverHTCut,
            double AplanarityCut, double MeffInclCut) {
  
  int index=0;
  if((int)Jet20PT.size() < NjetsCut) return false;
  
  //// pt cut ////
  // pt(j1)
  if(Jet20PT[0] < jetPTCuts[0]) return false;
  cutflow[index]++; index++;
  
  // case of 2j SR
  if(NjetsCut==2){
    //pt(j2)
    if(Jet20PT[1] < jetPTCuts[1]) return false;
    cutflow[index]++; index++;
    //dummy counter 
    cutflow[index]++; index++;
  }
  // case of 3j SR
  else if(NjetsCut==3){
    //pt(j2)
    if(Jet20PT[1] < jetPTCuts[1]) return false;
    cutflow[index]++; index++;
    //pt(j3)
    if(Jet20PT[2] < jetPTCuts[2]) return false;
    cutflow[index]++; index++;
  }
  // case of 4j SR
  else if(NjetsCut==4){
    if(Jet20PT[3] < jetPTCuts[3]) return false;
    cutflow[index]++; index++;
    //dummy counter            
    cutflow[index]++; index++;
  }
  // case of 5j SR
  else if(NjetsCut==5 && SRname==TString("5j-1700")){
    // pt(j4)              
    if(Jet20PT[3] < jetPTCuts[3]) return false;
    cutflow[index]++; index++;
    // pt(j5)                              
    if(Jet20PT[4] < jetPTCuts[4]) return false;
    cutflow[index]++; index++;
  }
  // case of 6j SR
  else if(NjetsCut>=5 && Jet20PT.size()>5){
    // pt(j6)           
    if(Jet20PT[5] < jetPTCuts[5]) return false;
    cutflow[index]++; index++;
    //dummy counter         
    cutflow[index]++; index++;
  }//else{
  //  cout <<"Error NjetsCut "<< NjetsCut <<" not defined!"<< endl;
  //}
  
  // EtaCut
  if(NjetsCut<=3){
    //eta1,2   
    for(int j=0;j<2;j++){
      if(fabs(Jet20ETA[j]) > jetEtaCuts) return false;
    }
    cutflow[index]++; index++;
  }
  else{
    //eta1,2,3,4,5,6           
    for(int j=0;j<NjetsCut;j++){
      if(fabs(Jet20ETA[j]) > jetEtaCuts) return false;
    }
    cutflow[index]++; index++;
  }

  // dPhiMin cut                                                     
  if(dPhiMin20 < dPhiMin20Cut) return false;
  cutflow[index]++; index++;

  // dPhiMin50 cut         
  if(NjetsCut>3  && dPhiMin50 < dPhiMin50Cut) return false;
  cutflow[index]++; index++;

  // MET/meff cut      
  double MeffNj=0;
  for(int j=0;j<NjetsCut;j++)
    MeffNj+= Jet20PT.at(j);
  double EtmissOverMeff = mET/MeffNj;
  if(EtmissOverMeff < EtmissOverMeffCut) return false;
  cutflow[index]++; index++;

  // MET/sqrt(HT) cut        
  if(EtmissOverHT < EtmissOverHTCut) return false;
  cutflow[index]++; index++;
  
  //Aplanarity cut                                                                     
  if(Aplanarity < AplanarityCut    ) return false;
  cutflow[index]++; index++;

  fillHistogram(SRname+TString("_Aplanarity"),100,0,0.5,Aplanarity, 1.0);

  //Meff(inclusive) cut                                                                   
  if(MeffIncl < MeffInclCut) return false;
  cutflow[index]++; index++;  

  return true;  
}
 

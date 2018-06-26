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
#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;


int kinematicsHHWH(){

  //Read input
  //TFile *fin=new TFile("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/f700/jmetk1/tag_1_delphes_events.root","READONLY");
  TChain *t = new TChain("Delphes");

    t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/HHmet/Validation/13processHHSample_f1200_k05/f1200k5_DATA/dihiggs*.root");
  

    //TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  //nevent = 20000; 
  Double_t totevent = nevent*1.;
  cout << "total events"  <<nevent << endl;
  t->SetMakeClass(1);

  //Setup Branches
  const Int_t    MaxnMET = 5000;
  Int_t          nMET = -1;
  Float_t        MET[MaxnMET];
  Float_t        MET_PHI[MaxnMET];
  TBranch        *met;
  TBranch        *met_phi;
  t->SetBranchAddress("MissingET", &nMET);
  t->SetBranchAddress("MissingET.MET", MET, &met);
  t->SetBranchAddress("MissingET.Phi", MET_PHI, &met_phi);

  const Int_t    MaxnJet = 5000;
  Int_t          nJET = -1;
  Float_t        JET_PT[MaxnJet];
  Float_t        JET_ETA[MaxnJet];
  Float_t        JET_PHI[MaxnJet];
  Float_t        JET_M[MaxnJet];
  UInt_t         JET_BTAG[MaxnJet];
  Int_t          JET_FLAVOR[MaxnJet];
  TLorentzVector all_vector[MaxnJet];
  TBranch        *jet_pt;
  TBranch        *jet_eta;
  TBranch        *jet_phi;
  TBranch        *jet_m;
  TBranch        *jet_btag;
  TBranch        *jet_flavor;
  t->SetBranchAddress("Jet", &nJET);
  t->SetBranchAddress("Jet.PT", JET_PT, &jet_pt);
  t->SetBranchAddress("Jet.Eta", JET_ETA, &jet_eta);
  t->SetBranchAddress("Jet.Phi", JET_PHI, &jet_phi);
  t->SetBranchAddress("Jet.Mass", JET_M, &jet_m);
  t->SetBranchAddress("Jet.BTag", JET_BTAG, &jet_btag);
  t->SetBranchAddress("Jet.Flavor", JET_FLAVOR, &jet_flavor);

  const Int_t     kMaxParticle = 5000;
  Int_t           nParticle = -1;
  Int_t           Particle_PID[kMaxParticle];   //[Particle_]                              
  Int_t           Particle_Status[kMaxParticle];   //[Particle_]                           
  Float_t         Particle_PT[kMaxParticle];   //[Particle_]                               
  Float_t         Particle_Eta[kMaxParticle];   //[Particle_]                              
  Float_t         Particle_Phi[kMaxParticle];   //[Particle_]                              
  Float_t         Particle_Mass[kMaxParticle];   //[Particle_]                             
  Int_t           PARTICLE_M1[kMaxParticle];
  Int_t           PARTICLE_M2[kMaxParticle];
  TBranch        *particle_pid;
  TBranch        *particle_status;
  TBranch        *particle_pt; 
  TBranch        *particle_eta; 
  TBranch        *particle_phi;
  TBranch        *particle_mass;
  TBranch        *particle_m1;
  TBranch        *particle_m2;
  t->SetBranchAddress("Particle", &nParticle);
  t->SetBranchAddress("Particle.PID", Particle_PID, &particle_pid);
  t->SetBranchAddress("Particle.Status", Particle_Status, &particle_status);
  t->SetBranchAddress("Particle.PT", Particle_PT, &particle_pt);
  t->SetBranchAddress("Particle.Eta", Particle_Eta, &particle_eta);
  t->SetBranchAddress("Particle.Phi", Particle_Phi, &particle_phi);
  t->SetBranchAddress("Particle.Mass", Particle_Mass, &particle_mass);
  t->SetBranchAddress("Particle.M1", PARTICLE_M1, &particle_m1);
  t->SetBranchAddress("Particle.M2", PARTICLE_M2, &particle_m2);

  int NDivision = nevent/10;
  int wcount, hcount;
  int totalww=0, totalhh=0, totalwh=0, totalw=0, totalh=0, totalnon=0;
  int totalOther=0;
  for(Long64_t jevent = 0; jevent < nevent; jevent++){
    if(jevent%NDivision ==0) cout << jevent << "/" << nevent <<endl;

    Long64_t ievent = t->LoadTree(jevent);

    met->GetEntry(ievent);
    met_phi->GetEntry(ievent);
    jet_pt->GetEntry(ievent);
    jet_phi->GetEntry(ievent);
    jet_eta->GetEntry(ievent);
    jet_m->GetEntry(ievent);
    jet_btag->GetEntry(ievent);
    jet_flavor->GetEntry(ievent);
    particle_pid->GetEntry(ievent);
    particle_status->GetEntry(ievent);
    particle_pt->GetEntry(ievent);
    particle_eta->GetEntry(ievent);
    particle_phi->GetEntry(ievent);
    particle_mass->GetEntry(ievent);
    particle_m1->GetEntry(ievent);
    particle_m2->GetEntry(ievent);

    //////////////                                                                          
    // build table                                                                           
    // WW, HH, WH                                                                               
    //////////////                                                                           
    wcount=0, hcount=0;
    vector<int> WhIndex;
    vector<int> ZhIndex;
    vector<int> Windex;
    vector<int> hindex;


    // loop over address of Wh and Zh particles     
    for(int i=0; i<nParticle; i++){
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==8880024) WhIndex.push_back(i);
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==8880023) ZhIndex.push_back(i);
    }

    for(int i=0; i<nParticle; i++){
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==24 && WhIndex.size()!=0){
        for(int j=0; j<WhIndex.size(); j++){
          if(PARTICLE_M1[i]==WhIndex[j] || PARTICLE_M2[i]==WhIndex[j]) {
	    wcount++;
	    Windex.push_back(i);
	  }
        }
      }
      if(Particle_Status[i]==3 && Particle_PID[i]==25 && ZhIndex.size()!=0){
        for(int j=0; j<ZhIndex.size(); j++){
          if(PARTICLE_M1[i]==ZhIndex[j] || PARTICLE_M2[i]==ZhIndex[j]) {
	    hcount++;
	    hindex.push_back(i);
	  }
        }
      }
    }

    //HH
    Double_t H_lead_Pt, H_sublead_Pt;
    if(wcount==0 && hcount==2){
      H_lead_Pt    = Particle_PT[hindex[0]];
      H_sublead_Pt = Particle_PT[hindex[1]];
      if(Particle_PT[hindex[1]]>H_lead_Pt){
	H_lead_Pt    = Particle_PT[hindex[1]];
	H_sublead_Pt = Particle_PT[hindex[0]];
	fillHistogram("leading higgs pT",50,0,1500,H_lead_Pt);
	fillHistogram("subleading higgs pT",50,0,1500,H_sublead_Pt);
	fillHistogram("MET from hh final states",30,0,1500,MET[0]);
      }
    }
    
    //WW
    Double_t W_lead_Pt, W_sublead_Pt;
    if(wcount==2 && hcount==0){
      W_lead_Pt    = Particle_PT[Windex[0]];
      W_sublead_Pt = Particle_PT[Windex[1]];
      if(Particle_PT[Windex[1]]>W_lead_Pt){
	W_lead_Pt    = Particle_PT[Windex[1]];
	W_sublead_Pt = Particle_PT[Windex[0]];
	fillHistogram("leading W pT",50,0,1500,W_lead_Pt);
	fillHistogram("subleading W pT",50,0,1500,W_sublead_Pt);
	fillHistogram("MET from WW final states",30,0,1500,MET[0]);
      }
    }
    
    //WH
    Double_t Single_W_Pt, Single_H_Pt;
    if(wcount==1 && hcount==1){
      Single_H_Pt = Particle_PT[hindex[0]];
      Single_W_Pt = Particle_PT[Windex[0]];
      fillHistogram("leading higgs pT(Wh)",50,0,1500,Single_H_Pt);
      fillHistogram("subleading W pT(Wh)",50,0,1500,Single_W_Pt);
      fillHistogram("MET from Wh final states",30,0,1500,MET[0]);
    }


    //counting final states
    if(wcount==2 && hcount==0) totalww++;
    else if(wcount==0 && hcount==2) totalhh++;
    else if(wcount==1 && hcount==1) totalwh++;
    else if(wcount==1 && hcount==0) totalw++;
    else if(wcount==0 && hcount==1) totalh++;
    else if(wcount==0 && hcount==0) totalnon++;
    else
      totalOther++;

    //cout << "event#" << ievent << "**" << "w count= " << wcount << " ; " << "h count= " << hcount << endl;
    if(wcount>2 || hcount>2){
      cout << "-------- event# --------" << ievent << endl;
      cout <<"**** wcount= " << wcount << " " << "**** hcount= " << hcount << endl;
      for(int j=0; j<nParticle; j++){
        if(Particle_Status[j]==3){
          cout << "Particle Adress= " << j << " ; "
               << "PID= "    << Particle_PID[j]    << " ; "
               << "Status= " << Particle_Status[j] << " ; "
               << "Mass= "   << Particle_Mass[j]   << " ; "
               << "Mother1 " << PARTICLE_M1[j]     << " ; "
               << "Mother2 " << PARTICLE_M2[j]     << endl;
	}
      }
    }

  }//End Event


  //Cross check a single event only contains final states of either WW Wh hh w, or h
  if(totalOther>0) cout << "Warning! " << "total Others= "<< totalOther << " "
                        <<"wcount= " << wcount << " " << "hcount= " << hcount << endl;

  if((totalww+ totalhh+ totalwh+ totalw+ totalh+ totalnon+ totalOther)/totevent != 1){
    cout << "ww count= " << totalww << " " << "or " << totalww/totevent << "; "
         << "hh count= " << totalhh << " " << "or " << totalhh/totevent << "; "
         << "wh count= " << totalwh << " " << "or " << totalwh/totevent << "; "
         << "w count=  " << totalw  << " " << "or " << totalw/totevent << "; "
         << "h count=  " << totalh  << " " << "or " << totalh/totevent << "; "
         << "zero h & w = " << totalnon << " "  << "or " << totalnon/totevent << "; "
         << "other count= " << totalOther << " " << "or " << totalOther/totevent << "; "
         << "sum all= " << (totalww+ totalhh+ totalwh+ totalw+ totalh+ totalnon+ totalOther)/totevent
         << endl;
  }

  cout << "totalww= " << totalww << ";   totalhh= " << totalhh << ";    "
       << "totalwh= " << totalwh << ";   totalw= "  << totalw  << ";    "
       << "totalh=  " << totalh  << endl;
  
  cout << "  " << endl;
  cout << "Sum all= " << (totalww+ totalhh+ totalwh+ totalw+ totalh+ totalnon+ totalOther) << endl;

  cout << "store histogram info ..." << endl;
  //Save Histograms                       
  int ipic = 0;
  int nhist = 14;
  TCanvas *c[nhist];
  TFile* fhistout=new TFile("hist_output.root","recreate");
  for( const auto& sm_pair : h1f_map ){
    TH1F* h1f = (TH1F*) sm_pair.second;
    h1f->Write();

    ipic +=1;
    c[ipic] = new TCanvas(Form("c%d",ipic));
    h1f->Draw();
    c[ipic]->SetLogy();
    c[ipic]->SaveAs();
  }
  fhistout->Close();

  

  return 0;
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

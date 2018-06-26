#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"

#include<iostream>
#include<string>
#include <algorithm>

using namespace std;

int kinematics_2D_local(){


  Int_t n[1]={1200};//,2000,2400};                                         
  Int_t ik[1]={5};
  Int_t fnum = 1;
  Int_t knum = 1;

  TFile* fhistout=new TFile("hist_output.root","recreate");

  //Canvas
  TCanvas *c;
  c = new TCanvas("c","c",900,900);
  TCanvas *c2;
  c2 = new TCanvas("c2","c2",900,900);  

  TH2 *hpt = new TH2D("higgs pT", "higgs pT", 15, 0, 1500, \
		                              15, 0, 1500);
  TH2 *heta  = new TH2D("higgs eta", "higgs eta", 30, -5, 5, \
                             		        30, -8, 8); 


  //Local Truth Particle Kinematic Study
  TChain *t = new TChain("Delphes");
  // t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/HHmet/Validation/13processHHSample_f1200_k05/f1200k5_DATA/dihiggs*.root");
  t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/HHmet/Validation/13processHHSample_f1200_k05/f1200k5_DATA/dihiggs*.root");
  Long64_t nevent = t->GetEntries();
  cout << "total events"  << nevent << endl;
  t->SetMakeClass(1);
  

  //Setup Branches      
  const Int_t    MaxnMET = 5000;
  Int_t          nMET = -1;
  Float_t        MET[MaxnMET];
  TBranch        *met;
  t->SetBranchAddress("MissingET", &nMET);
  t->SetBranchAddress("MissingET.MET", MET, &met);
  
  const Int_t     kMaxParticle = 5000;
  Int_t           nParticle = -1;
  Int_t           Particle_PID[kMaxParticle];
  Int_t           Particle_Status[kMaxParticle];
  Float_t         Particle_PT[kMaxParticle];
  Float_t         Particle_Eta[kMaxParticle];
  Float_t         Particle_Phi[kMaxParticle];
  Float_t         Particle_Mass[kMaxParticle];
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
    particle_pid->GetEntry(ievent);
    particle_status->GetEntry(ievent);
    particle_pt->GetEntry(ievent);
    particle_eta->GetEntry(ievent);
    particle_phi->GetEntry(ievent);
    particle_mass->GetEntry(ievent);
    particle_m1->GetEntry(ievent);
    particle_m2->GetEntry(ievent);
    
    //Build table              
    wcount=0, hcount=0;
    vector<int> WhIndex;
    vector<int> ZhIndex;
    vector<int> Windex;
    vector<int> hindex;
    //loop over address of Wh and Zh particles                            
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
    Double_t H_lead_Eta, H_sublead_Eta;
    if(wcount==0 && hcount==2){
      H_lead_Pt    = Particle_PT[hindex[0]];
      H_sublead_Pt = Particle_PT[hindex[1]];
      if(Particle_PT[hindex[1]]>H_lead_Pt){
	H_lead_Pt      = Particle_PT[hindex[1]];
	H_sublead_Pt   = Particle_PT[hindex[0]];
	H_lead_Eta     = Particle_Eta[hindex[1]]; 
	H_sublead_Eta  = Particle_Eta[hindex[0]];
	hpt->Fill(H_lead_Pt, H_sublead_Pt);
	heta->Fill(H_lead_Eta, H_sublead_Eta);        
      }
    }
    
    /*
    //WH                
    Double_t Single_W_Pt, Single_H_Pt;
    if(wcount==1 && hcount==1){
    Single_H_Pt = Particle_PT[hindex[0]];
    Single_W_Pt = Particle_PT[Windex[0]];
    }
    
    //WW     
    Double_t W_lead_Pt, W_sublead_Pt;
    if(wcount==2 && hcount==0){
    W_lead_Pt    = Particle_PT[Windex[0]];
    W_sublead_Pt = Particle_PT[Windex[1]];
    if(Particle_PT[Windex[1]]>W_lead_Pt){
    W_lead_Pt    = Particle_PT[Windex[1]];
    W_sublead_Pt = Particle_PT[Windex[0]];
    }
    }
    */
    
    
  }//End Event
  
  hpt->Write();
  hpt->SetStats(0);
  heta->Write();
  heta->SetStats(0);
  
  
 
  
  //Write Canvas
  c->cd();
  hpt->Draw("colz");
  hpt->GetXaxis()->SetTitle("Leading Higgs pT.[GeV]");
  hpt->GetYaxis()->SetTitle("Subleading Higgs pT.[GeV]");
  hpt->GetYaxis()->SetTitleOffset(1.55);
  c->SaveAs();
  
  c2->cd();
  heta->Draw("colz");
  heta->GetXaxis()->SetTitle("Leading Higgs eta.[GeV]");
  heta->GetYaxis()->SetTitle("Subleading Higgs eta.[GeV]");
  heta->GetYaxis()->SetTitleOffset(1.59);
  c2->SaveAs();


  fhistout->Close();




  return 0;
}

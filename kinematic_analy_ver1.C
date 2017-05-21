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

using namespace std;

int kinematic_analy_ver1() {
  // Set debug flag
  Int_t flag = 0;

  TFile *fin=new TFile("/Users/yuwahlim/Desktop/DataChallenge/Cut_analysis/hhdata/1000k1/tag_1_delphes_events.root","READONLY");
  TFile *fout=new TFile("kinematics.root","RECREATE");

  // HISTOGRAM SETUP FOR GenParticle
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",50, 0, 1500);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",50, 0, 1500);
  TH1D *dR1_bb = new TH1D("dR1_bb", "leading_dR", 50, 0, 5);
  TH1D *dR2_bb = new TH1D("dR2_bb", "subleading_dR", 50, 0, 5);

  // READ EVENT FROM DELPHES
  TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  cout << nevent << endl;
  t->SetMakeClass(1); 
  cout << "Total events to read = " << nevent << endl;

  //// GET BRANCH TO GenParticle ////
  Int_t nPTCL=-1;
  Int_t PTCL_PID[nevent];
  Int_t PTCL_M1[nevent];
  Float_t PTCL_PT[nevent];
  Float_t PTCL_ETA[nevent];
  Float_t PTCL_PHI[nevent];
  Float_t PTCL_M[nevent];

  t->SetBranchAddress("Particle", &nPTCL);
  TBranch *ptcl_pid = t->GetBranch("Particle.PID");
  ptcl_pid->SetAddress(&PTCL_PID);
  TBranch *ptcl_m1 = t->GetBranch("Particle.M1");
  ptcl_m1->SetAddress(&PTCL_M1);
  TBranch *ptcl_pt = t->GetBranch("Particle.PT");
  ptcl_pt->SetAddress(&PTCL_PT);
  TBranch *ptcl_eta = t->GetBranch("Particle.Eta");
  ptcl_eta->SetAddress(&PTCL_ETA);
  TBranch *ptcl_phi = t->GetBranch("Particle.Phi");
  ptcl_phi->SetAddress(&PTCL_PHI);
  TBranch *ptcl_m = t->GetBranch("Particle.Mass");
  ptcl_m->SetAddress(&PTCL_M);

  // event loop
  Int_t nevents = nevent;
  for(Long64_t ievent = 0; ievent < nevents; ievent++){    
    ptcl_pid->GetEntry(ievent);
    ptcl_m1->GetEntry(ievent);
    ptcl_pt->GetEntry(ievent);
    ptcl_eta->GetEntry(ievent);
    ptcl_phi->GetEntry(ievent);
    ptcl_m->GetEntry(ievent);
    
    cout << "running... event(" << ievent << ")" << endl;
    // GenParticle loop counting h and b
    Int_t h_num = 0, b_num = 0;
    for(int j = 0; j < nPTCL; j++){
      if(PTCL_PID[j] == 25) { h_num +=1; }
      if(PTCL_PID[j] == 5 || PTCL_PID[j] == -5) { b_num +=1; }
    } // end counting
    
    Int_t hindex[h_num], bindex[b_num]; 
    // get particles index
    Int_t index1 = -1, index2 = -1; 
    for(int j = 0; j < nPTCL; j++){
      // higgs
      if(PTCL_PID[j] == 25){
	index1 +=1; hindex[index1] = j;
	if(flag == 1){
	cout << " hindex " << hindex[index1]			\
	     << " index1 " << index1				\
	     << "  PID  "<< PTCL_PID[hindex[index1]] << endl;   
	}
      }
      // b-tag
      if(PTCL_PID[j] == 5 || PTCL_PID[j] == -5){
	index2 +=1; bindex[index2] = j;
	if(flag == 2){
	cout << " b_mother " << PTCL_M1[bindex[index2]]	\
	     <<" index2 " << index2				\
	     << "  PID  " << PTCL_PID[bindex[index2]] << endl;
	}
      }      
    }

    Int_t index3 = -1, index4 = -1;
    Int_t new_hindex[h_num], new_bindex[h_num];
    Double_t truth_hpt[h_num];
    Double_t leading_hpt, subleading_hpt;
    for(int j  = 0; j < h_num; j++){
      for(int i = 0; i < b_num; i++){
	if(hindex[j] == PTCL_M1[bindex[i]]){
	  index3 +=1;
	  new_hindex[index3] = hindex[j];
	  new_bindex[index3] = bindex[i];
	  if(flag == 3){
          cout << "show hindex " << hindex[j] << endl;    
          cout << "show bindex " << bindex[i] << endl;  
	  cout << "new_hindex " << new_hindex[index3]		\
	       << " new_bindex " << new_bindex[index3]		\
	       << " index3 " << index3 << endl; 
	  }
	}
      }
    }

    TLorentzVector lead_b_vecs[2], sub_b_vecs[2];
    Double_t delta_R1, delta_R2;
    if(PTCL_PT[new_hindex[0]] > PTCL_PT[new_hindex[2]]){
      leading_hpt = PTCL_PT[new_hindex[0]];
      subleading_hpt = PTCL_PT[new_hindex[2]];
      lead_b_vecs[0].SetPtEtaPhiM(PTCL_PT[new_bindex[0]],PTCL_ETA[new_bindex[0]], \
				   PTCL_PHI[new_bindex[0]],PTCL_M[new_bindex[0]]);
      lead_b_vecs[1].SetPtEtaPhiM(PTCL_PT[new_bindex[1]],PTCL_ETA[new_bindex[1]], \
				   PTCL_PHI[new_bindex[1]],PTCL_M[new_bindex[1]]);
      sub_b_vecs[0].SetPtEtaPhiM(PTCL_PT[new_bindex[2]],PTCL_ETA[new_bindex[2]],	\
				  PTCL_PHI[new_bindex[2]],PTCL_M[new_bindex[2]]);
      sub_b_vecs[1].SetPtEtaPhiM(PTCL_PT[new_bindex[3]],PTCL_ETA[new_bindex[3]],	\
				  PTCL_PHI[new_bindex[3]],PTCL_M[new_bindex[3]]);
      if(flag == 4){
        cout << "new bindex " <<  new_bindex[0] << " " << new_bindex[1] <<endl;
        cout << "b propertie " << PTCL_PT[new_bindex[1]] << " "         \
             << PTCL_ETA[new_bindex[1]] << " " << PTCL_PHI[new_bindex[1]] << " " \
             << PTCL_M[new_bindex[1]] << endl;
      }
    } else {
      leading_hpt = PTCL_PT[new_hindex[2]];
      subleading_hpt = PTCL_PT[new_hindex[1]];
      lead_b_vecs[0].SetPtEtaPhiM(PTCL_PT[new_bindex[2]],PTCL_ETA[new_bindex[2]], \
                                   PTCL_PHI[new_bindex[2]],PTCL_M[new_bindex[2]]);
      lead_b_vecs[1].SetPtEtaPhiM(PTCL_PT[new_bindex[3]],PTCL_ETA[new_bindex[3]], \
                                   PTCL_PHI[new_bindex[3]],PTCL_M[new_bindex[3]]);
      sub_b_vecs[0].SetPtEtaPhiM(PTCL_PT[new_bindex[0]],PTCL_ETA[new_bindex[0]], \
                                  PTCL_PHI[new_bindex[0]],PTCL_M[new_bindex[0]]);
      sub_b_vecs[1].SetPtEtaPhiM(PTCL_PT[new_bindex[1]],PTCL_ETA[new_bindex[1]], \
                                  PTCL_PHI[new_bindex[1]],PTCL_M[new_bindex[1]]);
      if(flag == 4){
	cout << "new bindex " <<  new_bindex[0] << " " << new_bindex[1] <<endl; 
	cout << "b propertie " << PTCL_PT[new_bindex[1]] << " "		\
	     << PTCL_ETA[new_bindex[1]] << " " << PTCL_PHI[new_bindex[1]] << " " \
	     << PTCL_M[new_bindex[1]] << endl;
      }
    }

    delta_R1 = abs(lead_b_vecs[0].DeltaR(lead_b_vecs[1]));
    delta_R2 = abs(sub_b_vecs[0].DeltaR(sub_b_vecs[1]));
    
    // Fill histogram
    h_leading_pt->Fill(leading_hpt);
    h_subleading_pt->Fill(subleading_hpt);
    dR1_bb->Fill(delta_R1);
    dR2_bb->Fill(delta_R2);
  } // end event loop

  fout->cd();
  h_leading_pt->Write();
  h_subleading_pt->Write();
  dR1_bb->Write();
  dR2_bb->Write();
  
  cout << "mission completed!" << endl;

  return 0;
}

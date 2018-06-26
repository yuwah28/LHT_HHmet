#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector2.h"
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
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

//global variables:
//Boost 
vector<TLorentzVector> akt10trimTLV;
vector<TLorentzVector> akt2trackTLV;
vector<UInt_t> akt2trackBTag;
vector<double> selectLargeRjetPt;
Double_t dEtahh;
Double_t Xhh;
Int_t eventnum;

bool passBoost(Long64_t* cutflow);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;


int hh_rso_analysis_done(){

  TChain *t = new TChain("Delphes");
  t->Add("/phys/groups/tev/scratch1/users/paklim/SMBackground/ttbarHadron/had*.root");
 
  Long64_t nevent = t->GetEntries();
  cout << "Total Event = " << nevent << endl;
  t->SetMakeClass(1);
  //----------------SETUP BRANCHES-----------------------                                                                             
  const Int_t    MaxnMET = 50;
  Int_t          nMET = -1;
  Float_t        MET[MaxnMET];
  Float_t        MET_PHI[MaxnMET];
  TBranch        *met;
  TBranch        *met_phi;
  t->SetBranchAddress("MissingET", &nMET);
  t->SetBranchAddress("MissingET.MET", MET, &met);
  t->SetBranchAddress("MissingET.Phi", MET_PHI, &met_phi);

  const Int_t    MaxnHT = 50;
  Int_t          nHT = -1;
  Float_t        ScalarHT[MaxnHT];
  TBranch        *scalar_ht;
  t->SetBranchAddress("ScalarHT", &nHT);
  t->SetBranchAddress("ScalarHT.HT", ScalarHT, &scalar_ht);

  const Int_t    MaxnJet = 1000;
  Int_t          nJET = -1;
  Float_t        JET_PT[MaxnJet];
  Float_t        JET_ETA[MaxnJet];
  Float_t        JET_PHI[MaxnJet];
  Float_t        JET_M[MaxnJet];
  UInt_t         JET_BTAG[MaxnJet];
  TLorentzVector all_vector[MaxnJet];
  TBranch        *jet_pt;
  TBranch        *jet_eta;
  TBranch        *jet_phi;
  TBranch        *jet_m;
  TBranch        *jet_btag;
  t->SetBranchAddress("Jet", &nJET);
  t->SetBranchAddress("Jet.PT", JET_PT, &jet_pt);
  t->SetBranchAddress("Jet.Eta", JET_ETA, &jet_eta);
  t->SetBranchAddress("Jet.Phi", JET_PHI, &jet_phi);
  t->SetBranchAddress("Jet.Mass", JET_M, &jet_m);
  t->SetBranchAddress("Jet.BTag", JET_BTAG, &jet_btag);

  const          Int_t MaxTrkNum = 10000;
  Int_t          nTRACK = -1;
  Int_t          TRACK_CHR[MaxTrkNum];
  Float_t        TRACK_PT[MaxTrkNum];
  Float_t        TRACK_ETA[MaxTrkNum];
  Float_t        TRACK_PHI[MaxTrkNum];
  TBranch        *track_pt;
  TBranch        *track_eta;
  TBranch        *track_phi ;
  TBranch        *track_charge;
  t->SetBranchAddress("Track", &nTRACK);
  t->SetBranchAddress("Track.PT", TRACK_PT, &track_pt);
  t->SetBranchAddress("Track.Eta", TRACK_ETA, &track_eta);
  t->SetBranchAddress("Track.Phi", TRACK_PHI, &track_phi);
  t->SetBranchAddress("Track.Charge", TRACK_CHR, &track_charge);

  const Int_t     kMaxJetAK10Trim = 3000;
  Int_t           nJetAK10Trim = nJetAK10Trim;
  Float_t         JetAK10Trim_PT[kMaxJetAK10Trim];   //[JetAK10Trim_]    
  Float_t         JetAK10Trim_Eta[kMaxJetAK10Trim];   //[JetAK10Trim_]       
  Float_t         JetAK10Trim_Phi[kMaxJetAK10Trim];   //[JetAK10Trim_]        
  Float_t         JetAK10Trim_Mass[kMaxJetAK10Trim];   //[JetAK10Trim_]    
  TBranch        *b_JetAK10Trim_PT;     
  TBranch        *b_JetAK10Trim_Eta;       
  TBranch        *b_JetAK10Trim_Phi;        
  TBranch        *b_JetAK10Trim_Mass;  
  t->SetBranchAddress("JetAK10Trim", &nJetAK10Trim);
  t->SetBranchAddress("JetAK10Trim.PT", JetAK10Trim_PT, &b_JetAK10Trim_PT);
  t->SetBranchAddress("JetAK10Trim.Eta", JetAK10Trim_Eta, &b_JetAK10Trim_Eta);
  t->SetBranchAddress("JetAK10Trim.Phi", JetAK10Trim_Phi, &b_JetAK10Trim_Phi);
  t->SetBranchAddress("JetAK10Trim.Mass", JetAK10Trim_Mass, &b_JetAK10Trim_Mass);

  const Int_t     kMaxJetAK2Track = 3000;
  Int_t           nJetAK2Track = nJetAK2Track;
  Float_t         JetAK2Track_PT[kMaxJetAK2Track];   //[JetAK2Track_] 
  Float_t         JetAK2Track_Eta[kMaxJetAK2Track];   //[JetAK2Track_]     
  Float_t         JetAK2Track_Phi[kMaxJetAK2Track];   //[JetAK2Track_]   
  Float_t         JetAK2Track_T[kMaxJetAK2Track];   //[JetAK2Track_]   
  Float_t         JetAK2Track_Mass[kMaxJetAK2Track];   //[JetAK2Track_]   
  UInt_t          JetAK2Track_BTag[kMaxJetAK2Track];   //[JetAK2Track_]    
  Int_t           JetAK2Track_Flavor[kMaxJetAK2Track];
  TBranch        *b_JetAK2Track_PT; 
  TBranch        *b_JetAK2Track_Eta;  
  TBranch        *b_JetAK2Track_Phi;   
  TBranch        *b_JetAK2Track_Mass;  
  TBranch        *b_JetAK2Track_BTag;  
  TBranch        *b_JetAK2Track_Flavor;
  t->SetBranchAddress("JetAK2Track", &nJetAK2Track);
  t->SetBranchAddress("JetAK2Track.PT", JetAK2Track_PT, &b_JetAK2Track_PT);
  t->SetBranchAddress("JetAK2Track.Eta", JetAK2Track_Eta, &b_JetAK2Track_Eta);
  t->SetBranchAddress("JetAK2Track.Phi", JetAK2Track_Phi, &b_JetAK2Track_Phi);
  t->SetBranchAddress("JetAK2Track.Mass", JetAK2Track_Mass, &b_JetAK2Track_Mass);
  t->SetBranchAddress("JetAK2Track.BTag", JetAK2Track_BTag, &b_JetAK2Track_BTag);
  t->SetBranchAddress("JetAK2Track.Flavor", JetAK2Track_Flavor, &b_JetAK2Track_Flavor);


  ///////// Start Cut Analysis //////////      
  Long64_t cut_flow[20]; std::fill_n(cut_flow, 20, 0);
  //nevent=500000;                               
  int NDivision = nevent/10;
  for(Long64_t jevent = 0; jevent < nevent; jevent++){
    
    
    if(jevent%NDivision ==0) cout << jevent << "/" << nevent <<endl;
    cut_flow[0]++;
    
    Long64_t ievent = t->LoadTree(jevent);
    met->GetEntry(ievent);
    met_phi->GetEntry(ievent);
    scalar_ht->GetEntry(ievent);
    jet_pt->GetEntry(ievent);
    jet_phi->GetEntry(ievent);
    jet_eta->GetEntry(ievent);
    jet_m->GetEntry(ievent);
    jet_btag->GetEntry(ievent);
    
    track_pt->GetEntry(ievent);
    track_eta->GetEntry(ievent);
    track_phi->GetEntry(ievent);
    track_charge->GetEntry(ievent);

    b_JetAK10Trim_PT->GetEntry(ievent);
    b_JetAK10Trim_Eta->GetEntry(ievent);
    b_JetAK10Trim_Phi->GetEntry(ievent);
    b_JetAK10Trim_Mass->GetEntry(ievent);
    b_JetAK2Track_PT->GetEntry(ievent);
    b_JetAK2Track_Eta->GetEntry(ievent);
    b_JetAK2Track_Phi->GetEntry(ievent);
    b_JetAK2Track_Mass->GetEntry(ievent);
    b_JetAK2Track_BTag->GetEntry(ievent);

    eventnum = ievent;
    
    //Boosted Analysis
					       
    //Large-R jet
    akt10trimTLV.clear();
    selectLargeRjetPt.clear();
    vector<int> LIndex;
    for(int j=0;j<nJetAK10Trim ;j++){
      if(JetAK10Trim_PT[j]>250 && fabs(JetAK10Trim_Eta[j])<2.0 && JetAK10Trim_Mass[j]>50){
        selectLargeRjetPt.push_back(JetAK10Trim_PT[j]);

	if(selectLargeRjetPt[0]>450){
	  LIndex.push_back(j);
	  TLorentzVector tlv;
	  tlv.SetPtEtaPhiM(JetAK10Trim_PT[j],JetAK10Trim_Eta[j],JetAK10Trim_Phi[j],JetAK10Trim_Mass[j]);
	  akt10trimTLV.push_back(tlv);
	  fillHistogram("akt10M",100,0,300,tlv.M());
       
	}
      }
    }
   

    Xhh=0;
    if(akt10trimTLV.size()>1){
      //Eta separation between hh 
      dEtahh = fabs(JetAK10Trim_Eta[LIndex[0]]-JetAK10Trim_Eta[LIndex[1]]);
      fillHistogram("Delta Etahh",50,-5,5,dEtahh);

      //Resonance mass
      Double_t mjLead = JetAK10Trim_Mass[LIndex[0]];
      Double_t mjSub  = JetAK10Trim_Mass[LIndex[1]];
      Xhh = sqrt( pow((mjLead-124.)/(0.1*mjLead),2)+pow((mjSub-115.)/(0.1*mjSub),2) );
      fillHistogram("hh mass windown",100,0,2000,Xhh);
    }
    
    //Small-R 0.2 Track-jet
    akt2trackTLV.clear();
    akt2trackBTag.clear();
    for(int j=0;j<nJetAK2Track ;j++){
      if(JetAK2Track_PT[j]>10  && fabs(JetAK2Track_Eta[j])<2.5){
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(JetAK2Track_PT[j],JetAK2Track_Eta[j],JetAK2Track_Phi[j],JetAK2Track_Mass[j]);
        akt2trackTLV.push_back(tlv);
        akt2trackBTag.push_back(JetAK2Track_BTag[j]);
      }
    }
    

    passBoost(&cut_flow[0]);

  }// end of events 

  cout << "store cut info ..." << endl;


  ofstream fileoutB("cutinfoB_local.dat");
  TString cutnameR[8]={"no selection", "large-R jet >=2", "|Delta_Eta_hh|<1.7",
		       "b-tagged_track-jets>=2", "Xhh<1.6", "2b-tagge track-jets",
		       "3b-tagge track-jets", "4b-tagge track-jets"};

  for(int j=0;j<8;j++) fileoutB << "[" << j << "] " << cutnameR[j] << ": " << cut_flow[j] <<endl;
  fileoutB.close();


  ///////// Save histogram ///////////                    
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
    //c[ipic]->SetLogy();        
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


bool passBoost(Long64_t* cut_flow){
  
  if(akt10trimTLV.size()<2) return false;
  cut_flow[1]++;
  
  if(dEtahh>=1.7) return false;
  cut_flow[2]++;
  
  //b-tagging: one or two ID track jets matched to large Rjet    
  int Nbtag_track[2]={0,0};
  if(akt10trimTLV.size()>1){
    for(int j=0; j<akt2trackTLV.size(); j++){
      double minR = 1.0;
      int match = -1;
      for(int i=0; i<2; i++){
	if(akt10trimTLV[i].DeltaR(akt2trackTLV[j]) <minR){
	  minR=akt10trimTLV[i].DeltaR(akt2trackTLV[j]);
	  match=i;
	}
      }
      if(match==-1) continue;
      if(Bool_t(akt2trackBTag[j] &(1<<1))==false) continue;
      Nbtag_track[match]++;
     
    }
  }
  
  if(Nbtag_track[0]<1 || Nbtag_track[1]<1 ) return false;
  cut_flow[3]++;

  if(Xhh>=1.6) return false;
  cut_flow[4]++;


  Int_t numbtag=0;
  if(Nbtag_track[0]>0 && Nbtag_track[1]>0){
    numbtag = Nbtag_track[0] + Nbtag_track[1];
  }

  //if(numbtag!=2) return false;
  if(numbtag==2) cut_flow[5]++;
  
 
  //if(numbtag!=3) return false; 
  if(numbtag==3) cut_flow[6]++;


  //if(Nbtag_track[0]!=2 || Nbtag_track[1]!=2) return false;
  if(Nbtag_track[0]==2 && Nbtag_track[1]==2) cut_flow[7]++;
  

  
}

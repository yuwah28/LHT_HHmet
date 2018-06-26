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
#include "TRandom3.h"
#include "TCanvas.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

bool passBoost(Long64_t* cutbinaryB);
bool passResolved(Long64_t* cutbinaryR);
bool passHybrid(Long64_t* cutflow, Long64_t* cutflow2, Long64_t* cutbinaryR, Long64_t* cutbinaryB, Long64_t* cutflowH);
int getBosonTag(TLorentzVector reco, vector<TLorentzVector>& truth, float coneR=1.0);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;
TRandom3 randGen;


//global variables:        
Int_t eventnum;
Int_t resolveCutSize;
Int_t boostCutSize;

//Resolved         
Int_t num_Lep;
Int_t num_JET;
Int_t num_bJET;
Int_t Nb_M;
Int_t Nb_L;
Int_t iso_flag;
Int_t dphi_flag;
Double_t ETmiss ;
Double_t hmass;
Double_t Wmass;
Int_t dRjj_flag;
Int_t dRbb_flag;
vector<TLorentzVector> selectJetsLV;

//Boosted
Double_t Lhmass;
Double_t LWmass;
vector<TLorentzVector> akt10trimTLV;
vector<TLorentzVector> akt2trackTLV;
vector<TLorentzVector> HiggsCandidateTLV;
vector<TLorentzVector> BosonCandidateTLV;

vector<int> Nakt10trimBosonTag;
vector<UInt_t> akt2trackBTag;
vector<UInt_t> akt10trimBosonTag;
vector<double> selectLargeRjetPt;
Int_t Windex;


int wh_had_HybridH(){

  randGen.SetSeed(2);

  TChain *t = new TChain("Delphes");

  //4 Processes Sample at ffnum kknum
  //t->Add("/phys/groups/tev/scratch1/users/paklim/Samples/4processSample_ffnumk5/ffnumk5_Hybrid_WH_DATA/Hdihiggs*.root");

  //13 Processes Sample at ffnum kknum
  //t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/HHmet/Validation/13processHHSample_ffnum_kknum/WH_Hybrid_ffnum/dihiggs*.root");

  //t->Add("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/ffnum/CHNAMEkknum/WHhybr.root");


  //t->Add("/phys/groups/tev/scratch1/users/paklim/HLSamples/WH_HH_met/WHffnumkknum.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/HLSamples/WH_HH_met/WHffnumkknuma.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/HLSamples/WH_HH_met/WHffnumkknumb.root");

  t->Add("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/\
ffnum/CHNAMEkknum/WHhybr.root");

  Long64_t nevent = t->GetEntries();
  cout << nevent << endl;
  t->SetMakeClass(1);

  //----------------SETUP BRANCHES----------------------- 
  const Int_t MaxnMET = 5000;
  Int_t nMET = -1;
  Float_t MET[MaxnMET];
  Float_t MET_PHI[MaxnMET];
  TBranch *met;
  TBranch *met_phi ;
  t->SetBranchAddress("MissingET", &nMET);
  t->SetBranchAddress("MissingET.MET", MET, &met);
  t->SetBranchAddress("MissingET.Phi", MET_PHI, &met_phi);

  const Int_t MaxnJet = 10000;
  Int_t nJET = -1;
  Float_t JET_PT[MaxnJet];
  Float_t JET_ETA[MaxnJet];
  Float_t JET_PHI[MaxnJet];
  Float_t JET_M[MaxnJet];
  UInt_t JET_BTAG[MaxnJet];
  TLorentzVector all_vector[MaxnJet];
  TBranch *jet_pt;
  TBranch *jet_eta;
  TBranch *jet_phi;
  TBranch *jet_m;
  TBranch *jet_btag;
  t->SetBranchAddress("Jet", &nJET);
  t->SetBranchAddress("Jet.PT",   JET_PT  , &jet_pt);
  t->SetBranchAddress("Jet.Eta",  JET_ETA , &jet_eta);
  t->SetBranchAddress("Jet.Phi",  JET_PHI , &jet_phi);
  t->SetBranchAddress("Jet.Mass", JET_M   , &jet_m);
  t->SetBranchAddress("Jet.BTag", JET_BTAG, &jet_btag);

  const Int_t MaxnELEC = 200;
  Int_t nElec = -1;
  Float_t ELEC_PT[MaxnELEC];
  Float_t ELEC_ETA[MaxnELEC];
  Float_t ELEC_PHI[MaxnELEC];
  Float_t ELEC_SUMPT[MaxnELEC];
  TBranch  *e_pt;
  TBranch  *e_eta;
  TBranch  *e_phi;
  TBranch  *e_sumpt;
  t->SetBranchAddress("Electron", &nElec);
  t->SetBranchAddress("Electron.PT",           ELEC_PT,    &e_pt);
  t->SetBranchAddress("Electron.Eta",          ELEC_ETA,   &e_eta);
  t->SetBranchAddress("Electron.Phi",          ELEC_PHI,   &e_phi);
  t->SetBranchAddress("Electron.SumPtCharged", ELEC_SUMPT, &e_sumpt);

  const Int_t MaxnMUON = 200;
  Int_t nMUON = -1;
  Float_t MUON_PT[MaxnMUON];
  Float_t MUON_ETA[MaxnMUON];
  Float_t MUON_PHI[MaxnMUON];
  Float_t MUON_SUMPT[MaxnMUON];
  TBranch *muon_pt;
  TBranch *muon_eta;
  TBranch *muon_phi;
  TBranch *muon_sumpt;
  t->SetBranchAddress("Muon", &nMUON);
  t->SetBranchAddress("Muon.PT",            MUON_PT,     &muon_pt);
  t->SetBranchAddress("Muon.Eta",           MUON_ETA,    &muon_eta);
  t->SetBranchAddress("Muon.Phi",           MUON_PHI,    &muon_phi);
  t->SetBranchAddress("Muon.SumPtCharged",  MUON_SUMPT,  &muon_sumpt);

  const          Int_t MaxTrkNum = 10000;
  Int_t          nTRACK = -1;
  Float_t        TRACK_PT[MaxTrkNum];
  Float_t        TRACK_ETA[MaxTrkNum];
  Float_t        TRACK_PHI[MaxTrkNum];
  Float_t        TRACK_CHR[MaxTrkNum];
  TBranch        *track_pt;
  TBranch        *track_eta;
  TBranch        *track_phi ;
  TBranch        *track_charge;
  t->SetBranchAddress("Track", &nTRACK);
  t->SetBranchAddress("Track.PT", TRACK_PT, &track_pt);
  t->SetBranchAddress("Track.Eta", TRACK_ETA, &track_eta);
  t->SetBranchAddress("Track.Phi", TRACK_PHI, &track_phi);
  t->SetBranchAddress("Track.Charge", &track_charge);

  const Int_t kMaxJetAK10Trim = 10000;
  Int_t nJetAK10Trim = nJetAK10Trim;
  Float_t         JetAK10Trim_PT[kMaxJetAK10Trim];
  Float_t         JetAK10Trim_Eta[kMaxJetAK10Trim];
  Float_t         JetAK10Trim_Phi[kMaxJetAK10Trim];
  Float_t         JetAK10Trim_Mass[kMaxJetAK10Trim];
  TBranch        *b_JetAK10Trim_PT;
  TBranch        *b_JetAK10Trim_Eta;
  TBranch        *b_JetAK10Trim_Phi;
  TBranch        *b_JetAK10Trim_Mass;
  t->SetBranchAddress("JetAK10Trim", &nJetAK10Trim);
  t->SetBranchAddress("JetAK10Trim.PT", JetAK10Trim_PT, &b_JetAK10Trim_PT);
  t->SetBranchAddress("JetAK10Trim.Eta", JetAK10Trim_Eta, &b_JetAK10Trim_Eta);
  t->SetBranchAddress("JetAK10Trim.Phi", JetAK10Trim_Phi, &b_JetAK10Trim_Phi);
  t->SetBranchAddress("JetAK10Trim.Mass", JetAK10Trim_Mass, &b_JetAK10Trim_Mass);

  const Int_t kMaxJetAK2Track = 1000;
  Int_t nJetAK2Track = nJetAK2Track;
  Float_t         JetAK2Track_PT[kMaxJetAK2Track];
  Float_t         JetAK2Track_Eta[kMaxJetAK2Track];
  Float_t         JetAK2Track_Phi[kMaxJetAK2Track];
  Float_t         JetAK2Track_T[kMaxJetAK2Track];
  Float_t         JetAK2Track_Mass[kMaxJetAK2Track];
  UInt_t          JetAK2Track_BTag[kMaxJetAK2Track];
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

  const Int_t     kMaxParticle = 5000;
  Int_t           nParticle = -1;
  Int_t           Particle_PID[kMaxParticle];
  Int_t           Particle_Status[kMaxParticle];
  Float_t         Particle_PT[kMaxParticle];
  Float_t         Particle_Eta[kMaxParticle];
  Float_t         Particle_Phi[kMaxParticle];
  Float_t         Particle_Mass[kMaxParticle];
  Int_t           Particle_M1[kMaxParticle];
  Int_t           Particle_M2[kMaxParticle];
  Int_t           Particle_D1[kMaxParticle];
  Int_t           Particle_D2[kMaxParticle];
  TBranch        *particle_pid;
  TBranch        *particle_status;
  TBranch        *particle_pt;
  TBranch        *particle_eta;
  TBranch        *particle_phi;
  TBranch        *particle_mass;
  TBranch        *particle_m1;
  TBranch        *particle_m2;
  TBranch        *particle_d1;
  TBranch        *particle_d2;
  t->SetBranchAddress("Particle", &nParticle);
  t->SetBranchAddress("Particle.PID", Particle_PID, &particle_pid);
  t->SetBranchAddress("Particle.Status", Particle_Status, &particle_status);
  t->SetBranchAddress("Particle.PT", Particle_PT, &particle_pt);
  t->SetBranchAddress("Particle.Eta", Particle_Eta, &particle_eta);
  t->SetBranchAddress("Particle.Phi", Particle_Phi, &particle_phi);
  t->SetBranchAddress("Particle.Mass", Particle_Mass, &particle_mass);
  t->SetBranchAddress("Particle.M1", Particle_M1, &particle_m1);
  t->SetBranchAddress("Particle.M2", Particle_M2, &particle_m2);
  t->SetBranchAddress("Particle.D1", Particle_D1, &particle_d1);
  t->SetBranchAddress("Particle.D2", Particle_D2, &particle_d2);



  //nevent=200000;
  //Truth counter parameters  
  int wcount, hcount;
  int totalww=0, totalhh=0, totalwh=0, totalw=0, totalh=0, totalnon=0;
  int totalOther=0;

  int totalWtag=0;
  int totalW=0;
  int totalH=0;
  int totalHtag=0;
  int zeroWtag=0;
  int oneWtag=0;
  int total_reco_Hbb=0;
  int total_truth_Hbb=0;
  int total_reco_Wjj=0;
  int total_truth_Wjj=0;

  //Define Cutflow Entries
  Long64_t resolvedCutFlow[20]; std::fill_n(resolvedCutFlow, 20, 0);
  Long64_t boostedCutFlow[20]; std::fill_n(boostedCutFlow, 20, 0);
  Long64_t combinedCutFlow[20]; std::fill_n(combinedCutFlow, 20, 0);

  //Long64_t cut_flow[20];  std::fill_n(cut_flow, 20, 0);
  //Long64_t cut_flow_boost[20];  std::fill_n(cut_flow_boost, 20, 0);
  int NDivision = nevent/10;
  // Loop over events
  for(Long64_t jevent = 0; jevent < nevent; jevent++){

    if(jevent%NDivision ==0) cout << jevent << "/" << nevent <<endl;
    //cut_flow[0]++;
    //cut_flow_boost[0]++;
  
    Long64_t ievent = t->LoadTree(jevent);
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

    muon_pt->GetEntry(ievent);
    muon_eta->GetEntry(ievent);
    muon_phi->GetEntry(ievent);
    muon_sumpt->GetEntry(ievent);

    track_pt->GetEntry(ievent);
    track_eta->GetEntry(ievent);
    track_phi->GetEntry(ievent);

    b_JetAK10Trim_PT->GetEntry(ievent);
    b_JetAK10Trim_Eta->GetEntry(ievent);
    b_JetAK10Trim_Phi->GetEntry(ievent);
    b_JetAK10Trim_Mass->GetEntry(ievent);

    b_JetAK2Track_PT->GetEntry(ievent);
    b_JetAK2Track_Eta->GetEntry(ievent);
    b_JetAK2Track_Phi->GetEntry(ievent);
    b_JetAK2Track_Mass->GetEntry(ievent);
    b_JetAK2Track_BTag->GetEntry(ievent);

    particle_pid->GetEntry(ievent);
    particle_status->GetEntry(ievent);
    particle_pt->GetEntry(ievent);
    particle_eta->GetEntry(ievent);
    particle_phi->GetEntry(ievent);
    particle_mass->GetEntry(ievent);
    particle_m1->GetEntry(ievent);
    particle_m2->GetEntry(ievent);
    particle_d1->GetEntry(ievent);
    particle_d2->GetEntry(ievent);


    //
    // build table for WW, HH, WH 
    //
    wcount=0, hcount=0;
    vector<int> WhIndex;
    vector<int> ZhIndex;
    vector<int> TWindex;
    vector<int> Thindex;
    // loop over address of Wh and Zh particles          
    for(int i=0; i<nParticle; i++){
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==8880024) WhIndex.push_back(i);
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==8880023) ZhIndex.push_back(i);
    }

    for(int i=0; i<nParticle; i++){
      if(Particle_Status[i]==3 && abs(Particle_PID[i])==24 && WhIndex.size()!=0){
        for(int j=0; j<WhIndex.size(); j++){
          if(Particle_M1[i]==WhIndex[j] || Particle_M2[i]==WhIndex[j]) {
            wcount++;
            TWindex.push_back(i);
          }
        }
      }
      if(Particle_Status[i]==3 && Particle_PID[i]==25 && ZhIndex.size()!=0){
        for(int j=0; j<ZhIndex.size(); j++){
          if(Particle_M1[i]==ZhIndex[j] || Particle_M2[i]==ZhIndex[j]) {
            hcount++;
            Thindex.push_back(i);
          }
        }
      }
    }

    //
    //Begin Truth Level Study 
    //
    // if(wcount==1 && hcount==1){
      //Begin analysis
      
      //cut_flow[0]++;
      //cut_flow_boost[0]++;

      //
      //Identify hadronic W/Z/H
      //
      
    vector<TLorentzVector> hadronicWTLV;
    vector<TLorentzVector> hadronicHiggsTLV, hadronicWZTLV, hadronicWZHTLV;
    int index=0;
    for(int j=0;j<100;j++){
      int PID = abs(Particle_PID[j]);
      int  D1 = Particle_D1[j];
      int  D2 = Particle_D2[j];
      //Select W,Z,H only   
      if(PID <23 || PID>25) continue;
      //Daughter must exist   
      if(D1<0) continue;
      int D1_PID = abs(Particle_PID[D1]);
      //Daughter must be u,d,c,s,b  (Hadronic)  
      if(D1_PID>5) continue;
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(Particle_PT[j],Particle_Eta[j],Particle_Phi[j],Particle_Mass[j]);
      hadronicWZHTLV.push_back(tmp);
      
      if(PID == 24){
	hadronicWTLV.push_back(tmp);
      }
      if(PID == 25){
	hadronicHiggsTLV.push_back(tmp);
      }
      else{
	//23: Zboson  24: Wboson  
	hadronicWZTLV.push_back(tmp);
      }
    }
    totalH += hadronicHiggsTLV.size();
    totalW += hadronicWTLV.size();
    
    
    // Truth Parton Jets 
    vector<TLorentzVector> truthWoffspringTLV;
    vector<TLorentzVector> truthHoffspringTLV;
    for(int j=0;j<nParticle;j++){
      
      if(Particle_Status[j]!=3) continue;
      int PID = abs(Particle_PID[j]);
      int M1 = Particle_M1[j];
      int M1_PID = abs(Particle_PID[M1]);
      
      //cout << "Particle info: " << "  PID= " << PID << "   M1= " << M1 << "   M1_PID= " << M1_PID << endl; 	
      if(PID>5) continue;
      
      if(M1_PID==24) { //if W
	TLorentzVector tmpw;
	tmpw.SetPtEtaPhiM(Particle_PT[j],Particle_Eta[j],Particle_Phi[j],Particle_Mass[j]);
	truthWoffspringTLV.push_back(tmpw);
      }
      
      if(PID==5){
	if(M1_PID==25) { //if h
	  TLorentzVector tmph;
	  tmph.SetPtEtaPhiM(Particle_PT[j],Particle_Eta[j],Particle_Phi[j],Particle_Mass[j]);
	  truthHoffspringTLV.push_back(tmph);
	}
      }
    }    
    
    
    Double_t dRjj;
    Double_t dPhijj;
    Double_t dEtajj;
    TLorentzVector WTLV;
    for(int j=0; j<truthWoffspringTLV.size(); j++){
      if(truthWoffspringTLV.size()<2) continue;
      dRjj = truthWoffspringTLV[0].DeltaR(truthWoffspringTLV[1]);
      dPhijj =  truthWoffspringTLV[0].DeltaPhi(truthWoffspringTLV[1]); 
      dEtajj =  sqrt(dRjj*dRjj-dPhijj*dPhijj);
      WTLV = truthWoffspringTLV[0]+truthWoffspringTLV[1];
    }
    if(truthWoffspringTLV.size()>1){
      fillHistogram("Truth W mass",100,0,800,WTLV.M());
      fillHistogram("Truth dR(j,j)",50,0,4,dRjj);
      //fillHistogram("Truth dEta(j,j)",50,0,4,dEtajj);
      //fillHistogram("Truth lead q pt",100,0,1500, truthWoffspringTLV[0].Pt());
      //fillHistogram("Truth sublead q pt",100,0,1500,truthWoffspringTLV[1].Pt());
      //fillHistogram("Truth lead inv mass (q)",100,0,40,truthWoffspringTLV[0].M());
      //fillHistogram("Truth sublead inv mass (q)",100,0,40,truthWoffspringTLV[1].M());
    }
    
    Double_t dRbb;
    Double_t dPhibb;
    Double_t dEtabb;
    TLorentzVector HTLV;
    for(int j=0; j<truthHoffspringTLV.size(); j++){
      if(truthHoffspringTLV.size()<2) continue;
      dRbb = truthHoffspringTLV[0].DeltaR(truthHoffspringTLV[1]);
      dPhibb = truthHoffspringTLV[0].DeltaPhi(truthHoffspringTLV[1]);
      dEtabb = sqrt(dRbb*dRbb-dPhibb*dPhibb);
      HTLV = truthHoffspringTLV[0]+truthHoffspringTLV[1];
    }
    if(truthHoffspringTLV.size()>1){
      fillHistogram("Truth H mass",100,0,800,HTLV.M());
      fillHistogram("Truth dR(b,b)",50,0,4,dRbb);
      fillHistogram("Truth dEta(b,b)",50,0,4,dEtabb);
      //fillHistogram("Truth lead b pt",100,0,1500, truthHoffspringTLV[0].Pt());
      //fillHistogram("Truth sublead b pt",100,0,1500,truthHoffspringTLV[1].Pt());
      //fillHistogram("Truth lead inv mass (b)",100,0,40,truthHoffspringTLV[0].M());
      //fillHistogram("Truth sublead inv mass (b)",100,0,40,truthHoffspringTLV[1].M());
    }
    
    if(hadronicHiggsTLV.size()>0 && hadronicWTLV.size()>0){
      fillHistogram("Wpt",100,0,1500,hadronicHiggsTLV[0].Pt());
      fillHistogram("Hpt",100,0,1500, hadronicWTLV[0].Pt());
      fillHistogram("Hpt/Wpt",100,0,5, hadronicHiggsTLV[0].Pt()/hadronicWTLV[0].Pt());
      //fillHistogram("W Eta",50,-5,5,hadronicHiggsTLV[0].Eta());
      //fillHistogram("H Eta",50,-5,5, hadronicWTLV[0].Eta());
      fillHistogram("dR(W,H)",100,0,8,hadronicHiggsTLV[0].DeltaR(hadronicWTLV[0]));
    }
    
    
    //
    // Baseline
    //
    
    //MET
    if(nMET<1) continue;
    ETmiss = MET[0];
    fillHistogram("MET",100,0,2000,ETmiss);
    
    //0l + 4-5 jets  
    num_Lep = 0;
    num_JET = 0;
    num_bJET = 0;
    for(int j=0; j<nJET; j++){
      if(j>=MaxnJet) break;
      if(JET_PT[j] > 30. && fabs(JET_ETA[j]) < 2.4) num_JET++;
	if((JET_PT[j] > 30. && fabs(JET_ETA[j]) < 2.4) && JET_BTAG[j]>0) num_bJET++;
    }
    //fillHistogram("N-CalJet",10,0,10,num_JET);
    //fillHistogram("N-bJet",10,0,10,num_bJET);
    for(int j=0; j<nMUON; j++){
      if(j>=MaxnJet) break;
      if(MUON_PT[j] > 10. && fabs(MUON_ETA[j]) < 2.4) num_Lep++;
    }
    for(int j=0; j<nElec; j++){
      if(j>=MaxnJet) break;
      if(ELEC_PT[j] > 10. && fabs(ELEC_ETA[j]) < 2.5) num_Lep++;
    }
    
    //Track Veto
    iso_flag = 0;
    for(int j=0; j<nElec; j++){
      if(j>=MaxnELEC) break;
      if( ELEC_PT[j]<5. || ELEC_PT[j]/ELEC_SUMPT[j]<5. ) continue;
      Double_t deltaphi = fmod( fabs(MET_PHI[0] - ELEC_PHI[j]), TMath::TwoPi() );
      if(deltaphi > TMath::Pi()) deltaphi = TMath::TwoPi() - deltaphi;
      //if(deltaphi<0) deltaphi += 4.*TMath::Pi();
      //if(deltaphi>2.*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      //if(deltaphi>TMath::Pi()) deltaphi = 2.*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*ELEC_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
      }
    for(int j=0; j<nMUON; j++){    
      if(j>=MaxnMUON) break;
      if( MUON_PT[j]<5. || MUON_PT[j]/MUON_SUMPT[j]<5. ) continue;
      Double_t deltaphi = fmod( fabs(MET_PHI[0] - MUON_PHI[j]), TMath::TwoPi() );
      if(deltaphi > TMath::Pi()) deltaphi = TMath::TwoPi() - deltaphi;
      //if(deltaphi<0) deltaphi += 4.*TMath::Pi();
      //if(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      //if(deltaphi>TMath::Pi()) deltaphi = 2.*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*MUON_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }
    
    TVector3 trkpj, trkpk;
    for(int j=0; j< nTRACK; j++){
      if( TRACK_PT[j]<10. ) continue;
      //pT scalar sum of all tracks (pT>0.5GeV & |eta|<2.5) within dR=0.3 around the parent track
      Double_t TRACK_SUMPT=0;
      trkpj.SetPtEtaPhi(TRACK_PT[j], TRACK_ETA[j], TRACK_PHI[j]);
      for(int k=0; k<nTRACK; k++){
	if(k==j) continue;
	if(TRACK_PT[k]<0.5 || fabs(TRACK_ETA[k]) >2.5) continue;
	trkpk.SetPtEtaPhi(TRACK_PT[k], TRACK_ETA[k], TRACK_PHI[k]);
	
	if( trkpj.DeltaR(trkpk) < 0.3)  TRACK_SUMPT+=TRACK_PT[k];
      }
      //it's better to use a denominator which guarantee to be non-zero
      if(TRACK_SUMPT/TRACK_PT[j] > 0.1 ) continue;
      Double_t deltaphi = fmod( fabs(MET_PHI[0] - TRACK_PHI[j]), TMath::TwoPi() );
      if(deltaphi > TMath::Pi()) deltaphi = TMath::TwoPi() - deltaphi;
      //if(deltaphi<0) deltaphi += 4.*TMath::Pi();
      //if(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      //if(deltaphi>TMath::Pi()) deltaphi = 2.*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*TRACK_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }
    
    //Delta Phi
    Double_t deltaphi[4];
    Int_t jetIndex[4]; std::fill_n(jetIndex, 4, -1);
    Int_t integer = 0;
    for(int j=0; j<nJET; j++){                                          
      if(j>=MaxnJet) break;
      if( JET_PT[j]>30. && fabs(JET_ETA[j] )<2.4){
	jetIndex[integer] = j;
	integer++;
	if(integer == 4) break;
      }
    }
    for(int j=0; j<4; j++){
      if(jetIndex[j] <0) continue;
      deltaphi[j] = fmod( fabs(MET_PHI[0] - JET_PHI[jetIndex[j]]), TMath::TwoPi() );
      if(deltaphi[j] > TMath::Pi()) deltaphi[j] = TMath::TwoPi() - deltaphi[j];
      //if(deltaphi[j]<0) deltaphi[j] += 4.*TMath::Pi();
      //if(deltaphi[j]>2.*TMath::Pi()) deltaphi[j] -= 2.*TMath::Pi();
      //if(deltaphi[j]>TMath::Pi()) deltaphi[j] = 2.*TMath::Pi() - deltaphi[j];
    }
    dphi_flag=0;
    if( (deltaphi[0]>0 && deltaphi[0]<=0.5) ||
	(deltaphi[1]>0 && deltaphi[1]<=0.5) ||
	(deltaphi[2]>0 && deltaphi[2]<=0.3) ||
	(deltaphi[3]>0 && deltaphi[3]<=0.3))         dphi_flag=1;
    
    
    //
    // Resolved
    //
    
    //Jet Selection
    selectJetsLV.clear();
    TLorentzVector lvec;
    vector<double> selectJetsPT, selectJetsPhi;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
	lvec.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	selectJetsLV.push_back(lvec);
	selectJetsPT.push_back(JET_PT[j]);
	selectJetsPhi.push_back(JET_PHI[j]);
      }
    }
    
    //btagging
    Nb_M=0;
    Nb_L=0;
    TLorentzVector wlvec;
    vector<TLorentzVector> selectWJetsLV;
    vector<int> mediumBindex, looseBindex, otherIndex;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
	if( Bool_t(JET_BTAG[j] & (1<<1)) ) mediumBindex.push_back(j);
	else if( Bool_t(JET_BTAG[j] & (1<<0)) ) looseBindex.push_back(j);
	else if(JET_BTAG[j]==0){
	  otherIndex.push_back(j); // jet fails b-tagging
	  wlvec.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	  selectWJetsLV.push_back(wlvec); // collecting possible W-jets
	}
      }
    }
    Nb_M = mediumBindex.size();
    Nb_L = Nb_M + looseBindex.size();
    
    
    vector<int> bIndex;
    // search 2 calo-jet that passed medium b-tag 
    for(unsigned int k=0; k<mediumBindex.size(); k++)
      bIndex.push_back(mediumBindex[k]);
    
    for(unsigned int k=0; k<looseBindex.size(); k++)
      bIndex.push_back(looseBindex[k]);
    
    
    // higgs mass
    hmass=0;
    dRbb_flag = -1;
    TLorentzVector higgsCandidateLV;
    vector<TLorentzVector> higgsBucketTLV;
    if(bIndex.size()>=2){
      int idx0=bIndex[0];
      int idx1=bIndex[1];
      TLorentzVector b0vec, b1vec;
      b0vec.SetPtEtaPhiM(JET_PT[idx0],JET_ETA[idx0],JET_PHI[idx0],JET_M[idx0]);
      b1vec.SetPtEtaPhiM(JET_PT[idx1],JET_ETA[idx1],JET_PHI[idx1],JET_M[idx1]);
      higgsCandidateLV = b0vec+b1vec;
      hmass = higgsCandidateLV.M();
      higgsBucketTLV.push_back(b0vec);
      higgsBucketTLV.push_back(b1vec);
      Double_t dRbb = b0vec.DeltaR(b1vec);
      if(dRbb <= 1.0) dRbb_flag=1;
      //cout << "dRbb= " << dRbb << endl;                                               
      fillHistogram("higgs mass window",100,0,300,hmass);
      fillHistogram("dR(b,b)",40,0,4,dRbb);
    }
    
    // Count num. of reco b-tag match to truth b
    int btag_count=0;
    for(int i=0; i<truthHoffspringTLV.size(); i++){
      double bconeR = 0.4;
      bool bmatch = false;
      for(int j=0; j<higgsBucketTLV.size(); j++){
	if(higgsBucketTLV.size()!=2) continue;
	if(higgsBucketTLV[j].DeltaR(truthHoffspringTLV[i])<bconeR){
	  bconeR = higgsBucketTLV[j].DeltaR(truthHoffspringTLV[i]);
	  bmatch = true;
	}
      }
      if(bmatch==false)continue;
      btag_count++;
    }
    
    total_reco_Hbb += btag_count;
    
    // Counting num. truth b from H                                     
    if(truthHoffspringTLV.size()==2){
      total_truth_Hbb += truthHoffspringTLV.size();
    }
    
    
    // W mass
    // select two leading Jets that are not b-tagged
    Wmass=0;
    dRjj_flag = -1;
    TLorentzVector WCandidateLV;
    vector<TLorentzVector> WBucketTLV;
    if(otherIndex.size()>1){
      Double_t dRjj_min;
      if(otherIndex.size()==2){
	TLorentzVector q_vec[2];
	for(int i=0; i<2; i++){
	  int idx = otherIndex[i];
	  q_vec[i].SetPtEtaPhiM(JET_PT[idx],JET_ETA[idx],JET_PHI[idx],JET_M[idx]);
	}
	dRjj_min = q_vec[0].DeltaR(q_vec[1]);
	WCandidateLV = q_vec[0] + q_vec[1];
	Wmass = WCandidateLV.M();
	WBucketTLV.push_back(q_vec[0]);
	WBucketTLV.push_back(q_vec[1]);	
      }
      
      if(otherIndex.size()==3){
	TLorentzVector q_vec[3];
	TLorentzVector q1vec, q2vec;
	for(int i=0; i<3; i++){
	  int idx = otherIndex[i];
	  q_vec[i].SetPtEtaPhiM(JET_PT[idx],JET_ETA[idx],JET_PHI[idx],JET_M[idx]);
	}
	dRjj_min = q_vec[0].DeltaR(q_vec[1]);
	WCandidateLV = q_vec[0] + q_vec[1];
	Wmass = WCandidateLV.M();
	q1vec = q_vec[0];
	q2vec = q_vec[1];
	for(int j=0; j<2; j++){
	  Double_t dRjj_new = q_vec[j].DeltaR(q_vec[2]);
	  if(dRjj_new < dRjj_min) {
	    dRjj_min = dRjj_new;
	    WCandidateLV = q_vec[j] + q_vec[2];
	    Wmass = WCandidateLV.M();
	    q1vec = q_vec[j];
	    q2vec = q_vec[2];
	  }
	}
	WBucketTLV.push_back(q1vec);
	WBucketTLV.push_back(q2vec);
      }
      
      if(otherIndex.size()==4){
	TLorentzVector q_vec[4];
	TLorentzVector q1vec, q2vec;
	for(int i=0; i<4; i++){
	  int idx = otherIndex[i];
	  q_vec[i].SetPtEtaPhiM(JET_PT[idx],JET_ETA[idx],JET_PHI[idx],JET_M[idx]);
	}
	dRjj_min = q_vec[0].DeltaR(q_vec[1]);
	WCandidateLV = q_vec[0] + q_vec[1];
	Wmass = WCandidateLV.M();
	q1vec = q_vec[0];
	q2vec = q_vec[1];
	for(int j=0; j<3; j++){
	  for(int k=2; k<4; k++){
	    if(j==k) continue;
	    Double_t dRjj_new = q_vec[j].DeltaR(q_vec[k]);
	    if(dRjj_new <dRjj_min) {
	      dRjj_min = dRjj_new;
	      WCandidateLV = q_vec[j] + q_vec[k];
	      Wmass = WCandidateLV.M();
	      q1vec =q_vec[j];
	      q2vec =q_vec[k];
	    }
	  }
	}
	WBucketTLV.push_back(q1vec);
	WBucketTLV.push_back(q2vec);
      }
      
      if(dRjj_min <= 1.0) dRjj_flag=1;
      fillHistogram("W mass window",100,1,800,Wmass);
      fillHistogram("dR(j,j)",40,0,4,dRjj_min);
      fillHistogram("Num. non-b-tagg",8,0,8,otherIndex.size());
    }
    
    // Count num. of reco q-tag match to truth q form W                             
    int qtag_count=0;
    for(int i=0; i<truthWoffspringTLV.size(); i++){
      double qconeR = 0.4;
      bool qmatch = false;
      for(int j=0; j<WBucketTLV.size(); j++){
	if(WBucketTLV.size()!=2) continue;
	if(WBucketTLV[j].DeltaR(truthWoffspringTLV[i])<qconeR){
	  qconeR = WBucketTLV[j].DeltaR(truthWoffspringTLV[i]);
	  qmatch = true;
	}
      }
      if(qmatch==false)continue;
      qtag_count++;
    }
    
    total_reco_Wjj += qtag_count;
    
    // Counting num. truth q from W                                                 
    if(truthWoffspringTLV.size()==2){
      total_truth_Wjj += truthWoffspringTLV.size();
    }
    
      
      
    //
    // Boosted
    //
    
    // Small R trackjet                 
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
    
    // Large R jet                                 
    vector<TLorentzVector> CandFjetTLV;
    vector<int> CandFIndex;
    for(int j=0;j<nJetAK10Trim ;j++){
      if(fabs(JetAK10Trim_Eta[j])<2.0 && JetAK10Trim_Mass[j]>50){
	TLorentzVector tlvPrelim;
	tlvPrelim.SetPtEtaPhiM(JetAK10Trim_PT[j],JetAK10Trim_Eta[j],JetAK10Trim_Phi[j],JetAK10Trim_Mass[j]);
	CandFjetTLV.push_back(tlvPrelim);
	CandFIndex.push_back(j);
      }
    }
    
    //Constructed trim LR-jet 
    akt10trimTLV.clear();
    akt10trimBosonTag.clear();
    vector<int> getFjetIndex;
    for(int k=0; k<CandFjetTLV.size(); k++){
      if(CandFjetTLV[0].Pt()<=200) continue;
      if(CandFjetTLV[1].Pt()<=200) continue;
      akt10trimTLV.push_back(CandFjetTLV[k]);
      getFjetIndex.push_back(CandFIndex[k]);
      //fillHistogram("Cand. LRjet",5,0,5,akt10trimTLV.size());
      //fillHistogram("LRjet Mass distribution",100,0,800,akt10trimTLV[k].M());
    }  
    
    //Get Boson Tag 
    for(int k=0; k<akt10trimTLV.size(); k++){
      akt10trimBosonTag.push_back( getBosonTag(akt10trimTLV[k], hadronicWZTLV) );
      fillHistogram("Boson Tag", 3,0,3, getBosonTag(akt10trimTLV[k], hadronicWZTLV));
    }
    
    //W & h Candidate 
    Lhmass=0;
    LWmass=0;
    HiggsCandidateTLV.clear();
    BosonCandidateTLV.clear();
    if(akt10trimTLV.size()>1){
      if(akt10trimTLV[0].M()>akt10trimTLV[1].M()){
	HiggsCandidateTLV.push_back(akt10trimTLV[0]);
	BosonCandidateTLV.push_back(akt10trimTLV[1]);
	Windex = getFjetIndex[1];
      }
      else if(akt10trimTLV[1].M()>akt10trimTLV[0].M()){
	HiggsCandidateTLV.push_back(akt10trimTLV[1]);
	BosonCandidateTLV.push_back(akt10trimTLV[0]);
	Windex = getFjetIndex[0];
      }
      
      Lhmass = HiggsCandidateTLV[0].M();
      LWmass = BosonCandidateTLV[0].M();
      
      fillHistogram("Higgs mass(boosted)",100,0,800, HiggsCandidateTLV[0].M());
      fillHistogram("W mass(boosted)",100,0,500, BosonCandidateTLV[0].M());
      fillHistogram("H PT", 100,0,2000, HiggsCandidateTLV[0].Pt());
      fillHistogram("W PT", 100,0,2000, BosonCandidateTLV[0].Pt());
    }
    
    
    
    eventnum=ievent;
    
    //
    //Call Functions 
    //
    
    //Hybrid resolved cut binary intialization
    vector<Long64_t> resolvedCut;
    resolvedCut.clear();
    for(int j=0;j<13;j++){
      resolvedCut.push_back(0);
    }
    resolveCutSize = resolvedCut.size();
    
    vector<Long64_t> boostedCut;
    boostedCut.clear();
    for(int j=0;j<13;j++){
      boostedCut.push_back(0);
    }
    boostCutSize = boostedCut.size();
    
    //passBoost(&cut_flow_boost[0]);
    //passResolved(&cut_flow[0]);
    passResolved(&resolvedCut[0]);
    passBoost(&boostedCut[0]);
    passHybrid(&resolvedCutFlow[0], &boostedCutFlow[0], &resolvedCut[0], &boostedCut[0], &combinedCutFlow[0]);
    
    //}// End Truth Level Study  
  }// End of Events
  cout << "store cut info ..." << endl;
  
  
  cout << "totalW= " << totalW << "  totalWtag= " << totalWtag << endl;
  cout << "totalH= " << totalH << "  totalHtag= " << totalHtag << endl;
  cout << "getBosonTag -> 0 = " << zeroWtag << endl;
  cout << "getBosonTag -> 1 = "<< oneWtag << endl;
  cout << "truth Hbb= " << total_truth_Hbb << "    reco Hbb= " << total_reco_Hbb << "   "
       << "reoc/trut= "<< total_reco_Hbb/total_truth_Hbb << endl;
  cout << "truth Wjj= " << total_truth_Wjj << "    reco Wjj= " << total_reco_Wjj << "   "
       << "reoc/trut= "<< total_reco_Wjj/total_truth_Wjj << endl;

  //Resolved
  ofstream fileoutR("cutinfoR_wh_had.dat");
  TString cutname[13]={"No Selection: ", 
		       "0l + (4~5) Jets: ",
		       "Track Veto: ",
		       "dPhi_{1,2}>0.5,dPhi_{3,4}>0.3 : ",
		       "dRjj < 1.0: ",
		       "Wjet mass in [70,100]: ", 
		       "",
		       "Hjet mass in [100,140]: ",
		       "Nb-tag=2: ",
		       "MET > 150: ",
		       "MET > 200 : ",
		       "MET > 300 : ",
		       "MET > 450 : "};

  for(int j=0;j<13;j++) fileoutR << "[" << j << "]: " << cutname[j] << resolvedCutFlow[j] << endl;
  fileoutR.close();

  //Boosted
  ofstream fileoutB("cutinfoB_wh_had.dat");
  TString cutnameB[13]={"No Selection: ",
                        "0l + rm(4~5) Jets cond: ",
                        "Track Veto: ",
                        "dPhi_{1,2}>0.5,dPhi_{3,4}>0.3 : ",
                        "2 LR-jet: ",
                        "Wjet mass in [70,100]: ",
                        "Boson tagging on Wjet: ",
                        "Hjet mass in [100,140]: ",
                        "Nb-tag=2: ",
			"MET > 150 : ",
			"MET > 200 : ",
                        "MET > 300 : ",
                        "MET > 450 : "};

  for(int j=0;j<13;j++) fileoutB << "[" << j << "]: " << cutnameB[j] << boostedCutFlow[j] << endl;
  fileoutB.close();

  //Combine                
  ofstream fileoutC("cutinfoH_wh_had.dat");
  TString cutnameH[13]={"no selection", "c1","c2", "c3","c4",
			"c5","c6","c7","c8",
			"c9","c10","c11","c12"};
  for(int j=0;j<13;j++) fileoutC << cutnameH[j] << "[" << j << "]: "<< combinedCutFlow[j] <<endl;
  fileoutC.close();


  //Analyze data storage                          
  Double_t fvar = 1000 ;
  Double_t kvar = 1.0 ;
   
  // Resolved 
  ofstream fileoutAr("resolved_whhad_acc_HL.dat", ios::app);
  fileoutAr << fvar << " " << kvar << " " << resolvedCutFlow[9] << " " << resolvedCutFlow[10] << " " << resolvedCutFlow[11] << " " << resolvedCutFlow[12] << "" << endl;
  fileoutAr.close();

  // Boost
  ofstream fileoutAb("boosted_whhad_acc_HL.dat", ios::app);
  fileoutAb << fvar << " " << kvar << " " << boostedCutFlow[9] << " " << boostedCutFlow[10] << " " << boostedCutFlow[11] << " " << boostedCutFlow[12] << "" << endl;
  fileoutAb.close();

  // Hybrid 
  ofstream fileoutAh("hybrid_whhad_acc_HL.dat", ios::app);
  fileoutAh << fvar << " " << kvar << " " << combinedCutFlow[9] << " " << combinedCutFlow[10] << " " << combinedCutFlow[11] << " " << combinedCutFlow[12] << "" << endl;
  fileoutAh.close();
  

  //Save Histogram   
  int ipic = 0;
  int nhist = 45;
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


int getBosonTag(TLorentzVector reco, vector<TLorentzVector>& truth, float coneR){
  //D2 tagger 
  // http://cds.cern.ch/record/2258132/files/ATLAS-CONF-2017-018.pdf Sec5

  float Signal = 0.5;
  float QCD=0.02;
  bool isMatch= false;
  
  for(int i=0; i< truth.size(); i++){
    if(reco.DeltaR(truth[i])< coneR) isMatch = true;
    
  }

  float efficiency = QCD;
  if(isMatch==true) efficiency = Signal;
  

  int BitNumber = 0;
  int BosonTag = 0;
  // randGen.Uniform() generate distribution in [0,1] 
  // Operand:= (A <= B) << Bit; if () is true we have 1, else 0
  // BosonTag = BosonTag(default=0) | Operand. If 0 | 0 we have 0 else we have 1.
  BosonTag |= (randGen.Uniform() <= efficiency)  << BitNumber;

  return BosonTag;
}


bool passBoost(Long64_t* cutbinaryB){

  if(num_Lep>0) return false;
  //if(num_JET<4 || num_JET>5 || num_Lep>0) return false;
  cutbinaryB[1]=1;

  if(iso_flag==1) return false;
  cutbinaryB[2]=1;

  if(dphi_flag==1) return false;
  cutbinaryB[3]=1;

  if(akt10trimTLV.size()<1) return false;
  cutbinaryB[4]=1;

  //if(BosonCandidateTLV[0].M()<70 || BosonCandidateTLV[0].M()>100) return false;
  if(LWmass<=70 || LWmass>=100) return false;
  cutbinaryB[5]=1;

  if(akt10trimBosonTag[Windex] != 1) return false;
  cutbinaryB[6]=1;

  int NBtag_H = 0;
  for(int j=0;j<akt2trackTLV.size();j++){
    double minR=1.0;
    int match=-1;
    if(HiggsCandidateTLV.size()==0)continue;
    if(HiggsCandidateTLV[0].DeltaR(akt2trackTLV[j]) <minR){
      minR=HiggsCandidateTLV[0].DeltaR(akt2trackTLV[j]);
      match=1;
    }
    if(match==-1) continue;
    if(akt2trackBTag[j]==0) continue;
    NBtag_H++;
  }
  fillHistogram("N b-tag tkj",6,0,6,NBtag_H);

  if(Lhmass<=100 || Lhmass>=140) return false;
  //if(HiggsCandidateTLV[0].M()<=100 || HiggsCandidateTLV[0].M()>=140) return false;
  cutbinaryB[7]=1;

  if(NBtag_H!=2) return false;
  cutbinaryB[8]=1;

  //dummy
  if(ETmiss<=150) return false;
  cutbinaryB[9]=1;

  //dummy
  if(ETmiss<=200) return false;
  cutbinaryB[10]=1;

  if(ETmiss<=300) return false;
  cutbinaryB[11]=1;

  if(ETmiss<=450) return false;
  cutbinaryB[12]=1;

  return true;
}


bool passResolved(Long64_t* cutbinaryR){

  if(num_JET<4 || num_JET>5 || num_Lep>0) return false;
  cutbinaryR[1]=1;

  if(iso_flag==1) return false;
  cutbinaryR[2]=1;

  if(dphi_flag==1) return false;
  cutbinaryR[3]=1;

  if(dRjj_flag==-1) return false;
  cutbinaryR[4]=1;

  if(Wmass<=70 || Wmass>=100) return false;
  cutbinaryR[5]=1;

  //dummy
  cutbinaryR[6]=1;

  if(hmass<=100 || hmass>=140) return false;
  cutbinaryR[7]=1;

  if(Nb_M<1 || Nb_L!=2) return false;
  cutbinaryR[8]=1;

  if(ETmiss<=150) return false;
  cutbinaryR[9]=1;

  if(ETmiss<=200) return false;
  cutbinaryR[10]=1;

  if(ETmiss<=300) return false;
  cutbinaryR[11]=1;

  if(ETmiss<=450) return false;
  cutbinaryR[12]=1;

  return true;
}
  
bool passHybrid(Long64_t* cut_flow,Long64_t* cut_flow2,Long64_t* cutbinaryR,Long64_t* cutbinaryB, Long64_t* cut_flowH){

  cut_flow[0]++;
  for(int i=1;i<resolveCutSize;i++){
    if(cutbinaryR[i]==0) continue;
    cut_flow[i]++;
  }

  cut_flow2[0]++;
  for(int i=1;i<boostCutSize;i++){
    if(cutbinaryB[i]==0) continue;
    cut_flow2[i]++;
  }

  cut_flowH[0]++;
  for(int i=1;i<resolveCutSize;i++){
    if(cutbinaryR[i]==0 && cutbinaryB[i]==0)continue;
    cut_flowH[i]++;
  }


  return true;
}

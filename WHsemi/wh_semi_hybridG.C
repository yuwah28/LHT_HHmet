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

bool passBoost(Long64_t* cutbinaryB);
bool passResolved(Long64_t* cutbinaryR);
bool passHybrid(Long64_t* cutflow, Long64_t* cutflow2, Long64_t* cutbinaryR, Long64_t* cutbinaryB, Long64_t* cutflowH);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;


//global variables:  
Int_t eventnum;
Int_t resolveCutSize;
Int_t boostCutSize;

//Resolved                                                                       
Int_t num_Lep;
Int_t Nb_M;
Int_t Nb_L;
Int_t iso_track;
Double_t hmass;
Double_t ETmiss;
Double_t MassT;
Double_t MassCT;
vector<TLorentzVector> selectJetsLV;

//Boost
vector<TLorentzVector> akt10trimTLV;
vector<TLorentzVector> akt2trackTLV;
vector<UInt_t> akt2trackBTag;


int wh_semi_hybridG(){

  TChain *t = new TChain("Delphes");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/WHmetReso/Validation/4processSample_ffnumk5/ffnumk5_Hybrid_WH_DATA/Hdihiggs*.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/HHmet/Validation/13processHHSample_ffnum_kknum/WH_Hybrid_ffnum/dihiggs*.root");

  t->Add("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/\
ffnum/CHNAMEkknum/WHhybr.root");

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

  /*
  const Int_t     kMaxParticle = 5000;
  Int_t           nParticle = -1;
  Int_t           Particle_PID[kMaxParticle];  
  Int_t           Particle_Status[kMaxParticle];  
  Float_t         Particle_PT[kMaxParticle];  
  Float_t         Particle_Eta[kMaxParticle];  
  Float_t         Particle_Phi[kMaxParticle];  
  Float_t         Particle_Mass[kMaxParticle];  
  TBranch        *particle_pid;  
  TBranch        *particle_status;
  TBranch        *particle_pt;  
  TBranch        *particle_eta;   
  TBranch        *particle_phi;   
  TBranch        *particle_mass;  
  t->SetBranchAddress("Particle", &nParticle);
  t->SetBranchAddress("Particle.PID", Particle_PID, &particle_pid);
  t->SetBranchAddress("Particle.Status", Particle_Status, &particle_status);
  t->SetBranchAddress("Particle.PT", Particle_PT, &particle_pt);
  t->SetBranchAddress("Particle.Eta", Particle_Eta, &particle_eta);
  t->SetBranchAddress("Particle.Phi", Particle_Phi, &particle_phi);
  t->SetBranchAddress("Particle.Mass", Particle_Mass, &particle_mass);
  */
  const          Int_t MaxnELEC = 1000;
  Int_t          nElec = -1;
  Int_t          ELEC_CHARGE[MaxnELEC];
  Float_t        ELEC_PT[MaxnELEC];
  Float_t        ELEC_ETA[MaxnELEC];
  Float_t        ELEC_PHI[MaxnELEC];
  Float_t        ELEC_SUMPT[MaxnELEC];
  TBranch        *e_charge;
  TBranch        *e_pt;
  TBranch        *e_eta;
  TBranch        *e_phi;
  TBranch        *e_sumpt;
  t->SetBranchAddress("Electron", &nElec);
  t->SetBranchAddress("Electron.Charge", ELEC_CHARGE, &e_charge);
  t->SetBranchAddress("Electron.PT", ELEC_PT, &e_pt);
  t->SetBranchAddress("Electron.Eta", ELEC_ETA, &e_eta);
  t->SetBranchAddress("Electron.Phi", ELEC_PHI, &e_phi);
  t->SetBranchAddress("Electron.SumPtCharged", ELEC_SUMPT, &e_sumpt);

  const          Int_t MaxnMUON = 1000;
  Int_t          nMUON = -1;
  Int_t          MUON_CHARGE[MaxnMUON];
  Float_t        MUON_PT[MaxnMUON];
  Float_t        MUON_ETA[MaxnMUON];
  Float_t        MUON_PHI[MaxnMUON];
  Float_t        MUON_SUMPT[MaxnMUON];
  TBranch        *muon_charge;
  TBranch        *muon_pt;
  TBranch        *muon_eta ;
  TBranch        *muon_phi;
  TBranch        *muon_sumpt;
  t->SetBranchAddress("Muon", &nMUON);
  t->SetBranchAddress("Muon.Charge",MUON_CHARGE, &muon_charge);
  t->SetBranchAddress("Muon.PT", MUON_PT, &muon_pt);
  t->SetBranchAddress("Muon.Eta", MUON_ETA, &muon_eta);
  t->SetBranchAddress("Muon.Phi", MUON_PHI, &muon_phi);
  t->SetBranchAddress("Muon.SumPtCharged", MUON_SUMPT, &muon_sumpt);  

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


  ///////// Start Cut Analysis //////////
  Long64_t resolvedCutFlow[20]; std::fill_n(resolvedCutFlow, 20, 0);
  Long64_t boostedCutFlow[20]; std::fill_n(boostedCutFlow, 20, 0);
  Long64_t combinedCutFlow[20]; std::fill_n(combinedCutFlow, 20, 0); 

  //Long64_t cut_flow[20]; std::fill_n(cut_flow, 20, 0);
  //Long64_t cut_flow_boost[20]; std::fill_n(cut_flow_boost, 20, 0);
  //nevent=200000;
  int NDivision = nevent/10;
  for(Long64_t jevent = 0; jevent < nevent; jevent++){
    
    if(jevent%NDivision ==0) cout << jevent << "/" << nevent <<endl;
    //cut_flow[0]++;
    //cut_flow_boost[0]++;
    
    Long64_t ievent = t->LoadTree(jevent);

    met->GetEntry(ievent);
    met_phi->GetEntry(ievent);
    scalar_ht->GetEntry(ievent);
    jet_pt->GetEntry(ievent);
    jet_phi->GetEntry(ievent);
    jet_eta->GetEntry(ievent);
    jet_m->GetEntry(ievent);
    jet_btag->GetEntry(ievent);
    jet_flavor->GetEntry(ievent);
    /*
    particle_pid->GetEntry(ievent);
    particle_status->GetEntry(ievent);
    particle_pt->GetEntry(ievent);
    particle_eta->GetEntry(ievent);
    particle_phi->GetEntry(ievent);
    particle_mass->GetEntry(ievent);
    */
    e_charge->GetEntry(ievent);
    e_pt->GetEntry(ievent);
    e_eta->GetEntry(ievent);
    e_phi->GetEntry(ievent);
    e_sumpt->GetEntry(ievent);
    muon_charge->GetEntry(ievent);
    muon_pt->GetEntry(ievent);
    muon_eta->GetEntry(ievent);
    muon_phi->GetEntry(ievent);
    muon_sumpt->GetEntry(ievent);
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


    //Lepton selection
    num_Lep=0;
    vector<int> mIndex,mIndex2;
    vector<int> eIndex,eIndex2;
    for(int i=0; i<nMUON; i++){
      if(MUON_PT[i]>25 && fabs(MUON_ETA[i])<2.1) mIndex.push_back(i);
    }
    if(mIndex.size()==1){
      for(int i=0; i<nMUON; i++){
	if(i!=mIndex[0] && (MUON_PT[i]>5 && MUON_PT[i]<25) && MUON_SUMPT[i]<5) mIndex2.push_back(i);
      }
    }
    for(int i=0; i<nElec; i++){
      if(ELEC_PT[i]>30 && fabs(ELEC_ETA[i])<1.44) eIndex.push_back(i);
    }
    if(eIndex.size()==1){
      for(int i=0; i<nElec; i++){
        if(i!=eIndex[0] && (ELEC_PT[i]>5 && ELEC_PT[i]<30) && ELEC_SUMPT[i]<5) eIndex2.push_back(i);
      }
    }
    
    num_Lep = mIndex.size()+eIndex.size()+mIndex2.size()+eIndex2.size();
    //cout << "lepton number =  " << num_Lep << endl;
    Int_t lepIndex=-1;
    Int_t lepCharge=0; //if the lepton candidate presence charge is either +1 or -1
    Double_t lepPT, lepPhi;
    if(mIndex.size()==1 && mIndex2.size()==0 && eIndex.size()==0) {
      lepIndex = mIndex[0];
      lepPT = MUON_PT[mIndex[0]];
      lepPhi = MUON_PHI[mIndex[0]];
      lepCharge = MUON_CHARGE[mIndex[0]];
    }
    else if (eIndex.size()==1 && eIndex2.size()==0 && mIndex.size()==0) {
      lepIndex = eIndex[0];
      lepPT = ELEC_PT[eIndex[0]];
      lepPhi = ELEC_PHI[eIndex[0]];
      lepCharge= ELEC_CHARGE[eIndex[0]];
    }

    //Transverse Mass
    if( (mIndex.size()==1 && mIndex2.size()==0 && eIndex.size()==0) || (eIndex.size()==1 && eIndex2.size()==0 && mIndex.size()==0)){
      Double_t dPhi = fmod( fabs(MET_PHI[0]-lepPhi), TMath::TwoPi() );
      if(dPhi>TMath::Pi()) dPhi = TMath::TwoPi()-dPhi;
      MassT = TMath::Sqrt(2.*lepPT*MET[0]*(1.-TMath::Cos(dPhi)));
      //fillHistogram("MassT",100,0,1000,MassT);
    }
    
    //Track Veto
    iso_track = 0;
    TVector3 trkpj, trkpk;
    for(int j=0; j<nTRACK; j++){
      if( TRACK_PT[j]<=10. ) continue;
      Double_t TRACK_SUMPT=0;  
  
      trkpj.SetPtEtaPhi(TRACK_PT[j], TRACK_ETA[j], TRACK_PHI[j]);
    
      for(int k=0; k<nTRACK; k++){
        if(k==j) continue;
        if(TRACK_PT[k]<=10. || fabs(TRACK_ETA[k])>=2.4) continue;
        trkpk.SetPtEtaPhi(TRACK_PT[k], TRACK_ETA[k], TRACK_PHI[k]);
        if( trkpj.DeltaR(trkpk) < 0.3) TRACK_SUMPT += TRACK_PT[k];
      }

      if(TRACK_SUMPT<6 && TRACK_SUMPT<0.1*TRACK_PT[j] && TRACK_CHR[j]*lepCharge==-1) iso_track=1;
    }
  
    /*
    for(int j=0; j<nTRACK; j++){
      if(TRACK_PT[j]>10. && fabs(TRACK_ETA[j])<2.4) {
	fillHistogram("Track PT",100,0,800,TRACK_PT[j]);
	fillHistogram("Track ETA",50,-4,4,TRACK_ETA[j]);
      }
    }
    */

    //Jet selection
    selectJetsLV.clear();
    TLorentzVector lvec;
    vector<double> selectJetsPT, selectJetsPhi;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
	lvec.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	selectJetsLV.push_back(lvec);
      }
    }
    
    Nb_M=0;
    Nb_L=0;
    vector<int> mediumBindex, looseBindex, otherIndex;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
        if( Bool_t(JET_BTAG[j] & (1<<1)) ) mediumBindex.push_back(j);
        else if( Bool_t(JET_BTAG[j] & (1<<0)) ) looseBindex.push_back(j);
        else otherIndex.push_back(j);
      }
    }
    Nb_M = mediumBindex.size();
    Nb_L = Nb_M + looseBindex.size();
    
    
    vector<int> bIndex;
    if(selectJetsLV.size()==2){
      for(unsigned int k=0; k<mediumBindex.size(); k++)
	bIndex.push_back(mediumBindex[k]);
    
      if(bIndex.size()==1){
	for(unsigned int k=0; k<looseBindex.size(); k++)
	  bIndex.push_back(looseBindex[k]);
      }
    }
    

    //higgs candidate mass
    hmass=0;
    TLorentzVector higgsCandidateLV;
    if(bIndex.size()==2){
      int idx0 = bIndex[0];
      int idx1 = bIndex[1];
      TLorentzVector b0vec, b1vec;
      b0vec.SetPtEtaPhiM(JET_PT[idx0],JET_ETA[idx0],JET_PHI[idx0],JET_M[idx0]);
      b1vec.SetPtEtaPhiM(JET_PT[idx1],JET_ETA[idx1],JET_PHI[idx1],JET_M[idx1]);

      higgsCandidateLV = b0vec+b1vec;
      hmass = higgsCandidateLV.M();
      fillHistogram("higgs mass-Resolved",100,0,800,hmass);
      fillHistogram("Leading pT",100,0,1500, selectJetsLV[bIndex[0]].Pt());
      fillHistogram("Subleading pT",100,0,1500, selectJetsLV[bIndex[1]].Pt());
        
      //CoTransver Mass 
      Double_t dPhibb = fmod( fabs(JET_PHI[idx0]-JET_PHI[idx1]), TMath::TwoPi() );
      if(dPhibb>TMath::Pi()) dPhibb = TMath::TwoPi()-dPhibb;

      MassCT = TMath::Sqrt(2.*JET_PT[idx0]*JET_PT[idx1]*(1.+TMath::Cos(dPhibb))); 
      //fillHistogram("MassCT",100,0,1000,MassCT);
    }    
    
    //MET
    if(nMET<1) continue;
    ETmiss = MET[0];
    fillHistogram("MET",100,0,1000,MET[0]);

    //Boosted Jet
    //LRjet
    akt10trimTLV.clear();
    for(int j=0;j<nJetAK10Trim ;j++){
      if(JetAK10Trim_PT[j]>200 && fabs(JetAK10Trim_Eta[j])<2.5){
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(JetAK10Trim_PT[j],JetAK10Trim_Eta[j],JetAK10Trim_Phi[j],JetAK10Trim_Mass[j]);
        akt10trimTLV.push_back(tlv);                                                                         
        fillHistogram("akt10M",100,0,300,tlv.M());
      }
    }
    
    //Small TrackJet
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


    eventnum = ievent;
    //Call Functions 
    //Hybrid resolved cut binary intialization 
    vector<Long64_t> resolvedCut;
    resolvedCut.clear();
    for(int j=0;j<14;j++){
      resolvedCut.push_back(0);
    }
    resolveCutSize = resolvedCut.size();

    vector<Long64_t> boostedCut;
    boostedCut.clear();
    for(int j=0;j<14;j++){
      boostedCut.push_back(0);
    }
    boostCutSize = boostedCut.size();

    passResolved(&resolvedCut[0]);
    passBoost(&boostedCut[0]);
    passHybrid(&resolvedCutFlow[0], &boostedCutFlow[0], &resolvedCut[0], &boostedCut[0], &combinedCutFlow[0]);


    //passBoost(&cut_flow_boost[0]);
    //passResolved(&cut_flow[0]);
    
  }// end of events

 
  cout << "store cut info ..." << endl;

  //Resolved
  ofstream fileoutR("cutinfoR_local.dat");
  TString cutnameR[12]={"no selection", 
		       "N(lepton)=1", 
		       "Isolated Track Veto",
		       "",
		       "N(b-tags)=2", 
		       "Mbbar in [90,150]", 
		       "MET>125",
		       "MT>150", 
		       "MCT>170", 
		       "MET>200",
		       "",
		       ""};
  for(int j=0;j<12;j++) fileoutR << "[" << j << "] " << cutnameR[j] << ": " << resolvedCutFlow[j] <<endl;
  fileoutR.close();
  
  //Boosted
  ofstream fileoutB("cutinfoB_local.dat");
  TString cutname[12]={"no selection", 
		       "N(lepton)=1", 
		       "Isolated Track Veto",
		       "large R-jets=1",
		       "Nb-tags Leading Large R-jet=2", 
		       "Higgs Cand. Mass in [90,150]", 
		       "MET>125",
		       "MT>150",
		       "", 
		       "MET>200", 
		       "MET>300", 
		       "MET>450"};
  for(int j=0;j<12;j++) fileoutB << "[" << j << "] " << cutname[j] << ": "<< boostedCutFlow[j] << endl;
  fileoutB.close();

  //Combine               
  ofstream fileoutC("cutinfoC_3process.dat");
  TString cutnameH[12]={"no selection", "c1","c2", "c3","c4",
			"c5","c6","c7","c8",
                        "c9","c10","c11"};
  for(int j=0;j<12;j++) fileoutC << cutnameH[j] << "[" << j << "]: "<< combinedCutFlow[j] << endl;
  fileoutC.close();

  //Analyze data storage
  Double_t fvar = fnum ;
  Double_t kvar = nknum ;
  
  // Resolved 
  ofstream fileoutAr("resolved_whsemi_acc_eff.dat", ios::app);
  fileoutAr << fvar << " " << kvar << " " << resolvedCutFlow[8] << " " << resolvedCutFlow[9] 
	    << "" << endl;
  fileoutAr.close();

  // Boost
  ofstream fileoutAb("boosted_whsemi_acc_eff.dat", ios::app);
  fileoutAb << fvar << " " << kvar << " " << boostedCutFlow[8] << " " << boostedCutFlow[9] 
	    << " " << boostedCutFlow[10] << " " << boostedCutFlow[11] << "" << endl;
  fileoutAb.close();

  // Hybrid  
  ofstream fileoutAh("hybrid_whsemi_acc_eff.dat", ios::app);
  fileoutAh << fvar << " " << kvar << " " << combinedCutFlow[8] << " " << combinedCutFlow[9] 
	    << " " << combinedCutFlow[10] << " " << combinedCutFlow[11] << "" << endl;
  fileoutAh.close();

  
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

bool passBoost (Long64_t* cutbinaryB){

  if(num_Lep!=1) return false;
  cutbinaryB[1]=1;

  if(iso_track==1) return false;
  cutbinaryB[2]=1;

  //1 LRjet
  if(akt10trimTLV.size()!=1) return false;
  cutbinaryB[3]=1;
 
  //b-tagging
  Double_t higgsMass = 0;
  int Nbtag_track = 0;
  if(akt10trimTLV.size()==1){
    for(int j=0; j<akt2trackTLV.size(); j++){
      double minR = 1.0;
      int match = -1;

      if(akt10trimTLV[0].DeltaR(akt2trackTLV[j]) <minR){
	minR=akt10trimTLV[0].DeltaR(akt2trackTLV[j]);
	match=1;        
      }
      if(match==-1) continue;
      if(akt2trackBTag[j]>0) {
	Nbtag_track++;
	higgsMass = akt10trimTLV[0].M();
	fillHistogram("higgs mass Boosted",100,0,800,higgsMass);
	fillHistogram("higgs pT",100,0,1500, akt10trimTLV[0].Pt());
      }
    }
  }
  //2 btag LRjet
  if(Nbtag_track!=2 ) return false;
  cutbinaryB[4]=1;
  
  if(higgsMass<100 || higgsMass>140) return false;
  cutbinaryB[5]=1;

  if(ETmiss<=125) return false;
  cutbinaryB[6]=1;

  if(MassT<=150) return false;
  cutbinaryB[7]=1;

  //dummy
  cutbinaryB[8]=1;

  if(ETmiss<=200) return false;
  cutbinaryB[9]=1;

  if(ETmiss<=300) return false;
  cutbinaryB[10]=1;

  if(ETmiss<=450) return false;
  cutbinaryB[11]=1;

  return true;
}

bool passResolved(Long64_t* cutbinaryR){

  if(num_Lep!=1) return false;
  cutbinaryR[1]=1;

  if(iso_track==1) return false;
  cutbinaryR[2]=1;

  //dummy
  cutbinaryR[3]=1;

  if(selectJetsLV.size()!=2 || Nb_M<1 || Nb_L<2) return false;
  cutbinaryR[4]=1;

  if(hmass<=90 || hmass>=150) return false;
  cutbinaryR[5]=1;

  if(ETmiss<=125) return false;
  cutbinaryR[6]=1;
  
  if(MassT<=150) return false;
  cutbinaryR[7]=1;

  if(MassCT<=170) return false;
  cutbinaryR[8]=1;

  if(ETmiss<=200) return false;
  cutbinaryR[9]=1;

  if(ETmiss<=300) return false;
  cutbinaryR[10]=1;

  if(ETmiss<=450) return false;
  cutbinaryR[11]=1;

  //dummy
  //cutbinaryR[10]=1;
  
  //dummy
  //cutbinaryR[11]=1;


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

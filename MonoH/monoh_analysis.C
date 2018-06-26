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
vector<TLorentzVector> forwardRjet;
vector<double> akt2trackPhi;
//vector<double> dR_SRbjLRj;
vector<double> leadingLargeRjetPt;
vector<double> dissociate_CentralJet_Pt_sum;
vector<UInt_t> akt2trackBTag;
Int_t iso_dR;
Int_t largeRjet_counter;
//Double_t dissociate_CentralJet_Pt_sum;

//Resolved
Int_t num_Lep;
Int_t num_JET;
Int_t num_btagSmallJet;
Int_t NRjet_central;
Int_t btagCentralCounter;
Double_t ETmiss ;
Int_t iso_flag;
Double_t dMmin;
Double_t dRmax;
Double_t PTmiss;
Double_t dPhiSmallRMin;
Double_t dPhiPtMiss;
Double_t higgsCandidateJet_leadPT;
Double_t delta_Rjhjh;
Double_t dPhihiggs;
Double_t SCALAR_HT;
Double_t HTcentral2;
Double_t HTcentral3;
Double_t HT_CandidateJets_extraJet;

bool passBoost(Long64_t* cutflow);
bool passSemiBoost(Long64_t* cutflow);
bool passResolved(Long64_t* cutflow);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;



int monoh_analysis_ver3(){  

  TFile *fin=new TFile("/phys/groups/tev/scratch4/users/paklim/LHT_Data_Backup/monoH/puremonoh_myatlas.root","READONLY");
  //TFile *fin=new TFile("/phys/groups/tev/scratch4/users/paklim/LHT_Data_Backup/save_scan_jets_met_13TeV/f1000/jmetk05/monoh.root","READONLY"); 

  //histogram
  TFile *fout=new TFile("signal3.root","RECREATE");
  Float_t bins[]={150,180,280,420,650};
  TH1F* tt_MET = new TH1F("tt_MET","tt",4,bins);
  //TH1F* Zj_MET = new TH1F("Zj_MET","Zjets",4,bins);

  TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  cout << "Total Event = " << nevent << endl;
  t->SetMakeClass(1);
  //nevent=10000;
  //----------------SETUP BRANCHES----------------------- 
  const Int_t    MaxnMET = 5000;
  Int_t          nMET = -1;
  Float_t        MET[MaxnMET];
  Float_t        MET_PHI[MaxnMET];
  TBranch        *met;
  TBranch        *met_phi;  
  t->SetBranchAddress("MissingET", &nMET);
  t->SetBranchAddress("MissingET.MET", MET, &met);
  t->SetBranchAddress("MissingET.Phi", MET_PHI, &met_phi);

  const Int_t    MaxnHT = 5000;
  Int_t          nHT = -1;
  Float_t        ScalarHT[MaxnHT];
  TBranch        *scalar_ht;
  t->SetBranchAddress("ScalarHT", &nHT);
  t->SetBranchAddress("ScalarHT.HT", ScalarHT, &scalar_ht);

  const Int_t    MaxnJet = 10000;
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

  const Int_t     kMaxParticle = 10000;
  Int_t           nParticle = -1;
  Int_t           Particle_PID[kMaxParticle];   //[Particle_]
  Int_t           Particle_Status[kMaxParticle];   //[Particle_]
  Float_t         Particle_PT[kMaxParticle];   //[Particle_]
  Float_t         Particle_Eta[kMaxParticle];   //[Particle_]
  Float_t         Particle_Phi[kMaxParticle];   //[Particle_]
  Float_t         Particle_Mass[kMaxParticle];   //[Particle_]
  TBranch        *b_Particle_PID;   //!
  TBranch        *b_Particle_Status;   //!
  TBranch        *b_Particle_PT;   //!
  TBranch        *b_Particle_Eta;   //!
  TBranch        *b_Particle_Phi;   //!
  TBranch        *b_Particle_Mass;   //!
  t->SetBranchAddress("Particle", &nParticle);
  t->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
  t->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
  t->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
  t->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
  t->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
  t->SetBranchAddress("Particle.Mass", Particle_Mass, &b_Particle_Mass);
  
  const Int_t     kMaxJetAK10Trim = 10000;
  Int_t           nJetAK10Trim = nJetAK10Trim;
  Float_t         JetAK10Trim_PT[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Eta[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Phi[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Mass[kMaxJetAK10Trim];   //[JetAK10Trim_]
  TBranch        *b_JetAK10Trim_PT;   //!
  TBranch        *b_JetAK10Trim_Eta;   //!
  TBranch        *b_JetAK10Trim_Phi;   //!
  TBranch        *b_JetAK10Trim_Mass;   //!
  t->SetBranchAddress("JetAK10Trim", &nJetAK10Trim);
  t->SetBranchAddress("JetAK10Trim.PT", JetAK10Trim_PT, &b_JetAK10Trim_PT);
  t->SetBranchAddress("JetAK10Trim.Eta", JetAK10Trim_Eta, &b_JetAK10Trim_Eta);
  t->SetBranchAddress("JetAK10Trim.Phi", JetAK10Trim_Phi, &b_JetAK10Trim_Phi);
  t->SetBranchAddress("JetAK10Trim.Mass", JetAK10Trim_Mass, &b_JetAK10Trim_Mass);
  
  const Int_t     kMaxJetAK2Track = 10000; 
  Int_t           nJetAK2Track = nJetAK2Track;
  Float_t         JetAK2Track_PT[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Eta[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Phi[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_T[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Mass[kMaxJetAK2Track];   //[JetAK2Track_]
  UInt_t          JetAK2Track_BTag[kMaxJetAK2Track];   //[JetAK2Track_]
  Int_t           JetAK2Track_Flavor[kMaxJetAK2Track];
  TBranch        *b_JetAK2Track_PT;   //!
  TBranch        *b_JetAK2Track_Eta;   //!
  TBranch        *b_JetAK2Track_Phi;   //!
  TBranch        *b_JetAK2Track_Mass;   //!
  TBranch        *b_JetAK2Track_BTag;   //!
  TBranch        *b_JetAK2Track_Flavor;
  t->SetBranchAddress("JetAK2Track", &nJetAK2Track);
  t->SetBranchAddress("JetAK2Track.PT", JetAK2Track_PT, &b_JetAK2Track_PT);
  t->SetBranchAddress("JetAK2Track.Eta", JetAK2Track_Eta, &b_JetAK2Track_Eta);
  t->SetBranchAddress("JetAK2Track.Phi", JetAK2Track_Phi, &b_JetAK2Track_Phi);
  t->SetBranchAddress("JetAK2Track.Mass", JetAK2Track_Mass, &b_JetAK2Track_Mass);
  t->SetBranchAddress("JetAK2Track.BTag", JetAK2Track_BTag, &b_JetAK2Track_BTag);
  t->SetBranchAddress("JetAK2Track.Flavor", JetAK2Track_Flavor, &b_JetAK2Track_Flavor);

  const          Int_t MaxnELEC = 5000;
  Int_t          nElec = -1;
  Float_t        ELEC_PT[MaxnELEC];
  Float_t        ELEC_ETA[MaxnELEC];
  Float_t        ELEC_PHI[MaxnELEC];
  Float_t        ELEC_SUMPT[MaxnELEC];
  TBranch        *e_pt;
  TBranch        *e_eta;
  TBranch        *e_phi;
  TBranch        *e_sumpt;
  t->SetBranchAddress("Electron", &nElec);
  t->SetBranchAddress("Electron.PT", ELEC_PT, &e_pt);
  t->SetBranchAddress("Electron.Eta", ELEC_ETA, &e_eta);
  t->SetBranchAddress("Electron.Phi", ELEC_PHI, &e_phi);
  t->SetBranchAddress("Electron.SumPtCharged", ELEC_SUMPT, &e_sumpt);

  const          Int_t MaxnMUON = 5000;
  Int_t          nMUON = -1;
  Float_t        MUON_PT[MaxnMUON];
  Float_t        MUON_ETA[MaxnMUON];
  Float_t        MUON_PHI[MaxnMUON];
  Float_t        MUON_SUMPT[MaxnMUON];
  TBranch        *muon_pt;
  TBranch        *muon_eta ;
  TBranch        *muon_phi;
  TBranch        *muon_sumpt;
  t->SetBranchAddress("Muon", &nMUON);
  t->SetBranchAddress("Muon.PT", MUON_PT, &muon_pt);
  t->SetBranchAddress("Muon.Eta", MUON_ETA, &muon_eta);
  t->SetBranchAddress("Muon.Phi", MUON_PHI, &muon_phi);
  t->SetBranchAddress("Muon.SumPtCharged", MUON_SUMPT, &muon_sumpt);
  
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

  Int_t myMETcount = 0;
  /////////////////////////////////
  // Start Cut Analysis          //
  /////////////////////////////////
  Long64_t cut_flow[20]; std::fill_n(cut_flow, 20, 0);
  Long64_t cut_flow_boost[20];  std::fill_n(cut_flow_boost, 20, 0);
  //nevent=1000;
  int NDivision = nevent/10;
  for(Long64_t ievent = 0; ievent < nevent; ievent++){

    if(ievent%NDivision ==0) cout << ievent << "/" << nevent <<endl; 
    cut_flow[0]++;
    cut_flow_boost[0]++; 

    met->GetEntry(ievent);
    met_phi->GetEntry(ievent);
    scalar_ht->GetEntry(ievent);    
    jet_pt->GetEntry(ievent);
    jet_phi->GetEntry(ievent);
    jet_eta->GetEntry(ievent);
    jet_m->GetEntry(ievent);
    jet_btag->GetEntry(ievent);
    jet_flavor->GetEntry(ievent);
    b_Particle_PID->GetEntry(ievent);
    b_Particle_Status->GetEntry(ievent);
    b_Particle_PT->GetEntry(ievent);
    b_Particle_Eta->GetEntry(ievent);
    b_Particle_Phi->GetEntry(ievent);
    b_Particle_Mass->GetEntry(ievent);

    b_JetAK10Trim_PT->GetEntry(ievent);
    b_JetAK10Trim_Eta->GetEntry(ievent);
    b_JetAK10Trim_Phi->GetEntry(ievent);
    b_JetAK10Trim_Mass->GetEntry(ievent);

    b_JetAK2Track_PT->GetEntry(ievent);
    b_JetAK2Track_Eta->GetEntry(ievent);
    b_JetAK2Track_Phi->GetEntry(ievent);
    b_JetAK2Track_Mass->GetEntry(ievent);
    b_JetAK2Track_BTag->GetEntry(ievent);

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

    //Truth study
    TLorentzVector genHiggsTLV[2];
    int index=0;
    for(int j=0;j<20;j++){
      if(Particle_PID[j] == 25){
	genHiggsTLV[index].SetPtEtaPhiM(Particle_PT[j],Particle_Eta[j],Particle_Phi[j],Particle_Mass[j]);
	index++;
      }
      if(index==2) break;
    }

    ////////// Resolved Region ///////////
    //////////////////       
    //lepton isolation          
    //////////////////       
    iso_flag = 0;
    for(int j=0; j<nMUON; j++){
      if(j>=MaxnMUON) break;
      if(MUON_PT[j]>7. && fabs(MUON_ETA[j])<2.7 && MUON_SUMPT[j]<0.15*MUON_PT[j]) iso_flag =1;
    }
    for(int j=0; j<nElec; j++){
      if(j>=MaxnELEC) break;
      if(ELEC_PT[j]>7. && fabs(ELEC_ETA[j])<2.47 && ELEC_SUMPT[j]<0.15*ELEC_PT[j]) iso_flag =1;
    }
    
    //////  
    // MET      
    //////          
    if(nMET<1) continue;
    ETmiss = MET[0];
    tt_MET->Fill(MET[0]);
    //Zj_MET->Fill(MET[0]); 
    if(MET[0]>150) { 
      myMETcount++;
      fillHistogram("My MET hist",4,150,800,MET[0]);
    }
    
    ///////////
    //Scalar HT
    ///////////
    SCALAR_HT = ScalarHT[0];
    
    /////////////////////////      
    //DeltaPhi on small Rjets       
    /////////////////////////       
    int smallJetCounter = 0;
    dPhiSmallRMin = TMath::TwoPi();
    for(unsigned int j=0; j< nJET; j++){
      if(JET_PT[j]>20 && fabs(JET_ETA[j])<2.5){
	smallJetCounter++;
	if(smallJetCounter>=3) break;
	Double_t dPhiSmallR = fmod( fabs(MET_PHI[0] - JET_PHI[j]), TMath::TwoPi() );
	if(dPhiSmallR > TMath::Pi()) dPhiSmallR = TMath::TwoPi() - dPhiSmallR;
	if(dPhiSmallR <  dPhiSmallRMin) dPhiSmallRMin = dPhiSmallR;
      }
    }
    
    /////////  
    //pT-miss  
    /////////
    Double_t sumPx = 0, sumPy = 0;
    for(int j=0; j <nTRACK; j++){
      if(TRACK_CHR[j]!=0 && TRACK_PT[j]>0.5 && abs(TRACK_ETA[j])<2.5){
	ROOT::Math::RhoEtaPhiVector vec_ct(TRACK_PT[j], TRACK_ETA[j], TRACK_PHI[j]);
        sumPx += vec_ct.X();
        sumPy += vec_ct.Y();
      }
    }
    TVector2 vec2PTmiss(-sumPx, -sumPy);
    Double_t PTmissPhi = vec2PTmiss.Phi();
    PTmiss = sumPx*sumPx + sumPy*sumPy;
    //DeltaPhi PTmiss and MET     
    dPhiPtMiss = fmod(fabs(MET_PHI[0] - PTmissPhi), TMath::TwoPi() );
    if(dPhiPtMiss > TMath::Pi()) dPhiPtMiss = TMath::TwoPi() - dPhiPtMiss;
    
    //////////////////////  
    //SmallRjet with b-tag  
    ////////////////////// 
    num_btagSmallJet = 0;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>20  && fabs(JET_ETA[j])<4.5 && Bool_t(JET_BTAG[j]&(1<<0)) ) num_btagSmallJet++;
    }

    /////////////
    //Central Jet
    /////////////
    vector<double> centralPt;
    vector<TLorentzVector> centralRjet;
    NRjet_central = 0;
    for(int j=0; j<nJET; j++){         
      if(JET_PT[j]>20 && fabs(JET_ETA[j])<2.5){
	NRjet_central++;
        TLorentzVector central_tlv;
        central_tlv.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
        centralRjet.push_back(central_tlv);
	centralPt.push_back(JET_PT[j]);
      }
    }
    //central bTag
    vector<double> btagCentralJet_Pt, blikeCentralJet_Pt, btagCentralJet_Phi, blikeCentralJet_Phi;
    vector<TLorentzVector> order_bTagCentral, order_blikeCentral;
    btagCentralCounter = 0;
    for(int j=0; j<nJET; j++){
      if(JET_PT[j]>20 && fabs(JET_ETA[j])<2.5 && Bool_t(JET_BTAG[j]&(1<<0))){
	btagCentralCounter++;        // counting btag central jet
	TLorentzVector btagCentralLV;
	btagCentralLV.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	order_bTagCentral.push_back(btagCentralLV);
	btagCentralJet_Pt.push_back(JET_PT[j]);
	btagCentralJet_Phi.push_back(JET_PHI[j]);
      }
      //central bLike candidate
      else if(JET_PT[j]>20 && fabs(JET_ETA[j])<2.5 && Bool_t(JET_BTAG[j]&(1<<0))==0){
	TLorentzVector blikeCentralLV;
	blikeCentralLV.SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	order_blikeCentral.push_back(blikeCentralLV);
	blikeCentralJet_Pt.push_back(JET_PT[j]);
	blikeCentralJet_Phi.push_back(JET_PHI[j]);	
      } 
    }
    
    /////////////////
    //Higgs candidate
    /////////////////
    //Double_t HT_CandidateJets_extraJet;
    TLorentzVector higgsCandidate;
    if(order_bTagCentral.size()>=2 && blikeCentralJet_Pt.size()>=0){
      // leading higg candidate jet
      higgsCandidateJet_leadPT = max(btagCentralJet_Pt[0],btagCentralJet_Pt[1]);
      fillHistogram("higgs Candidate Jet PT",100,0,1500,higgsCandidateJet_leadPT);  //histogram    

      // higgs reco
      higgsCandidate = order_bTagCentral[0]+order_bTagCentral[1];
      fillHistogram("reco higgs inv mass",100,0,800,higgsCandidate.M());

      // deltaR
      delta_Rjhjh=fabs(order_bTagCentral[0].DeltaR(order_bTagCentral[1]));
      
      // angle seperation
      Double_t higgsPhi = higgsCandidate.Phi(); 
      dPhihiggs = fmod(fabs(MET_PHI[0] - higgsPhi), TMath::TwoPi());
      if(dPhihiggs > TMath::Pi()) dPhihiggs = TMath::TwoPi() - dPhihiggs;
      
      //HT sum over higgs candidate
      if(blikeCentralJet_Pt.size()!=0){
        HT_CandidateJets_extraJet = btagCentralJet_Pt[0] + btagCentralJet_Pt[1] + blikeCentralJet_Pt[0];
      }
    }
    else if(order_bTagCentral.size()==1 && order_blikeCentral.size()>=1){
      // leading higg candidate jet 
      higgsCandidateJet_leadPT =  max(btagCentralJet_Pt[0],blikeCentralJet_Pt[0]);
      fillHistogram("higgs Candidate Jet PT",100,0,1500,higgsCandidateJet_leadPT);  //histogram 
      
      // higgs reco
      higgsCandidate = order_bTagCentral[0]+order_blikeCentral[0];
      fillHistogram("reco higgs inv mass",100,0,800,higgsCandidate.M());

      // deltaR
      delta_Rjhjh=fabs(order_bTagCentral[0].DeltaR(order_blikeCentral[0])); 
      
      // angle separation
      Double_t higgsPhi = higgsCandidate.Phi();
      dPhihiggs = fmod(fabs(MET_PHI[0] - higgsPhi), TMath::TwoPi());
      if(dPhihiggs > TMath::Pi()) dPhihiggs = TMath::TwoPi() - dPhihiggs;
      
      //HT sum over higgs b-tag candidate
      if(blikeCentralJet_Pt.size() > 1) {
	HT_CandidateJets_extraJet = btagCentralJet_Pt[0] + blikeCentralJet_Pt[0] + blikeCentralJet_Pt[1];
      }
    }
    
    //////////////////////
    //HT sum of central PT
    //////////////////////
    HTcentral2 = 0;
    if(centralPt.size()==2){
      HTcentral2 = centralPt[0]+centralPt[1];
    }
    HTcentral3 = 0;
    if (centralPt.size()>=3){
      HTcentral3 = centralPt[0]+centralPt[1]+centralPt[2];
    }
    
    
    /////// Boosted Region ////////
    ////////////////
    //Small-R0.2 Jet
    ////////////////
    akt2trackTLV.clear();
    akt2trackBTag.clear();
    for(int j=0;j<nJetAK2Track ;j++){
      if(JetAK2Track_PT[j]>10  && fabs(JetAK2Track_Eta[j])<2.5){
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(JetAK2Track_PT[j],JetAK2Track_Eta[j],JetAK2Track_Phi[j],JetAK2Track_Mass[j]);
        akt2trackTLV.push_back(tlv);
        akt2trackBTag.push_back(JetAK2Track_BTag[j]);
        akt2trackPhi.push_back(JetAK2Track_Phi[j]);
      }
    }

    //////////////       
    // Large-R Jet          
    //////////////    
    akt10trimTLV.clear();
    leadingLargeRjetPt.clear();
    //cout <<"SC " << nJetAK10Trim << endl;             
    largeRjet_counter = 0;
    for(int j=0;j<nJetAK10Trim ;j++){
      //if(akt10trimTLV.size()==0)break;
      if(JetAK10Trim_PT[j]>200 && fabs(JetAK10Trim_Eta[j])<2.0){
        largeRjet_counter++;
        leadingLargeRjetPt.push_back(JetAK10Trim_PT[j]);
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(JetAK10Trim_PT[j],JetAK10Trim_Eta[j],JetAK10Trim_Phi[j],JetAK10Trim_Mass[j]);
        akt10trimTLV.push_back(tlv);
        fillHistogram("akt10M",100,0,300,tlv.M());
      }
    }
    
    //Dissociate central b-tags 
    iso_dR = 0; 
    for(int j; j<order_bTagCentral.size(); j++){
      if(akt10trimTLV.size()!=0){
	for(int i; i<akt10trimTLV.size(); i++){
	  Double_t deltaR_SRbjLRj = fabs(order_bTagCentral[j].DeltaR(akt10trimTLV[i]));
	  if (deltaR_SRbjLRj>1) iso_dR++;
	}
      }
    }
    
    //HT-ratio selection
    Double_t SUM= 0;
    if(akt10trimTLV.size()!=0){
      for(int i; i<akt10trimTLV.size(); i++){
	for(int j; j<centralRjet.size(); j++){
	  Double_t deltaR_SRjLRj = fabs(centralRjet[j].DeltaR(akt10trimTLV[i]));
	  if(deltaR_SRjLRj>1) SUM += centralPt[j];
	}
	dissociate_CentralJet_Pt_sum.push_back(SUM);
      }
    }

    //////// Get histogram info /////////                                  
    fillHistogram(TString("num_electron_before"),3,0,3,nElec,1.);
    fillHistogram(TString("num_muon_before"),3,0,3,nMUON,1.);
    fillHistogram(TString("num_jet_before"),30,0,30,nJET,1.);
    fillHistogram(TString("MET_before"),4,150,800,MET[0],1.);
    int jetcount = 0;
    for(int i=0; i<nJET; i++){
      jetcount ++;
      fillHistogram(TString("jet_momentum"),100,0,3000,JET_PT[i],1.);
      fillHistogram(TString("jet_rapidity_before"),100,-6,6,JET_ETA[i],1.);
    }
    int bquarkcount = 0;
    for(int i=0; i<nParticle; i++){
      if(Particle_PID[i] == 5){
	bquarkcount++;
      }
    }
    fillHistogram(TString("num_bjet_before"),10,0,10,bquarkcount,1.);
    
    passResolved(&cut_flow[0]);
    passBoost(&cut_flow_boost[0]);
  
  } // End Of Events
                                           
  //////// output cut info //////////   
  cout << "store cut info ..." << endl;
  cout << "My MET count " << myMETcount << endl;
  //Resolved
  ofstream fileoutR("cutinfoR_local.dat");
  TString cutnameR[17]={"no selection", "No e or mu", "MET >150","pTmiss >30 (1b only)",
			"min[DeltaPhi(MET,jets)]>pi/9", "DeltaPhi(MET,pTmiss)>pi/2",
			"Num of central small Rjet>=2", "leading Higgs candidate small Rjet PT>45",
			"HT2jet>120 or HT3jet>150 (central)", "DeltaPhi(MET,higssPT)>2pi/3",
			"DeltaR(higgCanJet1,higgCanJet2)<1.8", "Veto on events with > 2btag", 
			"Sum of PT of two Higgs candidate jets and leading extra jet > 0.63xHT",
			"b-tagging: one or two small R jet",
			"MET>200","MET > 300","MET>450"};
 
  for(int j=0;j<17;j++) fileoutR << "[" << j << "] " << cutnameR[j] << ": " << cut_flow[j] <<endl; 
  fileoutR.close();
  //Boosted
  ofstream fileoutB("cutinfoB_local.dat");
  TString cutname[10]={"no selection", "No e or mu", "MET >500","pTmiss >30 (1b only)",
                       "min[DeltaPhi(MET,jets)]>pi/9","DeltaPhi(MET,pTmiss)>pi/2","Number of Large Rjet>=1",
		       "Veto on bjets not associate with Large Rjet", "HT ratio selection(<0.57)",
		       "b-tagging: one or two ID trackjets matched large Rjet"};
  for(int j=0;j<10;j++) fileoutB << "[" << j << "] " << cutname[j] << ": "<< cut_flow_boost[j] <<endl;
  fileoutB.close();
  
  /*
  //Analyze data storage
  Double_t fvar = fnum ;
  Double_t kvar = nknum ;

  // resolved
  ofstream fileoutAr("resolved_hh_acc_eff.dat", ios::app);
  fileoutAr << fvar << " " << kvar << " " << cut_flow[10] << " " << cut_flow[11] << " " << cut_flow[12] << " " << cut_flow[13] << "" << endl;
  fileoutAr.close();

  // Boost
  ofstream fileoutAb("boosted_hh_acc_eff.dat", ios::app);
  fileoutAb << fvar << " " << kvar << " " << cut_flow_boost[10] << " " << cut_flow_boost[11] << " " << cut_flow_boost[12] << " " << cut_flow_boost[13] << "" << endl;
  fileoutAb.close();
  */
  
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
    c[ipic]->SetLogy();
    c[ipic]->SaveAs();
  }
  fhistout->Close();

  //bkg hist
  fout->cd();
  tt_MET->Write();
  //Zj_MET->Write();
  TCanvas *cB;
  cB = new TCanvas();
  tt_MET->Draw();
  //Zj_MET->Draw();
  cB->SetLogy();
  cB->SaveAs();
 
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


bool passBoost (Long64_t* cut_flow_boost){

  if(iso_flag==1) return false;
  cut_flow_boost[1]++;


  if(ETmiss<=500) return false;
  cut_flow_boost[2]++;


  if(btagCentralCounter==1){
    if(PTmiss <=30) return false;
  }
  cut_flow_boost[3]++;


  if(dPhiSmallRMin <=TMath::Pi()/9.) return false;
  cut_flow_boost[4]++;


  if(dPhiPtMiss >=TMath::Pi()/2.) return false;
  cut_flow_boost[5]++;


  if(largeRjet_counter <1) return false;
  cut_flow_boost[6]++;

 
  if(iso_dR>0) return false;
  cut_flow_boost[7]++;

  for(int i=0; i<dissociate_CentralJet_Pt_sum.size(); i++){
    if(dissociate_CentralJet_Pt_sum[i] >0.57*(dissociate_CentralJet_Pt_sum[i]+leadingLargeRjetPt[0])) return false;
  }
  cut_flow_boost[8]++;

  
  //b-tagging: one or two ID track jets matched to large Rjet
  int Nbtag_track[2]={0,0};
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
    if(Bool_t(akt2trackBTag[j] &(1<<0))==false) continue;
    Nbtag_track[match]++;
  }
  if((Nbtag_track[0]+Nbtag_track[1])!=1 && (Nbtag_track[0]+Nbtag_track[1])!=2 ) return false;
  cut_flow_boost[9]++;
  
}


bool passResolved(Long64_t* cut_flow){
  
  if(iso_flag==1) return false;
  cut_flow[1]++;
  

  if(ETmiss <=150) return false;
  cut_flow[2]++;
 

  if(btagCentralCounter==1){
    if(PTmiss <=30) return false;
  }
  cut_flow[3]++;


  if(dPhiSmallRMin <=TMath::Pi()/9.) return false;
  cut_flow[4]++;


  if(dPhiPtMiss >=TMath::Pi()/2.) return false;
  cut_flow[5]++;
  

  if(NRjet_central <2) return false;
  cut_flow[6]++;


  if(higgsCandidateJet_leadPT <=45) return false;
  cut_flow[7]++;

  
  if(HTcentral2 <=120 && HTcentral3 <=150) return false;
  cut_flow[8]++;


  if(dPhihiggs <=TMath::TwoPi()/3.) return false;
  cut_flow[9]++;


  if(delta_Rjhjh >=1.8) return false;
  cut_flow[10]++;

  
  if(btagCentralCounter >2 ) return false;
  cut_flow[11]++;
   

  if(HT_CandidateJets_extraJet <=0.63*SCALAR_HT ) return false;
  cut_flow[12]++;


  if(num_btagSmallJet!=1 && num_btagSmallJet!=2) return false;
  cut_flow[13]++;

  //////////////////////
  //(11) MET > 200    //
  //////////////////////
  if(ETmiss<=200.) return false;
  cut_flow[14]++;
  
  //////////////////////
  //(12) MET > 350    //
  //////////////////////
  if(ETmiss<=350.) return false;
  cut_flow[15]++;
  
  //////////////////////
  //(13) MET > 500    //
  //////////////////////
  if(ETmiss<=500.) return false;
  cut_flow[16]++;
  
}


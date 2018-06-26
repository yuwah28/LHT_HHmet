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


//Functions Cells
bool passBoost(Long64_t* cutflow);
int getBosonTag(TLorentzVector reco, vector<TLorentzVector>& truth, float coneR=1.0);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;
TRandom3 randGen;

//Global Variables
//Boost
vector<TLorentzVector> akt10trimTLV;
vector<TLorentzVector> akt2trackTLV;
vector<TLorentzVector> HiggsCandidateTLV;
vector<TLorentzVector> BosonCandidateTLV;

vector<int> Nakt10trimBosonTag;
vector<int> akt10trimBosonTag;
vector<UInt_t> akt2trackBTag;
vector<double> selectLargeRjetPt;

Double_t rapidity_diff;
Double_t ETmiss;
Double_t dPhi;
Double_t mJJ;
Int_t Windex;
Int_t num_Lep;

Double_t PI = TMath::Pi();
Int_t eventnum;


//main
int whrso_analysis_ver4(){
 
  randGen.SetSeed(2);

  TChain *t = new TChain("Delphes");
  //t->Add("/phys/groups/tev/scratch4/users/paklim/LHT_Data_Backup/save_scan_jets_met_13TeV/f1200/jmetk05/cmswboost.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/WHmetReso/Validation/4processSample_f1200k5/Hdihiggs*.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/SMBackground/ttbarJets/BWHmet_Data/tj*.root");
  //t->Add("/phys/groups/tev/scratch1/users/paklim/AnalysisCodes/WHmetReso/Validation/4processSample_f1200k5/f1200k5_ATLASboost_DATA/Hdihiggs*.root");

  t->Add("/phys/groups/tev/scratch1/users/paklim/Samples/4processSample_f1200k5/f1200k5_Hybrid_WH_DATA/Hdihiggs*.root");

  Long64_t nevent = t->GetEntries();
  cout << nevent << endl;
  //nevent = 300000;
  t->SetMakeClass(1);


  //Root branches
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

  //Cutflow array 
  Long64_t cut_flow_boost[20];  std::fill_n(cut_flow_boost, 20, 0);
  int NDivision = nevent/10;
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
    //Begin analysis
    //
    cut_flow_boost[0]++;
    
    //Lepton Veto
    Int_t num_elec = 0;
    Int_t num_muon = 0;
    num_Lep = 0;
    for(int j=0; j<nMUON; j++){     
      if(MUON_PT[j]>7. && fabs(MUON_ETA[j])<2.5 && MUON_SUMPT[j]/MUON_PT[j]<0.15) num_muon++;
    }
    fillHistogram("Nmuon",5,0,5,num_muon);
    for(int j=0; j<nElec; j++){
      if(ELEC_PT[j]>7. && fabs(ELEC_ETA[j])<2.47 && ELEC_SUMPT[j]/ELEC_PT[j]<0.15) num_elec++;
    }
    fillHistogram("Nelectron",5,0,5,num_elec);
    num_Lep = num_muon + num_elec;
   
    //MET        
    if(nMET<1) continue;
    ETmiss = MET[0];
    
    //
    //Identify truth type event, jets+MET, h+MET, w+MET, hh+MET, wh+MET, ww+MET  
    //
    /*
    vector<int> genHiggsIndex, genWZIndex;
    for(int j=0;j<50;j++){
      if(Particle_Status[j]!=3) continue;
      int PID = abs(Particle_PID[j]);
      int M1  = Particle_M1[j];
      //Select W,Z,H                
      if(PID <23 || PID>25) continue;
      //Mother must exist and not Higgs (H->WW/H->ZZ are W/Z decayed from Higgs)
      if(M1<0 ) continue;
      int M1_PID = abs(Particle_PID[M1]);
      if(M1_PID==25) continue;
      if(PID==25) genHiggsIndex.push_back(j);
      else        genWZIndex.push_back(j);
    }
    if( (genHiggsIndex.size()+genWZIndex.size()) >2){
      cout <<"Event "<< ievent <<" NHiggs "<< genHiggsIndex.size() <<" N_WZ "<< genWZIndex.size() << endl;
      for(int j=0;j<50;j++){
	cout << j <<" "<< Particle_PID[j]<<" "<<Particle_Status[j] <<" M "<< Particle_M1[j] <<" D1 "<< Particle_D1[j]<<" "<<Particle_D2[j] <<" "<<Particle_PT[j] <<" "<<Particle_Eta[j] <<" "<<Particle_Phi[j] <<" "<< Particle_Mass[j] << endl;
      }
    }
    */

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

    //
    //Large R-Jets
    //
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
      if(CandFjetTLV[0].Pt()<=450) continue;
      if(CandFjetTLV[1].Pt()<=250) continue;
      akt10trimTLV.push_back(CandFjetTLV[k]);
      getFjetIndex.push_back(CandFIndex[k]);
      fillHistogram("Cand. LRjet",5,0,5,akt10trimTLV.size());
      fillHistogram("LRjet Mass distribution",100,0,800,akt10trimTLV[k].M());
    }
    /*
    if(hadronicWZTLV.size()!=0 && akt10trimTLV.size()!=0){
      cout << "DeltaR= " << akt10trimTLV[0].DeltaR(hadronicWZTLV[0]) << endl;
    }
    */

    //Get Boson Tag
    for(int k=0; k<akt10trimTLV.size(); k++){
      if(getBosonTag(akt10trimTLV[k], hadronicWZHTLV)==0) zeroWtag++;
      if(getBosonTag(akt10trimTLV[k], hadronicWZHTLV)==1) oneWtag++;
      akt10trimBosonTag.push_back( getBosonTag(akt10trimTLV[k], hadronicWZTLV) ); 
      fillHistogram("Boson Tag", 3,0,3, getBosonTag(akt10trimTLV[k], hadronicWZHTLV));
      /*
	for(int l=0; l<hadronicWZTLV.size(); l++){
	cout << ievent << "   *** Print DeltaR= " << akt10trimTLV[k].DeltaR(hadronicWZTLV[l]) << endl;
	}
      */
    }
  
    //
    //Wh candidacy    
    //
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
      fillHistogram("Higgs mass",100,0,1000, HiggsCandidateTLV[0].M());
      fillHistogram("W mass",100,0,1000, BosonCandidateTLV[0].M());
      //Rapidity difference  
      rapidity_diff = fabs(akt10trimTLV[0].Rapidity()-akt10trimTLV[1].Rapidity());
      fillHistogram("|delta_y12|",30,0,3,rapidity_diff);
      
      //cout << ievent << "  akt10Btag= "  << akt10trimBosonTag[Windex] << endl;

    }

   

    //Check W-tag
    int NRealWtag = 0;
    for(int k=0; k<hadronicWTLV.size(); k++){
      double wconeR = 1.0;
      int wmatch=-1;
      for(int j=0; j<BosonCandidateTLV.size(); j++){
        if(BosonCandidateTLV[j].DeltaR(hadronicWTLV[k])<wconeR){
          wconeR = BosonCandidateTLV[j].DeltaR(hadronicWTLV[k]);
          wmatch = 1;
        }
        if(wmatch==-1)continue;
        NRealWtag++;
      }
    }
    totalWtag += NRealWtag;
     
    //Check H b-tag
    int NRealHtag = 0;
    for(int k=0; k<hadronicHiggsTLV.size(); k++){
      double hconeR = 1.0;
      int hmatch=-1;
      for(int j=0; j<CandFjetTLV.size(); j++){
        if( CandFjetTLV[j].DeltaR(hadronicHiggsTLV[k])<hconeR){
          hconeR = CandFjetTLV[j].DeltaR(hadronicHiggsTLV[k]);
          hmatch = 1;
        }
        if(hmatch==-1)continue;
        NRealHtag++;
      }
    }
    totalHtag += NRealHtag;
       
    //
    //BTag TrackJet
    //
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
    
    //
    //Angle separation
    //
    if(HiggsCandidateTLV.size()!=0) {
      dPhi = fmod( fabs(MET_PHI[0]-HiggsCandidateTLV[0].Phi()),TMath::TwoPi() );
      if(dPhi>TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
      //cout << ievent << ":  angle separation= " << dPhi << endl;
    }
      
    //
    //Dijet Mass
    //
    mJJ = 0;
    if(HiggsCandidateTLV.size()!=0 && BosonCandidateTLV.size()!=0) {
      TLorentzVector DijetTLV = HiggsCandidateTLV[0]+BosonCandidateTLV[0];
      mJJ = DijetTLV.M();
    }
    //cout << ievent << ";  Dijet mass= " << mJJ << endl;
        

    eventnum=ievent;
    passBoost(&cut_flow_boost[0]);

  }//End Of Events
 
  cout << "totalW= " << totalW << "  totalWtag= " << totalWtag << endl;
  cout << "totalH= " << totalH << "  totalHtag= " << totalHtag << endl;
  cout << "getBosonTag -> 0 = " << zeroWtag << endl;
  cout << "getBosonTag -> 1 = "<< oneWtag << endl;
  
  //Output Analysis
  cout << "store cut info ..." << endl;

  //Boost
  ofstream fileoutB("cutinfoB_ttbar1M.dat");
  TString cutname[12]={"no selection", "Lepton veto",">1 large-R jets",
		       "Rapidity difference","MET veto",
		       "W Mass","W-tagging",
		       "Higgs Mass","Higgs b-tag=1",
		       "(SR1)Dijet Mass>1000",
		       "Higgs b=tag>1",
		       "(SR2)Dijet Mass>1000"};
  for(int j=0;j<12;j++) fileoutB << cutname[j] << "[" << j << "]: "<< cut_flow_boost[j] <<endl;
  fileoutB.close();

  
  //Save Histogram
  int ipic = 0;
  int nhist = 20;
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
    //cout << "cone size= " << coneR << endl;
    //cout << " Print Delta in getBosonTag func= " << reco.DeltaR(truth[i]) << endl;
  }
  float efficiency = QCD;
  if(isMatch==true) efficiency = Signal;
  
  int BitNumber = 0;
  int BosonTag = 0;
  // randGen.Uniform() generate distribution in [0,1]
  // Operand:= (A <= B) << Bit; if () is true we have 1, else 0
  // BosonTag = BosonTag(default=0) | Operand. If 0 | 0 we have 0 else we have 1.
  BosonTag |= (randGen.Uniform() <= efficiency)  << BitNumber;
  //double randnum = ((double) rand() / (RAND_MAX));
  //BosonTag |= (     randnum     <= efficiency)  << BitNumber;
  //cout << ((double) rand() / (RAND_MAX)) << endl;
  /*
  if( (randGen.Uniform() <= efficiency) && efficiency == Signal ){
    cout << "   Print Signal= " << Signal 
	 << "   Print BosonTag= " << BosonTag 
	 << "   and distribution= " << randGen.Uniform() 
	 << "   also print eff= "  << efficiency<< endl;
    
  }
  */
  return BosonTag;
}


bool passBoost(Long64_t* cut_flow){

  //Lepton
  if(num_Lep>0) return false;
  cut_flow[1]++;

  if(akt10trimTLV.size()<1) return false;
  cut_flow[2]++;
  
  if(rapidity_diff>=1.6) return false;
  cut_flow[3]++;
  
  if(dPhi>TMath::Pi()*120/180 && ETmiss>150) return false;
  cut_flow[4]++;

  //if(BosonCandidateTLV.size()==0) return false;
  //if(BosonCandidateTLV.size()==0 || (BosonCandidateTLV[0].M()<=50 || BosonCandidateTLV[0].M()>=70)) return false;  //VR-SR 
  if(BosonCandidateTLV[0].M()<75 || BosonCandidateTLV[0].M()>100) return false;     
  cut_flow[5]++;
  fillHistogram("W Candidate Mass",100,0,100,BosonCandidateTLV[0].M());

  //W tagging     
  if(akt10trimBosonTag[Windex] != 1) return false;
  cut_flow[6]++;

  //Higgs tagging
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
 
  if(HiggsCandidateTLV.size()==0 || (HiggsCandidateTLV[0].M()<=100 || HiggsCandidateTLV[0].M()>=145)) return false;
  cut_flow[7]++;
  fillHistogram("Higgs Candidate Mass",100,0,150,HiggsCandidateTLV[0].M());


  if(NBtag_H==1){
    //if(NBtag_H!=1) return false;
    cut_flow[8]++;

    if(mJJ<1000) return false;
    cut_flow[9]++;
  }

  if(NBtag_H>1){
    cut_flow[10]++;
    
    if(mJJ<1000) return false;
    cut_flow[11]++;
  }

  /*
  if(NBtag_H>0){
    if(NBtag_H!=1) return false;
    cut_flow[8]++;
    
    if(mJJ<1000) return false;
    cut_flow[9]++;
  }
  
  if(NBtag_H<2) return false;    
  cut_flow[10]++;
    
  if(mJJ<1000) return false;
  cut_flow[11]++; 
  */

  //if(NBtag_H==1)     cut_flow[9]++;
  //else if(NBtag_H>2) cut_flow[10]++;
  
  
  return true;
  
}


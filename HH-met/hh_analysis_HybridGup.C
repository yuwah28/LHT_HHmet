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
#include "TCanvas.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

//global variables:
Int_t eventnum;
//Boost
vector<TLorentzVector> akt10trimTLV;
vector<TLorentzVector> akt2trackTLV;
vector<UInt_t> akt2trackBTag;
Int_t resolveCutSize;
Int_t boostCutSize;

//Resolved
Int_t num_Lep;
Int_t num_JET;
Int_t Nb_T ;
Int_t Nb_M ;
Int_t Nb_L ;
double ETmiss ;
Int_t iso_flag;
Int_t dphi_flag;
Double_t dMmin;
Double_t dRmax;
Double_t averageM;

//Hybrid
//vector<Int_t> resolvedCut;
//vector<Int_t> boostCut;

//bool passBoost(Long64_t* cutflow);
bool passBoost(Long64_t* cutbinaryB);
bool passResolved(Long64_t* cutbinaryR);
bool passHybrid(Long64_t* cutflow, Long64_t* cutflow2, Long64_t* cutbinaryR, Long64_t* cutbinaryB, Long64_t* cutflowH);
void fillHistogram(TString hname, double nBin, double min, double max,
                   double value, double weight=1.0);
std::map<TString, TH1F*> h1f_map;

int hh_analysis_HybridGup(){
  

  TChain *t = new TChain("Delphes");

  //t->Add("/phys/groups/tev/scratch4/users/paklim/LHT_Data_Backup/save_scan_HH_met_13TeV/ffnum/hhmetkknum/cmsboost.root");

  //t->Add("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/ffnum/CHNAMEkknum/HHhybr.root");
  
  //t->Add("/phys/groups/tev/scratch1/users/paklim/HLSamples/3processes/WH_HH_met3/HHffnumkknumb.root");
 
  t->Add("/phys/groups/tev/scratch1/users/paklim/LHTDATA/save_scan_jets_met_13TeV/ffnum/CHNAMEkknum/HHhybr.root");

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
  
  const Int_t kMaxJetAK2Track = 10000; 
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
  

  /////////////////////////////////
  // Start Cut Analysis      
  Long64_t resolvedCutFlow[20]; std::fill_n(resolvedCutFlow, 20, 0);
  Long64_t boostedCutFlow[20]; std::fill_n(boostedCutFlow, 20, 0);
  Long64_t combinedCutFlow[20]; std::fill_n(combinedCutFlow, 20, 0);
  
  //nevent=10000;
  int NDivision = nevent/10;
 
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

  
    
    akt10trimTLV.clear();
     // Boosted Jet selection
    Double_t Lead_Fjet_mass, Sub_Fjet_mass;
    for(int j=0;j<nJetAK10Trim ;j++){
      if(JetAK10Trim_PT[j]>200 && fabs(JetAK10Trim_Eta[j])<2.5){
        TLorentzVector tlv;
        tlv.SetPtEtaPhiM(JetAK10Trim_PT[j],JetAK10Trim_Eta[j],JetAK10Trim_Phi[j],JetAK10Trim_Mass[j]);
        akt10trimTLV.push_back(tlv);
	Lead_Fjet_mass = akt10trimTLV[0].M();
        Sub_Fjet_mass  = akt10trimTLV[1].M();
        fillHistogram("lead_akt10 Mass", 100, 0, 400, Lead_Fjet_mass);
        fillHistogram("sub_akt10 Mass", 100, 0,400, Sub_Fjet_mass);
	fillHistogram("num akt10M",5,0,5,akt10trimTLV.size());
      }
    }

    //BTag TrackJet
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

    //Ol + 4-5 jets 
    num_Lep = 0;
    num_JET = 0;
    
    for(int j=0; j<nJET; j++){
      //if(j>=20) break;
      if(j>=MaxnJet) break;
      if(JET_PT[j] > 30. && fabs(JET_ETA[j]) < 2.4) num_JET++;
    }
    for(int j=0; j<nMUON; j++){
      //if(j>=10) break;
      if(j>=MaxnMUON) break;
      if(MUON_PT[j] > 10. && fabs(MUON_ETA[j]) < 2.4) num_Lep++;
    }
    for(int j=0; j<nElec; j++){
      //if(j>=10) break;
      if(j>=MaxnELEC) break;
      if(ELEC_PT[j] > 10. && fabs(ELEC_ETA[j]) < 2.5) num_Lep++;
    }
    
    
   
    //N(b,T) >= 2 
    vector<int> tightBindex,  mediumBindex, looseBindex, otherIndex;
    for(int j=0; j<nJET; j++){
      //if(j>=20) break;
      if(j>=MaxnJet) break;
      //With basic jet selection you want;
      //After that, you do exclusive jet index 
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
	if( Bool_t(JET_BTAG[j] & (1<<2)) ) tightBindex.push_back(j);
	else if( Bool_t(JET_BTAG[j] & (1<<1)) ) mediumBindex.push_back(j);
	else if( Bool_t(JET_BTAG[j] & (1<<0)) ) looseBindex.push_back(j);
	else otherIndex.push_back(j);
      }
    }
    Nb_T = tightBindex.size();
    Nb_M = Nb_T + mediumBindex.size();
    Nb_L = Nb_M + looseBindex.size();

    //Fill two Tight b-jet
    if(Nb_T>=2){
    Int_t tightB0 = tightBindex.at(0);
    Int_t tightB1 = tightBindex.at(1);
    }
   
    //MET > 150    
    if(nMET<1) continue;
    ETmiss = MET[0];
       
    //Track Veto    
    iso_flag = 0;
    for(int j=0; j<nElec; j++){
      //if(j>=10) break;
      if(j>=MaxnELEC) break;
      if( ELEC_PT[j]<5. || ELEC_PT[j]/ELEC_SUMPT[j]<5. ) continue;
      Double_t deltaphi = fmod( fabs(MET_PHI[0] - ELEC_PHI[j]), TMath::TwoPi() );
      if(deltaphi > TMath::Pi()) deltaphi = TMath::TwoPi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*ELEC_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1; 
    }
    for(int j=0; j<nMUON; j++){
      //if(j>=10) break;
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
      //it's better to use a denominator which guarantee to be non-zero.
      if(TRACK_SUMPT/TRACK_PT[j] > 0.1 ) continue;
      Double_t deltaphi = fmod( fabs(MET_PHI[0] - TRACK_PHI[j]), TMath::TwoPi() );
      if(deltaphi > TMath::Pi()) deltaphi = TMath::TwoPi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*TRACK_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }

    
    //Delta Phi   
    Double_t deltaphi[4];
    Int_t jetIndex[4]; std::fill_n(jetIndex, 4, -1);
    Int_t integer = 0;
    for(int j=0; j<nJET; j++){
      //if(j>=20) break;
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
    }

    dphi_flag=0;
    if( (deltaphi[0]>0 && deltaphi[0]<=0.5) || 
        (deltaphi[1]>0 && deltaphi[1]<=0.5) || 
        (deltaphi[2]>0 && deltaphi[2]<=0.3) ||
        (deltaphi[3]>0 && deltaphi[3]<=0.3))         dphi_flag=1;

    //if(ievent == 470) cout << ievent << "  flag dphi_flag= " << dphi_flag << endl;
    // Reco Higgs pair           
    // b-tag indentification                                              
    // BitNumber0:Loose 1:Medium 2:Tight
    TLorentzVector b_vector[4];
    vector<int> Bindex;

    for(int j=0; j< tightBindex.size(); j++) Bindex.push_back( tightBindex.at(j));
    for(int j=0; j< mediumBindex.size(); j++) Bindex.push_back( mediumBindex.at(j));
    for(int j=0; j< looseBindex.size(); j++) Bindex.push_back( looseBindex.at(j));
    for(int j=0; j< otherIndex.size(); j++) Bindex.push_back(otherIndex.at(j));
    for(int k=0;k< Bindex.size();k++){
      if(k>=4) break;
      int j = Bindex.at(k);
      b_vector[k].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
    }

    // Higgs Reco
    TLorentzVector aux1, aux2;
    Double_t delta_m[3], mean_m[3], delta_Ra[3], delta_Rb[3];
    Int_t z1=-1, z2=-1, z3=-1;
    for(int j=0; j<3; j++){
      z1 = j + 1;
      z2 = (z1+3) % 3 + 1;
      z3 = (z1+1) % 3 + 1;
      
      aux1=b_vector[0]+b_vector[z1];
      aux2=b_vector[z2]+b_vector[z3];
      delta_m[j]=fabs(aux1.M()-aux2.M());
      mean_m[j]=(aux1.M()+aux2.M())/2.0;
      delta_Ra[j]=fabs(b_vector[0].DeltaR(b_vector[z1]));
      delta_Rb[j]=fabs(b_vector[z2].DeltaR(b_vector[z3]));
    }
   
    //Delta mass 
    Double_t mini = delta_m[0];
    Int_t index_min = 0;
    if(delta_m[1]<mini){ 
      mini = delta_m[1]; index_min = 1; 
    }
    if(delta_m[2]<mini){
      mini = delta_m[2]; index_min = 2;
    }
    dMmin  = delta_m[index_min];

    //Delta Rmax  
    dRmax = TMath::Max(delta_Ra[index_min], delta_Rb[index_min]);
   
    //Higg inv mass
    averageM = mean_m[index_min];

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

    /*
    for(int j=0;j<14;j++){
      cout << "check resolvedCut array update= " << resolvedCut[j] << endl;
    }
    */
 
  } // End Of Events
  /////////////////////                                                       
  // output cut info //                                                       
  ///////////////////// 
  cout << "store cut info ..." << endl;

 
  ofstream fileout("cutinfoR_3process.dat");
  fileout << "No Selection : "    << " " << resolvedCutFlow[0] << endl;
  fileout << "Ol + (4~5) Jets : " << " " << resolvedCutFlow[1] << endl;
  fileout << "N_{b,T} >= 2 : "    << " " << resolvedCutFlow[2] << endl;
  fileout << "MET > 150 : "       << " " << resolvedCutFlow[3] << endl;
  fileout << "Track Veto : "      << " " << resolvedCutFlow[4] << endl;
  fileout << "#Delta#phi_{1,2}>0.5,#Delta#phi_{3,4}>0.3 : " << " " << resolvedCutFlow[5] << endl;
  fileout << "|#Deltam|<40GeV : "    << " " << resolvedCutFlow[6] << endl;
  fileout << "#DeltaR_{max}<2.2 : "  << " " << resolvedCutFlow[7] << endl;
  fileout << "100<#bar{m}<140GeV : " << " " << resolvedCutFlow[8] << endl;
  fileout << "3b+4b : " << " " << resolvedCutFlow[9] << endl;
  fileout << "4b : "    << " " << resolvedCutFlow[10] << endl;
  fileout << "MET > 200 : " << " " << resolvedCutFlow[11] << endl;
  fileout << "MET > 300 : " << " " << resolvedCutFlow[12] << endl;
  fileout << "MET > 450 : " << " " << resolvedCutFlow[13] << endl;
  fileout.close();

  //Boost
  ofstream fileoutB("cutinfoB_3process.dat");
  TString cutname[14]={"no selection", 
		       "0l+2akt10trim",
		       "N_btag>=2", 
		       "MET >150",
		       "Track veto",
                       "DeltaPhi",
		       "dummy",
		       "M1",
		       "M2",
		       "N1_btag>=2",
                       "N2_btag>=2",
		       "MET >200",
		       "MET > 300",
		       "MET>450"};
  for(int j=0;j<14;j++) fileoutB << cutname[j] << "[" << j << "]: "<< boostedCutFlow[j] <<endl;
  fileoutB.close();

  //Combine
  ofstream fileoutC("cutinfoC_3process.dat");
  TString cutnameH[14]={"no selection", "c1","c2", "c3","c4",
                       "c5","c6","c7","c8",
			"c9","c10","c11","c12","c13"};
  for(int j=0;j<14;j++) fileoutC << cutnameH[j] << "[" << j << "]: "<< combinedCutFlow[j] <<endl;
  fileoutC.close();

//Analyze data storage
  Double_t fvar = 1000 ;
  Double_t kvar = 1. ;

  
  // Resolved
  ofstream fileoutAr("resolved_hh_acc_HL.dat", ios::app);
  fileoutAr << fvar << " " << kvar << " " << resolvedCutFlow[10] << " " << resolvedCutFlow[11] << " " << resolvedCutFlow[12] << " " << resolvedCutFlow[13] << "" << endl;
  fileoutAr.close();

  // Boost
  ofstream fileoutAb("boosted_hh_acc_HL.dat", ios::app);
  fileoutAb << fvar << " " << kvar << " " << boostedCutFlow[10] << " " << boostedCutFlow[11] << " " << boostedCutFlow[12] << " " << boostedCutFlow[13] << "" << endl;
  fileoutAb.close();
  
  // Hybrid
  ofstream fileoutAh("hybrid_hh_acc_HL.dat", ios::app);
  fileoutAh << fvar << " " << kvar << " " << combinedCutFlow[10] << " " << combinedCutFlow[11] << " " << combinedCutFlow[12] << " " << combinedCutFlow[13] << "" << endl;
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

bool passBoost(Long64_t* cutbinaryB){
 
  if(num_Lep>0 || akt10trimTLV.size()<2) return false;
  //if(akt10trimTLV.size()<2) return false;
  cutbinaryB[1]=1;
  
  //Btag  each trackjet is assigned to the closet akknum0Trimjet
  int NBtag_T[2]={0,0}, NBtag_M[2]={0,0}, NBtag_L[2]={0,0};
  int Nbtag_track[2]={0,0};
  for(int j=0;j<akt2trackTLV.size();j++){
    double minR=1.0; 
    int match=-1;
    for(int i=0;i<2;i++){
      if(akt10trimTLV[i].DeltaR(akt2trackTLV[j]) <minR){
	minR=akt10trimTLV[i].DeltaR(akt2trackTLV[j]);
	match=i;
      }
    }
    if(match==-1) continue;
    if(akt2trackBTag[j]>0) Nbtag_track[match]++;
    /*
    if( Bool_t(akt2trackBTag[j] & (1<<0)) == false) continue;
    NBtag_L[match]++;
    if( Bool_t(akt2trackBTag[j] & (1<<1)) == false) continue;
    NBtag_M[match]++;
    if( Bool_t(akt2trackBTag[j] & (1<<2)) == false) continue;
    NBtag_T[match]++;
    */
  }

  // if( NBtag_T[0]<1 || NBtag_T[1]<1) return false;
  if(Nbtag_track[0]<1 || Nbtag_track[1]<1 ) return false;
  cutbinaryB[2]=1;
  
  if(ETmiss<=150) return false;
  cutbinaryB[3]=1;
  
  if(iso_flag==1) return false;
  cutbinaryB[4]=1;
  
  //cout << "eventnum=" << eventnum << ";   dphi cut #B= " << dphi_flag << endl;
  if(dphi_flag==1) return false;
  cutbinaryB[5]=1;
    
  //if(akt10trimTLV.size()<2) return false;
  cutbinaryB[6]=1;
  
  if(akt10trimTLV[0].M()<100 || akt10trimTLV[0].M()>140 ) return false;
  cutbinaryB[7]=1;
  
  if(akt10trimTLV[1].M()<100 || akt10trimTLV[1].M()>140 ) return false;
  cutbinaryB[8]=1;
  
  //if( NBtag_L[0] <2) return false;
  if(Nbtag_track[0]!=2) return false;
  cutbinaryB[9]=1;

  //if( NBtag_L[1] <2) return false;
  if(Nbtag_track[1]!=2) return false;
  cutbinaryB[10]=1;

  if(ETmiss<=200.) return false;
  cutbinaryB[11]=1;
  
  if(ETmiss<=300.) return false;
  cutbinaryB[12]=1;
  
  if(ETmiss<=450.) return false;
  cutbinaryB[13]=1;
  
  return true;
}

bool passResolved(Long64_t* cutbinaryR){
 
  if(num_JET<4 || num_JET>5 || num_Lep>0) return false;
  cutbinaryR[1]=1;
  
  //else cutbinaryR[1]=0;
  
  
  if(Nb_T<2) return false;
  cutbinaryR[2]=1;
 
  //cut_flow[2]++;
  
  if(ETmiss<=150) return false;
  cutbinaryR[3]=1;
  
  if(iso_flag==1) return false;
  cutbinaryR[4]=1;
  
  //cut_flow[4]++;
  //cout << "eventnum=" << eventnum << ";   dphi cut #R= " << dphi_flag << endl;
  if(dphi_flag==1) return false;
  cutbinaryR[5]=1;
  
  //cut_flow[5]++;
  
  if(dMmin  > 40.) return false;
  cutbinaryR[6]=1;
  
  // cut_flow[6]++;
  
  //(7) Delta Rmax  
  if(dRmax>=2.2) return false;
  cutbinaryR[7]=1;
  
  //cut_flow[7]++;
  
  //(8) Higg inv mass   
  if(averageM<=100. || averageM>=140.) return false;
  cutbinaryR[8]=1;
  
  //cut_flow[8]++;
  
  //(9) 3b+4b  
  if(Nb_M<3 || Nb_L<3) return false;
  cutbinaryR[9]=1;
  
  //cut_flow[9]++;
  
  //(10)  4b   
  if(Nb_M<3 || Nb_L<4) return false;
  cutbinaryR[10]=1;
  
  //cut_flow[10]++;
  
  //(11) MET > 200  
  if(ETmiss<=200.) return false;
  cutbinaryR[11]=1;
  
  //cut_flow[11]++;
  
  //(12) MET > 300  
  if(ETmiss<=300.) return false;
  cutbinaryR[12]=1;
  
  //cut_flow[12]++;
  
  //(13) MET > 450 
  if(ETmiss<=450.) return false;
  cutbinaryR[13]=1;
  
  //cut_flow[13]++;

  return true;
}


//flagE
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

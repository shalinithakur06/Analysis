#define FinalPlots_cxx
#include "FinalPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <fstream>
#include "muonSF.h"

using namespace std;

void FinalPlots::Loop()
{
   //Open the PU ratio file for pileup reweighting
   //TFile *puratiofile = TFile::Open("ratio_data_mc_pu_2017_10_22.root");
   //TH1F *puratiohisto = (TH1F*)puratiofile->Get("ratiohisto");
    
   if (fChain == 0) return;
   TFile *f = new TFile("TTbar_2017_10_16.root","RECREATE");

   TH1D *h_deltaR = new TH1D("h_deltaR","Jets Eta", 140, -7, 7);
   TH1D *h_nMuons = new TH1D("h_nMuons","Number of selected Muons", 20, 0, 20);
   TH1D *h_LeadMuonPt = new TH1D("h_LeadMuonPt","Muons Pt", 3000, 0, 3000);
   TH1D *h_LeadMuonEta = new TH1D("h_LeadMuonEta","Muons Eta", 100, -5, 5);
   TH1D *h_2ndMuonPt = new TH1D("h_2ndMuonPt","Muons Pt", 3000, 0, 3000);
   TH1D *h_2ndMuonEta = new TH1D("h_2ndMuonEta","Muons Eta", 100, -5, 5);

   TH1D *h_ZPt = new TH1D("h_ZPt","DiMuon Pt", 3000, 0, 3000);
   TH1D *h_ZMass = new TH1D("h_ZMass","DiMuon Pt", 3000, 0, 3000);
   TH1D *h_ZEta = new TH1D("h_ZEta","DiMuon Eta", 100, -5, 5);
   TH1D *h_ZPhi = new TH1D("h_ZPhi","Zs Phi", 68, -3.4, 3.4);

   TH1D *h_nJets = new TH1D("h_nJets","Number of selected Jets", 20, 0, 20);
   TH1D *h_JetPt = new TH1D("h_JetPt","Jets Pt", 3000, 0, 3000);
   TH1D *h_JetEta = new TH1D("h_JetEta","Jets Eta", 100, -5, 5);
   TH1D *h_JetPhi = new TH1D("h_JetPhi","Jets Phi", 68, -3.4, 3.4);
   TH1D *h_JetMass = new TH1D("h_JetMass","Jets Mass", 3000, 0, 3000);
   TH1D *h_JetSubJet = new TH1D("h_JetSubJet","Jets Mass", 500, 0, 1);
   ///preselection cuts
   TH1D *h_preSelnMuons = new TH1D("h_preSelnMuons","Number of selected Muons", 20, 0, 20);
   TH1D *h_preSelLeadMuonPt = new TH1D("h_preSelLeadMuonPt","Muons Pt", 500, 0, 500);
   TH1D *h_preSelLeadMuonEta = new TH1D("h_preSelLeadMuonEta","Muons Eta", 100, -5, 5);
   TH1D *h_preSel2ndMuonPt = new TH1D("h_preSel2ndMuonPt","Muons Pt", 500, 0, 500);
   TH1D *h_preSel2ndMuonEta = new TH1D("h_preSel2ndMuonEta","Muons Eta", 100, -5, 5);

   TH1D *h_preSelZPt = new TH1D("h_preSelZPt","DiMuon Pt", 500, 0, 500);
   TH1D *h_preSelZMass = new TH1D("h_preSelZMass","DiMuon Pt", 500, 0, 500);
   TH1D *h_preSelZEta = new TH1D("h_preSelZEta","DiMuon Eta", 100, -5, 5);
   TH1D *h_preSelZPhi = new TH1D("h_preSelZPhi","Zs Phi", 68, -3.4, 3.4);

   TH1D *h_preSelnJets = new TH1D("h_preSelnJets","Number of selected Jets", 20, 0, 20);
   TH1D *h_preSelJetPt = new TH1D("h_preSelJetPt","Jets Pt", 500, 0, 500);
   TH1D *h_preSelJetEta = new TH1D("h_preSelJetEta","Jets Eta", 100, -5, 5);
   TH1D *h_preSelJetPhi = new TH1D("h_preSelJetPhi","Jets Phi", 68, -3.4, 3.4);
   TH1D *h_preSelJetMass = new TH1D("h_preSelJetMass","Jets Mass", 500, 0, 500);
   TH1D *h_preSelJetSubJet = new TH1D("h_preSelJetSubJet","Jets Mass", 500, 0, 1);

   //tagging Z boson
   TH1D *h_ZTagLeadMuonPt = new TH1D("h_ZTagLeadMuonPt","Muons Pt", 500, 0, 500);
   TH1D *h_ZTagLeadMuonEta = new TH1D("h_ZTagLeadMuonEta","Muons Eta", 100, -5, 5);
   TH1D *h_ZTag2ndMuonPt = new TH1D("h_ZTag2ndMuonPt","Muons Pt", 500, 0, 500);
   TH1D *h_ZTag2ndMuonEta = new TH1D("h_ZTag2ndMuonEta","Muons Eta", 100, -5, 5);

   TH1D *h_ZTagZPt = new TH1D("h_ZTagZPt","DiMuon Pt", 500, 0, 500);
   TH1D *h_ZTagZMass = new TH1D("h_ZTagZMass","DiMuon Pt", 500, 0, 500);
   TH1D *h_ZTagZEta = new TH1D("h_ZTagZEta","DiMuon Eta", 100, -5, 5);
   TH1D *h_ZTagZPhi = new TH1D("h_ZTagZPhi","Zs Phi", 68, -3.4, 3.4);

   TH1D *h_ZTagnJets = new TH1D("h_ZTagnJets","Number of selected Jets", 20, 0, 20);
   TH1D *h_ZTagJetPt = new TH1D("h_ZTagJetPt","Jets Pt", 500, 0, 500);
   TH1D *h_ZTagJetEta = new TH1D("h_ZTagJetEta","Jets Eta", 100, -5, 5);
   TH1D *h_ZTagJetPhi = new TH1D("h_ZTagJetPhi","Jets Phi", 68, -3.4, 3.4);
   TH1D *h_ZTagJetMass = new TH1D("h_ZTagJetMass","Jets Mass", 500, 0, 500);
   TH1D *h_ZTagJetSubJet = new TH1D("h_ZTagJetSubJet","Jets Mass", 500, 0, 1);
  
   TH1D *h_nvtx = new TH1D("h_nvtx", "Number of Primary Vertices", 75, 0, 75);
   h_nvtx->Sumw2();
   TH1D *h_nGVtx = new TH1D("h_nGVtx", "Number of Good PV", 75, 0, 75);
   h_nGVtx->Sumw2();
   TH1D *h_nvtx_new = new TH1D("h_nvtx_new", "Number of Primary Vertices", 75, 0, 75);
   h_nvtx_new->Sumw2();
   TH1D *h_nGVtx_new = new TH1D("h_nGVtx_new", "Number of Good PV", 75, 0, 75);
   h_nGVtx_new->Sumw2();

   Long64_t nentries = fChain->GetEntries();
   //cout << "number of entries " << nentries << endl;
   //nentries = 100;
   bool debug(false); 
   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry % 100000 ==0) cout << "Processing Event Number " << jentry << endl;
     //if(debug) cout << "number of muons " << nMuons << endl;
     
     //double weight = 0.;
     //cout << puratiohisto->GetBinContent(1) << endl;
     //double puReweightingFactor = puratiohisto->GetBinContent(nGVtxs);
     //cout << "puReweightingFactor" << puReweightingFactor <<endl;

     //weight+=mc_EventWeight;
     //double wnorm=1.;

     double wnorm=PUwei;
     //cout << "wnorm =" << wnorm <<endl;
     //double efficiency=1.;
     //wnorm = puReweightingFactor;
     //double luminosity = 1.;
     //double cross_section = 1.; 
     //int totalNumEvsmc = 1.;
     //luminosity = 35867, cross_section = 6025.2, totalNumEvsmc = 73079437;   
     //wnorm = luminosity*cross_section/(totalNumEvsmc);
     h_nvtx_new->Fill(nVtxs,wnorm );
     h_nGVtx_new->Fill(nGVtxs, wnorm);

     if(nMuons < 2) continue;
     vector<TLorentzVector> nmu;
     vector<int> SelMuons;
     for(int i =0; i < nMuons; ++i){
       if(MuPt[i] < 35) continue;
       if( MuType[i] ==4) continue;
	else if( MuType[i] ==12) continue;  
       else if( MuType[i] ==36) continue;
       else if( MuType[i] ==40) continue;
       else if( MuType[i] ==44) continue;
       else if( MuType[i] ==68) continue;
       else if( MuType[i] ==76) continue;
       else if( MuType[i] ==96) continue;
       else if( MuType[i] ==100)continue; 
       else if( MuType[i] ==108)continue;
       TLorentzVector v;
       v.SetPtEtaPhiE(MuPt[i], MuEta[i], MuPhi[i], MuEn[i]);
       nmu.push_back(v);
     if(debug){
     cout << "Muon Pt/Eta/Phi/Charge " << MuPt[i]<< "\t" << MuEta[i]<< "\t" << MuPhi[i]<< "\t" << MuEn[i] << "\t" << MuCharge[i] << endl;
     }
     SelMuons.push_back(i);
     }
     if(nmu.size() < 2) continue;
     if(nmu.size() > 2) continue;
     if(nmu.size() ==2 && (MuCharge[0] == MuCharge[1])){
     continue;
     }
     TLorentzVector vz, vm1, vm2;
     int bestpair1(-1), bestpair2(-1); double dzmass(1000);
     if(nmu.size() > 2){
     //to calculate the best pair of Z
     for(size_t i=0; i < nmu.size(); ++i){
     for(size_t j=i+1; j < nmu.size(); ++j){
     vz=nmu[i] + nmu[j];
     if(debug){
     cout << "Combination " << i << ":" << j << "\t" << vz.M() << endl;
     }//matches debug condition.
     if( (vz.M() -90) < dzmass){
      bestpair1=i;
      bestpair2=j;
      dzmass=(vz.M() -90);
           }//if mass is close to 90 then select this pair.
         }//loop over muons-1
       }//loop over muons
     } //if more than two muons
     //if(dzmass <= 50) continue;

     //best pair of Z
     if(nmu.size() > 2){
     vm1=nmu[bestpair1];
     vm2=nmu[bestpair2];
     vz=vm1+vm2;
     if(debug){
     cout << "selected indices are : " << bestpair1 << ":" << bestpair2 <<"\t with mass : " << vz.M()  << endl;
     cout << "Lead Muon Pt/Eta/Phi/Energy " << vm1.Pt()<< "\t" << vm1.Eta()<< "\t" << vm1.Phi()<< "\t" << vm1.E()  << endl;
     cout << "2nd Muon  Pt/Eta/Phi/Energy " << vm2.Pt()<< "\t" << vm2.Eta()<< "\t" << vm2.Phi()<< "\t" << vm2.E()  << endl;
     cout << "Z  Pt/Eta/Phi/Mass " << vz.Pt()<< "\t" << vz.Rapidity()<< "\t" << vz.Phi()<< "\t" << vz.M()  << endl;
       }//matches debug;
     }//matches with more than two muon condition.
     else{
     vm1=nmu[0];
     vm2=nmu[1];
     vz=vm1+vm2;
     if(debug){
     cout << "Lead Muon Pt/Eta/Phi/Energy " << vm1.Pt()<< "\t" << vm1.Eta()<< "\t" << vm1.Phi()<< "\t" << vm1.E()  << endl;
     cout << "2nd  Muon Pt/Eta/Phi/Energy " << vm2.Pt()<< "\t" << vm2.Eta()<< "\t" << vm2.Phi()<< "\t" << vm2.E()  << endl;
     cout << "Z  Pt/Eta/Phi/Mass " << vz.Pt()<< "\t" << vz.Rapidity()<< "\t" << vz.Phi()<< "\t" << vz.M()  << endl;
       }//matches debug
     }//matches else condition
     //construct the Z candidate.
     ///fill the control plots
     //cout<<"nGVvtx wnorm "<<nGVtxs<<" "<<wnorm<<endl;
      if(vz.M()<60) continue;    
     //---------------------------------------------------//
     //get muon SF for 1st muon 
     //---------------------------------------------------//
     double lumi_BCDEF = 18.67; double lumi_GH = 18.17;	
     double lumi = lumi_BCDEF + lumi_GH;
     //trigger 	
     double muSFtrig_BCDEF 	= getMuonTrigSF(h2_trigSF_BCDEF, vm1.Eta(), vm1.Pt());
     double muSFtrig_GH 	= getMuonTrigSF(h2_trigSF_GH, vm1.Eta(), vm1.Pt());
     double muSFtrig 		= (muSFtrig_BCDEF*lumi_BCDEF + muSFtrig_GH*lumi_GH)/lumi; 
     //identification
     double muSFid_BCDEF 	= getMuonSF(h2_idSF_BCDEF, vm1.Eta(), vm1.Pt());
     double muSFid_GH 		= getMuonSF(h2_idSF_GH, vm1.Eta(), vm1.Pt());
     double muSFid 		= (muSFid_BCDEF*lumi_BCDEF + muSFid_GH*lumi_GH)/lumi; 
     //isolation 
     double muSFiso_BCDEF 	= getMuonSF(h2_isoSF_BCDEF, vm1.Eta(), vm1.Pt());
     double muSFiso_GH 		= getMuonSF(h2_isoSF_GH, vm1.Eta(), vm1.Pt());
     double muSFiso 		= (muSFiso_BCDEF*lumi_BCDEF + muSFiso_GH*lumi_GH)/lumi; 
     //tracking 
     double muSFtrack_BCDEF 	= getMuonTrackSF(tg_trackSF_BCDEF, vm1.Eta()); 
     double muSFtrack_GH 	= getMuonTrackSF(tg_trackSF_GH, vm1.Eta()); 
     double muSFtrack 		= (muSFtrack_BCDEF*lumi_BCDEF + muSFtrack_GH*lumi_GH)/lumi;
     
/*     cout<<"---------------------------"<<endl;
     cout<<"Muon PT = "<<vm1.Pt()<<endl;
     cout<<"Muon Eta = "<<vm1.Eta()<<endl;
     cout<<"\t"<<"muSFtrig_BCDEF = " 	<<muSFtrig_BCDEF<<endl;
     cout<<"\t"<<"muSFtrig_GH = " 		<<muSFtrig_GH<<endl;
     cout<<"\t"<<"muSFtrig = " 		<<muSFtrig<<endl;
     cout<<"\t"<<"muSFid_BCDEF = " 		<<muSFid_BCDEF<<endl;
     cout<<"\t"<<"muSFid_GH = " 		<<muSFid_GH<<endl;
     cout<<"\t"<<"muSFid = " 		<<muSFid<<endl;
     cout<<"\t"<<"muSFiso_BCDEF = " 	<<muSFiso_BCDEF<<endl;
     cout<<"\t"<<"muSFiso_GH = " 		<< muSFiso_GH <<endl;
     cout<<"\t"<<"muSFiso = " 		<< muSFiso <<endl;
     cout<<"\t"<<"muSFtrack_BCDEF = " 	<< muSFtrack_BCDEF <<endl;
     cout<<"\t"<<"muSFtrack_GH = " 		<<muSFtrack_GH  <<endl;
     cout<<"\t"<<"muSFtrack = " 		<< muSFtrack <<endl;
*/
     //combined SF
     double muSF =1.0;
     muSF = muSFtrig*muSFid*muSFiso*muSFtrack;	
     
     //---------------------------------------------------//
     //get muon SF for 2nd muon 
     //---------------------------------------------------//
     //trigger 	
     double muSFtrig_BCDEF2 	= getMuonTrigSF(h2_trigSF_BCDEF, vm2.Eta(), vm2.Pt());
     double muSFtrig_GH2 	= getMuonTrigSF(h2_trigSF_GH, vm2.Eta(), vm2.Pt());
     double muSFtrig2 		= (muSFtrig_BCDEF2*lumi_BCDEF + muSFtrig_GH2*lumi_GH)/lumi; 
     //identification
     double muSFid_BCDEF2 	= getMuonSF(h2_idSF_BCDEF, vm2.Eta(), vm2.Pt());
     double muSFid_GH2 		= getMuonSF(h2_idSF_GH, vm2.Eta(), vm2.Pt());
     double muSFid2 		= (muSFid_BCDEF2*lumi_BCDEF + muSFid_GH2*lumi_GH)/lumi; 
     //isolation 
     double muSFiso_BCDEF2 	= getMuonSF(h2_isoSF_BCDEF, vm2.Eta(), vm2.Pt());
     double muSFiso_GH2 		= getMuonSF(h2_isoSF_GH, vm2.Eta(), vm2.Pt());
     double muSFiso2 		= (muSFiso_BCDEF2*lumi_BCDEF + muSFiso_GH2*lumi_GH)/lumi; 
     //tracking 
     double muSFtrack_BCDEF2 	= getMuonTrackSF(tg_trackSF_BCDEF, vm2.Eta()); 
     double muSFtrack_GH2 	= getMuonTrackSF(tg_trackSF_GH, vm2.Eta()); 
     double muSFtrack2 		= (muSFtrack_BCDEF2*lumi_BCDEF + muSFtrack_GH2*lumi_GH)/lumi;
     
     //combined SF
     double muSF2 =1.0;
     //muSF2 = muSFtrig2*muSFid2*muSFiso2*muSFtrack2;	
       muSF2 = 0.98*muSFid2*muSFiso2*muSFtrack2;
     //apply SF

     wnorm =wnorm*muSF*muSF2; 
     /////////////////////////////////////////////////////////////////
     h_nvtx->Fill(nVtxs, wnorm);
     h_nGVtx->Fill(nGVtxs, wnorm);
     
     h_nMuons->Fill(nmu.size());
     h_LeadMuonPt-> Fill(vm1.Pt(), wnorm);
     h_LeadMuonEta->Fill(vm1.Eta(), wnorm);
     h_2ndMuonPt->  Fill(vm2.Pt(), wnorm);
     h_2ndMuonEta-> Fill(vm2.Eta(), wnorm);
     h_ZPt->Fill(vz.Pt(), wnorm);
     h_ZEta->Fill(vz.Rapidity(), wnorm);
     h_ZPhi->Fill(vz.Phi(), wnorm);
     h_ZMass->Fill(vz.M(), wnorm);
     h_nJets->Fill(nAK8Jets, wnorm);
     if(nAK8Jets){
     for(int ij=0; ij !=nAK8Jets; ++ij){
     if(debug){
     //cout << "Jet Pt/Eta/Phi/Mass/tau3 : " << aK8JetPt[ij] << "\t" << aK8JetEta[ij] << "\t" <<aK8JetPhi[ij] << "\t" <<aK8JetPrunedMass[ij]      << "\t" << aK8Jet_tau3[ij] << endl;
     }
     h_JetPt ->Fill(aK8JetPt[ij], wnorm);
     h_JetEta->Fill(aK8JetEta[ij], wnorm);
     h_JetPhi->Fill(aK8JetPhi[ij], wnorm);
     h_JetMass->Fill(aK8JetPrunedMass[ij], wnorm);
     h_JetSubJet->Fill(aK8Jet_tau3[ij], wnorm);
       }//loop over ak8 jets
     } //if there are ak8 jets
     else continue;
     if(debug) cout << "##############################################" << endl;
     //cut on 
     if(vz.Pt() <=100) continue;
     //applying cuts on the AK8Jets.
     vector<int> selectedJets;
     if(nAK8Jets){
     for(int ij=0; ij !=nAK8Jets; ++ij){
     if( aK8JetPt[ij] <= 100) continue;
     if( fabs(aK8JetEta[ij]) >= 2.4) continue;
     if(aK8JetPrunedMass[ij] <=40) continue; 
     //apply the deltaR cut here.
     TLorentzVector vj;
     vj.SetPtEtaPhiE(aK8JetPt[ij], aK8JetEta[ij], aK8JetPhi[ij], aK8JetEn[ij]);
     if(debug){
     cout << "Jet Pt/Eta/Phi/Mass/tau3 : " << aK8JetPt[ij] << "\t" << aK8JetEta[ij] << "\t" <<aK8JetPhi[ij] << "\t" <<aK8JetPrunedMass[ij]      << "\t" << aK8Jet_tau3[ij] << endl;
     cout << "Lead Muon Pt/Eta/Phi/Charge " << vm1.Pt()<< "\t" << vm1.Eta()<< "\t" << vm1.Phi()<< "\t" << vm1.E()  << endl;
     cout << "2nd Muon Pt/Eta/Phi/Charge " << vm2.Pt()<< "\t" << vm2.Eta()<< "\t" << vm2.Phi()<< "\t" << vm2.E()  << endl;
     cout << "Z  Pt/Eta/Phi/Mass " << vz.Pt()<< "\t" << vz.Rapidity()<< "\t" << vz.Phi()<< "\t" << vz.M()  << endl;
     cout << "Delta Eta/Phi/deltaR with muon1 : " << aK8JetEta[ij] - vm1.Eta() << "\t" << aK8JetPhi[ij] - vm1.Phi() << "\t" << vj.DeltaR(vm1) << endl;
     cout << "Delta Eta/Phi/deltaR with muon2 : " << aK8JetEta[ij] - vm2.Eta() << "\t" << aK8JetPhi[ij] - vm2.Phi() << "\t" << vj.DeltaR(vm2) << endl;
     }
     if( fabs(vj.DeltaR(vm1)) < 0.8 || fabs(vj.DeltaR(vm2)) < 0.8) continue;
     h_deltaR->Fill(vj.DeltaR(vm1), wnorm);
     h_deltaR->Fill(vj.DeltaR(vm2), wnorm);
     selectedJets.push_back(ij);
       } 
     }
     if(selectedJets.size() ==0) continue;
     ///
     h_preSelnMuons->Fill(nmu.size(), wnorm);
     h_preSelnJets->Fill(selectedJets.size(), wnorm);
     h_preSelLeadMuonPt-> Fill(vm1.Pt(), wnorm);
     h_preSelLeadMuonEta->Fill(vm1.Eta(), wnorm);
     h_preSel2ndMuonPt->  Fill(vm2.Pt(), wnorm);
     h_preSel2ndMuonEta-> Fill(vm2.Eta(), wnorm);
     h_preSelZPt->Fill(vz.Pt(), wnorm);
     h_preSelZEta->Fill(vz.Rapidity(), wnorm);
     h_preSelZPhi->Fill(vz.Phi(), wnorm);
     h_preSelZMass->Fill(vz.M(), wnorm);
     ///for jets
     int fillonce(1), ZtagJets(0);
     for(size_t ijet=0; ijet < selectedJets.size(); ++ijet){
     int iijet = selectedJets[ijet];
     if(aK8JetPrunedMass[iijet] < 40)
     cout << "surviving Jet Pt Eta Phi En " << JetPt[iijet] << "\t" << JetEta[iijet] << "\t" << JetPhi[iijet] << "\t" 
     << JetEn[iijet] << "\t" << aK8JetPrunedMass[iijet] << endl;
     h_preSelJetPt ->Fill(aK8JetPt[iijet], wnorm);
     h_preSelJetEta->Fill(aK8JetEta[iijet], wnorm);
     h_preSelJetPhi->Fill(aK8JetPhi[iijet], wnorm);
     h_preSelJetMass->Fill(aK8JetPrunedMass[iijet], wnorm);
     h_preSelJetSubJet->Fill(aK8Jet_tau3[iijet], wnorm);
     //tagging Z boson
     if(vz.M() > 200){
     if(aK8JetPt[iijet] > 200){
     if(aK8JetPrunedMass[iijet] >= 70 && aK8JetPrunedMass[iijet] <=110){ 
     if(fabs(aK8Jet_tau3[iijet]) < 0.5){
     if(fillonce){
     h_ZTagLeadMuonPt-> Fill(vm1.Pt(), wnorm);
     h_ZTagLeadMuonEta->Fill(vm1.Eta(), wnorm);
     h_ZTag2ndMuonPt->  Fill(vm2.Pt(), wnorm);
     h_ZTag2ndMuonEta-> Fill(vm2.Eta(), wnorm);
     h_ZTagZPt->Fill(vz.Pt(), wnorm);
     h_ZTagZEta->Fill(vz.Rapidity(), wnorm);
     h_ZTagZPhi->Fill(vz.Phi(), wnorm);
     h_ZTagZMass->Fill(vz.M(), wnorm);
     fillonce=0;
     }
     h_ZTagJetPt->Fill(JetPt[iijet], wnorm);
     h_ZTagJetEta->Fill(JetEta[iijet], wnorm);
     h_ZTagJetPhi->Fill(JetPhi[iijet], wnorm);
     h_ZTagJetMass->Fill(aK8JetPrunedMass[iijet], wnorm);
     h_ZTagJetSubJet->Fill(aK8Jet_tau3[iijet], wnorm);
     ++ZtagJets;
            }//subjetiness cut.
          }//pruned mass cut.
        }//ak8JetPt
      }// Zmass cut.
     
     }
     if(ZtagJets) h_ZTagnJets ->Fill(ZtagJets, wnorm);
     if(debug) cout << "**************************************" << endl;  
          }
     f->Write();
     f->Close();
     }

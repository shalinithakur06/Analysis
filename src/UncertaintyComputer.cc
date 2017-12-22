#include "interface/UncertaintyComputer.hh"
#include <iostream>
#include <iomanip>
#include "TRandom3.h"
#include <stack/Roch/RoccoR.cc>

ClassImp(UncertaintyComputer)

using namespace std;
RoccoR  rc("stack/Roch/rcdata.2016.v3");
double UncertaintyComputer::muPtWithRochCorr(const MyMuon *mu, bool isData, double u1, double u2, int s, int m){
  //double genPt  = mu.Genp4.pt();  
  double charge = mu->charge;
  double pt 	= mu->p4.pt();
  double eta 	= mu->p4.eta();
  double phi 	= mu->p4.phi();
  int nl 	= mu->nTrackerLayers;
  double dataSF = rc.kScaleDT(charge, pt, eta, phi, s, m); 
  double mcSF 	= rc.kScaleFromGenMC(charge, pt, eta, phi, nl, u1, u2, s, m); 
  double SF = 1.0; 
  if(isData)SF = dataSF;
  else SF = mcSF;
  return SF*pt;
}

double UncertaintyComputer::metWithJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jer){
  double metX = MET.p4.px();
  double metY = MET.p4.py();
  for(size_t i = 0; i < j->size(); i++){
    int j_ind = j->at(i);
    double gen_pt = vJ[j_ind].Genp4.pt();
    double jet_pt = vJ[j_ind].p4.pt();
    double sigmaJER = vJ[j_ind].resolution;
    //apply JER uncert, scaling
    double delR = DeltaR(vJ[j_ind].Genp4, vJ[j_ind].p4);
    double rCone = 0.4;
    if(gen_pt> 0 && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
      //if(gen_pt <= 0) continue;
      MyLorentzVector rawJet = vJ[j_ind].p4; 
      metX += rawJet.px();
      metY += rawJet.py();
      double jet_pt = vJ[j_ind].p4.pt();
      double SF = getJERSF(vJ[j_ind].p4.eta(), jer);
      double ptscale = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt);
      rawJet *= ptscale;
      metX -= rawJet.px();
      metY -= rawJet.py();
    }
  } 
  return sqrt(metX*metX + metY*metY);
}

double UncertaintyComputer::getJERSF(double eta, int jer){
  double SF = 1.0;
  for(size_t i = 0; i < 13; i++){
    if(TMath::Abs(eta) >= JEREtaMap[i] && TMath::Abs(eta) < JEREtaMap[i+1]){
      if(jer == 0)SF = JERSF[i];
      else if (jer == 1) SF = JERSFUp[i];
      else if(jer == -1) SF = JERSFDown[i];  
    }
  }
  return SF;
}

double UncertaintyComputer::metWithJES(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes){
  double metX = MET.p4.px(); 
  double metY = MET.p4.py(); 
  for(size_t i = 0; i < j->size(); i++){ 
    int j_ind = j->at(i); 
    metX -= (vJ[j_ind].p4.px()*(vJ[j_ind].JECUncertainty*double(jes)));
    metY -= (vJ[j_ind].p4.py()*(vJ[j_ind].JECUncertainty*double(jes)));
  }
  return sqrt(metX*metX + metY*metY);
}

double UncertaintyComputer::metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes, int jer) 
{ 
  double metX = MET.p4.px(); 
  double metY = MET.p4.py(); 
  //get JER uncert.
  for(size_t i = 0; i < j->size(); i++){ 
    int j_ind = j->at(i); 
    double gen_pt = vJ[j_ind].Genp4.pt(); 
    double jet_pt = vJ[j_ind].p4.pt();
    double sigmaJER = vJ[j_ind].resolution;
    //apply JER uncert, scaling
    double delR = DeltaR(vJ[j_ind].Genp4, vJ[j_ind].p4);
    double rCone = 0.4;
    if(gen_pt> 0 && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
    //if(gen_pt <= 0) continue; 
    MyLorentzVector rawJet = vJ[j_ind].p4; 
    metX += rawJet.px(); 
    metY += rawJet.py(); 
    double jet_pt = vJ[j_ind].p4.pt(); 
    double SF = getJERSF(vJ[j_ind].p4.eta(), jer); 
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    double ptscale = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt); 
    rawJet *= ptscale; 
    metX -= rawJet.px(); 
    metY -= rawJet.py();
    } 
  }  
  //get JES unc.
  for(size_t i = 0; i < j->size(); i++){  
    int j_ind = j->at(i);  
    metX -= (vJ[j_ind].p4.px()*(vJ[j_ind].JECUncertainty*double(jes))); 
    metY -= (vJ[j_ind].p4.py()*(vJ[j_ind].JECUncertainty*double(jes))); 
  } 
  return sqrt(metX*metX + metY*metY); 
}

double UncertaintyComputer::metWithUncl(const vector<MyJet> & vJ, vector<int> *j, const vector<MyMuon> &vMu, vector<int> *m, const vector<MyElectron> &vEle, vector<int> *el, MyMET MET, int unc)
{
  double metX = MET.p4.px(); 
  double metY = MET.p4.py(); 

  //remove jets
  for(size_t i = 0; i < j->size(); i++){  
    int j_ind = j->at(i);
    metX += vJ[j_ind].p4.px();
    metY += vJ[j_ind].p4.py();

  }
  //remove leptons
  for(size_t i = 0; i < m->size(); i++){   
    int m_ind = m->at(i); 
    metX += vMu[m_ind].p4.px();
    metY += vMu[m_ind].p4.py();
  }
  for(size_t i = 0; i < el->size(); i++){
    int e_ind = el->at(i);
    metX += vEle[e_ind].p4.px();
    metY += vEle[e_ind].p4.py();
  }

  metX *= (1 + double(unc)*0.1); //vary by 10%
  metY *= (1 + double(unc)*0.1);

  //Re add objects
  for(size_t i = 0; i < m->size(); i++){    
    int m_ind = m->at(i);  
    metX -= vMu[m_ind].p4.px(); 
    metY -= vMu[m_ind].p4.py(); 
  } 
  for(size_t i = 0; i < el->size(); i++){ 
    int e_ind = el->at(i); 
    metX -= vEle[e_ind].p4.px(); 
    metY -= vEle[e_ind].p4.py(); 
  } 
  for(size_t i = 0; i < j->size(); i++){   
    int j_ind = j->at(i); 
    metX -= vJ[j_ind].p4.px(); 
    metY -= vJ[j_ind].p4.py(); 
 
  }
  return sqrt(metX*metX + metY*metY);
}

double UncertaintyComputer::jetPtWithJESJER(MyJet jet, int jes, int jer){

  double gen_pt = jet.Genp4.pt();  
  double jet_pt = jet.p4.pt();
  double sigmaJER = jet.resolution ;
  //apply JES uncert scaling 
  jet_pt *= (1+(jet.JECUncertainty*double(jes)));
  //apply JER uncert, scaling
  double delR = DeltaR(jet.Genp4, jet.p4);
  double rCone = 0.4;
  if(gen_pt> 0 && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
  //if(gen_pt > 0){
    double SF = getJERSF(jet.p4.eta(), jer);
    double ptscale = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt);
    jet_pt *= ptscale;
  }
  return jet_pt;
}

double UncertaintyComputer::EffUncOnSV(MyJet jet)
{
  
  double Uncert = 0.0;
  if(abs(jet.partonFlavour) >= 3 && abs(jet.partonFlavour) <= 5)
    {
      Uncert = sveffunc->getUncC(jet.p4.pt());
    }
  else
    {
      Uncert = sveffunc->getUncL(jet.p4.pt(), jet.p4.eta());
    }

  return Uncert;
}
double UncertaintyComputer::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}

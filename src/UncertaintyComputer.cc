#include "interface/UncertaintyComputer.hh"
#include <iostream>
#include <iomanip>
#include "TRandom3.h"

ClassImp(UncertaintyComputer)

using namespace std;
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

//bottom mistagging, by event re-weighting 
double UncertaintyComputer::getBTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet){
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double pMC = 1.0; 
  pMC = btsf->getBTagPmc(h2_qTagEff_Num, h2_qTagEff_Denom, jet.p4.eta(), jet.p4.pt(), csv);
  return pMC;
}
double UncertaintyComputer::getBTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, int scale){
  double pData = 1.0;
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double eta = jet.p4.eta();
  double pt = jet.p4.pt();
  int flavor = abs(jet.partonFlavour);
  if(scale == 0) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,0);
  else if(scale == 1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,1);
  else if(scale == -1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,-1);
  return pData;
}

double UncertaintyComputer::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}


#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TTimeStamp.h"
#include "Math/VectorUtil.h"

#include "interface/Reader.h"
#include "interface/ObjectSelector.hh"
#include "interface/MomentumVec.h"
#include "interface/LumiReweighting.h"
#include "interface/UncertaintyComputer.hh"
#include "interface/HistogramPlotter.hh"


//---------------------------------------------------//
// Main Class
//---------------------------------------------------//
class Analyzer : public ObjectSelector, HistogramPlotter
{
public :
  Analyzer() : ObjectSelector(), HistogramPlotter()
  {
    DRMIN_JET = 0.4;
    DRMIN_ELE = 0.4;
    //---------------------------------------------------//
    //Pileup reweigting 
    //---------------------------------------------------//
    //PU info for Data: 
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
    //PU info for MC:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information
    LumiWeights_ = reweight::LumiReWeighting("stack/lumiRewgt/trueInTimePU_mcDY.root","stack/lumiRewgt/trueMinBiasPU_dataMu.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    
    //---------------------------------------------------//
    //MC cross sections at 13 TeV 
    //---------------------------------------------------//
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    //evtDBS= event at Data Base Server i.e in DAS (https://cmsweb.cern.ch/das/).
    xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62079400;
    //xss["DY1JetsToLL"]       =  1016;          evtDBS["DY1JetsToLL"]       =  62627174;
    xss["DY2JetsToLL"]       =  331.3;         evtDBS["DY2JetsToLL"]       =  19970551;
    xss["DY3JetsToLL"]       =  96.6;          evtDBS["DY3JetsToLL"]       =  5856110;
    xss["DY4JetsToLL"]       =  51.4;          evtDBS["DY4JetsToLL"]       =  4197868;
    xss["DYJetsToLL"]        =  4895;          evtDBS["DYJetsToLL"]        =  48103700;
    xss["DYJetsToLLamcatnlo"]=  6025.2;        evtDBS["DYJetsToLLamcatnlo"]=  115626000;
    //xss["DYJetsToLLamcatnlo"]=  6025.2;      evtDBS["DYJetsToLLamcatnlo"]=  122055388;
    //xss["DYJetsToLL"]        =  4895;          evtDBS["DYJetsToLL"]        =  49144274;
    xss["QCD_Pt-15to20_Mu"]  =  3819570;       evtDBS["QCD_Pt-15to20_Mu"]  =  4141251;
    xss["QCD_Pt-20to30_Mu"]  =  2960198;       evtDBS["QCD_Pt-20to30_Mu"]  =  31475157;
    xss["QCD_Pt-30to50_Mu"]  =  1652471;       evtDBS["QCD_Pt-30to50_Mu"]  =  29954815;
    xss["QCD_Pt-50to80_Mu"]  =  437504;        evtDBS["QCD_Pt-50to80_Mu"]  =  19806915;
    xss["QCD_Pt-80to120_Mu"] =  106033;        evtDBS["QCD_Pt-80to120_Mu"] =  13786971;
    xss["QCD_Pt-120to170_Mu"]=  25190;         evtDBS["QCD_Pt-120to170_Mu"]=  8042721;
    xss["QCD_Pt-170to300_Mu"]=  8654;          evtDBS["QCD_Pt-170to300_Mu"]=  7947159;
    xss["QCD_Pt-300to470_Mu"]=  797;           evtDBS["QCD_Pt-300to470_Mu"]=  7937590;
    xss["QCD_Pt-470to600_Mu"]=  79 ;           evtDBS["QCD_Pt-470to600_Mu"]=  18976018;
    xss["QCD_Pt-600to800_Mu"]=  25.09;         evtDBS["QCD_Pt-600to800_Mu"]=  9981311;
    xss["QCD_Pt-800to1000_Mu"]= 4.70;          evtDBS["QCD_Pt-800to1000_Mu"]=  19767439;
    xss["QCD_Pt-1000toInf_Mu"]= 1.62;          evtDBS["QCD_Pt-1000toInf_Mu"]=  13400031;

    xss["QCD_Pt-15to20_EM"]  =  254600;        evtDBS["QCD_Pt-15to20_EM"]  =  5652601;
    xss["QCD_Pt-20to30_EM"]  =  5352960;       evtDBS["QCD_Pt-20to30_EM"]  =  9218954;
    xss["QCD_Pt-30to50_EM"]  =  9928000;       evtDBS["QCD_Pt-30to50_EM"]  =  4730195;
    xss["QCD_Pt-50to80_EM"]  =  2890800;       evtDBS["QCD_Pt-50to80_EM"]  =  22337070;
    xss["QCD_Pt-80to120_EM"] =  350000;        evtDBS["QCD_Pt-80to120_EM"] =  35841783;
    xss["QCD_Pt-120to170_EM"]=  62964;         evtDBS["QCD_Pt-120to170_EM"]=  35817281;
    xss["QCD_Pt-170to300_EM"]=  18810;         evtDBS["QCD_Pt-170to300_EM"]=  11540163;
    xss["QCD_Pt-300toInf_EM"]=  1350;          evtDBS["QCD_Pt-300toInf_EM"]=  7373633;

    xss["TTJetsP"]           =  831.76;        evtDBS["TTJetsP"]           =  77081156;   
    xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  44813600;
    //xss["W1JetsToLNu"]       =  9493;          evtDBS["W1JetsToLNu"]       =  45367044;
    xss["W2JetsToLNu"]       =  3120;          evtDBS["W2JetsToLNu"]       =  29878415;
    xss["W3JetsToLNu"]       =  942.3;         evtDBS["W3JetsToLNu"]       =  19798117;
    xss["W4JetsToLNu"]       =  524.2;         evtDBS["W4JetsToLNu"]       =  9170576;
    xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29181900;
    //xss["WJetsToLNu"]        =  50690;         evtDBS["WJetsToLNu"]        =  29705748;
    xss["WW"]                =  118.7;         evtDBS["WW"]                =  994012;
    xss["WZ"]                =  46.74;         evtDBS["WZ"]                =  1000000;
    xss["ZZ"]                =  17.72;         evtDBS["ZZ"]                =  990064; 
    //xss["sampCode_"]         =  1;           evtDBS["sampCode_"]         =  1; 
    xss["sampCode_"]         =   0.00427;      evtDBS["sampCode_"]         =  10000; 
    //xss["Signal sample"]
    //xss["TestGenSim_Mu_250"]        =  0.00427;   evtDBS["TestGenSim_Mu_250"]        = 10000;
    xss["ExcitedLepton_MuMuZ_250"]  =  0.00427;   evtDBS["ExcitedLepton_MuMuZ_250"]  = 198200;
    xss["ExcitedLepton_MuMuZ_1500"]  = 0.00427;   evtDBS["ExcitedLepton_MuMuZ_1500"] = 190900;
    xss["ExcitedLepton_MuMuZ_2000"]  = 0.00427;   evtDBS["ExcitedLepton_MuMuZ_2000"] = 200000;
    xss["ExcitedLepton_MuMuZ_4000"]  = 0.00427;   evtDBS["ExcitedLepton_MuMuZ_4000"] = 192300;
    //xss["TestLHEGeneration_Mu_6000"]    =  0.1662;      evtDBS["TestLHEGeneration_Mu_6000"] =  10000;

    //xss["TestLHEGeneration_Ele_2000"]    =  0.1662;   evtDBS["TestLHEGeneration_Ele_2000"]  =  9280;
    //xss["TestLHEGeneration_Ele_4000"]    =  0.1662;   evtDBS["TestLHEGeneration_Ele_4000"]  =  7536;
    //xss["TestLHEGeneration_Ele_6000"]    =  0.1662;   evtDBS["TestLHEGeneration_Ele_6000"]  =  8848;

    //Lumis(inverse pb) of single muon DATA at 13TeV
  };
  ~Analyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url,  string myKey="PFlow", string evtType="data");
  void CutFlowProcessor(TString url,  string myKey="PFlow", TString cutflowType="base", TFile *outFile_=0);
  //void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();
  float reweightHEPNUPWJets(int hepNUP);
  float reweightHEPNUPDYJets(int hepNUP);
private :
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_;   //pileup syst up
  reweight::PoissonMeanShifter PShiftDown_; //pileup syst down 
  std::map<string, double> xss;
  std::map<string, double> evtDBS;
  std::map<string, double> muSF;
  std::map<string, double> eleSF;
  /*
  BTagCalibrationReader readCSVfile(const std::string &tagger, const std::string &filename);
  */
  ofstream outfile_;
  Double_t getMuonSF(TH2D *h2, double eta, double pt);
  Double_t getMuonTrigSF(TH2D *h2, double eta, double pt);
  Double_t getMuonTrackSF(TGraphAsymmErrors *tg, double eta);
  Double_t getEleSF(TH2D *h2, double etaSC, double pt);
  Double_t getEleTrigSF(TH2D *h2, double pt, double etaSC);
  Double_t getEleHeep1SF(TGraphAsymmErrors *tg, double pt);
  Double_t getEleHeep2SF(TGraphAsymmErrors *tg, double eta);
  double deltaPhi12(double phi1, double phi2);
  double phi0to2pi(double phi);
};

float Analyzer::reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  if(nJets==0)      return 2.11;
  else if(nJets==1) return 0.23;
  else if(nJets==2) return 0.119;
  else if(nJets==3) return 0.0562;
  else if(nJets>=4) return 0.0671;
  else return 1 ;
}

float Analyzer::reweightHEPNUPDYJets(int hepNUP){

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.120;
  else if(nJets==1) return 0.0164;
  else if(nJets==2) return 0.0168;
  else if(nJets==3) return 0.0167;
  else if(nJets>=4) return 0.0128;
  else return 1 ;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
Double_t Analyzer::getMuonSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=120){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(120);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}
Double_t Analyzer::getMuonTrigSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}

Double_t Analyzer::getMuonTrackSF(TGraphAsymmErrors *tg, double eta){
 
  Double_t *eta_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();

  double SF = 1.0;
  // eta < eta_array[0]
  if(abs(eta)<eta_array[0]) SF = sf_array[0];
  
  // eta_array[0]<eta<eta_array[n_points -1]
  for(Int_t i = 0; i < n_points-1; i++){
    if(abs(eta) >= eta_array[i] && abs(eta) < eta_array[i+1]) SF = sf_array[i+1];
  }
  // eta > eta_array[n_points -]
  if(abs(eta)>eta_array[n_points-1]) SF = sf_array[n_points -1];
  return SF;
}

Double_t Analyzer::getEleSF(TH2D *h2, double etaSC, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <500
  //for Pt >500, we use SF of Pt = 500
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}
Double_t Analyzer::getEleHeep1SF(TGraphAsymmErrors *tg, double pt){
  Double_t *pt_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();
  double SF = 1.0;
if(pt < pt_array[0]) SF = sf_array[0];
for(Int_t i = 0; i < n_points-1; i++){
  if(pt >= pt_array[i] && pt < pt_array[i+1]) SF = sf_array[i+1];
  }
  if(pt > pt_array[n_points-1]) SF = sf_array[n_points-1];
  return SF;
}

Double_t Analyzer::getEleHeep2SF(TGraphAsymmErrors *tg, double etaSC){

  Double_t *eta_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();

double SF = 1.0;
if(abs(etaSC)<eta_array[0]) SF = sf_array[0];
for(Int_t i = 0; i < n_points-1; i++){
  if(abs(etaSC) >= eta_array[i] && abs(etaSC) < eta_array[i+1]) SF = sf_array[i+1];
  }
  if(abs(etaSC)>eta_array[n_points-1]) SF = sf_array[n_points -1];
return SF;
}

Double_t Analyzer::getEleTrigSF(TH2D *h2, double pt, double etaSC){
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=200){
    Int_t binX = xaxis->FindBin(pt);
    Int_t binY = yaxis->FindBin(etaSC);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(200);
    Int_t binY = yaxis->FindBin(etaSC);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}

double phi0to2pi(double phi){
    double pi = 3.1415926535;
    while (phi >= 2.*pi) phi -= 2.*pi;
    while (phi < 0.) phi += 2.*pi;
    return phi;
}

double deltaPhi12(double phi1_, double phi2_){
    // build the delta Phi angle between the two vectors
    double pi = 3.1415926535;
    double phi1 = phi0to2pi(phi1_);
    double phi2 = phi0to2pi(phi2_);
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2*pi - dPhi)) ? 2*pi - dPhi : dPhi;
    return dPhi;
}

//---------------------------------------------------//
//muon scale factors from 2D histograms 
//---------------------------------------------------//      
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2 
//Trigger SF
TFile *f_trigSF_BCDEF 	= new TFile("stack/muonSF/triggreSF_BCDEF.root");
TFile *f_trigSF_GH 		= new TFile("stack/muonSF/triggreSF_GH.root");
TH2D *h2_trigSF_BCDEF 	= (TH2D*)f_trigSF_BCDEF->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
TH2D *h2_trigSF_GH 		= (TH2D*)f_trigSF_GH->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
//Identification SF
TFile *f_idSF_BCDEF 		= new TFile("stack/muonSF/idSF_BCDEF.root");
TFile *f_idSF_GH 		= new TFile("stack/muonSF/idSF_GH.root");
TH2D *h2_idSF_BCDEF 		= (TH2D*)f_idSF_BCDEF->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
TH2D *h2_idSF_GH 		= (TH2D*)f_idSF_GH->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
//Isolation SF
TFile *f_isoSF_BCDEF 		= new TFile("stack/muonSF/isoSF_BCDEF.root");
TFile *f_isoSF_GH 		= new TFile("stack/muonSF/isoSF_GH.root");
TH2D *h2_isoSF_BCDEF 		= (TH2D*)f_isoSF_BCDEF->Get("TightISO_MediumID_pt_eta/abseta_pt_ratio");
TH2D *h2_isoSF_GH 		= (TH2D*)f_isoSF_GH->Get("TightISO_MediumID_pt_eta/abseta_pt_ratio");
//Tracking SF
TFile *f_trackSF_BCDEF 	= new TFile("stack/muonSF/trackingSF_BCDEF.root");
TFile *f_trackSF_GH 		= new TFile("stack/muonSF/trackingSF_GH.root");
TGraphAsymmErrors *tg_trackSF_BCDEF 	= (TGraphAsymmErrors*)f_trackSF_BCDEF->Get("ratio_eff_aeta_dr030e030_corr");
TGraphAsymmErrors *tg_trackSF_GH 	= (TGraphAsymmErrors*)f_trackSF_GH->Get("ratio_eff_aeta_dr030e030_corr");

//---------------------------------------------------//
//Electron scale factors from 2D histograms 
//---------------------------------------------------//      
//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
//Reconstruction SF
TFile *f_ele_recoSF 	  	= new TFile("stack/eleSF/ele_recoSF.root");
TH2D *h2_ele_recoSF 		= (TH2D*)f_ele_recoSF->Get("EGamma_SF2D");
//Identification SF
TFile *f_ele_veto_idSF 	= new TFile("stack/eleSF/ele_veto_idSF.root");
TH2D *h2_ele_veto_idSF 	= (TH2D*)f_ele_veto_idSF->Get("EGamma_SF2D");
TFile *f_ele_loose_idSF 	= new TFile("stack/eleSF/ele_loose_idSF.root");
TH2D *h2_ele_loose_idSF 	= (TH2D*)f_ele_loose_idSF->Get("EGamma_SF2D");
TFile *f_ele_medium_idSF 	= new TFile("stack/eleSF/ele_medium_idSF.root");
TH2D *h2_ele_medium_idSF 	= (TH2D*)f_ele_medium_idSF->Get("EGamma_SF2D");
TFile *f_ele_tight_idSF 	= new TFile("stack/eleSF/ele_tight_idSF.root");
TH2D *h2_ele_tight_idSF 	= (TH2D*)f_ele_tight_idSF->Get("EGamma_SF2D");
TFile *f_ele_heep_SF_EE     = new TFile("stack/eleSF/HEEP_SF.root");
TFile *f_ele_heep_SF_EB     = new TFile("stack/eleSF/HEEP_SF.root");
TFile *f_ele_heep_SF            = new TFile("stack/eleSF/HEEP_SF.root");
TGraphAsymmErrors *tg_heep_SF_EE      = (TGraphAsymmErrors*)f_ele_heep_SF_EE->Get("SF_Et_Endcap");
TGraphAsymmErrors *tg_heep_SF_EB      = (TGraphAsymmErrors*)f_ele_heep_SF_EB->Get("SF_Et_Barrel");
TGraphAsymmErrors *tg_heep_SF             = (TGraphAsymmErrors*)f_ele_heep_SF->Get("SF_Eta");

//Trigger scale factors
//https://indico.cern.ch/event/604912/
//TFile *f_ele_trigSF 		= new TFile("stack/eleSF/ele_trigSF_Run2016All_v1.root");
//TH2D *h2_ele_trigSF 		= (TH2D*)f_ele_trigSF->Get("Ele27_WPTight_Gsf");

TFile *f_ele_trigSF 		= new TFile("stack/eleSF/MW_2ndleg_EGM2D.root");
TH2D *h2_ele_trigSF 		= (TH2D*)f_ele_trigSF->Get("EGamma_SF2D");

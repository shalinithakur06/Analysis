
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
    //Backgrounds
    //only 97.3% proccessed, hence less number of events for DYJetsToLL_M50
    //xss["DYJetsToLL_M50"]         = 5758.4;              evtDBS["DYJetsToLL_M50"]         = 122055388; 
    xss["DYJetsToLL_M50"]         = 5758.4;              evtDBS["DYJetsToLL_M50"]         = 0.6648* 118692000; 
    xss["DYJetsToLL_M100to200"]   = 226.6 ;              evtDBS["DYJetsToLL_M100to200"]   = 0.6439* 1083606;
    xss["DYJetsToLL_M200to400"]   = 7.77 ;               evtDBS["DYJetsToLL_M200to400"]   = 0.5688* 2925885;    
    xss["DYJetsToLL_M400to500"]   = 0.423;               evtDBS["DYJetsToLL_M400to500"]   = 0.5208* 287262;
    xss["DYJetsToLL_M500to700"]   = 0.240;               evtDBS["DYJetsToLL_M500to700"]   = 0.5099* 280940;  
    xss["DYJetsToLL_M700to800"]   = 0.0350;              evtDBS["DYJetsToLL_M700to800"]   = 276235;//no neg weight
    xss["DYJetsToLL_M800to1000"]  = 0.0300;              evtDBS["DYJetsToLL_M800to1000"]  = 271768;//no neg weight
    xss["DYJetsToLL_M1000to1500"] = 0.0160;              evtDBS["DYJetsToLL_M1000to1500"] = 258620;//no neg weight
    xss["DYJetsToLL_M1500to2000"] = 0.0020;              evtDBS["DYJetsToLL_M1500to2000"] = 258625;//no neg weight
    xss["DYJetsToLL_M2000to3000"] = 0.00054;             evtDBS["DYJetsToLL_M2000to3000"] = 255342;//no neg weight
    xss["TT"] 			  = 831.76;              evtDBS["TT"] 		          = 77081156;
    //only 97.7% proccessed, hence less number of events for WJetsToLNu
    //xss["WJetsToLNu"]		  = 50690;               evtDBS["WJetsToLNu"]	          = 29705748;
    xss["WJetsToLNu"]		  = 50690;               evtDBS["WJetsToLNu"]	          = 29017600;
    xss["WW"]			  = 118.7;               evtDBS["WW"]		          = 994012;
    xss["WZ"]			  = 46.74;	         evtDBS["WZ"]		          = 1000000;
    xss["ZZ"] 			  = 17.72;               evtDBS["ZZ"] 		          = 990064;
    //Muon signal
    xss["ExLepMuMuZ_M250"]	  = 0.00427;             evtDBS["ExLepMuMuZ_M250"]	  = 198200;
    xss["ExLepMuMuZ_M1500"]	  = 0.0004267;           evtDBS["ExLepMuMuZ_M1500"]	  = 190900;
    xss["ExLepMuMuZ_M2000"] 	  = 0.00002021;          evtDBS["ExLepMuMuZ_M2000"]       = 200000;
    xss["ExLepMuMuZ_M2500"] 	  = 0.00006755;          evtDBS["ExLepMuMuZ_M2500"]       = 200000;
    xss["ExLepMuMuZ_M4000"]       = 0.000003209;         evtDBS["ExLepMuMuZ_M4000"]       = 192300;
    //Electron signal
    xss["ExLepEEZ_M250"] 	  = 0.004407;            evtDBS["ExLepEEZ_M250"] 	  = 182100;
    xss["ExLepEEZ_M2000"] 	  = 0.0001908;	         evtDBS["ExLepEEZ_M2000"] 	  = 189000;
    xss["ExLepEEZ_M2500"]	  = 0.00005731;          evtDBS["ExLepEEZ_M2500"]	  = 187400;
    xss["ExLepEEZ_M4000"]	  = 0.000003897;         evtDBS["ExLepEEZ_M4000"]	  = 199800;
    //Dummy sample
    xss["sampCode_"]    	  =  1;  	         evtDBS["sampCode_"]              =  1; 

  };
  ~Analyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url,  string myKey="PFlow", string evtType="data");
  void CutFlowProcessor(TString url,  string myKey="PFlow", TString cutflowType="base", TFile *outFile_=0);
  //void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();
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

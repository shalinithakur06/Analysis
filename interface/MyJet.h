#ifndef __MYJET_H__
#define __MYJET_H__

#include "TROOT.h"
#include <map>
#include <vector>
#include <string>

#include "MomentumVec.h"

class MyJet 
{
 public:
  MyJet();
  ~MyJet();
  void   Reset();
  
  std::string GetName() { 
      std::string name(jetName);
      return name; 
  };
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016 
 
  ///basic
  MyLorentzVector Genp4;      //4 vector of matched genJet
  double jetCharge;
  ///std::string jetName;
  TString jetName;
  MyLorentzVector p4;         // 4 vector of jet.
  int parton_id;
  int parton_mother_id;
  int partonFlavour;
  int hadronFlavour;
  Point3D vertex;
  
  ///ids 
  double neutralHadronEnergyFraction;   
  double neutralEmEnergyFraction; 
  int NumConst; //Number of Constituents
  double muonEnergyFraction; //Muon Fraction
  double chargedHadronEnergyFraction;
  int chargedMultiplicity;
  double chargedEmEnergyFraction; 
  int neutralMultiplicity; //Number of Neutral Particles
  bool jetIDLoose; 
 
  ///Btag, JEC & SV
  std::map<std::string, double>bDiscriminator;
  std::map<std::string, double>JECs;
  double JECUncertainty;
  std::vector<MyLorentzVector>SVP4; 
  std::vector<double>SVflightDistance;
  std::vector<double>SVflightDistanceErr;
  std::vector<float>SVNChi2;
  
  //Jet resolution and scale factors
  double scaleFactor;
  double resolution;
  double triggerJet_pt;
  int quality;

 private :

}; 
#endif

#include "interface/MyJet.h"

MyJet::MyJet():
  
  ///basic  
  jetCharge(0),
  jetName(""),
  parton_id(0),
  parton_mother_id(0),
  partonFlavour(0),
  hadronFlavour(0),

  ///ids
  neutralHadronEnergyFraction(-9.),
  neutralEmEnergyFraction(-9.),
  NumConst(0.),
  muonEnergyFraction(-9.),
  chargedHadronEnergyFraction(-9.),
  chargedMultiplicity(0.),
  chargedEmEnergyFraction(-9.),
  neutralMultiplicity(0.),
  jetIDLoose(false),
  
  ///JEC
  JECUncertainty(-9.),
  scaleFactor(0.0),
  resolution(0.0),
  triggerJet_pt(-9.),
  quality(0)
{
}

MyJet::~MyJet() 
{
}

void MyJet::Reset()
{
  ///basic
  Genp4.SetCoordinates(0.0, 0.0, 0.0, 0.0);
  jetCharge = 0;
  jetName = "";
  p4.SetCoordinates(0.0, 0.0, 0.0, 0.0);
  parton_id = 0;
  parton_mother_id = 0;
  partonFlavour = 0;
  hadronFlavour = 0;
  vertex.SetCoordinates(-9.0,-9.0,-9.0);

  ///ids
  neutralHadronEnergyFraction = -9.;
  neutralEmEnergyFraction = -9.;
  NumConst = 0.;
  muonEnergyFraction = -9.;
  chargedHadronEnergyFraction = -9.;
  chargedMultiplicity = 0.;
  chargedEmEnergyFraction = -9.;
  neutralMultiplicity = 0.;
  jetIDLoose = false;
  
  ///btag, JEC, SV
  bDiscriminator.clear();
  JECs.clear();
  JECUncertainty = -9.;
  SVP4.clear();
  SVflightDistance.clear();
  SVflightDistanceErr.clear();
  SVNChi2.clear();

  scaleFactor = 0.0;
  resolution = 0.0;
  triggerJet_pt = 0.;
  quality = 0;

}

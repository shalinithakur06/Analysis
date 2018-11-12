#include "interface/MyElectron.h"

MyElectron::MyElectron():

  ///basic
  charge(0), 
  gen_id(0),
  gen_mother_id(0),
  //name(""),
  eleName(""),
  eleSCEta(0),

  ///sel
  sigmaIetaIeta(-999),
  dEtaInSeed(-999),
  dPhiIn(-999),
  hadOverEm(-999),
  iEminusiP(-999),
  nInnerHits(-999),
  nInnerLostHits(-999),
  isPassConVeto(true),
  isEcalDriven(true),
  energy5x5(-999),
  energy2x5(-999),
  eleRho(-999),
  eleTrkPt(-999),
  GsfEleEmHadD1IsoRhoCut(-999),

  ///ids
  isEE(-1),
  isEB(-1),
  
  ///iso
  ChHadIso(999.),  
  PhotonIso(999.),   
  NeuHadIso(999.),   
  PileupIso(999.), 
  relCombPFIsoEA(999.),
  D0(999.),
  Dz(999.),
  trigger_ele_pt(0),
  quality(0),
  passConversionVeto(true)
{
}

MyElectron::~MyElectron()
{
}


void MyElectron::Reset()
{
  ///basic
  charge = 0; 
  gen_id = 0;
  gen_mother_id = 0;
 // name = "";
  eleName = "";
  p4.SetCoordinates(0.0, 0.0, 0.0, 0.0);
  eleSCEta = 0;
  vertex.SetCoordinates(0.0, 0.0, 0.0);
 
  ///sel
  sigmaIetaIeta = -999.;
  dEtaInSeed = -999.;
  dPhiIn = -999.;
  hadOverEm = -999.;
  iEminusiP = -999.;
  nInnerHits = -999.;
  nInnerLostHits=-999;
  isPassConVeto = true;
  isEcalDriven = true;
  energy5x5 = -999;
  energy2x5 = -999;
  eleRho = -999;
  eleTrkPt= -999;
  GsfEleEmHadD1IsoRhoCut= -999;

  ///ids
  isEE = -1;
  isEB = -1;
  eidWPs.clear();
  
  ///iso
  ChHadIso = 999.;   
  PhotonIso = 999.;    
  NeuHadIso = 999.;    
  PileupIso = 999.;  
  relCombPFIsoEA = 999.; 
  D0 = 999.;
  Dz = 999.;
  trigger_ele_pt = 0;
  passConversionVeto = true;
  quality = 0;
  }

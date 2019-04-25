#include "interface/ObjectSelector.hh"
#include <iostream>
#include <iomanip>
ClassImp(ObjectSelector)

using namespace std;

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
//Electron ID: veto
bool ObjectSelector::cutBasedElectronID_Summer16_80X_V1_veto(const MyElectron *e)
{
  bool passID = false;
  //barrel
  if(abs(e->eleSCEta) <= 1.479 
     && e->sigmaIetaIeta 	< 0.0115	 
     && abs(e->dEtaInSeed) 	< 0.00749	
     && abs(e->dPhiIn) 		< 0.228	
     && e->hadOverEm 		< 0.356	
     && e->relCombPFIsoEA 	< 0.175	
     && abs(e->iEminusiP) 	< 0.299	
     && e->nInnerHits       	<= 2
     && e->passConversionVeto  
    )passID = true;

  //endcap 
  if(abs(e->eleSCEta) > 1.479 
     && e->sigmaIetaIeta 	< 0.037	
     && abs(e->dEtaInSeed) 	< 0.00895	 
     && abs(e->dPhiIn) 		< 0.213	
     && e->hadOverEm 		< 0.211	
     && e->relCombPFIsoEA 	< 0.159	
     && abs(e->iEminusiP) 	< 0.15	
     && e->nInnerHits       	<= 3
     && e->passConversionVeto  
     )passID = true;
  return passID;
}

//Electron ID: loose
bool ObjectSelector::cutBasedElectronID_Summer16_80X_V1_loose(const MyElectron *e)
{
  bool passID = false;
  //barrel
  if(abs(e->eleSCEta) <= 1.479 
     && e->sigmaIetaIeta 	< 0.011	 
     && abs(e->dEtaInSeed) 	< 0.00477	
     && abs(e->dPhiIn) 		< 0.222	
     && e->hadOverEm 		< 0.298	
     && e->relCombPFIsoEA 	< 0.0994	
     && abs(e->iEminusiP) 	< 0.241	
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
    )passID = true;

  //endcap 
  if(abs(e->eleSCEta) > 1.479 
     && e->sigmaIetaIeta 	< 0.0314	
     && abs(e->dEtaInSeed) 	< 0.00868	 
     && abs(e->dPhiIn) 		< 0.213	
     && e->hadOverEm 		< 0.101	
     && e->relCombPFIsoEA 	< 0.107	
     && abs(e->iEminusiP) 	< 0.14	
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
     )passID = true;
  return passID;
}

//Electron ID: medium
bool ObjectSelector::cutBasedElectronID_Summer16_80X_V1_medium(const MyElectron *e)
{
  bool passID = false;
  //barrel
  if(abs(e->eleSCEta) <= 1.479 
     && e->sigmaIetaIeta 	< 0.00998	 
     && abs(e->dEtaInSeed) 	< 0.00311	
     && abs(e->dPhiIn) 		< 0.103	
     && e->hadOverEm 		< 0.253	
     ///&& e->relCombPFIsoEA 	< 0.0695	
     && abs(e->iEminusiP) 	< 0.134	
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
    )passID = true;

  //endcap
  if(abs(e->eleSCEta) > 1.479 
     && e->sigmaIetaIeta 	< 0.0298	
     && abs(e->dEtaInSeed) 	< 0.00609	 
     && abs(e->dPhiIn) 		< 0.045	
     && e->hadOverEm 		< 0.0878	
     ///&& e->relCombPFIsoEA 	< 0.0821	
     && abs(e->iEminusiP) 	< 0.13	
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
     )passID = true;
  return passID;
}

//Electron ID: tight
bool ObjectSelector::cutBasedElectronID_Summer16_80X_V1_tight(const MyElectron *e)
{
  bool passID = false;
  //barrel
  if(abs(e->eleSCEta) <= 1.479 
     && e->sigmaIetaIeta 	< 0.00998 
     //&& abs(e->dEtaInSeed) 	< 0.00308
     && abs(e->dPhiIn) 		< 0.0816 
     && e->hadOverEm 		< 0.0414 
     && e->relCombPFIsoEA 	< 0.0588
     && abs(e->iEminusiP) 	< 0.0129 
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
    )passID = true;
  //endcap
  if(abs(e->eleSCEta) > 1.479 
     && e->sigmaIetaIeta 	< 0.0292 
     && abs(e->dEtaInSeed) 	< 0.00605 
     && abs(e->dPhiIn) 		< 0.0394
     && e->hadOverEm 		< 0.0641
     && e->relCombPFIsoEA 	< 0.0571
     && abs(e->iEminusiP) 	< 0.0129
     && e->nInnerHits       	<= 1
     && e->passConversionVeto  
     )passID = true;
  return passID;
}
//Electron HEEP ID
  bool ObjectSelector::heepElectronID_HEEPV70(const MyElectron *e, MyVertex & vertex)
{
  bool passID = false;
  float energy2x5Overenergy5x5 = e->energy2x5/e->energy5x5;
  //for barrel
  if(abs(e->eleSCEta)                <=1.444
     && e->isEcalDriven              == 1 
     && abs(e->dEtaInSeed)           < 0.004 
     && abs(e->dPhiIn)               < 0.06 
     && e->hadOverEm                 < 1/e->p4.E() + 0.05 
     && energy2x5Overenergy5x5       > 0.94 
     && e->GsfEleEmHadD1IsoRhoCut    < 2+0.03*e->p4.pt()+0.28*e->eleRho 
     && e->eleTrkPt                  < 5  
     && e->nInnerHits                <=1   //Inner Lost Hits
     && abs(e->D0)                   < 0.02
     )passID = true;
  //endcap
  double HadDepth = 0.0;
  if(e->p4.E() < 50) HadDepth = 2.5 +0.28*e->eleRho;
  else HadDepth = 2.5+0.03*(e->p4.E()-50) +0.28*e->eleRho; 
  if(abs(e->eleSCEta) > 1.566 && abs(e->eleSCEta) < 2.5
     && e->isEcalDriven                == 1
     && abs(e->dEtaInSeed)             < 0.006
     && abs(e->dPhiIn)                 < 0.06 
     && e->hadOverEm                   < 5/e->p4.E() + 0.05 
     && e->sigmaIetaIeta               <0.03 
     && e->GsfEleEmHadD1IsoRhoCut      < HadDepth
     && e->eleTrkPt                    < 5
     && e->nInnerHits                  <=1
     && abs(e->D0)                     < 0.05
     )passID = true;
  return passID; 
}

void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  for(unsigned int i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    double ePt     	   = TMath::Abs(e->p4.pt());
    double eEta     	   = TMath::Abs(e->p4.eta());
    //bool passID = cutBasedElectronID_Summer16_80X_V1_loose(e);
    ///bool passID = cutBasedElectronID_Summer16_80X_V1_medium(e);
    //bool passID = cutBasedElectronID_Summer16_80X_V1_tight(e);
    bool passID = heepElectronID_HEEPV70(e, vertex); 
    if(passID && ePt >35 && eEta <2.5 ){e_i->push_back(i);}
  }
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
bool ObjectSelector::isHighPtMuon(const MyMuon * m){
  bool isHighPt(false);
  isHighPt = m->isGlobalMuon &&
	  m->nMuonHits > 0 && 
	  m->nMatchedStations >1 &&
          m->bestMuPtErr/m->bestMuPtTrack < 0.3 &&
          m->nPixelHits > 0 &&
          m->nTrackerLayers > 5;
return isHighPt;
}
void ObjectSelector::preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isData){
  for( int i=0;i< (int) vM.size();i++){
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt = TMath::Abs(m->p4.pt());
    bool passID = isHighPtMuon(m);
    double iso = m->pfRelIso;
    if(passID && mPt > 35 && mEta < 2.4 && 
		    abs(m->D0) < 0.2 && abs(m->Dz) < 0.5 && iso < 0.15){ 
      m_i->push_back(i);
    }
  }
}

void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes, int jer){
  for(unsigned int i=0;i<vJ.size();i++){
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    double jetPt      = jetPtWithJESJER(vJ[i], jes, jer); 
    double neutralHadEnFrac = jet->neutralHadronEnergyFraction;
    double neutralEmEnFrac = jet->neutralEmEnergyFraction;
    double chargedHadEnFrac = jet->chargedHadronEnergyFraction;
    double chargedEmFrac = jet->chargedEmEnergyFraction;
    if(jetEta < 2.4){
      if(jetPt > 170 
        && neutralHadEnFrac < 0.90
        && neutralEmEnFrac  < 0.90
	&& jet->NumConst > 1
        && chargedHadEnFrac > 0
	&& chargedEmFrac < 0.99
	&& jet->chargedMultiplicity > 0){
        j_i->push_back(i);
      }
    }
    if(jetEta >= 2.4 && jetEta < 2.5){
      if(jetPt > 170 
        && neutralHadEnFrac < 0.90
        && neutralEmEnFrac  < 0.90
	&& jet->NumConst > 1){
        j_i->push_back(i);
      }
    }
  }
}


bool ObjectSelector::looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM, bool isPFlow){
  bool looseVeto(false);
  for(int i=0;i< (int)vM.size();i++){
    if( i==selectedMuon ){continue;}
    const MyMuon * m = &vM[i];
    bool isGlobalMuon = m->isGlobalMuon; 
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mRelIso  = m->pfRelIso;
    if(! isGlobalMuon) continue;
    if(mEta<2.4  && mPt> 15 && mRelIso < 0.20 ){ looseVeto = true; }
  }
  return looseVeto;
    
}

bool ObjectSelector::looseElectronVeto(unsigned long selectedElectron, const vector<MyElectron> & vE, MyVertex & vertex, bool isPFlow){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  bool looseVeto(false);
  for(unsigned long i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    //double eEta     	   = TMath::Abs(e->p4.eta());
    double ePt     	   = TMath::Abs(e->p4.pt());
    double d0      	   = fabs(e->D0);
    double zvertex   	   = vertex.XYZ.z();
    double zelectron 	   = e->vertex.z();
    double dz 		   = fabs(zvertex - zelectron);
    if( i==selectedElectron) continue; 
    bool passID = cutBasedElectronID_Summer16_80X_V1_veto(e);
    if(passID && ePt >15 && d0 < 0.05 && dz < 0.1){looseVeto = true;}
  }
  return looseVeto;
}


void ObjectSelector::ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR ){
  for(size_t i = 0; i < e_old->size(); i++){
    int iele = (*e_old)[i];
    double delR2Mu = 5.0;
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vE[iele].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }
    if(delR2Mu > DR) e_new->push_back(iele);
  }
}

void ObjectSelector::JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE,vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, double DR){
  for(size_t i = 0; i < j_old->size(); i++){
    int ijet = (*j_old)[i];

    double delR2Mu = 5.0, delR2Ele = 5.0;
    
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vJ[ijet].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }
    for(size_t k = 0; k < el->size(); k++){
      int iele = (*el)[k];
      double delR = DeltaR(vJ[ijet].p4, vE[iele].p4);
      if(delR < delR2Ele)delR2Ele = delR;
    }
    if(delR2Mu > DR && delR2Ele > DR )
    {
        j_new->push_back(ijet);
    }
    }
}

double ObjectSelector::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}


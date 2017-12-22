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
     //&& e->relCombPFIsoEA 	< 0.175	
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
     //&& e->relCombPFIsoEA 	< 0.159	
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
     && abs(e->dEtaInSeed) 	< 0.00308
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

void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  for(unsigned int i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];

    double ePt     	   = TMath::Abs(e->p4.pt());
    double eEta     	   = TMath::Abs(e->p4.eta());
    double d0      	   = fabs(e->D0);
    double zvertex   	   = vertex.XYZ.z();
    double zelectron 	   = e->vertex.z();
    double dz 		   = fabs(zvertex - zelectron);
    //bool passID = cutBasedElectronID_Summer16_80X_V1_loose(e);
    bool passID = cutBasedElectronID_Summer16_80X_V1_medium(e);
    //bool passID = cutBasedElectronID_Summer16_80X_V1_tight(e);
    
    if(passID && ePt >30 && eEta <2.5 && d0 < 0.05 && dz < 0.2){e_i->push_back(i);}
  }
}

void ObjectSelector::preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isData, double random_u1, double random_u2, int err_member, int err_set){
  for( int i=0;i< (int) vM.size();i++){
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    ///double mPt      = TMath::Abs(m->p4.pt());
    double mD0      = fabs(m->D0);
    //double mRelIso  = m->pfRelIso;
    double mPt   = muPtWithRochCorr(m, isData, random_u1, random_u2, err_set, err_member); 
    /*
    bool isGlobalMuon = m->isGlobalMuon; 
    bool isPFMuon = m->isPFMuon; 
    bool passId = (isGlobalMuon && isPFMuon && m->nMuonHits >=1 
		   && m->nPixelHits >= 1 && m->nMatchedStations >=2 
		   && m->nTrackerLayers >= 6 && m->normChi2 < 10); 
    */
    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ && mD0 < M_D0_MAX_&& dz < ZMAX_){ 
    //if(passId && mPt > 20 && mEta < 2.1 && mRelIso <0.50){ 
      m_i->push_back(i);
    }
  }
}

void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes, int jer){
 
  for(unsigned int i=0;i<vJ.size();i++){
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    //double jetPt      = TMath::Abs(jet->p4.pt());
    double jetPt      = jetPtWithJESJER(vJ[i], jes, jer); 
    double neutralHadEnFrac = jet->neutralHadronEnergyFraction;
    double neutralEmEnFrac = jet->neutralEmEnergyFraction;
    double chargedHadEnFrac = jet->chargedHadronEnergyFraction;

    ///double pujetid    = int(jet->puIDMVALoose);
    ///if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_ && pujetid==1.0  )
    if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_
      && neutralHadEnFrac < 0.99
      && neutralEmEnFrac  < 0.99
      && chargedHadEnFrac > 0
      ){
      j_i->push_back(i);
    }
  }
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
bool ObjectSelector::isMediumMuon(const MyMuon * m, bool isPFlow){
  
  bool isMedium(false);
  bool goodGlob = m->isGlobalMuon && 
	  m->normChi2 <3 && 
	  m->chi2LocalPosition < 12 && 
	  m->trkKink < 20; 
  bool isLooseMuon = m->isPFMuon && 
          (m->isGlobalMuon || m->isTrackerMuon);
  isMedium =  isLooseMuon &&  
	    m->validFraction > 0.8 && 
	    m->segmentCompatibility >(goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}

bool ObjectSelector::isMediumMuonGH(const MyMuon * m, bool isPFlow){
  
  bool isMedium(false);
  bool goodGlob = m->isGlobalMuon && 
	  m->normChi2 <3 && 
	  m->chi2LocalPosition < 12 && 
	  m->trkKink < 20; 
  bool isLooseMuon = m->isPFMuon && 
          (m->isGlobalMuon || m->isTrackerMuon);
  isMedium =  isLooseMuon &&  
	    m->validFraction > 0.49 && 
	    m->segmentCompatibility >(goodGlob ? 0.303 : 0.451); 
  return isMedium; 
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
    if(isMediumMuon(m, isPFlow) && mEta<LOOSE_M_ETA_MAX_  && mPt> LOOSE_M_PT_MIN_ && mRelIso < LOOSE_M_RELISO_MAX_ ){ looseVeto = true; }
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


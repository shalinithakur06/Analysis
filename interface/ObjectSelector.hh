#ifndef _objectselector_h_
#define _objectselector_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

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
#include "TTimeStamp.h"
#include <exception>

#ifdef _STANDALONE
#include "Reader.h"
#else
#include "interface/Reader.h"
#endif
#include "interface/UncertaintyComputer.hh"
#endif

class ObjectSelector : public UncertaintyComputer
{
public : 
  ObjectSelector() : UncertaintyComputer()
  {
    setDefaultSelection();
  }
  virtual ~ObjectSelector(){}
  ///~ObjectSelector(){}

  void setDefaultSelection(){
    defaultSelection_=true;

    //Vertex //////
    ZMAX_    = 0.5; // initally was using 0.2
    ///////////////

    // Jets ///////////////////////////////////////////
    JET_PT_MIN_        = 25;  
    JET_ETA_MAX_       = 2.4; 
    JET_EMF_MIN_       = 0.01;
    JET_LEPTON_DRMIN_  = 0.5;
    //////////////////////////////////////////////////

    // electron //////////////////////////////////////
    E_ETA_MAX_          = 2.1;
    E_PT_MIN_           = 25;
    E_D0_MAX_           = 0.05;
    //////////////////////////////////////////////////


    // muon //////////////////////////////////////////
    M_RELISO_MAX_  = 0.15; //sync with muon POG
    //https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2  
    M_PT_MIN_      = 25;  
    M_ETA_MAX_     = 2.1; 
    M_D0_MAX_      = 0.2; 
    
    LOOSE_M_RELISO_MAX_ = 0.25;
    LOOSE_M_ETA_MAX_    = 2.4;
    LOOSE_M_PT_MIN_     = 15;
    ///////////////////////////////////////////////////
  }

  // preselection of objects
  void preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow=false);
  void preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isData=false, double random_u1=0.5, double random_u2 =0.4, int err_set=0, int err_member=0);
  void preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes=0, int jer=0);

  //Loose Lepton veto
  bool looseElectronVeto(unsigned long selectedElectron, const vector<MyElectron> & vE, MyVertex & vertex, bool isPFlow=false);
  bool looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM, bool isPFlow=false);
  
  //Medium muon ID
  bool isMediumMuon(const MyMuon * m, bool isPFlow);
  bool isMediumMuonGH(const MyMuon * m, bool isPFlow);
  
  //egmGsfElectronIDs: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  bool cutBasedElectronID_Summer16_80X_V1_veto(const MyElectron *e); 
  bool cutBasedElectronID_Summer16_80X_V1_loose(const MyElectron *e); 
  bool cutBasedElectronID_Summer16_80X_V1_medium(const MyElectron *e); 
  bool cutBasedElectronID_Summer16_80X_V1_tight(const MyElectron *e); 

  // object cleaning
  void ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR );

  void JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE, vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, double DR);
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV); 
    
  
private :
  // Vertex
  double ZMAX_;

  // jet
  double JET_PT_MIN_,JET_ETA_MAX_,JET_LTK_PT_MIN_, JET_BTAGGING_, JET_EMF_MIN_, JET_LTK_,JET_LEPTON_DRMIN_  ;
  
  // electron
  double E_RELISO_MAX_, E_ETA_MAX_, E_D0_MAX_, E_PT_MIN_, E_ET_MIN_,RHO_AEFF_E_;
  double LOOSE_E_RELISO_MAX_, LOOSE_E_ETA_MAX_, LOOSE_E_D0_MAX_, LOOSE_E_ET_MIN_;
  
  
  // muon
  double M_RELISO_MAX_, M_ETA_MAX_, M_D0_MAX_, M_PT_MIN_,RHO_AEFF_M_;
  double LOOSE_M_RELISO_MAX_, LOOSE_M_ETA_MAX_, LOOSE_M_D0_MAX_, LOOSE_M_PT_MIN_;
  
  int HPS_ISO_;
  bool defaultSelection_;

  ClassDef(ObjectSelector, 1)
};
#endif

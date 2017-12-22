#ifndef _uncertaintycomputer_h_
#define _uncertaintycomputer_h_

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
#include "interface/SVEffUnc.hh"

#endif


//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
const double JEREtaMap[14] = {0., 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0}; 
const double JERSF[13] = {1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16};
const double JERSFUp[13] = {1.109+0.008 , 1.138+0.013, 1.114+0.013, 1.123+0.024, 1.084+0.011, 1.082+0.035, 1.140+0.047, 1.067+0.053, 1.177+0.041, 1.364+0.039, 1.857+0.071, 1.328+0.022, 1.16+0.029};
const double JERSFDown[13] = {1.109-0.008 , 1.138-0.013, 1.114-0.013, 1.123-0.024, 1.084-0.011, 1.082-0.035, 1.140-0.047, 1.067-0.053, 1.177-0.041, 1.364-0.039, 1.857-0.071, 1.328-0.022, 1.16-0.029};

class UncertaintyComputer{

public :
  UncertaintyComputer()
  {
    sveffunc = new SVEffUnc();
  }

   virtual ~UncertaintyComputer(){
   ///~UncertaintyComputer(){
     delete sveffunc;
  }
  double muPtWithRochCorr(const MyMuon *mu, bool isData=false, double u1=0.5, double u2=0.4, int s=0, int m=0); 
  double metWithJES(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0);
  double metWithJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jer=0);
  double metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0, int jer=0);
  double metWithUncl(const vector<MyJet> & vJ, vector<int> *j, const vector<MyMuon> &vMu, vector<int> *m, const vector<MyElectron> &vEle, vector<int> *el, MyMET MET, int unc=0);
  double getJERSF(double eta, int jer=0);
  double jetPtWithJESJER(MyJet jet, int jes=0, int jer=0); 
  void  openCSVfile(const std::string &filename); 
  double EffUncOnSV(MyJet jet);
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV);
  
private :
  enum BVariation{kNo = 0, kDown = 1, kUp = 2};
  SVEffUnc* sveffunc;

  ClassDef(UncertaintyComputer, 1)
};
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 


double totLumi = 35.9;
TString inFileDir="stack_2016Data_20181104_Mu_sys";
//TString histname="mlZ_max_sig1500";
bool isMuChannel = true;
bool isEleChannel = false;
TFile* fData    = TFile::Open(inFileDir+"/all_muData.root");
//bkg
TFile* fVV      = TFile::Open(inFileDir+"/all_VV.root");
TFile* fDY      = TFile::Open(inFileDir+"/all_DY.root");
TFile* fWJ      = TFile::Open(inFileDir+"/all_WJets.root");
TFile* fTT      = TFile::Open(inFileDir+"/all_TT.root");

//signal
TFile *fLstar250       = TFile::Open(inFileDir+"/all_ExLepMuMuZ_M250.root");
TFile *fLstar1500      = TFile::Open(inFileDir+"/all_ExLepMuMuZ_M1500.root");
TFile *fLstar2000      = TFile::Open(inFileDir+"/all_ExLepMuMuZ_M2000.root");
TFile *fLstar2500      = TFile::Open(inFileDir+"/all_ExLepMuMuZ_M2500.root");
TFile *fLstar4000      = TFile::Open(inFileDir+"/all_ExLepMuMuZ_M4000.root");

//----------------------------------------//
//Variuos functions
//----------------------------------------//
//Read histos from input file. Return empty hist if the hist does not exist. 
TH1F* readHisto(TFile *inFile, TString histPath, TString inHistName, double sf, TString process){
  TH1F* hist;
  if(!(inFile->Get(histPath+inHistName))){
    hist = (TH1F*)(fTT->Get("base/Iso/ControlP/mll"));//initialise an empty hist
    hist->Reset();
  }else hist = (TH1F*)(inFile->Get(histPath+inHistName));
  hist->Scale(sf);
  cout<<setw(10)<<process<<setw(15)<<histPath<<setw(15)<<inHistName<<setw(15)<<hist->Integral()<<endl;
  return hist;
}  

//get normalised uncertainity
double getSysUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown){
  return 1 + max(fabs(hUp->Integral() - hCentral->Integral()), fabs(hCentral->Integral() - hDown->Integral()))/hCentral->Integral();
}

//get statistical uncertainity
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
  return statUnc;
}
//----------------------------------------//
//function to make data card for each mass
//----------------------------------------//
void makeOneDataCard(TFile *fLstar, int mass=250, TString label="lstar250", TString histname = "mlZ_max_sig1500"){
  cout<<" ======> mass point: "<<mass<<endl;
  //ttbar
  double sf_ttbar = 1; 
  TH1F* hTTbar = readHisto(fTT, "base/Iso/ZTag/", histname, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JESUp = readHisto(fTT, "JESPlus/Iso/ZTag/", histname, sf_ttbar, "TTbar");
  TH1F* hTTbar_JERUp = readHisto(fTT, "JERPlus/Iso/ZTag/", histname, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JESDown = readHisto(fTT, "JESMinus/Iso/ZTag/", histname, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JERDown = readHisto(fTT, "JERMinus/Iso/ZTag/", histname, sf_ttbar, "TTbar"); 
  
  //w+jets
  double sf_wjet = 1;
  TH1F* hWJet = readHisto(fWJ, "base/Iso/ZTag/", histname, sf_wjet, "WJets"); 
  TH1F* hWJet_JESUp = readHisto(fWJ, "JESPlus/Iso/ZTag/", histname, sf_wjet, "WJets"); 
  TH1F* hWJet_JERUp = readHisto(fWJ, "JERPlus/Iso/ZTag/", histname, sf_wjet, "WJets"); 
  TH1F* hWJet_JESDown = readHisto(fWJ, "JESMinus/Iso/ZTag/", histname, sf_wjet, "WJets"); 
  TH1F* hWJet_JERDown = readHisto(fWJ, "JERMinus/Iso/ZTag/", histname, sf_wjet, "WJets"); 
  
  //Z+Jets
  double sf_dyjet = 1;
  TH1F* hDYJet = readHisto(fDY, "base/Iso/ZTag/", histname, sf_wjet, "DYJets"); 
  TH1F* hDYJet_JESUp = readHisto(fDY, "JESPlus/Iso/ZTag/", histname, sf_wjet, "DYJets"); 
  TH1F* hDYJet_JERUp = readHisto(fDY, "JERPlus/Iso/ZTag/", histname, sf_wjet, "DYJets"); 
  TH1F* hDYJet_JESDown = readHisto(fDY, "JESMinus/Iso/ZTag/", histname, sf_wjet, "DYJets"); 
  TH1F* hDYJet_JERDown = readHisto(fDY, "JERMinus/Iso/ZTag/", histname, sf_wjet, "DYJets"); 
  
  //Dibosons
  double sf_vv = 1;
  TH1F* hVV = readHisto(fVV, "base/Iso/ZTag/", histname, sf_vv, "VV"); 
  TH1F* hVV_JESUp = readHisto(fVV, "JESPlus/Iso/ZTag/", histname, sf_vv, "VV"); 
  TH1F* hVV_JERUp = readHisto(fVV, "JERPlus/Iso/ZTag/", histname, sf_vv, "VV"); 
  TH1F* hVV_JESDown = readHisto(fVV, "JESMinus/Iso/ZTag/", histname, sf_vv, "VV"); 
  TH1F* hVV_JERDown = readHisto(fVV, "JERMinus/Iso/ZTag/", histname, sf_vv, "VV"); 
  //Data
  double sf_data = 1; //should be 1, always
  TH1F* hData = readHisto(fData, "base/Iso/ZTag/", histname, sf_data, "Data"); 
  //lstar
  double sf_lstar = 1; 
  TH1F* hLstar = readHisto(fLstar, "base/Iso/ZTag/", histname, sf_lstar, "Signal"); 
  TH1F* hLstar_JESUp = readHisto(fLstar, "JESPlus/Iso/ZTag/", histname, sf_lstar, "Signal"); 
  TH1F* hLstar_JERUp = readHisto(fLstar, "JERPlus/Iso/ZTag/", histname, sf_lstar, "Signal"); 
  TH1F* hLstar_JESDown = readHisto(fLstar, "JESMinus/Iso/ZTag/", histname, sf_lstar, "Signal"); 
  TH1F* hLstar_JERDown = readHisto(fLstar, "JERMinus/Iso/ZTag/", histname, sf_lstar, "Signal"); 
  //open input template data card of 8 TeV
  ifstream in;
  char* c = new char[1000];
  in.open("template_datacard_count.txt");
  //create output data card for 13 TeV
  string outDataCard = "datacard_llstar_llZ_llq_13TeV_M250.txt";
  string histname_str(histname);
  if(isMuChannel) outDataCard = "datacard_llstar_llZ_llq_13TeV_mu_M%d.txt"; 
  else outDataCard = "datacard_llstar_llZ_llq_13TeV_ele_M%d.txt";
  ofstream out(Form(outDataCard.c_str(), mass));
  out.precision(8);

  time_t secs=time(0);
  tm *t=localtime(&secs);
  while (in.good()){
    in.getline(c,1000,'\n');
    if (in.good()){
      string line(c);
      if(line.find("Date")!=string::npos){
        string day = string(Form("%d",t->tm_mday));
        string month = string(Form("%d",t->tm_mon+1));
        string year = string(Form("%d",t->tm_year+1900));
        line.replace( line.find("XXX") , 3 , day+"/"+month+"/"+year);
        out << line << endl;
      }
      else if(line.find("Description")!=string::npos){
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        line.replace( line.find("ZZZ") , 3 , string(Form("%f", totLumi)) ); 
        out << line << endl;
      }
      else if(line.find("Observation")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%.0f", hData->Integral())));
        out << line << endl;
      }
      else if(line.find("process")!=string::npos && line.find("lstar")!=string::npos){
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        out << line << endl;
      }
      else if(line.find("rate")!=string::npos){
        string rate = "rate               ";  
        string space = "     ";
        out << rate ;
	out << space << hLstar->Integral()
            << space << hTTbar->Integral()
            << space << hWJet->Integral()
            << space << hDYJet->Integral()
            << space << hVV->Integral()
            << endl;
      }
      else if(line.find("CMS_stat_lstar")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(hLstar,  0))));  
        out << line << endl;
      } 
      else if(line.find("CMS_stat_tt")!=string::npos){  
	line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(hTTbar,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_wjet")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(hWJet,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_zjet")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(hDYJet,  0))));  
        out << line << endl; 
      }
      else if(line.find("CMS_stat_vv")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(hVV,  0))));   
        out << line << endl;  
      }
     else if(line.find("CMS_scale_j")!=string::npos){
        float JESUnc_lstar = (hLstar->Integral() > 0) ? getSysUnc(hLstar, hLstar_JESUp, hLstar_JESDown) : 1.00;
        line.replace( line.find("LLLL") , 4 , string(Form("%.3f", JESUnc_lstar)) );

        float JESUnc_ttbar = (hTTbar->Integral() > 0) ? getSysUnc(hTTbar, hTTbar_JESUp, hTTbar_JESDown) : 1.00;
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", JESUnc_ttbar)) );

        float JESUnc_wjet = (hWJet->Integral() > 0) ? getSysUnc(hWJet, hWJet_JESUp, hWJet_JESDown) : 1.00;
        line.replace( line.find("WWWW") , 4 , string(Form("%.3f", JESUnc_wjet)) );

        float JESUnc_zjet = (hDYJet->Integral() > 0) ? getSysUnc(hDYJet, hDYJet_JESUp, hDYJet_JESDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", JESUnc_zjet)) );

        float JESUnc_vv = (hVV->Integral() > 0) ? getSysUnc(hVV, hVV_JESUp, hVV_JESDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", JESUnc_vv)) );
        out << line << endl;
      }
     else if(line.find("CMS_res_j")!=string::npos){
        float JERUnc_hLstar = (hLstar->Integral() > 0) ? getSysUnc(hLstar, hLstar_JERUp, hLstar_JERDown) : 1.00;
        line.replace( line.find("LLLL") , 4 , string(Form("%.3f", JERUnc_hLstar)) );

        float JERUnc_ttbar = (hTTbar->Integral() > 0) ? getSysUnc(hTTbar, hTTbar_JERUp, hTTbar_JERDown) : 1.00;
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", JERUnc_ttbar)) );

        float JERUnc_wjet = (hWJet->Integral() > 0) ? getSysUnc(hWJet, hWJet_JERUp, hWJet_JERDown) : 1.00;
        line.replace( line.find("WWWW") , 4 , string(Form("%.3f", JERUnc_wjet)) );

        float JERUnc_zjet = (hDYJet->Integral() > 0) ? getSysUnc(hDYJet, hDYJet_JERUp, hDYJet_JERDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", JERUnc_zjet)) );

        float JERUnc_vv = (hVV->Integral() > 0) ? getSysUnc(hVV, hVV_JERUp, hVV_JERDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", JERUnc_vv)) );
        out << line << endl;
      }
      else{ //default without changes
        out << line << endl;
      }
    }
  } 
  out.close();
  in.close();
}

//----------------------------------------//
//make data card for all masses
//----------------------------------------//
void makeExLepDataCard_13TeV_count(){
  makeOneDataCard(fLstar250, 250,  "lstar250", "mlZ_max_sig250");
  makeOneDataCard(fLstar1500, 1500,  "lstar1500", "mlZ_max_sig1500");
  makeOneDataCard(fLstar2000, 2000,  "lstar2000", "mlZ_max_sig2000");
  makeOneDataCard(fLstar2500, 2500,  "lstar2500", "mlZ_max_sig2500");
  makeOneDataCard(fLstar4000, 4000,  "lstar4000", "mlZ_max_sig4000");
}

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 

/////////////////////////// USERS INPUT ///////////////////////////
///INPUT FILES
//data
double totLumi = 35.9;
TString inFileDir="stack_20180224_Mu_sys_2";
//TString inFileDir="stack_20171223_Mu_sys";
TString histname="mjj_max";
bool isMuChannel = true;
bool isEleChannel = false;

TFile* fData    = TFile::Open(inFileDir+"/all_muData.root");
//OUTPUT FILE
TString outShapeFile ="ShapesExcitedLepton";
TFile *fout = new TFile(outShapeFile+"_mu_"+histname+"_13TeV.root", "RECREATE");

//bkg
TFile* fVV      = TFile::Open(inFileDir+"/all_VV.root");
TFile* fDY      = TFile::Open(inFileDir+"/all_DY.root");
TFile* fWJ      = TFile::Open(inFileDir+"/all_WJets.root");
TFile* fTT      = TFile::Open(inFileDir+"/all_TTJetsP.root");
TFile* fQCD      = TFile::Open(inFileDir+"/all_QCD.root");

//signal
TFile *fWH250      = TFile::Open(inFileDir+"/all_signal250.root");
TFile *fWH1500      = TFile::Open(inFileDir+"/all_signal1500.root");
TFile *fWH2000      = TFile::Open(inFileDir+"/all_signal2000.root");
TFile *fWH4000      = TFile::Open(inFileDir+"/all_signal4000.root");


//////////////////////////////////////////////////////////////////
//----------------------------------------//
//Variuos functions
//----------------------------------------//
//Read histos from input file. Return empty hist if the hist does not exist. Write to another file.
TH1F* readWriteHisto(TFile *inFile, TString histPath, TString inHistName, double sf, TFile *outFile, TString outHistName, bool isWrite = false){
  TH1F* hist;
  if(!(inFile->Get(histPath+inHistName))){
    hist = (TH1F*)(fTT->Get(histPath+inHistName))->Clone(outHistName);
    hist->Add(hist, -1);
  }else hist = (TH1F*)(inFile->Get(histPath+inHistName))->Clone(outHistName);
  hist->Scale(sf);
  cout<<inHistName<<"\t"<<outHistName<<"\t"<<"\t"<<hist->Integral()<<endl;
  ///hist->GetXaxis()->SetRangeUser(0, 200);
  if(isWrite){
    fout->cd();
    hist->Write(outHistName);
  }
  return hist;
}  

//get normalised uncertainity
double getBTagUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown){
  return 1 + max(fabs(hUp->Integral() - hCentral->Integral()), fabs(hCentral->Integral() - hDown->Integral()))/hCentral->Integral();
}

//get statistical uncertainity
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
  return statUnc;
}
//----------------------------------------//
//Global variables
//----------------------------------------//
//ttbar
double sf_ttbar = 1; 
TH1F* ttbar = readWriteHisto(fTT, "base/Iso/ZTag/", histname, sf_ttbar, fout, "ttbar", true);
TH1F* ttbar_JESUp = readWriteHisto(fTT, "JESPlus/Iso/ZTag/", histname, sf_ttbar, fout, "ttbar_JESUp", true);
TH1F* ttbar_JERUp = readWriteHisto(fTT, "JERPlus/Iso/ZTag/", histname, sf_ttbar, fout, "ttbar_JERUp", true);
TH1F* ttbar_JESDown = readWriteHisto(fTT, "JESMinus/Iso/ZTag/", histname, sf_ttbar, fout, "ttbar_JESDown", true);
TH1F* ttbar_JERDown = readWriteHisto(fTT, "JERMinus/Iso/ZTag/", histname, sf_ttbar, fout, "ttbar_JERDown", true);

//ttll
double sf_ttll = 0;
TH1F* ttll = readWriteHisto(fTT, "base/Iso/ZTag/", histname, sf_ttll, fout, "ttll", true);
TH1F* ttll_JESUp = readWriteHisto(fTT, "JESPlus/Iso/ZTag/", histname, sf_ttll, fout, "ttll_JESUp", true);
TH1F* ttll_JESDown = readWriteHisto(fTT, "JESMinus/Iso/ZTag/", histname, sf_ttll, fout, "ttll_JESDown", true);

//w+jets
double sf_wjet = 1;
TH1F* wjet = readWriteHisto(fWJ, "base/Iso/ZTag/", histname, sf_wjet, fout, "wjet", true);
TH1F* wjet_JESUp = readWriteHisto(fWJ, "JESPlus/Iso/ZTag/", histname, sf_wjet, fout, "wjet_JESUp", true);
TH1F* wjet_JERUp = readWriteHisto(fWJ, "JERPlus/Iso/ZTag/", histname, sf_wjet, fout, "wjet_JERUp", true);
TH1F* wjet_JESDown = readWriteHisto(fWJ, "JESMinus/Iso/ZTag/", histname, sf_wjet, fout, "wjet_JESDown", true);
TH1F* wjet_JERDown = readWriteHisto(fWJ, "JERMinus/Iso/ZTag/", histname, sf_wjet, fout, "wjet_JERDown", true);

//Z+Jets
double sf_zjet = 1;
TH1F* zjet = readWriteHisto(fDY, "base/Iso/ZTag/", histname, sf_zjet, fout, "zjet", true);
TH1F* zjet_JESUp = readWriteHisto(fDY, "JESPlus/Iso/ZTag/", histname, sf_zjet, fout, "zjet_JESUp", true);
TH1F* zjet_JERUp = readWriteHisto(fDY, "JERPlus/Iso/ZTag/", histname, sf_zjet, fout, "zjet_JERUp", true);
TH1F* zjet_JESDown = readWriteHisto(fDY, "JESMinus/Iso/ZTag/", histname, sf_zjet, fout, "zjet_JESDown", true);
TH1F* zjet_JERDown = readWriteHisto(fDY, "JERMinus/Iso/ZTag/", histname, sf_zjet, fout, "zjet_JERDown", true);

//Dibosons
double sf_diboson = 1;
TH1F* diboson = readWriteHisto(fVV, "base/Iso/ZTag/", histname, sf_diboson, fout, "diboson", true);
TH1F* diboson_JESUp = readWriteHisto(fVV, "JESPlus/Iso/ZTag/", histname, sf_diboson, fout, "diboson_JESUp", true);
TH1F* diboson_JERUp = readWriteHisto(fVV, "JERPlus/Iso/ZTag/", histname, sf_diboson, fout, "diboson_JERUp", true);
TH1F* diboson_JESDown = readWriteHisto(fVV, "JESMinus/Iso/ZTag/", histname, sf_diboson, fout, "diboson_JESDown", true);
TH1F* diboson_JERDown = readWriteHisto(fVV, "JERMinus/Iso/ZTag/", histname, sf_diboson, fout, "diboson_JERDown", true);

//QCD MC
double sf_qcd = 1;
TH1F* qcd = readWriteHisto(fQCD, "base/Iso/ZTag/", histname, sf_qcd, fout, "qcd", true);
TH1F* qcd_JESUp = readWriteHisto(fQCD, "JESPlus/Iso/ZTag/", histname, sf_qcd, fout, "qcd_JESUp", true);
TH1F* qcd_JERUp = readWriteHisto(fQCD, "JERPlus/Iso/ZTag/", histname, sf_qcd, fout, "qcd_JERUp", true);
TH1F* qcd_JESDown = readWriteHisto(fQCD, "JESMinus/Iso/ZTag/", histname, sf_qcd, fout, "qcd_JESDown", true);
TH1F* qcd_JERDown = readWriteHisto(fQCD, "JERMinus/Iso/ZTag/", histname, sf_qcd, fout, "qcd_JERDown", true);

//Data
double sf_data = 1; //should be 1, always
TH1F* data_obs = readWriteHisto(fData, "base/Iso/ZTag/", histname, sf_data, fout, "data_obs", true);

//----------------------------------------//
//function to make data card for each mass
//----------------------------------------//
void makeOneDataCard(TFile *fWH, int mass=250, TString label="WH250", TString label2="HH250", TString histname = "mjj"){
  cout<<" mass point: "<<mass<<endl;
  //wh
  double sf_wh = 1; 
  TH1F* wh = readWriteHisto(fWH, "base/Iso/ZTag/", histname, sf_wh, fout, label, true);
  TH1F* wh_JESUp = readWriteHisto(fWH, "JESPlus/Iso/ZTag/", histname, sf_wh, fout, label+"_JESUp", true);
  TH1F* wh_JERUp = readWriteHisto(fWH, "JERPlus/Iso/ZTag/", histname, sf_wh, fout, label+"_JERUp", true);
  TH1F* wh_JESDown = readWriteHisto(fWH, "JESMinus/Iso/ZTag/", histname, sf_wh, fout, label+"_JESDown", true);
  TH1F* wh_JERDown = readWriteHisto(fWH, "JERMinus/Iso/ZTag/", histname, sf_wh, fout, label+"_JERDown", true);
  
  //open input template data card of 8 TeV
  ifstream in;
  char* c = new char[1000];
  in.open("template_datacard_shape.txt");
  //create output data card for 13 TeV
  string outDataCard = "datacard_exlep_13TeV_mH.txt";
  string histname_str(histname);
  if(isMuChannel) outDataCard = "datacard_exlep_mu_"+histname_str+"_13TeV_mH%d.txt"; 
  else outDataCard = "datacard_exlep_ele_"+histname_str+"_13TeV_mH%d.txt";
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
      else if(line.find("shapes")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(outShapeFile+"_mu_"+histname+"_13TeV"));
        out << line << endl;
      }
      else if(line.find("Observation")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%.0f", data_obs->Integral())));
        out << line << endl;
      }
      else if(line.find("process")!=string::npos && line.find("WH")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%d", mass)) );
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        out << line << endl;
      }
      else if(line.find("rate")!=string::npos){
        string rate = "rate               ";  
        string space = "     ";
        out << rate ;
        out << space << 0*wh->Integral()
	    << space << wh->Integral()
            << space << ttbar->Integral()
            << space << wjet->Integral()
            << space << zjet->Integral()
            << space << diboson->Integral()
            << space << qcd->Integral()
            << endl;
      }
      else if(line.find("CMS_stat_wh")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(wh,  0))));  
        out << line << endl;
      } 
      else if(line.find("CMS_stat_tt")!=string::npos){  
	line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(ttbar,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_wjet")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(wjet,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_zjet")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(zjet,  0))));  
        out << line << endl; 
      }
      else if(line.find("CMS_stat_vv")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(diboson,  0))));   
        out << line << endl;  
      }
      else if(line.find("CMS_stat_qcd")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", getStatUnc(qcd,  0))));   
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
void makeExLepDataCard_13TeV(){
  makeOneDataCard(fWH250,  250,   "WH250",   "HH250",  histname);
  makeOneDataCard(fWH1500, 1500,  "WH1500", "HH1500",  histname);
  makeOneDataCard(fWH2000, 2000,  "WH2000", "HH2000",  histname);
  makeOneDataCard(fWH4000, 4000,  "WH4000", "HH4000",  histname);
  fout->Close();
}

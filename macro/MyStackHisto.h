#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

using namespace std;

//R. K. Verma
//Sat Jul 14 14:47:08 IST 2018

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY_M50.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile *fSig     = TFile::Open("all_ExLep_M250.root");

TString baseDir = "baseLowMET";
class MyStackHisto{
  public : 
	//get histogram from root file. Return empty hist, if the hist does not exist.
	TH1F* getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName);
        //decorate histogram
        TH1F* decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit);
	//function to stack histos
        void stackHisto(TFile *inRootFile, TString lable, TString baseDir, TString isoDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, TLegend* leg);
	//function to add histograms
	TH1F*  addHistoForUnc(TString baseDir, TString isoDir, TString histDir, TString histName);
	//Up/down error in the unc band
	double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus);
	double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus);
	//unc graph
	TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, bool isFullGraph = false, bool isRatioGraph = false);
        TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );
  private :
	int dont_use ;
};

//--------------------------------------------//
//define various functions
//--------------------------------------------//
TH1F*  MyStackHisto:: getHisto(TFile *inRootFile, TString baseDir, TString isoDir, TString histDir, TString histName){
  TH1F* hist;
  TString fullPath = baseDir+isoDir+histDir+histName;
  //TString fullPath = string(inRootFile->GetName())+"/"+baseDir+isoDir+histDir+histName;
  string exception_msg ("The histogram path, "+fullPath+", does not exist"); 
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"WARNING:"<<e<<endl;
  }
  try{
    if(!(fTT->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"\033[01;31mERROR: \033[0m"<<e<< endl;
    exit(0);
  }
  if(!(inRootFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(inRootFile->Get(fullPath))->Clone(histName);
  return hist;
}

TH1F* MyStackHisto:: decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit){
  hist->SetTitle(myTit);
  hist->GetXaxis()->SetTitle(xTit);
  hist->GetYaxis()->SetTitle(yTit);
  hist->GetYaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetTitleOffset(1.00);
  hist->GetYaxis()->SetTitleSize(0.10);   
  hist->GetXaxis()->SetTitleSize(0.10);
  hist->GetXaxis()->SetLabelSize(0.10);   
  hist->GetYaxis()->SetLabelSize(0.10);   
  hist->GetXaxis()->SetTickLength(0.05);
  hist->GetXaxis()->SetNdivisions(10);
  hist->GetYaxis()->SetNdivisions(5);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTickLength(0.04);
  return hist;
}

void  MyStackHisto:: stackHisto(TFile *inRootFile, TString lable, TString baseDir, TString isoDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, TLegend* leg){
  TH1F* hist = getHisto(inRootFile, baseDir, isoDir, histDir, histName);
  hist->Scale(scale);  
  hist->SetFillColor(color);
  leg->AddEntry(hist,lable,"F");
  MuptStack->Add(hist);
  hMC->Add(hist);
}

double MyStackHisto:: errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus){
  double errUp = sqrt(pow(fabs(hJESPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hJERPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bTagPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2));
  return errUp;
}

double MyStackHisto:: errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus){
  double errDown =sqrt(pow(fabs(hCentral->GetBinContent(iBin+1) - hJESMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - hJERMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bTagMinus->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2));
  return errDown;
}

TGraphAsymmErrors * MyStackHisto:: UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, bool isFullGraph = false, bool isRatioGraph = false){
  TGraphAsymmErrors *gr;
  int n1 = hCentral->GetNbinsX(); 
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  //cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++){
    if(isFullGraph){
    Yval[i]   = hCentral->GetBinContent(i+1);
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus); 
    }
    if(isRatioGraph){
    Yval[i]   = 1;
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus); 
    //cout<<"bin = "<<i<<endl;
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
    errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
    errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
    }
    Xval[i]   = hCentral->GetBinCenter(i+1);
    XerrorU[i]= hCentral->GetBinWidth(i+1)/2;
    XerrorD[i]= hCentral->GetBinWidth(i+1)/2;
  }
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
}

TH1F*  MyStackHisto:: addHistoForUnc(TString baseDir, TString isoDir, TString histDir, TString histName){
  TH1F* hVV =   	getHisto(fVV,   baseDir, isoDir, histDir, histName);
  TH1F* hDY =   	getHisto(fDY,    baseDir, isoDir, histDir, histName);
  TH1F* hWJ =   	getHisto(fWJ,   baseDir, isoDir, histDir, histName);
  TH1F* hTT =   	getHisto(fTT,   baseDir, isoDir, histDir, histName);
  TH1F* hAll = (TH1F*)hVV->Clone("hAllMC");
  hAll->Add(hDY);
  hAll->Add(hWJ);
  hAll->Add(hTT);
  return hAll;
}

TPaveText * MyStackHisto:: paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}

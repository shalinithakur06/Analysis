#include <iostream>
#include <fstream>
#include <iomanip>

//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY_M50.root");
//TFile* fDY	= TFile::Open("all_DY_Pt.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fBkg	= TFile::Open("all_Bkg.root");

TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TString fullPath = histPath+"/"+dir+"/"+histName;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get("base/Iso/ZTag/mll")->Clone(histName));
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

//get statistical uncertainity
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = fabs(sError);
  return statUnc;
}

string doubleToStr(double val){
     ostringstream convert;
     string result("");
     convert <<std::setprecision(3)<<val;
     result = convert.str();
  return result;
}

double getSysError( TH1F * h_JESPlus, TH1F * h_base, TH1F * h_JESMinus, TH1F * h_JERPlus, TH1F * h_JERMinus, TH1F * h_bTagPlus, TH1F * h_bTagMinus){
  double uncJES = pow(TMath::Max(fabs(h_JESPlus->Integral() - h_base->Integral()), fabs(h_base->Integral() - h_JESMinus->Integral())), 2);
  double uncJER = pow(TMath::Max(fabs(h_JERPlus->Integral() - h_base->Integral()), fabs(h_base->Integral() - h_JERMinus->Integral())), 2);
  double uncBTag = pow(TMath::Max(fabs(h_bTagPlus->Integral() - h_base->Integral()), fabs(h_base->Integral() - h_bTagMinus->Integral())), 2);
  double unc = sqrt(uncJES +uncJER  +uncBTag);
  return unc ;
}

string getBkgYield(TFile *file, TString histName){
  TH1F * hBase = (TH1F*)getHisto(file, "base/Iso", "ZTag", histName);
  double yield = hBase->Integral();
  //stat err
  double stat_err = getStatUnc(hBase, 0.0);
  //sys error
  TH1F * h_JESPlus  = (TH1F*)getHisto(file, "JESPlus/Iso", "ZTag", histName);
  TH1F * h_JESMinus = (TH1F*)getHisto(file, "JESMinus/Iso", "ZTag", histName);
  TH1F * h_JERPlus  = (TH1F*)getHisto(file, "JERPlus/Iso", "ZTag", histName);
  TH1F * h_JERMinus = (TH1F*)getHisto(file, "JERMinus/Iso", "ZTag", histName);
  TH1F * h_bTagPlus = (TH1F*)getHisto(file, "bTagPlus/Iso", "ZTag", histName);
  TH1F * h_bTagMinus= (TH1F*)getHisto(file, "bTagMinus/Iso", "ZTag", histName);
  double sys_err = 0.0;
  sys_err = getSysError(h_JESPlus, hBase, h_JESMinus, h_JERPlus, h_JERMinus, h_bTagPlus, h_bTagMinus);
  //retrun yield + stat_err + sys_err
  return doubleToStr(yield)+"$\\pm$"+doubleToStr(stat_err)+"$\\pm$"+doubleToStr(sys_err);
}

string getSigYield(TString sigFile, TString histName){
  TFile *file     = TFile::Open(sigFile);
  TH1F * hBase = (TH1F*)getHisto(file, "base/Iso", "ZTag", histName);
  double yield = hBase->Integral();
  //stat err
  double stat_err = getStatUnc(hBase, 0.0);
  //sys error
  TH1F * h_JESPlus  = (TH1F*)getHisto(file, "JESPlus/Iso", "ZTag", histName);
  TH1F * h_JESMinus = (TH1F*)getHisto(file, "JESMinus/Iso", "ZTag", histName);
  TH1F * h_JERPlus  = (TH1F*)getHisto(file, "JERPlus/Iso", "ZTag", histName);
  TH1F * h_JERMinus = (TH1F*)getHisto(file, "JERMinus/Iso", "ZTag", histName);
  TH1F * h_bTagPlus = (TH1F*)getHisto(file, "bTagPlus/Iso", "ZTag", histName);
  TH1F * h_bTagMinus= (TH1F*)getHisto(file, "bTagMinus/Iso", "ZTag", histName);
  double sys_err = 0.0;
  sys_err = getSysError(h_JESPlus, hBase, h_JESMinus, h_JERPlus, h_JERMinus, h_bTagPlus, h_bTagMinus);
  //retrun yield + stat_err + sys_err
  return doubleToStr(yield)+"$\\pm$"+doubleToStr(stat_err)+"$\\pm$"+doubleToStr(sys_err);
}

void MyLCutTable(){  
  ofstream outFile; 
  outFile.open("lCutTable.tex"); 
  outFile<<"\\documentclass[]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\textheight = 8.in"<<endl;  
  outFile<<"\\textwidth 7.0in"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<""<<endl;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outFile<<"\\begin{table}"<<endl;
  outFile<<"\\begin{center}"<<endl;  
  //outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\scriptsize{"<<endl;  
  outFile<<"\\begin{tabular}{ |r|r|r|r|r|r|r|}"<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"L-Cut & "<<"$M_{l^*}$ & "<<"$N_{sig}\\pm stat \\pm sys $ &"<<"$N_{t\\bar{t}}\\pm stat \\pm sys $ & "<<"$N_{DY}\\pm stat \\pm sys $ & "<<"$N_{VV}\\pm stat \\pm sys $ &"<<"$N_{Bkg}\\pm stat \\pm sys$\\\\"<<endl;
  outFile<<"(X, Y) & "<<"(GeV) & "<<"&"<<"& "<<"&"<<"&"<<"\\\\"<<endl;
  outFile<<"\\hline "<<endl; 

  vector<string>lCutVec;
  lCutVec.push_back("$440, 300$");
  lCutVec.push_back("$450, 560$");
  lCutVec.push_back("$700, 900$");
  lCutVec.push_back("$950, 1080$");
  lCutVec.push_back("$1200, 1370$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  lCutVec.push_back("$1300, 1700$");
  
  vector<string>sigMassVec;
  sigMassVec.push_back("250");
  sigMassVec.push_back("500");
  sigMassVec.push_back("750");
  sigMassVec.push_back("1000");
  sigMassVec.push_back("1250");
  sigMassVec.push_back("1500");
  sigMassVec.push_back("1750");
  sigMassVec.push_back("2000");
  sigMassVec.push_back("2500");
  sigMassVec.push_back("3000");
  sigMassVec.push_back("3500");
  sigMassVec.push_back("4000");
  sigMassVec.push_back("4500");
  sigMassVec.push_back("5000");
  
  for(int i =0; i<lCutVec.size(); i++){
    TString sigFile = "all_ExLep_M"+sigMassVec[i]+".root";
    TString hist = "mlZ_max_sig"+sigMassVec[i];
    outFile<<lCutVec[i]<<" & "<<sigMassVec[i]<<" & "<<
	  getSigYield(sigFile, hist)<<" & "<<
	  getBkgYield(fTT, hist)<<" & "<<
	  getBkgYield(fDY, hist)<<" & "<<
	  getBkgYield(fVV, hist)<<" & "<<
	  getBkgYield(fBkg, hist)<<" \\\\"<<endl;
  }
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  string chName = "muon";
  if(isEleChannel) chName = "electron";
  outFile<<"\\caption{Event yields after different L-cut $(M_{lZ}^{max} > X, M_{lZ}^{min} < Y)$ for "+chName+" channel.}"<<endl;
  outFile<<"\\end{table}"<<endl;

  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 

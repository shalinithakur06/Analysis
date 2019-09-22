#include <iostream>
#include <fstream>
#include <iomanip>

//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
//TFile* fDY	= TFile::Open("all_DY_Pt.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
//TFile *fSig     = TFile::Open("all_ExLepEEZ_M250.root");
TFile *fSig250     = TFile::Open("all_ExLep_M250.root");
TFile *fSig500     = TFile::Open("all_ExLep_M500.root");
TFile *fSig750     = TFile::Open("all_ExLep_M750.root");
TFile *fSig1000     = TFile::Open("all_ExLep_M1000.root");
TFile *fSig1250     = TFile::Open("all_ExLep_M1250.root");
TFile *fSig1500     = TFile::Open("all_ExLep_M1500.root");
TFile *fSig1750     = TFile::Open("all_ExLep_M1750.root");
TFile *fSig2000     = TFile::Open("all_ExLep_M2000.root");
TFile *fSig2500     = TFile::Open("all_ExLep_M2500.root");
TFile *fSig3000     = TFile::Open("all_ExLep_M3000.root");
TFile *fSig3500     = TFile::Open("all_ExLep_M3500.root");
TFile *fSig4000     = TFile::Open("all_ExLep_M4000.root");
TFile *fSig4500     = TFile::Open("all_ExLep_M4500.root");
TFile *fSig5000     = TFile::Open("all_ExLep_M5000.root");

string doubleToStr(double val){
     ostringstream convert;
     string result("");
     convert <<std::setprecision(4)<<val;
     result = convert.str();
  return result;
}

string makeCutFlowRow(TFile *file, string rowLable){
  TH1F * hCutFlow = (TH1F*)file->Get("base/Iso/cutflow");
  TH1F * hTotalEvents = (TH1F*)file->Get("base/totalYield");
  double totalEvents = hTotalEvents->Integral();
  cout<<rowLable<<"\t"<<totalEvents<<endl;
  string row = rowLable;
  for(int i=1; i<hCutFlow->GetSize()-1; i++){
    if(i<=5){
    //if(hCutFlow->GetBinContent(i)!=0){
      double preYield = 1.0; //previous yield
      if (i==1) preYield = totalEvents; 
      else preYield = hCutFlow->GetBinContent(i-1);
      double currYield = 1.0; //current yueld
      currYield = hCutFlow->GetBinContent(i);
      double eff = 100*currYield/preYield;
      row += " & "+ doubleToStr(currYield)+" ("+doubleToStr(eff)+")";
    }
  }
  row += "\\\\";
  //cout<<row<<endl;
  return row;
}
string makeCutFlowRowRatio(TFile *fData, string rowLableBkg, string rowLableRatio, bool isRatio=false){
  TH1F * hCutFlowData = (TH1F*)fData->Get("base/Iso/cutflow");
  //add all bkgs
  TH1F * hEventsTT = (TH1F*)fTT->Get("base/totalYield");
  TH1F * hEventsDY = (TH1F*)fDY->Get("base/totalYield");
  TH1F * hEventsWJ = (TH1F*)fWJ->Get("base/totalYield");
  TH1F * hEventsVV = (TH1F*)fVV->Get("base/totalYield");
  TH1F * hEventsBkg = (TH1F*)hEventsTT->Clone();
  hEventsBkg->Add(hEventsDY);
  hEventsBkg->Add(hEventsWJ);
  hEventsBkg->Add(hEventsVV);
  
  TH1F * hCutFlowTT = (TH1F*)fTT->Get("base/Iso/cutflow");
  TH1F * hCutFlowDY = (TH1F*)fDY->Get("base/Iso/cutflow");
  TH1F * hCutFlowWJ = (TH1F*)fWJ->Get("base/Iso/cutflow");
  TH1F * hCutFlowVV = (TH1F*)fVV->Get("base/Iso/cutflow");
  TH1F * hCutFlowBkg = (TH1F*)hCutFlowTT->Clone();
  hCutFlowBkg->Add(hCutFlowDY);
  hCutFlowBkg->Add(hCutFlowWJ);
  hCutFlowBkg->Add(hCutFlowVV);

  //print row for all bkgs
  double totalEvents = hEventsBkg->Integral();
  string rowBkg = rowLableBkg;
  for(int i=1; i<hCutFlowBkg->GetSize()-1; i++){
    if(i<=5){
    //if(hCutFlowBkg->GetBinContent(i)!=0){
      double preYield = 1.0; //previous yield
      if (i==1) preYield = totalEvents; 
      else preYield = hCutFlowBkg->GetBinContent(i-1);
      double currYield = 1.0; //current yueld
      currYield = hCutFlowBkg->GetBinContent(i);
      double eff = 100*currYield/preYield;
      rowBkg += " & "+ doubleToStr(currYield)+" ("+doubleToStr(eff)+")";
    }
  }
  rowBkg += "\\\\";

  //print row for data/MC ratio
  string rowRatio = rowLableRatio;
  for(int i=1; i<hCutFlowBkg->GetSize()-1; i++){
    if(i<=5){
      double yieldBkg = 1.0; //current yueld
      double yieldData = 1.0; //current yueld
      yieldBkg = hCutFlowBkg->GetBinContent(i);
      yieldData = hCutFlowData->GetBinContent(i);
      rowRatio += " & "+ doubleToStr(yieldData/yieldBkg);
    }
  }
  rowRatio += "\\\\";

  //return 
  if(isRatio){
    cout<<rowRatio<<endl;
    return rowRatio;
  }
  else{
    //cout<<rowBkg<<endl;
    return rowBkg;
  }
}

void MyCutFlowTable(){  
  ofstream outFile; 
  outFile.open("cutflowTable.tex"); 
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
  outFile<<"\\begin{tabular}{ |c|c|c|c|c|c|}"<<endl;  
  outFile<<"\\hline "<<endl; 
  if(isMuChannel) outFile<<" & "<<"Muon trigger & "<<"ControlSel & "<<"b-jet Veto & "<<"PreSel & "<<"ZTag \\\\"<<endl;
  if(isEleChannel) outFile<<" & "<<"Electron trigger & "<<"ControlSel & "<<"b-jet Veto & "<<"PreSel & "<<"ZTag \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"Process & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<makeCutFlowRow(fSig250 , "$M_{l^*} $ = 250  GeV")<<endl;
  outFile<<makeCutFlowRow(fSig500 , "$M_{l^*} $ = 500  GeV")<<endl;
  outFile<<makeCutFlowRow(fSig750 , "$M_{l^*} $ = 750  GeV")<<endl;
  outFile<<makeCutFlowRow(fSig1000, "$M_{l^*} $ = 1000 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig1250, "$M_{l^*} $ = 1250 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig1500, "$M_{l^*} $ = 1500 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig1750, "$M_{l^*} $ = 1750 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig2000, "$M_{l^*} $ = 2000 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig2500, "$M_{l^*} $ = 2500 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig3000, "$M_{l^*} $ = 3000 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig3500, "$M_{l^*} $ = 3500 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig4000, "$M_{l^*} $ = 4000 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig4500, "$M_{l^*} $ = 4500 GeV")<<endl;
  outFile<<makeCutFlowRow(fSig5000, "$M_{l^*} $ = 5000 GeV")<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<makeCutFlowRow(fTT, "$t\\bar{t}$ + jets")<<endl;
  outFile<<makeCutFlowRow(fDY, "$Z/\\gamma$+ jets")<<endl;
  outFile<<makeCutFlowRow(fVV, "$VV$")<<endl;
  outFile<<makeCutFlowRow(fWJ, "$W$ + jets")<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<makeCutFlowRowRatio(fData, "All Bkg", "Data/MC", false)<<endl;
  outFile<<makeCutFlowRow(fData, "Data")<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<makeCutFlowRowRatio(fData, "All Bkg", "Data/MC", true)<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  string chName = "muon";
  if(isEleChannel) chName = "electron";
  outFile<<"\\caption{Cutflow table for "+chName+" channel. The $\\epsilon$ is the efficiency $(N_{after}/N_{before})$ and $N_{evt} (=N_{after})$ is the number of survived events after different cuts.}"<<endl;
  outFile<<"\\end{table}"<<endl;

  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 

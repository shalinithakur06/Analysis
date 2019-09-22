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

using namespace std;
 
///////////////////////////////////////////  
//CHANNEL
bool isMuChannel = false;
bool isEleChannel = true;

//INPUT FILES
TFile* fData = TFile::Open("all_EleData.root");
//if(isMuChannel)fData  = TFile::Open("all_muData.root");

TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fST	= TFile::Open("all_ST.root");
TFile* fTT	= TFile::Open("all_TTJetsP.root");
TFile *fSig     = TFile::Open("all_Hplus120.root");
double sf_ttbar = 1.00218; 
  
//SAVE HISTOS ON DISK
bool isSaveHisto = false;
///////////////////////////////////////////  

//*-------------------------------
//* get histo from root files
//*-------------------------------
TH1F* getHisto(TFile *histFile, TString dirBase, TString dirIso, TString dirBTag, TString histName, double sf){
  TH1F* hist;
  TString histPath("");
  if(histName=="cutflow") histPath = dirBase+"/"+dirIso+"/"+histName;
  else  histPath = dirBase+"/"+dirIso+"/"+dirBTag+"/"+histName;

  if(!(histFile->Get(histPath))){
    hist = (TH1F*)(fTT->Get(histPath))->Clone(histName);
    hist->Add(hist, -1);
  }else hist = (TH1F*)(histFile->Get(histPath))->Clone(histName);
  hist->Scale(sf);
  return hist;
}

//*-------------------------------
//* function for (Data -nonQCDBkg)
//*-------------------------------
TH1F * dataMCdiff(TString dirIso, TString dirBTag, TString histname, TString xaxis_title, int bin, bool log=false, bool axisrange=false, double xmin=0, double xmax=10){
  //hMC = all Bkg MC samples
  TH1F* hMC =  (TH1F*) getHisto(fVV, "baseLowMET", dirIso, dirBTag, histname, 1)->Clone("hMC");
  hMC->Add(getHisto(fDY, "baseLowMET", dirIso, dirBTag, histname, 1));
  hMC->Add(getHisto(fST, "baseLowMET", dirIso, dirBTag, histname, 1));
  hMC->Add(getHisto(fWJ, "baseLowMET", dirIso, dirBTag, histname, 1));
  hMC->Add(getHisto(fTT, "baseLowMET", dirIso, dirBTag, histname, sf_ttbar));
  
  //DATA
  TH1F* hData = (TH1F*)getHisto(fData, "baseLowMET", dirIso, dirBTag, histname, 1)->Clone("data");
  //Data - non QCD Bkg
  TH1F *hDiff = (TH1F*)hData->Clone("hDiff");
  hDiff->Add(hMC, -1);
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(1.0);
  //hDiff->GetYaxis()->SetRangeUser(0, 2);
  hDiff->GetXaxis()->SetRangeUser(xmin, xmax);
  hDiff->SetTitle("");
  //hDiff->GetXaxis()->SetTitle(xaxis_title); 
  hDiff->GetYaxis()->SetTitleOffset(0.65);
  hDiff->GetYaxis()->SetTitle("Data - nonQCDBkg"); hDiff->GetYaxis()->CenterTitle();
  hDiff->GetYaxis()->SetTitleSize(0.08); hDiff->GetXaxis()->SetTitleSize(0.08);
  hDiff->GetXaxis()->SetLabelSize(0.06); hDiff->GetXaxis()->LabelsOption("u"); // extra
  hDiff->GetYaxis()->SetLabelSize(0.06); hDiff->GetYaxis()->LabelsOption("u"); // extra
  ///hDiff->Draw("e1"); // use "P" or "AP"
  //hDiff->Draw("E same");
  return hDiff;
}

//*-------------------------------
//* Overlap of 2 (Data -nonQCDBkg)
//*-------------------------------

void dataMCdiffOverlap(TString dirIso, TString dirNoniso, TString dirBTag, TString histname, TString xaxis_title, double xmin=0, double xmax = 100){
  gStyle->SetOptStat(0);
  //TCanvas *c1 = new TCanvas("ddd", "aaa");
  //c1->Divide(2, 2);
  //c1->cd(postion);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  //legend box
  TLegend* leg = new TLegend(0.6318792,0.6261504,0.8012081,0.9198861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.1);
  
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.11,0.9354,0.90,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.08);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  TText *text = pt->AddText(dirBTag+": #sqrt{s} = 13 TeV,    34.94 fb^{-1}; ");
  //TText *text = pt->AddText(dir+":  CMS Preliminary,    #sqrt{s} = 13 TeV,    35.45 fb^{-1}; ");
  text->SetTextAlign(11);
  
  //pave text channel box
  TPaveText *ch = new TPaveText(0.953,0.9154898,0.6210067,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetTextSize(0.08);
  ch->SetBorderSize(1);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  
  //pave text histname
  TPaveText *hLable = new TPaveText(0.20,0.5254898,0.3210067,0.6262187,"brNDC");
  hLable->SetFillColor(19);
  hLable->SetFillStyle(0);
  hLable->SetLineColor(0);
  hLable->SetTextSize(0.08);
  hLable->SetBorderSize(1);
  hLable->AddText(xaxis_title);

  //data-MC from isolated region
  TH1F *hDiff = dataMCdiff(dirIso, dirBTag, histname, xaxis_title, 1, true, true, xmin, xmax);
  leg->AddEntry(hDiff,"qcd iso","P");
  hDiff->SetMarkerColor(kRed);
  hDiff->SetLineColor(kRed);
  hDiff->Scale(1/hDiff->Integral());
  hDiff->Draw("e1"); 
  //data-MC from non-isolated region
  TH1F *hDiff_NonIso = dataMCdiff(dirNoniso, dirBTag, histname, xaxis_title, 1, true, true, xmin, xmax);
  leg->AddEntry(hDiff_NonIso,"qcd non-iso","P");
  hDiff_NonIso->SetMarkerColor(kGreen);
  hDiff_NonIso->SetLineColor(kGreen);
  hDiff_NonIso->Scale(1/hDiff_NonIso->Integral());
  hDiff_NonIso->Draw("SAME");  
  pt->Draw();
  leg->Draw();
  ch->Draw();
  hLable->Draw();
  TF1 *baseLine = new TF1("baseLine","0", -100, 2000); 
  baseLine->SetLineColor(kBlue);
  baseLine->Draw("SAME");
  //c1->Update();
  
}

void qcd_sf(){
  TH1F *hDiff = dataMCdiff("Iso", " ", "cutflow", "cutflow", 1, true, true, 0, 15);
  TH1F *hDiff_NonIso = dataMCdiff("NonIso", " ", "cutflow", "cutflow", 1, true, true, 0, 15);
  int nbin= hDiff->GetSize();
  TString steps[15] = {"muon trig","= 1 muon","0 electron","muon SF","RelIso", "#geq 4 jets","#slash{E}_{T} #geq 20GeV", "MT >20", "#geq 2 b-jets","BTag SF", "fit converg","kfJetPt >25","dRJets <0.2",""};
  hDiff->Divide(hDiff_NonIso);
  std::cout << std::setprecision(4);
  for(int i=1; i<17; i++){
    double sf = hDiff->GetBinContent(i); 
    double sf_err = hDiff->GetBinError(i);
    cout<<sf<<"\t"<<" +- "<<sf_err<<"\t"<<i<<"\t"<<steps[i-1]<<endl;
  }
  hDiff->Draw("e1");
  hDiff->SetTitle("");
  hDiff->GetYaxis()->SetTitle("QCD scale factors");
}

//void qcd_shape_kfit(TString dir="KinFit"){
void qcd_shape_btag(TString dir="BTag"){
  TCanvas *c1 = new TCanvas(dir, dir);
  c1->Divide(2, 3);
  if(isMuChannel){
    c1->cd(1);
    dataMCdiffOverlap("Iso", "NonIso", dir, "pt_mu", "Pt^{#mu}", 0, 300);
    c1->cd(2);
    dataMCdiffOverlap("Iso", "NonIso", dir, "eta_mu", "#eta^{#mu}", -10, 10);
  }
  if(isEleChannel){
    c1->cd(1);
    dataMCdiffOverlap("Iso", "NonIso", dir, "pt_ele", "Pt^{e}", 0, 300);
    c1->cd(2);
    dataMCdiffOverlap("Iso", "NonIso", dir, "eta_ele", "#eta^{e}", -10, 10);
  }
  c1->cd(3);
  dataMCdiffOverlap("Iso", "NonIso", dir, "pt_jet", "Pt^{jets}", 0, 300);
  c1->cd(4);
  dataMCdiffOverlap("Iso", "NonIso", dir, "eta_jet", "#eta^{jets}", -10, 10);
  c1->cd(5);
  dataMCdiffOverlap("Iso", "NonIso", dir, "final_pt_met", "MET", 0, 50);
  c1->cd(6);
  dataMCdiffOverlap("Iso", "NonIso", dir, "wmt", "MT", 0, 200);
  /*
  c1->cd(7);
  dataMCdiffOverlap("/nvtx", "nvtx", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(8);
  dataMCdiffOverlap("/rhoAll", "#rho", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(9);
  if(dir=="BTag") dataMCdiffOverlap("/mjj", "mjj", "Iso/"+dir, "NonIso/"+dir, dir);
  else dataMCdiffOverlap("/mjj_kfit", "mjj_kfit", "Iso/"+dir, "NonIso/"+dir, dir);
  */
  c1->Update();
  if(isSaveHisto){
    TString outFile("$PWD/");
    //outFile += histname;
    if(isMuChannel)outFile += "mu_"+dir+".png";
    if(isEleChannel)outFile += "ele_"+dir+".png";
    c1->SaveAs(outFile);
    c1->Close();
  }
  //dataMCdiffOverlap("/wmt", "MT","Iso/"+dir, "NonIso/"+dir, dir);
}


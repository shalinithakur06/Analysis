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
 
//*-------------------------------
//* function to add histograms
//*-------------------------------
 void addHisto(TString path, TString filename, TString lable, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, TH1F* hMC){
  TFile* f2 = TFile::Open(path+filename);
  TH1F* h2_base = (TH1F*)(f2->Get("baseLowMET/"+histname))->Clone("h2_base");
  h2_base->Scale(scale);  
  h2_base->SetFillColor(color);
  if(axisrange) h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  hMC->Add(h2_base);
}

//*-------------------------------
//* function for (Data -nonQCDBkg)
//*-------------------------------
TH1F * dataMCdiff(TString histname, TString xaxis_title, int bin, bool log=false, bool axisrange=false, double xmin=0, double xmax=10){
  //path to the histo root files
  TString inFile("$PWD/");
  TFile* f = TFile::Open(inFile+"all_VV.root");
  TH1F* h1_base;
  if(!(f->Get("baseLowMET/"+histname))){
    TFile* f_tmp = TFile::Open(inFile+"all_TTJetsP.root");
    h1_base = (TH1F*)(f_tmp->Get("baseLowMET/"+histname))->Clone("h1_base");
    h1_base->Add(h1_base, -1);
  }else h1_base = (TH1F*)(f->Get("baseLowMET/"+histname))->Clone("h1_base");
  
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  addHisto(inFile, "all_DY.root", "DY", histname, 9, 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_ST.root", "ST", histname, 800 , 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_WJets.root", "WJets", histname, 6 , 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_TTJetsP.root","ttbar", histname, 433, 1, axisrange, xmin, xmax,  hMC);
  //DATA
  TFile* data_ = TFile::Open(inFile+"all_muData.root");
  //TFile* data_ = TFile::Open(inFile+"all_EleData.root");
  TH1F* data = (TH1F*)(data_->Get("baseLowMET/"+histname))->Clone("data");
  //Data - non QCD Bkg
  TH1F *hDiff = (TH1F*)data->Clone("hDiff");
  hDiff->Add(hMC, -1);
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(1.0);
  //hDiff->GetYaxis()->SetRangeUser(0, 2);
  //hDiff->GetXaxis()->SetRangeUser(xmin, xmax);
  hDiff->SetTitle("");
  hDiff->GetXaxis()->SetTitle(xaxis_title);
  hDiff->GetYaxis()->SetTitle("Data - nonQCDBkg");
  ///hDiff->Draw("e1"); // use "P" or "AP"
  //hDiff->Draw("E same");
  return hDiff;
}

//*-------------------------------
//* Overlap of 2 (Data -nonQCDBkg)
//*-------------------------------

void dataMCdiffOverlap(TString histname, TString xaxis_title, TString dir_iso, TString dir_Noniso, TString dir){
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
  leg->SetTextSize(0.03);
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.11,0.9354,0.90,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.04);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  TText *text = pt->AddText(dir+":  CMS Preliminary,    #sqrt{s} = 13 TeV,    35.45 fb^{-1}; ");
  text->SetTextAlign(11);
  //pave text channel box
  TPaveText *ch = new TPaveText(0.953,0.9154898,0.6210067,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetTextSize(0.045);
  ch->SetBorderSize(1);
  ch->AddText("#mu + jets");
  //data-MC from isolated region
  TH1F *hDiff = dataMCdiff(dir_iso+histname, xaxis_title, 1, true, true, 0, 10);
  leg->AddEntry(hDiff,"qcd iso, MET<20 GeV","P");
  hDiff->SetMarkerColor(kRed);
  hDiff->SetLineColor(kRed);
  hDiff->Scale(1/hDiff->Integral());
  hDiff->Draw("e1"); 
  //data-MC from non-isolated region
  TH1F *hDiff_NonIso = dataMCdiff(dir_Noniso+histname, xaxis_title, 1, true, true, 0, 10);
  leg->AddEntry(hDiff_NonIso,"qcd non-iso, MET<20 GeV","P");
  hDiff_NonIso->SetMarkerColor(kGreen);
  hDiff_NonIso->SetLineColor(kGreen);
  hDiff_NonIso->Scale(1/hDiff_NonIso->Integral());
  hDiff_NonIso->Draw("SAME");  
  pt->Draw();
  leg->Draw();
  ch->Draw();
  TF1 *baseLine = new TF1("baseLine","0", -100, 2000); 
  baseLine->SetLineColor(kBlack);
  baseLine->Draw("SAME");
  //c1->Update();
  
}

void qcd_sf(){
  TH1F *hDiff = dataMCdiff("Iso/cutflow", "cutflow", 1, true, true, 0, 15);
  TH1F *hDiff_NonIso = dataMCdiff("NonIso/cutflow", "cutflow", 1, true, true, 0, 15);
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

void qcd_shape_btag(TString dir="KinFit"){
//void qcd_shape_btag(TString dir="BTag"){
  TCanvas *c1 = new TCanvas(dir, dir);
  c1->Divide(2, 3);
  c1->cd(1);
  dataMCdiffOverlap("/pt_mu", "Pt^{#mu}", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(2);
  dataMCdiffOverlap("/eta_mu"," #eta^{#mu}","Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(3);
  dataMCdiffOverlap("/pt_jet","Pt^{jets}", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(4);
  dataMCdiffOverlap("/eta_jet", "#eta^{jet}", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(5);
  dataMCdiffOverlap("/final_pt_met", "MET", "Iso/"+dir, "NonIso/"+dir, dir);
  c1->cd(6);
  dataMCdiffOverlap("/wmt", "MT", "Iso/"+dir, "NonIso/"+dir, dir);
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
  bool savePDF = false;
  if(savePDF){
    TString outFile("/home/ravindra/lxplus/StackHisto/MuPlots/QCD/"+dir);
    //outFile += histname;
    outFile += "_mu"+dir+".pdf";
    c1->SaveAs(outFile);
    c1->Close();
  }
  //dataMCdiffOverlap("/wmt", "MT","Iso/"+dir, "NonIso/"+dir, dir);
}

void qcd_shape_check(TString dir="BTag"){
  TCanvas *c1 = new TCanvas(dir, dir);
  c1->Divide(1, 1);
  c1->cd(1);
  dataMCdiffOverlap("/pt_mu", "Pt^{#mu}", "Iso/"+dir, "NonIso/"+dir, dir);
}

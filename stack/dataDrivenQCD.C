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
  TH1F* h2_base = (TH1F*)(f2->Get("base/"+histname))->Clone("h2_base");
  h2_base->Scale(scale);  
  h2_base->SetFillColor(color);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  hMC->Add(h2_base);
}

//*-------------------------------
//* function for (Data -nonQCDBkg)
//*-------------------------------
TH1F * dataMCdiff(TString histname, TString xaxis_title, int bin, bool log=false, bool axisrange=false, double xmin=0, double xmax=10){

  //path to the histo root files
  TString inFile("$PWD/stack_20170919_SR_Ele_lowMET/");
  TFile* f1 = TFile::Open(inFile+"all_VV.root");
  TH1F* h1_base = (TH1F*)(f1->Get("base/"+histname))->Clone("h1_base");
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  hMC->Reset(); 

  //addHisto(inFile, "all_Hplus.root", "Hplus", histname, 2, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //addHisto(inFile, "all_QCD.root", "QCD", histname, 3, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  addHisto(inFile, "all_DY.root", "DY", histname, 9, 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_ST.root", "ST", histname, 800 , 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_WJets.root", "WJets", histname, 6 , 1, axisrange, xmin, xmax,  hMC);
  addHisto(inFile, "all_TTJetsP.root","ttbar", histname, 433, 1.0134, axisrange, xmin, xmax,  hMC);
  //DATA
  //TFile* data_ = TFile::Open(inFile+"all_muData.root");
  TFile* data_ = TFile::Open(inFile+"all_EleData.root");
  TH1F* data = (TH1F*)(data_->Get("base/"+histname))->Clone("data");
  
  //Data - non QCD Bkg
  TH1F *hDiff = (TH1F*)data->Clone("hDiff");
  int nbin= data->GetSize();
  double diff_=0.0; 
  double sigma=0.0;
  double binC_data= 0.0;
  double binC_mc = 0.0;
  
  //vector<double> sf;
  for(int i=1; i<nbin; i++){
    binC_data = data->GetBinContent(i);
    binC_mc = hMC->GetBinContent(i);
    if(binC_mc !=0&& binC_data!=0){
      diff_ = binC_data - binC_mc;
      sigma = sqrt(binC_data) + sqrt(binC_mc);
      hDiff->SetBinContent(i, diff_);
      hDiff->SetBinError(i, sigma);
      //if(qcdSF) sf.push_back(binC_data -binC_mc);
    }
  }
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(1.3);
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
void dataMCdiffOverlap(TString histname, TString dir_iso, TString dir_Noniso){
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  c1->cd(1);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  //legend box
  TLegend* leg = new TLegend(0.7818792,0.6261504,0.9312081,0.9198861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  //pave text box
  TPaveText *pt = new TPaveText(0.15,0.94,0.93,0.97, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.04);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  TText *text = pt->AddText("CMS Preliminary,    #sqrt{s} = 13 TeV,    34.31 fb^{-1}");
  text->SetTextAlign(11);
  //data-MC from isolated region
  TH1F *hDiff = dataMCdiff(dir_iso+histname, histname, 1, true, true, 0, 10);
  leg->AddEntry(hDiff,"Isolated","P");
  hDiff->SetMarkerColor(kRed);
  hDiff->SetLineColor(kRed);
  hDiff->Scale(1/hDiff->Integral());
  hDiff->Draw(); 
  
  //data-MC from non-isolated region
  TH1F *hDiff_NonIso = dataMCdiff(dir_Noniso+histname, histname, 1, true, true, 0, 10);
  leg->AddEntry(hDiff_NonIso,"Non-Isolated","P");
  hDiff_NonIso->SetMarkerColor(kGreen);
  hDiff_NonIso->SetLineColor(kGreen);
  hDiff_NonIso->Scale(1/hDiff_NonIso->Integral());
  hDiff_NonIso->Draw("SAME");  
  pt->Draw();
  leg->Draw();
  //c1->Update();
  TF1 *baseLine = new TF1("baseLine","0", -100, 2000); 
  baseLine->SetLineColor(kBlack);
  baseLine->Draw("SAME");
}

void qcd_shape_btag(TString dir="BTag/"){
  dataMCdiffOverlap("pt_ele", "Iso/"+dir, "NonIso/"+dir);
  dataMCdiffOverlap("eta_ele", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("pt_mu", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("eta_mu", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("pt_jet", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("eta_jet", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("nvtx", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("rhoAll", "Iso/"+dir, "NonIso/"+dir);
  //dataMCdiffOverlap("wmt", "Iso/"+dir, "NonIso/"+dir);
}

void qcd_shape_kf(TString dir="KinFit/"){
  dataMCdiffOverlap("pt_ele", "Iso/"+dir, "NonIso/"+dir);
  dataMCdiffOverlap("eta_ele", "Iso/"+dir, "NonIso/"+dir);
  dataMCdiffOverlap("pt_jet", "Iso/"+dir, "NonIso/"+dir);
  dataMCdiffOverlap("eta_jet", "Iso/"+dir, "NonIso/"+dir);
  dataMCdiffOverlap("final_pt_met", "Iso/"+dir, "NonIso/"+dir);
}
  

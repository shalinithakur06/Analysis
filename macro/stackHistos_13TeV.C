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
bool isMuChannel = true;
bool isEleChannel = false;

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY_M50.root");
//TFile* fDY	= TFile::Open("all_DY_Pt.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
//TFile *fSig     = TFile::Open("all_ExLepEEZ_M250.root");
TFile *fSig     = TFile::Open("all_ExLep_M250.root");

//SAVE HISTOS ON DISK
//bool isSaveHisto = true;
bool isSaveHisto = true;
///////////////////////////////////////////  

//--------------------------------------------//
//stack histos
//--------------------------------------------//
void example_stack(TString dir, TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=true, bool axisrange=false, double xmin=0, double xmax=10, bool label=false, double unc=false){
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  //TCanvas *c1 = new TCanvas();
  TCanvas *c1 = new TCanvas("canvas_"+histname, "canvas_"+histname, 800, 900);
  const float xpad[2] = {0.0,1.0};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  if(ratio){
    c1->Divide(1,2, 0, 0); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    gPad->SetLogy(log);
  }
  
  //Legends
  TLegend* leg = new TLegend(0.7518792,0.5061504,0.8512081,0.8898861,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.05);
  //Data
  TH1F* data = getHisto(fData, "base/Iso", dir, histname);
  //TH1F* data = (TH1F*)(fData->Get("base/Iso/"+dir+"/"+histname))->Clone("data");
  //data->SetBinErrorOption(TH1::kPoisson);
  data->SetMarkerStyle(20); data->SetMarkerSize(0.8);
  if(axisrange) data->GetXaxis()->SetRangeUser(xmin,xmax);
  data->SetFillColor(kBlack);
  //data->GetYaxis()->SetTitleOffset(1.35);
  ///data->Rebin(bin);
  data->SetTitle("");
  if(dir=="") data->SetAxisRange(1.0, 1.0e09 ,"Y");
  //else data->SetAxisRange(1.0, 1.0e6 ,"Y");
  if(drawdata)data->Draw("E");
  if(drawdata)leg->AddEntry(data,"Data","PE");
  data->GetYaxis()->SetTitle("Events");
  //data->GetYaxis()->SetRangeUser(1, 3000);
  data->GetYaxis()->SetRangeUser(1, 1.4* data->GetMaximum());
  //data->GetXaxis()->SetTitle(xaxis_title);
  data->GetYaxis()->SetTitleOffset(1.10);
  data->GetYaxis()->SetNdivisions(5);
  if(!label) data->GetXaxis()->SetNdivisions(5);
  data->GetXaxis()->SetTitleOffset(1.00);
  data->GetYaxis()->SetTitleSize(0.07);   data->GetXaxis()->SetTitleSize(0.07);
  data->GetXaxis()->SetLabelSize(0.07);   
  data->GetYaxis()->SetLabelSize(0.07);   
  data->GetXaxis()->SetTickLength(0.05);
  data->GetYaxis()->SetTickLength(0.04);
 
  //VV is the base histo
  TH1F* h1_base = getHisto(fVV, "base/Iso", dir, histname);
  h1_base->SetFillColor(kOrange+1);
  if(axisrange) h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  leg->AddEntry(h1_base,"VV","F");
  //Define stacked histo
  THStack* MuptStack = new THStack("MuptStack","");
  MuptStack->Add(h1_base);
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  stackHisto(fWJ, "W+ jets", dir, histname, kViolet + 1 , 1.0, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fTT,"t#bar{t}", dir, histname, kAzure +10, 1.0, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fDY, "Z/#gamma+jets", dir, histname, kSpring -9, 1.0, axisrange, xmin, xmax, MuptStack, hMC, leg);
  //Signal 
  TH1F* hSig = getHisto(fSig, "base/Iso", dir, histname);
  hSig->Scale(10);
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  
  //Draw stacked histo
  if(drawdata){
    //gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    MuptStack->Draw("HISTSAME");
  }
  else {
  MuptStack->Draw("HIST");
  //gPad->SetBottomMargin(0.20);
  gPad->SetTopMargin(0.10);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  }
  MuptStack->SetTitle(""); 
  //MuptStack->SetMaximum(1.1*(MuptStack->GetMaximum()));
  //MuptStack->SetMaximum(1.2*(hSig->GetMaximum()));
  MuptStack->GetYaxis()->SetTitle("Events");
  //MuptStack->GetYaxis()->SetRangeUser(1, 100);
  MuptStack->GetYaxis()->SetRangeUser(1, 1.4* MuptStack->GetMaximum());
  //MuptStack->GetXaxis()->SetTitle(xaxis_title);
  MuptStack->GetYaxis()->SetTitleOffset(1.10);
  MuptStack->GetXaxis()->SetTitleOffset(1.00);
  MuptStack->GetYaxis()->SetTitleSize(0.07);   
  MuptStack->GetXaxis()->SetTitleSize(0.07);
  MuptStack->GetXaxis()->SetLabelSize(0.07);   
  MuptStack->GetYaxis()->SetLabelSize(0.07);   
  MuptStack->GetXaxis()->SetTickLength(0.05);
  MuptStack->GetYaxis()->SetTickLength(0.04);
  
  if(axisrange) MuptStack->GetXaxis()->SetRangeUser(xmin,xmax);
  MuptStack->GetXaxis()->SetTitle(xaxis_title);
  MuptStack->GetYaxis()->SetTitle("Events");
  if(!label) MuptStack->GetXaxis()->SetNdivisions(5);
  MuptStack->SetMinimum(1.0);
  MuptStack->Draw("HISTSAME");
  if(drawdata)data->Draw("SAME"); 
  
  if(drawsignal)hSig->Draw("HISTSAME"); //  
  leg->AddEntry(hSig, "10*Signal","L"); 
  leg->Draw();
  
  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //signal
  //TPaveText *cct = paveText(0.1513423,0.7754898,0.4010067,0.8962187, 0, 19, 1, 0, 132);
  TPaveText *cct = paveText(0.492462,0.768605,0.742462,0.889247, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.049);
  //cct->AddText("M_{H^{+}} = 120 GeV");
  cct->AddText("M_{l^{*}} = 250 GeV");
    ///cct->AddText("Br(t#rightarrow H^{+}b) = 0.1");
  
  //hist name
  TPaveText *hLable = paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.080);
  hLable->AddText(xaxis_title);
  
  //channel
  TPaveText *ch = paveText(0.853,0.9154898,0.9010067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.08);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = paveText(0.09,0.9354,0.88,0.9362, 0, 19, 1, 0, 132);
  if(drawdata) pt->SetTextSize(0.059);
  else pt->SetTextSize(0.039);
  TText *text = pt->AddText(dir+": #sqrt{s} = 13 TeV, 41.9 fb^{-1}");
  //TText *text = pt->AddText(dir+": CMS Preliminary, #sqrt{s} = 13 TeV, 35.12 fb^{-1}");
  text->SetTextAlign(11);
  pt->Draw();
  if(drawsignal) cct->Draw();
  ch->Draw();
  //hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(ratio){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.5); //gPad->SetGridy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    //gPad->SetPad(xpad[0],0.05,xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)data->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(data);
    if(!label) hRatio->GetYaxis()->SetNdivisions(5);
    //hRatio->Add(hMC, -1);
    hRatio->Divide(hMC); hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(0.8);
    hRatio->SetMarkerColor(kBlack); hRatio->SetLineColor(kBlack); hRatio->GetYaxis()->SetRangeUser(0., 2.0);
    //hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetXaxis()->SetTickLength(0.13);
    hRatio->GetYaxis()->SetTickLength(0.04);
    hRatio->GetXaxis()->SetTitle(xaxis_title);
    hRatio->GetYaxis()->SetTitleOffset(0.450);
    hRatio->GetXaxis()->SetTitleOffset(1.00);
    hRatio->GetYaxis()->SetTitle("#frac{Data}{Bkg}"); hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitleSize(0.15); hRatio->GetXaxis()->SetTitleSize(0.20);
    hRatio->GetXaxis()->SetLabelSize(0.15); 
    hRatio->GetYaxis()->SetLabelSize(0.13); 
  //lable x-axis, for cutflow
  if(label){
    TString steps[7]={"muon trig","= 2 muon", "ControlP","bJetVeto","preSel", "Z-Tag",""}; 
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      hRatio->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      //hRatio->GetXaxis()->SetLabelSize(0.1);
      hRatio->GetXaxis()->LabelsOption("u");
    }
  }
    hRatio->Draw("E"); // use "P" or "AP"
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    TString outFile("$PWD/");
    outFile += dir+"/"+histname;
    if(isMuChannel) outFile += "_mu"+dir+".png";
    if(isEleChannel) outFile += "_ele"+dir+".png";
    c1->SaveAs(outFile);
    ///c1->Close();
  }
}

TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  if(!(histFile->Get(histPath+"/"+dir+"/"+histName))){
    hist = (TH1F*)(fTT->Get(histPath+"/"+dir+"/"+histName))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(histPath+"/"+dir+"/"+histName))->Clone(histName);
  return hist;
}
void stackHisto(TFile *filename, TString lable, TString dir, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, THStack* MuptStack, TH1F* hMC, TLegend* leg){
  TH1F* h2_base = getHisto(filename, "base/Iso", dir, histname);
  //h2_base->Draw();
  h2_base->Scale(scale);  
  h2_base->SetFillColor(color);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  leg->AddEntry(h2_base,lable,"F");
  MuptStack->Add(h2_base);
  hMC->Add(h2_base);
}



TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}
void example_stack_all(){
  example_stack("", "cutflow", "cutflow", 1, true, true, true, true, true, 0, 10, true, false);
  //example_stack("final_RelIso_mu", "RelIso of muons", 1, true, true, true, false, false, 0, 1, false, dir);
}

void example_stack_ControlP(TString dir="ControlP"){
//void example_stack_PreSel(TString dir="bJetVeto"){
//void example_stack_PreSel(TString dir="PreSel"){
//void example_stack_ZTag(TString dir="ZTag"){

    example_stack("BTag","multi_bjet","N^{b-jet}",1,true,true, true,false,true,       1,      5,  false, false);
    /*
  if(isMuChannel){
    example_stack(dir,"pt_1stMu","Pt^{#mu_{1}} [GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"pt_2ndMu","Pt^{#mu_{2}} [GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"mll","M_{#mu_{1}#mu_{2}} [GeV]",1,true,true,true,true, true,               50,    1000.0, false, false);
    example_stack(dir,"eta_1stMu","#eta^{#mu_{1}}", 1, true,true,true,true,true,       -3.0,   4.5,  false, false);
    example_stack(dir,"eta_2ndMu","#eta^{#mu_{2}}", 1, true,true,true,true,true,       -3.0,   4.5,  false, false);
    example_stack(dir,"multi_mu","N^{#mu}",1,true,true, true,false,true,       1,      5,  false, false);
    example_stack(dir,"pt_Z","Pt_{#mu_{1}#mu_{2}}[GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"eta_Z","#eta^{#mu_{1}#mu_{2}}", 1, true,true,true,true, true, -3.0,   4.5,  false, false);
  }
  else{
    example_stack(dir,"pt_1stEle","Pt^{e_{1}} [GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"pt_2ndEle","Pt^{e_{2}} [GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"mll","M_{e_{1}e_{2}} [GeV]",1,true,true,true,true, true,               50,    1000.0, false, false);
    example_stack(dir,"eta_1stEle","#eta^{e_{1}}", 1, true,true,true,true,true,       -3.0,   4.5,  false, false);
    example_stack(dir,"eta_2ndEle","#eta^{e_{2}}", 1, true,true,true,true,true,       -3.0,   4.5,  false, false);
    example_stack(dir,"multi_Ele","N^{e}",1,true,true, true,false,true,       1,      5,  false, false);
    example_stack(dir,"pt_Z","Pt_{e_{1}e_{2}}[GeV]", 1, true,true,true,true, true,   0.0,    1000.0, false, false);
    example_stack(dir,"eta_Z","#eta^{e_{1}e_{2}}", 1, true,true,true,true, true, -3.0,   4.5,  false, false);
  }
  example_stack(dir,"pt_jet","Pt^{jets} [GeV]", 1, true,true,true,true, true, 170,    1000,  false, false);
  example_stack(dir,"eta_jet","#eta^{jets}", 1, true,true,true,true,true,     -3.0,   4.5,  false, false);
  example_stack(dir,"final_multi_jet","N^{jets}",1,true,true,true,false,true,   3,      15,  false, false);
  example_stack(dir,"nvtx","N^{vertex}",1,true,true,true,true, true,           0.0,    60.0,  false, false);
  example_stack(dir,"rhoAll","#rho",1,true,true,true,true, true,               0.0,    70.0,  false, false);
  if(dir=="ZTag"){
   if(isMuChannel){
     example_stack(dir,"mlZ_max","M_{#mu Z}^{max}",1,true,false,false,false, true,                 100,    1000.0, false, false);
     example_stack(dir,"mlZ_min","M_{#mu Z}^{min}",1,true,false,false,false, true,                 0.0,    1000.0, false, false);
   }
   else{
     example_stack(dir,"mlZ_max","M_{e Z}^{max}",1,true,false,false,false, true,                 100,    1000.0, false, false);
     example_stack(dir,"mlZ_min","M_{e Z}^{min}",1,true,false,false,false, true,                 0.0,    1000.0, false, false);
   }
  }
  */
}


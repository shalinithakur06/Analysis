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
//else(isEleChannel)fData = TFile::Open("all_EleData.root");

TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fQCD	= TFile::Open("all_QCD.root");
TFile* fST	= TFile::Open("all_ST.root");
TFile* fTT	= TFile::Open("all_TTJetsP.root");
TFile *fSig     = TFile::Open("all_Hplus120.root");
double sf_ttbar = 1.00218; 

//USER'S INPUT FOR DATA DRIVEN QCD 
bool isDataDrivenQCD = false;
double qcd_sf_btag =  2.204 ; 
double qcd_sf_kfit = 2.380 ;
double qcd_sf_btag_err = 0.1416;
double qcd_sf_kfit_err = 0.2842;
TFile *f_QCD_dd = new TFile("all_QCD_dd.root","RECREATE");
  
//SAVE HISTOS ON DISK
bool isSaveHisto = false;
///////////////////////////////////////////  

//--------------------------------------------//
//various functions
//--------------------------------------------//
void stackHisto(TFile *filename, TString lable, TString dir, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, THStack* MuptStack, TH1F* hMC, TLegend* leg);
//function to add histograms
TH1F* addHistoForUnc(TString dir, TString histname, TString sys, bool isDataDrivenQCD = false, double qcd_sf=1);
TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );
//get histogram from root file. Return empty hist, if the hist does not exit.
TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName);

//--------------------------------------------//
//stack histos
//--------------------------------------------//
void example_stack(TString dir, TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10, bool label=false, double unc=false){
  //Pad
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  if(ratio){
    c1->Divide(1,2); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    gPad->SetLogy(log);
  }
  
  //Legends
  TLegend* leg = new TLegend(0.7518792,0.6061504,0.9312081,0.8898861,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.03);
  //Data
  TH1F* data = getHisto(fData, "base/Iso", dir, histname);
  //TH1F* data = (TH1F*)(fData->Get("base/Iso/"+dir+"/"+histname))->Clone("data");
  //data->SetBinErrorOption(TH1::kPoisson);
  data->SetMarkerStyle(20); data->SetMarkerSize(0.8);
  if(axisrange) data->GetXaxis()->SetRangeUser(xmin,xmax);
  data->SetFillColor(kBlack);
  data->GetYaxis()->SetTitleOffset(1.35);
  data->GetYaxis()->SetTitle("Events");
  data->GetXaxis()->SetTitle(xaxis_title);
  ///data->Rebin(bin);
  data->SetTitle("");
  if(label) data->SetAxisRange(1.0, 1.0e9 ,"Y");
  else data->SetAxisRange(1.0, 1.0e6 ,"Y");
  if(drawdata)data->Draw("E"); 
  if(drawdata)leg->AddEntry(data,"Data","PE"); 
  
  //VV is the base histo
  TH1F* h1_base = getHisto(fVV, "base/Iso", dir, histname);
  h1_base->SetFillColor(13);
  if(axisrange) h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  leg->AddEntry(h1_base,"VV","F");
  //Define stacked histo
  THStack* MuptStack = new THStack("MuptStack","");
  MuptStack->Add(h1_base);
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  
  stackHisto(fQCD, "QCD", dir, histname, 3, 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fDY, "Z+jets", dir, histname, 9, 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fST, "Single t", dir, histname, 800 , 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fWJ, "W+ jets", dir, histname, 6 , 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fTT,"t#bar{t}", dir, histname, 433, sf_ttbar, axisrange, xmin, xmax, MuptStack, hMC, leg);

  //Signal 
  TH1F* hSig = getHisto(fSig, "base/Iso", dir, histname);
  //TH1F* hSig = (TH1F*)(fSig->Get("base/Iso/"+dirhistname))->Clone("Signal");
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  
  //Draw stacked histo
  if(drawdata) MuptStack->Draw("HISTSAME");
  else MuptStack->Draw("HIST");
  MuptStack->SetTitle(""); 
  //MuptStack->SetMaximum(1.1*(MuptStack->GetMaximum()));
  MuptStack->SetMaximum(1.2*(hSig->GetMaximum()));
  
  if(axisrange) MuptStack->GetXaxis()->SetRangeUser(xmin,xmax);
  MuptStack->GetXaxis()->SetTitle(xaxis_title);
  MuptStack->GetYaxis()->SetTitle("Events");
  MuptStack->GetYaxis()->SetTitleOffset(1.45);
  MuptStack->SetMinimum(1.0);
  MuptStack->Draw("HISTSAME");
  if(drawdata)data->Draw("SAME"); 
  
  //hSigBkg = Bkg MC+ signal MC
  TH1F* hSigBkg = (TH1F*)hMC->Clone("hSigBkg"); 
  hSigBkg->Add(hSig); hSigBkg->SetLineStyle(2);
  hSigBkg->SetLineWidth(3); hSigBkg->SetFillColor(0);
  if(drawsignal)hSig->Draw("HISTSAME"); // only for hSig histogram 
  //if(drawsignal)hSigBkg->Draw("HISTSAME"); //  
  leg->AddEntry(hSig, "Signal","L"); 
  leg->Draw();
  
  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //signal
  TPaveText *cct = paveText(0.1513423,0.7754898,0.4010067,0.8962187, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.049);
  cct->AddText("M_{H^{+}} = 120 GeV");
  ///cct->AddText("Br(t#rightarrow H^{+}b) = 0.1");
  
  //hist name
  TPaveText *hLable = paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.080);
  hLable->AddText(xaxis_title);
  
  //channel
  TPaveText *ch = paveText(0.953,0.9154898,0.6210067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.08);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = paveText(0.09,0.9354,0.88,0.9362, 0, 19, 1, 0, 132);
  if(drawdata) pt->SetTextSize(0.059);
  else pt->SetTextSize(0.039);
  TText *text = pt->AddText(dir+": CMS Preliminary, #sqrt{s} = 13 TeV, 34.94 fb^{-1}");
  text->SetTextAlign(11);
  pt->Draw();
  if(drawsignal) cct->Draw();
  ch->Draw();
  hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(ratio){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.3); gPad->SetGridy();
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    //gPad->SetPad(xpad[0],0.05,xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)data->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(data);
    //hRatio->Add(hMC, -1);
    hRatio->Divide(hMC); hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(0.8);
    hRatio->SetMarkerColor(kBlack); hRatio->SetLineColor(kBlack); hRatio->GetYaxis()->SetRangeUser(0, 2);
    //hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetXaxis()->SetTitle(xaxis_title); hRatio->GetYaxis()->SetTitleOffset(0.45);
    hRatio->GetYaxis()->SetTitle("Data/Bkg"); hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitleSize(0.1); hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetLabelSize(0.12); hRatio->GetXaxis()->LabelsOption("u"); // extra
    hRatio->GetYaxis()->SetLabelSize(0.08); hRatio->GetXaxis()->LabelsOption("u"); // extra
  //lable x-axis, for cutflow
  if(label){
    TString steps[15]={"electron trig","= 1 electron","0 muon","electron SF","RelIso < 0.8", "#geq 4 jets","#slash{E}_{T} #geq 20 GeV", "MT >20 GeV", "#geq 2 b-jets","BTag SF", "fit converges","kfJetPt >25 GeV","dRJets <0.2",""}; 
    //if(isMuChannel) steps[15] = {"muon trig","= 1 muon","0 electron","muon SF","RelIso", "#geq 4 jets","#slash{E}_{T} #geq 20 GeV", "MT >20 GeV", "#geq 2 b-jets","BTag SF", "fit converges","kfJetPt >25 GeV","dRJets <0.2",""};
    //if(isEleChannel) steps[15] = {"electron trig","= 1 electron","0 muon","electron SF","RelIso < 0.8", "#geq 4 jets","#slash{E}_{T} #geq 20 GeV", "MT >20 GeV", "#geq 2 b-jets","BTag SF", "fit converges","kfJetPt >25 GeV","dRJets <0.2",""};
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      hRatio->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      hRatio->GetXaxis()->SetLabelSize(0.1);
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
    c1->Close();
  }
}

TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  if(!(histFile->Get(histPath+"/"+dir+"/"+histName))){
    hist = (TH1F*)(fTT->Get(histPath+"/"+dir+"/"+histName))->Clone(histName);
    hist->Add(hist, -1);
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
  example_stack("", "cutflow", "cutflow", 1, true, true, true, true, true, 0, 15, true, true);
  //example_stack("final_RelIso_mu", "RelIso of muons", 1, true, true, true, false, false, 0, 1, false, dir);
}

//void example_stack(TString dir, TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10, bool label=false, double unc=false){

void example_stack_btag(TString dir="BTag"){//dir=KinFit, BTag
//void example_stack_kfit(TString dir="KinFit"){//dir=KinFit, BTag
 example_stack(dir,"pt_mu","Pt^{#mu}[GeV]", 1, true,true,true,true, true,   0.0,    500.0, false, true);
 example_stack(dir,"eta_mu","#eta^{#mu}", 1, true,true,true,true,true,       -2.5,   3.5,  false, true);
 example_stack(dir,"pt_jet","Pt^{jets}[GeV]", 1, true,true,true,true, true,  0.0,    500,  false, true);
 example_stack(dir,"eta_jet","#eta^{jets}", 1, true,true,true,true,true,     -2.5,   3.5,  false, true);
 example_stack(dir,"final_multi_jet","N^{jets}",1,true,true,true,true,true,   3,      15,  false, true);
 example_stack(dir,"final_pt_met","MET", 1, true,true,true,true, true,        0.0,    500, false, true);
 example_stack(dir,"nvtx","N^{vertex}",1,true,true,true,true, true,           0.0,    70.0,  false, true);
 example_stack(dir,"rhoAll","#rho",1,true,true,true,true, true,               0.0,    70.0,  false, true);
 example_stack(dir,"wmt","MT[GeV]",1,true,true,true,true,true,                0,      300,   false, true);
} 


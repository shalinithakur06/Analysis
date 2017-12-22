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

void stackHisto(TString path, TString filename, TString lable, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, THStack* MuptStack, TH1F* hMC, TH1F* hSig, TLegend* leg){
  TFile* f2 = TFile::Open(path+filename);
  if(f2 == 0) return;
  if(f2->IsZombie()){f2->Close(); return;}
  TH1F* h2_base = (TH1F*)(f2->Get("base/Iso/"+histname))->Clone("h2_base");
  h2_base->Scale(scale);  
  h2_base->SetFillColor(color);
  if(axisrange){
    h2_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  leg->AddEntry(h2_base,lable,"F");
  MuptStack->Add(h2_base);
  hMC->Add(h2_base);
  hSig->Add(h2_base);
}


void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10, TString dir="BTag", bool label=false)
{
  ////////////////////////////////////////////////////////////
  // 
  // Basic decorations for stacked plots
  //
  ////////////////////////////////////////////////////////////
  
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  const float xpad[2] = {0.,1};
  //const float ypad[4] = {0.,0.3,0.3,1.0};
  const float ypad[4] = {0.,0.2351916,0.2351916,0.98};
  if(ratio){
    c1->Divide(1,2);
    c1->cd(1);
    //gPad->SetBottomMargin(0);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    //p->SetGridx();
    //p->SetGridy();
    gPad->SetLogy(log);
  }
  
  //TLegend* leg = new TLegend(0.7818792,0.6261504,0.9312081,0.9198861,NULL,"brNDC");
  TLegend* leg = new TLegend(0.7518792,0.6061504,0.9312081,0.8898861,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);

  TPaveText *pt = new TPaveText(0.09,0.9354,0.88,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetLineColor(0);
  pt->SetTextFont(132);

  ////////////////////////////////////////////////////////////
  // 
  // Stack all the other MC histos on top of the base histo
  //
  ////////////////////////////////////////////////////////////

  TString inFile("$PWD/stack_20170919_SR_Ele_highMET/");
  
  //VV is the base histo
  //TFile* f1 = TFile::Open(inFile+"WW_Merged.root");
  TFile* f1 = TFile::Open(inFile+"all_VV.root");
  if(f1 == 0) return;
  if(f1->IsZombie()){f1->Close(); return;}
  TH1F* h1_base = (TH1F*)(f1->Get("base/Iso/"+histname))->Clone("h1_base");
  double scale_factor = 1.0;
  h1_base->Scale(scale_factor);  
  h1_base->SetFillColor(13);
  if(axisrange){
    h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  leg->AddEntry(h1_base,"VV","F");
  //Define stacked histo
  THStack* MuptStack = new THStack("MuptStack","");
  MuptStack->Add(h1_base);
  
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  hMC->Reset();  
  
  //hSig = Bkg MC+ signal MC
  TH1F* hSig = (TH1F*)h1_base->Clone("hSig"); 
  hSig->Reset(); 
  hSig->SetLineColor(kCyan);
  hSig->Add(h1_base);
  hSig->SetLineColor(kMagenta+0);
  hSig->SetLineStyle(2);
  hSig->SetLineWidth(3);
  hSig->SetFillColor(0);
  
  //---------------------------
  // QCD from Data
  //---------------------------
  TFile* fVV = TFile::Open(inFile+"all_VV.root");
  TFile* fDY = TFile::Open(inFile+"all_DY.root");
  TFile* fST = TFile::Open(inFile+"all_ST.root");
  TFile* fWJ = TFile::Open(inFile+"all_WJets.root");
  TFile* fTT = TFile::Open(inFile+"all_TTJetsP.root");
  TFile* fDa = TFile::Open(inFile+"all_EleData.root");
  
  TH1F* hVV = (TH1F*)(fVV->Get("base/NonIso/"+histname));
  TH1F* hDY = (TH1F*)(fDY->Get("base/NonIso/"+histname));
  TH1F* hST = (TH1F*)(fST->Get("base/NonIso/"+histname));
  TH1F* hWJ = (TH1F*)(fWJ->Get("base/NonIso/"+histname));
  TH1F* hTT = (TH1F*)(fTT->Get("base/NonIso/"+histname));
  TH1F* hDa = (TH1F*)(fDa->Get("base/NonIso/"+histname));
  
  TH1F* hOtherMC = (TH1F*)hVV->Clone("hOtherMC"); 
  hOtherMC->Add(hDY); 
  hOtherMC->Add(hST); 
  hOtherMC->Add(hWJ);
  hTT->Scale(1.013/0.9865);
  hOtherMC->Add(hTT); 
  TH1F* hQCD = (TH1F*)hDa->Clone("QCD from Data"); 
  hQCD->Add(hOtherMC, -1);
  hQCD->Scale(1.21);
  if(axisrange){
    hQCD->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  hQCD->SetFillColor(kGreen);
  //hQCD->Draw();
  
  //---------------------------
  // QCD and other Bkg from MC
  //---------------------------
  //COLOR SCHEME: VIBGYOR : V=kViolet, I=kIndigo, B=kBlue, G=kGreen, Y=kYellow, O=kOrange, Red= kRed    
  //stackHisto(inFile, "all_Hplus.root", "Hplus", histname, 2, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //data-driven QCD
  bool dataDrivenQCD = false;
  if(dataDrivenQCD){
    leg->AddEntry(hQCD,"QCD","F");
    MuptStack->Add(hQCD);
    hMC->Add(hQCD);
  }else stackHisto(inFile, "all_QCD.root", "QCD", histname, 3, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_DY.root", "DY", histname, 9, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_ST.root", "ST", histname, 800 , 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_WJets.root", "WJets", histname, 6 , 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  //stackHisto(inFile, "all_TTJetsM.root","ttbar", histname, 433, 1, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  stackHisto(inFile, "all_TTJetsP.root","TTJets", histname, 433, 1.00/0.9865, axisrange, xmin, xmax, MuptStack, hMC, hSig,leg);
  
  ////////////////////////////////////////////////////////////
  //
  // Plot DATA points on the stacked MC histos
  //
  ////////////////////////////////////////////////////////////

  TFile* data_ = TFile::Open(inFile+"all_EleData.root");
  if(data_ == 0) return;
  if(data_->IsZombie()){data_->Close(); return;}
  //TH1F* data = (data_->Get("base/BTag/"+histname))->Clone("data");
  TH1F* data = (TH1F*)(data_->Get("base/Iso/"+histname))->Clone("data");
  ///data->SetBinErrorOption(TH1::kPoisson); // added this line for Asymmetric Error Bars for Poisson Event Counts
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.1);
  data->Rebin(bin);
  if(axisrange){
    data->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  data->SetFillColor(kBlack);
  
  //lable x-axis, for cutflow
  if(label){
    TString steps[15] = {"no cut","ele trig","= 1 electron","0 muon","electron SF","RelIso <0.08", "#geq 4 jets","#slash{E}_{T} #geq 20GeV", "#geq 2 b-jets","BTag SF","MT >0 GeV","Top Pt Rwt","KinFit noCut","KinFit cut",""};
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      data->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      data->GetXaxis()->SetLabelSize(0.06);
      data->GetXaxis()->LabelsOption("u");
    }
  }
    
  
  if(drawdata)leg->AddEntry(data,"Data","PE"); 
  if(drawsignal)leg->AddEntry(hSig,"Sig+bkgs","L"); // only for hSig histogram
  TPaveText *cct = new TPaveText(0.2013423,0.7754898,0.4010067,0.8962187,"brNDC");
  cct->SetFillColor(19);
  cct->SetFillStyle(0);
  cct->SetLineColor(0);
  cct->SetBorderSize(1);
  cct->AddText("M_{H^{+}} = 120 GeV");
  ///cct->AddText("Br(t#rightarrow H^{+}b) = 0.1");
  TPaveText *ch = new TPaveText(0.953,0.9154898,0.6210067,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetBorderSize(1);
  ch->AddText("e + jets");
  pt->SetTextSize(0.059);
  TText *text = pt->AddText(dir+":   CMS Preliminary,    #sqrt{s} = 13 TeV,    35.32 fb^{-1}; ");
  //TText *text = pt->AddText("#sqrt{s} = 13 TeV");
  text->SetTextAlign(11);
  if(histname.Contains("Pre_RelIso") || histname.Contains("cutflow")) {
   data->SetAxisRange(1.0, 5.0e9 ,"Y");
  }else{
   data->SetAxisRange(1.0, 1.0e6 ,"Y");
  }

  data->GetYaxis()->SetTitleOffset(1.45);
  data->GetYaxis()->SetTitle("Events");
  data->GetXaxis()->SetTitle(xaxis_title);
  if(drawdata){
    data->Draw("E"); // not for dijet mass
    MuptStack->Draw("HISTSAME"); // no histsame for dijet calculation
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->SetMinimum(1.0);
    MuptStack->SetMaximum(1.1*(MuptStack->GetMaximum()));
    if(axisrange){
      MuptStack->GetXaxis()->SetRangeUser(xmin,xmax);
    }
    //stack->GetXaxis()->SetLimits(xmin, xmax);
    //stack->Draw("H"); // make sure it's redrawn
    MuptStack->Draw("HIST");
    MuptStack->GetXaxis()->SetTitle(xaxis_title);
    MuptStack->GetYaxis()->SetTitle("Events");
    MuptStack->GetYaxis()->SetTitleOffset(1.45);
    MuptStack->Draw("HISTSAME");
  }
  if(drawsignal)hSig->Draw("HISTSAME"); // only for hSig histogram 
  //  if(drawsignal)h7_base->Draw("HISTSAME"); // for sv_mass plot only
  if(drawdata)data->Draw("ESAME"); // not for dijet mass 
  MuptStack->SetMinimum(1.0);
  leg->Draw();
  pt->Draw();
  if(drawsignal) cct->Draw();
  ch->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  ////////////////////////////////////////////////////////////
  //
  // DATA - Bkg
  //
  ////////////////////////////////////////////////////////////

  if(ratio){
    c1->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3);
    gPad->SetGridy();
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    //gPad->SetPad(xpad[0],0.05,xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)data->Clone("hRatio");
    //hRatio->Reset();
    //hRatio->Add(data);
    //hRatio->Add(hMC, -1);
    hRatio->Divide(hMC);
    /*
    int nbin= data->GetSize();
    double div=0.0; 
    double sigma=0.0;
    double a1 =0.0;
    double binC_data= 0.0;
    double binC_mc = 0.0;

    for(int i=1; i<nbin; i++){
      binC_data = data->GetBinContent(i);
      binC_mc = hMC->GetBinContent(i);
      if(binC_mc !=0&& binC_data!=0){
        div = binC_data/binC_mc;
        a1 = sqrt(1.0/binC_data + 1.0/binC_mc);
        sigma = div*a1;
        hRatio->SetBinContent(i, div);
        hRatio->SetBinError(i, sigma);
      }
    }
    */
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.1);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetRangeUser(0, 2);
    //hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetXaxis()->SetTitle(xaxis_title);
    hRatio->GetYaxis()->SetTitleOffset(0.55);
    hRatio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetLabelSize(0.12); // 0.1
    hRatio->GetXaxis()->LabelsOption("u"); // extra
    hRatio->GetYaxis()->SetLabelSize(0.06);
    hRatio->Draw("e1"); // use "P" or "AP"
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    //hRatio->Draw("E same");
    c1->Update();
  }

  bool savePDF = false;
  if(savePDF){
    TString outFile("/home/ravindra/lxplus/StackHisto/ElePlots/");
    outFile += histname;
    outFile += "_ele"+dir+".pdf";
    c1->SaveAs(outFile);
    c1->Close();
  }
}

void example_stack_all(TString dir = ""){
  example_stack("cutflow", "cutflow", 1, true, true, true, false, true, 1, 15, dir, true);
  example_stack("pre_RelIso_ele", "RelIso of muons", 1, true, true, true, false, true, 0, 1, dir);
}

void example_stack_btag(TString dir="BTag"){//dir=KinFit, BTag
 example_stack(dir+"/pt_ele","Pt^{#e}[GeV]", 1, true,true,true,false, true,   0.0,    500.0, dir);
 example_stack(dir+"/eta_ele","#eta^{#e}", 1, true,true,true,false,true,       -3.0,   4.0,  dir);
 example_stack(dir+"/pt_jet","Pt^{jets}[GeV]", 1, true,true,true,false, true,  0.0,    500,  dir);
 example_stack(dir+"/eta_jet","#eta^{jets}", 1, true,true,true,false,true,     -3.0,   4.0,  dir);
 example_stack(dir+"/final_multi_jet","N^{jets}",1,true,true,true,false,true,   3,      15,  dir);
 example_stack(dir+"/final_pt_met","MET", 1, true,true,true,false, true,        0.0,    500, dir);
 //example_stack(dir+"/pfCISV","pfCombinedInclusiveSecondaryVertexV2BJetTags", 1, true,true,true,false,false,0.0,1.0);
 //example_stack(dir+"/pfCMVA","pfCombinedMVAV2BJetTags", 1, true,true,true,false,false,0.0,1.0);
 //example_stack(dir+"/pfCCvsL","pfCombinedCvsLJetTags", 1, true,true,true,false,false,0.0,1.0);
 //example_stack(dir+"/pfCCvsB","pfCombinedCvsBJetTags", 1, true,true,true,false,false,0.0,1.0);
 //example_stack(dir+"/mjj","m^{jj}",1,true,true,true,false, true,                0.0,    400.0);
 example_stack(dir+"/nvtx","N^{vertex}",1,true,true,true,false, true,           0.0,    70.0, dir);
 example_stack(dir+"/rhoAll","#rho",1,true,true,true,false, true,               0.0,    70.0, dir);
 example_stack(dir+"/wmt","MT[GeV]",1,true,true,true,false,true,                0,      300,  dir);
 //example_stack(dir+"/CSVL_count","N^{bjet}",1,true,true, true,false,true,       1,      10);
} 
  

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
//INPUT FILES
TFile* fData  = TFile::Open("all_muData.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fQCD	= TFile::Open("all_QCD.root");
TFile* fST	= TFile::Open("all_ST.root");
TFile* fTT	= TFile::Open("all_TTJetsP.root");
TFile *fSig   = TFile::Open("all_Hplus120.root");

//USER'S INPUT FOR DATA DRIVEN QCD 
bool isDataDrivenQCD = true;
double qcd_sf_btag = 2.618 ; 
double qcd_sf_kfit = 2.279 ;
double qcd_sf_btag_err = 0.20;
double qcd_sf_kfit_err = 0.34;
TFile *f_QCD_dd = new TFile("all_QCD_dd.root","RECREATE");
  
//SAVE HISTOS ON DISK
////bool isSaveHisto = true;
bool isSaveHisto = false;
////////////////////////////////////////////////////////

//--------------------------------------------//
//various functions
//--------------------------------------------//
//UP error in the unc band
double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus, TH1F* hQCD_dd, double qcd_sf_err);
//down error in the unc band
double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus, TH1F* hQCD_dd, double qcd_sf_err);
//unc graph
TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F* hQCD_dd, double qcd_sf_err, bool isFullGraph = false, bool isRatioGraph = false);
//function to stack histos
void stackHisto(TFile *filename, TString lable, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, THStack* MuptStack, TH1F* hMC, TLegend* leg);
//qcd from data
TH1F* getDataDrivenQCD(TString histname, double qcd_sf=1);
//function to add histograms
TH1F* addHistoForUnc(TString histname, TString sys, bool isDataDrivenQCD = false, double qcd_sf=1);
TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );
//get histogram from root file. Return empty hist, if the hist does not exit.
TH1F* getHisto(TFile *histFile, TString histPath, TString histName);

//--------------------------------------------//
//stack histos
//--------------------------------------------//
void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10, TString dir="BTag", bool label=false, double unc=false){
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
  cout<<"AAA"<<endl; 
  //Data
  TH1F* data = getHisto(fData, "base/Iso/", histname);
  //TH1F* data = (TH1F*)(fData->Get("base/Iso/"+histname))->Clone("data");
  //data->SetBinErrorOption(TH1::kPoisson);
  data->SetMarkerStyle(20); data->SetMarkerSize(0.8);
  if(axisrange) data->GetXaxis()->SetRangeUser(xmin,xmax);
  data->SetFillColor(kBlack);
  data->GetYaxis()->SetTitleOffset(1.35);
  data->GetYaxis()->SetTitle("Events");
  data->GetXaxis()->SetTitle(xaxis_title);
  ///data->Rebin(bin);
  data->SetTitle("");
  data->SetAxisRange(1.0, 1.0e6 ,"Y");
  if(drawdata)data->Draw("E"); 
  if(drawdata)leg->AddEntry(data,"Data","PE"); 
  
  //VV is the base histo
  TH1F* h1_base = getHisto(fVV, "base/Iso/", histname);
  h1_base->SetFillColor(13);
  if(axisrange) h1_base->GetXaxis()->SetRangeUser(xmin,xmax);
  leg->AddEntry(h1_base,"VV","F");
  //Define stacked histo
  THStack* MuptStack = new THStack("MuptStack","");
  MuptStack->Add(h1_base);
  //hMC = all Bkg MC samples
  TH1F* hMC = (TH1F*)h1_base->Clone("hMC");
  
  //---------------------------
  // QCD from Data
  //---------------------------
  TH1F * hQCD_dd;
  if(isDataDrivenQCD){
    if(dir=="BTag") hQCD_dd = getDataDrivenQCD(histname, qcd_sf_btag);
    if(dir=="KinFit") hQCD_dd = getDataDrivenQCD(histname, qcd_sf_kfit);
    if(axisrange)hQCD_dd->GetXaxis()->SetRangeUser(xmin,xmax);
    hQCD_dd->SetFillColor(kGreen);
    cout<<hQCD_dd->GetEntries()<<endl;
    cout<<hQCD_dd->Integral()<<endl;
    //hQCD->Draw();
    f_QCD_dd->cd();
    hQCD_dd->Write();
    leg->AddEntry(hQCD_dd,"QCD","F");
    MuptStack->Add(hQCD_dd);
    hMC->Add(hQCD_dd);
  }
  else stackHisto(fQCD, "QCD", histname, 3, 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fDY, "Z+jets", histname, 9, 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fST, "Single t", histname, 800 , 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fWJ, "W+ jets", histname, 6 , 1, axisrange, xmin, xmax, MuptStack, hMC, leg);
  stackHisto(fTT,"t#bar{t}", histname, 433, 1, axisrange, xmin, xmax, MuptStack, hMC, leg);

  //Signal 
  TH1F* hSig = (TH1F*)(fSig->Get("base/Iso/"+histname))->Clone("Signal");
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
  
  
  //-------------------------------------///
  //unc band
  //-------------------------------------///
  double qcd_sf = 1;
  double qcd_sf_err = 1;
  if(dir=="BTag"){ 
    qcd_sf = qcd_sf_btag;
    qcd_sf_err = qcd_sf_btag_err;
  }
  if(dir=="KinFit"){
    qcd_sf = qcd_sf_kfit;
    qcd_sf_err = qcd_sf_kfit_err;
  }
  if(unc){
  TGraphAsymmErrors *UncBand;
  UncBand = UNCGRAPH(addHistoForUnc(histname, "base", isDataDrivenQCD, qcd_sf),
      	    addHistoForUnc(histname, "JESPlus"),
      	    addHistoForUnc(histname, "JESMinus"),
      	    addHistoForUnc(histname, "JERPlus"),
      	    addHistoForUnc(histname, "JERMinus"),
      	    addHistoForUnc(histname, "bTagPlus"),
      	    addHistoForUnc(histname, "bTagMinus"),
	    hQCD_dd, qcd_sf_err, true, false);
  UncBand->SetFillColor(1);
  UncBand->SetFillStyle(3017);
  UncBand->Draw(" E2 same");
  leg->AddEntry(UncBand, "Unc","F"); 
  }
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
  //channel
  TPaveText *ch = paveText(0.953,0.9154898,0.6210067,0.9762187, 0, 19, 1, 0, 132);
  ch->AddText("#mu + jets");
  //CMS prili
  TPaveText *pt = paveText(0.09,0.9354,0.88,0.9362, 0, 19, 1, 0, 132);
  if(drawdata) pt->SetTextSize(0.059);
  else pt->SetTextSize(0.039);
  TText *text = pt->AddText(dir+":   CMS Preliminary,    #sqrt{s} = 13 TeV,    35.45 fb^{-1}");
  text->SetTextAlign(11);
  pt->Draw();
  if(drawsignal) cct->Draw();
  ch->Draw();
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
    hRatio->GetYaxis()->SetLabelSize(0.06);
  //lable x-axis, for cutflow
  if(label){
    TString steps[15] = {"muon trig","= 1 muon","0 electron","muon SF","RelIso", "#geq 4 jets","#slash{E}_{T} #geq 20 GeV", "MT >20 GeV", "#geq 2 b-jets","BTag SF", "fit converges","kfJetPt >25 GeV","dRJets <0.2",""};
    const size_t nsteps = sizeof(steps)/sizeof(TString);
    for(int istep=0; istep<nsteps; istep++ ){
      hRatio->GetXaxis()->SetBinLabel(istep+1, steps[istep]);
      hRatio->GetXaxis()->SetLabelSize(0.1);
      hRatio->GetXaxis()->LabelsOption("u");
    }
  }
    hRatio->Draw("E"); // use "P" or "AP"
    //unc band
    if(unc){
    TGraphAsymmErrors *UncBand_Ratio;
    UncBand_Ratio = UNCGRAPH(addHistoForUnc(histname, "base", isDataDrivenQCD, qcd_sf),
      	    addHistoForUnc(histname, "JESPlus"),
      	    addHistoForUnc(histname, "JESMinus"),
      	    addHistoForUnc(histname, "JERPlus"),
      	    addHistoForUnc(histname, "JERMinus"),
      	    addHistoForUnc(histname, "bTagPlus"),
      	    addHistoForUnc(histname, "bTagMinus"),
	    hQCD_dd, qcd_sf_err, false, true);
    UncBand_Ratio->SetFillColor(20);
    UncBand_Ratio->SetFillStyle(3001);
    UncBand_Ratio->Draw(" E2 same");
    }
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    TString outFile("$PWD/");
    outFile += histname;
    outFile += "_mu"+dir+".png";
    c1->SaveAs(outFile);
    c1->Close();
  }
}

double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus, TH1F *hQCD_dd, double qcd_sf_err){
  double errUp = sqrt(pow(fabs(hJESPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hJERPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(fabs(bTagPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+ pow(qcd_sf_err*hQCD_dd->GetBinContent(iBin+1),2));
  return errUp;
}

double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus, TH1F *hQCD_dd, double qcd_sf_err){
  double errDown =sqrt(pow(fabs(hCentral->GetBinContent(iBin+1) - hJESMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - hJERMinus->GetBinContent(iBin+1)),2) + 
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bTagMinus->GetBinContent(iBin+1)),2) + 
		  pow(hCentral->GetBinError(iBin+1),2)+pow(qcd_sf_err*hQCD_dd->GetBinContent(iBin+1),2));
  return errDown;
}

TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, TH1F* hQCD_dd, double qcd_sf_err, bool isFullGraph = false, bool isRatioGraph = false){
  TGraphAsymmErrors *gr;
  int n1 = hCentral->GetNbinsX(); 
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++){
    if(isFullGraph){
    Yval[i]   = hCentral->GetBinContent(i+1);
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus, hQCD_dd, qcd_sf_err); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus, hQCD_dd, qcd_sf_err); 
    }
    if(isRatioGraph){
    Yval[i]   = 1;
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus, hQCD_dd, qcd_sf_err); 
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus, hQCD_dd, qcd_sf_err); 
    errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
    errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
    }
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<errorD[i]<<endl;
    Xval[i]   = hCentral->GetBinCenter(i+1);
    XerrorU[i]= hCentral->GetBinWidth(i+1)/2;
    XerrorD[i]= hCentral->GetBinWidth(i+1)/2;
  }
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
} 

TH1F* getHisto(TFile *histFile, TString histPath, TString histName){
  TH1F* hist; 
  if(!(histFile->Get(histPath+histName))){
    hist = (TH1F*)(fTT->Get(histPath+histName))->Clone(histName);
    hist->Add(hist, -1);
  }else hist = (TH1F*)(histFile->Get(histPath+histName))->Clone(histName);
  return hist;
}

void stackHisto(TFile *filename, TString lable, TString histname, int color, double scale, bool axisrange, double xmin, double xmax, THStack* MuptStack, TH1F* hMC, TLegend* leg){
  TH1F* h2_base = getHisto(filename, "base/Iso/", histname);
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


TH1F* getDataDrivenQCD(TString histname, double qcd_sf=1){
  TH1F* hVV = getHisto(fVV, "base/NonIso/", histname); 
  TH1F* hDY = getHisto(fDY, "base/NonIso/", histname); 
  TH1F* hWJ = getHisto(fWJ, "base/NonIso/", histname); 
  TH1F* hST = getHisto(fST, "base/NonIso/", histname); 
  TH1F* hTT = getHisto(fTT, "base/NonIso/", histname); 
  TH1F* hData = getHisto(fData, "base/NonIso/", histname); 
  TH1F* hOtherMC = (TH1F*)hVV->Clone("hOtherMC"); 
  hOtherMC->Add(hDY); 
  hOtherMC->Add(hST); 
  hOtherMC->Add(hWJ); 
  hOtherMC->Add(hTT); 
  TH1F* hQCD = (TH1F*)hData->Clone(histname); 
  hQCD->Add(hOtherMC, -1);
  hQCD->Scale(qcd_sf);
  return hQCD;
}

TH1F* addHistoForUnc(TString histname, TString sys, bool isDataDrivenQCD = false, double qcd_sf=1){
  TH1F* hVV = getHisto(fVV, sys+"/Iso/", histname); 
  TH1F* hDY = getHisto(fDY, sys+"/Iso/", histname); 
  TH1F* hQCD_mc = getHisto(fQCD, sys+"/Iso/", histname); 
  TH1F* hWJ = getHisto(fWJ, sys+"/Iso/", histname); 
  TH1F* hST = getHisto(fST, sys+"/Iso/", histname); 
  TH1F* hTT = getHisto(fTT, sys+"/Iso/", histname); 
  TH1F* hAll = (TH1F*)hVV->Clone("hAllMC");
  hAll->Add(hDY);
  hAll->Add(hWJ);
  hAll->Add(hST);
  hAll->Add(hTT);
  TH1F* hQCD_dd = getDataDrivenQCD(histname, qcd_sf); 
  if(isDataDrivenQCD) hAll->Add(hQCD_dd);
  else hAll->Add(hQCD_mc);
  return hAll;
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
  example_stack("cutflow", "cutflow", 1, true, true, true, true, true, 0, 15, "", true, false);
  //example_stack("final_RelIso_mu", "RelIso of muons", 1, true, true, true, false, false, 0, 1, false, dir);
}

//void example_stack(TString histname, TString xaxis_title, int bin, bool log=false, bool drawdata=true, bool ratio=false, bool drawsignal=false, bool axisrange=false, double xmin=0, double xmax=10, TString dir="BTag", bool label=false, double unc=false){
//void example_stack_btag(TString dir="BTag"){//dir=KinFit, BTag
void example_stack_kfit(TString dir="KinFit"){//dir=KinFit, BTag
 //example_stack(dir+"/pt_mu","Pt^{#mu}[GeV]", 1, true,true,true,true, true,   0.0,    500.0, dir, false, true);
 //example_stack(dir+"/eta_mu","#eta^{#mu}", 1, true,true,true,true,true,       -2.5,   3.5,  dir, false, true);
 //example_stack(dir+"/pt_jet","Pt^{jets}[GeV]", 1, true,true,true,true, true,  0.0,    500,  dir, false, true);
 //example_stack(dir+"/eta_jet","#eta^{jets}", 1, true,true,true,true,true,     -2.5,   3.5,  dir, false, true);
 //example_stack(dir+"/final_multi_jet","N^{jets}",1,true,true,true,true,true,   3,      15,  dir, false, true);
 //example_stack(dir+"/final_pt_met","MET", 1, true,true,true,true, true,        0.0,    500, dir, false, true);
 if(dir=="BTag"){
   example_stack(dir+"/mjj","m^{jj}",1,true,true,true,true, true,                 0.0,    400.0, dir, false, true);
   example_stack(dir+"/CSVL_count","N^{bjet}",1,true,true, true,true,true,       1,      10, dir, false, true);
   example_stack(dir+"/pfCCvsL","pfCombinedCvsLJetTags",1,true,true, true,true,true,   -1.5, 2, dir, false, true);
   example_stack(dir+"/pfCCvsB","pfCombinedCvsBJetTags",1,true,true, true,true,true,   -1.5, 2, dir, false, true);
 }
 if(dir=="KinFit"){
   example_stack(dir+"/mjj_kfit","m^{jj}",1,true,true,true, true, true,      0.0,    400.0, dir, false, true);
   //example_stack(dir+"/mjj_kfit","m^{jj}",1,true,false,false, true, true,       0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_CTagL","mjj_kfit_CTagL",1,true,false,false, true, true,  0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_noCTagL","mjj_kfit_noCTagL",1,true,false,false, true, true, 0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_CTagM","mjj_kfit_CTagM",1,true,false,false, true, true,   0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_noCTagM","mjj_kfit_noCTagM",1,true,false,false, true, true,  0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_CTagT","mjj_kfit_CTagT",1,true,false,false, true, true,   0.0,    200.0, dir, false, false);
   //example_stack(dir+"/mjj_kfit_noCTagT","mjj_kfit_noCTagT",1,true,false,false, true, true, 0.0,    200.0, dir, false, false);
 }
 ///example_stack(dir+"/nvtx","N^{vertex}",1,true,true,true,true, true,           0.0,    70.0, dir, false, true);
 ///example_stack(dir+"/rhoAll","#rho",1,true,true,true,true, true,               0.0,    70.0, dir, false, true);
 ///example_stack(dir+"/wmt","MT[GeV]",1,true,true,true,true,true,                0,      165,  dir, false, true);
} 


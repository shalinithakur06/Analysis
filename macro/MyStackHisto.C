#include "MyStackHisto.h"

///////////////////////////////////////////  
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//SAVE HISTOS ON DISK
bool isSaveHisto = true;
///////////////////////////////////////////  

void plotStackedHisto(TString baseDir, TString isoDir, TString histDir, TString histName, TString xTitle,bool isData=false, bool isSig=false, double xmin=0, double xmax=10, double unc=false){
  MyStackHisto MyStackHisto_;
  string hist_name (histName);
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  TCanvas *c1 = new TCanvas(histName, histName, 700, 700);
  //TCanvas *c1 = new TCanvas("c1", "Data_MC", 400, 600);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.30,0.30,0.98};
  if(isData){
    c1->Divide(1, 2, 0, 0); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    if(isData) gPad->SetLogy(true);
    if(hist_name.find("mjj") != string::npos) gPad->SetLogy(false);
  }

  //-------------------------------
  //Legends
  //-------------------------------
  TLegend* leg = new TLegend(0.7518792,0.3061504,0.9512081,0.8798861,NULL,"brNDC");
  //TLegend* leg = new TLegend(0.7618792,0.3061504,0.9712081,0.8798861,NULL,"brNDC");
  if(hist_name.find("pt") != string::npos || hist_name.find("mt") != string::npos || hist_name.find("Fit") != string::npos ||hist_name.find("RelIso") != string::npos){
    leg = new TLegend(0.6018792,0.6061504,0.9912081,0.8898861,NULL,"brNDC");
    leg->SetNColumns(2);
  }
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.05);
  
  //-------------------------------
  // stack MC Bkg histo
  //-------------------------------
  THStack* hStack = new THStack("hStack","");
  TH1F* hVV = MyStackHisto_.getHisto(fVV, baseDir, isoDir, histDir, histName);
  TH1F* hMC = (TH1F*)hVV->Clone("hMC");
  int col_depth =0;
  hVV->SetFillColor(kOrange + col_depth);
  leg->AddEntry(hVV,"VV","F");
  hStack->Add(hVV);
  MyStackHisto_.stackHisto(fWJ, "W+ jets", baseDir, isoDir, histDir, histName, kViolet +col_depth , 1,   hStack, hMC, leg);
  MyStackHisto_.stackHisto(fTT,"t#bar{t} + jets", baseDir, isoDir, histDir, histName, kCyan+col_depth, 1,   hStack, hMC, leg);
  MyStackHisto_.stackHisto(fDY, "Z/#gamma+jets", baseDir, isoDir, histDir, histName, kSpring+col_depth, 1,   hStack, hMC, leg);
  
  // trim the histDir string
  std::string histDir_str;
  std::string histDir_orig(histDir);
  std::remove_copy(histDir_orig.begin(), histDir_orig.end(), std::back_inserter(histDir_str), '/');
  TString histDir_(histDir_str);

  gPad->SetTopMargin(0.10);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hStack->Draw("HIST");
  hStack->SetMinimum(1.0);
  hStack ->GetXaxis()->SetRangeUser(xmin, xmax);
  //cout<<hStack->GetMaximum()<<endl;
  if(isData){
    hStack->GetYaxis()->SetTitleOffset(0.70);
    hStack->GetYaxis()->SetTitleSize(0.10);   
    hStack->GetYaxis()->SetLabelSize(0.07);   
    hStack->GetYaxis()->SetTickLength(0.04); 
    hStack->GetYaxis()->SetTitle("Events");
    hStack->GetXaxis()->SetTitleOffset(1.20);
  }
  else{
  hStack->GetYaxis()->SetTitle("Events");
  hStack->GetXaxis()->SetTitle(xTitle);
  hStack->GetXaxis()->SetTitleSize(0.07);
  hStack->GetXaxis()->SetLabelSize(0.05);   
  hStack->GetXaxis()->SetTickLength(0.05); 
  //hStack->GetYaxis()->SetNdivisions(5);
  hStack->GetYaxis()->SetTitleSize(0.08);   
  hStack->GetYaxis()->SetLabelSize(0.05);   
  hStack->GetYaxis()->SetTickLength(0.04); 
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.15);
  hStack->GetYaxis()->SetTitleOffset(0.80);
  hStack->GetXaxis()->SetTitleOffset(0.90);
  }

  //-------------------------------------///
  //unc band
  //-------------------------------------///
  if(unc){
  TGraphAsymmErrors *UncBand;
  UncBand = MyStackHisto_.UNCGRAPH(
            MyStackHisto_.addHistoForUnc("base/", 	 isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JESPlus/",      isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JESMinus/",     isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JERPlus/",      isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JERMinus/",     isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("bTagPlus/",     isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("bTagMinus/",    isoDir, histDir, histName),
	    true, false);
  UncBand->SetFillColor(17);
  UncBand->SetFillStyle(3008);
  UncBand->Draw(" E2 same");
  leg->AddEntry(UncBand, "Unc","F"); 
  }

  //-------------------------------
  //Data
  //-------------------------------
  TH1F* hData = MyStackHisto_.getHisto(fData, baseDir, isoDir, histDir, histName);;
  ///MyStackHisto_.decorateHisto(hData, "", xTitle, "Events");
  hData->SetFillColor(kBlack);
  hData->SetMarkerStyle(20); hData->SetMarkerSize(1.2);
  if(isData)hData->Draw("Esame"); 
  if(isData)leg->AddEntry(hData,"Data","PE"); 

  //-------------------------------
  //Signal 
  //-------------------------------
  TH1F* hSig = MyStackHisto_.getHisto(fSig, baseDir, isoDir, histDir, histName);
  if(histDir_str.find("PtbJet") != string::npos)
    hSig = MyStackHisto_.getHisto(fSig, baseDir, isoDir, histDir, histName);
  ///MyStackHisto_.decorateHisto(hSig, "", xTitle, "Events");
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  if(isSig)hSig->Draw("HISTSAME"); 
  leg->AddEntry(hSig, "Signal","L"); 
  //leg->AddEntry(hSig, "#splitline{Signal}{M_{H^{+}} = 120 GeV}","L"); 
  leg->Draw();
  double yMax = 0;
  if(hData->GetMaximum() > hSig->GetMaximum()) yMax = hData->GetMaximum();
  else yMax = hSig->GetMaximum();
  if(yMax < hMC->GetMaximum()) yMax = hMC->GetMaximum();

  if(isData) hStack->SetMaximum(4.0*hStack->GetMaximum());
  else hStack->SetMaximum(1.1*yMax);
  TPaveText *cct = MyStackHisto_.paveText(0.35,0.8454,0.45,0.8462, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.07);
  cct->AddText("M_{l^{*}} = 250 GeV");

  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //hist name
  TPaveText *hLable = MyStackHisto_.paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.07);
  hLable->AddText(xTitle);
  
  //channel
  TPaveText *ch = MyStackHisto_.paveText(0.823,0.9154898,0.9210067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.10);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = MyStackHisto_.paveText(0.01,0.9554,0.82,0.9562, 0, 19, 1, 0, 132);
  if(isData) pt->SetTextSize(0.080);
  else pt->SetTextSize(0.05);
  pt->AddText(histDir_+": 35.9 fb^{-1} (13 TeV)");
  //TText *text = pt->AddText(histDir+": CMS Preliminary, #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  pt->Draw();
  if(isSig) cct->Draw();
  ch->Draw();
  //hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(isData){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.5); //gPad->SetGridy();
    if(histDir=="") gPad->SetBottomMargin(0.55);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(hData);
    hRatio->Divide(hMC); 
    MyStackHisto_.decorateHisto(hRatio, "", xTitle, "#frac{Data}{Bkg}");
    hRatio->SetFillColor(kBlack);
    if(baseDir=="baseLowMET/") hRatio->GetYaxis()->SetRangeUser(0, 2);
    else hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetYaxis()->SetTitleOffset(0.40);
    hRatio->GetXaxis()->SetTitleOffset(0.90);
    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.2);
    hRatio->GetYaxis()->SetTitleSize(0.15); 
    hRatio->GetXaxis()->SetTitleSize(0.23);
    hRatio->GetXaxis()->SetLabelSize(0.15); 
    hRatio->GetYaxis()->SetLabelSize(0.10); 
    if(hist_name.find("mjj") != string::npos){
      hRatio->GetXaxis()->SetTitleSize(0.15); 
      hRatio->GetXaxis()->SetTitleOffset(1.40);
    }
    //lable x-axis, for cutflow
    if(histName=="cutflow"){
      vector<string >cut_label;
      if(isEleChannel){
        cut_label.push_back("Ele trigger");
      }
      if(isMuChannel){
        cut_label.push_back("Mu trigger");
      }
      cut_label.push_back("Control Sel");
      cut_label.push_back("b-jet veto");
      cut_label.push_back("PreSel");
      cut_label.push_back("ZTag Sel");
      for(int istep=0; istep<cut_label.size(); istep++ ){
       hRatio->GetXaxis()->SetBinLabel(istep+1, cut_label[istep].c_str());
      }
      hRatio->GetXaxis()->LabelsOption("v");
      hRatio->GetXaxis()->SetTickLength(0.08); 
      hRatio->GetXaxis()->SetLabelOffset(0.03);
      hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    }
    //unc band
    hRatio->Draw("E"); // use "P" or "AP"
    if(unc){
    TGraphAsymmErrors *UncBand_Ratio;
    UncBand_Ratio = MyStackHisto_.UNCGRAPH(
	    MyStackHisto_.addHistoForUnc("base/", 	isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JESPlus/",     isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JESMinus/",    isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JERPlus/",     isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("JERMinus/",    isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("bTagPlus/",    isoDir, histDir, histName),
      	    MyStackHisto_.addHistoForUnc("bTagMinus/",   isoDir, histDir, histName),
	    false, true);
    UncBand_Ratio->SetFillColor(17);
    //UncBand_Ratio->SetFillStyle(19);
    UncBand_Ratio->Draw("E2 same");
    }
    hRatio->Draw("E same"); // use "P" or "AP"
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    mkdir(histDir_, S_IRWXU);
    TString outFile("$PWD/");
    outFile += histDir_+"/"+histName;
    if(isMuChannel) outFile += "_mu"+histDir_+".pdf";
    if(isEleChannel) outFile += "_ele"+histDir_+".pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}

void MyStackHisto(){

  TString histDir="ControlP/"; // ControlP/, PreSel/, ZTag/;
  TString baseDir = "base/"; // base/;

  TString isoDir = "Iso/";
  bool isDataMjj=false;
  //flags
  bool isData = true;
  bool isSig = true;
  bool isUnc = true;
  plotStackedHisto(baseDir, isoDir, "", "cutflow","cutflow", isData,  isSig,  0.5, 7.5, false);
  plotStackedHisto(baseDir, isoDir, histDir, "pt_jet", "Pt^{jets} [GeV]", isData, isSig,  100, 1000.0, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "eta_jet", "#eta^{jets}", isData, isSig,  -3.0, 5.0, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "final_multi_jet", "N^{jets}", isData, isSig,  0, 7 , isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "nvtx", "N^{vertex}", isData, isSig,  0, 70, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "rhoAll", "#rho", isData, isSig,  0, 70, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "ak8Pmass", "Pruned mass of jets", isData, isSig,  0.0, 500, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "ak8Tau21", "Jet-subjetiness (#tau_{21})", isData, isSig,  0, 1.5, isUnc);
  if(isMuChannel){
  plotStackedHisto(baseDir, isoDir, histDir, "dR1", "#Delta R (jets, #mu_{1})", isData, isSig,  0, 8, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "dR2", "#Delta R (jets, #mu_{2})", isData, isSig,  0, 8, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_1stMu", "Pt^{#mu_{1}} [GeV]", isData, isSig,  0.0, 1000.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_1stMu", "#eta^{#mu_{1}}", isData, isSig,  -3.0, 5.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_2ndMu", "Pt^{#mu_{2}} [GeV]", isData, isSig,  0.0, 500.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_2ndMu", "#eta^{#mu_{2}}", isData, isSig,  -3.0, 5.0, isUnc);
    //plotStackedHisto(baseDir, isoDir, histDir, "multi_mu", "N^{#mu}", isData, isSig,  1, 5, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "mll", "M^{#mu_{1}#mu_{2}} [GeV]", isData, isSig,  50, 1000, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_Z", "Pt^{#mu_{1}#mu_{2}}[GeV]", isData, isSig,  0, 1000, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_Z", "#eta^{#mu_{1}#mu_{2}}", isData, isSig,  -3.0, 4.5, isUnc);
    if(histDir=="ZTag/"){
      plotStackedHisto(baseDir, isoDir, histDir, "mlZ_max","M^{#mu Z}_{max}", isDataMjj, isSig,  200, 3000, isUnc);
      plotStackedHisto(baseDir, isoDir, histDir, "mlZ_min","M^{#mu Z}_{min}", isDataMjj, isSig,  50, 1000, isUnc);
    }
  }
  if(isEleChannel){
  plotStackedHisto(baseDir, isoDir, histDir, "dR1", "#Delta R (jets, e_{1})", isData, isSig,  0, 8, isUnc);
  plotStackedHisto(baseDir, isoDir, histDir, "dR2", "#Delta R (jets, e_{2})", isData, isSig, 0, 8, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_1stEle", "Pt^{e_{1}} [GeV]", isData, isSig,  0.0, 1000.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_1stEle", "#eta^{e_{1}}", isData, isSig,  -3.0, 5.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_2ndEle", "Pt^{e_{2}} [GeV]", isData, isSig,  0.0, 500.0, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_2ndEle", "#eta^{e_{2}}", isData, isSig,  -3.0, 5.0, isUnc);
    //plotStackedHisto(baseDir, isoDir, histDir, "multi_Ele", "N^{e}", isData, isSig,  1, 5, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "mll", "M^{e_{1}e_{2}} [GeV]", isData, isSig,  50, 1000, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "pt_Z", "Pt^{e_{1}e_{2}}[GeV]", isData, isSig,  0, 1000, isUnc);
    plotStackedHisto(baseDir, isoDir, histDir, "eta_Z", "#eta^{e_{1}e_{2}}", isData, isSig,  -3.0, 4.5, isUnc);
    if(histDir=="ZTag/"){
      plotStackedHisto(baseDir, isoDir, histDir, "mlZ_max","M^{eZ}_{max}", isDataMjj, isSig,  200, 3000, isUnc);
      plotStackedHisto(baseDir, isoDir, histDir, "mlZ_min","M^{eZ}_{min}", isDataMjj, isSig,  50, 1000, isUnc);
    }
  }

  plotStackedHisto(baseDir, isoDir, "ControlP/", "pfCISV","b-discriminator", isData, isSig,  0.0, 1.5, isUnc);
  plotStackedHisto(baseDir, isoDir, "BTag/", "multi_bjet","N^{b-jet}", isData, isSig,  0.5, 3.5, isUnc);
}

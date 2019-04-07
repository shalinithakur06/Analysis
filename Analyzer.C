
///////////////////////
// Electron Channel
///////////////////////

#include "Analyzer.h"
#include <map>

using namespace std;
void Analyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  //check if the input file is MC or Data  
  Reader *evR_;  
  evR_ = new Reader();
  TFile *f_ = TFile::Open(url);
  int nEntries = evR_->AssignEventTreeFrom(f_);
  MyEvent *ev_;
  ev_ = evR_->GetNewEvent(1);

  CutFlowProcessor(url, myKey, "base", outFile_);
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
 if(!ev_->isData){ 
    CutFlowProcessor(url, myKey, "JESPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JESMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "bTagMinus", 	outFile_);
  }
  outFile_->Write(); 
  outFile_->Close();
  f_->Close();
  delete f_;
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void Analyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, TFile *outFile_){
  int input_count_PreSel = 0;
  int input_count_ZTag = 0;
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  int jes = 0, jer = 0, metuc = 0, bScale = 0;

  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("bTagPlus"))bScale = 1;
  else if (cutflowType.Contains("bTagMinus"))bScale = -1; 
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35860;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  fillHisto(outFile_, cutflowType, "", "totalEvents", 10, 0, 10000000000, initialEvents, 1 );
  MyEvent *ev;
  
  //---------------------------------------------------//
  //BTag SF: read CSV file for SF, 2D histos for eff 
  //---------------------------------------------------//      
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  const std::string & bTagCSVfile 	= "stack/CSVv2_Moriond17_B_H.csv";
  const std::string & bTagName 		= "CSVv2";
  const std::string & bTagSys 		= "central"; 
  if(bScale==1) const std::string &bTagSys 		= "up"; 
  if(bScale==-1)const std::string &bTagSys 		= "down"; 
  const std::vector<std::string> & otherSysTypes = {"up", "down"};
  //b-quark
  BTagCalibrationReader readBTagCSV_bT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_B);
  //c-quark
  BTagCalibrationReader readBTagCSV_cT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "comb", bTagSys, otherSysTypes, BTagEntry::FLAV_C);
  //other(light) quarks and gluon
  BTagCalibrationReader readBTagCSV_lT= readCSV(bTagCSVfile, bTagName, BTagEntry::OP_TIGHT,
    	      "incl", bTagSys, otherSysTypes, BTagEntry::FLAV_UDSG);
  
  //getBTagEffHistos(f);
  TString histPath("myMiniTreeProducer/Jets/");
  TH2D* h2_BTagEff_Denom_b 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_b"));
  TH2D* h2_BTagEff_Denom_c 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_c"));
  TH2D* h2_BTagEff_Denom_udsg 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Denom_udsg")); 
  TH2D* h2_BTagEff_Num_bT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_bT"));
  TH2D* h2_BTagEff_Num_cT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_cT"));
  TH2D* h2_BTagEff_Num_udsgT 		= (TH2D*)(f->Get(histPath+"h2_BTagEff_Num_udsgT")); 
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  double kfCount = 0;
  double n_negEvt = 0.0;
  double n_posEvt = 0.0;
  double n_noCharge = 0.0;
  double n_oppCharge = 0.0;
  double n_sameCharge = 0.0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%10000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
    //if(i > 20000) break;
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      TString sampleName_(sampleName);
      if(sampleName_.Contains("DYJetsToLL")){
        double sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight;
        evtWeight *= ev->sampleInfo.gen_weight;
        if(ev->sampleInfo.gen_weight ==1) n_posEvt++;
        else n_negEvt++;
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      }	  
      //lumi weight
      else {
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight; 
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      } 
      //pileup weight
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
        float npu = pu[0];
        double weightPU = LumiWeights_.weight(npu);
        evtWeight *= weightPU;  
        fillHisto(outFile_, cutflowType, "", "puWeight", 1000, 0, 100, weightPU, 1 );
	//////cout<<"evtWeight="<< weightPU<<endl;
      }
      if(i==0){
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
	fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents*sampleWeight);
      }
    } 
    else{ 
      if(i==0)fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents);
    }
    
    //---------------------------------------------------//
    //apply electron triggers
    //---------------------------------------------------//
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_DoubleEle33_CaloIdL") != string::npos) {
        passTrig = true;
      }
    }
    if(!passTrig) continue; 
    double nCutPass = 1.0;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
   
    //---------------------------------------------------//
    //get all objets e.g. leptons, jets, vertices etc.
    //---------------------------------------------------//
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    //preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    //clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    int nEle = e_init.size();
    double pri_vtxs = Vertices[0].totVtx;
    if(nEle < 2)continue;
    int e1 = e_init[0];
    int e2 = e_init[1];
    //fill histo for relative isolation of both electrons
    double e1RelIso = pfElectrons[e1].relCombPFIsoEA;
    double e2RelIso = pfElectrons[e2].relCombPFIsoEA;
    fillHisto(outFile_, cutflowType, "", "e1RelIso", 100, 0, 1, e1RelIso, evtWeight);
    fillHisto(outFile_, cutflowType, "", "e2RelIso", 100, 0, 1, e2RelIso, evtWeight);
    
    //get Medium ID for first electron
    //veto first or 2nd electron, if they are fake
    ///if(looseMuonVeto( e1, pfElectrons, isPFlow) ) continue;
    ///if(looseMuonVeto( e2, pfElectrons, isPFlow) ) continue;
    //both electrons should have opposite charge
    //charge selection
    n_noCharge++ ;
    if(pfElectrons[e1].charge == pfElectrons[e2].charge) n_sameCharge++;
    if(pfElectrons[e1].charge != pfElectrons[e2].charge) n_oppCharge++;
    //both muons should have opposite charge
    //if(pfElectrons[e1].charge == pfElectrons[e2].charge) continue; // to be consistent with other AN
     
    //events should not have any electron
    ///if(looseElectronVeto(-1, pfElectrons, Vertices[0], isPFlow)) continue;
    int count_ele = e_init.size();
    
    //---------------------------------------------------//
    //apply Electron SF to eventWeights 
    //---------------------------------------------------//
    //Reco, ID, trigger, HEEP
    double eleSF1 =0;
    double ele_recoSF1       = getEleSF(h2_ele_recoSF, pfElectrons[e1].eleSCEta, pfElectrons[e1].p4.pt());
    //This is cut-based ID, we are using Heep ID
    //double ele_medium_idSF1  = getEleSF(h2_ele_medium_idSF, pfElectrons[e1].eleSCEta, pfElectrons[e1].p4.pt());
    double ele_trigSF1       = getEleTrigSF(h2_ele_trigSF, pfElectrons[e1].eleSCEta, pfElectrons[e1].p4.pt());
    double ele_heep_SF1      = getEleHeep2SF(tg_heep_SF, pfElectrons[e1].eleSCEta);
    eleSF1 = ele_recoSF1*ele_trigSF1*ele_heep_SF1;  

    double eleSF2 =0;
    double ele_recoSF2       = getEleSF(h2_ele_recoSF, pfElectrons[e2].eleSCEta, pfElectrons[e2].p4.pt());
    //double ele_medium_idSF2  = getEleSF(h2_ele_medium_idSF, pfElectrons[e2].eleSCEta, pfElectrons[e2].p4.pt());
    double ele_trigSF2       = getEleTrigSF(h2_ele_trigSF, pfElectrons[e2].eleSCEta,  pfElectrons[e2].p4.pt());
    double ele_heep_SF2      = getEleHeep2SF(tg_heep_SF, pfElectrons[e2].eleSCEta);
    eleSF2 = ele_recoSF2*ele_trigSF2*ele_heep_SF2;

    //Scale factors are applied on MC only.
    double eleSF = 1.0;
    if(!ev->isData) eleSF = eleSF1*eleSF2;
    evtWeight *= eleSF;
    fillHisto(outFile_, cutflowType, "", "eleSF", 1000, 0, 100, eleSF, 1 );
    //---------------------------------------------------//
    // Iso(<0.15) and Non-iso(>0.15) region 
    //---------------------------------------------------//
    bool isofound = false;
    string cutflowType_(cutflowType);
    cutflowType_ = cutflowType+"/Iso";
    if(e1RelIso > 0.08) continue;
    if(e2RelIso > 0.08) continue;
    
    //---------------------------------------------------//
    //get 4 vector for Z boson
    //---------------------------------------------------//
    MyLorentzVector vZ = pfElectrons[e1].p4 + pfElectrons[e2].p4;
    if(vZ.mass() < 60) continue;    
    double dR1 = 0.0;
    double dR2 = 0.0;
    int count_jets = j_final.size();
    //---------------------------------------------------//
    //Fill histos with for Control Plots
    //---------------------------------------------------//
    //fill histos for muon
    double elePt1 = pfElectrons[e1].p4.pt();
    double elePt2 = pfElectrons[e2].p4.pt(); 
    if(j_final.size()==0) continue;
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e2].p4);
      fillHisto(outFile_, cutflowType_, "ControlP","dR1", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","dR2", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","dR", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","dR", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","pt_jet", 500, 0, 10000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","ak8Pmass", 500, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "ControlP","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","multi_Ele",  15, 0.5, 15.5, count_ele, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_1stEle", 500, 0, 10000, elePt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_2ndEle", 500, 0, 10000, elePt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_1stEle", 50, -5, 5, pfElectrons[e1].p4.eta(), evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ControlP", "ptEle1_ptEle2", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "ControlP","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ControlP","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","mll", 500, 0, 10000, vZ.M(), evtWeight );
    
    
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "ControlP","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "ControlP","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );

    //---------------------------------------------------//
    //apply B-tagging
    //---------------------------------------------------//
    double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
    double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
    int count_CSVT_SF = 0; 
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
      pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
      fillHisto(outFile_, cutflowType_, "ControlP", "pfCISV", 100, -2, 2, pfCISV, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP", "pfCMVA", 100, -2, 2, pfCMVA, evtWeight );
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
      if(pfCISV > 0.9535){
        count_CSVT_SF++; 
        double jetPt = jetPtWithJESJER(pfJets[ijet], jes, jer);
        fillHisto(outFile_, cutflowType_, "BTag", "pt_bjet", 100, 0, 1000, jetPt, evtWeight );
        fillHisto(outFile_, cutflowType_, "BTag", "eta_bjet", 50, -5, 5, pfJets[ijet].p4.eta(), evtWeight );
      }
    }
    fillHisto(outFile_, cutflowType_, "BTag","multi_bjet",  15, 0.5, 15.5, count_CSVT_SF, evtWeight );
    if(count_CSVT_SF !=0) continue; // reject events if there is any b-jet

    //Apply b-tag SF
    //https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a)%20Event%20reweighting%20using%20scal
    double pmc_btag = 1.0;
    double pdata_btag = 1.0;
    double bTagWt = 1.0; 
    if(!ev->isData){
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
        int ind_jet = j_final[ijet];
        double pMC_ = 1.0;
        double pData_ = 1.0;
        //b-quark
        if(abs(pfJets[ind_jet].partonFlavour) ==5){
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_bT, h2_BTagEff_Denom_b, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_bT, h2_BTagEff_Num_bT, h2_BTagEff_Denom_b, pfJets[ind_jet],bScale);
        }
        //c-quark
        else if(abs(pfJets[ind_jet].partonFlavour) ==4){ 
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_cT, h2_BTagEff_Denom_c, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_cT, h2_BTagEff_Num_cT, h2_BTagEff_Denom_c, pfJets[ind_jet],bScale);
        }
        //other quarks and gluon
        else{ 
          pMC_ = getBTagPmcSys(h2_BTagEff_Num_udsgT, h2_BTagEff_Denom_udsg, pfJets[ind_jet]); 
          pData_ = getBTagPdataSys(readBTagCSV_lT, h2_BTagEff_Num_udsgT, h2_BTagEff_Denom_udsg, pfJets[ind_jet], bScale); 
        }
        pmc_btag = pmc_btag*pMC_;
        pdata_btag = pdata_btag*pData_;
      }
    bTagWt = pdata_btag/pmc_btag;
    }
    evtWeight *= bTagWt;
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight);
    fillHisto(outFile_, cutflowType, "", "bTagWeight", 100, 0, 2, bTagWt, 1);

    //---------------------------------------------------//
    //Fill histos with pre-selection
    //---------------------------------------------------//
    bool isPreSel = false;
    if(vZ.Pt() < 100) continue;    
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      if(jetPt <= 100) continue;    
      if(fabs(pfJets[ind_jet].p4.eta()) >= 2.5) continue;    
      if(pfJets[ind_jet].ak8Pmass <= 40) continue;    
      dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e2].p4);
      if(dR1 < 0.8 || dR2 < 0.8) continue;    
      isPreSel = true;
      fillHisto(outFile_, cutflowType_, "PreSel","dR1", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR2", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","pt_jet", 500, 0, 10000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","ak8Pmass", 500, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
    }
    if(isPreSel){
    fillHisto(outFile_, cutflowType_, "PreSel","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///if(muonPt1 <100) continue;    
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "PreSel","multi_Ele",  15, 0.5, 15.5, count_ele, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_1stEle", 500, 0, 5000, elePt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_2ndEle", 500, 0, 10000, elePt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_1stEle", 50, -5, 5, pfElectrons[e1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto2D(outFile_, cutflowType_,"PreSel", "ptEle1_ptEle", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "PreSel","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "PreSel","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","mll", 500, 0, 10000, vZ.M(), evtWeight );
    
    
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "PreSel","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "PreSel","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    input_count_PreSel++;
    if(input_count_PreSel%100==0)
    cout << "input count after PreSel: "<< input_count_PreSel << endl;
    }

    //---------------------------------------------------//
    // fill histo after ZTag:
    //---------------------------------------------------//
    // tag a jet as Z-boson
    vector<size_t> allZjet;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double ak8Pmass_ = pfJets[ind_jet].ak8Pmass;
      double ak8Tau21 = pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1;
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      if(jetPt <= 100) continue;    
      if(fabs(pfJets[ind_jet].p4.eta()) >= 2.5) continue;    
      if(pfJets[ind_jet].ak8Pmass <= 40) continue;    
      dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e2].p4);
      if(dR1 < 0.8 || dR2 < 0.8) continue;    
      if(vZ.M()> 200 && jetPt >200 && ak8Pmass_ > 70 && ak8Pmass_ < 110){ 
	  fillHisto(outFile_, cutflowType_, "ZTag", "totalEvt", 20, -1, 1, 0, evtWeight);
	if(ak8Tau21<0.45) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.45, evtWeight);
	if(ak8Tau21<0.46) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.46, evtWeight);
	if(ak8Tau21<0.47) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.47, evtWeight);
	if(ak8Tau21<0.48) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.48, evtWeight);
	if(ak8Tau21<0.49) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.49, evtWeight);
	if(ak8Tau21<0.50) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.50, evtWeight);
	if(ak8Tau21<0.51) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.51, evtWeight);
	if(ak8Tau21<0.52) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.52, evtWeight);
	if(ak8Tau21<0.53) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.53, evtWeight);
	if(ak8Tau21<0.54) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.54, evtWeight);
	if(ak8Tau21<0.55) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.55, evtWeight);
	if(ak8Tau21<0.56) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.56, evtWeight);
	if(ak8Tau21<0.57) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.57, evtWeight);
	if(ak8Tau21<0.58) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.58, evtWeight);
	if(ak8Tau21<0.59) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.59, evtWeight);
	if(ak8Tau21<0.60) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.60, evtWeight);
	if(ak8Tau21<0.61) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.61, evtWeight);
	if(ak8Tau21<0.62) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.62, evtWeight);
	if(ak8Tau21<0.63) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.63, evtWeight);
	if(ak8Tau21<0.64) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.64, evtWeight);
	if(ak8Tau21<0.65) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.65, evtWeight);
	if(ak8Tau21<0.66) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.66, evtWeight);
	if(ak8Tau21<0.67) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.67, evtWeight);
	if(ak8Tau21<0.68) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.68, evtWeight);
	if(ak8Tau21<0.69) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.69, evtWeight);
	if(ak8Tau21<0.70) 
	  fillHisto(outFile_, cutflowType_, "ZTag", "passedEvt", 26, 0.445, 0.705, 0.70, evtWeight);
	if(ak8Tau21 < 0.60) allZjet.push_back(ijet);
      }
    }

    //apply tau21 scale factor
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging#2016_scale_factors_and_correctio
    if(!ev->isData) evtWeight *= 1.11; 

    if(allZjet.size()==0) continue;
    double tau21 = pfJets[j_final[allZjet[0]]].ak8Tau2/pfJets[j_final[allZjet[0]]].ak8Tau1;
    fillHisto2D(outFile_, cutflowType_,"", "tau21_mll", 500, 0, 1, tau21, 50, 0, 1000,  vZ.M(), 1);
    /*////////////////////////
    if(vZ.M() > 200 && tau21 < 0.5) cutflowType_ = cutflowType_+"/RegionA";
    if(vZ.M() < 200 && tau21 < 0.5) cutflowType_ = cutflowType_+"/RegionB";
    if(vZ.M() < 200 && tau21 > 0.5) cutflowType_ = cutflowType_+"/RegionC";
    if(vZ.M() < 200 && tau21 > 0.5) cutflowType_ = cutflowType_+"/RegionC";
    */////////////////////
    MyLorentzVector vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[e1].p4;
    MyLorentzVector vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[e2].p4;
    
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "ZTag","multi_Ele",  15, 0.5, 15.5, count_ele, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_1stEle", 500, 0, 10000, elePt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_2ndEle", 500, 0, 10000, elePt2, evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ZTag", "ptEle1_ptEle", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "ZTag","eta_1stEle", 50, -5, 5, pfElectrons[e1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ZTag","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_Z",  50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","phi_Z",  50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","mll", 500, 0, 10000, vZ.M(), evtWeight );

    //Section 5.1 of https://drive.google.com/drive/folders/1e8PiEjw7sWXJn7-ET8TaK81QOgFHV0_g
    double mlZmax = 0.0;
    double mlZmin = 0.0;
    if(vZmax.M() >= vZmin.M()){
      mlZmax = vZmax.M();
      mlZmin = vZmin.M();
    }
    else{
      mlZmax = vZmin.M();
      mlZmin = vZmax.M();
    }
    fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min", 500, 0, 10000, mlZmin, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max", 500, 0, 10000, mlZmax, evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ZTag", "mlZmin_mlZmax",500, 0, 10000, mlZmin, 500, 0, 10000, mlZmax, 1);

    //L-cut
    vector<string> sigMass;
    vector<double> lCutMin;
    vector<double> lCutMax;
    sigMass.push_back("250");   lCutMax.push_back(440);  lCutMin.push_back(300);
    sigMass.push_back("500");   lCutMax.push_back(450);  lCutMin.push_back(560);
    sigMass.push_back("750");   lCutMax.push_back(700);  lCutMin.push_back(900);
    sigMass.push_back("1000");  lCutMax.push_back(950);  lCutMin.push_back(1080);
    sigMass.push_back("1250");  lCutMax.push_back(1200); lCutMin.push_back(1370);
    sigMass.push_back("1500");  lCutMax.push_back(1300); lCutMin.push_back(1700);
    sigMass.push_back("1750");  lCutMax.push_back(1300); lCutMin.push_back(1950);
    sigMass.push_back("2000");  lCutMax.push_back(1300); lCutMin.push_back(2200);
    sigMass.push_back("2500");  lCutMax.push_back(1300); lCutMin.push_back(2700);
    sigMass.push_back("3000");  lCutMax.push_back(1300); lCutMin.push_back(3200);
    sigMass.push_back("3500");  lCutMax.push_back(1300); lCutMin.push_back(3700);
    sigMass.push_back("4000");  lCutMax.push_back(1300); lCutMin.push_back(4200);
    sigMass.push_back("4500");  lCutMax.push_back(1300); lCutMin.push_back(4700);
    sigMass.push_back("5000");  lCutMax.push_back(1300); lCutMin.push_back(5200);
    for(int l =0; l<sigMass.size(); l++){
      double max = lCutMax[l];
      double min = lCutMin[l];
      TString mass = sigMass[l];
      if(mlZmax > max && mlZmin < min){
      fillHisto(outFile_, cutflowType_, "ZTag1","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag1","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max + 50) && mlZmin < (min - 50) ){
      fillHisto(outFile_, cutflowType_, "ZTag2","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag2","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max + 100) && mlZmin < (min - 100) ){
      fillHisto(outFile_, cutflowType_, "ZTag3","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag3","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max - 50) && mlZmin < (min + 50) ){
      fillHisto(outFile_, cutflowType_, "ZTag4","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag4","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max - 100) && mlZmin < (min + 100) ){
      fillHisto(outFile_, cutflowType_, "ZTag5","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag5","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max) && mlZmin < (min - 50) ){
      fillHisto(outFile_, cutflowType_, "ZTag6","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag6","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max) && mlZmin < (min - 100) ){
      fillHisto(outFile_, cutflowType_, "ZTag7","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag7","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max + 50) && mlZmin < (min) ){
      fillHisto(outFile_, cutflowType_, "ZTag8","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag8","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
      if(mlZmax > (max + 100) && mlZmin < (min) ){
      fillHisto(outFile_, cutflowType_, "ZTag9","mlZ_min_sig"+mass, 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag9","mlZ_max_sig"+mass, 500, 0, 10000, mlZmax, evtWeight );
      }
    }
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e2].p4);
      fillHisto(outFile_, cutflowType_, "ZTag","dR1", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","dR2", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","dR", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","dR", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","pt_jet", 500, 0, 10000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","ak8Pmass", 500, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );

    }
    fillHisto(outFile_, cutflowType_, "ZTag","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "ZTag","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "ZTag","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    input_count_ZTag++;
    if(input_count_ZTag%10==0)
    cout << "input count after ZTag: "<< input_count_ZTag << endl;
    //if(i > 200000) break;*/
  }//event loop
  cout<<"Total events  = "<<nEntries<<endl;
  cout<<"Total events with negative weight = "<<n_negEvt<<endl;
  cout<<"Total events with positive weight = "<<n_posEvt<<endl;
  double effective_evt = (n_posEvt-n_negEvt);
  double amcnlo_weight = 1.0;
  if(effective_evt !=0) amcnlo_weight = effective_evt/nEntries;
  fillHisto(outFile_, cutflowType, "", "amcnlo_weight", 10, 0, 1, amcnlo_weight,1);
  fillHisto(outFile_, cutflowType, "", "noCharge", 10, -2, 2, 0, n_noCharge);
  fillHisto(outFile_, cutflowType, "", "oppCharge", 10, -2, 2, -1, n_oppCharge);
  fillHisto(outFile_, cutflowType, "", "sameCharge", 10, -2, 2, 1, n_sameCharge);
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("DYJetsToLL_M50_EleMC_20190117_Ntuple_8.root", "PF", ""); 
  ///CutFlowAnalysis("TT_EleMC_20190117_Ntuple_7.root", "PF", "");
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/sthakur/ntuple_for2016Data_EleMC_20190117/EleMC_20190117/TT_EleMC_20190117/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TT_EleMC_20190117/190117_103502/0000/TT_EleMC_20190117_Ntuple_1.root", "PF", "");

  //====================================
  //condor submission
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 

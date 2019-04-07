
///////////////////////
// Muon Channel
///////////////////////

#include "Analyzer.h"
#include <map>

using namespace std;
void Analyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel(9);
  
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
  double n_oneTrig = 0.0;
  double n_twoTrig = 0.0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
    //if(i > 1000) break; 
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
    //apply muon triggers
    //---------------------------------------------------//
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6
    bool passTrig = false;
    bool passOneTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_IsoMu24") != string::npos || 
		      trig[it].find("HLT_IsoTkMu") != string::npos) {
        passTrig = true;
      }
    } 
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_IsoMu24") != string::npos) passOneTrig = true;
    } 
    if(passTrig)n_twoTrig++;
    if(passOneTrig) n_oneTrig++;
    
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
    int nMuon = m_init.size();
    double pri_vtxs = Vertices[0].totVtx;
    if(nMuon < 2)continue;
    int m1 = m_init[0];
    int m2 = m_init[1];
    double m1RelIso = pfMuons[m1].pfRelIso;
    double m2RelIso = pfMuons[m2].pfRelIso;
    fillHisto(outFile_, cutflowType, "", "m1RelIso", 100, 0, 1, m1RelIso, evtWeight);
    fillHisto(outFile_, cutflowType, "", "m2RelIso", 100, 0, 1, m2RelIso, evtWeight);
    
    //charge selection
    n_noCharge++ ;
    if(pfMuons[m1].charge == pfMuons[m2].charge) n_sameCharge++;
    if(pfMuons[m1].charge != pfMuons[m2].charge) n_oppCharge++;
    //both muons should have opposite charge
    if(pfMuons[m1].charge == pfMuons[m2].charge) continue;

    //veto first or 2nd muon, if they are fake
    ///if(looseMuonVeto( m1, pfMuons, isPFlow) ) continue;
    ///if(looseMuonVeto( m2, pfMuons, isPFlow) ) continue;
    int count_muon = m_init.size();
    
    //---------------------------------------------------//
    //apply muon SF to eventWeights 
    //---------------------------------------------------//
    double lumi_BCDEF = 19.14; double lumi_GH = 16.23;	
    double lumi = lumi_BCDEF + lumi_GH;
    //get muon scale factor for fist muon
    //trigger 	
    double muSFtrig_BCDEF1 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFtrig_GH1 	= getMuonTrigSF(h2_trigSF_GH, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFtrig1 		= (muSFtrig_BCDEF1*lumi_BCDEF + muSFtrig_GH1*lumi_GH)/lumi; 

    //identification
    double muSFid_BCDEF1 	= getMuonSF(h2_idSF_BCDEF, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFid_GH1 		= getMuonSF(h2_idSF_GH, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFid1 		= (muSFid_BCDEF1*lumi_BCDEF + muSFid_GH1*lumi_GH)/lumi; 
    //isolation 
    double muSFiso_BCDEF1 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFiso_GH1 		= getMuonSF(h2_isoSF_GH, pfMuons[m1].p4.eta(), pfMuons[m1].p4.pt());
    double muSFiso1 		= (muSFiso_BCDEF1*lumi_BCDEF + muSFiso_GH1*lumi_GH)/lumi; 
    //tracking 
    double muSFtrack_BCDEF1 	= getMuonTrackSF(tg_trackSF_BCDEF, pfMuons[m1].p4.eta()); 
    double muSFtrack_GH1 	= getMuonTrackSF(tg_trackSF_GH, pfMuons[m1].p4.eta()); 
    double muSFtrack1 		= (muSFtrack_BCDEF1*lumi_BCDEF + muSFtrack_GH1*lumi_GH)/lumi;
    double muSF1 = muSFtrig1*muSFid1*muSFiso1*muSFtrack1;	
    
    //get muon scale factor for 2nd muon
    //trigger 	
    double muSFtrig_BCDEF2 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFtrig_GH2 	= getMuonTrigSF(h2_trigSF_GH, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFtrig2 		= (muSFtrig_BCDEF2*lumi_BCDEF + muSFtrig_GH2*lumi_GH)/lumi; 
    //identification
    double muSFid_BCDEF2 	= getMuonSF(h2_idSF_BCDEF, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFid_GH2 		= getMuonSF(h2_idSF_GH, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFid2 		= (muSFid_BCDEF2*lumi_BCDEF + muSFid_GH2*lumi_GH)/lumi; 
    //isolation 
    double muSFiso_BCDEF2 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFiso_GH2 		= getMuonSF(h2_isoSF_GH, pfMuons[m2].p4.eta(), pfMuons[m2].p4.pt());
    double muSFiso2 		= (muSFiso_BCDEF2*lumi_BCDEF + muSFiso_GH2*lumi_GH)/lumi; 
    //tracking 
    double muSFtrack_BCDEF2 	= getMuonTrackSF(tg_trackSF_BCDEF, pfMuons[m2].p4.eta()); 
    double muSFtrack_GH2 	= getMuonTrackSF(tg_trackSF_GH, pfMuons[m2].p4.eta()); 
    double muSFtrack2 		= (muSFtrack_BCDEF2*lumi_BCDEF + muSFtrack_GH2*lumi_GH)/lumi;
    double muSF2 = muSFtrig2*muSFid2*muSFiso2*muSFtrack2;	
    
    //combined SF
    double muSF =1.0;
    if(!ev->isData) muSF = muSF1*muSF2;
    evtWeight *= muSF;
    fillHisto(outFile_, cutflowType, "", "muonSF", 1000, 0, 100, muSF, 1 );
 
    //---------------------------------------------------//
    // Iso(<0.15) and Non-iso(>0.15) region 
    //---------------------------------------------------//
    bool isofound = false;
    string cutflowType_(cutflowType);
    cutflowType_ = cutflowType+"/Iso";
    if(m1RelIso > 0.15) continue;
    if(m2RelIso > 0.15) continue;
    
    //---------------------------------------------------//
    //get 4 vector for Z boson
    //---------------------------------------------------//
    MyLorentzVector vZ = pfMuons[m1].p4 + pfMuons[m2].p4;
    if(vZ.mass() < 60) continue;    
    double dR1 = 0.0;
    double dR2 = 0.0;
    int count_jets = j_final.size();
    //---------------------------------------------------//
    //Fill histos with for Control Plots
    //---------------------------------------------------//
    //fill histos for muon
    double muonPt1 = pfMuons[m1].p4.pt();
    double muonPt2 = pfMuons[m2].p4.pt(); 
    if(j_final.size()==0) continue;
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
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
    fillHisto(outFile_, cutflowType_, "ControlP","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_1stMu", 500, 0, 10000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_2ndMu", 500, 0, 10000, muonPt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ControlP", "ptMu1_ptMu2", 500, 0, 10000, muonPt1,500, 0, 10000, muonPt2, 1);
    fillHisto(outFile_, cutflowType_, "ControlP","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","m1RelIso", 100, 0, 1, m1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","m2RelIso", 100, 0, 1, m2RelIso, evtWeight );
   
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
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
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
    //fillHisto(outFile_, cutflowType_, "PreSel","pfJets_size", 15, 0, 15, pfJets.size(), evtWeight );
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///if(muonPt1 <100) continue;    
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "PreSel","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_1stMu", 500, 0, 10000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_2ndMu", 500, 0, 10000, muonPt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto2D(outFile_, cutflowType_,"PreSel", "ptMu1_ptMu2", 500, 0, 10000, muonPt1,500, 0, 10000, muonPt2, 1);
    fillHisto(outFile_, cutflowType_, "PreSel","m1RelIso", 100, 0, 1, m1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","m2RelIso", 100, 0, 1, m2RelIso, evtWeight );
   
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
    if(input_count_PreSel%10==0)
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
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
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
    if(allZjet.size()==0) continue;
    //apply tau21 scale factor
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging#2016_scale_factors_and_correctio
    if(!ev->isData) evtWeight *= 1.11; 

    MyLorentzVector vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[m1].p4;
    MyLorentzVector vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[m2].p4;
    
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "ZTag","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_1stMu", 500, 0, 10000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_2ndMu", 500, 0, 10000, muonPt2, evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ZTag", "ptMu1_ptMu2", 500, 0, 10000, muonPt1,500, 0, 10000, muonPt2, 1);
    fillHisto(outFile_, cutflowType_, "ZTag","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","m1RelIso", 100, 0, 1, m1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","m2RelIso", 100, 0, 1, m2RelIso, evtWeight );
   
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
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
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
    //if(i > 200000) break;
  }//event loop
  cout<<"Total events  = "<<nEntries<<endl;
  cout<<"Total events with negative weight = "<<n_negEvt<<endl;
  cout<<"Total events with positive weight = "<<n_posEvt<<endl;

  double effective_evt = (n_posEvt-n_negEvt);
  double amcnlo_weight = 1.0;
  if(effective_evt !=0) amcnlo_weight = effective_evt/nEntries;
  fillHisto(outFile_, cutflowType, "", "amcnlo_weight", 10, 0, 1, amcnlo_weight, 1 );
  fillHisto(outFile_, cutflowType, "", "noCharge", 10, -2, 2, 0, n_noCharge);
  fillHisto(outFile_, cutflowType, "", "oppCharge", 10, -2, 2, -1, n_oppCharge);
  fillHisto(outFile_, cutflowType, "", "sameCharge", 10, -2, 2, 1, n_sameCharge);
  fillHisto(outFile_, cutflowType, "", "oneTrig", 10, -2, 4, 1, n_oneTrig);
  fillHisto(outFile_, cutflowType, "", "twoTrig", 10, -2, 4, 2, n_twoTrig);
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJetsP_MuMC_20171104_Ntuple_1.root", "PF", ""); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");

  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/sthakur/ntuple_for2016Data_MuMC_20190117/MuMC_20190117/DYJetsToLL_M50_MuMC_20190117/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M50_MuMC_20190117/190117_091524/0000/DYJetsToLL_M50_MuMC_20190117_Ntuple_1.root" , "PF", "");

  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/sthakur/ntuple_for2016Data_MuMC_20190117/MuMC_20190117/TT_MuMC_20190117/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TT_MuMC_20190117/190117_092153/0000/TT_MuMC_20190117_Ntuple_92.root", "PF", "MuRunB2v2_MuData_20190117_Ntuple_87");
   //====================================
  //condor submission

  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 

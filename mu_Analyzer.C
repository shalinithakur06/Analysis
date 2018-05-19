
///////////////////////
// Muon Channel
///////////////////////

#include "Analyzer.h"
#include <map>

using namespace std;
void Analyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  //TString Filename_ = outFile+evtType+"_8Feb_signal_6000.root";
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  
  //check if the input file is MC or Data  
  Reader *evR_;  
  evR_ = new Reader();
  TFile *f_ = TFile::Open(url);
  int nEntries = evR_->AssignEventTreeFrom(f_);
  MyEvent *ev_;
  ev_ = evR_->GetNewEvent(1);

  CutFlowProcessor(url, myKey, "base", outFile_);
  //CutFlowProcessor(url, myKey, "baseLowMET", outFile_);
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  if(!ev_->isData){ 
    CutFlowProcessor(url, myKey, "JESPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JESMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERMinus", 	outFile_);

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
  int jes = 0, jer = 0, metuc = 0, bscale = 0, minMET =20, minMT =0;
  //to estimate unc in the data-driven qcd 
  bool isLowMET = false, isIso20 = false;

  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  //to estimate unc in the data-driven qcd 
  else if (cutflowType.Contains("LowMET"))isLowMET = true;
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35328;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  fillHisto(outFile_, cutflowType, "", "totalEvents", 10, 0, 10000000000, initialEvents, 1 );
  MyEvent *ev;
  
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  double kfCount = 0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
  
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      //k-factor weight (along with lumi weight) 
      if(sampleName.find("WJetsToLNu") != string::npos || sampleName.find("W1JetsToLNu") != string::npos || sampleName.find("W2JetsToLNu") != string::npos || sampleName.find("W3JetsToLNu") != string::npos || sampleName.find("W4JetsToLNu") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        double weightK = reweightHEPNUPWJets(hepNUP) * (lumiTotal/1000.0);
        if(i < 1){
        }
        evtWeight *= weightK;  
      }
      else if(sampleName.find("DYJetsToLL") != string::npos || sampleName.find("DY1JetsToLL") != string::npos || sampleName.find("DY2JetsToLL") != string::npos || sampleName.find("DY3JetsToLL") != string::npos || sampleName.find("DY4JetsToLL") != string::npos){
        int hepNUP = ev->sampleInfo.hepNUP;
        std::vector<int> hepIDUP = ev->sampleInfo.hepIDUP;
        std::vector<int> hepISTUP = ev->sampleInfo.hepISTUP;
        int countZ = 0;
        for(size_t p=0; p<hepIDUP.size(); p++){
          if(hepIDUP[p]==23 && hepISTUP[p]==2)
            countZ = countZ + 1;
        }
        if(countZ==0) hepNUP = hepNUP+1;
        double weightK = reweightHEPNUPDYJets(hepNUP) * (lumiTotal/1000.0);
        evtWeight *= weightK;  
      }
      //lumi weight
      else {
      double sampleWeight(1.0);
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      evtWeight *= sampleWeight; 
      fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      if(i < 1){
        }
      }
      //pileup weight
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
      float npu = pu[0];
      double weightPU = LumiWeights_.weight(npu);
      evtWeight *= weightPU;  
      }
    } 
    
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_Mu17_TrkIsoVVL") != string::npos) {
        passTrig = true;
      }
    }
    //if(!passTrig) continue;

    double nCutPass = 1.0;
    double nCutPass_NonIso = 1.0;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
   
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
    double u1 	= gRandom->Rndm();//used for rochester corrections
    double u2 	= gRandom->Rndm();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData, u1, u2, 0, 0);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    //clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    int nMuon = m_init.size();
    double pri_vtxs = Vertices[0].totVtx;
    if(nMuon != 2)continue;
    //fill histo for relative isolation of both muons
    fillHisto(outFile_, cutflowType+"/Iso", "", "RelIso", 100, 0, 1, pfMuons[m_init[0]].pfRelIso, evtWeight);
    fillHisto(outFile_, cutflowType+"/Iso", "", "RelIso", 100, 0, 1, pfMuons[m_init[1]].pfRelIso, evtWeight);
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "RelIso", 100, 0, 1, pfMuons[m_init[0]].pfRelIso, evtWeight);
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "RelIso", 100, 0, 1, pfMuons[m_init[1]].pfRelIso, evtWeight);
    
    //get Medium ID for first muon
    int m1 = m_init[0];
    bool passID1 = false;
    string input_file(url);
    if(input_file.find("RunG") != string::npos 
		    ||input_file.find("RunH") != string::npos)
	    passID1 = isMediumMuonGH(&pfMuons[m1], isPFlow);
    else
	    passID1 = isMediumMuon(&pfMuons[m1], isPFlow);
    
    //get Medium ID for 2nd muon
    int m2 = m_init[1];
    bool passID2 = false;
    if(input_file.find("RunG") != string::npos 
		    ||input_file.find("RunH") != string::npos)
	    passID2 = isMediumMuonGH(&pfMuons[m2], isPFlow);
    else
	    passID2 = isMediumMuon(&pfMuons[m2], isPFlow);
    
    //Make sure that both muon satisfy medium ID
    if(!passID1) continue;
    if(!passID2) continue;
    //veto first or 2nd muon, if they are fake
    ///if(looseMuonVeto( m1, pfMuons, isPFlow) ) continue;
    ///if(looseMuonVeto( m2, pfMuons, isPFlow) ) continue;
    //both muons should have opposite charge
    if(pfMuons[m1].charge == pfMuons[m2].charge) continue;

    nCutPass++; 
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
     
    //events should not have any electron
    ///if(looseElectronVeto(-1, pfElectrons, Vertices[0], isPFlow)) continue;
    ///nCutPass++; 
    ///nCutPass_NonIso++;
    ///fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
    int count_muon = m_init.size();
    ///int muCharge = pfMuons[m_i].charge;
    
    //---------------------------------------------------//
    //apply muon SF to eventWeights 
    //---------------------------------------------------//
    double lumi_BCDEF = 19.09; double lumi_GH = 16.22;	
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
    if(!ev->isData){
      muSF = muSF1*muSF2;	
    }
    evtWeight *= muSF;
    nCutPass++;
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
 
    //---------------------------------------------------//
    // Iso(<0.15) and Non-iso(>0.15) region 
    //---------------------------------------------------//
    bool noisofound = false;
    bool isofound = false;
    double tmp_iso = pfMuons[m1].pfRelIso;
    fillHisto(outFile_, cutflowType, "", "RelIso_mu", 100, 0, 1, tmp_iso, evtWeight );
    string cutflowType_(cutflowType);
    cutflowType_ = cutflowType+"/Iso";
    //if(tmp_iso <= 0.15) cutflowType_ = cutflowType+"/Iso";
    //if(tmp_iso > 0.15 && tmp_iso <= 0.4) cutflowType_ = cutflowType+"/NonIso";
    //if(tmp_iso > 0.4) continue;
    if(tmp_iso > 0.15) continue;
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    double mRelIso = pfMuons[m1].pfRelIso;
    
    //---------------------------------------------------//
    //get 4 vector for Z boson
    //---------------------------------------------------//
    MyLorentzVector vZ = pfMuons[m1].p4 + pfMuons[m2].p4;
    if(vZ.mass() < 60) continue;    
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    double dR1 = 0.0;
    double dR2 = 0.0;
    int count_jets = j_final.size();
    //---------------------------------------------------//
    //Fill histos with for Control Plots
    //---------------------------------------------------//
    //fill histos for muon
    double muonPt1 = muPtWithRochCorr(&pfMuons[m1], ev->isData, u1, u2, 0, 0);
    double muonPt2 = muPtWithRochCorr(&pfMuons[m2], ev->isData, u1, u2, 0, 0);
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
      fillHisto(outFile_, cutflowType_, "ControlP","pt_jet", 100, 0, 1000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","ak8Pmass", 500, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "ControlP","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_1stMu", 500, 0, 5000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt_2ndMu", 500, 0, 5000, muonPt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","final_RelIso_mu", 100, 0, 1, mRelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ControlP","pt_Z",  500, 0, 5000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","mjj", 200, 0, 1000, vZ.M(), evtWeight );
    
    
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "ControlP","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "ControlP","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }

    //---------------------------------------------------//
    //Fill histos with pre-selection
    //---------------------------------------------------//
    if(vZ.Pt() < 100) continue;    
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      if(jetPt <= 100) continue;    
      if(fabs(pfJets[ind_jet].p4.eta()) >= 2.4) continue;    
      if(pfJets[ind_jet].ak8Pmass <= 40) continue;    
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
      if(dR1 < 0.8 || dR2 < 0.8) continue;    
      fillHisto(outFile_, cutflowType_, "PreSel","dR1", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR2", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","pt_jet", 100, 0, 1000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","ak8Pmass", 500, 0, 5000, pfJets[ind_jet].ak8Pmass, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","ak8Tau21", 50, 0, 5, pfJets[ind_jet].ak8Tau2/pfJets[ind_jet].ak8Tau1, evtWeight );
    }
    
    fillHisto(outFile_, cutflowType_, "PreSel","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///if(muonPt1 <100) continue;    
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "PreSel","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_1stMu", 500, 0, 5000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_2ndMu", 500, 0, 5000, muonPt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","final_RelIso_mu", 100, 0, 1, mRelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "PreSel","pt_Z",  500, 0, 5000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_Z", 50, -5, -5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","phi_Z", 50, -5, -5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","mjj", 200, 0, 1000, vZ.M(), evtWeight );
    
    
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "PreSel","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "PreSel","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    input_count_PreSel++;
    if(input_count_PreSel%10==0)
    cout << "input count after PreSel: "<< input_count_PreSel << endl;

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
      if(vZ.M()> 200 && jetPt >200 && ak8Pmass_ > 70 && ak8Pmass_ < 110 && ak8Tau21 < 0.5){ 
	allZjet.push_back(ijet);
      }
    }
    
    if(allZjet.size()!=1) continue;
    MyLorentzVector vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[m1].p4;
    MyLorentzVector vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfMuons[m2].p4;
    
    //MyLorentzVector vZmax = pfJets[j_final[0]].p4 + pfJets[j_final[1]].p4 + pfMuons[m1].p4;
    //MyLorentzVector vZmin = pfJets[j_final[0]].p4 + pfJets[j_final[1]].p4 + pfMuons[m2].p4;
    
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "ZTag","multi_mu",  15, 0.5, 15.5, count_muon, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_1stMu", 500, 0, 5000, muonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_2ndMu", 500, 0, 5000, muonPt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_1stMu", 50, -5, 5, pfMuons[m1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_2ndMu", 50, -5, 5, pfMuons[m2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","final_RelIso_mu", 100, 0, 1, mRelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ZTag","pt_Z",  500, 0, 5000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_Z",  500, 0, 5000, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","phi_Z",  500, 0, 5000, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","mjj", 200, 0, 1000, vZ.M(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","mjj_max", 200, 0, 10000, vZmax.M(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","mjj_min", 200, 0, 10000, vZmin.M(), evtWeight );
    
    //fill histos for jets
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      dR1 = DeltaR(pfJets[ind_jet].p4, pfMuons[m1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfMuons[m2].p4);
      fillHisto(outFile_, cutflowType_, "PreSel","dR1", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR2", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR1, evtWeight );
      fillHisto(outFile_, cutflowType_, "PreSel","dR", 100, 0, 10, dR2, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","pt_jet", 100, 0, 1000, jetPt, evtWeight );
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
    //if(i > 2000) break;
  }//event loop
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJetsP_MuMC_20171104_Ntuple_1.root", "PF", ""); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  CutFlowAnalysis("All_Mu2000_Ntuple.root", "PF", ""); 

  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/sthakur/ntuple_MuMC_20180224/MuMC_20180224/DYJetsToLL_MuMC_20180224/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_MuMC_20180224/180224_185201/0000/DYJetsToLL_MuMC_20180224_Ntuple_1.root", "PF", "");

  //====================================
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 

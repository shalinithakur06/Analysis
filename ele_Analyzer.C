
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
  //CutFlowProcessor(url, myKey, "baseLowMET", outFile_);
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  /*
  if(!ev_->isData){ 
    CutFlowProcessor(url, myKey, "JESPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JESMinus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERPlus", 	outFile_);
    CutFlowProcessor(url, myKey, "JERMinus", 	outFile_);

  }
  */
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
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35381;
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
  double n_negEvt = 0.0;
  double n_posEvt = 0.0;
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
      }
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
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData, u1, u2, 0, 0);
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
    ///if(pfElectrons[e1].charge == pfElectrons[e2].charge) continue; // to be consistent with other AN
     
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
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
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
    fillHisto2D(outFile_, cutflowType_,"ControlP", "ptMu1_ptEle", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "ControlP","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ControlP","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","mjj", 500, 0, 10000, vZ.M(), evtWeight );
    
    
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
      dR1 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e1].p4);
      dR2 = DeltaR(pfJets[ind_jet].p4, pfElectrons[e2].p4);
      if(dR1 < 0.8 || dR2 < 0.8) continue;    
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
    fillHisto(outFile_, cutflowType_, "PreSel","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pfJets_size", 15, 0, 15, pfJets.size(), evtWeight );
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    ///if(muonPt1 <100) continue;    
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "PreSel","multi_Ele",  15, 0.5, 15.5, count_ele, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_1stEle", 500, 0, 5000, elePt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","pt_2ndEle", 500, 0, 10000, elePt2, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_1stEle", 50, -5, 5, pfElectrons[e1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto2D(outFile_, cutflowType_,"PreSel", "ptMu1_ptEle", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "PreSel","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "PreSel","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "PreSel","mjj", 500, 0, 10000, vZ.M(), evtWeight );
    
    
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
    
    if(allZjet.size()==0) continue;
    MyLorentzVector vZmax =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[e1].p4;
    MyLorentzVector vZmin =  pfJets[j_final[allZjet[0]]].p4 + pfElectrons[e2].p4;
    
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    //fill histos for muon
    fillHisto(outFile_, cutflowType_, "ZTag","multi_Ele",  15, 0.5, 15.5, count_ele, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_1stEle", 500, 0, 10000, elePt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","pt_2ndEle", 500, 0, 10000, elePt2, evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ZTag", "ptMu1_ptEle", 500, 0, 10000, elePt1,500, 0, 10000, elePt2, 1);
    fillHisto(outFile_, cutflowType_, "ZTag","eta_1stEle", 50, -5, 5, pfElectrons[e1].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_2ndEle", 50, -5, 5, pfElectrons[e2].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","e1RelIso", 100, 0, 1, e1RelIso, evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","e2RelIso", 100, 0, 1, e2RelIso, evtWeight );
   
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ZTag","pt_Z",  500, 0, 10000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","eta_Z",  500, 0, 10000, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ZTag","phi_Z",  500, 0, 10000, vZ.Phi(), evtWeight );
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
    double mlZmax_low  = 0.0;
    double mlZmax_high = 0.0;
    double mlZmin_low  = 0.0;
    double mlZmin_high = 0.0;
    if(mlZmax > 440 && mlZmin < 300){
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_sig250", 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max_sig250", 500, 0, 10000, mlZmax, evtWeight );
    }
    if(mlZmax > 1300 && mlZmin < 1700){
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_sig1500", 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max_sig1500", 500, 0, 10000, mlZmax, evtWeight );
    }
    if(mlZmax > 1300 && mlZmin < 2200){
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_sig2000", 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max_sig2000", 500, 0, 10000, mlZmax, evtWeight );
    }
    if(mlZmax > 1300 && mlZmin < 2700){
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_sig2500", 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max_sig2500", 500, 0, 10000, mlZmax, evtWeight );
    }
    if(mlZmax > 1300 && mlZmin < 4200){
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_min_sig4000", 500, 0, 10000, mlZmin, evtWeight );
      fillHisto(outFile_, cutflowType_, "ZTag","mlZ_max_sig4000", 500, 0, 10000, mlZmax, evtWeight );
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
    //if(i > 200000) break;
  }//event loop
  cout<<"Total events  = "<<nEntries<<endl;
  cout<<"Total events with negative weight = "<<n_negEvt<<endl;
  cout<<"Total events with positive weight = "<<n_posEvt<<endl;
  double effective_evt = (n_posEvt-n_negEvt)/nEntries;
  double amcnlo_weght = 1.0;
  if(effective_evt !=0) amcnlo_weght = 1/effective_evt;
  fillHisto(outFile_, cutflowType, "", "amcnlo_weght", 10, 0, 2, amcnlo_weght, 1 );
 
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJetsP_MuMC_20171104_Ntuple_1.root", "PF", ""); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  CutFlowAnalysis("outFile_.root", "PF", ""); 
   //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/sthakur/ntuple_EleMC_20180505/EleMC_20180505/DYJetsToLLamcatnlo_EleMC_20180505/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLLamcatnlo_EleMC_20180505/180505_150933/0000/outFileEle_2000_1.root", "PF", "");

  //====================================
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 

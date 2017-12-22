
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
  //to estimate unc in the data-driven qcd 
  //CutFlowProcessor(url, myKey, "baseIso20HighMET", outFile_);
  //CutFlowProcessor(url, myKey, "baseIso20LowMET", outFile_);
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
  int input_count = 0;
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
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  else if (cutflowType.Contains("bTagMinus"))bscale = -1; 
  //to estimate unc in the data-driven qcd 
  else if (cutflowType.Contains("baseIso")){ 
    if (cutflowType.Contains("Iso20HighMET"))isIso20 = true; 
    if (cutflowType.Contains("Iso20LowMET")){isIso20 = true; isLowMET= true;}
  }
  else if (cutflowType.Contains("LowMET"))isLowMET = true;
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 34698+758;
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
        if(i < 1){
        }
      }
      //lumi weight
      else {
      double sampleWeight(1.0);
      sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
      evtWeight *= sampleWeight; 
      fillHisto(outFile_, cutflowType, "", "MuPt30JetPt25MET20MT20", 10, 0, 1000, sampleWeight, 1 );
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
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("HLT_IsoMu24") != string::npos) {
        passTrig = true;
      }
    }
    if(!passTrig){
    //cout << "not satisfying trigger" << endl;
      continue;
    }
    double nCutPass = 1.0;
    double nCutPass_NonIso = 1.0;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
   
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

    // preselect objects 
    vector<int> m_init; m_init.clear();
    double u1 	= gRandom->Rndm();//used for rochester corrections
    double u2 	= gRandom->Rndm();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData, u1, u2, 0, 0);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    if(m_init.size() > 0){
      int m_i_noiso = m_init[0];
      double mRelIso_no_iso = pfMuons[m_i_noiso].pfRelIso;
      fillHisto(outFile_, cutflowType+"/Iso", "", "RelIso", 100, 0, 1, mRelIso_no_iso, evtWeight);
      fillHisto(outFile_, cutflowType+"/NonIso", "", "RelIso", 100, 0, 1, mRelIso_no_iso, evtWeight);
    }
    int nMuon = m_init.size();
    double pri_vtxs = Vertices[0].totVtx;
    if(nMuon != 1)continue;
    //check if 0th muon has mediumMuon ID
    int m_i = m_init[0];
    bool passID = false;
    string input_file(url);
    if(input_file.find("RunG") != string::npos 
		    ||input_file.find("RunH") != string::npos)
	    passID = isMediumMuonGH(&pfMuons[m_i], isPFlow);
    else
	    passID = isMediumMuon(&pfMuons[m_i], isPFlow);
    if(!passID) continue;
    
    //veto 0th muon, if other muons are stroger than the 0th.
    //we veto 0th if it is a fake muon
    if(looseMuonVeto( m_i, pfMuons, isPFlow) ) continue;
    nCutPass = 2; 
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
     
    //events should not have any electron
    if(looseElectronVeto(-1, pfElectrons, Vertices[0], isPFlow)) continue;
    nCutPass = 3; 
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
    int count_muon = m_init.size();
    ///int muCharge = pfMuons[m_i].charge;
    
    //---------------------------------------------------//
    //apply muon SF to eventWeights 
    //---------------------------------------------------//
    double lumi_BCDEF = 18.85+0.337; double lumi_GH = 15.84+0.420;	
    double lumi = lumi_BCDEF + lumi_GH;
    //trigger 	
    double muSFtrig_BCDEF 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFtrig_GH 		= getMuonTrigSF(h2_trigSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFtrig 		= (muSFtrig_BCDEF*lumi_BCDEF + muSFtrig_GH*lumi_GH)/lumi; 
    //identification
    double muSFid_BCDEF 	= getMuonSF(h2_idSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFid_GH 		= getMuonSF(h2_idSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFid 		= (muSFid_BCDEF*lumi_BCDEF + muSFid_GH*lumi_GH)/lumi; 
    //isolation 
    double muSFiso_BCDEF 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFiso_GH 		= getMuonSF(h2_isoSF_GH, pfMuons[m_i].p4.eta(), pfMuons[m_i].p4.pt());
    double muSFiso 		= (muSFiso_BCDEF*lumi_BCDEF + muSFiso_GH*lumi_GH)/lumi; 
    //tracking 
    double muSFtrack_BCDEF 	= getMuonTrackSF(tg_trackSF_BCDEF, pfMuons[m_i].p4.eta()); 
    double muSFtrack_GH 	= getMuonTrackSF(tg_trackSF_GH, pfMuons[m_i].p4.eta()); 
    double muSFtrack 		= (muSFtrack_BCDEF*lumi_BCDEF + muSFtrack_GH*lumi_GH)/lumi;

    //combined SF
    double muSF =1.0;
    if(!ev->isData){
      muSF = muSFtrig*muSFid*muSFiso*muSFtrack;	
    }
    evtWeight *= muSF;
    nCutPass = 4;
    nCutPass_NonIso++;
    fillHisto(outFile_, cutflowType+"/Iso", "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    fillHisto(outFile_, cutflowType+"/NonIso", "", "cutflow", 20, 0.5, 20.5, nCutPass_NonIso, evtWeight );
    //---------------------------------------------------//
    // Iso(<0.12) and Non-iso(>0.12) region 
    //---------------------------------------------------//
    bool noisofound = false;
    bool isofound = false;
    double tmp_iso = pfMuons[m_i].pfRelIso;
    fillHisto(outFile_, cutflowType, "", "RelIso_mu", 100, 0, 1, tmp_iso, evtWeight );
    string cutflowType_(cutflowType);
    if(isIso20){
      if(tmp_iso <= 0.20) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.20 && tmp_iso <= 0.4) cutflowType_ = cutflowType+"/NonIso";
    }
    else{
      if(tmp_iso <= 0.15) cutflowType_ = cutflowType+"/Iso";
      if(tmp_iso > 0.15 && tmp_iso <= 0.4) cutflowType_ = cutflowType+"/NonIso";
    }
    if(tmp_iso > 0.4) continue;
    nCutPass = 5;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    double mRelIso = pfMuons[m_i].pfRelIso;
    //double muonPt = pfMuons[m_i].p4.pt();
    double muonPt = muPtWithRochCorr(&pfMuons[m_i], ev->isData, u1, u2, 0, 0);
    
    //---------------------------------------------------//
    // Apply Jet Selection
    //---------------------------------------------------//
    int count_jets = j_final.size();
    if(count_jets < 4)continue;  // events should have 4 or more jets
    nCutPass = 6;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    //apply MET selection   
    //---------------------------------------------------//
    double metPt = 0; 
    metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    //if(metuc==0)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    //else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    
    if(isLowMET){
      if(metPt > minMET) continue;  
    }
    else if(metPt < minMET) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
    nCutPass = 7;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //Transverse mass b/w lepton and MET
    double deltaPhi(0);
    //leptonPt = TMath::Abs(pfMuons[m_i].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_i].p4, met.p4);
    double mt = sqrt (  2*muonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    if(mt < minMT) continue; // mt cuts
    nCutPass = 8;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
    
    //---------------------------------------------------//
    // add set of plots after BTag:
    //---------------------------------------------------//
    //fillHisto("pt_mu", cutflowType_+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
    fillHisto(outFile_, cutflowType_, "BTag","pt_mu", 50, 0, 500, muonPt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","eta_mu", 50, -5, 5, pfMuons[m_i].p4.eta(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","phi_mu", 50, -5, 5, pfMuons[m_i].p4.phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_RelIso_mu", 100, 0, 1, mRelIso, evtWeight );
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      fillHisto(outFile_, cutflowType_, "BTag","pt_jet", 50, 0, 500, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "BTag","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "BTag","final_multi_jet", 15, 0, 15, count_jets, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","final_pt_met", 50, 0, 500, metPt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","wmt", 50, 0, 500, mt, evtWeight );
    fillHisto(outFile_, cutflowType_, "BTag","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "BTag","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    input_count++;
    if(input_count%10==0)
    cout << "input count iso: "<< input_count << endl;
    //if(i > 2000) break;
  }//event loop
  cout<<"kfCount = "<<kfCount<<endl;
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJetsP_MuMC_20171104_Ntuple_1.root", "PF", ""); 
  CutFlowAnalysis("outFile_.root", "PF", ""); 
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/", "PF", "");
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_kfitM_20171115/MuMC_20171115/TTJetsP_MuMC_20171115/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTJetsP_MuMC_20171115/171115_113943/0000/TTJetsP_MuMC_20171115_Ntuple_99.root", "PF", "");
  
  //====================================
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  //====================================
} 

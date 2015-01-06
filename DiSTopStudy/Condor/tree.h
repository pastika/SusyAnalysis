
   // Declaration of leaf types
   UInt_t          t_EvtRun;
   UInt_t          t_EvtLS;
   UInt_t          t_EvtEvent;
   Int_t           t_NVertices;
   Double_t        t_PUWeight;
   Double_t        t_PFMetPx;
   Double_t        t_PFMetPy;
   Double_t        t_PFMetE;
   Double_t        t_PFMetPhi;
   Double_t        t_PFMetSignificance;
   Int_t           NJets;
   vector<double>  *t_PFJetPt;
   vector<double>  *t_PFJetEta;
   vector<double>  *t_PFJetPhi;
   vector<double>  *t_PFJetE;
   vector<double>  *t_PFJetBTag;
   Int_t           t_NJetsPt30Eta2p5;
   Int_t           t_NJetsPt30Eta5p0;
   Int_t           t_NJetsPt50Eta2p5;
   Int_t           t_NJetsPt50Eta5p0;
   Double_t        t_PFht;
   Double_t        t_PFmht;
   Double_t        t_PFmphi;
   Int_t           t_NVetoMuon;
   Int_t           t_NVetoEle;
   vector<double>  *t_PFMuonPt;
   vector<double>  *t_PFMuonEta;
   vector<double>  *t_PFMuonPhi;
   vector<double>  *t_PFMuonE;
   vector<double>  *t_PFMuonCh;
   vector<double>  *t_PFMuonChiSq;
   vector<double>  *t_PFMuonValidMuonHits;
   vector<double>  *t_PFMuonMatchedStations;
   vector<double>  *t_PFMuonValidPixelHits;
   vector<double>  *t_PFMuonTrackerLayerMeasured;
   vector<double>  *t_PFMuonCharHadEt;
   vector<double>  *t_PFMuonNeutHadEt;
   vector<double>  *t_PFMuonPhotEt;
   vector<double>  *t_PFMuonSumPUPt;
   vector<double>  *t_PFMuonRelIso;
   vector<double>  *t_PFMuonCharHadEt03;
   vector<double>  *t_PFMuonNeutHadEt03;
   vector<double>  *t_PFMuonPhotEt03;
   vector<double>  *t_PFMuonSumPUPt03;
   vector<double>  *t_PFMuonRelIso03;
   vector<double>  *t_PFMuonDirIso01;
   vector<double>  *t_PFMuonDirIso02;
   vector<double>  *t_PFMuonDirIso03;
   vector<double>  *t_PFMuonDirIso04;
   vector<double>  *t_PFMuonDirIso05;
   vector<bool>    *t_PFMuonID;
   vector<bool>    *t_PFMuonVtxAss;
   vector<bool>    *t_PFMuonIsolated;
   vector<bool>    *t_PFMuonIDVtxIso;
   vector<bool>    *t_PFMuonIsWdau;
   vector<bool>    *t_PFMuonIsBdau;
   vector<double>  *t_PFMuonGenMuPt;
   vector<double>  *t_PFMuonGenMuEta;
   vector<double>  *t_PFMuonGenMuPhi;
   vector<double>  *t_PFElecPt;
   vector<double>  *t_PFElecEta;
   vector<double>  *t_PFElecPhi;
   vector<double>  *t_PFElecE;
   vector<double>  *t_PFElecCh;
   vector<bool>    *t_PFEleIsEB;
   vector<bool>    *t_PFEleIsEE;
   vector<double>  *t_PFEleDEtaIn;
   vector<double>  *t_PFEleDPhiIn;
   vector<double>  *t_PFEleSigIEtaIEta;
   vector<double>  *t_PFEleHOE;
   vector<double>  *t_PFEleD0Vtx;
   vector<double>  *t_PFEleZVtx;
   vector<double>  *t_PFEleIsoCH;
   vector<double>  *t_PFEleIsoEM;
   vector<double>  *t_PFEleIsoNH;
   vector<double>  *t_PFEleIsoRho;
   vector<double>  *t_PFEleIsoEffA;
   vector<double>  *t_PFEleRelIso;
   vector<bool>    *t_PFEleID;
   vector<bool>    *t_PFEleIsolated;
   vector<bool>    *t_PFEleIDIso;
   vector<double>  *t_PFEleDirIso01;
   vector<double>  *t_PFEleDirIso02;
   vector<double>  *t_PFEleDirIso03;
   vector<double>  *t_PFEleDirIso04;
   vector<double>  *t_PFEleDirIso05;
   vector<bool>    *t_PFEleIsWdau;
   vector<bool>    *t_PFEleIsBdau;
   vector<double>  *t_PFEleGenElePt;
   vector<double>  *t_PFEleGenEleEta;
   vector<double>  *t_PFEleGenElePhi;
   vector<int>     *t_genPartPdgId;
   vector<int>     *t_genPartStatus;
   vector<double>  *t_genPartPt;
   vector<double>  *t_genPartEta;
   vector<double>  *t_genPartPhi;
   vector<double>  *t_genPartE;
   vector<int>     *t_genPartDecayMode;
   Double_t        t_fatJet1_m123;
   Double_t        t_fatJet1_pt123;
   Double_t        t_fatJet1_eta123;
   Double_t        t_fatJet1_phi123;
   Double_t        t_fatJet1_m12;
   Double_t        t_fatJet1_m13;
   Double_t        t_fatJet1_m23;
   Double_t        t_fatJet2_m123;
   Double_t        t_fatJet2_pt123;
   Double_t        t_fatJet2_eta123;
   Double_t        t_fatJet2_phi123;
   Double_t        t_fatJet2_m12;
   Double_t        t_fatJet2_m13;
   Double_t        t_fatJet2_m23;
   UInt_t          t_nFatJet;
   UInt_t          t_fatJet1_cond1;
   UInt_t          t_fatJet1_cond2;
   UInt_t          t_fatJet1_cond3;
   UInt_t          t_fatJet2_cond1;
   UInt_t          t_fatJet2_cond2;
   UInt_t          t_fatJet2_cond3;
   vector<vector<double> > *t_fatJetSubJetPt;
   vector<vector<double> > *t_fatJetSubJetEta;
   vector<vector<double> > *t_fatJetSubJetPhi;
   vector<vector<double> > *t_fatJetSubJetE;
   Double_t        t_DTWMass1;
   Double_t        t_DTWMass2;
   Double_t        t_DTMassChiSqr1;
   Double_t        t_DTMassChiSqr2;
   Double_t        t_DTTopMass1;
   Double_t        t_DTTopMass2;
   Double_t        t_DTTotalSystemMass;
   Double_t        t_DTTopPt1;
   Double_t        t_DTTopPt2;
   Double_t        t_DTTotalSystemPt;
   Double_t        t_DTTopEta1;
   Double_t        t_DTTopEta2;
   Double_t        t_DTTotalSystemEta;
   Double_t        t_DTTopPhi1;
   Double_t        t_DTTopPhi2;
   Double_t        t_DTTotalSystemPhi;
   Int_t           t_DTindexB1;
   Int_t           t_DTindexB2;
   Int_t           t_DTindexW11;
   Int_t           t_DTindexW12;
   Int_t           t_DTindexW21;
   Int_t           t_DTindexW22;
   Double_t        t_DTdPhi1;
   Double_t        t_DTdPhi2;
   Double_t        t_DTTopMu1;
   Double_t        t_DTTopMu2;

   // List of branches
   TBranch        *b_t_EvtRun;   //!
   TBranch        *b_t_EvtLS;   //!
   TBranch        *b_t_EvtEvent;   //!
   TBranch        *b_t_NVertices;   //!
   TBranch        *b_t_PUWeight;   //!
   TBranch        *b_t_PFMetPx;   //!
   TBranch        *b_t_PFMetPy;   //!
   TBranch        *b_t_PFMetE;   //!
   TBranch        *b_t_PFMetPhi;   //!
   TBranch        *b_t_PFMetSignificance;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_t_PFJetPt;   //!
   TBranch        *b_t_PFJetEta;   //!
   TBranch        *b_t_PFJetPhi;   //!
   TBranch        *b_t_PFJetE;   //!
   TBranch        *b_t_PFJetBTag;   //!
   TBranch        *b_t_NJetsPt30Eta2p5;   //!
   TBranch        *b_t_NJetsPt30Eta5p0;   //!
   TBranch        *b_t_NJetsPt50Eta2p5;   //!
   TBranch        *b_t_NJetsPt50Eta5p0;   //!
   TBranch        *b_t_PFht;   //!
   TBranch        *b_t_PFmht;   //!
   TBranch        *b_t_PFmphi;   //!
   TBranch        *b_t_NVetoMuon;   //!
   TBranch        *b_t_NVetoEle;   //!
   TBranch        *b_t_PFMuonPt;   //!
   TBranch        *b_t_PFMuonEta;   //!
   TBranch        *b_t_PFMuonPhi;   //!
   TBranch        *b_t_PFMuonE;   //!
   TBranch        *b_t_PFMuonCh;   //!
   TBranch        *b_t_PFMuonChiSq;   //!
   TBranch        *b_t_PFMuonValidMuonHits;   //!
   TBranch        *b_t_PFMuonMatchedStations;   //!
   TBranch        *b_t_PFMuonValidPixelHits;   //!
   TBranch        *b_t_PFMuonTrackerLayerMeasured;   //!
   TBranch        *b_t_PFMuonCharHadEt;   //!
   TBranch        *b_t_PFMuonNeutHadEt;   //!
   TBranch        *b_t_PFMuonPhotEt;   //!
   TBranch        *b_t_PFMuonSumPUPt;   //!
   TBranch        *b_t_PFMuonRelIso;   //!
   TBranch        *b_t_PFMuonCharHadEt03;   //!
   TBranch        *b_t_PFMuonNeutHadEt03;   //!
   TBranch        *b_t_PFMuonPhotEt03;   //!
   TBranch        *b_t_PFMuonSumPUPt03;   //!
   TBranch        *b_t_PFMuonRelIso03;   //!
   TBranch        *b_t_PFMuonDirIso01;   //!
   TBranch        *b_t_PFMuonDirIso02;   //!
   TBranch        *b_t_PFMuonDirIso03;   //!
   TBranch        *b_t_PFMuonDirIso04;   //!
   TBranch        *b_t_PFMuonDirIso05;   //!
   TBranch        *b_t_PFMuonID;   //!
   TBranch        *b_t_PFMuonVtxAss;   //!
   TBranch        *b_t_PFMuonIsolated;   //!
   TBranch        *b_t_PFMuonIDVtxIso;   //!
   TBranch        *b_t_PFMuonIsWdau;   //!
   TBranch        *b_t_PFMuonIsBdau;   //!
   TBranch        *b_t_PFMuonGenMuPt;   //!
   TBranch        *b_t_PFMuonGenMuEta;   //!
   TBranch        *b_t_PFMuonGenMuPhi;   //!
   TBranch        *b_t_PFElecPt;   //!
   TBranch        *b_t_PFElecEta;   //!
   TBranch        *b_t_PFElecPhi;   //!
   TBranch        *b_t_PFElecE;   //!
   TBranch        *b_t_PFElecCh;   //!
   TBranch        *b_t_PFEleIsEB;   //!
   TBranch        *b_t_PFEleIsEE;   //!
   TBranch        *b_t_PFEleDEtaIn;   //!
   TBranch        *b_t_PFEleDPhiIn;   //!
   TBranch        *b_t_PFEleSigIEtaIEta;   //!
   TBranch        *b_t_PFEleHOE;   //!
   TBranch        *b_t_PFEleD0Vtx;   //!
   TBranch        *b_t_PFEleZVtx;   //!
   TBranch        *b_t_PFEleIsoCH;   //!
   TBranch        *b_t_PFEleIsoEM;   //!
   TBranch        *b_t_PFEleIsoNH;   //!
   TBranch        *b_t_PFEleIsoRho;   //!
   TBranch        *b_t_PFEleIsoEffA;   //!
   TBranch        *b_t_PFEleRelIso;   //!
   TBranch        *b_t_PFEleID;   //!
   TBranch        *b_t_PFEleIsolated;   //!
   TBranch        *b_t_PFEleIDIso;   //!
   TBranch        *b_t_PFEleDirIso01;   //!
   TBranch        *b_t_PFEleDirIso02;   //!
   TBranch        *b_t_PFEleDirIso03;   //!
   TBranch        *b_t_PFEleDirIso04;   //!
   TBranch        *b_t_PFEleDirIso05;   //!
   TBranch        *b_t_PFEleIsWdau;   //!
   TBranch        *b_t_PFEleIsBdau;   //!
   TBranch        *b_t_PFEleGenElePt;   //!
   TBranch        *b_t_PFEleGenEleEta;   //!
   TBranch        *b_t_PFEleGenElePhi;   //!
   TBranch        *b_t_genPartPdgId;   //!
   TBranch        *b_t_genPartStatus;   //!
   TBranch        *b_t_genPartPt;   //!
   TBranch        *b_t_genPartEta;   //!
   TBranch        *b_t_genPartPhi;   //!
   TBranch        *b_t_genPartE;   //!
   TBranch        *b_t_genPartDecayMode;   //!
   TBranch        *b_fatJet1_m123;   //!
   TBranch        *b_fatJet1_pt123;   //!
   TBranch        *b_fatJet1_eta123;   //!
   TBranch        *b_fatJet1_phi123;   //!
   TBranch        *b_fatJet1_m12;   //!
   TBranch        *b_fatJet1_m13;   //!
   TBranch        *b_fatJet1_m23;   //!
   TBranch        *b_fatJet2_m123;   //!
   TBranch        *b_fatJet2_pt123;   //!
   TBranch        *b_fatJet2_eta123;   //!
   TBranch        *b_fatJet2_phi123;   //!
   TBranch        *b_fatJet2_m12;   //!
   TBranch        *b_fatJet2_m13;   //!
   TBranch        *b_fatJet2_m23;   //!
   TBranch        *b_t_nFatJet;   //!
   TBranch        *b_t_fatJet1_cond1;   //!
   TBranch        *b_t_fatJet1_cond2;   //!
   TBranch        *b_t_fatJet1_cond3;   //!
   TBranch        *b_t_fatJet2_cond1;   //!
   TBranch        *b_t_fatJet2_cond2;   //!
   TBranch        *b_t_fatJet2_cond3;   //!
   TBranch        *b_t_fatJetSubJetPt;   //!
   TBranch        *b_t_fatJetSubJetEta;   //!
   TBranch        *b_t_fatJetSubJetPhi;   //!
   TBranch        *b_t_fatJetSubJetE;   //!
   TBranch        *b_t_DTWMass1;   //!
   TBranch        *b_t_DTWMass2;   //!
   TBranch        *b_t_DTMassChiSqr1;   //!
   TBranch        *b_t_DTMassChiSqr2;   //!
   TBranch        *b_t_DTTopMass1;   //!
   TBranch        *b_t_DTTopMass2;   //!
   TBranch        *b_t_DTTotalSystemMass;   //!
   TBranch        *b_t_DTTopPt1;   //!
   TBranch        *b_t_DTTopPt2;   //!
   TBranch        *b_t_DTTotalSystemPt;   //!
   TBranch        *b_t_DTTopEta1;   //!
   TBranch        *b_t_DTTopEta2;   //!
   TBranch        *b_t_DTTotalSystemEta;   //!
   TBranch        *b_t_DTTopPhi1;   //!
   TBranch        *b_t_DTTopPhi2;   //!
   TBranch        *b_t_DTTotalSystemPhi;   //!
   TBranch        *b_t_DTindexB1;   //!
   TBranch        *b_t_DTindexB2;   //!
   TBranch        *b_t_DTindexW11;   //!
   TBranch        *b_t_DTindexW12;   //!
   TBranch        *b_t_DTindexW21;   //!
   TBranch        *b_t_DTindexW22;   //!
   TBranch        *b_t_DTdPhi1;   //!
   TBranch        *b_t_DTdPhi2;   //!
   TBranch        *b_t_DTTopMu1;   //!
   TBranch        *b_t_DTTopMu2;   //!


void Init(TTree *fChain)
{
   
   t_PFJetPt = 0;
   t_PFJetEta = 0;
   t_PFJetPhi = 0;
   t_PFJetE = 0;
   t_PFJetBTag = 0;
   t_PFMuonPt = 0;
   t_PFMuonEta = 0;
   t_PFMuonPhi = 0;
   t_PFMuonE = 0;
   t_PFMuonCh = 0;
   t_PFMuonChiSq = 0;
   t_PFMuonValidMuonHits = 0;
   t_PFMuonMatchedStations = 0;
   t_PFMuonValidPixelHits = 0;
   t_PFMuonTrackerLayerMeasured = 0;
   t_PFMuonCharHadEt = 0;
   t_PFMuonNeutHadEt = 0;
   t_PFMuonPhotEt = 0;
   t_PFMuonSumPUPt = 0;
   t_PFMuonRelIso = 0;
   t_PFMuonCharHadEt03 = 0;
   t_PFMuonNeutHadEt03 = 0;
   t_PFMuonPhotEt03 = 0;
   t_PFMuonSumPUPt03 = 0;
   t_PFMuonRelIso03 = 0;
   t_PFMuonDirIso01 = 0;
   t_PFMuonDirIso02 = 0;
   t_PFMuonDirIso03 = 0;
   t_PFMuonDirIso04 = 0;
   t_PFMuonDirIso05 = 0;
   t_PFMuonID = 0;
   t_PFMuonVtxAss = 0;
   t_PFMuonIsolated = 0;
   t_PFMuonIDVtxIso = 0;
   t_PFMuonIsWdau = 0;
   t_PFMuonIsBdau = 0;
   t_PFMuonGenMuPt = 0;
   t_PFMuonGenMuEta = 0;
   t_PFMuonGenMuPhi = 0;
   t_PFElecPt = 0;
   t_PFElecEta = 0;
   t_PFElecPhi = 0;
   t_PFElecE = 0;
   t_PFElecCh = 0;
   t_PFEleIsEB = 0;
   t_PFEleIsEE = 0;
   t_PFEleDEtaIn = 0;
   t_PFEleDPhiIn = 0;
   t_PFEleSigIEtaIEta = 0;
   t_PFEleHOE = 0;
   t_PFEleD0Vtx = 0;
   t_PFEleZVtx = 0;
   t_PFEleIsoCH = 0;
   t_PFEleIsoEM = 0;
   t_PFEleIsoNH = 0;
   t_PFEleIsoRho = 0;
   t_PFEleIsoEffA = 0;
   t_PFEleRelIso = 0;
   t_PFEleID = 0;
   t_PFEleIsolated = 0;
   t_PFEleIDIso = 0;
   t_PFEleDirIso01 = 0;
   t_PFEleDirIso02 = 0;
   t_PFEleDirIso03 = 0;
   t_PFEleDirIso04 = 0;
   t_PFEleDirIso05 = 0;
   t_PFEleIsWdau = 0;
   t_PFEleIsBdau = 0;
   t_PFEleGenElePt = 0;
   t_PFEleGenEleEta = 0;
   t_PFEleGenElePhi = 0;
   t_genPartPdgId = 0;
   t_genPartStatus = 0;
   t_genPartPt = 0;
   t_genPartEta = 0;
   t_genPartPhi = 0;
   t_genPartE = 0;
   t_genPartDecayMode = 0;
   t_fatJetSubJetPt = 0;
   t_fatJetSubJetEta = 0;
   t_fatJetSubJetPhi = 0;
   t_fatJetSubJetE = 0;
  
   fChain->SetBranchAddress("t_EvtRun", &t_EvtRun, &b_t_EvtRun);
   fChain->SetBranchAddress("t_EvtLS", &t_EvtLS, &b_t_EvtLS);
   fChain->SetBranchAddress("t_EvtEvent", &t_EvtEvent, &b_t_EvtEvent);
   fChain->SetBranchAddress("t_NVertices", &t_NVertices, &b_t_NVertices);
   fChain->SetBranchAddress("t_PUWeight", &t_PUWeight, &b_t_PUWeight);
   fChain->SetBranchAddress("t_PFMetPx", &t_PFMetPx, &b_t_PFMetPx);
   fChain->SetBranchAddress("t_PFMetPy", &t_PFMetPy, &b_t_PFMetPy);
   fChain->SetBranchAddress("t_PFMetE", &t_PFMetE, &b_t_PFMetE);
   fChain->SetBranchAddress("t_PFMetPhi", &t_PFMetPhi, &b_t_PFMetPhi);
   fChain->SetBranchAddress("t_PFMetSignificance", &t_PFMetSignificance, &b_t_PFMetSignificance);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("t_PFJetPt", &t_PFJetPt, &b_t_PFJetPt);
   fChain->SetBranchAddress("t_PFJetEta", &t_PFJetEta, &b_t_PFJetEta);
   fChain->SetBranchAddress("t_PFJetPhi", &t_PFJetPhi, &b_t_PFJetPhi);
   fChain->SetBranchAddress("t_PFJetE", &t_PFJetE, &b_t_PFJetE);
   fChain->SetBranchAddress("t_PFJetBTag", &t_PFJetBTag, &b_t_PFJetBTag);
   fChain->SetBranchAddress("t_NJetsPt30Eta2p5", &t_NJetsPt30Eta2p5, &b_t_NJetsPt30Eta2p5);
   fChain->SetBranchAddress("t_NJetsPt30Eta5p0", &t_NJetsPt30Eta5p0, &b_t_NJetsPt30Eta5p0);
   fChain->SetBranchAddress("t_NJetsPt50Eta2p5", &t_NJetsPt50Eta2p5, &b_t_NJetsPt50Eta2p5);
   fChain->SetBranchAddress("t_NJetsPt50Eta5p0", &t_NJetsPt50Eta5p0, &b_t_NJetsPt50Eta5p0);
   fChain->SetBranchAddress("t_PFht", &t_PFht, &b_t_PFht);
   fChain->SetBranchAddress("t_PFmht", &t_PFmht, &b_t_PFmht);
   fChain->SetBranchAddress("t_PFmphi", &t_PFmphi, &b_t_PFmphi);
   fChain->SetBranchAddress("t_NVetoMuon", &t_NVetoMuon, &b_t_NVetoMuon);
   fChain->SetBranchAddress("t_NVetoEle", &t_NVetoEle, &b_t_NVetoEle);
   fChain->SetBranchAddress("t_PFMuonPt", &t_PFMuonPt, &b_t_PFMuonPt);
   fChain->SetBranchAddress("t_PFMuonEta", &t_PFMuonEta, &b_t_PFMuonEta);
   fChain->SetBranchAddress("t_PFMuonPhi", &t_PFMuonPhi, &b_t_PFMuonPhi);
   fChain->SetBranchAddress("t_PFMuonE", &t_PFMuonE, &b_t_PFMuonE);
   fChain->SetBranchAddress("t_PFMuonCh", &t_PFMuonCh, &b_t_PFMuonCh);
   fChain->SetBranchAddress("t_PFMuonChiSq", &t_PFMuonChiSq, &b_t_PFMuonChiSq);
   fChain->SetBranchAddress("t_PFMuonValidMuonHits", &t_PFMuonValidMuonHits, &b_t_PFMuonValidMuonHits);
   fChain->SetBranchAddress("t_PFMuonMatchedStations", &t_PFMuonMatchedStations, &b_t_PFMuonMatchedStations);
   fChain->SetBranchAddress("t_PFMuonValidPixelHits", &t_PFMuonValidPixelHits, &b_t_PFMuonValidPixelHits);
   fChain->SetBranchAddress("t_PFMuonTrackerLayerMeasured", &t_PFMuonTrackerLayerMeasured, &b_t_PFMuonTrackerLayerMeasured);
   fChain->SetBranchAddress("t_PFMuonCharHadEt", &t_PFMuonCharHadEt, &b_t_PFMuonCharHadEt);
   fChain->SetBranchAddress("t_PFMuonNeutHadEt", &t_PFMuonNeutHadEt, &b_t_PFMuonNeutHadEt);
   fChain->SetBranchAddress("t_PFMuonPhotEt", &t_PFMuonPhotEt, &b_t_PFMuonPhotEt);
   fChain->SetBranchAddress("t_PFMuonSumPUPt", &t_PFMuonSumPUPt, &b_t_PFMuonSumPUPt);
   fChain->SetBranchAddress("t_PFMuonRelIso", &t_PFMuonRelIso, &b_t_PFMuonRelIso);
   fChain->SetBranchAddress("t_PFMuonCharHadEt03", &t_PFMuonCharHadEt03, &b_t_PFMuonCharHadEt03);
   fChain->SetBranchAddress("t_PFMuonNeutHadEt03", &t_PFMuonNeutHadEt03, &b_t_PFMuonNeutHadEt03);
   fChain->SetBranchAddress("t_PFMuonPhotEt03", &t_PFMuonPhotEt03, &b_t_PFMuonPhotEt03);
   fChain->SetBranchAddress("t_PFMuonSumPUPt03", &t_PFMuonSumPUPt03, &b_t_PFMuonSumPUPt03);
   fChain->SetBranchAddress("t_PFMuonRelIso03", &t_PFMuonRelIso03, &b_t_PFMuonRelIso03);
   fChain->SetBranchAddress("t_PFMuonDirIso01", &t_PFMuonDirIso01, &b_t_PFMuonDirIso01);
   fChain->SetBranchAddress("t_PFMuonDirIso02", &t_PFMuonDirIso02, &b_t_PFMuonDirIso02);
   fChain->SetBranchAddress("t_PFMuonDirIso03", &t_PFMuonDirIso03, &b_t_PFMuonDirIso03);
   fChain->SetBranchAddress("t_PFMuonDirIso04", &t_PFMuonDirIso04, &b_t_PFMuonDirIso04);
   fChain->SetBranchAddress("t_PFMuonDirIso05", &t_PFMuonDirIso05, &b_t_PFMuonDirIso05);
   fChain->SetBranchAddress("t_PFMuonID", &t_PFMuonID, &b_t_PFMuonID);
   fChain->SetBranchAddress("t_PFMuonVtxAss", &t_PFMuonVtxAss, &b_t_PFMuonVtxAss);
   fChain->SetBranchAddress("t_PFMuonIsolated", &t_PFMuonIsolated, &b_t_PFMuonIsolated);
   fChain->SetBranchAddress("t_PFMuonIDVtxIso", &t_PFMuonIDVtxIso, &b_t_PFMuonIDVtxIso);
   fChain->SetBranchAddress("t_PFMuonIsWdau", &t_PFMuonIsWdau, &b_t_PFMuonIsWdau);
   fChain->SetBranchAddress("t_PFMuonIsBdau", &t_PFMuonIsBdau, &b_t_PFMuonIsBdau);
   fChain->SetBranchAddress("t_PFMuonGenMuPt", &t_PFMuonGenMuPt, &b_t_PFMuonGenMuPt);
   fChain->SetBranchAddress("t_PFMuonGenMuEta", &t_PFMuonGenMuEta, &b_t_PFMuonGenMuEta);
   fChain->SetBranchAddress("t_PFMuonGenMuPhi", &t_PFMuonGenMuPhi, &b_t_PFMuonGenMuPhi);
   fChain->SetBranchAddress("t_PFElecPt", &t_PFElecPt, &b_t_PFElecPt);
   fChain->SetBranchAddress("t_PFElecEta", &t_PFElecEta, &b_t_PFElecEta);
   fChain->SetBranchAddress("t_PFElecPhi", &t_PFElecPhi, &b_t_PFElecPhi);
   fChain->SetBranchAddress("t_PFElecE", &t_PFElecE, &b_t_PFElecE);
   fChain->SetBranchAddress("t_PFElecCh", &t_PFElecCh, &b_t_PFElecCh);
   fChain->SetBranchAddress("t_PFEleIsEB", &t_PFEleIsEB, &b_t_PFEleIsEB);
   fChain->SetBranchAddress("t_PFEleIsEE", &t_PFEleIsEE, &b_t_PFEleIsEE);
   fChain->SetBranchAddress("t_PFEleDEtaIn", &t_PFEleDEtaIn, &b_t_PFEleDEtaIn);
   fChain->SetBranchAddress("t_PFEleDPhiIn", &t_PFEleDPhiIn, &b_t_PFEleDPhiIn);
   fChain->SetBranchAddress("t_PFEleSigIEtaIEta", &t_PFEleSigIEtaIEta, &b_t_PFEleSigIEtaIEta);
   fChain->SetBranchAddress("t_PFEleHOE", &t_PFEleHOE, &b_t_PFEleHOE);
   fChain->SetBranchAddress("t_PFEleD0Vtx", &t_PFEleD0Vtx, &b_t_PFEleD0Vtx);
   fChain->SetBranchAddress("t_PFEleZVtx", &t_PFEleZVtx, &b_t_PFEleZVtx);
   fChain->SetBranchAddress("t_PFEleIsoCH", &t_PFEleIsoCH, &b_t_PFEleIsoCH);
   fChain->SetBranchAddress("t_PFEleIsoEM", &t_PFEleIsoEM, &b_t_PFEleIsoEM);
   fChain->SetBranchAddress("t_PFEleIsoNH", &t_PFEleIsoNH, &b_t_PFEleIsoNH);
   fChain->SetBranchAddress("t_PFEleIsoRho", &t_PFEleIsoRho, &b_t_PFEleIsoRho);
   fChain->SetBranchAddress("t_PFEleIsoEffA", &t_PFEleIsoEffA, &b_t_PFEleIsoEffA);
   fChain->SetBranchAddress("t_PFEleRelIso", &t_PFEleRelIso, &b_t_PFEleRelIso);
   fChain->SetBranchAddress("t_PFEleID", &t_PFEleID, &b_t_PFEleID);
   fChain->SetBranchAddress("t_PFEleIsolated", &t_PFEleIsolated, &b_t_PFEleIsolated);
   fChain->SetBranchAddress("t_PFEleIDIso", &t_PFEleIDIso, &b_t_PFEleIDIso);
   fChain->SetBranchAddress("t_PFEleDirIso01", &t_PFEleDirIso01, &b_t_PFEleDirIso01);
   fChain->SetBranchAddress("t_PFEleDirIso02", &t_PFEleDirIso02, &b_t_PFEleDirIso02);
   fChain->SetBranchAddress("t_PFEleDirIso03", &t_PFEleDirIso03, &b_t_PFEleDirIso03);
   fChain->SetBranchAddress("t_PFEleDirIso04", &t_PFEleDirIso04, &b_t_PFEleDirIso04);
   fChain->SetBranchAddress("t_PFEleDirIso05", &t_PFEleDirIso05, &b_t_PFEleDirIso05);
   fChain->SetBranchAddress("t_PFEleIsWdau", &t_PFEleIsWdau, &b_t_PFEleIsWdau);
   fChain->SetBranchAddress("t_PFEleIsBdau", &t_PFEleIsBdau, &b_t_PFEleIsBdau);
   fChain->SetBranchAddress("t_PFEleGenElePt", &t_PFEleGenElePt, &b_t_PFEleGenElePt);
   fChain->SetBranchAddress("t_PFEleGenEleEta", &t_PFEleGenEleEta, &b_t_PFEleGenEleEta);
   fChain->SetBranchAddress("t_PFEleGenElePhi", &t_PFEleGenElePhi, &b_t_PFEleGenElePhi);
   fChain->SetBranchAddress("t_genPartPdgId", &t_genPartPdgId, &b_t_genPartPdgId);
   fChain->SetBranchAddress("t_genPartStatus", &t_genPartStatus, &b_t_genPartStatus);
   fChain->SetBranchAddress("t_genPartPt", &t_genPartPt, &b_t_genPartPt);
   fChain->SetBranchAddress("t_genPartEta", &t_genPartEta, &b_t_genPartEta);
   fChain->SetBranchAddress("t_genPartPhi", &t_genPartPhi, &b_t_genPartPhi);
   fChain->SetBranchAddress("t_genPartE", &t_genPartE, &b_t_genPartE);
   fChain->SetBranchAddress("t_genPartDecayMode", &t_genPartDecayMode, &b_t_genPartDecayMode);
   fChain->SetBranchAddress("t_fatJet1_m123", &t_fatJet1_m123, &b_fatJet1_m123);
   fChain->SetBranchAddress("t_fatJet1_pt123", &t_fatJet1_pt123, &b_fatJet1_pt123);
   fChain->SetBranchAddress("t_fatJet1_eta123", &t_fatJet1_eta123, &b_fatJet1_eta123);
   fChain->SetBranchAddress("t_fatJet1_phi123", &t_fatJet1_phi123, &b_fatJet1_phi123);
   fChain->SetBranchAddress("t_fatJet1_m12", &t_fatJet1_m12, &b_fatJet1_m12);
   fChain->SetBranchAddress("t_fatJet1_m13", &t_fatJet1_m13, &b_fatJet1_m13);
   fChain->SetBranchAddress("t_fatJet1_m23", &t_fatJet1_m23, &b_fatJet1_m23);
   fChain->SetBranchAddress("t_fatJet2_m123", &t_fatJet2_m123, &b_fatJet2_m123);
   fChain->SetBranchAddress("t_fatJet2_pt123", &t_fatJet2_pt123, &b_fatJet2_pt123);
   fChain->SetBranchAddress("t_fatJet2_eta123", &t_fatJet2_eta123, &b_fatJet2_eta123);
   fChain->SetBranchAddress("t_fatJet2_phi123", &t_fatJet2_phi123, &b_fatJet2_phi123);
   fChain->SetBranchAddress("t_fatJet2_m12", &t_fatJet2_m12, &b_fatJet2_m12);
   fChain->SetBranchAddress("t_fatJet2_m13", &t_fatJet2_m13, &b_fatJet2_m13);
   fChain->SetBranchAddress("t_fatJet2_m23", &t_fatJet2_m23, &b_fatJet2_m23);
   fChain->SetBranchAddress("t_nFatJet", &t_nFatJet, &b_t_nFatJet);
   fChain->SetBranchAddress("t_fatJet1_cond1", &t_fatJet1_cond1, &b_t_fatJet1_cond1);
   fChain->SetBranchAddress("t_fatJet1_cond2", &t_fatJet1_cond2, &b_t_fatJet1_cond2);
   fChain->SetBranchAddress("t_fatJet1_cond3", &t_fatJet1_cond3, &b_t_fatJet1_cond3);
   fChain->SetBranchAddress("t_fatJet2_cond1", &t_fatJet2_cond1, &b_t_fatJet2_cond1);
   fChain->SetBranchAddress("t_fatJet2_cond2", &t_fatJet2_cond2, &b_t_fatJet2_cond2);
   fChain->SetBranchAddress("t_fatJet2_cond3", &t_fatJet2_cond3, &b_t_fatJet2_cond3);
   fChain->SetBranchAddress("t_fatJetSubJetPt", &t_fatJetSubJetPt, &b_t_fatJetSubJetPt);
   fChain->SetBranchAddress("t_fatJetSubJetEta", &t_fatJetSubJetEta, &b_t_fatJetSubJetEta);
   fChain->SetBranchAddress("t_fatJetSubJetPhi", &t_fatJetSubJetPhi, &b_t_fatJetSubJetPhi);
   fChain->SetBranchAddress("t_fatJetSubJetE", &t_fatJetSubJetE, &b_t_fatJetSubJetE);
   fChain->SetBranchAddress("t_DTWMass1", &t_DTWMass1, &b_t_DTWMass1);
   fChain->SetBranchAddress("t_DTWMass2", &t_DTWMass2, &b_t_DTWMass2);
   fChain->SetBranchAddress("t_DTMassChiSqr1", &t_DTMassChiSqr1, &b_t_DTMassChiSqr1);
   fChain->SetBranchAddress("t_DTMassChiSqr2", &t_DTMassChiSqr2, &b_t_DTMassChiSqr2);
   fChain->SetBranchAddress("t_DTTopMass1", &t_DTTopMass1, &b_t_DTTopMass1);
   fChain->SetBranchAddress("t_DTTopMass2", &t_DTTopMass2, &b_t_DTTopMass2);
   fChain->SetBranchAddress("t_DTTotalSystemMass", &t_DTTotalSystemMass, &b_t_DTTotalSystemMass);
   fChain->SetBranchAddress("t_DTTopPt1", &t_DTTopPt1, &b_t_DTTopPt1);
   fChain->SetBranchAddress("t_DTTopPt2", &t_DTTopPt2, &b_t_DTTopPt2);
   fChain->SetBranchAddress("t_DTTotalSystemPt", &t_DTTotalSystemPt, &b_t_DTTotalSystemPt);
   fChain->SetBranchAddress("t_DTTopEta1", &t_DTTopEta1, &b_t_DTTopEta1);
   fChain->SetBranchAddress("t_DTTopEta2", &t_DTTopEta2, &b_t_DTTopEta2);
   fChain->SetBranchAddress("t_DTTotalSystemEta", &t_DTTotalSystemEta, &b_t_DTTotalSystemEta);
   fChain->SetBranchAddress("t_DTTopPhi1", &t_DTTopPhi1, &b_t_DTTopPhi1);
   fChain->SetBranchAddress("t_DTTopPhi2", &t_DTTopPhi2, &b_t_DTTopPhi2);
   fChain->SetBranchAddress("t_DTTotalSystemPhi", &t_DTTotalSystemPhi, &b_t_DTTotalSystemPhi);
   fChain->SetBranchAddress("t_DTindexB1", &t_DTindexB1, &b_t_DTindexB1);
   fChain->SetBranchAddress("t_DTindexB2", &t_DTindexB2, &b_t_DTindexB2);
   fChain->SetBranchAddress("t_DTindexW11", &t_DTindexW11, &b_t_DTindexW11);
   fChain->SetBranchAddress("t_DTindexW12", &t_DTindexW12, &b_t_DTindexW12);
   fChain->SetBranchAddress("t_DTindexW21", &t_DTindexW21, &b_t_DTindexW21);
   fChain->SetBranchAddress("t_DTindexW22", &t_DTindexW22, &b_t_DTindexW22);
   fChain->SetBranchAddress("t_DTdPhi1", &t_DTdPhi1, &b_t_DTdPhi1);
   fChain->SetBranchAddress("t_DTdPhi2", &t_DTdPhi2, &b_t_DTdPhi2);
   fChain->SetBranchAddress("t_DTTopMu1", &t_DTTopMu1, &b_t_DTTopMu1);
   fChain->SetBranchAddress("t_DTTopMu2", &t_DTTopMu2, &b_t_DTTopMu2);
}


#include "FileHeader.h"
#include "tree.h"

Char_t cTemp[500];
Char_t cTemp2[300];
Char_t cTemp3[300];

TH1F *hDTWDist[2][2][10];
TH1F *hHTWDist[2][2][10];


Double_t iAllEvents = 0;
Double_t iEventsHT = 0;
Double_t iEventsDT =0;


Char_t *cHOrder[2]={"Leading","SubLeading"};
Char_t *cHVar[]={"Mass","PT","Eta","Phi","Met","MetSignificance","BG1","MVA","BG2","BG3"};
Char_t *cBGSel[2]={"Signal","BG"};

Float_t fHVar1Min[]={0,0,-3.0,-3.5,0,0,0,-1,0,0};
Float_t fHVar1Max[]={1000,500,3.0,3.5,500,1000,1000,1,200,1000};
Int_t iHVar1Bins[]={25,25,15,20,25,50,50,40,40,50};


Char_t *cHTagger[] = {"Combinotorial Tagger","HepMC Tagger"};
Char_t *cHTaggerS[] = {"DT","HT"};

void MakeIndividualPlots(Int_t iS){
  cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  "<<iS<<endl;

  iAllEvents = 0;
  iEventsHT = 0;
  iEventsDT =0;

  
  TMVA::Reader *reader = new TMVA::Reader();
  


  Float_t f_DTWMass1,f_DTMassChiSqr1,f_DTMassChiSqr2,f_DTTopEta1,f_DTdPhi1,f_DTTopMu1,fNJets,f_PFMetSignificance,f_PFMetE;
  reader->AddVariable( "t_DTWMass1", &f_DTWMass1 );
  reader->AddVariable( "t_DTMassChiSqr1",           &f_DTMassChiSqr1);
  reader->AddVariable( "t_DTMassChiSqr2",           &f_DTMassChiSqr2);
  reader->AddVariable( "t_DTTopEta1",               &f_DTTopEta1 );
  reader->AddVariable( "t_DTdPhi1",                 &f_DTdPhi1);
  reader->AddVariable( "t_DTTopMu1",                &f_DTTopMu1 );
  reader->AddVariable( "NJets",                &fNJets);
  reader->AddVariable( "t_PFMetSignificance",                &f_PFMetSignificance);
  reader->AddVariable( "t_PFMetE",                &f_PFMetE);
  
 
  reader->AddSpectator( "t_fatJet1_m123",                &(Float_t)t_fatJet1_m123);
  reader->AddSpectator( "t_fatJet1_m12",                &(Float_t)t_fatJet1_m12);
  reader->AddSpectator( "t_fatJet1_m13",                &(Float_t)t_fatJet1_m13);
  reader->AddSpectator( "t_fatJet1_m23",                &(Float_t)t_fatJet1_m23);
  reader->AddSpectator( "t_fatJet2_m123",                &(Float_t)t_fatJet2_m123);
  reader->AddSpectator( "t_fatJet2_m12",                &(Float_t)t_fatJet2_m12);
  reader->AddSpectator( "t_fatJet2_m13",                &(Float_t)t_fatJet2_m13);
  reader->AddSpectator( "t_fatJet2_m23",                &(Float_t)t_fatJet2_m23);
   


  reader->BookMVA("BDT","/uscms/home/vasurang/nobackup/Run2012/SUSYSTop/MVA/TMVA-v4.1.2/test/weights/TMVAClassification_BDT.weights.xml");


  
  for(Int_t i=0;i<2;i++){
    for(Int_t j=0;j<9;j++){
      for(Int_t l=0;l<2;l++){
	sprintf(cTemp,"WDist_%s_%s_%s",cBGSel[l],cHOrder[i],cHVar[j]);
	hDTWDist[l][i][j]=new TH1F(cTemp,"",iHVar1Bins[j],fHVar1Min[j],fHVar1Max[j]);
	
	sprintf(cTemp,"TDist_%s_%s_%s",cBGSel[l],cHOrder[i],cHVar[j]);
	hHTWDist[l][i][j]=new TH1F(cTemp,"",iHVar1Bins[j],fHVar1Min[j],fHVar1Max[j]);
      }
    }
  }

  

  for(Int_t iF=0;iF<iFileLen[iS];iF++){
      sprintf(cTemp,"%s/%s",cFileDir[iS],cFileNames[iS][iF]);
      cout<<cTemp<<endl;
      
      TFile *file = new TFile(cTemp);
      TTree *EventTree  = (TTree*)file->Get("lostLeptonTree/tree");
      Init(EventTree);
      cout<<"Processing :  "<<EventTree->GetEntries()<<endl;
      TH1F *hAllEventsT = (TH1F*)file->Get("allEvents");
      iAllEvents = 	iAllEvents+ (Double_t)hAllEventsT->Integral();

      for(Long_t l1=0;l1<EventTree->GetEntries();l1++){
      //      for(Long_t l1=0;l1<1000;l1++){
	   
	EventTree->GetEntry(l1);

	iAllEvents++;
	
	if(t_PFMetE<170)
	  continue;
	f_PFMetSignificance = t_PFMetSignificance ;
	if(f_PFMetSignificance<50)
	  continue;

	if(t_DTTotalSystemMass>0){

	   f_DTWMass1 = t_DTWMass1;
	  f_DTMassChiSqr1 =t_DTMassChiSqr1;
	  f_DTMassChiSqr2 = t_DTMassChiSqr2;
	  f_DTTopEta1  =f_DTTopEta1;
	  f_DTdPhi1 = t_DTdPhi1;
	  f_DTTopMu1 = 	t_DTTopMu1;
	  fNJets = (Float_t)NJets;
	  f_PFMetSignificance = t_PFMetSignificance ;
	  f_PFMetE = t_PFMetE;
	  
	  
	  if(reader->EvaluateMVA("BDT")<0)
	    continue;

	  hDTWDist[0][0][7]->Fill(reader->EvaluateMVA("BDT"));
	  hDTWDist[0][1][7]->Fill(reader->EvaluateMVA("BDT"));
	  
	  
	  iEventsDT++;

	  hDTWDist[0][0][0]->Fill(t_DTTopMass1);
	  hDTWDist[0][1][0]->Fill(t_DTTopMass2);
	  
	  hDTWDist[0][0][1]->Fill(t_DTTopPt1);
	  hDTWDist[0][1][1]->Fill(t_DTTopPt2);
	  
	  hDTWDist[0][0][2]->Fill(t_DTTopEta1);
	  hDTWDist[0][1][2]->Fill(t_DTTopEta2);
	  
	  hDTWDist[0][0][3]->Fill(t_DTTopPhi1);
	  hDTWDist[0][1][3]->Fill(t_DTTopPhi2);
	  
	  hDTWDist[0][0][4]->Fill(t_PFMetE);
	  hDTWDist[0][1][4]->Fill(t_PFMetE);
	  
	  hDTWDist[0][0][5]->Fill(t_PFMetSignificance);
	  hDTWDist[0][1][5]->Fill(t_PFMetSignificance);
	  
	  
	  hDTWDist[0][0][6]->Fill(t_DTWMass1);
	  hDTWDist[0][1][6]->Fill(t_DTWMass2);
	  
	 
	  
	}
	
	Bool_t bPresent = kFALSE;
	
	
	for(Int_t iJets =0;iJets<NJets;iJets++){
	  if((*t_PFJetBTag)[iJets]>0.25){
	    if( deltaR((*t_PFJetEta)[iJets],(*t_PFJetPhi)[iJets],t_fatJet1_eta123,t_fatJet1_phi123)>0.5)
	      bPresent = kTRUE;
	    
	  }
	  
	}
	
	
	if(t_fatJet1_pt123>200&&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3&&bPresent){
	  
	  iEventsHT++;
	  hHTWDist[0][0][0]->Fill(t_fatJet1_m123);
	  hHTWDist[0][1][0]->Fill(t_fatJet2_m123);
	  
	  hHTWDist[0][0][1]->Fill(t_fatJet1_pt123);
	  hHTWDist[0][1][1]->Fill(t_fatJet2_pt123);
	  
	  hHTWDist[0][0][2]->Fill(t_fatJet1_eta123);
	  hHTWDist[0][1][2]->Fill(t_fatJet2_eta123);
	  
	  hHTWDist[0][0][3]->Fill(t_fatJet1_phi123);
	  hHTWDist[0][1][3]->Fill(t_fatJet1_phi123);
	  
	  hHTWDist[0][0][4]->Fill(t_PFMetE);
	  hHTWDist[0][1][4]->Fill(t_PFMetE);
	  
	  hHTWDist[0][0][5]->Fill(t_PFMetSignificance);
	  hHTWDist[0][1][5]->Fill(t_PFMetSignificance);
	  
	  
	  hHTWDist[0][0][6]->Fill(t_DTWMass1);
	  hHTWDist[0][1][6]->Fill(t_DTWMass2);
	  
	  hHTWDist[0][0][8]->Fill(TMath::Min(TMath::Min(t_fatJet1_m12,t_fatJet1_m13),t_fatJet1_m23));
	  hHTWDist[0][1][8]->Fill(TMath::Min(TMath::Min(t_fatJet2_m12,t_fatJet2_m13),t_fatJet2_m23));
	} 	   
      }
  }


  sprintf(cTemp,"SensitivityStudy/BasicDist_%s.root",cFileStr[iS]);
  TFile *fOutput = new TFile(cTemp,"recreate");
  
  
  for(Int_t i=0;i<2;i++){
    for(Int_t j=0;j<9;j++){
      for(Int_t l=0;l<2;l++){			
	
	hDTWDist[l][i][j]->Scale(15000.0*dCrossSection[iS]/iAllEvents );
	hHTWDist[l][i][j]->Scale(15000.0*dCrossSection[iS]/iAllEvents);
	
	

	hDTWDist[l][i][j]->Write();
	hHTWDist[l][i][j]->Write();
      }
    }
  }


  fOutput->Close();
}





double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return deta*deta + dphi*dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return TMath::Sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}


double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > TMath::Pi()) 
    result -= (2.*TMath::Pi());
  while (result <= -1.*TMath::Pi()) 
    result += 2.*TMath::Pi();
  return result;
}

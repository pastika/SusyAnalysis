#include "tree.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

Char_t *cFileNames[]={"outTest_T2tt_mStop550_mLSP250_NOCUTS.root","test/Background/outTest_TTJets_TuneZ2star_8TeV-madgraph-tauola_1.root","test/ZInv/temp.root"};
Char_t *cFileNamesLegend[]={"T2tt_mStop550_mLSP250","TJets_TuneZ2star_8TeV-madgraph-tauola","ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph"};
Char_t *cFileNamesH[]={"Signal","TTJets","ZInv"};

Int_t iMarker[] = {24,25,26,27,28,20};
Int_t iMarkerColor[] = {1,2,4,6,8,10,12};

Int_t iFiles = 3;

TH1F *hDTWDist[2][2][10][10];
TH1F *hHTWDist[2][2][10][10];


Double_t iAllEvents[]={0,0,0};
Double_t iEventsHT[]={0,0,0};
Double_t iEventsDT[]={0,0,0};


Char_t *cHOrder[2]={"Leading","SubLeading"};
Char_t *cHVar[]={"Mass","PT","Eta","Phi","Met","MetSignificance","BG1","MVA","BG2","BG3"};
Char_t *cBGSel[2]={"Signal","BG"};

Float_t fHVar1Min[]={0,0,-3.0,-3.5,0,0,0,-1,0,0};
Float_t fHVar1Max[]={1000,500,3.0,3.5,500,1000,1000,1,200,1000};
Int_t iHVar1Bins[]={25,25,15,20,25,50,50,40,40,50};


Char_t cTemp[300];


Char_t *cHTagger[] = {"Combinotorial Tagger","HepMC Tagger"};
Char_t *cHTaggerS[] = {"DT","HT"};

void CompareBG(){

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
      for(Int_t k=0;k<iFiles;k++){
	for(Int_t l=0;l<2;l++){
	  sprintf(cTemp,"WDist_%s_%s_%s_%s",cBGSel[l],cHOrder[i],cHVar[j],cFileNamesH[k]);
	  hDTWDist[l][i][j][k]=new TH1F(cTemp,"",iHVar1Bins[j],fHVar1Min[j],fHVar1Max[j]);

	  sprintf(cTemp,"TDist_%s_%s_%s_%s",cBGSel[l],cHOrder[i],cHVar[j],cFileNamesH[k]);
	  hHTWDist[l][i][j][k]=new TH1F(cTemp,"",iHVar1Bins[j],fHVar1Min[j],fHVar1Max[j]);
	}
      }
    }
  }

  for(Int_t iF = 0; iF<iFiles;iF++){
    TFile  *fPlot = new TFile(cFileNames[iF]);
    fPlot->ls();

    TTree *EventTree  = (TTree*)fPlot->Get("lostLeptonTree/tree");
    init(EventTree);
    //    EventTree->Print();
    cout<<"Processing :  "<<EventTree->GetEntries()<<endl;

    for(Long_t l1=0;l1<EventTree->GetEntries();l1++){
    //      for(Long_t l1=0;l1<10000;l1++){
      iAllEvents[iF]++;

    
      
      EventTree->GetEntry(l1);

     if(t_PFMetE<100)
	continue;
      if(t_DTTotalSystemMass>0){
	iEventsDT[iF]++;
	hDTWDist[0][0][0][iF]->Fill(t_DTTopMass1);
	hDTWDist[0][1][0][iF]->Fill(t_DTTopMass2);
	
	hDTWDist[0][0][1][iF]->Fill(t_DTTopPt1);
	hDTWDist[0][1][1][iF]->Fill(t_DTTopPt2);
	
	hDTWDist[0][0][2][iF]->Fill(t_DTTopEta1);
	hDTWDist[0][1][2][iF]->Fill(t_DTTopEta2);
	
	hDTWDist[0][0][3][iF]->Fill(t_DTTopPhi1);
	hDTWDist[0][1][3][iF]->Fill(t_DTTopPhi2);
	
	hDTWDist[0][0][4][iF]->Fill(t_PFMetE);
	hDTWDist[0][1][4][iF]->Fill(t_PFMetE);
	
	hDTWDist[0][0][5][iF]->Fill(t_PFMetSignificance);
	hDTWDist[0][1][5][iF]->Fill(t_PFMetSignificance);
	

	hDTWDist[0][0][6][iF]->Fill(t_DTWMass1);
	hDTWDist[0][1][6][iF]->Fill(t_DTWMass2);

	f_DTWMass1 = t_DTWMass1;
	f_DTMassChiSqr1 =t_DTMassChiSqr1;
	f_DTMassChiSqr2 = t_DTMassChiSqr2;
	f_DTTopEta1  =f_DTTopEta1;
	f_DTdPhi1 = t_DTdPhi1;
	f_DTTopMu1 = 	t_DTTopMu1;
	fNJets = (Float_t)NJets;
	f_PFMetSignificance = t_PFMetSignificance ;
	f_PFMetE = t_PFMetE;
	
     
	hDTWDist[0][0][7][iF]->Fill(reader->EvaluateMVA("BDT"));
	hDTWDist[0][1][7][iF]->Fill(reader->EvaluateMVA("BDT"));

      }
      
      
      Bool_t bPresent = kFALSE;


      for(Int_t iJets =0;iJets<NJets;iJets++){
	if((*t_PFJetBTag)[iJets]>0.25){
	  if( deltaR((*t_PFJetEta)[iJets],(*t_PFJetPhi)[iJets],t_fatJet1_eta123,t_fatJet1_phi123)>0.5)
	  bPresent = kTRUE;

	}
	
      }
      
      
      if(t_fatJet1_pt123>200&&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3&&bPresent){

	iEventsHT[iF]++;
	hHTWDist[0][0][0][iF]->Fill(t_fatJet1_m123);
	hHTWDist[0][1][0][iF]->Fill(t_fatJet2_m123);
	
	hHTWDist[0][0][1][iF]->Fill(t_fatJet1_pt123);
	hHTWDist[0][1][1][iF]->Fill(t_fatJet2_pt123);
	
	hHTWDist[0][0][2][iF]->Fill(t_fatJet1_eta123);
	hHTWDist[0][1][2][iF]->Fill(t_fatJet2_eta123);
	
	hHTWDist[0][0][3][iF]->Fill(t_fatJet1_phi123);
	hHTWDist[0][1][3][iF]->Fill(t_fatJet1_phi123);
	
	hHTWDist[0][0][4][iF]->Fill(t_PFMetE);
	hHTWDist[0][1][4][iF]->Fill(t_PFMetE);
	
	hHTWDist[0][0][5][iF]->Fill(t_PFMetSignificance);
	hHTWDist[0][1][5][iF]->Fill(t_PFMetSignificance);
	

	hHTWDist[0][0][6][iF]->Fill(t_DTWMass1);
	hHTWDist[0][1][6][iF]->Fill(t_DTWMass2);
	
	hHTWDist[0][0][8][iF]->Fill(TMath::Min(TMath::Min(t_fatJet1_m12,t_fatJet1_m13),t_fatJet1_m23));
	hHTWDist[0][1][8][iF]->Fill(TMath::Min(TMath::Min(t_fatJet2_m12,t_fatJet2_m13),t_fatJet2_m23));
      }
      
    }

    fPlot->Close();
  }


  cout<< iEventsHT[0]/iAllEvents[0]<<" : "<<iEventsHT[1]/iAllEvents[1]<<" : "<<iEventsHT[2]/iAllEvents[2]<<endl;
  cout<< iEventsDT[0]/iAllEvents[0]<<" : "<<iEventsDT[1]/iAllEvents[1]<<" : "<<iEventsDT[2]/iAllEvents[2]<<endl;
    
  cout<< iEventsHT[0]/TMath::Sqrt(iEventsHT[1])<<" : "<<iEventsHT[0]/TMath::Sqrt(iEventsHT[2])<<endl;
  cout<< iEventsDT[0]/TMath::Sqrt(iEventsDT[1])<<" : "<<endl;;
  cout<<iEventsDT[0]/TMath::Sqrt(iEventsDT[2])<<endl;


  for(Int_t i=0;i<2;i++){
    for(Int_t j=0;j<8;j++){
      if(j<8)
	plotOverlayDist(hDTWDist[0][i][j],j,0,i,0);
       if(j!=6&&j!=7)
	 plotOverlayDist(hHTWDist[0][i][j],j,1,i,0);
      
    }
  }

}



void plotOverlayDist(TH1F *hTemp[10],Int_t iVar,Int_t iType,Int_t iOrder,Int_t iSel){
    for(Int_t iF = 0; iF<iFiles;iF++){
      hTemp[iF]->Sumw2();
      hTemp[iF]->Scale(1.0/hTemp[iF]->GetMaximum());
      hTemp[iF]->SetMarkerStyle( iMarker[iF]);
      hTemp[iF]->SetMarkerColor( iMarkerColor[iF]);
      hTemp[iF]->SetLineColor( iMarkerColor[iF]);
      hTemp[iF]->GetXaxis()->SetTitle(cHVar[iVar]);
      hTemp[iF]->GetYaxis()->SetRangeUser(0.0,1.5);
      if(iF==0){
	hTemp[iF]->Draw("P");
      }else{
	hTemp[iF]->Draw("Same P");
      }
      
    }
    
    TLatex *la1 = new TLatex(); 
    la1->SetNDC();
    la1->SetTextFont(42);
    la1->SetTextSize(0.045);
    la1->DrawLatex(0.5,0.96,cHTagger[iType]);


    sprintf(cTemp,"BGCompareDist_%s_%s_%s_%s.gif",cBGSel[iSel],cHOrder[iOrder],cHTaggerS[iType],cHVar[iVar]);
    c1->Print(cTemp,"gif");

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

#include "tree.h"



Char_t *cString="T2tt_mStop550_mLSP250";
Char_t *cFileName="outTest_T2tt_mStop550_mLSP250_NOCUTS.root";

//Char_t *cString="TTJets_TuneZ2star_8TeV-madgraph-tauola";
//Char_t *cFileName="test/Background/outTest_TTJets_TuneZ2star_8TeV-madgraph-tauola.root";


TH1F *hWMass[2][3];
TH1F *hTopMass[5][2][3];
TH1F *hDiSTopSystem[4][3];
TH1F *hTopGlobal[3][3];
Char_t cTemp[200];

Char_t *cTempKin[]={"Mass","Pt","Eta","Phi","deltaPhiMetWrsBJet"};
Char_t *cTempGlobal[]={"MET","WMassChiSquared","TopMassChiSquared"};
Char_t *cTempType[]={"All","WithPtCut","WithHepTagPtCut"};
Char_t *cTempTop[]={"Leading Top","Sub-leading Top"};

Char_t *cTempLegend1[]={"All DT events","DT Lead Top Pt>200","HT With FatJet Pt>200"};
Char_t *cTempLegend2[]={"All events","All DT Events","HT With FatJet Pt>200"};

Int_t iMarkerType[] = {20,24,25};
Int_t iMarkerColor[] = {1,2,4};

Float_t fMinx =0;
Float_t fMaxx =1000;

ofstream DQHtml;

Char_t cTemp2[300];


void MakePlots( ){

  for(Int_t j=0;j<3;j++){
    sprintf(cTemp,"%s_Top%sDT%s",cString,cTempGlobal[0],cTempType[j]);
    cout<<cTemp<<endl;
    hTopGlobal[0][j]=new TH1F(cTemp,"",200,0,2000);

    sprintf(cTemp,"%s_Top%sDT%s",cString,cTempGlobal[1],cTempType[j]);
    cout<<cTemp<<endl;
    hTopGlobal[1][j]=new TH1F(cTemp,"",100,0,1000);

    sprintf(cTemp,"%s_Top%sDT%s",cString,cTempGlobal[2],cTempType[j]);
    cout<<cTemp<<endl;
    hTopGlobal[2][j]=new TH1F(cTemp,"",100,0,1000);

    sprintf(cTemp,"%s_DiStop%sDT%s",cString,cTempKin[0],cTempType[j]);
    cout<<cTemp<<endl;
    hDiSTopSystem[0][j]=new TH1F(cTemp,"",100,0,2000);
    
    sprintf(cTemp,"%s_DiStop%sDT%s%d",cString,cTempKin[1],cTempType[j]);
    cout<<cTemp<<endl;
    hDiSTopSystem[1][j]=new TH1F(cTemp,"",30,0,600);
    
    sprintf(cTemp,"%s_DiStop%sDT%s%d",cString,cTempKin[2],cTempType[j]);
    cout<<cTemp<<endl;
    hDiSTopSystem[2][j]=new TH1F(cTemp,"",50,-5,5);
    
    sprintf(cTemp,"%s_DiStop%sDT%s%d",cString,cTempKin[3],cTempType[j]);
    cout<<cTemp<<endl;
    hDiSTopSystem[3][j]=new TH1F(cTemp,"",40,-4,4);
    
    

  }

  for(Int_t i=0;i<2;i++){
    
    for(Int_t j=0;j<3;j++){

      cout<<i<<" : "<<j<<endl;
      sprintf(cTemp,"%s_WMassDT%s%d",cString,cTempType[j],i+1);
      cout<<cTemp<<endl;
      hWMass[i][j]=new TH1F(cTemp,"",200,0,2000);
      
      sprintf(cTemp,"%s_Top%sDT%s%d",cString,cTempKin[0],cTempType[j],i+1);
      cout<<cTemp<<endl;
      hTopMass[0][i][j]=new TH1F(cTemp,"",200,0,2000);
      
      sprintf(cTemp,"%s_Top%sDT%s%d",cString,cTempKin[1],cTempType[j],i+1);
      cout<<cTemp<<endl;
      hTopMass[1][i][j]=new TH1F(cTemp,"",200,0,2000);

      sprintf(cTemp,"%s_Top%sDT%s%d",cString,cTempKin[2],cTempType[j],i+1);
      cout<<cTemp<<endl;
      hTopMass[2][i][j]=new TH1F(cTemp,"",50,-5,5);
      
      sprintf(cTemp,"%s_Top%sDT%s%d",cString,cTempKin[3],cTempType[j],i+1);
      cout<<cTemp<<endl;
      hTopMass[3][i][j]=new TH1F(cTemp,"",40,-4,4);

      sprintf(cTemp,"%s_Top%sDT%s%d",cString,cTempKin[4],cTempType[j],i+1);
      cout<<cTemp<<endl;
      hTopMass[4][i][j]=new TH1F(cTemp,"",70,-7,7);
    }
  }
  
  gSystem->mkdir(cString);
  sprintf(cTemp,"%s/Plots.html",cString);
  DQHtml.open(cTemp);
  
  

   TFile  *fPlot = new TFile(cFileName);
   fPlot->ls();

   TTree *EventTree  = (TTree*)fPlot->Get("lostLeptonTree/tree");
   init(EventTree);
   EventTree->Print();
   cout<<"Processing :  "<<EventTree->GetEntries()<<endl;

   Double_t dEvents = 0;
   Double_t dEventsDT = 0;
   Double_t dEventsPt200DT = 0;
   Double_t dEventsPt200HT = 0;
   Double_t dEventsPt200HTAndDT = 0;
   Double_t dEventsHTAndDT = 0;

   

   for(Long_t l=0;l<EventTree->GetEntries();l++){
     // for(Long_t l=0;l<10000;l++){
     dEvents++;
     
     EventTree->GetEntry(l);
     

     hTopGlobal[0][0]->Fill(t_PFMetE);


     //     if(t_PFMetE)
     //         continue;

     //Make plots if Dan's has fund some thing in the events;
     if(t_DTTotalSystemMass>0){
       dEventsDT++;

       hDiSTopSystem[0][0]->Fill(t_DTTotalSystemMass);
       hDiSTopSystem[1][0]->Fill(t_DTTotalSystemPt);
       hDiSTopSystem[2][0]->Fill(t_DTTotalSystemEta);
       hDiSTopSystem[3][0]->Fill(t_DTTotalSystemPhi);

       hTopGlobal[0][1]->Fill(t_PFMetE);
    
       hTopGlobal[1][0]->Fill(TMath::Sqrt(t_DTMassChiSqr1));
       hTopGlobal[2][0]->Fill(TMath::Sqrt(t_DTMassChiSqr2));
       
       hWMass[0][0]->Fill(t_DTWMass1);
       hWMass[1][0]->Fill(t_DTWMass2);
       
       hTopMass[0][0][0]->Fill(t_DTTopMass1);
       hTopMass[0][1][0]->Fill(t_DTTopMass2);
       
       hTopMass[1][0][0]->Fill(t_DTTopPt1);
       hTopMass[1][1][0]->Fill(t_DTTopPt2);
       
       hTopMass[2][0][0]->Fill(t_DTTopEta1);
       hTopMass[2][1][0]->Fill(t_DTTopEta2);
       
       hTopMass[3][0][0]->Fill(t_DTTopPhi1);
       hTopMass[3][1][0]->Fill(t_DTTopPhi2);
       
       hTopMass[4][0][0]->Fill(t_DTdPhi1);
       hTopMass[4][1][0]->Fill(t_DTdPhi2);
       
       if(t_DTTopPt1>200){
	 
	 dEventsPt200DT++;

	 hDiSTopSystem[0][1]->Fill(t_DTTotalSystemMass);
	 hDiSTopSystem[1][1]->Fill(t_DTTotalSystemPt);
	 hDiSTopSystem[2][1]->Fill(t_DTTotalSystemEta);
	 hDiSTopSystem[3][1]->Fill(t_DTTotalSystemPhi);
	 
	 
	 hWMass[0][1]->Fill(t_DTWMass1);
	 hWMass[1][1]->Fill(t_DTWMass2);
	 
	 
	 hTopMass[0][0][1]->Fill(t_DTTopMass1);
	 hTopMass[0][1][1]->Fill(t_DTTopMass2);
	 
	 hTopMass[1][0][1]->Fill(t_DTTopPt1);
	 hTopMass[1][1][1]->Fill(t_DTTopPt2);
	 
	 hTopMass[2][0][1]->Fill(t_DTTopEta1);
	 hTopMass[2][1][1]->Fill(t_DTTopEta2);
	 
	 hTopMass[3][0][1]->Fill(t_DTTopPhi1);
	 hTopMass[3][1][1]->Fill(t_DTTopPhi2);
	 
	 hTopMass[4][0][1]->Fill(t_DTdPhi1);
	 hTopMass[4][1][1]->Fill(t_DTdPhi2);
	 
	 hTopGlobal[1][1]->Fill(TMath::Sqrt(t_DTMassChiSqr1));
	 hTopGlobal[2][1]->Fill(TMath::Sqrt(t_DTMassChiSqr2));
	 
	 
       }
       
       if(t_fatJet1_pt123>200&&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3&&t_DTTopPt1>200){
	 hWMass[0][2]->Fill(t_DTWMass1);
	 hWMass[1][2]->Fill(t_DTWMass2);
	 
	  
	 hTopMass[0][0][2]->Fill(t_DTTopMass1);
	 hTopMass[0][1][2]->Fill(t_DTTopMass2);
	 
	 hTopMass[1][0][2]->Fill(t_DTTopPt1);
	 hTopMass[1][1][2]->Fill(t_DTTopPt2);
	 
	 hTopMass[2][0][2]->Fill(t_DTTopEta1);
	 hTopMass[2][1][2]->Fill(t_DTTopEta2);
	 
	 hTopMass[3][0][2]->Fill(t_DTTopPhi1);
	 hTopMass[3][1][2]->Fill(t_DTTopPhi2);
	 
	 hTopMass[4][0][2]->Fill(t_DTdPhi1);
	 hTopMass[4][1][2]->Fill(t_DTdPhi2);
	 
	 hTopGlobal[1][2]->Fill(TMath::Sqrt(t_DTMassChiSqr1));
	 hTopGlobal[2][2]->Fill(TMath::Sqrt(t_DTMassChiSqr2));

	 hDiSTopSystem[0][2]->Fill(t_DTTotalSystemMass);
	 hDiSTopSystem[1][2]->Fill(t_DTTotalSystemPt);
	 hDiSTopSystem[2][2]->Fill(t_DTTotalSystemEta);
	 hDiSTopSystem[3][2]->Fill(t_DTTotalSystemPhi);

       }
     }
     
     if(t_fatJet1_pt123>200&&t_DTTotalSystemMass>0&&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3){
       dEventsHTAndDT++;
     }
     
     
     if(t_fatJet1_pt123>200&&  t_DTTopPt1>200 &&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3){
       dEventsPt200HTAndDT++;
     }
     
     if(t_fatJet1_pt123>200 &&t_fatJet1_m123>147.3&&t_fatJet1_m123<197.3){
       hTopGlobal[0][2]->Fill(t_PFMetE);
       dEventsPt200HT++;
     }
   }
 
   cout<<"Making Plots"<<endl;
  
   DQHtml << "<H2 style=\"background-color: #BBBBBB;\">" << endl;
   DQHtml <<cString<<endl;
   DQHtml << "</H2>" << endl;

   DQHtml << "<table cellspacing=\"\" cellpadding=\"\" border=\"\" align=\"center\" >" << endl;
   DQHtml << "" << endl;
   DQHtml << "" << endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Total Number of events : "<<dEvents<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Number of events that pass Dan's tagger requirements : "<<dEventsDT<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Number of events that pass Dan's tagger requirements and leading top pT > 200 : "<<dEventsPt200DT<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Number of events that pass HepTopTagger requirements : "<<dEventsPt200HT<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Number of events that pass Dan's tagger requirements and HepTopTagger requirements : "<<dEventsHTAndDT<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> Number of events that pass Dan's tagger requirements and leading top pT > 200  and HepTopTagger requirements : "<<dEventsPt200HTAndDT<<"</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;

   DQHtml << "</table>" << endl;
   
   DQHtml << "<table cellspacing=\"\" cellpadding=\"\" border=\"\" align=\"center\" >" << endl;
   DQHtml << "" << endl;
   DQHtml << "" << endl;
   DQHtml<<"<TR>"<< endl;
   DQHtml<<"<TD  align='center' colspan=2> In the plots DT stands for Dan's tagger and HT stands for HepTopTagger</TD>"<< endl;
   DQHtml<<"</TR>"<< endl;
   DQHtml << "</table>" << endl;

   fPlot->Close();

   
   DQHtml << "<H3 style=\"background-color: #BBBBBB;\">" << endl;
   DQHtml <<"Di-stop system"<<endl;
   DQHtml << "</H3>" << endl;
   
   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_DiTopMass_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger Di-Top - Mass ");
   makePlots(cTemp, hDiSTopSystem[0],0);
   
   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_DiTopPt_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger Di-Top - Pt ");
   makePlots(cTemp, hDiSTopSystem[1],1);
   
   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_DiTopEta_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger Di-Top - Eta ");
   makePlots(cTemp, hDiSTopSystem[2],2);
   
   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_DiTopPhi_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger Di-Top - Phi ");
   makePlots(cTemp, hDiSTopSystem[3],3);

   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_WMassChiSquared_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger W ChiSquared ");
   makePlots2(cTemp,  hTopGlobal[1],1);

   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_TopMassChiSquared_CompareWithHT");
   sprintf(cTemp2,"Dan's Tagger Top ChiSquared ");
   makePlots2(cTemp,  hTopGlobal[2], 2);
   
   DQHtml << "<H3 style=\"background-color: #BBBBBB;\">" << endl;
   DQHtml <<"Event variables"<<endl;
   DQHtml << "</H3>" << endl;
   

   fMinx =0;
   fMaxx =1000;
   sprintf(cTemp,"DanTagger_MET_CompareWithHT");
   sprintf(cTemp2,"MET ");
   makePlots3(cTemp,  hTopGlobal[0],0);

   for(Int_t i=0;i<2;i++){

     DQHtml << "<H3 style=\"background-color: #BBBBBB;\">" << endl;
     DQHtml <<cTempTop[i]<<endl;
     DQHtml << "</H3>" << endl;

     
     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_WMass%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger W - Mass ");
     makePlots(cTemp,hWMass[i]);
     
   
     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_TopMass%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger Top - Mass ");
     makePlots(cTemp, hTopMass[0][i],0);

     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_TopPt%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger Top - Pt ");
     makePlots(cTemp, hTopMass[1][i],1);

     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_TopEta%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger Top - Eta ");
     makePlots(cTemp, hTopMass[2][i],2);

     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_TopPhi%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger Top - Phi ");
     makePlots(cTemp, hTopMass[3][i],3);
   
     
     fMinx =0;
     fMaxx =1000;
     sprintf(cTemp,"DanTagger_DeltaPhiMetWithBJet%i_CompareWithHT",i);
     sprintf(cTemp2,"Dan's Tagger DeltaPhi Missing Et with respect to b-tagged jet");
     makePlots(cTemp, hTopMass[4][i],4);
     
   }



   DQHtml << "</BODY>" << endl;
   DQHtml << "</HTML>" << endl;
   DQHtml.close();

   sprintf(cTemp,"Plots_%s.root",cString);
   TFile *fOut = new TFile(cTemp,"recreate");
   
   for(Int_t i=0;i<2;i++){
     for(Int_t j=0;j<3;j++){
       hWMass[i][j]->Write();

       for(Int_t k=0;k<3;k++){
	   hTopMass[k][i][j]->Write();
       }
     }
   }

   fOut->Close();

   cout<<"All done"<<endl;
}



void makePlots(Char_t *cName,TH1F *hTemp[3],Int_t iKin=0){
  
  TLatex *la1 = new TLatex(); 
  la1->SetNDC();
  la1->SetTextFont(42);
  la1->SetTextSize(0.045);
  
  
  TLegend *legend =new TLegend(0.6,0.9,0.9,0.75);

  for(Int_t i=0;i<3;i++){
    hTemp[i]->SetMarkerColor(iMarkerColor[i]);
    hTemp[i]->SetLineColor(iMarkerColor[i]);
    hTemp[i]->SetMarkerStyle(iMarkerType[i]);
    hTemp[i]->Sumw2();

    if(i==0){
      hTemp[i]->Draw("P");
    }else{
      hTemp[i]->Draw("same P");
    }
    
    hTemp[i]->GetXaxis()->SetTitle(cTempKin[iKin]);
    
    TLegendEntry *leTempN= legend->AddEntry(hTemp[i],cTempLegend1[i],"P");
    leTempN->SetTextColor(iMarkerColor[i]);
    leTempN->SetTextFont(42);
    leTempN->SetTextSize(0.03);
    
  }

  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  Char_t cTemp[500];
  
  la1->DrawLatex(0.5,0.96,cString);


  DQHtml << "<table cellspacing=\"\" cellpadding=\"\" border=\"\" align=\"center\" >" << endl;
  DQHtml << "" << endl;
  DQHtml << "" << endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center' colspan=2>"<<cTemp2<<"</TD>"<< endl;
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;
  
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s.gif",cString,cName);
  c1->Print(cTemp,"gif");

  sprintf(cTemp,"%s.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;
  sprintf(cTemp,"%s/%s.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s.pdf",cName);

  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;
  sprintf(cTemp,"%s/%s.C",cString,cName);
  c1->Print(cTemp,"C");
  
  sprintf(cTemp,"%s.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;

  c1->SetLogy();

  DQHtml<<"</TD>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s_LogScale.gif",cString,cName);
  c1->Print(cTemp,"gif");
  sprintf(cTemp,"%s_LogScale.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;

  sprintf(cTemp,"%s/%s_LogScale.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s_LogScale.pdf",cString,cName);
  
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;

  
  sprintf(cTemp,"%s/%s_LogScale.C",cString,cName);
  c1->Print(cTemp,"C");
  sprintf(cTemp,"%s_LogScale.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;
  
  c1->SetLogy(kFALSE);

  DQHtml<<"</TR>"<< endl;
  DQHtml << "</table>" << endl;
  
}




void makePlots2(Char_t *cName,TH1F *hTemp[3],Int_t iKin=0){
  
  TLatex *la1 = new TLatex(); 
  la1->SetNDC();
  la1->SetTextFont(42);
  la1->SetTextSize(0.045);
  
  
  TLegend *legend =new TLegend(0.6,0.9,0.9,0.75);

  for(Int_t i=0;i<3;i++){
    hTemp[i]->SetMarkerColor(iMarkerColor[i]);
    hTemp[i]->SetLineColor(iMarkerColor[i]);
    hTemp[i]->SetMarkerStyle(iMarkerType[i]);
    hTemp[i]->Sumw2();

    if(i==0){
      hTemp[i]->Draw("P");
    }else{
      hTemp[i]->Draw("same P");
    }
    
    hTemp[i]->GetXaxis()->SetTitle(cTempGlobal[iKin]);
    
    TLegendEntry *leTempN= legend->AddEntry(hTemp[i],cTempLegend1[i],"P");
    leTempN->SetTextColor(iMarkerColor[i]);
    leTempN->SetTextFont(42);
    leTempN->SetTextSize(0.03);
    
  }

  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  Char_t cTemp[500];
  
  la1->DrawLatex(0.5,0.96,cString);


  DQHtml << "<table cellspacing=\"\" cellpadding=\"\" border=\"\" align=\"center\" >" << endl;
  DQHtml << "" << endl;
  DQHtml << "" << endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center' colspan=2>"<<cTemp2<<"</TD>"<< endl;
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;
  
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s.gif",cString,cName);
  c1->Print(cTemp,"gif");

  sprintf(cTemp,"%s.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;
  sprintf(cTemp,"%s/%s.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s.pdf",cName);

  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;
  sprintf(cTemp,"%s/%s.C",cString,cName);
  c1->Print(cTemp,"C");
  
  sprintf(cTemp,"%s.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;

  c1->SetLogy();

  DQHtml<<"</TD>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s_LogScale.gif",cString,cName);
  c1->Print(cTemp,"gif");
  sprintf(cTemp,"%s_LogScale.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;

  sprintf(cTemp,"%s/%s_LogScale.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s_LogScale.pdf",cString,cName);
  
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;

  
  sprintf(cTemp,"%s/%s_LogScale.C",cString,cName);
  c1->Print(cTemp,"C");
  sprintf(cTemp,"%s_LogScale.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;
  
  c1->SetLogy(kFALSE);

  DQHtml<<"</TR>"<< endl;
  DQHtml << "</table>" << endl;
  
}



void makePlots3(Char_t *cName,TH1F *hTemp[3],Int_t iKin=0){
  
  TLatex *la1 = new TLatex(); 
  la1->SetNDC();
  la1->SetTextFont(42);
  la1->SetTextSize(0.045);
  
  
  TLegend *legend =new TLegend(0.6,0.9,0.9,0.75);

  for(Int_t i=0;i<3;i++){
    hTemp[i]->SetMarkerColor(iMarkerColor[i]);
    hTemp[i]->SetLineColor(iMarkerColor[i]);
    hTemp[i]->SetMarkerStyle(iMarkerType[i]);
    hTemp[i]->Sumw2();

    if(i==0){
      hTemp[i]->Draw("P");
    }else{
      hTemp[i]->Draw("same P");
    }
    
    hTemp[i]->GetXaxis()->SetTitle(cTempGlobal[iKin]);
    
    TLegendEntry *leTempN= legend->AddEntry(hTemp[i],cTempLegend2[i],"P");
    leTempN->SetTextColor(iMarkerColor[i]);
    leTempN->SetTextFont(42);
    leTempN->SetTextSize(0.03);
    
  }

  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  Char_t cTemp[500];
  
  la1->DrawLatex(0.5,0.96,cString);


  DQHtml << "<table cellspacing=\"\" cellpadding=\"\" border=\"\" align=\"center\" >" << endl;
  DQHtml << "" << endl;
  DQHtml << "" << endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center' colspan=2>"<<cTemp2<<"</TD>"<< endl;
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;
  
  DQHtml<<"</TR>"<< endl;
  DQHtml<<"<TR>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s.gif",cString,cName);
  c1->Print(cTemp,"gif");

  sprintf(cTemp,"%s.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;
  sprintf(cTemp,"%s/%s.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s.pdf",cName);

  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;
  sprintf(cTemp,"%s/%s.C",cString,cName);
  c1->Print(cTemp,"C");
  
  sprintf(cTemp,"%s.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;

  c1->SetLogy();

  DQHtml<<"</TD>"<< endl;
  DQHtml<<"<TD  align='center'>"<< endl;

  sprintf(cTemp,"%s/%s_LogScale.gif",cString,cName);
  c1->Print(cTemp,"gif");
  sprintf(cTemp,"%s_LogScale.gif",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"<img align='center' id=\""<<cTemp<<"\" name=\""<<cTemp<<"\" src=\""<<cTemp<<"\" alt=\"Click\" width='100' height='100' > " <<endl;
  DQHtml<< "</A>" << endl;
  
  DQHtml<< "<BR>" << endl;

  sprintf(cTemp,"%s/%s_LogScale.pdf",cString,cName);
  c1->Print(cTemp,"pdf");
  sprintf(cTemp,"%s_LogScale.pdf",cString,cName);
  
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(pdf)" <<endl;
  DQHtml<< "</A>" << endl;

  
  sprintf(cTemp,"%s/%s_LogScale.C",cString,cName);
  c1->Print(cTemp,"C");
  sprintf(cTemp,"%s_LogScale.C",cName);
  DQHtml<< "<A href=\""<<cTemp<<"\">"<<endl;
  DQHtml <<"(.C)" <<endl;
  DQHtml<< "</A>" << endl;
  
  c1->SetLogy(kFALSE);

  DQHtml<<"</TR>"<< endl;
  DQHtml << "</table>" << endl;
  
}

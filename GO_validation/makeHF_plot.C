#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include <iomanip>
#include <utility>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TCut.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <cstdio>
#include <cstring> 
#include <fstream>
#include <TTreeReader.h>



using namespace std;


void makeHF_plot(){
  
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetOptStat(0);  
  gStyle->SetStatStyle(1001);                                                                        
  gStyle->SetOptFit(1110);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelFont(132,"xyz");
  gStyle->SetLabelSize(0.04,"xyz"); 
  gStyle->SetTitleFont(132,"xyz");  
  gStyle->SetTitleStyle(1001);                                                                                  
  gStyle->SetTitleSize(0.04,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz"); 
  gStyle->SetNdivisions(510,"xyz");
  gStyle->SetStatX(0.995);		
  gStyle->SetStatY(0.995); 

  

  
  

  TH1::SetDefaultSumw2(true);
  TH1D::SetDefaultSumw2();

  Int_t nbin=200;
  Int_t axisMin=0;
  Int_t axisMax_hf=10000;
  Int_t axisMax_zdc=2e6;


  TH1F*myHist1 = new TH1F("HiHF", "", nbin, axisMin, axisMax_hf);
  TH1F*myHist2 = new TH1F("HiZDC", "", 20000, axisMin, axisMax_zdc);
  TH2F*myHist3 = new TH2F("HiHF_vs_ZDC_bin200", "", 200, axisMin, axisMax_hf, 2000, axisMin, axisMax_zdc);
  TH2F*myHist4 = new TH2F("HiHF_vs_ZDC_bin1000", "", 1000, axisMin, axisMax_hf, 20000, axisMin, axisMax_zdc);
  //TH2F*myHist5 = new TH2F("HiHF_vs_ZDC_bin10000", "", 10000, axisMin, axisMax_hf, 200000, axisMin, axisMax_zdc);

  TH1F*myHist1_evtsel = new TH1F("HiHF_evtsel", "", nbin, axisMin, axisMax_hf);
  TH1F*myHist2_evtsel = new TH1F("HiZDC_evtsel", "", 20000, axisMin, axisMax_zdc);
  TH2F*myHist3_evtsel = new TH2F("HiHF_vs_ZDC_bin200_evtsel", "", 200, axisMin, axisMax_hf, 2000, axisMin, axisMax_zdc);
  TH2F*myHist4_evtsel = new TH2F("HiHF_vs_ZDC_bin1000_evtsel", "", 1000, axisMin, axisMax_hf, 20000, axisMin, axisMax_zdc);
  //TH2F*myHist5_evtsel = new TH2F("HiHF_vs_ZDC_bin10000_evtsel", "", 10000, axisMin, axisMax_hf, 200000, axisMin, axisMax_zdc);



  TH1F * hist_hibin = new TH1F("CentBin", "CentBin", 100, 0, 100);

  const int Ncent= 7;
  Int_t cent[Ncent+1]={0, 10, 20, 40, 60, 80, 100, 200};
  const char*label_cent[Ncent]={"cent05","cent510","cent1020","cent2030","cent3040","cent4050","cent50100"};

  /*TH1F*h[Ncent];
  for(int ihist=0; ihist<Ncent; ihist++){
    h[ihist]=new TH1F(Form("HF_%s", label_cent[ihist]), Form("HF_%s", label_cent[ihist]), nbin, axisMin, axisMax);
    }*/



 
  const TString input_file = "HIForest_RawPrime0_all.txt";
  TFile * fout = new TFile("Pileup_out_v2.root","RECREATE");

  TString Str;
  ifstream fpr(Form("%s",input_file.Data()), ios::in);
  if(!fpr.is_open()){
    cout << "List of input files not found!" << endl;
    return;
  }

  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(fpr, file_chain))
    {
      file_name_vector.push_back(file_chain);
    }
  TChain *t = new TChain("hiEvtAnalyzer/HiTree");
  TChain *tskimanalysis = new TChain("skimanalysis/HltTree");
  TChain *thltanalysis = new TChain("hltanalysis/HltTree");

  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      t->Add(*listIterator);
      tskimanalysis->Add(*listIterator);
      thltanalysis->Add(*listIterator);
    }

  t->AddFriend(tskimanalysis);
  t->AddFriend(thltanalysis);
  
      
  Int_t CentBin, hiBin, pclusterCompatibilityFilter, pprimaryVertexFilter, pphfCoincFilterPF3Th5, HLT_HIMinimumBiasHF1ANDZDC1nOR_v4;
  Float_t hiHF_pf, hiZDC;
  UInt_t run;
  t->SetBranchAddress("run", &run);      
  t->SetBranchAddress("hiHF_pf", &hiHF_pf);
  t->SetBranchAddress("hiZDC", &hiZDC);
  //t->SetBranchAddress("CentBin", &CentBin);
  t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  t->SetBranchAddress("pphfCoincFilterPF3Th5",&pphfCoincFilterPF3Th5);
  t->SetBranchAddress("HLT_HIMinimumBiasHF1ANDZDC1nOR_v4", &HLT_HIMinimumBiasHF1ANDZDC1nOR_v4);

  t->SetBranchStatus("*", 0);
  for (const auto& p : {"run", "hiHF_pf", "hiZDC", "hiBin", "pclusterCompatibilityFilter", "pprimaryVertexFilter", "pphfCoincFilterPF3Th5", "HLT_HIMinimumBiasHF1ANDZDC1nOR_v4" })
    //for (const auto& p : {"CentBin" })
    t->SetBranchStatus(p, 1);  

  Int_t nevent =t->GetEntries();
      
  for (int ievt =0; ievt <nevent ; ievt++){
    if(ievt%100000 == 0) cout<<"Processing data event: " << ievt << " / " << nevent << endl;

    t->GetEntry(ievt);

    //hist_hibin->Fill(CentBin);


    
    myHist1->Fill(hiHF_pf);
    myHist2->Fill(hiZDC);
    myHist3->Fill(hiHF_pf, hiZDC);
    myHist4->Fill(hiHF_pf, hiZDC);
    //myHist5->Fill(hiHF_pf, hiZDC);


    if ((HLT_HIMinimumBiasHF1ANDZDC1nOR_v4 ==1)&& (pprimaryVertexFilter ==1) && (pclusterCompatibilityFilter ==1) && (pphfCoincFilterPF3Th5 > 0)){
    myHist1_evtsel->Fill(hiHF_pf);
    myHist2_evtsel->Fill(hiZDC);
    myHist3_evtsel->Fill(hiHF_pf, hiZDC);
    myHist4_evtsel->Fill(hiHF_pf, hiZDC);
    //myHist5_evtsel->Fill(hiHF_pf, hiZDC);
    }

    
    /*for(int icent=0; icent<Ncent; icent++){
	if (hiBin > cent[icent] && hiBin<cent[icent+1]){
	  h[icent]->Fill(hiHF_pf);	  
	}
	}*/
    //}
    
  }//---event loop
  
  delete t;
  //inFile->Close();
  //delete inFile;
  
  
  Double_t maxBinContent[Ncent];
  Double_t BinCenter[Ncent];
  TLine * line[7];

/*TCanvas *c1 = new TCanvas("c1", "c1",339,84,931,692);
  gStyle->SetOptStat(0);
  c1->Range(-2104.596,-1.690284,9673.534,6.820934);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogy();
  c1->SetLeftMargin(0.1786868);
  c1->SetRightMargin(0.1420883);
  c1->SetTopMargin(0.0640625);
  c1->SetBottomMargin(0.1359375);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  int maxBin = -1;
  
  c1->cd();

 
  myHist1->SetLineColor(1);
  myHist1->SetLineWidth(2);
  myHist1->SetMarkerStyle(20);
  myHist1->GetXaxis()->SetTitle("HF E_{T} (GeV)");
  myHist1->GetXaxis()->SetRange(1,160);
  myHist1->GetXaxis()->CenterTitle(true);
  myHist1->GetXaxis()->SetLabelFont(132);
  myHist1->GetXaxis()->SetLabelOffset(0.01);
  myHist1->GetXaxis()->SetLabelSize(0.04);
  myHist1->GetXaxis()->SetTitleSize(0.05);
  myHist1->GetXaxis()->SetTitleOffset(1);
  myHist1->GetXaxis()->SetTitleFont(132);
  myHist1->GetYaxis()->SetTitle("Events");
  myHist1->GetYaxis()->CenterTitle(true);
  myHist1->GetYaxis()->SetLabelFont(132);
  myHist1->GetYaxis()->SetLabelOffset(0.01);
  myHist1->GetYaxis()->SetLabelSize(0.04);
  myHist1->GetYaxis()->SetTitleSize(0.05);
  myHist1->GetYaxis()->SetTitleFont(132);
  myHist1->Draw("E1");

  TLatex *   tex = new TLatex(5413.629,2562545,"PbPb #sqrt{S_{NN}} = 5.36 TeV");
  tex->SetTextFont(132);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(-88.74797,2267129,"CMS");
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(621.2362,2337626,"Internal");
  tex->SetTextFont(52);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(587.7885,14.78892,"HLT_HIMinimumBiasHF1ANDZDC1nOR_v4");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(602.9768,6.81819,"PrimaryVertexFilter");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(602.9768,3.143414,"ClusterCompatibilityFilter");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(618.1652,1.37631,"pfConicFilterPF3Th5");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(602.9768,29.68776,"Run = 1111");
  tex->SetTextColor(2);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(4112.978,1.776478,"0-5%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(3423.886,1.535551,"5-10%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(2468.828,1.566999,"10-20%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(1670.931,1.776478,"20-30%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(1126.911,1.821621,"30-40%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(752.1413,1.915378,"40-50%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(292.7463,1.867911,"50-100%");
  tex->SetTextColor(2);
  tex->SetTextFont(132);
  tex->SetTextSize(0.03);
  tex->SetTextAngle(90);
  tex->SetLineWidth(2);
  tex->Draw();  
  
  for(int icnt=0; icnt<Ncent; icnt++){
    for (int iBin = 1; iBin <= h[icnt]->GetNbinsX(); ++iBin) {
      if (h[icnt]->GetBinContent(iBin) > 0) {
	maxBin = iBin;
	break;
      }
    }
    //  std::cout<<"maxBin="<<maxBin<<std::endl;
    maxBinContent[icnt] = myHist1->GetBinContent(maxBin);  
    BinCenter[icnt] = myHist1->GetXaxis()->GetBinCenter(maxBin);
    
    line[icnt] = new TLine(BinCenter[icnt], 0.0,  BinCenter[icnt] , maxBinContent[icnt]);
    line[icnt]->Draw("SAME");
    line[icnt]->SetLineColor(2);
    line[icnt]->SetLineWidth(2);
    line[icnt]->SetLineStyle(2);
}
*/  

  
  //c->SaveAs("hist_HF.root");
  
  //std::cout<<"maxBinXValue ="<<maxBinXValue<<std::endl;
  
  fout->cd();
  
//c1->Write();
  myHist1->Write();
  myHist2->Write();
  myHist3->Write();
  myHist4->Write();
  //myHist5->Write();

  myHist1_evtsel->Write();
  myHist2_evtsel->Write();
  myHist3_evtsel->Write();
  myHist4_evtsel->Write();
  //myHist5_evtsel->Write();  
//for(int j =0; j<Ncent; j++){
//  h[j]->Write();
//}
  
}

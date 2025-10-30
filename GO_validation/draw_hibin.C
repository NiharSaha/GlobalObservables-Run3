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


void draw_hibin(){
  
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

  

  
  

  TH1F::SetDefaultSumw2(true);
  TH1F::SetDefaultSumw2();




  //TH1F * hist_hibin_dfinder = new TH1F("CentBin_dfinder", "CentBin_dfinder", 200, 0, 200);
  TH1F * hist_hibin_forest = new TH1F("CentBin_forest", "CentBin_forest", 200, 0, 200);
  //TH1F * hist_hiHF = new TH1F("hiHf_pf", "hiHf_pf", 200, 0, 10000);
  
  const TString input_file = "2024A_v1.txt";
  TFile * fout = new TFile("Draw_hiBin_2024A_Oct29.root","RECREATE");

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
  //TChain *t = new TChain("Dfinder/ntDkpi");
  //TChain *tskimanalysis = new TChain("skimanalysis/HltTree");
  TChain *t = new TChain("hiEvtAnalyzer/HiTree");

  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      //t->Add(*listIterator);
      t->Add(*listIterator);
      //thltanalysis->Add(*listIterator);
    }

  //t->AddFriend(tevt);
  //t->AddFriend(thltanalysis);
  
      
  Int_t CentBin, hiBin, pclusterCompatibilityFilter, pprimaryVertexFilter, pphfCoincFilterPF3Th5, HLT_HIMinimumBiasHF1ANDZDC1nOR_v4;
  Float_t hiHF_pf, hiZDC;
  UInt_t run;
  
  //t->SetBranchAddress("CentBin", &CentBin);
  t->SetBranchAddress("hiBin", &hiBin);
  //t->SetBranchAddress("hiHF_pf", &hiHF_pf);
  

  t->SetBranchStatus("*", 0);
  
  for (const auto& p : {"hiBin"})
    t->SetBranchStatus(p, 1);  

  Int_t nevent =t->GetEntries();
      
  for (int ievt =0; ievt <nevent ; ievt++){
    if(ievt%100000 == 0) cout<<"Processing data event: " << ievt << " / " << nevent << endl;

    t->GetEntry(ievt);

    //hist_hibin_dfinder->Fill(CentBin);
    hist_hibin_forest->Fill(hiBin);
    //hist_hiHF->Fill(hiHF_pf);

    
  }//---event loop
  
  delete t;
  
  

  
  fout->cd();
  //hist_hibin_dfinder->Write();
  hist_hibin_forest->Write();
  //hist_hiHF->Write();
  fout->Close();
  
}

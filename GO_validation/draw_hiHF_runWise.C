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


void draw_hiHF_runWise(){
  
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





  std::map<UInt_t, TH1F*> hiHF_hists_wSel;
  std::map<UInt_t, TH1F*> hiHF_hists_noSel;
  
  const TString input_file = "Forest_2024B_v1.txt";
  TFile * fout = new TFile("Draw_hiHF_evtsel_Runwise_oct22_new.root","RECREATE");

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
  TChain *thltanalysis = new TChain("hltanalysis/HltTree");
  TChain *tskimanalysis = new TChain("skimanalysis/HltTree");

  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      t->Add(*listIterator);
      thltanalysis->Add(*listIterator);
      tskimanalysis->Add(*listIterator);
    }

  t->AddFriend(thltanalysis);
  t->AddFriend(tskimanalysis);

  
      
  Int_t CentBin, hiBin, pclusterCompatibilityFilter, pprimaryVertexFilter, pphfCoincFilterPF3Th5, HLT_HIMinimumBiasHF1ANDZDC1nOR_v4;
  Float_t hiHF_pf, hiZDC;
  UInt_t run;

  std::map<std::string, int> varI;
  std::vector<std::string> hlt_trigger_names;
  bool HLT = true;
  
  if (HLT) {
    std::string pattern = "HLT_HIMinimumBiasHF1AND";
    TObjArray* all_branches = thltanalysis->GetListOfBranches();
    TIter next(all_branches); TObject* obj;
    while ((obj = next())) {
      std::string branch_name = obj->GetName();
      if (branch_name.rfind(pattern, 0) == 0 && branch_name.find("Prescale") == std::string::npos) {
        hlt_trigger_names.push_back(branch_name);
        std::cout << "  Found HLT branch for logic check: " << branch_name << std::endl;
      }
    }
  }


  
  t->SetBranchStatus("*", 0); // Start by deactivating everything for performance
  
  t->SetBranchStatus("run", 1);
  t->SetBranchStatus("hiHF_pf", 1);
  t->SetBranchStatus("hiBin", 1);
  t->SetBranchStatus("pprimaryVertexFilter", 1);
  t->SetBranchStatus("pclusterCompatibilityFilter", 1);
  t->SetBranchStatus("pphfCoincFilterPF3Th5", 1);


  if (HLT) {
    t->SetBranchStatus("HLT_HIMinimumBias*", 1);
    std::cout << "[INFO] Activated all branches matching 'HLT_HIMinimumBias*'" << std::endl;
  }
  
  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress("hiHF_pf", &hiHF_pf);
  t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  t->SetBranchAddress("pphfCoincFilterPF3Th5", &pphfCoincFilterPF3Th5);


  t->SetBranchStatus("HLT_HIMinimumBiasHF1AND*", 1);
  for (const auto& name : hlt_trigger_names) {
    t->SetBranchAddress(name.c_str(), &(varI[name]));
  }

  
  

  Int_t nevent =t->GetEntries();
      
  for (int ievt =0; ievt <nevent ; ievt++){
    if(ievt%100000 == 0) cout<<"Processing data event: " << ievt << " / " << nevent << endl;

    t->GetEntry(ievt);

  
    if (hiHF_hists_wSel.find(run) == hiHF_hists_wSel.end()) {

      std::cout << "Found new run: " << run << ". Creating new histogram." << std::endl;
      TString hist_name_wSel = Form("hiHF_pf_run_%u_wSel", run);
      TString hist_title_wSel = Form("hiHF_pf (w/ Selections) for Run %u;hiHF_pf;Entries", run);
      hiHF_hists_wSel[run] = new TH1F(hist_name_wSel, hist_title_wSel, 200, 0, 10000);

      TString hist_name_noSel = Form("hiHF_pf_run_%u_noSel", run);
      TString hist_title_noSel = Form("hiHF_pf (No Selections) for Run %u;hiHF_pf;Entries", run);
      hiHF_hists_noSel[run] = new TH1F(hist_name_noSel, hist_title_noSel, 200, 0, 10000);

    }


    hiHF_hists_noSel[run]->Fill(hiHF_pf);
	
    bool hlt_passed = !HLT; // Default to true if HLT is off
    if (HLT) {
      for (const auto& name : hlt_trigger_names) {
	if (varI.count(name) && varI.at(name) > 0) {
	  hlt_passed = true;
	  break; // Found a passing trigger, no need to check more                                                                                                            
	}
      }
    }
    
    bool pass = (hlt_passed && pprimaryVertexFilter >0 && pclusterCompatibilityFilter > 0 && pphfCoincFilterPF3Th5 > 0);

    if(pass){
      hiHF_hists_wSel[run]->Fill(hiHF_pf);
    }


    

    
    
  }//---event loop
  
  delete t;

  std::map<UInt_t, TH1F*> hiHF_hists_wSel_norm;
  for (auto const& [run_num, hist_raw] : hiHF_hists_wSel) {
    
    TString norm_name = Form("%s_norm", hist_raw->GetName());
    TH1F* hist_norm = (TH1F*)hist_raw->Clone(norm_name);
    hist_norm->SetTitle(Form("hiHF_pf Normalized Shape for Run %u;hiHF_pf;Normalized Entries", run_num));

    if(hist_norm->Integral() > 0) {
      hist_norm->Scale(1.0 / hist_norm->Integral());
    }
    
    hiHF_hists_wSel_norm[run_num] = hist_norm;
  }

  std::map<UInt_t, TH1F*> hiHF_hists_noSel_norm;
  for (auto const& [run_num, hist_raw_noSel] : hiHF_hists_noSel) {
    TString norm_name = Form("%s_norm", hist_raw_noSel->GetName());
    TH1F* hist_norm = (TH1F*)hist_raw_noSel->Clone(norm_name);
    hist_norm->SetTitle(Form("hiHF_pf (No Sel) Normalized Shape for Run %u;hiHF_pf;Normalized Entries", run_num));
 
    if(hist_norm->Integral() > 0) {
      hist_norm->Scale(1.0 / hist_norm->Integral());
    }
    hiHF_hists_noSel_norm[run_num] = hist_norm;
  }

  
  fout->cd();

  std::cout << "Writing raw (with selection) histograms to file..." << std::endl;
  for (auto const& [run_num, hist] : hiHF_hists_wSel) {
    hist->Write();
  }
  std::cout << "Writing normalized (with selection) histograms to file..." << std::endl;
  for (auto const& [run_num, hist] : hiHF_hists_wSel_norm) {
    hist->Write();
  }

  std::cout << "Writing raw (no selection) histograms to file..." << std::endl;
  for (auto const& [run_num, hist] : hiHF_hists_noSel) {
    hist->Write();
  }
  std::cout << "Writing normalized (no selection) histograms to file..." << std::endl;
  for (auto const& [run_num, hist] : hiHF_hists_noSel_norm) {
    hist->Write();
  }
  

  
  fout->Close();
  
}

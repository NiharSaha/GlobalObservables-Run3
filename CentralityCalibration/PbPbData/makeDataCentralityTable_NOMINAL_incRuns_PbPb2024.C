#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <numeric>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"


void makeDataCentralityTable_NOMINAL_incRuns_PbPb2024(
                                          const TString input_file = "HIForest_RawPrime0_all.txt",
					  // --- MODIFIED: HLT_trg is no longer used, pattern is hardcoded below ---
                                          // const char* HLT_trg = "HLT_HIMinimumBiasHF1ANDZDC1nOR_v4",
                                          const bool HLT = true,
                                          const char* variable = "hiHF_pf",
                                          // --- MODIFIED: RUN variable is no longer needed for inclusive analysis ---
                                          // const int RUN = 387973,
                                          const char*RawPrime = "RawPrime0",
                                          const char* date = "Oct1_incRun_woPU", // Updated date for new version
                                          const char* CoinFilter = "pphfCoincFilterPF3Th5",
                                          const double threshold = 100.0,
                                          const char* label = "Nominal",
                                          const size_t nbins = 200
                                          )
{
  // Constant parameters
  const auto mcXscale = 1.012;
  const auto threshold_Min = 1000.0;
  const auto thresholdMax = 4000.0;

  //Tag
  const char* testMessage = Form("HIPhysics%s_%s", RawPrime, date);
  const char* tag = "CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal" ;
  const std::string outputTag = Form("2024Run_HYDMC_xSF%0.2f_%s_Threshold%.0f_%s_Normalisation%.0f_%.0f_%s", mcXscale, CoinFilter, threshold, label, threshold_Min, thresholdMax, testMessage);


  // Process data
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


  TFile *outFile = new TFile(Form("CentralityTable_HFtowers200_DataPbPb_usingMC_%s.root", outputTag.c_str()),"recreate");
  TNtuple * nt = new TNtuple("nt","","value");


  TH1F*hfData1 = new TH1F("hfData1","hfData1", 200,0, 10000);
  // TH1F*hfData2 = new TH1F("hfData2","hfData2", 100,0, 10000);
  TH1F*hfMc1 = new TH1F("hfMc1","hf mc", 200, 0, 10000);
  TH1F*hfMc2 = new TH1F("hfMc2","hf mc", 200, 0, 10000);
  TH1F* hfCombined = new TH1F("hfCombined","hf_combined", 200,0, 10000);

  TEfficiency*dataEff1 = new TEfficiency("dataEff1", "dataEff1", 1000, 0, 1000);
  // TEfficiency*dataEff2 = new TEfficiency("dataEff2", "dataEff2", 1000, 0, 1000);
  TEfficiency*mcEff= new TEfficiency("mcEff", "mcEff", 1000, 0, 1000);

  TH2D* h_zdc_vs_hihf = new TH2D("h_zdc_vs_hihf", "ZDC vs hiHF_pf;E_{T}^{HF, PF};ZDC E_{sum}", 200, 0, 10000, 250, 0, 800e3);
  TH1D* h_hf_before_cut = new TH1D("h_hf_before_cut", "hiHF_pf (Before Pileup Cut);E_{T}^{HF, PF};Events", 200, 0, 10000);
  TH1D* h_hf_after_cut = new TH1D("h_hf_after_cut", "hiHF_pf (After Pileup Cut);E_{T}^{HF, PF};Events", 200, 0, 10000);

  const int runNum = 1;
  CentralityBins*bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);



  std::vector<std::string> hlt_trigger_names;
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

  // --- Step 2: Activate all necessary branches for reading ---
  t->SetBranchStatus("*", 0); // Start by deactivating everything for performance

  // Activate physics and filter branches by name
  t->SetBranchStatus("run", 1);
  t->SetBranchStatus(variable, 1);
  t->SetBranchStatus("vz", 1);
  t->SetBranchStatus("hiZDC", 1);
  t->SetBranchStatus("hiBin", 1);
  t->SetBranchStatus("pprimaryVertexFilter", 1);
  t->SetBranchStatus("pclusterCompatibilityFilter", 1);
  t->SetBranchStatus(CoinFilter, 1);

  // --- MODIFIED: Activate all HLT branches using the wildcard, as you suggested ---
  if (HLT) {
    t->SetBranchStatus("HLT_HIMinimumBias*", 1);
    std::cout << "[INFO] Activated all branches matching 'HLT_HIMinimumBias*'" << std::endl;
  }


  // --- Step 3: Set branch addresses for all activated branches ---
  std::map<std::string, int> varI;
  std::map<std::string, float> varF;
  UInt_t run;

  // Set addresses for variables we will read
  t->SetBranchAddress("run", &run);
  t->SetBranchAddress(variable, &(varF[variable]));
  t->SetBranchAddress("vz", &(varF["vz"]));
  t->SetBranchAddress("hiZDC", &(varF["hiZDC"]));
  t->SetBranchAddress("hiBin", &(varI["hiBin"]));
  t->SetBranchAddress("pprimaryVertexFilter", &(varI["pprimaryVertexFilter"]));
  t->SetBranchAddress("pclusterCompatibilityFilter", &(varI["pclusterCompatibilityFilter"]));
  t->SetBranchAddress(CoinFilter, &(varI[CoinFilter]));

  // Set addresses for all the HLT trigger branches we found (excluding Prescale)
  for (const auto& name : hlt_trigger_names) {
      t->SetBranchAddress(name.c_str(), &(varI[name]));
  }
  

  
  std::vector<std::pair<float, bool>> values, hfdata;
  std::vector<std::pair<int, bool>> hibin;
  TH1D::SetDefaultSumw2();
  
  auto Nevents = t->GetEntries();
  double mcYscale_data(0);

  std::cout<<"Total Number of events = "<<Nevents<<std::endl;

  // --- Pileup Cut Parameters ---
  const double zdc_cut_intercept = 750e3;
  const double hihf_cut_intercept = 9000;
  const double pileup_slope = -zdc_cut_intercept / hihf_cut_intercept;

  for(Long64_t iev = 0; iev < Nevents; iev++) {

    if(iev%100000 == 0) cout<<"Processing data event: " << iev << " / " << Nevents << endl;

    t->GetEntry(iev);
    const float current_hihf = varF.at("hiHF_pf");
    const float current_zdc = varF.at("hiZDC");

    h_zdc_vs_hihf->Fill(current_hihf, current_zdc);

    const double y_cut_line = pileup_slope * current_hihf + zdc_cut_intercept;
    const bool pass_pileup_cut = (current_zdc <= y_cut_line);

    // --- MODIFIED: Manually loop through triggers to check for a pass ---
    bool hlt_passed = !HLT; // Default to true if HLT is off
    if (HLT) {
      for (const auto& name : hlt_trigger_names) {
	if (varI.count(name) && varI.at(name) > 0) {
	  hlt_passed = true;
	  break; // Found a passing trigger, no need to check more
	}
      }
    }
    
        
        
    bool pass_original_filters = (hlt_passed &&
                                  varI.at("pprimaryVertexFilter") > 0 &&
                                  varI.at("pclusterCompatibilityFilter") > 0 &&
                                  varI.at(CoinFilter) > 0);
    
    
    if (pass_original_filters) {
      h_hf_before_cut->Fill(current_hihf);
    }

    //const bool pass = pass_original_filters && pass_pileup_cut;
    const bool pass = pass_original_filters;

    if (pass) {
      h_hf_after_cut->Fill(current_hihf);
      const auto& parameter = current_hihf;

      // Treat all runs inclusively
      hfdata.push_back({parameter, true});
      hibin.push_back({varI.at("hiBin"), true});

      // Fill one histogram and perform calculations for all runs
      hfData1->Fill(parameter);
      if (parameter > threshold) {
        values.push_back({parameter, false});
      }
      if (parameter > threshold_Min && parameter < thresholdMax) {
        mcYscale_data += 1;
      }
    }
  }//data events loop


  TCanvas* c_pileup_cut = new TCanvas("c_pileup_cut", "2D Data with Linear Pileup Cut", 900, 700);
  h_zdc_vs_hihf->Draw("COLZ");

  TLine* cutLine = new TLine(0, zdc_cut_intercept, hihf_cut_intercept, 0);
  cutLine->SetLineColor(kRed);
  cutLine->SetLineWidth(3);
  cutLine->SetLineStyle(2); // Dashed line
  cutLine->Draw("SAME");

  // =====================  Process MC ======================
   TFile inputMCfile("/eos/cms/store/group/phys_heavyions/nsaha/GO2024/2024PbPbRun3/forest_2024Run3_HYD2024_TuneCELLO_official_21072025/Hydjet_MinBias_TuneCELLO_5p36TeV_pythia8/HiForest_2024Run3_HYD2024_TuneCELLO_official_21072025/250721_164051/0000/HiForestMiniAOD_out_combined.root","READ");

  if (!inputMCfile.IsOpen()) throw std::logic_error("MC file was not found!");

  const auto& tmc = inputMCfile.Get<TTree>("hiEvtAnalyzer/HiTree");
  const auto& tskimanalysismc = inputMCfile.Get<TTree>("skimanalysis/HltTree");
  const auto& thltanalysismc = inputMCfile.Get<TTree>("hltanalysis/HltTree");

  tmc->AddFriend(tskimanalysismc);
  tmc->AddFriend(thltanalysismc);

  std::map<std::string, int> mcVarI;
  std::map<std::string, float> mcVarF;

  std::vector<std::string> mcIntBranches = {"pprimaryVertexFilter", "pclusterCompatibilityFilter", std::string(CoinFilter)};
  for (const std::string& p : mcIntBranches) tmc->SetBranchAddress(p.c_str(), &(mcVarI[p]));

  for (const std::string& p : std::vector<std::string>{std::string(variable), "vz"})
    tmc->SetBranchAddress(p.c_str(), &(mcVarF[p]));

  tmc->SetBranchStatus("*", 0);

  std::vector<std::string> mcActiveBranches = {std::string(variable), "pprimaryVertexFilter", "pclusterCompatibilityFilter", std::string(CoinFilter), "vz"};
  for (const std::string& p : mcActiveBranches) tmc->SetBranchStatus(p.c_str(), 1);

  Nevents = tmc->GetEntries();
  double mcYscale_mc(0);
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    if(iev%5000 == 0) cout<<"Processing mc event: " << iev << " / " << Nevents << endl;
    tmc->GetEntry(iev);
    const auto parameter = mcVarF.at(variable) * mcXscale;

    const bool pass = (mcVarI.at("pprimaryVertexFilter")>0 && mcVarI.at("pclusterCompatibilityFilter")>0 && mcVarI.at(CoinFilter)>0);

    if (pass) {
      hfMc1->Fill(parameter);
      if (parameter>threshold_Min && parameter<thresholdMax)
        mcYscale_mc += 1;
    }
    if (parameter <= threshold)
      values.push_back({parameter, true});
    hfMc2->Fill(parameter);
    mcEff->Fill(pass, parameter);
  } //end of mc loop
  inputMCfile.Close();

  // Scale MC
  const auto mcYscale = mcYscale_data / mcYscale_mc;
  std::cout<<"[INFO] mcYscale = "<< mcYscale << std::endl;
  hfMc1->Scale(mcYscale);
  hfMc2->Scale(mcYscale);
  for (const auto& v : values) {
    nt->Fill(v.first, v.second ? mcYscale : 1.);
    hfCombined->Fill(v.first, v.second ? mcYscale : 1.);
  }
  const auto totEff = hfData1->Integral() / hfCombined->Integral();

  const auto passed = values.size();
  double totalXsec(0);
  for(const auto& v : values)
    totalXsec += v.second ? mcYscale : 1.;

  std::cout << std::endl;
  std::cout << "Selected events = " << passed << std::endl;
  std::cout << "Selected weighed events = " << totalXsec << std::endl;
  std::cout << "Total efficiency = " << totEff << std::endl;

  // Create text file
  ofstream txtfile(Form("output_DataPbPb_usingMC_hiHF_%s.txt", outputTag.c_str()));
  txtfile << "Tag name: " << tag << endl;
  txtfile << "Number of events = " << passed << std::endl;
  txtfile << "Number of weighted events = " << totalXsec << std::endl;
  txtfile << "Total efficiency = " << totEff << std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << "Using MC to correct for peripheral events. Threshold = " << threshold<<". mc X scale factor = "<< mcXscale<< std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << "hiHF based cuts are: " << std::endl;
  txtfile << "(";

  // Store bin boundaries
  const auto size = values.size();
  std::sort(values.begin(), values.end());
  std::vector<double> binboundaries(nbins+1);
  binboundaries[0] = 0.;
  binboundaries[nbins] = values[size-1].first;
  std::cout << "Events per bin = " << totalXsec/nbins << std::endl;
  double integral(0);
  size_t currentbin(1);
  for(const auto& v : values) {
    const auto& val = v.first;
    integral += v.second ? mcYscale : 1.;
    const auto sum = currentbin*(totalXsec/nbins);
        if(integral > sum) {
          binboundaries[currentbin] = val >= 0 ? val : 0;
          currentbin++;
        }
  }
  std::cout << currentbin << std::endl;
  for(size_t i = 0; i < nbins; i++)
    txtfile << binboundaries[i] << ", ";
  txtfile << binboundaries[nbins] << ")" << std::endl;
  txtfile << std::endl;
  txtfile<<"-------------------------------------"<<endl;

  txtfile<<"# Bin BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    int ii = nbins-i;
    bins->table_[i].bin_edge = binboundaries[ii-1];
    txtfile << i << " " << bins->table_[i].bin_edge << " " << endl;
  }
  txtfile << endl;
  txtfile<<"-------------------------------------"<<endl;
  txtfile.close();


  // Check bin boundaries
  int newbin, oldbin;

  TTree t1("CentralityBin","CentralityBin");

  t1.Branch("newBin",&newbin,"newBin/I");
  t1.Branch("oldBin",&oldbin,"oldBin/I");

  for (size_t i=0; i<hfdata.size(); i++) {
    newbin = 199;
    for(size_t b = 0; b < 200; ++b){
      if(hfdata[i].first >= binboundaries[199-b]){
        newbin = b;
        break;
      }
    }
    oldbin = hibin[i].first;
    t1.Fill();
  }

  outFile->cd();
  TDirectory *dir = outFile->mkdir(tag);
  dir->cd();
  bins->Write();
  nt->Write();
  hfData1->Write();
  hfMc1->Write();
  hfMc2->Write();
  mcEff->Write();
  dataEff1->Write();
  hfCombined->Write();

  h_zdc_vs_hihf->Write();
  c_pileup_cut->Write();
  h_hf_before_cut->Write();
  h_hf_after_cut->Write();

  t1.Write();
  outFile->Close();
}

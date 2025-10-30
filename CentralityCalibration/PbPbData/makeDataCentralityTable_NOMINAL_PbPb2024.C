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

void makeDataCentralityTable_NOMINAL_PbPb2024(
				     const TString input_file = "Forest_2024B_v1.txt",
				     const bool HLT = true,
				     const char* HLT_trg = "HLT_HIMinimumBiasHF1ANDZDC1nOR_v4",
				     const char* variable = "hiHF",
				     const int RUN = 388401, 
				     const char*RawPrime = "RawPrime0_usinghiHF", 
				     const char* date = "Oct29",
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

  TH1F*hfData1 = new TH1F("hfData1","hfData1", 100,0, 10000);
  TH1F*hfData2 = new TH1F("hfData2","hfData2", 100,0, 10000);
  TH1F*hfMc1 = new TH1F("hfMc1","hf mc", 100, 0, 10000);
  TH1F*hfMc2 = new TH1F("hfMc2","hf mc", 100, 0, 10000);
  TH1F* hfCombined = new TH1F("hfCombined","hf_combined", 100,0, 10000);
  TEfficiency*dataEff1 = new TEfficiency("dataEff1", "dataEff1", 1000, 0, 1000);
  TEfficiency*dataEff2 = new TEfficiency("dataEff2", "dataEff2", 1000, 0, 1000);
  TEfficiency*mcEff= new TEfficiency("mcEff", "mcEff", 1000, 0, 1000);


  const int runNum = 1;
  CentralityBins*bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);


  UInt_t run;
  t->SetBranchAddress("run", &run);
  t->SetBranchStatus("*", 0);
  
  std::map<std::string, int> varI;
  std::map<std::string, float> varF;
  
  std::vector<std::string> intBranches = {"pprimaryVertexFilter", "pclusterCompatibilityFilter", std::string(CoinFilter), "hiBin"};
  if (HLT) intBranches.insert(intBranches.begin(), std::string(HLT_trg));
  for (const std::string& p : intBranches) t->SetBranchAddress(p.c_str(), &(varI[p]));
  
  for (const std::string& p : std::vector<std::string>{std::string(variable), "vz"})
    t->SetBranchAddress(p.c_str(), &(varF[p]));
  
  std::vector<std::string> activeBranches = {"run", std::string(variable), "pprimaryVertexFilter", "pclusterCompatibilityFilter", std::string(CoinFilter), "vz", "hiBin"};
  if (HLT) activeBranches.insert(activeBranches.begin() + 2, std::string(HLT_trg));
  for (const std::string& p : activeBranches) t->SetBranchStatus(p.c_str(), 1);
  
  
  
  std::vector<std::pair<float, bool>> values, hfdata;
  std::vector<std::pair<int, bool>> hibin;
  TH1D::SetDefaultSumw2();

  std::array<std::array<size_t, 500>, 4> numMinHFTowerV{};

  auto Nevents = t->GetEntries();
  double mcYscale_data(0);

  std::cout<<"Total Number of events = "<<Nevents<<std::endl;
 
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    
    if(iev%100000 == 0) cout<<"Processing data event: " << iev << " / " << Nevents << endl;

    t->GetEntry(iev);

    
    const auto& parameter = varF.at(variable);
    const auto& numMinHFTower = varI.at(CoinFilter);

    

    bool pass;
    if (HLT) {pass= (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && varI.at(CoinFilter)>0);}
    else {pass= (varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && varI.at(CoinFilter)>0);}
 
    
    if (pass) {
      hfdata.push_back({parameter, (run == RUN)});
      hibin.push_back({varI.at("hiBin"), (run == RUN)});
      if (run == RUN) {
        hfData1->Fill(parameter);
        if (parameter > threshold)
          values.push_back({parameter, false});
        if (parameter>threshold_Min && parameter<thresholdMax)
          mcYscale_data += 1;
      }
      else if (run != RUN)
        hfData2->Fill(parameter);
    }

    if (pass){
      ((run == RUN) ? dataEff1 : dataEff2)->Fill(pass, parameter);
    }


    
 }//data events loop

  std::cout<<"mcYscale_data="<<mcYscale_data<<std::endl;

  // Process MC
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
  //txtfile << "Input tree: " << inFileName << endl;
  txtfile << "Tag name: " << tag << endl;
  txtfile << "Number of events = " << passed << std::endl;
  txtfile << "Number of weighted events = " << totalXsec << std::endl;
  txtfile << "Total efficiency = " << totEff << std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << "Using MC to correct for peripheral events. Threshold = " << threshold<<". mc X scale factor = "<< mcXscale<< std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << variable <<"based cuts are: " << std::endl;
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
	  //std::cout << "current bin = " << currentbin << " ; integral = " << integral << " ; sum = " << sum << std::endl;
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
  TTree t1("CentralityBin_wRun","CentralityBin_wRun");
  TTree t2("CentralityBin_woRun","CentralityBin_woRun");

  for (auto& t : {&t1, &t2}){
    t->Branch("newBin",&newbin,"newBin/I");
    t->Branch("oldBin",&oldbin,"oldBin/I");
  }
  for (size_t i=0; i<hfdata.size(); i++) {
    newbin = 199;
    for(size_t b = 0; b < 200; ++b){
      if(hfdata[i].first >= binboundaries[199-b]){
        newbin = b;
        break;
      }
    }
    oldbin = hibin[i].first;
    
    (hfdata[i].second ? t1 : t2).Fill();
  }

  outFile->cd();
  TDirectory *dir = outFile->mkdir(tag);
  dir->cd();
  bins->Write();
  nt->Write();
  hfData1->Write();
  hfData2->Write();
  hfMc1->Write();
  hfMc2->Write();
  mcEff->Write();
  dataEff1->Write();
  dataEff2->Write();
  hfCombined->Write();

  t1.Write();
  t2.Write();
  outFile->Close();

  }


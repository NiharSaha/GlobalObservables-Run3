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

void makeDataCentralityTable_NOMINAL_PbPb2023(
					      const TString input_file = "forest_run387973_PbPb2023.txt",
					      const char * inputMC_file = "/eos/cms/store/group/phys_heavyions/nsaha/GO2024/2024PbPbRun3/forest_2024Run3_HYD2024_PFcand_27112024/Hydjet2024_v2/HiForest_2024Run3_HYD2024_new_PFcand_27112024/241127_054424/0000/HiForestMiniAOD_2024Run3_HYD2024_PFcand_27112024_out_combined.root"
					      const char* HLT_trg = "HLT_HIMinimumBiasHF1ANDZDC1nOR_v4",
					      const int RUN = 387973, 
					      const char*RawPrime = "RawPrime0", 
					      const char* CoinFilter = "pphfCoincFilter2Th4",
					      const double threshold = 100.0,
					      const char* label = "Nominal", 
					      const size_t nbins = 200
					      )
{
  // Constant parameters
  const auto mcXscale = 1.008;
  const auto threshold_Min = 1000.0;
  const auto thresholdMax = 4000.0;

  //Tag
  //const char* testMessage = Form("run%d_HIPhysics%s",RUN, RawPrime);
  const char* testMessage = Form("HIPhysics%s", RawPrime);
  //const char* tag = Form("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_%s", label) ;
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
  //TNtuple nt("nt","","value");
  // TH1F*hfData1 = new TH1F("hfData1",Form("hf data run == %d",RUN), 100,0, 10000);
  //TH1F*hfData2 = new TH1F("hfData2",Form("hf data run != %d",RUN), 100,0, 10000);

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
  //CentralityBins * bins = new CentralityBins();
  bins->table_.reserve(nbins);



  UInt_t run;
  t->SetBranchAddress("run", &run);
  std::map<std::string, int> varI, mcVarI;
  //const char* numMinHFTowerLbl = Form("pphfCoincFilter2Th%d", hfCoinThr);
  for (const auto& p : {HLT_trg, "pprimaryVertexFilter", "pclusterCompatibilityFilter", CoinFilter, "hiBin"})
    t->SetBranchAddress(p, &(varI[p]));
  std::map<std::string, float> varF, mcVarF;
  for (const auto& p : {"hiHF", "vz"})
    t->SetBranchAddress(p, &(varF[p]));
  t->SetBranchStatus("*", 0);
  for (const auto& p : {"run", "hiHF", HLT_trg, "pprimaryVertexFilter", "pclusterCompatibilityFilter", CoinFilter, "vz", "hiBin"})
    t->SetBranchStatus(p, 1);

  
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

    
    const auto& parameter = varF.at("hiHF");
    const auto& numMinHFTower = varI.at(CoinFilter);

    
    //const bool pass = (varI.at("HLT_HIMinimumBiasHF1AND_v1")>0 ); //with only HLT trigger!!
    const bool pass = (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && varI.at(CoinFilter)>0);
  
    //std::cout<<"***************"<<std::endl;    
    //To check with only HLT_HIZeroBias_v4
    //const bool pass = (varI.at("HLT_HIZeroBias_v4")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && numMinHFTower>=hfCoinN);

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
    if (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {

      //To check with only HLT_HIZeroBias_v4
      //if (varI.at("HLT_HIZeroBias_v4")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {
      //const bool passTight = (pass && varI.at(numMinHFTowerLbl)>=(hfCoinN+1));
      ((run == RUN) ? dataEff1 : dataEff2)->Fill(pass, parameter);
    }

    if (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {

    //if (varI.at("HLT_HIZeroBias_v4")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {
      size_t idx = (parameter > 1000.)*2 + (run == RUN);
      if ((parameter > 1000.) || (parameter < 140.))//(parameter > 25. && parameter < 60.))//(parameter > 60. && parameter < 140.))
        numMinHFTowerV[idx][numMinHFTower] += 1;
    }
    


    //****W/O RUN number *******
    /*    if (pass) {
      hfdata.push_back({parameter, ""});
      hibin.push_back({varI.at("hiBin"), ""});
      //if (run == RUN) {
      hfData1->Fill(parameter);
      if (parameter > threshold)
	values.push_back({parameter, false});
      if (parameter>threshold_Min && parameter<thresholdMax)
	mcYscale_data += 1;
      //}
      //else if (run != RUN)
      // hfData2->Fill(parameter);
    }
    

  //=*********** Noise study *************

  if (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {

      //To check with only HLT_HIZeroBias_v4
      //if (varI.at("HLT_HIZeroBias_v4")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {
      //const bool passTight = (pass && varI.at(numMinHFTowerLbl)>=(hfCoinN+1));
      ((run == RUN) ? dataEff1 : dataEff2)->Fill(pass, parameter);
    }

    if (varI.at(HLT_trg)>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {

    //if (varI.at("HLT_HIZeroBias_v4")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {
    //  size_t idx = (parameter > 1000.)*2 + (run == RUN);
      size_t idx = (parameter > 1000.)*2 ;
      if ((parameter > 1000.) || (parameter < 140.))//(parameter > 25. && parameter < 60.))//(parameter > 60. && parameter < 140.))
        numMinHFTowerV[idx][numMinHFTower] += 1;
    }
    //=**************** End noise study ***********

    */
    
 }//data events loop
  //inFile.Close();

  // Compute noise study
  //TH1F*hfNoise1 = new TH1F("hfNoise1",Form("hf noise data run == %d",RUN), 60, 0, 60);
  //TH1F*hfNoise2 = new TH1F("hfNoise2",Form("hf noise data run != %d",RUN), 60, 0, 60);

  TH1F*hfNoise1 = new TH1F("hfNoise1","hfNoise1", 60, 0, 60);
  TH1F*hfNoise2 = new TH1F("hfNoise2","hfNoise2", 60, 0, 60);
  std::array<std::array<double, 200>, 4> nEvt{};
  for (size_t i=0; i<101; i++)
    for (size_t j=0; j<4; j++)
      nEvt[j][i] = std::accumulate(numMinHFTowerV[j].begin()+i, numMinHFTowerV[j].end(), 0);
  for (size_t i=0; i<60; i++) {
    hfNoise1->SetBinContent(i+1, (nEvt[0][i] / nEvt[1][i])*(nEvt[3][i] / nEvt[2][i]));
    hfNoise2->SetBinContent(i+1, nEvt[1][i] / nEvt[1][20]);
  }

  // ----------  Process MC -------------

  //official MC forest file

  TFile inputMCfile(inputMC_file,"READ");	
  if (!inputMCfile.IsOpen()) throw std::logic_error("MC file was not found!");
  const auto& tmc = inputMCfile.Get<TTree>("hiEvtAnalyzer/HiTree");
  const auto& tskimanalysismc = inputMCfile.Get<TTree>("skimanalysis/HltTree");
  const auto& thltanalysismc = inputMCfile.Get<TTree>("hltanalysis/HltTree");
  tmc->AddFriend(tskimanalysismc);
  tmc->AddFriend(thltanalysismc);

  for (const auto& p : { "pprimaryVertexFilter", "pclusterCompatibilityFilter", CoinFilter})
    //for (const auto& p : { "pprimaryVertexFilter", "pclusterCompatibilityFilter", "pphfCoincFilter2Th4"})

  tmc->SetBranchAddress(p, &(mcVarI[p]));
  for (const auto& p : {"hiHF", "vz"})
    tmc->SetBranchAddress(p, &(mcVarF[p]));
  tmc->SetBranchStatus("*", 0);
  for (const auto& p : {"hiHF", "pprimaryVertexFilter", "pclusterCompatibilityFilter", CoinFilter, "vz"})
    //for (const auto& p : {"hiHF", "pprimaryVertexFilter", "pclusterCompatibilityFilter", "pphfCoincFilter2Th4", "vz"})
    tmc->SetBranchStatus(p, 1);

  

  Nevents = tmc->GetEntries();
  double mcYscale_mc(0);
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    if(iev%5000 == 0) cout<<"Processing mc event: " << iev << " / " << Nevents << endl;
    tmc->GetEntry(iev);
    const auto parameter = mcVarF.at("hiHF") * mcXscale;

    const bool pass = (mcVarI.at("pprimaryVertexFilter")>0 && mcVarI.at("pclusterCompatibilityFilter")>0 && mcVarI.at(CoinFilter)>0);
    //const bool pass = (mcVarI.at("pprimaryVertexFilter")>0 && mcVarI.at("pclusterCompatibilityFilter")>0 && mcVarI.at("pphfCoincFilter2Th4")>0);


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


  
  //-------------- Scale MC --------------
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

  // ----------  Create text file where hiHF bin edges (table) are saved -------------

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
  txtfile << "hiHF based cuts are: " << std::endl;
  txtfile << "(";

  //------------------------ Store bin boundaries ----------------------
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


  //-------- Store histograms ------------
  outFile->cd();
  //const auto& dir = outFile.mkdir(tag.c_str());
  TDirectory *dir = outFile->mkdir(tag);
  dir->cd();
  bins->Write();
  nt->Write();
  hfData1->Write();
  hfData2->Write();
  //c->Write();
  hfMc1->Write();
  hfMc2->Write();
  mcEff->Write();
  dataEff1->Write();
  dataEff2->Write();
  hfNoise1->Write();
  hfNoise2->Write();
  hfCombined->Write();
  outFile->Close();


  
  //--------------- Check bin boundaries (w/ a& w/o calibration) --------------
  int newbin, oldbin;
  TFile outf(Form("compare_centralitybins_%s.root", outputTag.c_str()),"recreate");
  //TTree t1(Form("anaCentrality_%d", RUN),"analysis level centrality");
  //TTree t2(Form("anaCentrality_not%d", RUN),"analysis level centrality");

  TTree t1("anaCentrality","analysis level centrality");
  TTree t2("anaCentrality_","analysis level centrality");
  for (auto& t : {&t1, &t2}) {
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
  t1.Write();
  t2.Write();
  outf.Close();
  }


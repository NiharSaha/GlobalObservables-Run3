#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;


void makeMCCentralityTable_PbPb2024(
			   const char * tag = "CentralityTable_HFtowers200_HydjetCello_v1401x0_official_MC2024",
			   const char* date = "Oct29",
			   const char* label = "HFtowers",
			   const char * variable = "hiHF"

			   )
{
  
  TH1::SetDefaultSumw2();

  int nbins = 200;
  //TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/2023PbPbRun3/forest_Run3_HYD_RECODEBUG_09102023/GENSIM_HydjetnoPU_CMSSW_13_2_0_pre3_28Jun2023/HiForest_Run3_HYD_RECODEBUG_09102023/231009_185152/0000/Forest_Run3_HYD_RECODEBUG_merged.root";

  TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2024/2024PbPbRun3/forest_2024Run3_HYD2024_TuneCELLO_official_21072025/Hydjet_MinBias_TuneCELLO_5p36TeV_pythia8/HiForest_2024Run3_HYD2024_TuneCELLO_official_21072025/250721_164051/0000/HiForestMiniAOD_out_combined.root";
  TString outFileName = Form("CentralityTable_%s200_HydjetCello_official_MC2024_%s.root",label, date);
  ofstream txtfile(Form("output_HydjetCello_using%s_official_MC2024_%s.txt", variable, date));


  TH1F*hfMc = new TH1F("hfMc","hf mc", 100, 0, 10000);
  
  TFile *inFile = TFile::Open(inFileName);

  TTree *t = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *t_skim = (TTree*)inFile->Get("skimanalysis/HltTree");
  t->AddFriend(t_skim);


  
  const int runNum = 1;
  CentralityBins * bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);
  

  txtfile << "Input tree: " << inFileName << endl;
  
  double binboundaries[nbins+1];
  vector<float> values;

  float var;
  int hiBin, pprimaryVertexFilter, pclusterCompatibilityFilter, pphfCoincFilterPF3Th5;
  
  t->SetBranchAddress(variable, &var);
  t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress("pprimaryVertexFilter",	&pprimaryVertexFilter);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  t->SetBranchAddress("pphfCoincFilterPF3Th5", &pphfCoincFilterPF3Th5);

  t->SetBranchStatus("*", 0);
  for (const auto& p : {variable, "hiBin", "pprimaryVertexFilter", "pclusterCompatibilityFilter", "pphfCoincFilterPF3Th5"})
    t->SetBranchStatus(p, 1);

  
  /*bool binB = label.compare("b") == 0;
  bool binNpart = label.compare("Npart") == 0;
  bool binNcoll = label.compare("Ncoll") == 0;
  bool binNhard = label.compare("Nhard") == 0;*/
  //  bool binHF = label.compare("hiHF") == 1;
  /*bool binHFplus = label.compare("HFtowersPlus") == 0;
  bool binHFminus = label.compare("HFtowersMinus") == 0;
  bool binHFplusTrunc = label.compare("HFtowersPlusTrunc") == 0;
  bool binHFminusTrunc = label.compare("HFtowersMinusTrunc") == 0;
  bool binNpix = label.compare("PixelHits") == 0;
  bool binNpixTrks = label.compare("PixelTracks") == 0;
  bool binNtrks = label.compare("Tracks") == 0;*/


    vector<int> hibins;
  const bool pass = (pprimaryVertexFilter>0 && pclusterCompatibilityFilter>0 && pphfCoincFilterPF3Th5>0);
  
  unsigned int Nevents = t->GetEntries();
  txtfile << "Number of events = " << Nevents << endl << endl;
  
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if(iev%100000 == 0) cout<<"Processing event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);



    float parameter = -1;
    if (pprimaryVertexFilter>0 && pclusterCompatibilityFilter>0 && pphfCoincFilterPF3Th5>0){
      parameter = var;
      values.push_back(parameter);
      hibins.push_back(hiBin);
      hfMc->Fill(parameter);
      }
  }
  
  sort(values.begin(),values.end());
  
  txtfile << "-------------------------------------" << endl;
  txtfile << label << " based cuts are: " << endl;
  txtfile << "(";
  
  int size = values.size();
  binboundaries[nbins] = values[size-1];
  
  for(int i = 0; i < nbins; i++) {
    int entry = (int)(i*(size/((float)(nbins))));
    if(entry < 0 || i == 0) binboundaries[i] = 0;
    else binboundaries[i] = values[entry];
  }
  
  for(int i = 0; i < nbins; i++) {
    if(binboundaries[i] < 0) binboundaries[i] = 0;
    txtfile << binboundaries[i] << ", ";
  }
  txtfile << binboundaries[nbins] << ")" << endl << "-------------------------------------" << endl;
  
  txtfile<<"-------------------------------------"<<endl;
  txtfile<<"# Bin  BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    
    int ii = nbins-i;
    bins->table_[i].bin_edge = binboundaries[ii-1];
    
    txtfile << i << " " << binboundaries[ii-1] << " " <<endl;
  }
  txtfile<<"-------------------------------------"<<endl;
  txtfile.close();



  TFile * outFile = new TFile(outFileName,"recreate");
  TDirectory* dir = outFile->mkdir(tag);
  dir->cd();  
  
  // Check bin boundaries
  int newbin, oldbin;

  TTree* t1 = new TTree("anaCentrality","anaCentrality");
  t1->Branch("newBin", &newbin, "newBin/I");
  t1->Branch("oldBin", &oldbin, "oldBin/I");
  

  for (size_t i=0; i<values.size(); i++) {
    newbin = 199;
    for(size_t b = 0; b < 200; ++b){
      if(values[i] >= binboundaries[199-b]){
        newbin = b;
        break;
      }
    }
    oldbin = hibins[i];
    t1->Fill();

  }


  bins->Write();
  hfMc->Write();
  t1->Write();
  bins->Delete();
  outFile->Close();

  
}

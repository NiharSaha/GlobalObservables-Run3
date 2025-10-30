#include "TSystem.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TClonesArray.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
//#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include "TLegend.h"
#include <fstream>
#include "TLatex.h"

using namespace std;


void plot_hibin( TFile*file11,TFile*file12, TFile*file13) {
  gStyle->SetOptStat(0);


  TDirectory * dir11 = (TDirectory*)file11->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal");
  TDirectory * dir12 = (TDirectory*)file12->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal");
  TDirectory * dir13 = (TDirectory*)file13->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal");
  TTree* tree1 = (TTree*)dir11->Get("CentralityBin_wRun");
  TTree* tree2 = (TTree*)dir12->Get("CentralityBin_wRun");
  TTree* tree3 = (TTree*)dir13->Get("CentralityBin_wRun");

  TTreeReader reader1(tree1);
  TTreeReaderValue<Int_t> branchValue1(reader1, "newBin");
  TTreeReader reader2(tree2);
  TTreeReaderValue<Int_t> branchValue2(reader2, "newBin");
  TTreeReader reader3(tree3);
  TTreeReaderValue<Int_t> branchValue3(reader3, "newBin");
  
  TH1D* hibin_nominal = new TH1D("hiBin_nominal", "hiBin_nominal", 100, 0, 200);
  TH1D* hibin_down = new TH1D("hiBin_down", "hiBin_down", 100, 0, 200);
  TH1D* hibin_up = new TH1D("hiBin_up", "hiBin_up", 100, 0, 200); 
  
  while (reader1.Next()) {hibin_nominal->Fill(*branchValue1);}
  while (reader2.Next()) {hibin_down->Fill(*branchValue2);}
  while (reader3.Next()) {hibin_up->Fill(*branchValue3);}
   
  
  TCanvas *canvas = new TCanvas("canvas", "canvas",442,165,857,688);
   gStyle->SetOptStat(0);
   canvas->Range(-41.90494,-1885.814,224.0152,11467.79);
   canvas->SetFillColor(0);
   canvas->SetBorderMode(0);
   canvas->SetBorderSize(2);
   canvas->SetLeftMargin(0.1575847);
   canvas->SetRightMargin(0.06774668);
   canvas->SetBottomMargin(0.1412214);
   canvas->SetFrameBorderMode(0);
   canvas->SetFrameBorderMode(0);


   hibin_down->SetTitle("");
   hibin_down->SetLineColor(30);
   hibin_down->SetMarkerColor(30);
   hibin_down->SetMarkerStyle(29);
   hibin_down->SetMarkerSize(1.4);
   hibin_down->GetXaxis()->SetTitle("hiBin");

   hibin_down->GetYaxis()->SetTitle("Events");
   hibin_down->GetXaxis()->SetRange(1,100);
   hibin_down->GetXaxis()->CenterTitle(true);
   hibin_down->GetXaxis()->SetLabelFont(42);
   hibin_down->GetXaxis()->SetTitleSize(0.05);
   hibin_down->GetXaxis()->SetTitleOffset(1);
   hibin_down->GetXaxis()->SetTitleFont(132);
   hibin_down->GetYaxis()->CenterTitle(true);
   hibin_down->GetYaxis()->SetLabelFont(42);
   hibin_down->GetYaxis()->SetTitleSize(0.05);
   hibin_down->GetYaxis()->SetTitleFont(132);
   
   hibin_down->Draw("E1");

   hibin_up->SetTitle("");
   hibin_up->SetLineColor(9);
   hibin_up->SetMarkerColor(9);
   hibin_up->SetMarkerStyle(22);
   hibin_up->SetMarkerSize(1.1);
   //hibin_up->GetXaxis()->SetTitle("newBin");
   hibin_up->GetXaxis()->SetLabelFont(42);
   hibin_up->GetXaxis()->SetTitleOffset(1);
   hibin_up->GetXaxis()->SetTitleFont(42);
   hibin_up->GetYaxis()->SetLabelFont(42);
   hibin_up->GetYaxis()->SetTitleFont(42);
   
   hibin_up->Draw("ESAME");

   hibin_nominal->SetTitle("");
   hibin_nominal->SetLineColor(2);
   hibin_nominal->SetMarkerColor(2);
   hibin_nominal->SetMarkerStyle(20);
   hibin_nominal->SetMarkerSize(1.2);
   //hibin_nominal->GetXaxis()->SetTitle("hiBin");

   hibin_nominal->Draw("ESAME");

   
   TLegend *leg = new TLegend(0.1944035,0.2461832,0.5891016,0.3721374,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

   TLegendEntry *entry=leg->AddEntry("hibin_nominal","Threshold =100 (Nominal)","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hibin_up","Threshold = 300 (Syst UP)","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(9);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(9);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hibin_down","Threshold = 50 (Syst DOWN)","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(30);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(30);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.4);
   entry->SetTextFont(42);
   leg->Draw();
   
   TLatex *   tex = new TLatex(139.0304,10422.94,"PbPb #sqrt{S_{NN}} = 5.36 TeV (2024)");
   tex->SetTextFont(132);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(-0.7832696,10244.56,"#bf{CMS} #it{Preliminary}");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   /*  tex = new TLatex(36.78254,3362.414,Form("Run = %d", RUN));
   tex->SetTextColor(2);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   */
  
  //canvas->SaveAs("histogram.png");
  

  // file1->Close();
}




//void plot_hiHF( TFile* file2, const int RUN) {
void plot_hiHF( TFile* file2) {
  gStyle->SetOptStat(0);

  TDirectory * dir = (TDirectory*)file2->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal");
  TH1F* hist_data = (TH1F*)dir->Get("hfData1");
  TH1F* hist_MC = (TH1F*)dir->Get("hfMc1");
  TH1F* hist_combined = (TH1F*)dir->Get("hfCombined");


  TCanvas* canvas = new TCanvas("canvas", "Canvas Title", 185,180,825,720);
  canvas->Range(0,0,1,1);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetFrameBorderMode(0);

  TPad* pad1 = new TPad("pad1", "Upper Pad", 0.0, 0.3, 1.0, 1.0);
  pad1->Draw();
  pad1->cd();
  pad1->Range(-875.0001,-0.5672784,7875,5.924375);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetLogy();
  pad1->SetBottomMargin(0.005235602);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderMode(0);
  
  hist_data->SetTitle("");
  hist_data->SetLineColor(1);
  hist_data->GetXaxis()->SetRangeUser(0, 7000);
  hist_data->SetMarkerStyle(20);
  hist_data->GetXaxis()->SetLabelFont(42);
  hist_data->GetXaxis()->SetTitleOffset(1);
  hist_data->GetXaxis()->SetTitleFont(42);
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->CenterTitle(true);
  hist_data->GetYaxis()->SetLabelFont(42);
  hist_data->GetYaxis()->SetLabelSize(0.05);
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(0.86);
  hist_data->GetYaxis()->SetTitleFont(132);
  hist_data->Draw("E1");
  

  hist_MC->SetTitle("");
  hist_MC->SetLineColor(2);
  hist_MC->SetLineWidth(2);
  hist_MC->SetMarkerColor(2);
  hist_MC->SetMarkerSize(1.5);
  hist_MC->GetXaxis()->SetRangeUser(0, 7000);
  hist_MC->GetXaxis()->SetLabelFont(42);
  hist_MC->GetXaxis()->SetTitleOffset(1);
  hist_MC->GetXaxis()->SetTitleFont(42);
  hist_MC->GetYaxis()->SetLabelFont(42);
  hist_MC->GetYaxis()->SetTitleFont(42);
  hist_MC->Draw("HISTSAME");
  
  hist_combined->SetTitle("");
  hist_combined->SetLineColor(4);
  hist_combined->SetMarkerColor(4);
  hist_combined->SetMarkerStyle(25);
  hist_combined->SetMarkerSize(1.4);
  hist_combined->GetXaxis()->SetRangeUser(0, 7000);
  hist_combined->GetXaxis()->SetLabelFont(42);
  hist_combined->GetXaxis()->SetTitleOffset(1);
  hist_combined->GetXaxis()->SetTitleFont(42);
  hist_combined->GetYaxis()->SetLabelFont(42);
  hist_combined->GetYaxis()->SetTitleFont(42);
  hist_combined->Draw("E1SAME");
   
 
   
   TLegend *legend = new TLegend(0.5893074,0.6690647,0.8918591,0.8787256,NULL,"brNDC");
   legend->SetBorderSize(0);
   legend->SetLineColor(1);
   legend->SetLineStyle(1);
   legend->SetLineWidth(1);
   legend->SetFillColor(0);
   legend->SetFillStyle(1001);
   legend->AddEntry(hist_data, "2024 RawPrime", "lp");
   legend->AddEntry(hist_MC, "Hydjet MC x 1.012", "lp");
   legend->AddEntry(hist_combined, "Combined", "lp");
   legend->Draw();


   TLatex *   tex = new TLatex(4834.295,253498.6,"PbPb #sqrt{S_{NN}} = 5.36 TeV");
   tex->SetTextFont(132);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(-45.71694,261408.1,"CMS Preliminary");
   tex->SetTextFont(52);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(602.825,10.8071,"HLT_HIMinimumBiasHF1ANDZDC1nOR_v4");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
   tex = new TLatex(602.825,4.396605,"PrimaryVertexFilter");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(581.5613,2.079979,"ClusterCompatibilityFilter");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(613.4568,0.9840117,"pfConicFilterPF3Th5");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   /*tex = new TLatex(2527.19,49145,Form("Run = %d", RUN));
   tex->SetTextColor(2);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();*/
      tex = new TLatex(2505.92,18213.2,"Threshold = 100");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04110997);
   tex->SetLineWidth(2);
   tex->Draw();


   TLine *l = new TLine(100.0,0.0,100.0, 2*hist_MC->GetMaximum());
   l->SetLineColor(2);
   l->SetLineStyle(2);
   l->SetLineWidth(2);
   l->Draw();

   pad1->Modified();
   //canvas->cd();
   
   
   canvas->cd();
   TPad* pad2 = new TPad("pad2", "Lower Pad", 0.0, 0.0, 1.0, 0.3);
   pad2->SetTopMargin(0.02);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-1222.445,-1.035714,11262.52,2.033929);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetTopMargin(0.01105294);
   pad2->SetBottomMargin(0.3374055);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);
   
   TH1D* ratio = (TH1D*)hist_data->Clone("ratio");
   ratio->Divide(hist_MC);
   ratio->SetLineColor(kBlack);
   ratio->GetYaxis()->SetRangeUser(0.0, 2.0); // Set the Y-axis range as needed
   ratio->GetXaxis()->SetRangeUser(0.0, 7000.0);
   ratio->SetTitle("");
   ratio->SetMinimum(0);
   ratio->SetMaximum(2);
   ratio->SetLineWidth(0);
   ratio->SetMarkerStyle(20);
   ratio->GetXaxis()->SetTitle("HF E_{T} (GeV)");
   ratio->GetXaxis()->CenterTitle(true);
   ratio->GetXaxis()->SetLabelFont(42);
   ratio->GetXaxis()->SetLabelSize(0.1);
   ratio->GetXaxis()->SetTitleSize(0.13);
   ratio->GetXaxis()->SetTitleOffset(1.01);
   ratio->GetXaxis()->SetTitleFont(132);
   ratio->GetYaxis()->SetTitle("Data/MC");
   ratio->GetYaxis()->CenterTitle(true);
   ratio->GetYaxis()->SetNdivisions(509);
   ratio->GetYaxis()->SetLabelFont(42);
   ratio->GetYaxis()->SetLabelSize(0.08);
   ratio->GetYaxis()->SetTitleSize(0.09);
   ratio->GetYaxis()->SetTitleOffset(0.44);
   ratio->GetYaxis()->SetTitleFont(132);
   ratio->Draw("E1");
    
   TLine *line = new TLine(0.0,1.0,7000.0,1.0);
   line->SetLineColor(2);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();
   
   pad2->Modified();
   canvas->cd();
   canvas->Modified();
   //canvas->cd();
   //canvas->SetSelected(canvas);
   
    //canvas->SaveAs("comparison_with_ratio.png");


    //file2->Close();
}

void plot_all_PbPb2024(){

    const int RUN = 388401;
    
    //Nominal
    TString inFile1 = "CentralityTable_HFtowers200_DataPbPb_usingMC_2024Run_HYDMC_xSF1.01_pphfCoincFilterPF3Th5_Threshold100_Nominal_Normalisation1000_4000_HIPhysicsRawPrime0_usinghiHF_Oct29.root";
    //Sys Down
    TString inFile2 = "CentralityTable_HFtowers200_DataPbPb_usingMC_2024Run_HYDMC_xSF1.01_pphfCoincFilterPF3Th5_Threshold50_SysDown_Normalisation1000_4000_HIPhysicsRawPrime0_usinghiHF_Oct29.root";
    //Sys Up
    TString inFile3 = "CentralityTable_HFtowers200_DataPbPb_usingMC_2024Run_HYDMC_xSF1.01_pphfCoincFilterPF3Th5_Threshold300_SysUp_Normalisation1000_4000_HIPhysicsRawPrime0_usinghiHF_Oct29.root";




    TFile *file1 = TFile::Open(inFile1);
    TFile *file2 = TFile::Open(inFile1);
    TFile *file3 = TFile::Open(inFile1);



    plot_hibin(file1, file2, file3);    

    //plot_hiHF(file1);
   



}



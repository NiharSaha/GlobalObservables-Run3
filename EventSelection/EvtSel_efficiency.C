#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

// --- Main function for the analysis ---
void EvtSel_efficiency() 
{
    using std::cout;
    using std::endl;
    using std::ifstream;
    using std::string;

    // =================================================================
    // ### 1. File Loading and TTree Setup
    // =================================================================

    const TString input_file = "forest_2024B.txt";
    TFile * fout = new TFile("EvtSelEff_Oct29.root","RECREATE");

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
    TChain * thltanalysis = new TChain("hltanalysis/HltTree");

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

    // =================================================================
    // ### 2. Branch and Histogram Setup (Expanded)
    // =================================================================
    
    
    std::vector<std::string> hlt_trigger_names;
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
    
    std::map<std::string, Int_t> hltTriggerVars;
    for (const auto& name : hlt_trigger_names) {
      hltTriggerVars[name] = 0; // Initialize
    }
    
    //--For different coincFilters list here --------
    std::vector<std::string> filterNames = {
        "pphfCoincFilterPF2Th3", "pphfCoincFilterPF2Th4", "pphfCoincFilterPF2Th5",
        "pphfCoincFilterPF3Th3", "pphfCoincFilterPF3Th4", "pphfCoincFilterPF3Th5",
        "pphfCoincFilterPF4Th3", "pphfCoincFilterPF4Th4", "pphfCoincFilterPF4Th5"
    };

    std::map<std::string, Int_t> filterVars;
    for (const auto& name : filterNames) {
        filterVars[name] = 0; // Initialize
    }



    // --- Basic variables ---
    Int_t hiBin, pclusterCompatibilityFilter, pprimaryVertexFilter, HLT;
    Float_t hiHF_pf;
    

    // --- Set Branch Addresses ---
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("hiHF_pf", &hiHF_pf);
    t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
    t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);


    for (const auto& name : hlt_trigger_names) {
      t->SetBranchAddress(name.c_str(), &hltTriggerVars[name]);
    }
    for (const auto& name : filterNames) {
      t->SetBranchAddress(name.c_str(), &filterVars[name]);
    }
    
    // --- Set Branch Status ---
    t->SetBranchStatus("*", 0);
    std::vector<const char*> branchesToEnable = {
      "hiBin", "hiHF_pf", "pprimaryVertexFilter", "pclusterCompatibilityFilter"
    };

    for (const auto& name : hlt_trigger_names) {
      branchesToEnable.push_back(name.c_str());
    }
    
    for (const auto& name : filterNames) {
      branchesToEnable.push_back(name.c_str());
    }
    
    for (const auto& branch : branchesToEnable) {
        t->SetBranchStatus(branch, 1);
    }

    // --- Define Binning (CHANGE THESE to match your needs) ---
    int nBins_hiBin = 200; double min_hiBin = 0; double max_hiBin = 200;
    int nBins_hiHF = 200; double min_hiHF = 0; double max_hiHF = 200;

    // --- Denominator Histograms (for 'pass' selection) ---
    TH1D* h_hiBin_total = new TH1D("h_hiBin_total", "hiBin (pass);hiBin;Events", nBins_hiBin, min_hiBin, max_hiBin);
    TH1D* h_hiHF_total = new TH1D("h_hiHF_total", "hiHF_pf (pass);hiHF_pf;Events", nBins_hiHF, min_hiHF, max_hiHF);
    h_hiBin_total->Sumw2();
    h_hiHF_total->Sumw2();

    // --- Numerator Histograms (Maps for all filters) ---
    std::map<std::string, TH1D*> h_hiBin_pass;
    std::map<std::string, TH1D*> h_hiHF_pass;

    for (const auto& name : filterNames) {
        // hiBin numerator
        TString hname_bin = TString::Format("h_hiBin_pass_%s", name.c_str());
        TString htitle_bin = TString::Format("hiBin (pass && %s);hiBin;Events", name.c_str());
        h_hiBin_pass[name] = new TH1D(hname_bin, htitle_bin, nBins_hiBin, min_hiBin, max_hiBin);
        h_hiBin_pass[name]->Sumw2();

        // hiHF numerator
        TString hname_hf = TString::Format("h_hiHF_pass_%s", name.c_str());
        TString htitle_hf = TString::Format("hiHF_pf (pass && %s);hiHF_pf;Events", name.c_str());
        h_hiHF_pass[name] = new TH1D(hname_hf, htitle_hf, nBins_hiHF, min_hiHF, max_hiHF);
        h_hiHF_pass[name]->Sumw2();
    }
    
    cout << "Histograms defined for " << filterNames.size() << " filters." << endl;

    // =================================================================
    // ### 3. Run Event Loop
    // =================================================================

    Int_t nevent = t->GetEntries();
    cout << "Total events to process: " << nevent << endl;

    for (int ievt = 0; ievt < nevent; ievt++) {
        if (ievt % 100000 == 0) cout << "Processing data event: " << ievt << " / " << nevent << endl;

        t->GetEntry(ievt);

	bool passHLT = false;
        for (const auto& name : hlt_trigger_names) {
	  if (hltTriggerVars[name]) {
	    passHLT = true;
	    break; // Found one that passed, no need to check others
	  }
        }

	bool pass = passHLT && pprimaryVertexFilter && pclusterCompatibilityFilter;

        if (pass) {
            // 1. Fill Denominator histograms
            h_hiBin_total->Fill(hiBin);
            h_hiHF_total->Fill(hiHF_pf);

            // 2. Loop over all filters and fill numerators
            for (const auto& name : filterNames) {
                if (filterVars[name]) { // Check if the filter variable is true
                    h_hiBin_pass[name]->Fill(hiBin);
                    h_hiHF_pass[name]->Fill(hiHF_pf);
                }
            }
        }
    } //--- End of event loop ---

    cout << "Event loop finished. Calculating, plotting, and saving..." << endl;

    // =================================================================
    // ### 4. Calculate, Plot, and Save (Expanded)
    // =================================================================

    // --- Define maps for efficiency histograms ---
    std::map<std::string, TH1D*> h_eff_hiBin;
    std::map<std::string, TH1D*> h_eff_hiHF;

    // --- Define styles for 12 plots ---
    std::vector<int> colors = {
        kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange+7,
        kRed-7, kBlue-7, kGreen-3, kMagenta-7, kCyan-7, kBlack
    };
    std::vector<int> markers = {
        kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown,
        kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown,
        kFullStar, kOpenStar, kFullCross, kOpenCross
    };
    
    // --- Create Canvases and Legends ---
    gStyle->SetOptStat(0); // Turn off stat box
    TCanvas* c_hiBin = new TCanvas("c_hiBin", "hiBin Efficiency Comparison", 1000, 800);
    c_hiBin->cd();
    gPad->SetLeftMargin(0.12);
    // Adjust legend position and size for many entries
    TLegend* leg_hiBin = new TLegend(0.55, 0.15, 0.9, 0.55);
    leg_hiBin->SetNColumns(2); // Use 2 columns
    leg_hiBin->SetBorderSize(0);

    TCanvas* c_hiHF = new TCanvas("c_hiHF", "hiHF Efficiency Comparison", 1000, 800);
    c_hiHF->cd();
    gPad->SetLeftMargin(0.12);
    TLegend* leg_hiHF = new TLegend(0.55, 0.15, 0.9, 0.55);
    leg_hiHF->SetNColumns(2);
    leg_hiHF->SetBorderSize(0);

    // --- Loop over filters to calculate, style, and draw ---
    for (size_t i = 0; i < filterNames.size(); ++i) {
        const auto& name = filterNames[i];

        // --- Calculate Efficiencies ---
        TString eff_hname_bin = TString::Format("h_eff_hiBin_%s", name.c_str());
        h_eff_hiBin[name] = (TH1D*)h_hiBin_pass[name]->Clone(eff_hname_bin);
        h_eff_hiBin[name]->SetTitle("Event Selection Efficiency vs. hiBin");
        h_eff_hiBin[name]->Divide(h_hiBin_total);

        TString eff_hname_hf = TString::Format("h_eff_hiHF_%s", name.c_str());
        h_eff_hiHF[name] = (TH1D*)h_hiHF_pass[name]->Clone(eff_hname_hf);
        h_eff_hiHF[name]->SetTitle("Event Selection Efficiency vs. hiHF_pf");
        h_eff_hiHF[name]->Divide(h_hiHF_total);

        // --- Style ---
        h_eff_hiBin[name]->SetLineColor(colors[i]);
        h_eff_hiBin[name]->SetMarkerColor(colors[i]);
        h_eff_hiBin[name]->SetMarkerStyle(markers[i]);
        
        h_eff_hiHF[name]->SetLineColor(colors[i]);
        h_eff_hiHF[name]->SetMarkerColor(colors[i]);
        h_eff_hiHF[name]->SetMarkerStyle(markers[i]);

        // --- Draw hiBin plot ---
        c_hiBin->cd();
        TString drawOpt = "E1 SAME";
        if (i == 0) { // Draw first plot to set axes
            h_eff_hiBin[name]->GetYaxis()->SetRangeUser(0.0, 1.1);
            h_eff_hiBin[name]->GetYaxis()->SetTitle("Efficiency");
            drawOpt = "E1";
        }
        h_eff_hiBin[name]->Draw(drawOpt);
        leg_hiBin->AddEntry(h_eff_hiBin[name], name.c_str(), "lep");

        // --- Draw hiHF plot ---
        c_hiHF->cd();
        drawOpt = "E1 SAME";
        if (i == 0) { // Draw first plot to set axes
            h_eff_hiHF[name]->GetYaxis()->SetRangeUser(0.0, 1.1);
            h_eff_hiHF[name]->GetYaxis()->SetTitle("Efficiency");
            drawOpt = "E1";
        }
        h_eff_hiHF[name]->Draw(drawOpt);
        leg_hiHF->AddEntry(h_eff_hiHF[name], name.c_str(), "lep");
    }

    // --- Finalize and Save Plots ---
    c_hiBin->cd();
    leg_hiBin->Draw();
    //c_hiBin->SaveAs("efficiency_vs_hiBin_all.pdf");

    c_hiHF->cd();
    leg_hiHF->Draw();
    //c_hiHF->SaveAs("efficiency_vs_hiHF_all.pdf");
    
    // --- Save everything to the ROOT file ---
    fout->cd();    
    // Write totals
    h_hiBin_total->Write();
    h_hiHF_total->Write();

    // Write all numerators and efficiencies
    for (const auto& name : filterNames) {
        h_hiBin_pass[name]->Write();
        h_hiHF_pass[name]->Write();
        h_eff_hiBin[name]->Write();
        h_eff_hiHF[name]->Write();
    }

    // Write canvases
    c_hiBin->Write();
    c_hiHF->Write();

    // --- Close the file ---
    fout->Close();

    cout << "All histograms and plots saved to " << fout->GetName() << endl;
}

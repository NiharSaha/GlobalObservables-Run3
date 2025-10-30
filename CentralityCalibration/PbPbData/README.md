Instructions to be followed:
1. You need to perform the calibration inside CMSSW, so make sure to use the relevant CMSSW version.
2. You need to change the arguments of the function according to PbPb run for the specific year. There will be changes to the HLT trigger, Run number, and coincFilter, along with your forest.txt file.
3. You should also consider changing the tag according to the run for clarity.
   For example:
   ```CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal```
   
   ```HFtowers200``` -> Tower Based
   
   ```periHYDJETshape``` -> Peripheral (below threshold) is taken from HYD MC
   
   ```run3v140x01``` -> contains Run3, CMSSW version
5. The output of the 2024 calibration will be two files: one is a root file containing all the histograms (for hiHF) and tree branch (for hiBin), and another txt file containing the calibration table.

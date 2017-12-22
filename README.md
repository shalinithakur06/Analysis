# Analysis
   
### Download and compile the package  ###  
* source /cvmfs/cms.cern.ch/cmsset_default.sh
* cmsrel CMSSW_8_0_28
* cd CMSSW_8_0_28/src
* cmsenv
* https://github.com/shalinithakur06/Analysis.git 
* cd Analysis/src
* make clean 
* make
* cd .. 

### Compile and run, in one go ### 
* root -l 'runMe.C("Analyzer")'

### Submit condor batch jobs  ###

* cd condor
* ./RunCond_TIFR.sh ntupleT2Paths.txt

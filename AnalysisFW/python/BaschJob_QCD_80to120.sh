# Lxplus Batch Job Script

cd /afs/cern.ch/user/b/bchazinq/work/CMSSW_5_3_32/src/
eval `scramv1 runtime  -sh`
cd 2011-jet-inclusivecrosssection-ntupleproduction-optimized/AnalysisFW/python/
cmsRun OpenDataTreeProducerOptimized_mcPAT_2011_cfg.py

#execute: 
# bsub -q 1nw -J job1 < BaschJob_QCD_80to120.sh  


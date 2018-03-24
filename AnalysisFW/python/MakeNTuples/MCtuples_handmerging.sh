# Helper script to merge tuples 

cwd=$(pwd)

# Environment
echo "CMSSW environment"
cd /afs/cern.ch/user/b/bchazinq/work/CMSSW_5_3_32
eval `scramv1 runtime -sh`
cd -

echo "EOS environment"
source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh

## Mount EOS
##mkdir eos
##eosmount eos/
##cd eos/cms/store/group/phys_smp/mhaapale/Jet/crab_OpenDataTree_all/160808_155332/0000


# Top-level directory of the N-tuples outputs
cd /eos/user/b/bchazinq/QCDPt15to30_00000

# Merge files
MC=output_*.root
#DATA=OpenDataTree_data_*.root
NUM=$(ls -l $MC | wc -l)
echo "Merging $NUM files..." 

hadd -f tuples_QCDPt15to30_00000.root $MC
#hadd -f tuples_QCDPt15to30_00001.root $MC
#hadd -f tuples_MC30-50.root $MC

cd $cwd
#eosumount eos/

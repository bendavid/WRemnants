combinetf_dir=$1
working_dir=$2

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $combinetf_dir
eval `scram runtime -sh`
popd

pushd $working_dir

shift 2
combinetf.py $@

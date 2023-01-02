export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ):$PYTHONPATH"
export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/narf:$PYTHONPATH"

WREM_BASE=""
if [[ "$PWD" == *"/WRemnants"* ]]; then
    WREM_BASE="$PWD"
elif [[ "$HOSTNAME" == *"mit.edu"* ]]; then
    WREM_BASE="${PWD}/WRemnants"
fi
echo "Created environment variable WREM_BASE=${WREM_BASE}"

COMBINE_STUDIES=""
if [[ "$HOSTNAME" == *"lxplus8s10.cern.ch"* ]]; then
    COMBINE_STUDIES="/scratch/${USER}/CombineStudies/"
#elif [[ "$HOSTNAME" == *"mit.edu"* ]]; then
#    COMBINE_STUDIES="/data/submit/cms/store/${USER}/CombineStudies/" # FIXME: EDIT PATH
else
    echo "Unknown hostname ${HOSTNAME}. Please update ${BASH_SOURCE}"
    exit 0
fi

# folders to store datacards and root files with histograms for combine
mkdir -pv ${COMBINE_STUDIES}/WMass/
mkdir -pv ${COMBINE_STUDIES}/ZMassWLike/
echo "Created environment variable COMBINE_STUDIES=${COMBINE_STUDIES}"

# web page folder to store plots
PLOTS="/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis/"
mkdir -pv ${PLOTS}
echo "Created environment variable PLOTS=${PLOTS}"

echo "Creating links to output folders"
ln -sv ${PLOTS} ${WREM_BASE}/scripts/analysisTools/plots
ln -sv ${COMBINE_STUDIES} ${WREM_BASE}/CombineStudies

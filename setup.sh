export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ):$PYTHONPATH"
export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/narf:$PYTHONPATH"

export WREM_BASE=""
if [[ "$PWD" == *"/WRemnants"* ]]; then
    WREM_BASE="$PWD"
else
    WREM_BASE="${PWD}/WRemnants"
fi
echo "Created environment variable WREM_BASE=${WREM_BASE}"

export COMBINE_STUDIES=""
if [[ "$HOSTNAME" == *"lxplus8s10.cern.ch"* ]]; then
    COMBINE_STUDIES="/scratch/${USER}/CombineStudies/"
elif [[ "$HOSTNAME" == *"mit.edu"* ]]; then
    COMBINE_STUDIES="/data/submit/${USER}/CombineStudies/"
elif [[ "$HOSTNAME" == *"cmsanalysis.pi.infn.it"* ]]; then
    COMBINE_STUDIES="/scratchnvme/${USER}/CombineStudies/"
else
    echo "Unknown hostname ${HOSTNAME}. Please update ${BASH_SOURCE}"
    exit 0
fi

# folders to store datacards and root files with histograms for combine
mkdir -pv ${COMBINE_STUDIES}/WMass/
mkdir -pv ${COMBINE_STUDIES}/ZMassWLike/
echo "Created environment variable COMBINE_STUDIES=${COMBINE_STUDIES}"

# web page folder to store plots
export PLOTS="/eos/user/${USER:0:1}/${USER}/www/WMassAnalysis/"
mkdir -pv ${PLOTS}
echo "Created environment variable PLOTS=${PLOTS}"

echo "Creating links to output folders (if not existing)"

plotLink="${WREM_BASE}/scripts/analysisTools/plots"
if [ -L ${plotLink} ] && [ -e ${plotLink} ]; then
    echo ">>> Existing link ${plotLink} --> ${PLOTS}"
else
    ln -sv ${PLOTS} ${plotLink}
fi

combineLink="${WREM_BASE}/CombineStudies"
if [ -L ${combineLink} ] && [ -e ${combineLink} ]; then
    echo ">>> Existing link ${combineLink} --> ${COMBINE_STUDIES}"
else
    ln -sv ${COMBINE_STUDIES} ${combineLink}
fi


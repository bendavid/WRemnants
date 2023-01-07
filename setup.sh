export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ):$PYTHONPATH"
export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/narf:$PYTHONPATH"

export WREM_BASE=""
if [[ "$PWD" == *"/WRemnants"* ]]; then
    WREM_BASE="$PWD"
else
    WREM_BASE="${PWD}/WRemnants"
fi
echo "Created environment variable WREM_BASE=${WREM_BASE}"

hostNameNotFound=false
export COMBINE_STUDIES="${WREM_BASE}/CombineStudies"
if [[ "$HOSTNAME" == *"lxplus8s10.cern.ch"* ]]; then
    COMBINE_STUDIES="/scratch/${USER}/CombineStudies/"
elif [[ "$HOSTNAME" == *"mit.edu"* ]]; then
    COMBINE_STUDIES="/data/submit/${USER}/CombineStudies/"
elif [[ "$HOSTNAME" == *"cmsanalysis.pi.infn.it"* ]]; then
    COMBINE_STUDIES="/scratchnvme/${USER}/CombineStudies/"
else
    echo ">>> Unknown hostname ${HOSTNAME}. Please update ${BASH_SOURCE} if you need a specific folder."
    echo ">>> \$COMBINE_STUDIES will point to default local folder ${COMBINE_STUDIES}"
    hostNameNotFound=true
fi
echo "Created environment variable COMBINE_STUDIES=${COMBINE_STUDIES}"

# folders to store datacards and root files with histograms for combine
mkdir -pv ${COMBINE_STUDIES}/WMass/
mkdir -pv ${COMBINE_STUDIES}/ZMassWLike/

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

# create link if $COMBINE_STUDIES was not a local folder
if [[ "$hostNameNotFound" == false ]]; then
    combineLink="${WREM_BASE}/CombineStudies"
    if [ -L ${combineLink} ] && [ -e ${combineLink} ]; then
        echo ">>> Existing link ${combineLink} --> ${COMBINE_STUDIES}"
    else
        ln -sv ${COMBINE_STUDIES} ${combineLink}
    fi
fi


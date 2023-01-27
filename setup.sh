export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ):$PYTHONPATH"
export PYTHONPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/narf:$PYTHONPATH"

export WREM_BASE=""
if [[ "$PWD" == *"/WRemnants"* ]]; then
    WREM_BASE="$PWD"
else
    WREM_BASE="${PWD}/WRemnants"
fi
echo "Created environment variable WREM_BASE=${WREM_BASE}"

# utility variables pointing to specific folders in the filesystem
export COMBINE_STUDIES="${WREM_BASE}/scripts/combine/"
echo "Created environment variable COMBINE_STUDIES=${COMBINE_STUDIES}"

export PLOTS="${WREM_BASE}/scripts/analysisTools/"
echo "Created environment variable PLOTS=${PLOTS}"

WREM_BASE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export WREM_BASE=$(readlink -f "$WREM_BASE")

source ${WREM_BASE}/narf/setup.sh

export PYTHONPATH="${WREM_BASE}:$PYTHONPATH"

echo "Created environment variable WREM_BASE=${WREM_BASE}"

# utility variables pointing to specific folders in the filesystem
export COMBINE_STUDIES="${WREM_BASE}/scripts/combine/"
echo "Created environment variable COMBINE_STUDIES=${COMBINE_STUDIES}"

export PLOTS="${WREM_BASE}/scripts/analysisTools/"
echo "Created environment variable PLOTS=${PLOTS}"

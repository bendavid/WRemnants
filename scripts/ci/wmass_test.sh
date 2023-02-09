SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

pushd $SCRIPT_DIR/../..
. ./setup.sh

python3 scripts/histmakers/mw_with_mu_eta_pt.py -j 8 --maxFiles 2

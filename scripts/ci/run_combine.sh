#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: run_combine.sh <combinetf_dir> <combine_cards>"
	exit 1
fi

combinetf_dir=$1

cards=()
for card in ${@:2}; do
	echo $card
	cards+=( `basename $card` )
done

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $combinetf_dir
eval `scram runtime -sh`
popd

working_dir=`dirname $2`
card_name=${cards[0]}

pushd $working_dir

if [[ $# -gt 2 ]]; then
	card_name=${working_dir}.txt
	combineCards.py $cards > $card_name
fi

outfile=${card_name/txt/hdf5}
text2hdf5.py --X-allow-no-signal $card_name
combinetf.py --doImpacts --binByBinStat -t -1 $outfile
popd
python3 scripts/combine/printImpacts -f $working_dir/$outfile

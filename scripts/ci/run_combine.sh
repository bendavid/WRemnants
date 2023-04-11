#!/bin/bash

if [[ $# -lt 3 ]]; then
	echo "Requires at least three arguments: run_combine.sh <combinetf_dir> <mode> <combine_cards>"
	exit 1
fi

combinetf_dir=$1

mode=$2

cards=()
for card in ${@:3}; do
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

if [ "$mode" == "mass" ]; then
	text2hdf5.py --X-allow-no-signal "$card_name"
	combinetf.py --doImpacts --binByBinStat -t -1 "$outfile"
elif [ "$mode" == "unfold" ]; then
	text2hdf5.py --X-allow-no-background --maskedChan=xnorm "$card_name"
	combinetf.py --doImpacts --binByBinStat -t 1 --saveHists --computeHistErrors --correlateXsecStat "$outfile"
fi

#!/bin/bash

if [[ $# -lt 4 ]]; then
	echo "Requires at least four arguments: run_combine.sh <combinetf_dir> <mode> <working_dir> <combine_cards>"
	exit 1
fi

combinetf_dir=$1
mode=$2
working_dir=$3

cards=()
for card in ${@:4}; do
	echo $card
	cards+=( $card )
done

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $combinetf_dir
eval `scram runtime -sh`
popd

pushd $working_dir

card_name=${working_dir}.txt
echo "combineCards.py ${cards[@]} > $card_name"
combineCards.py ${cards[@]} > $card_name

outfile=${card_name/txt/hdf5}

if [ "$mode" == "mass" ]; then
	echo text2hdf5.py --X-allow-no-signal "$card_name"
	text2hdf5.py --X-allow-no-signal "$card_name"
	echo combinetf.py --doImpacts --binByBinStat -t -1 "$outfile"
	combinetf.py --doImpacts --binByBinStat -t -1 "$outfile"
elif [ "$mode" == "unfold" ]; then
	echo text2hdf5.py --X-allow-no-background --maskedChan=xnorm "$card_name"
	text2hdf5.py --X-allow-no-background --maskedChan=xnorm "$card_name"
	echo combinetf.py --doImpacts --binByBinStat -t -1 --correlateXsecStat "$outfile"
	combinetf.py --doImpacts --binByBinStat -t -1 --correlateXsecStat "$outfile"
fi
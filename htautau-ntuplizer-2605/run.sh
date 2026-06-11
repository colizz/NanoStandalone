#!/bin/bash

set -e

INPUT_FILE_LIST=$1
OUTPUT_EOS_PATH_FATJET=$2
OUTPUT_EOS_PATH_BOOSTEDTAU=$3
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

echo "INPUT_FILE_LIST: $INPUT_FILE_LIST"

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
    source /cvmfs/cms.cern.ch/cmsset_default.sh
fi

export RELEASE=CMSSW_15_0_19
if [ -r "$RELEASE/src" ]; then
    echo "release $RELEASE already exists"
else
    scram p CMSSW "$RELEASE"
fi

cd "$RELEASE/src"
eval `scram runtime -sh`
cd ../..

IFS='|' read -r -a INPUT_FILES <<< "$INPUT_FILE_LIST"

FATJET_SELECTION_RULE="max_pt"
case "$OUTPUT_EOS_PATH_FATJET" in
    *htautau_pt200_file*.root)
        FATJET_SELECTION_RULE="gen_matched"
        ;;
esac
echo "FatJet object selection rule: ${FATJET_SELECTION_RULE}"

BOOSTEDTAU_SELECTION_RULE="max_pt_higgs"
case "$OUTPUT_EOS_PATH_BOOSTEDTAU" in
    *htautau_pt200_file*.root)
        BOOSTEDTAU_SELECTION_RULE="gen_matched"
        ;;
esac
echo "BoostedTau object selection rule: ${BOOSTEDTAU_SELECTION_RULE}"

FATJET_OUTPUTS=()
BOOSTEDTAU_OUTPUTS=()

for i in "${!INPUT_FILES[@]}"; do
    INPUT_FILE="${INPUT_FILES[$i]}"
    FATJET_OUTPUT="output_fatjet_${i}.root"
    BOOSTEDTAU_OUTPUT="output_boostedtau_${i}.root"

    echo "========================================="
    echo "Processing input ${i}: ${INPUT_FILE}"
    echo "FatJet output: ${FATJET_OUTPUT}"
    echo "BoostedTau output: ${BOOSTEDTAU_OUTPUT}"
    echo "========================================="

    root -l -b -q "$SCRIPT_DIR/nanoAOD_to_htautau_fatjet.C++(\"$INPUT_FILE\", \"$FATJET_OUTPUT\", \"$FATJET_SELECTION_RULE\")"
    root -l -b -q "$SCRIPT_DIR/nanoAOD_to_htautau_boostedTau.C++(\"$INPUT_FILE\", \"$BOOSTEDTAU_OUTPUT\", \"$BOOSTEDTAU_SELECTION_RULE\")"

    FATJET_OUTPUTS+=("$FATJET_OUTPUT")
    BOOSTEDTAU_OUTPUTS+=("$BOOSTEDTAU_OUTPUT")
done

hadd -f "output_fatjet.root" "${FATJET_OUTPUTS[@]}"
hadd -f "output_boostedtau.root" "${BOOSTEDTAU_OUTPUTS[@]}"

xrdcp -f "output_fatjet.root" "$OUTPUT_EOS_PATH_FATJET"
xrdcp -f "output_boostedtau.root" "$OUTPUT_EOS_PATH_BOOSTEDTAU"

touch dummy.cc

#!/bin/bash

set -e

INPUT_FILE_LIST=$1
OUTPUT_EOS_PATH_FATJET=$2
MAX_EVENTS=${3:--1}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

echo "INPUT_FILE_LIST: $INPUT_FILE_LIST"
echo "OUTPUT_EOS_PATH_FATJET: $OUTPUT_EOS_PATH_FATJET"
echo "MAX_EVENTS: $MAX_EVENTS"

if [ -z "$INPUT_FILE_LIST" ] || [ -z "$OUTPUT_EOS_PATH_FATJET" ]; then
    echo "Usage: $0 INPUT_FILE_LIST OUTPUT_EOS_PATH_FATJET [MAX_EVENTS]" >&2
    exit 1
fi

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

if [ ! -r "$SCRIPT_DIR/NanoTuples_nanov15_custom.tgz" ]; then
    echo "Missing NanoTuples package: $SCRIPT_DIR/NanoTuples_nanov15_custom.tgz" >&2
    exit 1
fi

tar -xzf "$SCRIPT_DIR/NanoTuples_nanov15_custom.tgz"

models=(
    "InclParticleTransformer-MD/ak8/V01/model.onnx"
    "InclParticleTransformer-MD/ak8/V02-HidLayer/model_embed.onnx"
)

for path in "${models[@]}"; do
    target="$CMSSW_BASE/src/PhysicsTools/NanoTuples/data/$path"
    if [ -s "$target" ]; then
        echo "Model already exists: $target"
        continue
    fi
    mkdir -p "$(dirname "$target")"
    echo "Downloading model: $path"
    wget -q "https://coli.web.cern.ch/coli/repo/NanoTuples_data/$path" -O "$target"
done

scram b -j8
cd ../..

IFS='|' read -r -a INPUT_FILES <<< "$INPUT_FILE_LIST"

FATJET_SELECTION_RULE="max_pt"
case "$OUTPUT_EOS_PATH_FATJET" in
    *htautau_pt200_file*.root)
        FATJET_SELECTION_RULE="gen_matched"
        ;;
esac
echo "FatJet object selection rule: ${FATJET_SELECTION_RULE}"

FATJET_OUTPUTS=()

for i in "${!INPUT_FILES[@]}"; do
    INPUT_FILE="${INPUT_FILES[$i]}"
    CMS_INPUT_FILE="$INPUT_FILE"
    case "$INPUT_FILE" in
        file:*|root://*|/store/*)
            ;;
        /*)
            CMS_INPUT_FILE="file:$INPUT_FILE"
            ;;
    esac
    NANO_OUTPUT="custom_nano_${i}.root"
    FATJET_OUTPUT="output_fatjet_${i}.root"
    CMSDRIVER_CFG="custom_nano_${i}_cfg.py"

    echo "========================================="
    echo "Processing MiniAOD input ${i}: ${INPUT_FILE}"
    echo "cmsRun filein: ${CMS_INPUT_FILE}"
    echo "Custom NanoAOD output: ${NANO_OUTPUT}"
    echo "FatJet output: ${FATJET_OUTPUT}"
    echo "========================================="

    cmsDriver.py \
        --python_filename "$CMSDRIVER_CFG" \
        --eventcontent NANOAODSIM \
        --customise PhysicsTools/NanoTuples/nanoTuples_cff.nanoTuples_customizeMC \
        --datatier NANOAODSIM \
        --fileout "file:$NANO_OUTPUT" \
        --conditions 150X_mcRun3_2024_realistic_v2 \
        --step NANO \
        --scenario pp \
        --filein "$CMS_INPUT_FILE" \
        --era Run3_2024 \
        --mc \
        -n "$MAX_EVENTS" \
        --no_exec

    cmsRun "$CMSDRIVER_CFG"

    root -l -b -q "$SCRIPT_DIR/nanoAOD_to_htautau_fatjet_custom_tagger.C++(\"$NANO_OUTPUT\", \"$FATJET_OUTPUT\", \"$FATJET_SELECTION_RULE\")"

    FATJET_OUTPUTS+=("$FATJET_OUTPUT")
done

hadd -f "output_fatjet.root" "${FATJET_OUTPUTS[@]}"

case "$OUTPUT_EOS_PATH_FATJET" in
    root://*)
        xrdcp -f "output_fatjet.root" "$OUTPUT_EOS_PATH_FATJET"
        ;;
    *)
        cp "output_fatjet.root" "$OUTPUT_EOS_PATH_FATJET"
        ;;
esac

touch dummy.cc

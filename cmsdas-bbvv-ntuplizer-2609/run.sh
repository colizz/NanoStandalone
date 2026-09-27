#!/usr/bin/env bash

set -euo pipefail

INPUT_FILE_LIST=$1
OUTPUT_EOS_PATH=$2
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# Same worker setup as htautau-ntuplizer-2605. Override only for local ROOT tests.
if [[ ${CMSDAS_SKIP_CMSSW:-0} != 1 ]]; then
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export RELEASE=CMSSW_15_0_19
    if [[ ! -d $RELEASE/src ]]; then
        scram p CMSSW "$RELEASE"
    fi
    cd "$RELEASE/src"
    eval "$(scram runtime -sh)"
    cd ../..
fi

IFS='|' read -r -a INPUT_FILES <<< "$INPUT_FILE_LIST"
[[ ${#INPUT_FILES[@]} -gt 0 ]]
OUTPUTS=()

export CMSDAS_MACRO="$SCRIPT_DIR/nanoAOD_to_cmsdas_bbvv.C"

for i in "${!INPUT_FILES[@]}"; do
    input=${INPUT_FILES[$i]}
    [[ -n $input ]]
    
    if [[ $input == /store/* ]]; then
        input="${XROOTD_REDIRECTOR:-root://cms-xrd-global.cern.ch/}/${input}"
    fi

    export NANO_INPUT="$input" NANO_OUTPUT="output_${i}.root"

    root -l -b <<'ROOT'
if (gSystem->CompileMacro(gSystem->Getenv("CMSDAS_MACRO"), "kO") == 0) {
    gSystem->Exit(1);
}

nanoAOD_to_cmsdas_bbvv(gSystem->Getenv("NANO_INPUT"), gSystem->Getenv("NANO_OUTPUT"));
gSystem->Exit(0);
ROOT

    [[ -s $NANO_OUTPUT ]]
    OUTPUTS+=("$NANO_OUTPUT")
done

hadd -f output.root "${OUTPUTS[@]}"

# EOS transfer follows the reference workflow. Local destinations support testing.
if [[ $OUTPUT_EOS_PATH == root://* ]]; then
    xrdcp -f output.root "$OUTPUT_EOS_PATH"
else
    cp output.root "$OUTPUT_EOS_PATH"
fi

touch dummy.cc

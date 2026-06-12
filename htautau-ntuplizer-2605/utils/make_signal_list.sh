#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_dir="$(cd "${script_dir}/.." && pwd)"

input_dir="/eos/cms/store/cmst3/group/vhcc/sfTuples/ext/GluGluH-Hto2Tau_Bin-PT-200_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8/24MiniAODv6"
output="${base_dir}/samples/htautau_pt200_mini.txt"

find "${input_dir}" \
    -maxdepth 1 -name 'miniv6*.root' -type f -printf '%T@ %p\n' \
| sort -n \
| cut -d' ' -f2- \
| awk '
{
    sub("^/eos/cms", "root://eoscms.cern.ch/", $0)
    files[n++] = $0
}
END {
    for (i = 0; i < n; i += 10) {
        line = files[i]
        for (j = i + 1; j < i + 10 && j < n; j++) {
            line = line "|" files[j]
        }
        printf "%s htautau_pt200_file%d\n", line, i/10
    }
}' > "${output}"

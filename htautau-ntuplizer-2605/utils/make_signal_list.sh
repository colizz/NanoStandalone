#!/bin/bash

set -e

find /eos/cms/store/cmst3/group/vhcc/sfTuples/ext/GluGluH-Hto2Tau_Bin-PT-200_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-pythia8/24MiniAODv6 \
    -maxdepth 1 -name 'nano*.root' -type f -printf '%T@ %p\n' \
| sort -n \
| cut -d' ' -f2- \
| awk '
{
    files[n++] = $0
}
END {
    for (i = 0; i < n; i += 10) {
        line = files[i]
        for (j = i + 1; j < i + 10 && j < n; j++) {
            line = line "," files[j]
        }
        printf "%s htautau_pt200_file%d\n", line, i/10
    }
}' > htautau_pt200.txt
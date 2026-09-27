# CMSDAS bbVV NanoAOD ntuplizer (2609)

The layout and worker workflow follow `htautau-ntuplizer-2605`: the macro and
`run.sh` are at the top level, with `jdl/`, `samples/`, and `log/` directories.

## Input lists and submission

Use `samples/hhbbvv.txt` and `samples/ttbar.txt` as single-column input lists,
with one NanoAOD path per line. For example:

```text
/store/.../file1.root
/store/.../file2.root
```

Full XRootD URLs are also supported. Each row creates one job. A row may also
contain multiple paths joined by `|`; the worker processes and merges that
group. Do not add a second output-name column, duplicate files, or blank rows.
Condor generates output names automatically using the sample name, Cluster,
and Process identifiers.

Run these commands from this directory with a valid X509 proxy:

```bash
mkdir -p log
condor_submit jdl/submit_hhbbvv.jdl
condor_submit jdl/submit_ttbar.jdl
```

## Worker workflow and output

The worker loads the CMS environment, creates or reuses `CMSSW_15_0_19`, and
loads its ROOT environment with `scram runtime`. It then processes each input
file with the compiled `nanoAOD_to_cmsdas_bbvv.C` macro, merges the job outputs
with `hadd`, and uploads the merged file with `xrdcp`.

Both JDL files use this output directory:

```text
root://eoscms.cern.ch//store/cmst3/group/vhcc/sfTuples/trees_sf/20260927_NanoV15_standalone_cmsdas_bbvv_ntuples
```

Output names are `ntuple_hhbbvv_<Cluster>_<Process>.root` and
`ntuple_ttbar_<Cluster>_<Process>.root`.
Ensure that the destination directory exists and is writable before production.
No remote directory was created and no remote upload was performed here.

A successful upload creates `dummy.cc`; Condor returns this marker and the job
logs. Failures return a nonzero exit code. As in the reference workflow,
`xrdcp -f` overwrites an existing output with the same name. A fresh submission gets a new Cluster identifier and new filenames. Merge only
one successful output per input row; do not combine duplicate inputs across
resubmitted batches. Use a separate `EOSDIR` for each production batch if needed.

For `/store/...` inputs, the default redirector is
`root://cms-xrd-global.cern.ch/`. Override it with `XROOTD_REDIRECTOR` if needed.

## Selection and ntuple contents

The macro uses `Jet_btagUParTAK4B > 0.1272` (2024 medium WP) and requires at least
two selected b-tagged jets and two selected leptons. No bbWW truth filter is
applied to the HH 2B2V signal.

- Jets: `pt > 15 GeV`, `|eta| < 2.5`.
- Electrons: `pt > 10 GeV`, `|eta| < 2.5`, `pfRelIso03_all < 0.12`.
- Muons: `pt > 10 GeV`, `|eta| < 2.4`, `pfRelIso04_all < 0.25`.

There are no additional ID, trigger, opposite-charge, or overlap-removal cuts.
The main `tree` contains the same 15 branches and types as the reference
ntuples. The aligned `audit` tree stores event UID and `genWeight`; auxiliary
objects retain cutflow, weight sums, input provenance, and selection settings.
No cross-section, luminosity, scale-factor, or systematic corrections are
applied. NanoAOD isolation and CMS b-tagging differ from the Delphes definitions.

## Local validation

For local tests, load an available ROOT environment and set
`CMSDAS_SKIP_CMSSW=1`. Run `run.sh` from a separate working directory, providing
an input file (or a `|`-joined list) and a local output path. This skips only the
CMSSW setup; it still runs the macro, `hadd`, and local output copy.

The wrapper produced 78 selected events from the provided 313-event NanoAOD.
All output branches agree with the previously validated implementation.
Remote EOS upload, target CMSSW initialization, and Condor scheduling have not
been tested.

# AVITI run summariser
`aviti_run_summariser.py` aggregates quality and performance metrics from demultiplexed [Element AVITI](https://www.elementbiosciences.com/products/aviti) sequencing runs processed by [Bases2Fastq](https://docs.elembio.io/docs/bases2fastq/). 

Given a parent directory containing one or more analysed run folders, the script parses `RunStats.json`, `RunParameters.json`, `Metrics.csv`, `UnassignedSequences.csv`, and per-sample `*_stats.json` files, then compiles two summary TSVs: a run-level summary spreadsheet and a sample-level summary spreadsheet. This provides a single point of reference for tracking run quality, yield, and assignment rates across multiple sequencing runs without needing to open individual QC reports or manually copy metrics.

This repository contains the python script `aviti_run_summariser.py`, explanations of parsed data and their source files, and guidance on running the script. 

> The AVITI run summariser is set up (via a crontab job) to run every 24 hours at midnight. Using `.processed_runs`, the script will only run on (bases2fastq)-processed AVITI run directories that have been generated since execution of the script had run. It will append new rows to the run_summary and sample_summary TSV files for each new AVITI run directory and samples within said directory, respectively.

---

## Requirements
 - Python ≥ 3.10
- No external dependencies (standard library only)
 
## Usage
If run manually:
```bash
python aviti_run_summariser.py -i /path/to/ANALYSED_RUNS -o /path/to/output_dir [--full]


-i: The input directory should contain one or more AVITI run subdirectories, each produced by Bases2Fastq. A valid run directory is identified by the presence of `RunStats.json` at its root and a completion flag in 'info/Bases2Fastq.log'
-o: The directory path where all outputs will be saved. Any directory that does not exist at runtime will be automatically created by the script.
--full: (optionally) Reprocess all runs and overwrite existing TSV files. Without this flag, only new (unprocessed) runs are appended to the existing outputs. Whether a run is new (unprocessed) or old (processed) is determined by the presence of the run directory name in a .processed_runs file within the output directory.
```

**The script will skip any subdirectories missing this file.**
 
## Outputs
| File | Description |
|------|-------------|
| `run_summary.tsv` | One row per run. Combines run configuration (from `RunParameters.json`), run-level quality metrics (from `RunStats.json` and `Metrics.csv`), and unassigned polony counts. |
| `sample_summary.tsv` | One row per sample. Includes sample-level yield, quality, and mismatch metrics parsed from `{SampleName}_stats.json`, along with run and project identifiers. |
 
---

## Column Reference
Element Bioscience definitions of the statistics reported in both files (below) are also outlined on their [website](https://docs.elembio.io/docs/elembio-cloud/run-output/).
### `run_summary.tsv`

| Column | Source | Description |
|--------|--------|-------------|
| `RunDirName` | Directory name | The on-disk directory name for the run (e.g. `20260401_AV251604_3571_OP_ANALYSIS`). Included as a fallback identifier since it typically encodes date, instrument, and run number. |
| `RunStats_RunName` | `RunStats.json` | The human-assigned run name set when configuring the sequencing run on the instrument. |
| `RunParameters_Date` | `RunParameters.json` | ISO 8601 timestamp of when the run started on the instrument. |
| `RunParameters_InstrumentName` | `RunParameters.json` | The instrument identifier (e.g. `AV251604`). |
| `RunParameters_RunDescription` | `RunParameters.json` | Free-text description entered at run setup. |
| `NumSamplesFound` | Derived | Count of `*_stats.json` files found under `Samples/` for this run. As such, this number represents all samples, and does not take into account replicates across lanes (e.g., if 48 samples were run on lane 1 and the same 48 on lane 2, NumSamplesFound == 96) |
| `RunParameters_ThroughputSelection` | `RunParameters.json` | `High` or `Standard`. Determines the target polony density on the flow cell. |
| `RunParameters_KitConfiguration` | `RunParameters.json` | The sequencing kit used, expressed as total available cycles (e.g. `150Cycles`, `300Cycles`). Constrains the maximum combined read length (R1+R2+I1+I2). |
| `RunParameters_PreparationWorkflow` | `RunParameters.json` | Library prep origin. `ThirdParty` = non-Element kits & conversion kits used. |
| `RunParameters_ChemistryVersion` | `RunParameters.json` | The sequencing chemistry version (e.g. `CloudbreakFS`). |
| `RunParameters_LowDiversity` | `RunParameters.json` | `True`/`False`. Whether the run was flagged as low-diversity at setup. When `True`, base-calling is adjusted for libraries with unbalanced base composition (e.g. amplicon pools). |
| `RunParameters_LibraryType` | `RunParameters.json` | `Linear` or `Circular`. Indicates the library topology. |
| `RunParameters_PlatformVersion` | `RunParameters.json` | The AVITI software/firmware version at the time of sequencing. |
| `RunParameters_Cycles_R1` | `RunParameters.json` | Number of cycles allocated to Read 1. This is the *requested* cycle count; actual usable read length after trimming may differ (see `RunStats_MeanReadLength`). |
| `RunParameters_Cycles_R2` | `RunParameters.json` | Number of cycles allocated to Read 2. |
| `RunParameters_Cycles_I1` | `RunParameters.json` | Number of cycles allocated to Index Read 1. |
| `RunParameters_Cycles_I2` | `RunParameters.json` | Number of cycles allocated to Index Read 2. |
| `RunStats_TotalYield` | `RunStats.json` | Total yield in gigabases from all reads (assigned + unassigned). The difference from `AssignedYield` represents data lost to unrecognised index combinations. |
| `Metrics_TotalYield (gb)_{lane}` | `Metrics.csv` | Total yield in gigabases. Per-lane as above. |
| `RunStats_AssignedYield` | `RunStats.json` | Total yield in gigabases from reads successfully assigned to a sample. This is the usable data output. |
| `Metrics_AssignedYield (gb)_{lane}` | `Metrics.csv` | Assigned yield in gigabases. Per-lane as above. |
| `RunStats_NumPolonies` | `RunStats.json` → `Lanes[]` | Total polony count (summed across lanes). The AVITI equivalent of Illumina's cluster count, or the raw number of sequenceable molecules on the flow cell. Also output per-lane as `_L1`, `_L2`. |
| `Unassigned_TotalCount_1+2` | `UnassignedSequences.csv` | Sum of all unassigned polony counts from the combined-lane (`1+2`) rows. Quantifies the absolute scale of data that could not be assigned to any sample. |
| `RunStats_PercentAssignedReads` | `RunStats.json` → `Lanes[]` | Percentage of reads assigned to a sample based on index matching. Low values (e.g. < 90%) suggest index issues, contamination, or incorrect run manifest entries. Also output per-lane. |
| `Metrics_PercentAssigned_{lane}` | `Metrics.csv` | Percentage of reads assigned. Per-lane as above. |
| `RunStats_PercentQ30` | `RunStats.json` → `Lanes[]` | Percentage of base calls at Q30+ (≤1 in 1,000 error probability). Excludes filtered reads and no-calls. Also output per-lane. |
| `Metrics_PercentQ30_{lane}` | `Metrics.csv` | Percentage of Q30+ base calls from the Bases2Fastq CSV summary. Suffixed `_1+2` (combined), `_1`, and `_2` for per-lane breakdown. |
| `RunStats_PercentQ40` | `RunStats.json` → `Lanes[]` | Percentage of base calls at Q40+ (≤1 in 10,000 error probability). Excludes filtered reads and no-calls. Also output per-lane. |
| `Metrics_PercentQ40_{lane}` | `Metrics.csv` | Percentage of Q40+ base calls. Per-lane as above. |
| `RunStats_QualityScoreMean` | `RunStats.json` → `Lanes[]` | Arithmetic mean Q score across all base calls. A single-number quality summary; >35 is ideal. Also output per-lane. |
| `Metrics_QualityScoreMean_{lane}` | `Metrics.csv` | Mean Q score. Per-lane as above. |
| `RunStats_PercentMismatch` | `RunStats.json` → `Lanes[]` | Percentage of assigned polonies with at least one index mismatch. High rates may indicate index hopping, synthesis errors, or phasing issues. Also output per-lane. |
| `Metrics_PercentMismatch_{lane}` | `Metrics.csv` | Index mismatch rate. Per-lane as above. |
| `RunStats_PercentMismatchI1` | `RunStats.json` → `Lanes[]` | Mismatch rate for Index 1 specifically. Comparing I1 vs I2 helps diagnose asymmetric index quality issues. Also output per-lane. |
| `RunStats_PercentMismatchI2` | `RunStats.json` → `Lanes[]` | Mismatch rate for Index 2. Also output per-lane. |
| `RunStats_PercentMismatchI*_{lane}` | `RunStats.json` → `Lanes[]` | Mismatch rate for Index 1 or 2, per-lane. |
| `Metrics_PercentUnexpectedIndexPairs_{lane}` | `Metrics.csv` | Unexpected index pair rate. Also output per-lane. |
| `RunStats_MeanReadLength` | `RunStats.json` → `Lanes[].Reads[]` | Average read length in bases after adapter trimming (averaged across R1/R2 per lane, then polony-weighted across lanes). If substantially shorter than requested cycles, many reads had adapters trimmed (i.e. had short inserts). As we do not typically conduct on-board trimming, this number should equal half the value in RunParameters_KitConfiguration. Also output per-lane. |
| `RunStats_MeanReadLength_{lane}` | `RunStats.json` → `Lanes[].Reads[]` | Average read length in bases after adapter trimming. Per-lane as above. |
| `RunStats_PercentReadsTrimmed` | `RunStats.json` → `Lanes[].Reads[]` | Percentage of reads with adapter sequence removed. High values are expected for short-insert libraries; unexpected high trimming on long-insert libraries suggests a library prep issue. As we do not typically conduct on-board trimming, this column should be blank/zero for all runs. Also output per-lane. |
| `RunStats_FlowCellID` | `RunStats.json` | The flow cell serial number. Ties back to the physical consumable. |
| `RunStats_AnalysisVersion` | `RunStats.json` | The Bases2Fastq version used for demultiplexing. |

---


### `sample_summary.tsv`

| Column | Source | Description |
|--------|--------|-------------|
| `RunDirName` | Directory name | The on-disk run directory name. Links this sample back to its parent run. |
| `RunStats_RunName` | `RunStats.json` | The run name from `RunStats.json`, providing a human-readable run identifier for each sample row. |
| `ProjectName` | Directory structure | The project name inferred from `Samples/{ProjectName}/{SampleName}/`. Empty if samples are stored directly under `Samples/` without a project subdirectory. |
| `SampleStats_SampleName` | `{SampleName}_stats.json` | The alphanumeric sample identifier as recorded in the stats file. |
| `SampleStats_SampleNumber` | `{SampleName}_stats.json` | The numeric sample identifier assigned by Bases2Fastq during demultiplexing. |
| `SampleStats_NumPolonies` | `{SampleName}_stats.json` | Total number of polonies assigned to this sample. The primary measure of sequencing depth for the sample. |
| `SampleStats_Yield` | `{SampleName}_stats.json` | Total bases for this sample in gigabases. Combined with `NumPolonies`, this gives the effective read length (Yield / NumPolonies × 10⁹). |
| `SampleStats_MeanReadLength` | `{SampleName}_stats.json` | Average read length in bases after adapter trimming for this sample. Shorter than the requested cycle count indicates adapter contamination or short insert sizes in this sample specifically. |
| `SampleStats_PercentQ30` | `{SampleName}_stats.json` | Percentage of base calls at Q30+ for this sample. Excludes filtered reads and no-calls. Compare against the run-level `PercentQ30` to identify samples with unusually low or high quality. |
| `SampleStats_PercentQ40` | `{SampleName}_stats.json` | Percentage of base calls at Q40+ for this sample. |
| `SampleStats_QualityScoreMean` | `{SampleName}_stats.json` | Mean Q score for this sample. Samples with substantially lower values than the run average may indicate library-specific issues (e.g. degraded input, poor library complexity). |
| `SampleStats_PercentMismatch` | `{SampleName}_stats.json` | Percentage of polonies assigned to this sample with an index mismatch. A sample with a notably higher mismatch rate than others in the same run may have index synthesis issues or be receiving misassigned reads from a similar index. |
| `SampleStats_PercentReadsTrimmed` | `{SampleName}_stats.json` | Percentage of reads trimmed for this sample. Comparing across samples within a run highlights which libraries had shorter inserts or more adapter read-through. |

## Authorship
This script was written by Dan Parsons @NHMUK.

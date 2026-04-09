#!/usr/bin/env python3
"""
aviti_run_summariser.py

Recursively scans an input directory containing AVITI (Bases2Fastq) sequencing
run subdirectories and compiles run-level and sample-level summary TSV files.

Usage:
    python aviti_run_summariser.py -i /path/to/ANALYSED_RUNS -o /path/to/output_dir

Outputs:
    - run_summary.tsv        : one row per run (RunStats.json + Metrics.csv + UnassignedSequences.csv)
    - sample_summary.tsv     : one row per sample (SampleName_stats.json)

Author: Dan Parsons @NHMUK / with Claude
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────

def safe_json_load(path: Path) -> dict | None:
    """Load a JSON file, returning None on failure."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("Could not read %s: %s", path, exc)
        return None


def find_run_dirs(root: Path) -> list[Path]:
    """
    Identify AVITI run directories under *root*.
    A valid run dir must contain RunStats.json at its top level.
    """
    run_dirs = sorted(
        d for d in root.iterdir()
        if d.is_dir() and (d / "RunStats.json").is_file()
    )
    if not run_dirs:
        log.warning("No run directories (containing RunStats.json) found under %s", root)
    return run_dirs


# ──────────────────────────────────────────────
# Metrics.csv parsing
# ──────────────────────────────────────────────

METRICS_CSV_FIELDS = [
    "PercentQ30",
    "PercentQ40",
    "QualityScoreMean",
    "PercentAssigned",
    "PercentMismatch",
    "PercentUnexpectedIndexPairs",
    "AssignedYield (gb)",
    "TotalYield (gb)",
]


def parse_metrics_csv(path: Path) -> dict:
    """
    Parse Metrics.csv and return a flat dict with keys like:
        Metrics_PercentQ30_1+2, Metrics_PercentQ30_1, Metrics_PercentQ30_2, ...
    Uses utf-8-sig encoding to handle BOM (byte order mark).
    """
    result = {}
    if not path.is_file():
        log.warning("Metrics.csv not found: %s", path)
        return result
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames:
                reader.fieldnames = [f.strip() for f in reader.fieldnames]
                log.debug("Metrics.csv headers: %s", reader.fieldnames)
            else:
                log.warning("Metrics.csv has no headers: %s", path)
                return result
            for row in reader:
                metric = row.get("Metric", "").strip()
                lane = row.get("Lane", "").strip()
                value = row.get("Value", "").strip()
                if metric and lane:
                    if metric == "PercentQ50":
                        continue
                    key = f"Metrics_{metric}_{lane}"
                    result[key] = value
    except Exception as exc:
        log.warning("Could not read %s: %s", path, exc)
    return result


# ──────────────────────────────────────────────
# UnassignedSequences.csv parsing
# ──────────────────────────────────────────────

def parse_unassigned_csv(path: Path) -> int:
    """
    Parse UnassignedSequences.csv and return the total unassigned polony
    count from the combined lane rows (Lane == '1+2').
    """
    total = 0
    if not path.is_file():
        log.warning("UnassignedSequences.csv not found: %s", path)
        return total
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            reader.fieldnames = [f.strip() for f in reader.fieldnames]
            for row in reader:
                lane = row.get("Lane", "").strip()
                if lane == "1+2":
                    try:
                        total += int(row.get("Count", 0))
                    except (ValueError, TypeError):
                        pass
    except OSError as exc:
        log.warning("Could not read %s: %s", path, exc)
    return total


# ──────────────────────────────────────────────
# RunParameters.json parsing
# ──────────────────────────────────────────────

RUNPARAMS_FIELDS = [
    "Date",
    "InstrumentName",
    "RunDescription",
    "ThroughputSelection",
    "KitConfiguration",
    "PreparationWorkflow",
    "ChemistryVersion",
    "LowDiversity",
    "LibraryType",
    "PlatformVersion",
]

# Cycles sub-object keys to flatten
RUNPARAMS_CYCLES = ["R1", "R2", "I1", "I2"]


def parse_runparameters(path: Path) -> dict:
    """
    Parse RunParameters.json and return a flat dict with keys prefixed
    RunParameters_.
    """
    result = {}
    data = safe_json_load(path)
    if data is None:
        log.warning("RunParameters.json not found or unreadable: %s", path)
        return result
    for field in RUNPARAMS_FIELDS:
        result[f"RunParameters_{field}"] = data.get(field, "")
    # Flatten Cycles sub-object
    cycles = data.get("Cycles", {})
    for key in RUNPARAMS_CYCLES:
        result[f"RunParameters_Cycles_{key}"] = cycles.get(key, "")
    return result


# ──────────────────────────────────────────────
# RunStats.json parsing
# ──────────────────────────────────────────────

# Fields that are reliably at the top level of RunStats.json
RUNSTATS_TOP_LEVEL_FIELDS = [
    "RunName",
    "FlowCellID",
    "AnalysisVersion",
    "AssignedYield",
    "TotalYield",
]

# Fields that may be top-level OR only present inside Lanes[].
# If top-level, use directly; otherwise aggregate across lanes
# (polony-weighted average for percentages, sum for counts/yield).
RUNSTATS_LANE_FIELDS = [
    "NumPolonies",
    "NumPoloniesBeforeTrimming",
    "PercentAssignedReads",
    "PercentQ30",
    "PercentQ40",
    "QualityScoreMean",
    "PercentMismatch",
    "PercentMismatchI1",
    "PercentMismatchI2",
    "PercentUnexpectedIndexPairs",
]

# Fields nested inside Lanes[].Reads[] — need to be aggregated up
# from the Reads array first, then across lanes.
RUNSTATS_READS_FIELDS = [
    "MeanReadLength",
    "PercentReadsTrimmed",
]

# Fields to sum across lanes (counts / yield); everything else is
# weighted-averaged by NumPolonies.
_SUM_FIELDS = {"NumPolonies", "NumPoloniesBeforeTrimming"}


def _aggregate_lane_field(lanes: list[dict], field: str) -> str:
    """
    Aggregate a field across lanes.
    - For count fields (NumPolonies, etc.): sum.
    - For percentage / mean fields: polony-weighted average.
    Returns "" if the field is absent from all lanes.
    """
    values = []
    weights = []
    for lane in lanes:
        val = lane.get(field)
        if val is None:
            continue
        values.append(val)
        weights.append(lane.get("NumPolonies", 0))

    if not values:
        return ""

    if field in _SUM_FIELDS:
        return sum(values)

    # Weighted average
    total_weight = sum(weights)
    if total_weight == 0:
        return sum(values) / len(values)
    return sum(v * w for v, w in zip(values, weights)) / total_weight


def _get_reads_field_for_lane(lane: dict, field: str):
    """
    Extract a field from Lanes[].Reads[] and average across reads
    (e.g. R1 and R2) for a single lane.
    Returns None if not found.
    """
    reads = lane.get("Reads", [])
    values = []
    for read in reads:
        val = read.get(field)
        if val is not None:
            values.append(val)
    if not values:
        return None
    return sum(values) / len(values)


def _aggregate_reads_field(lanes: list[dict], field: str) -> str:
    """
    Aggregate a Reads-level field across lanes (polony-weighted average
    of the per-lane read-averaged values).
    """
    values = []
    weights = []
    for lane in lanes:
        val = _get_reads_field_for_lane(lane, field)
        if val is None:
            continue
        values.append(val)
        weights.append(lane.get("NumPolonies", 0))

    if not values:
        return ""

    total_weight = sum(weights)
    if total_weight == 0:
        return sum(values) / len(values)
    return sum(v * w for v, w in zip(values, weights)) / total_weight


def _recursive_search(obj, field: str):
    """
    Recursively search a nested dict/list structure for a key.
    Returns the first value found, or None.
    """
    if isinstance(obj, dict):
        if field in obj:
            return obj[field]
        for v in obj.values():
            result = _recursive_search(v, field)
            if result is not None:
                return result
    elif isinstance(obj, list):
        for item in obj:
            result = _recursive_search(item, field)
            if result is not None:
                return result
    return None


def parse_runstats(data: dict) -> dict:
    """
    Extract selected fields from RunStats.json, prefixed with RunStats_.

    Fields may live at different nesting levels:
      - Top-level: RunName, FlowCellID, AssignedYield, TotalYield, etc.
      - Lanes[]:   PercentQ30, QualityScoreMean, NumPolonies, etc.
      - Lanes[].Reads[]: MeanReadLength, PercentReadsTrimmed
      - Possibly elsewhere: PercentBelowFilterThreshold

    Per-lane values are also emitted as RunStats_{field}_L{lane_number}.
    """
    result = {}

    # ── Top-level fields ──
    for field in RUNSTATS_TOP_LEVEL_FIELDS:
        result[f"RunStats_{field}"] = data.get(field, "")

    lanes = data.get("Lanes", [])

    # ── Direct lane-level fields ──
    for field in RUNSTATS_LANE_FIELDS:
        # Try top-level first
        top_val = data.get(field)
        if top_val is not None:
            result[f"RunStats_{field}"] = top_val
        else:
            result[f"RunStats_{field}"] = _aggregate_lane_field(lanes, field)

        # Per-lane values
        for lane in lanes:
            lane_num = lane.get("Lane", "?")
            lane_val = lane.get(field)
            if lane_val is not None:
                result[f"RunStats_{field}_L{lane_num}"] = lane_val

    # ── Reads-level fields (Lanes[].Reads[].field) ──
    for field in RUNSTATS_READS_FIELDS:
        # Try top-level first
        top_val = data.get(field)
        if top_val is not None:
            result[f"RunStats_{field}"] = top_val
        else:
            result[f"RunStats_{field}"] = _aggregate_reads_field(lanes, field)

        # Per-lane values (averaged across reads within each lane)
        for lane in lanes:
            lane_num = lane.get("Lane", "?")
            lane_val = _get_reads_field_for_lane(lane, field)
            if lane_val is not None:
                result[f"RunStats_{field}_L{lane_num}"] = lane_val

    # ── PercentBelowFilterThreshold — location varies, use recursive search ──
    top_val = data.get("PercentBelowFilterThreshold")
    if top_val is not None:
        result["RunStats_PercentBelowFilterThreshold"] = top_val
    else:
        # Try each lane
        for lane in lanes:
            lane_num = lane.get("Lane", "?")
            lane_val = lane.get("PercentBelowFilterThreshold")
            if lane_val is not None:
                result[f"RunStats_PercentBelowFilterThreshold_L{lane_num}"] = lane_val
        # If still no aggregate, search recursively for a top-level-ish value
        if "RunStats_PercentBelowFilterThreshold" not in result:
            found = _recursive_search(data, "PercentBelowFilterThreshold")
            if found is not None:
                result["RunStats_PercentBelowFilterThreshold"] = found
            else:
                result["RunStats_PercentBelowFilterThreshold"] = ""

    return result


# ──────────────────────────────────────────────
# Sample-level stats parsing
# ──────────────────────────────────────────────

# Fields at the top level of {SampleName}_stats.json
SAMPLE_TOP_FIELDS = [
    "SampleName",
    "SampleNumber",
    "NumPolonies",
    "Yield",
]

# Fields inside Occurrences[] — aggregated (polony-weighted) across occurrences
SAMPLE_OCC_FIELDS = [
    "PercentQ30",
    "PercentQ40",
    "QualityScoreMean",
    "PercentMismatch",
]

# Fields inside Occurrences[].Reads[] — averaged across reads, then across occurrences
SAMPLE_READS_FIELDS = [
    "MeanReadLength",
    "PercentReadsTrimmed",
]

# Explicit column order for sample_summary.tsv
SAMPLE_SUMMARY_COLUMNS = [
    "RunDirName",
    "RunStats_RunName",
    "ProjectName",
    "SampleStats_SampleName",
    "SampleStats_SampleNumber",
    "SampleStats_NumPolonies",
    "SampleStats_Yield",
    "SampleStats_MeanReadLength",
    "SampleStats_PercentQ30",
    "SampleStats_PercentQ40",
    "SampleStats_QualityScoreMean",
    "SampleStats_PercentMismatch",
    "SampleStats_PercentReadsTrimmed",
]


def find_sample_stats(samples_dir: Path) -> list[tuple[str, Path]]:
    """
    Recursively find all *_stats.json files under Samples/.
    Returns list of (project_name_or_empty, path) tuples.
    """
    results = []
    if not samples_dir.is_dir():
        return results
    for stats_file in sorted(samples_dir.rglob("*_stats.json")):
        # Determine project: Samples/<Project>/<Sample>/<file> vs Samples/<Sample>/<file>
        rel = stats_file.relative_to(samples_dir)
        parts = rel.parts  # e.g. ('iBOLP1', 'POLNB2852-19', 'POLNB2852-19_stats.json')
        if len(parts) == 3:
            project = parts[0]
        else:
            project = ""
        results.append((project, stats_file))
    return results


def _aggregate_occurrences(occurrences: list[dict], field: str):
    """Polony-weighted average of a field across Occurrences[]."""
    values, weights = [], []
    for occ in occurrences:
        val = occ.get(field)
        if val is not None:
            values.append(val)
            weights.append(occ.get("NumPolonies", 0))
    if not values:
        return ""
    total_w = sum(weights)
    if total_w == 0:
        return sum(values) / len(values)
    return sum(v * w for v, w in zip(values, weights)) / total_w


def _aggregate_occ_reads_field(occurrences: list[dict], field: str):
    """
    For fields nested in Occurrences[].Reads[] (e.g. MeanReadLength):
    average across reads within each occurrence, then polony-weighted
    average across occurrences.
    """
    values, weights = [], []
    for occ in occurrences:
        reads = occ.get("Reads", [])
        read_vals = [r[field] for r in reads if r.get(field) is not None]
        if read_vals:
            values.append(sum(read_vals) / len(read_vals))
            weights.append(occ.get("NumPolonies", 0))
    if not values:
        return ""
    total_w = sum(weights)
    if total_w == 0:
        return sum(values) / len(values)
    return sum(v * w for v, w in zip(values, weights)) / total_w


def parse_sample_stats(data: dict, project: str, run_name: str, run_dir_name: str) -> dict:
    """
    Extract selected fields from a sample stats JSON.

    Fields may live at different nesting levels:
      - Top-level: SampleName, NumPolonies, Yield, etc.
      - Occurrences[]: PercentQ30, QualityScoreMean, PercentMismatch, etc.
      - Occurrences[].Reads[]: MeanReadLength, PercentReadsTrimmed
    """
    result = {
        "RunDirName": run_dir_name,
        "RunStats_RunName": run_name,
        "ProjectName": project,
    }

    # Top-level fields
    for field in SAMPLE_TOP_FIELDS:
        result[f"SampleStats_{field}"] = data.get(field, "")

    occurrences = data.get("Occurrences", [])

    # Occurrence-level fields (check top-level first, then aggregate)
    for field in SAMPLE_OCC_FIELDS:
        top_val = data.get(field)
        if top_val is not None:
            result[f"SampleStats_{field}"] = top_val
        else:
            result[f"SampleStats_{field}"] = _aggregate_occurrences(occurrences, field)

    # Reads-level fields (check top-level first, then aggregate)
    for field in SAMPLE_READS_FIELDS:
        top_val = data.get(field)
        if top_val is not None:
            result[f"SampleStats_{field}"] = top_val
        else:
            result[f"SampleStats_{field}"] = _aggregate_occ_reads_field(occurrences, field)

    return result


# ──────────────────────────────────────────────
# Per-run aggregation
# ──────────────────────────────────────────────

def process_run(run_dir: Path) -> tuple[dict, list[dict]]:
    """
    Process a single run directory.
    Returns (run_row_dict, [sample_row_dicts]).
    """
    run_dir_name = run_dir.name
    log.info("Processing run: %s", run_dir_name)

    # ── RunStats.json ──
    runstats_data = safe_json_load(run_dir / "RunStats.json")
    if runstats_data is None:
        log.error("Skipping %s – RunStats.json unreadable", run_dir_name)
        return {}, []
    run_row = {"RunDirName": run_dir_name}
    run_row.update(parse_runstats(runstats_data))

    run_name = runstats_data.get("RunName", run_dir_name)

    # ── RunParameters.json (inserted early so columns appear after RunStats_RunName) ──
    runparams = parse_runparameters(run_dir / "RunParameters.json")
    # Insert RunParameters columns right after RunStats_RunName by rebuilding the dict
    early_keys = ["RunDirName", "RunStats_RunName"]
    reordered = {}
    for k in early_keys:
        if k in run_row:
            reordered[k] = run_row.pop(k)
    reordered.update(runparams)
    reordered.update(run_row)
    run_row = reordered

    # ── Metrics.csv ──
    metrics = parse_metrics_csv(run_dir / "Metrics.csv")
    log.info("Metrics dict has %d keys for %s", len(metrics), run_dir_name)
    run_row.update(metrics)

    # ── UnassignedSequences.csv ──
    unassigned_total = parse_unassigned_csv(run_dir / "UnassignedSequences.csv")
    run_row["Unassigned_TotalCount_1+2"] = unassigned_total

    # ── Sample stats ──
    sample_rows = []
    for project, stats_path in find_sample_stats(run_dir / "Samples"):
        sdata = safe_json_load(stats_path)
        if sdata is None:
            continue
        sample_rows.append(parse_sample_stats(sdata, project, run_name, run_dir_name))

    run_row["NumSamplesFound"] = len(sample_rows)

    return run_row, sample_rows


# ──────────────────────────────────────────────
# TSV writing
# ──────────────────────────────────────────────

# Explicit column order for run_summary.tsv.
# Columns not in this list are excluded from the output.
RUN_SUMMARY_COLUMNS = [
    "RunDirName",
    "RunStats_RunName",
    "RunParameters_Date",
    "RunParameters_InstrumentName",
    "RunParameters_RunDescription",
    "NumSamplesFound",
    "RunParameters_ThroughputSelection",
    "RunParameters_KitConfiguration",
    "RunParameters_PreparationWorkflow",
    "RunParameters_ChemistryVersion",
    "RunParameters_LowDiversity",
    "RunParameters_LibraryType",
    "RunParameters_PlatformVersion",
    "RunParameters_Cycles_R1",
    "RunParameters_Cycles_R2",
    "RunParameters_Cycles_I1",
    "RunParameters_Cycles_I2",
    "RunStats_TotalYield",
    "Metrics_TotalYield (gb)_1",
    "Metrics_TotalYield (gb)_2",
    "RunStats_AssignedYield",
    "Metrics_AssignedYield (gb)_1",
    "Metrics_AssignedYield (gb)_2",
    "RunStats_NumPolonies",
    "RunStats_NumPolonies_L1",
    "RunStats_NumPolonies_L2",
    "Unassigned_TotalCount_1+2",
    "RunStats_PercentAssignedReads",
    "RunStats_PercentAssignedReads_L1",
    "RunStats_PercentAssignedReads_L2",
    "RunStats_PercentQ30",
    "RunStats_PercentQ30_L1",
    "RunStats_PercentQ30_L2",
    "RunStats_PercentQ40",
    "RunStats_PercentQ40_L1",
    "RunStats_PercentQ40_L2",
    "RunStats_QualityScoreMean",
    "RunStats_QualityScoreMean_L1",
    "RunStats_QualityScoreMean_L2",
    "RunStats_PercentMismatch",
    "RunStats_PercentMismatch_L1",
    "RunStats_PercentMismatch_L2",
    "RunStats_PercentMismatchI1",
    "RunStats_PercentMismatchI2",
    "RunStats_PercentMismatchI1_L1",
    "RunStats_PercentMismatchI1_L2",
    "RunStats_PercentMismatchI2_L1",
    "RunStats_PercentMismatchI2_L2",
    "RunStats_PercentUnexpectedIndexPairs",
    "RunStats_PercentUnexpectedIndexPairs_L1",
    "RunStats_PercentUnexpectedIndexPairs_L2",
    "RunStats_MeanReadLength",
    "RunStats_MeanReadLength_L1",
    "RunStats_MeanReadLength_L2",
    "RunStats_PercentReadsTrimmed",
    "RunStats_FlowCellID",
    "RunStats_AnalysisVersion",
]


def write_tsv(rows: list[dict], path: Path, columns: list[str] | None = None):
    """
    Write a list of dicts to a TSV.

    If *columns* is provided, only those columns are written in that order.
    Otherwise, the union of all keys (insertion-ordered) is used.
    """
    if not rows:
        log.warning("No data to write to %s", path)
        return

    if columns is not None:
        # Deduplicate while preserving order
        seen_cols = set()
        fieldnames = []
        for c in columns:
            if c not in seen_cols:
                fieldnames.append(c)
                seen_cols.add(c)
    else:
        # Auto-detect: preserve insertion order from the first row,
        # then add any extra keys from subsequent rows
        fieldnames = list(rows[0].keys())
        seen = set(fieldnames)
        for row in rows[1:]:
            for k in row:
                if k not in seen:
                    fieldnames.append(k)
                    seen.add(k)

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    log.info("Wrote %d rows (%d columns) to %s", len(rows), len(fieldnames), path)


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate AVITI (Bases2Fastq) sequencing run statistics into summary TSVs."
    )
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        type=Path,
        help="Root directory containing AVITI run subdirectories.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        type=Path,
        help="Output directory for summary TSV files.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        log.error("Input directory does not exist: %s", input_dir)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    run_dirs = find_run_dirs(input_dir)
    log.info("Found %d run directories", len(run_dirs))

    all_run_rows = []
    all_sample_rows = []

    for rd in run_dirs:
        run_row, sample_rows = process_run(rd)
        if run_row:
            all_run_rows.append(run_row)
        all_sample_rows.extend(sample_rows)

    # Write outputs
    run_out = output_dir / "run_summary.tsv"
    sample_out = output_dir / "sample_summary.tsv"

    write_tsv(all_run_rows, run_out, columns=RUN_SUMMARY_COLUMNS)
    write_tsv(all_sample_rows, sample_out, columns=SAMPLE_SUMMARY_COLUMNS)

    log.info(
        "Done. %d runs, %d samples across all runs.",
        len(all_run_rows),
        len(all_sample_rows),
    )


if __name__ == "__main__":
    main()

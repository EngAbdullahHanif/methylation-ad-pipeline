#!/usr/bin/env python3
"""Run the methylation stages with dependency-aware caching and peak memory measurements."""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(".workflow-cache")
R_ENV = ("renv.lock", ".Rprofile", "renv/settings.json", "renv/activate.R")
RAW = "data/raw/GSE212682"
GEO_MATRIX = "data/metadata/GSE212682_series_matrix.txt.gz"
EPIC_MASK = "data/reference/EPIC.hg19.MASK_general.20180909.tsv.gz"
REPLICATION_DATA = "data/raw/GSE59685"


def results(names):
    return tuple(f"results/{name}" for name in names.split())


STAGES = {
    "download": {
        "command": ("Rscript", "scripts/00_download_data.R"),
        "inputs": (f"{RAW}/GSE212682_RAW.tar",),
        "outputs": (f"{RAW}/GSE212682_RAW.tar", f"{RAW}/idats/*.idat.gz", "data/metadata/GSE212682_series.txt"),
    },
    "qc": {
        "previous": "download", "command": ("Rscript", "scripts/01_import_and_qc.R"),
        "inputs": (GEO_MATRIX,),
        "outputs": ("data/metadata/sample_sheet_clean.tsv",) + results(
            "rg_set_raw.rds detection_p.rds qc_summary.rds sex_check.rds qc_summary.tsv sex_check.tsv "
            "qc_mean_detection_p.pdf qc_intensity.pdf"
        ),
    },
    "normalization": {
        "previous": "qc", "command": ("Rscript", "scripts/02_normalisation.R"),
        "inputs": (EPIC_MASK,),
        "outputs": results("methyl_set_noob.rds methyl_set_filtered.rds beta_filtered.rds probe_universe_before_mask.rds "
                           "normalisation_summary.tsv normalisation_beta_density.pdf probe_filter_counts.tsv m2_summary.tsv"),
    },
    "differential": {
        "previous": "normalization", "command": ("Rscript", "scripts/03_exploration_and_differential.R"),
        "outputs": results(
            "batch_overlap.tsv model_design.tsv m3_model_note.txt m3_variable_summary.tsv m3_results.rds "
            "pca_scores.tsv pca_variance.tsv pca_by_batch.pdf pca_by_diagnosis.pdf pc_variable_associations.tsv "
            "differential_methylation.rds differential_methylation.tsv differential_top.tsv differential_summary.tsv "
            "differential_sensitivity.rds differential_sensitivity_summary.tsv differential_sensitivity_top.tsv "
            "differential_qq.pdf differential_volcano.pdf dmr_summary.tsv dmr_ranges.tsv dmr_ranges.rds "
            "probe_mask_summary.tsv probe_mask_comparison.tsv"
        ),
    },
    "enrichment": {
        "previous": "differential", "command": ("Rscript", "scripts/06_enrichment.R"),
        "outputs": results("enrichment_results.rds enrichment_results.tsv enrichment_summary.tsv enrichment_top_pathways.pdf"),
    },
    "leave-one-batch-out": {
        "previous": "differential",
        "command": ("Rscript", "scripts/04_leave_one_batch_out.R"),
        "outputs": results("leave_one_batch_out_summary.tsv leave_one_batch_out_effects.tsv leave_one_batch_out.pdf"),
    },
    "report": {
        "previous": "enrichment",
        "requires": ("leave-one-batch-out",),
        "command": ("Rscript", "-e", 'rmarkdown::render("report/analysis_report.Rmd", output_dir = "report")'),
        "sources": ("report/analysis_report.Rmd",),
        "outputs": ("report/analysis_report.html",),
    },
    "replication": {
        "previous": "differential", "command": ("Rscript", "scripts/05_replication.R"),
        "sources": ("scripts/prepare_replication.py",),
        "inputs": (f"{REPLICATION_DATA}/GSE59685_series_matrix.txt.gz", f"{REPLICATION_DATA}/GSE59685_betas.csv.gz"),
        "outputs": results("replication_GSE59685_summary.tsv replication_GSE59685_effects.tsv replication_GSE59685.pdf"),
    },
}


def fingerprints(patterns):
    state = {}
    for pattern in patterns:
        files = sorted(Path().glob(pattern))
        state[pattern] = {str(path): [path.stat().st_size, path.stat().st_mtime_ns] for path in files if path.is_file()}
    return state


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature(stage, stamps):
    sources = stage.get("sources", ()) + R_ENV
    if stage["command"][1].endswith(".R"):
        sources += (stage["command"][1],)
    state = {
        "version": 1,
        "command": list(stage["command"]),
        "sources": {path: digest(path) for path in sources},
        "inputs": fingerprints(stage.get("inputs", ())),
        "upstream": stamps.get(stage.get("previous")),
    }
    if stage.get("requires"):
        state["requires"] = {name: stamps[name] for name in stage["requires"]}
    return state


def stamp_key(stamp):
    return hashlib.sha256(json.dumps(stamp, sort_keys=True).encode()).hexdigest()


def run_stage(name, command):
    started = time.perf_counter()
    with subprocess.Popen(command) as process:
        _, status, resource_usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    elapsed = time.perf_counter() - started
    rss_mib = resource_usage.ru_maxrss / (1024 ** (2 if sys.platform == "darwin" else 1))
    print(f"{name}: {elapsed:.2f} s, peak RSS {rss_mib:.1f} MiB", flush=True)
    with (CACHE / "benchmarks.tsv").open("a", newline="") as stream:
        if stream.tell() == 0:
            stream.write("stage\telapsed_seconds\tpeak_rss_mib\texit_code\n")
        stream.write(f"{name}\t{elapsed:.2f}\t{rss_mib:.1f}\t{process.returncode}\n")
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command)


def dependencies(stage):
    previous = (stage["previous"],) if "previous" in stage else ()
    return previous + stage.get("requires", ())


def stage_order(target):
    order = []

    def visit(name):
        if name not in order:
            for dependency in dependencies(STAGES[name]):
                visit(dependency)
            order.append(name)

    visit(target)
    return order


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=STAGES, default="report")
    parser.add_argument("--force-from", choices=STAGES, help="Rerun this stage and downstream stages up to the target")
    args = parser.parse_args()
    order = stage_order(args.target)
    if args.force_from and args.force_from not in order:
        parser.error("--force-from must be on the requested target's dependency path")
    os.chdir(ROOT)
    CACHE.mkdir(exist_ok=True)
    with (CACHE / "lock").open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.exit(1, "Another workflow run is active in this checkout.\n")
        stamps = {}
        forced_stages = set()
        for name in order:
            stage = STAGES[name]
            path = CACHE / f"{name}.json"
            forced = args.force_from == name or any(dep in forced_stages for dep in dependencies(stage))
            if forced:
                forced_stages.add(name)
            current = {"signature": signature(stage, stamps), "outputs": fingerprints(stage["outputs"])}
            try:
                cached = json.loads(path.read_text()) if path.exists() else None
            except json.JSONDecodeError:
                cached = None
            if not forced and cached == current and all(current["outputs"].values()):
                print(f"{name}: cached", flush=True)
            else:
                path.unlink(missing_ok=True)
                print(f"{name}: running", flush=True)
                run_stage(name, stage["command"])
                current = {"signature": signature(stage, stamps), "outputs": fingerprints(stage["outputs"])}
                missing = [pattern for pattern, files in current["outputs"].items() if not files]
                if missing:
                    parser.exit(1, f"{name} did not create expected outputs: {', '.join(missing)}\n")
                temporary = path.with_suffix(".tmp")
                temporary.write_text(json.dumps(current, sort_keys=True) + "\n")
                temporary.replace(path)
            stamps[name] = stamp_key(current)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode if error.returncode > 0 else 128 - error.returncode)

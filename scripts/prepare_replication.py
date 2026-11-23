#!/usr/bin/env python3
"""Extract selected frontal-cortex samples and probes from the public GSE59685 files."""

import argparse
import csv
import gzip
from pathlib import Path
import subprocess


GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59685"
DATA = Path("data/raw/GSE59685")


def download(path, url):
    if path.is_file():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    subprocess.run(
        ["curl", "--fail", "--location", "--retry", "3", "--continue-at", "-", "--output", str(partial), url],
        check=True,
    )
    partial.replace(path)


def frontal_samples(path):
    with gzip.open(path, "rt", newline="") as stream:
        rows = csv.reader((line for line in stream if line.startswith(("!Sample_geo_accession", "!Sample_characteristics_ch1"))),
                          delimiter="\t")
        accession = next(rows)
        if accession[0] != "!Sample_geo_accession":
            raise ValueError("Expected a GEO sample accession row before sample characteristics")
        fields = list(rows)
    if any(len(row) != len(accession) for row in fields):
        raise ValueError("GEO sample metadata has inconsistent row widths")
    samples = []
    for column, sample_id in enumerate(accession[1:], start=1):
        characteristics = dict(row[column].split(": ", 1) for row in fields if ": " in row[column])
        label = characteristics["ad.disease.status"]
        if characteristics["source tissue"] != "frontal cortex" or label not in {"AD", "C"}:
            continue
        barcode = characteristics["barcode"]
        samples.append({
            "sample_id": sample_id, "donor_id": characteristics["subjectid"], "barcode": barcode,
            "diagnosis": "Control" if label == "C" else label, "age": float(characteristics["age.brain"]),
            "sex": characteristics["Sex"].lower(), "chip": barcode.split("_", 1)[0],
        })
    if not samples or any(len({sample[field] for sample in samples}) != len(samples)
                          for field in ("sample_id", "donor_id", "barcode")):
        raise ValueError("Expected one frontal-cortex array per eligible donor")
    return samples


def extract_betas(path, candidates, samples, output):
    found = set()
    with gzip.open(path, "rb") as stream, output.open("w", newline="") as out:
        barcode_line = next(line for line in stream if line.strip() and not line.startswith(b"#"))
        barcodes = next(csv.reader([barcode_line.decode("ascii")]))
        sample_ids = next(csv.reader([next(stream).decode("ascii")]))
        original_ids = next(csv.reader([next(stream).decode("ascii")]))
        if len(set(barcodes)) != len(barcodes) or len(barcodes) != len(sample_ids) or len(barcodes) != len(original_ids):
            raise ValueError("Unexpected barcode or accession header in the GSE59685 beta matrix")
        positions = {barcode: i for i, barcode in enumerate(barcodes)}
        columns = [positions[sample["barcode"]] for sample in samples]
        if any(sample_ids[i] != sample["sample_id"] for i, sample in zip(columns, samples)):
            raise ValueError("GEO sample IDs and barcodes do not agree with the beta matrix")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["CpG"] + [sample["sample_id"] for sample in samples])
        for line in stream:
            probe = line.split(b",", 1)[0].strip(b'"\r\n')
            if probe not in candidates:
                continue
            if probe in found:
                raise ValueError(f"Duplicate beta row for {probe.decode()}")
            values = next(csv.reader([line.decode("ascii")]))
            if len(values) != len(barcodes):
                raise ValueError(f"Wrong column count on beta row {probe.decode()}")
            writer.writerow([probe.decode("ascii")] + [values[i] for i in columns])
            found.add(probe)
    return len(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("discovery_top", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    matrix = DATA / "GSE59685_series_matrix.txt.gz"
    beta = DATA / "GSE59685_betas.csv.gz"
    download(matrix, f"{GEO}/matrix/{matrix.name}")
    download(beta, f"{GEO}/suppl/{beta.name}")
    with args.discovery_top.open(newline="") as stream:
        candidates = {row["CpG"].encode("ascii") for row in csv.DictReader(stream, delimiter="\t")}
    if not candidates:
        parser.error("The discovery table contains no CpGs to compare")
    samples = frontal_samples(matrix)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "samples.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=samples[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)
    count = extract_betas(beta, candidates, samples, args.output_dir / "beta.tsv")
    print(f"Found {count} of {len(candidates)} candidate probes across {len(samples)} eligible frontal-cortex arrays")


if __name__ == "__main__":
    main()

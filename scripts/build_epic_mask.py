#!/usr/bin/env python3
"""Extract the published EPIC v1 hg19 general-mask IDs from the archived manifest."""

import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path

MANIFEST_SHA256 = "81fbb7def0759acb60385eff870f838cbbbb636ead4d4dcc2391d3ad6fb8b72e"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    digest = hashlib.sha256()
    with args.manifest.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != MANIFEST_SHA256:
        parser.error("Manifest checksum differs from the archived 20180909 EPIC hg19 release")

    with gzip.open(args.manifest, "rt", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        masked = sorted({row["probeID"] for row in rows if row["MASK_general"] == "TRUE"})
    if not masked:
        parser.error("The manifest did not return any MASK_general probes")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as stream:
                stream.write("CpG\n" + "\n".join(masked) + "\n")
    print(f"Wrote {len(masked):,} excluded probe IDs to {args.output}")


if __name__ == "__main__":
    main()

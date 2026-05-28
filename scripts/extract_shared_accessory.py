#!/usr/bin/env python3

import csv
import sys
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract shared accessory orthogroups from OrthoFinder "
            "Orthogroups.GeneCount.tsv.\n"
            "Shared accessory is defined here as non-core and present in >=2 strains."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to Orthogroups.GeneCount.tsv"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)

    if not input_path.exists():
        sys.stderr.write(f"ERROR: Input file not found: {input_path}\n")
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    gene_count_out = outdir / "shared_accessory_gene_count.tsv"
    presence_absence_out = outdir / "shared_accessory_presence_absence.tsv"
    summary_out = outdir / "shared_accessory_summary.tsv"

    n_all = 0
    n_core = 0
    n_singleton = 0
    n_shared_accessory = 0

    with input_path.open("r", newline="") as infile, \
         gene_count_out.open("w", newline="") as gc_outfile, \
         presence_absence_out.open("w", newline="") as pa_outfile, \
         summary_out.open("w", newline="") as sum_outfile:

        reader = csv.reader(infile, delimiter="\t")
        gc_writer = csv.writer(gc_outfile, delimiter="\t")
        pa_writer = csv.writer(pa_outfile, delimiter="\t")
        sum_writer = csv.writer(sum_outfile, delimiter="\t")

        header = next(reader)
        if len(header) < 4:
            sys.stderr.write("ERROR: Unexpected header format in input file.\n")
            sys.exit(1)

        # 1列目: Orthogroup, 最後: Total
        if header[-1] != "Total":
            sys.stderr.write(
                "ERROR: Expected last column to be 'Total'. "
                f"Found '{header[-1]}' instead.\n"
            )
            sys.exit(1)

        strain_names = header[1:-1]
        n_strains_total = len(strain_names)

        gc_writer.writerow(
            ["Orthogroup", "n_strains_present", "total_gene_count", "category"] + strain_names
        )
        pa_writer.writerow(
            ["Orthogroup", "n_strains_present", "category"] + strain_names
        )
        sum_writer.writerow(
            ["Orthogroup", "n_strains_present", "total_gene_count", "category"]
        )

        for row in reader:
            if not row:
                continue

            orthogroup = row[0]

            try:
                counts = [int(x) for x in row[1:-1]]   # Total列を除く
                total_from_file = int(row[-1])
            except ValueError:
                sys.stderr.write(f"ERROR: Non-integer value found in row: {orthogroup}\n")
                sys.exit(1)

            if len(counts) != n_strains_total:
                sys.stderr.write(
                    f"ERROR: Column number mismatch in row: {orthogroup}\n"
                )
                sys.exit(1)

            n_all += 1

            presence = [1 if x > 0 else 0 for x in counts]
            n_present = sum(presence)
            total_gene_count = sum(counts)

            if total_gene_count != total_from_file:
                sys.stderr.write(
                    f"WARNING: Total mismatch in {orthogroup}: "
                    f"sum(counts)={total_gene_count}, file_total={total_from_file}\n"
                )

            if n_present == n_strains_total:
                category = "core"
                n_core += 1
            elif n_present == 1:
                category = "singleton"
                n_singleton += 1
            elif 2 <= n_present < n_strains_total:
                category = "shared_accessory"
                n_shared_accessory += 1
            else:
                category = "absent_or_unexpected"

            if category == "shared_accessory":
                gc_writer.writerow(
                    [orthogroup, n_present, total_gene_count, category] + counts
                )
                pa_writer.writerow(
                    [orthogroup, n_present, category] + presence
                )
                sum_writer.writerow(
                    [orthogroup, n_present, total_gene_count, category]
                )

    sys.stderr.write("Finished.\n")
    sys.stderr.write(f"Input file: {input_path}\n")
    sys.stderr.write(f"Output directory: {outdir}\n")
    sys.stderr.write(f"Total orthogroups: {n_all}\n")
    sys.stderr.write(f"Core orthogroups: {n_core}\n")
    sys.stderr.write(f"Singleton orthogroups: {n_singleton}\n")
    sys.stderr.write(f"Shared accessory orthogroups: {n_shared_accessory}\n")
    sys.stderr.write(f"Wrote: {gene_count_out}\n")
    sys.stderr.write(f"Wrote: {presence_absence_out}\n")
    sys.stderr.write(f"Wrote: {summary_out}\n")


if __name__ == "__main__":
    main()


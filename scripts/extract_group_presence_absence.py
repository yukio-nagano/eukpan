#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract orthogroups based on focal-group and non-group "
            "presence thresholds from OrthoFinder Orthogroups.GeneCount.tsv, "
            "and output a presence/absence table."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to Orthogroups.GeneCount.tsv"
    )
    parser.add_argument(
        "-g", "--group-file",
        required=True,
        help="Text file listing focal group sample names, one per line"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )
    parser.add_argument(
        "--group-min-pct",
        type=float,
        default=0.0,
        help="Minimum allowed presence percentage in focal group (default: 0)"
    )
    parser.add_argument(
        "--group-max-pct",
        type=float,
        default=100.0,
        help="Maximum allowed presence percentage in focal group (default: 100)"
    )
    parser.add_argument(
        "--nongroup-min-pct",
        type=float,
        default=0.0,
        help="Minimum allowed presence percentage in non-group (default: 0)"
    )
    parser.add_argument(
        "--nongroup-max-pct",
        type=float,
        default=100.0,
        help="Maximum allowed presence percentage in non-group (default: 100)"
    )
    return parser.parse_args()


def read_sample_list(path: Path):
    samples = []
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if s:
                samples.append(s)
    return samples


def main():
    args = parse_args()

    input_path = Path(args.input)
    group_path = Path(args.group_file)
    output_path = Path(args.output)

    if not input_path.exists():
        sys.stderr.write(f"ERROR: Input file not found: {input_path}\n")
        sys.exit(1)

    if not group_path.exists():
        sys.stderr.write(f"ERROR: Group file not found: {group_path}\n")
        sys.exit(1)

    group_samples = read_sample_list(group_path)
    group_set = set(group_samples)

    if len(group_samples) == 0:
        sys.stderr.write("ERROR: Group file is empty.\n")
        sys.exit(1)

    with input_path.open("r", newline="") as infile, output_path.open("w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        header = next(reader)
        if len(header) < 4:
            sys.stderr.write("ERROR: Unexpected header format.\n")
            sys.exit(1)

        if header[-1] != "Total":
            sys.stderr.write(
                f"ERROR: Expected last column to be 'Total', but found '{header[-1]}'.\n"
            )
            sys.exit(1)

        all_samples = header[1:-1]

        missing = [s for s in group_samples if s not in all_samples]
        if missing:
            sys.stderr.write("ERROR: The following group samples are not in the input table:\n")
            for s in missing:
                sys.stderr.write(f"  {s}\n")
            sys.exit(1)

        nongroup_samples = [s for s in all_samples if s not in group_set]
        if len(nongroup_samples) == 0:
            sys.stderr.write("ERROR: No non-group samples remain.\n")
            sys.exit(1)

        sample_to_idx = {sample: i for i, sample in enumerate(all_samples)}
        group_idx = [sample_to_idx[s] for s in group_samples]
        nongroup_idx = [sample_to_idx[s] for s in nongroup_samples]

        writer.writerow(
            [
                "Orthogroup",
                "group_present_n",
                "group_total_n",
                "group_present_pct",
                "nongroup_present_n",
                "nongroup_total_n",
                "nongroup_present_pct",
            ] + all_samples
        )

        n_total = 0
        n_pass = 0

        for row in reader:
            if not row:
                continue

            n_total += 1
            orthogroup = row[0]

            try:
                counts = [int(x) for x in row[1:-1]]
            except ValueError:
                sys.stderr.write(f"ERROR: Non-integer count in row: {orthogroup}\n")
                sys.exit(1)

            if len(counts) != len(all_samples):
                sys.stderr.write(f"ERROR: Column number mismatch in row: {orthogroup}\n")
                sys.exit(1)

            presence = [1 if x > 0 else 0 for x in counts]

            group_present_n = sum(presence[i] for i in group_idx)
            group_total_n = len(group_idx)
            group_present_pct = 100.0 * group_present_n / group_total_n

            nongroup_present_n = sum(presence[i] for i in nongroup_idx)
            nongroup_total_n = len(nongroup_idx)
            nongroup_present_pct = 100.0 * nongroup_present_n / nongroup_total_n

            if (
                args.group_min_pct <= group_present_pct <= args.group_max_pct and
                args.nongroup_min_pct <= nongroup_present_pct <= args.nongroup_max_pct
            ):
                n_pass += 1
                writer.writerow(
                    [
                        orthogroup,
                        group_present_n,
                        group_total_n,
                        f"{group_present_pct:.2f}",
                        nongroup_present_n,
                        nongroup_total_n,
                        f"{nongroup_present_pct:.2f}",
                    ] + presence
                )

    sys.stderr.write("Finished.\n")
    sys.stderr.write(f"Input: {input_path}\n")
    sys.stderr.write(f"Group file: {group_path}\n")
    sys.stderr.write(f"Output: {output_path}\n")
    sys.stderr.write(f"group_min_pct: {args.group_min_pct}\n")
    sys.stderr.write(f"group_max_pct: {args.group_max_pct}\n")
    sys.stderr.write(f"nongroup_min_pct: {args.nongroup_min_pct}\n")
    sys.stderr.write(f"nongroup_max_pct: {args.nongroup_max_pct}\n")
    sys.stderr.write(f"Total orthogroups: {n_total}\n")
    sys.stderr.write(f"Matched orthogroups: {n_pass}\n")


if __name__ == "__main__":
    main()

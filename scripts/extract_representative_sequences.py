#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_fasta(path):
    records = []
    header = None
    seq_lines = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines)))

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Extract one representative sequence (longest) for each Orthogroup."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file (e.g. A_rare_nonA_common_presence_absence.tsv)"
    )
    parser.add_argument(
        "-s", "--seqdir",
        required=True,
        help="Orthogroup_Sequences directory"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    seqdir = Path(args.seqdir)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not seqdir.exists():
        raise FileNotFoundError(f"Sequence directory not found: {seqdir}")

    orthogroups = []
    with input_path.open() as fh:
        header = next(fh)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            og = line.split("\t")[0]
            orthogroups.append(og)

    n_found = 0
    n_missing = 0

    with output_path.open("w") as out:
        for og in orthogroups:
            fasta_path = seqdir / f"{og}.fa"
            if not fasta_path.exists():
                print(f"WARNING: missing sequence file for {og}")
                n_missing += 1
                continue

            records = parse_fasta(fasta_path)
            if not records:
                print(f"WARNING: empty sequence file for {og}")
                n_missing += 1
                continue

            rep_header, rep_seq = max(records, key=lambda x: len(x[1]))

            out.write(f">{og}|{rep_header}\n")
            for i in range(0, len(rep_seq), 80):
                out.write(rep_seq[i:i+80] + "\n")

            n_found += 1

    print(f"Orthogroups in input: {len(orthogroups)}")
    print(f"Representative sequences written: {n_found}")
    print(f"Missing/empty orthogroups: {n_missing}")
    print(f"Output FASTA: {output_path}")


if __name__ == "__main__":
    main()

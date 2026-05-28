#!/usr/bin/env bash

set -euo pipefail

PIPELINE_VERSION="0.1.0"

usage() {
    cat << EOF
Usage:
  bash eukpan.sh \\
    --gff-dir <GFF_DIR> \\
    --genome-dir <GENOME_DIR> \\
    --outdir <OUTDIR> \\
    --threads <THREADS> \\
    [--orthofinder-env <CONDA_ENV_NAME>] \\
    [--scripts-dir <SCRIPTS_DIR>]

Example:
  bash eukpan.sh \\
    --gff-dir gff \\
    --genome-dir assembled_genomes \\
    --outdir pangenome_results \\
    --threads 16 \\
    --orthofinder-env orthofinder \\
    --scripts-dir scripts
EOF
}

GFF_DIR=""
GENOME_DIR=""
OUTDIR=""
THREADS=""
ORTHOFINDER_ENV="orthofinder"
SCRIPTS_DIR="scripts"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gff-dir)
            GFF_DIR="$2"
            shift 2
            ;;
        --genome-dir)
            GENOME_DIR="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --orthofinder-env)
            ORTHOFINDER_ENV="$2"
            shift 2
            ;;
        --scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$GFF_DIR" || -z "$GENOME_DIR" || -z "$OUTDIR" || -z "$THREADS" ]]; then
    echo "ERROR: Missing required arguments." >&2
    usage
    exit 1
fi

if [[ ! -d "$GFF_DIR" ]]; then
    echo "ERROR: GFF/GTF directory not found: $GFF_DIR" >&2
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "ERROR: Genome directory not found: $GENOME_DIR" >&2
    exit 1
fi

if [[ ! -d "$SCRIPTS_DIR" ]]; then
    echo "ERROR: Scripts directory not found: $SCRIPTS_DIR" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
fi

if [[ "$THREADS" -lt 1 ]]; then
    echo "ERROR: --threads must be greater than or equal to 1." >&2
    exit 1
fi

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/logs"

PIPELINE_LOG="$OUTDIR/logs/pipeline.log"

exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "=============================================="
echo "EukPan: Eukaryotic microbial pangenome pipeline"
echo "Version: $PIPELINE_VERSION"
echo "=============================================="
echo "GFF/GTF directory:  $GFF_DIR"
echo "Genome directory:   $GENOME_DIR"
echo "Output directory:   $OUTDIR"
echo "Threads:            $THREADS"
echo "OrthoFinder env:    $ORTHOFINDER_ENV"
echo "Scripts directory:  $SCRIPTS_DIR"
echo "Pipeline log:       $PIPELINE_LOG"
echo "=============================================="

echo
echo "[Initial cleanup] Removing old pipeline-generated files in the output directory..."

rm -rf "$OUTDIR/fixed_gff"
rm -rf "$OUTDIR/longest_gff"
rm -rf "$OUTDIR/proteins"
rm -rf "$OUTDIR/proteomes"
rm -rf "$OUTDIR/idmap"
rm -rf "$OUTDIR/phylo"
rm -rf "$OUTDIR/accessory_results"
rm -rf "$OUTDIR/accessory_analysis"

mkdir -p "$OUTDIR/fixed_gff"
mkdir -p "$OUTDIR/longest_gff"
mkdir -p "$OUTDIR/proteins"
mkdir -p "$OUTDIR/proteomes"
mkdir -p "$OUTDIR/idmap"
mkdir -p "$OUTDIR/phylo"
mkdir -p "$OUTDIR/phylo/OG_seqs"
mkdir -p "$OUTDIR/phylo/ids"
mkdir -p "$OUTDIR/phylo/aln"
mkdir -p "$OUTDIR/phylo/trimmed_auto1"
mkdir -p "$OUTDIR/phylo/trimmed_auto1_taxa"
mkdir -p "$OUTDIR/phylo/trimmed_auto1_short_removed"
mkdir -p "$OUTDIR/phylo/trimmed_auto1_sequence_removed"
mkdir -p "$OUTDIR/accessory_results"
mkdir -p "$OUTDIR/accessory_analysis"

echo "Old pipeline-generated files were removed."

echo
echo "[Step 0] Checking required commands in the current environment..."

REQUIRED_COMMANDS=(
    agat
    agat_convert_sp_gxf2gxf.pl
    agat_sp_keep_longest_isoform.pl
    gffread
    seqkit
    mafft
    parallel
    trimal
    awk
    sed
    grep
    sort
    uniq
    diff
    python3
    conda
)

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found in current environment: $cmd" >&2
        exit 1
    fi
done

echo "Required preprocessing and alignment commands were found."

echo
echo "[Step 0b] Checking optional command: Inkscape"

if ! command -v inkscape >/dev/null 2>&1; then
    echo "WARNING: Inkscape was not found. EMF export in analyze_shared_accessory.py may fail." >&2
    echo "WARNING: PNG and SVG outputs should still be generated." >&2
else
    echo "Inkscape was found."
fi

echo
echo "[Step 0c] Checking Python helper scripts..."

REQUIRED_SCRIPTS=(
    extract_shared_accessory.py
    analyze_shared_accessory.py
    extract_group_presence_absence.py
    extract_representative_sequences.py
)

for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [[ ! -f "$SCRIPTS_DIR/$script" ]]; then
        echo "ERROR: Required Python script not found: $SCRIPTS_DIR/$script" >&2
        exit 1
    fi
done

echo "Required Python helper scripts were found."

echo
echo "[Step 0d] Checking OrthoFinder in conda environment: ${ORTHOFINDER_ENV}"

if ! conda run -n "$ORTHOFINDER_ENV" orthofinder --version >/dev/null 2>&1; then
    echo "ERROR: OrthoFinder could not be executed in conda environment: ${ORTHOFINDER_ENV}" >&2
    echo "Please check with:" >&2
    echo "  conda activate ${ORTHOFINDER_ENV}" >&2
    echo "  orthofinder --version" >&2
    exit 1
fi

echo "OrthoFinder was found in conda environment: ${ORTHOFINDER_ENV}"

echo
echo "[Step 0e] Software versions used directly by this pipeline"

echo
echo "AGAT:"
agat --version 2>&1 || true

echo
echo "gffread:"
gffread --version 2>&1 | head -n 3 || true

echo
echo "seqkit:"
seqkit version || true

echo
echo "MAFFT:"
mafft --version 2>&1 | head -n 2 || true

echo
echo "GNU parallel:"
parallel --version | head -n 2 || true

echo
echo "trimAl:"
trimal --version 2>&1 | head -n 3 || true

echo
echo "Python:"
python3 --version || true

if command -v inkscape >/dev/null 2>&1; then
    echo
    echo "Inkscape:"
    inkscape --version || true
fi

echo
echo "OrthoFinder:"
conda run -n "$ORTHOFINDER_ENV" orthofinder --version 2>&1 || true

echo
echo "[Step 1] Checking input files..."

shopt -s nullglob

GFF_FILES=( "$GFF_DIR"/*.gff "$GFF_DIR"/*.gff3 "$GFF_DIR"/*.gtf )
GENOME_FILES=( "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta )

if [[ ${#GFF_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .gff, .gff3, or .gtf files found in: $GFF_DIR" >&2
    exit 1
fi

if [[ ${#GENOME_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .fna, .fa, or .fasta files found in: $GENOME_DIR" >&2
    exit 1
fi

echo "Found ${#GFF_FILES[@]} GFF/GFF3/GTF files."
echo "Found ${#GENOME_FILES[@]} genome FASTA files."

echo
echo "[Step 2] Converting GFF/GFF3/GTF files with AGAT..."

for f in "${GFF_FILES[@]}"; do
    base=$(basename "$f")
    name="${base%.*}"

    echo "  Converting annotation file: $name"

    agat_convert_sp_gxf2gxf.pl \
        -g "$f" \
        -o "$OUTDIR/fixed_gff/${name}.fixed.gff3" \
        > "$OUTDIR/logs/${name}.agat_convert.log" \
        2> "$OUTDIR/logs/${name}.agat_convert.err"
done

echo
echo "[Step 3] Keeping longest isoforms with AGAT..."

for f in "$OUTDIR/fixed_gff"/*.fixed.gff3; do
    [ -e "$f" ] || continue

    name=$(basename "$f" .fixed.gff3)

    echo "  Keeping longest isoforms: $name"

    agat_sp_keep_longest_isoform.pl \
        --gff "$f" \
        -o "$OUTDIR/longest_gff/${name}.longest.gff3" \
        > "$OUTDIR/logs/${name}.agat_longest.log" \
        2> "$OUTDIR/logs/${name}.agat_longest.err"
done

echo
echo "[Step 4] Cleaning temporary AGAT files if present..."

rm -rf agat_*

echo
echo "[Step 5] Extracting protein sequences with gffread..."

for f in "$OUTDIR/longest_gff"/*.longest.gff3; do
    [ -e "$f" ] || continue

    name=$(basename "$f" .longest.gff3)

    genome=""

    for ext in fna fa fasta; do
        candidate="$GENOME_DIR/${name}.${ext}"
        if [[ -f "$candidate" ]]; then
            genome="$candidate"
            break
        fi
    done

    if [[ -z "$genome" ]]; then
        echo "WARNING: Genome FASTA not found for sample: $name. Skipping." >&2
        echo "WARNING: Expected one of: ${GENOME_DIR}/${name}.fna, .fa, or .fasta" >&2
        continue
    fi

    echo "  Extracting proteins: $name"

    gffread "$f" \
        -g "$genome" \
        -y "$OUTDIR/proteins/${name}.protein.fa" \
        > "$OUTDIR/logs/${name}.gffread.log" \
        2> "$OUTDIR/logs/${name}.gffread.err"
done

echo
echo "[Step 6] Simplifying protein FASTA headers..."

PROTEIN_FILES=( "$OUTDIR/proteins"/*.protein.fa )

if [[ ${#PROTEIN_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No protein FASTA files were generated in: $OUTDIR/proteins" >&2
    exit 1
fi

for f in "${PROTEIN_FILES[@]}"; do
    [ -e "$f" ] || continue

    s=$(basename "$f" .protein.fa)

    echo "  Simplifying headers: $s"

    rm -f "$OUTDIR/idmap/${s}.idmap.tsv"

    awk -v S="$s" -v OUTDIR="$OUTDIR" '
    BEGIN{n=0}
    /^>/{
        n++
        old=substr($0,2)
        new=sprintf("%s_%06d", S, n)
        print new "\t" old >> (OUTDIR "/idmap/" S ".idmap.tsv")
        print ">" new
        next
    }
    {print}
    ' "$f" > "$OUTDIR/proteomes/${s}.fa"
done

echo
echo "[Step 7] Checking proteome FASTA files before OrthoFinder..."

PROTEOME_FILES=( "$OUTDIR/proteomes"/*.fa )

if [[ ${#PROTEOME_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No proteome FASTA files found in: $OUTDIR/proteomes" >&2
    exit 1
fi

SAMPLE_COUNT=${#PROTEOME_FILES[@]}

echo "Found ${SAMPLE_COUNT} proteome FASTA files."

echo
echo "[Step 8] Running OrthoFinder in conda environment: ${ORTHOFINDER_ENV}"

conda run -n "$ORTHOFINDER_ENV" orthofinder \
    -f "$OUTDIR/proteomes" \
    -t "$THREADS" \
    -a "$THREADS" \
    > "$OUTDIR/logs/orthofinder.log" \
    2> "$OUTDIR/logs/orthofinder.err"

echo
echo "[Step 9] Detecting the latest OrthoFinder results directory..."

ORTHOFINDER_BASE="$OUTDIR/proteomes/OrthoFinder"

if [[ ! -d "$ORTHOFINDER_BASE" ]]; then
    echo "ERROR: OrthoFinder base directory not found: $ORTHOFINDER_BASE" >&2
    exit 1
fi

ORTHOFINDER_RESULTS=$(find "$ORTHOFINDER_BASE" -maxdepth 1 -type d -name "Results_*" | sort | tail -n 1)

if [[ -z "$ORTHOFINDER_RESULTS" || ! -d "$ORTHOFINDER_RESULTS" ]]; then
    echo "ERROR: Could not detect OrthoFinder Results_* directory." >&2
    exit 1
fi

echo "Detected OrthoFinder results directory:"
echo "  $ORTHOFINDER_RESULTS"

ORTHOGROUP_DIR="$ORTHOFINDER_RESULTS/Orthogroups"
ORTHO_TSV="$ORTHOGROUP_DIR/Orthogroups.tsv"
SCO_LIST="$ORTHOGROUP_DIR/Orthogroups_SingleCopyOrthologues.txt"
GENECOUNT_TSV="$ORTHOGROUP_DIR/Orthogroups.GeneCount.tsv"
OG_SEQ_DIR="$ORTHOFINDER_RESULTS/Orthogroup_Sequences"

if [[ ! -f "$ORTHO_TSV" ]]; then
    echo "ERROR: Orthogroups.tsv not found: $ORTHO_TSV" >&2
    exit 1
fi

if [[ ! -f "$SCO_LIST" ]]; then
    echo "ERROR: Orthogroups_SingleCopyOrthologues.txt not found: $SCO_LIST" >&2
    exit 1
fi

if [[ ! -f "$GENECOUNT_TSV" ]]; then
    echo "ERROR: Orthogroups.GeneCount.tsv not found: $GENECOUNT_TSV" >&2
    exit 1
fi

if [[ ! -d "$OG_SEQ_DIR" ]]; then
    echo "ERROR: Orthogroup_Sequences directory not found: $OG_SEQ_DIR" >&2
    exit 1
fi

SCO_COUNT=$(wc -l < "$SCO_LIST")
echo "Number of single-copy orthologues: $SCO_COUNT"

if [[ "$SCO_COUNT" -eq 0 ]]; then
    echo "ERROR: No single-copy orthologues were detected." >&2
    exit 1
fi

echo
echo "[Step 10] Preparing single-copy core orthologue sequences..."

echo "  Concatenating all proteomes into one FASTA file."

cat "$OUTDIR"/proteomes/*.fa > "$OUTDIR/phylo/all_proteins.fa"

echo "  Extracting single-copy core orthologue sequences."
echo "  seqkit messages will be written to:"
echo "    $OUTDIR/logs/seqkit_extract_orthogroups.log"

SEQKIT_LOG="$OUTDIR/logs/seqkit_extract_orthogroups.log"
: > "$SEQKIT_LOG"

while read -r og; do
    [[ -z "$og" ]] && continue

    grep -P "^${og}\t" "$ORTHO_TSV" \
    | tr '\t' '\n' \
    | tail -n +2 \
    | sed '/^$/d' \
    | sed 's/, /\n/g' \
    | sed 's/ /\n/g' \
    | sed '/^$/d' > "$OUTDIR/phylo/ids/${og}.ids"

    seqkit grep \
        -f "$OUTDIR/phylo/ids/${og}.ids" \
        "$OUTDIR/phylo/all_proteins.fa" \
        > "$OUTDIR/phylo/OG_seqs/${og}.fa" \
        2>> "$SEQKIT_LOG"
done < "$SCO_LIST"

echo
echo "[Step 11] Checking sequence counts in single-copy orthologue FASTA files..."

BAD_OG_COUNT=0

for f in "$OUTDIR"/phylo/OG_seqs/*.fa; do
    [ -e "$f" ] || continue

    n=$(grep -c '^>' "$f")

    if [[ "$n" -ne "$SAMPLE_COUNT" ]]; then
        echo "ERROR: Unexpected sequence count in $(basename "$f"): $n, expected $SAMPLE_COUNT" >&2
        BAD_OG_COUNT=$((BAD_OG_COUNT + 1))
    fi
done

if [[ "$BAD_OG_COUNT" -ne 0 ]]; then
    echo "ERROR: Some single-copy orthologue FASTA files do not contain exactly one sequence per sample." >&2
    echo "ERROR: Please check ID extraction, FASTA headers, and OrthoFinder outputs." >&2
    exit 1
fi

echo "All single-copy orthologue FASTA files contain exactly one sequence per sample."

echo
echo "[Step 12] Running MAFFT for all single-copy orthologues..."
echo "MAFFT stdout/stderr will be written to:"
echo "  $OUTDIR/logs/mafft.log"

: > "$OUTDIR/logs/mafft.log"

parallel -j "$THREADS" \
    'base=$(basename {} .fa); mafft --auto {} > '"$OUTDIR"'/phylo/aln/${base}.aln.fa 2>> '"$OUTDIR"'/logs/mafft.log' \
    ::: "$OUTDIR"/phylo/OG_seqs/*.fa

echo "MAFFT finished."

echo
echo "[Step 13] Trimming alignments with trimAl automated1..."

for f in "$OUTDIR"/phylo/aln/*.aln.fa; do
    [ -e "$f" ] || continue

    base=$(basename "$f" .aln.fa)

    trimal \
        -in "$f" \
        -out "$OUTDIR/phylo/trimmed_auto1/${base}.trim.fa" \
        -automated1 \
        > "$OUTDIR/logs/${base}.trimal.log" \
        2> "$OUTDIR/logs/${base}.trimal.err"
done

echo
echo "[Step 14] Checking sequence counts after trimAl..."

SEQ_LOSS_LOG="$OUTDIR/phylo/alignments_removed_due_to_sequence_loss.tsv"
echo -e "alignment\tobserved_sequence_count\texpected_sequence_count\treason" > "$SEQ_LOSS_LOG"

SEQ_LOSS_COUNT=0

for f in "$OUTDIR"/phylo/trimmed_auto1/*.trim.fa; do
    [ -e "$f" ] || continue

    n=$(grep -c '^>' "$f")

    if [[ "$n" -ne "$SAMPLE_COUNT" ]]; then
        echo "WARNING: Sequence count changed after trimming in $(basename "$f"): $n, expected $SAMPLE_COUNT" >&2
        echo "WARNING: This alignment will be excluded from concatenation." >&2

        echo -e "$(basename "$f")\t$n\t$SAMPLE_COUNT\tsequence_removed_by_trimAl_or_all_gap_sequence" >> "$SEQ_LOSS_LOG"

        mv "$f" "$OUTDIR/phylo/trimmed_auto1_sequence_removed/"
        SEQ_LOSS_COUNT=$((SEQ_LOSS_COUNT + 1))
    fi
done

if [[ "$SEQ_LOSS_COUNT" -eq 0 ]]; then
    echo "No alignments lost sequences after trimAl."
    rm -f "$SEQ_LOSS_LOG"
else
    echo "Alignments with sequence loss after trimAl were found and excluded."
    echo "Number of excluded alignments: $SEQ_LOSS_COUNT"
    echo "List of excluded alignments:"
    echo "  $SEQ_LOSS_LOG"
fi

N_AFTER_SEQ_LOSS_FILTER=$(find "$OUTDIR/phylo/trimmed_auto1" -maxdepth 1 -type f -name "*.trim.fa" | wc -l)
echo "Number of trimmed alignments remaining after sequence-loss filtering: $N_AFTER_SEQ_LOSS_FILTER"

if [[ "$N_AFTER_SEQ_LOSS_FILTER" -eq 0 ]]; then
    echo "ERROR: No trimmed alignments remain after sequence-loss filtering." >&2
    exit 1
fi

echo
echo "[Step 15] Checking and moving extremely short alignments..."

SHORT_LOG="$OUTDIR/phylo/short_alignments_removed.tsv"
echo -e "alignment\tlength" > "$SHORT_LOG"

SHORT_COUNT=0

for f in "$OUTDIR"/phylo/trimmed_auto1/*.trim.fa; do
    [ -e "$f" ] || continue

    len=$(awk '
    /^>/{
        if(seqlen>max) max=seqlen
        seqlen=0
        next
    }
    {
        seqlen+=length($0)
    }
    END{
        if(seqlen>max) max=seqlen
        print max+0
    }' "$f")

    if [[ "$len" -lt 30 ]]; then
        echo -e "$(basename "$f")\t$len" >> "$SHORT_LOG"
        mv "$f" "$OUTDIR/phylo/trimmed_auto1_short_removed/"
        SHORT_COUNT=$((SHORT_COUNT + 1))
    fi
done

if [[ "$SHORT_COUNT" -eq 0 ]]; then
    echo "No extremely short alignments were found."
    rm -f "$SHORT_LOG"
else
    echo "Short alignments were found and moved to:"
    echo "  $OUTDIR/phylo/trimmed_auto1_short_removed"
    echo "List of removed alignments:"
    echo "  $SHORT_LOG"
fi

N_TRIMMED=$(find "$OUTDIR/phylo/trimmed_auto1" -maxdepth 1 -type f -name "*.trim.fa" | wc -l)
echo "Number of trimmed alignments used for concatenation: $N_TRIMMED"

if [[ "$N_TRIMMED" -eq 0 ]]; then
    echo "ERROR: No trimmed alignments remain for concatenation." >&2
    exit 1
fi

echo
echo "[Step 16] Rewriting FASTA headers to taxon names..."

rm -f "$OUTDIR"/phylo/trimmed_auto1_taxa/*.trim.fa

for f in "$OUTDIR"/phylo/trimmed_auto1/*.trim.fa; do
    [ -e "$f" ] || continue

    base=$(basename "$f")

    awk '
    /^>/{
        id=substr($0,2)
        sub(/_[0-9]+$/, "", id)
        print ">" id
        next
    }
    {print}
    ' "$f" > "$OUTDIR/phylo/trimmed_auto1_taxa/$base"
done

echo
echo "[Step 17] Checking duplicate taxa in alignments..."

DUP_COUNT=0

for f in "$OUTDIR"/phylo/trimmed_auto1_taxa/*.trim.fa; do
    [ -e "$f" ] || continue

    dup=$(grep '^>' "$f" | sed 's/^>//' | awk '{print $1}' | sort | uniq -d)

    if [[ -n "$dup" ]]; then
        echo "ERROR: Duplicate taxon found in $f" >&2
        echo "$dup" >&2
        DUP_COUNT=$((DUP_COUNT + 1))
    fi
done

if [[ "$DUP_COUNT" -ne 0 ]]; then
    echo "ERROR: Duplicate taxa were detected. This is not valid for single-copy core alignment." >&2
    exit 1
fi

echo "No duplicate taxa were detected."

echo
echo "[Step 18] Checking taxa sets across all alignments..."

EXPECTED_IDS="$OUTDIR/phylo/expected_taxa.ids"
CURRENT_IDS="$OUTDIR/phylo/current.ids"
MISSING_IDS="$OUTDIR/phylo/missing.ids"
EXTRA_IDS="$OUTDIR/phylo/extra.ids"
BAD_TAXA_LOG="$OUTDIR/phylo/bad_taxa_set.tsv"

find "$OUTDIR/proteomes" -maxdepth 1 -type f -name "*.fa" \
    -exec basename {} .fa \; \
    | sort > "$EXPECTED_IDS"

EXPECTED_N=$(wc -l < "$EXPECTED_IDS")

echo "Expected number of taxa: $EXPECTED_N"
echo "Expected taxa list:"
cat "$EXPECTED_IDS"

echo -e "alignment\tobserved_taxa_count\texpected_taxa_count\tmissing_taxa\textra_taxa" > "$BAD_TAXA_LOG"

BAD_TAXA_COUNT=0

for f in "$OUTDIR"/phylo/trimmed_auto1_taxa/*.trim.fa; do
    [ -e "$f" ] || continue

    grep '^>' "$f" | sed 's/^>//' | awk '{print $1}' | sort > "$CURRENT_IDS"

    OBSERVED_N=$(wc -l < "$CURRENT_IDS")

    if ! diff -q "$EXPECTED_IDS" "$CURRENT_IDS" >/dev/null; then
        comm -23 "$EXPECTED_IDS" "$CURRENT_IDS" | paste -sd "," - > "$MISSING_IDS"
        comm -13 "$EXPECTED_IDS" "$CURRENT_IDS" | paste -sd "," - > "$EXTRA_IDS"

        missing=$(cat "$MISSING_IDS")
        extra=$(cat "$EXTRA_IDS")

        echo -e "$(basename "$f")\t${OBSERVED_N}\t${EXPECTED_N}\t${missing}\t${extra}" >> "$BAD_TAXA_LOG"

        echo "ERROR: Taxa set differs in $(basename "$f")" >&2
        echo "ERROR: Observed taxa count: ${OBSERVED_N}; expected: ${EXPECTED_N}" >&2
        echo "ERROR: Missing taxa: ${missing}" >&2
        echo "ERROR: Extra taxa: ${extra}" >&2

        BAD_TAXA_COUNT=$((BAD_TAXA_COUNT + 1))
    fi
done

rm -f "$CURRENT_IDS" "$MISSING_IDS" "$EXTRA_IDS"

if [[ "$BAD_TAXA_COUNT" -ne 0 ]]; then
    echo "ERROR: Some alignments have inconsistent taxa sets." >&2
    echo "ERROR: See detailed log:" >&2
    echo "  $BAD_TAXA_LOG" >&2
    exit 1
fi

rm -f "$BAD_TAXA_LOG"

echo "All alignments have the expected taxa set."

echo
echo "[Step 19] Concatenating core alignments..."

python3 - "$OUTDIR/phylo/trimmed_auto1_taxa" "$OUTDIR/phylo/concat_all.fa" "$OUTDIR/phylo/concat_partitions.txt" << 'PYCODE'
import sys
from pathlib import Path
from collections import OrderedDict

indir = Path(sys.argv[1])
output = Path(sys.argv[2])
partitions_out = Path(sys.argv[3])

files = sorted(indir.glob("*.trim.fa"))
if not files:
    sys.exit("ERROR: No alignment files found for concatenation.")

seqs = OrderedDict()
partitions = []
start = 1

for f in files:
    local = OrderedDict()
    current_id = None

    with f.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                if current_id in local:
                    raise ValueError(f"Duplicate taxon in {f}: {current_id}")
                local[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(f"Sequence line found before header in {f}")
                local[current_id].append(line)

    local = OrderedDict((k, "".join(v)) for k, v in local.items())

    if not local:
        raise ValueError(f"Empty alignment: {f}")

    lengths = {len(v) for v in local.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent sequence lengths in {f}: {sorted(lengths)}")

    aln_len = next(iter(lengths))

    if not seqs:
        for k in local:
            seqs[k] = []
    else:
        if set(local) != set(seqs):
            missing = sorted(set(seqs) - set(local))
            extra = sorted(set(local) - set(seqs))
            raise ValueError(f"{f}: missing={missing}, extra={extra}")

    for k in seqs:
        seqs[k].append(local[k])

    end = start + aln_len - 1
    gene_name = f.name
    partitions.append((gene_name, start, end))
    start = end + 1

with output.open("w") as out:
    for taxon, parts in seqs.items():
        s = "".join(parts)
        out.write(f">{taxon}\n")
        for i in range(0, len(s), 80):
            out.write(s[i:i+80] + "\n")

with partitions_out.open("w") as out:
    for gene_name, start, end in partitions:
        out.write(f"LG, {gene_name} = {start}-{end}\n")

print(f"Wrote {output}")
print(f"Wrote {partitions_out}")
print(f"Taxa: {len(seqs)}")
print(f"Sites: {sum(end - start + 1 for _, start, end in partitions)}")
print(f"Loci: {len(partitions)}")
PYCODE

echo
echo "[Step 20] Starting accessory genome analysis."

echo "Shared accessory orthogroups are defined as non-core orthogroups present in at least two strains."

python3 "$SCRIPTS_DIR/extract_shared_accessory.py" \
    -i "$GENECOUNT_TSV" \
    -o "$OUTDIR/accessory_results"

SHARED_ACCESSORY_TSV="$OUTDIR/accessory_results/shared_accessory_presence_absence.tsv"

if [[ -f "$SHARED_ACCESSORY_TSV" ]]; then
    SHARED_ACCESSORY_COUNT=$(( $(wc -l < "$SHARED_ACCESSORY_TSV") - 1 ))
else
    SHARED_ACCESSORY_COUNT="NA"
fi

echo "Number of shared accessory orthogroups: $SHARED_ACCESSORY_COUNT"

echo
echo "[Step 21] Visualizing shared accessory genome patterns."

python3 "$SCRIPTS_DIR/analyze_shared_accessory.py" \
    -i "$OUTDIR/accessory_results/shared_accessory_presence_absence.tsv" \
    -o "$OUTDIR/accessory_analysis"

echo
echo "=============================================="
echo "Pipeline finished successfully."
echo "Version: $PIPELINE_VERSION"
echo "=============================================="

echo
echo "Key analysis summary:"
echo

echo "  Number of input proteomes:"
echo "    $SAMPLE_COUNT"

echo
echo "  Number of single-copy orthologues:"
echo "    $SCO_COUNT"

echo
echo "  Number of trimmed alignments remaining after sequence-loss filtering:"
echo "    $N_AFTER_SEQ_LOSS_FILTER"

echo
echo "  Number of trimmed alignments used for concatenation:"
echo "    $N_TRIMMED"

echo
echo "  Number of shared accessory orthogroups:"
echo "    $SHARED_ACCESSORY_COUNT"

echo
echo "Important output files and directories:"
echo

echo "1. Main pipeline log"
echo "   File:"
echo "     $PIPELINE_LOG"
echo "   Meaning:"
echo "     Complete log of the pipeline screen output."

echo
echo "2. Simplified proteome FASTA files"
echo "   Directory:"
echo "     $OUTDIR/proteomes"
echo "   Meaning:"
echo "     Protein FASTA files with simplified headers used as OrthoFinder input."

echo
echo "3. ID mapping files"
echo "   Directory:"
echo "     $OUTDIR/idmap"
echo "   Meaning:"
echo "     Mapping tables between simplified protein IDs and original protein headers."

echo
echo "4. OrthoFinder results"
echo "   Directory:"
echo "     $ORTHOFINDER_RESULTS"
echo "   Meaning:"
echo "     Main OrthoFinder output directory."

echo
echo "5. Orthogroups table"
echo "   File:"
echo "     $ORTHO_TSV"
echo "   Meaning:"
echo "     Orthogroup membership table listing genes assigned to each orthogroup."

echo
echo "6. Orthogroup gene count table"
echo "   File:"
echo "     $GENECOUNT_TSV"
echo "   Meaning:"
echo "     Gene count matrix for orthogroups across samples. This is the main input for accessory genome and group-specific analyses."

echo
echo "7. Single-copy orthologue list"
echo "   File:"
echo "     $SCO_LIST"
echo "   Meaning:"
echo "     List of single-copy core orthogroups used for core genome phylogenetic analysis."

echo
echo "8. Single-copy orthologue FASTA files"
echo "   Directory:"
echo "     $OUTDIR/phylo/OG_seqs"
echo "   Meaning:"
echo "     Protein FASTA files for each single-copy core orthogroup."

echo
echo "9. Trimmed alignments"
echo "   Directory:"
echo "     $OUTDIR/phylo/trimmed_auto1"
echo "   Meaning:"
echo "     trimAl-trimmed alignments used for concatenation."

echo
echo "10. Alignments removed due to sequence loss"
echo "    Directory:"
echo "      $OUTDIR/phylo/trimmed_auto1_sequence_removed"
echo "    Meaning:"
echo "      Alignments excluded because trimAl removed one or more sequences, typically because a sequence became composed only of gaps."

echo
echo "11. Removed short alignments"
echo "    Directory:"
echo "      $OUTDIR/phylo/trimmed_auto1_short_removed"
echo "    Meaning:"
echo "      Alignments shorter than 30 amino acids after trimming. These were excluded from concatenation."

echo
echo "12. Concatenated core protein alignment"
echo "    File:"
echo "      $OUTDIR/phylo/concat_all.fa"
echo "    Meaning:"
echo "      Concatenated alignment of single-copy core orthologues. This file can be used for SplitsTree, ModelTest-NG, and IQ-TREE."

echo
echo "13. Partition file"
echo "    File:"
echo "      $OUTDIR/phylo/concat_partitions.txt"
echo "    Meaning:"
echo "      Partition information corresponding to the concatenated core protein alignment."

echo
echo "14. Shared accessory presence/absence table"
echo "    File:"
echo "      $OUTDIR/accessory_results/shared_accessory_presence_absence.tsv"
echo "    Meaning:"
echo "      Presence/absence matrix of non-core orthogroups present in at least two samples."

echo
echo "15. Accessory genome visualizations"
echo "    Directory:"
echo "      $OUTDIR/accessory_analysis"
echo "    Meaning:"
echo "      Visualization and clustering results based on shared accessory orthogroup presence/absence patterns."
echo
echo "    Main output files:"
echo "      01_presence_distribution.png"
echo "        Distribution of shared accessory orthogroups by the number of strains in which they are present."
echo
echo "      02_presence_absence_heatmap_reordered.png"
echo "        Binary presence/absence heatmap reordered according to accessory genome similarity."
echo
echo "      03_jaccard_distance_matrix_reordered.png"
echo "        Reordered Jaccard distance matrix based on shared accessory orthogroup profiles."
echo
echo "      04_hierarchical_clustering.png"
echo "        Hierarchical clustering dendrogram based on Jaccard distances among strains."
echo
echo "      05_dendrogram_heatmap_combined.png"
echo "        Combined hierarchical clustering dendrogram and binary presence/absence heatmap."

echo
echo "16. MAFFT log"
echo "    File:"
echo "      $OUTDIR/logs/mafft.log"
echo "    Meaning:"
echo "      Combined MAFFT log. MAFFT messages are not printed to the screen."

echo
echo "17. seqkit extraction log"
echo "    File:"
echo "      $OUTDIR/logs/seqkit_extract_orthogroups.log"
echo "    Meaning:"
echo "      seqkit messages from single-copy orthologue sequence extraction."

echo
echo "Downstream analysis examples:"
echo

echo "A. IQ-TREE version check:"
echo "  iqtree3 --version"
echo "  # If iqtree3 is unavailable, try:"
echo "  iqtree2 --version"
echo "  iqtree --version"

echo
echo "B. Standard IQ-TREE analysis using ModelFinder:"
echo "  iqtree3 \\"
echo "    -s $OUTDIR/phylo/concat_all.fa \\"
echo "    -st AA \\"
echo "    -m MFP \\"
echo "    -T AUTO \\"
echo "    -B 1000 \\"
echo "    --alrt 1000 \\"
echo "    --prefix $OUTDIR/phylo/iqtree_concat"

echo
echo "C. If IQ-TREE model search is too slow, first run ModelTest-NG and then run IQ-TREE with a fixed model."

echo
echo "ModelTest-NG version check:"
echo "  modeltest-ng --version"

echo
echo "IQ-TREE version check:"
echo "  iqtree3 --version"

echo
echo "ModelTest-NG example:"
echo "  modeltest-ng \\"
echo "    -i $OUTDIR/phylo/concat_all.fa \\"
echo "    -d aa \\"
echo "    -p $THREADS"

echo
echo "Then run IQ-TREE with the selected fixed model. For example, if LG+G4 is selected:"
echo "  iqtree3 \\"
echo "    -s $OUTDIR/phylo/concat_all.fa \\"
echo "    -st AA \\"
echo "    -m LG+G4 \\"
echo "    -T AUTO \\"
echo "    -B 1000 \\"
echo "    --alrt 1000 \\"
echo "    --prefix $OUTDIR/phylo/iqtree_concat_LG_G4"

echo
echo "Group-specific orthogroup extraction:"
echo
echo "First, prepare GROUP_SAMPLES.txt as a plain text file containing one sample name per line."
echo "The sample names must match the proteome FASTA basenames used in OrthoFinder."
echo
echo "Example GROUP_SAMPLES.txt:"
echo "  BCC7051"
echo "  BL18"
echo "  CBS_466.91"
echo
echo "Example command to create GROUP_SAMPLES.txt:"
echo "  cat > GROUP_SAMPLES.txt << 'EOF'"
echo "  BCC7051"
echo "  BL18"
echo "  CBS_466.91"
echo "  EOF"
echo
echo "Group-specific orthogroup extraction example:"
echo "  python3 $SCRIPTS_DIR/extract_group_presence_absence.py \\"
echo "    -i $GENECOUNT_TSV \\"
echo "    -g GROUP_SAMPLES.txt \\"
echo "    -o group_biased_presence_absence.tsv \\"
echo "    --group-min-pct 90 \\"
echo "    --group-max-pct 100 \\"
echo "    --nongroup-min-pct 0 \\"
echo "    --nongroup-max-pct 10"
echo
echo "This example extracts orthogroups present in 90-100% of the focal group and 0-10% of the non-group samples."

echo
echo "Representative sequence extraction example:"
echo "  python3 $SCRIPTS_DIR/extract_representative_sequences.py \\"
echo "    -i group_biased_presence_absence.tsv \\"
echo "    -s $OG_SEQ_DIR \\"
echo "    -o group_biased_representatives.fa"

echo
echo "Done."

#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Perform PCA on an accessory orthogroup presence/absence matrix "
            "and export R-ready files for visualization."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input shared_accessory_presence_absence.tsv"
    )

    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "-g", "--group-samples",
        default=None,
        help=(
            "Optional file containing focal group sample names, one per line. "
            "If provided, samples are labeled as group or non-group."
        )
    )

    parser.add_argument(
        "--group-name",
        default="Group_A",
        help="Name of the focal group used in metadata. Default: Group_A"
    )

    parser.add_argument(
        "--nongroup-name",
        default="Non_group_A",
        help="Name of non-focal samples used in metadata. Default: Non_group_A"
    )

    parser.add_argument(
        "--standardize",
        action="store_true",
        help=(
            "Standardize each orthogroup column to unit variance before PCA. "
            "Default is centered-only PCA, which is usually suitable for binary "
            "presence/absence data."
        )
    )

    return parser.parse_args()


def read_group_samples(group_file):
    """
    Read focal group sample names from a plain text file.

    The file should contain one sample name per line.
    Empty lines and lines beginning with '#' are ignored.
    """
    group_samples = set()

    with open(group_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                group_samples.add(line)

    return group_samples


def load_presence_absence(input_path):
    """
    Load shared accessory orthogroup presence/absence table.

    Expected input format:
        Orthogroup    n_strains_present    category    strain1    strain2 ...

    The strain columns should contain 0/1 values.

    Output matrix:
        rows = samples/strains
        columns = orthogroups
        values = 0/1
    """
    df = pd.read_csv(input_path, sep="\t")

    required_cols = {"Orthogroup", "n_strains_present", "category"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            "Input file must contain columns: "
            "Orthogroup, n_strains_present, category"
        )

    strain_cols = [
        c for c in df.columns
        if c not in ["Orthogroup", "n_strains_present", "category", "Total"]
    ]

    if len(strain_cols) < 3:
        raise ValueError("At least three strains are required for PCA.")

    pa = df[strain_cols].copy()
    pa = pa.astype(int)

    # Convert from orthogroup x strain to strain x orthogroup.
    matrix = pa.T
    matrix.index.name = "Sample"
    matrix.columns = df["Orthogroup"].astype(str).values

    return df, matrix


def run_pca(matrix, standardize=False):
    """
    Run PCA using numpy singular value decomposition.

    Input:
        matrix:
            rows = samples
            columns = orthogroups
            values = 0/1

    By default, columns are centered but not scaled.
    For binary presence/absence data, centered-only PCA is often a reasonable
    first choice because very rare and very common orthogroups are not
    over-weighted by scaling.

    If --standardize is used, columns are centered and scaled to unit variance.
    """
    X = matrix.values.astype(float)

    # Center columns.
    col_means = X.mean(axis=0)
    X_centered = X - col_means

    if standardize:
        col_sd = X_centered.std(axis=0, ddof=1)
        col_sd[col_sd == 0] = 1.0
        X_used = X_centered / col_sd
    else:
        X_used = X_centered

    # Remove zero-variance columns.
    variances = X_used.var(axis=0, ddof=1)
    keep = variances > 0

    X_used = X_used[:, keep]
    kept_features = matrix.columns[keep]

    if X_used.shape[1] < 2:
        raise ValueError("Fewer than two variable orthogroups remain after filtering.")

    # PCA by SVD.
    U, S, Vt = np.linalg.svd(X_used, full_matrices=False)

    eigenvalues = (S ** 2) / (X_used.shape[0] - 1)
    variance_ratio = eigenvalues / eigenvalues.sum()

    scores_array = U * S
    loadings_array = Vt.T

    n_pcs = scores_array.shape[1]
    pc_names = [f"PC{i}" for i in range(1, n_pcs + 1)]

    scores = pd.DataFrame(
        scores_array,
        index=matrix.index,
        columns=pc_names
    )
    scores.index.name = "Sample"

    loadings = pd.DataFrame(
        loadings_array,
        index=kept_features,
        columns=pc_names
    )
    loadings.index.name = "Orthogroup"

    variance_df = pd.DataFrame({
        "PC": pc_names,
        "eigenvalue": eigenvalues,
        "variance_explained": variance_ratio,
        "variance_explained_percent": variance_ratio * 100.0,
        "cumulative_variance_explained": np.cumsum(variance_ratio),
        "cumulative_variance_explained_percent": np.cumsum(variance_ratio) * 100.0
    })

    return scores, loadings, variance_df, kept_features


def make_metadata(samples, group_samples=None, group_name="Group_A", nongroup_name="Non_group_A"):
    """
    Create sample metadata for R visualization.
    """
    metadata = pd.DataFrame({"Sample": samples})

    if group_samples is None:
        metadata["Group"] = "All_samples"
    else:
        metadata["Group"] = metadata["Sample"].apply(
            lambda x: group_name if x in group_samples else nongroup_name
        )

    return metadata


def write_r_script(outdir):
    """
    Write an R script for non-interactive PCA visualization using ggplot2.
    """
    r_script = outdir / "plot_accessory_pca.R"

    r_code = r'''#!/usr/bin/env Rscript

# PCA visualization for EukPan accessory orthogroup presence/absence analysis

scores <- read.delim("accessory_pca_scores.tsv", header = TRUE, sep = "\t", check.names = FALSE)
variance <- read.delim("accessory_pca_variance.tsv", header = TRUE, sep = "\t", check.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required. Install it with: install.packages('ggplot2')")
}

library(ggplot2)

pc1_var <- round(variance$variance_explained_percent[variance$PC == "PC1"], 2)
pc2_var <- round(variance$variance_explained_percent[variance$PC == "PC2"], 2)

p <- ggplot(scores, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of accessory orthogroup presence/absence",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Group"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("accessory_pca_PC1_PC2.pdf", p, width = 7, height = 6)
ggsave("accessory_pca_PC1_PC2.png", p, width = 7, height = 6, dpi = 300)

message("Saved:")
message("  accessory_pca_PC1_PC2.pdf")
message("  accessory_pca_PC1_PC2.png")
'''

    with open(r_script, "w") as f:
        f.write(r_code)

    return r_script


def write_interactive_r_instructions(outdir):
    """
    Write an R text file containing commands for interactive visualization
    using identify().
    """
    instruction_file = outdir / "interactive_R_identify_commands.txt"

    text = r'''# Interactive PCA visualization in R
# Run these commands after starting R in the accessory_pca output directory.

scores <- read.delim("accessory_pca_scores.tsv", header = TRUE, sep = "\t", check.names = FALSE)
variance <- read.delim("accessory_pca_variance.tsv", header = TRUE, sep = "\t", check.names = FALSE)

pc1.var <- round(variance$variance_explained_percent[variance$PC == "PC1"], 2)
pc2.var <- round(variance$variance_explained_percent[variance$PC == "PC2"], 2)

group.factor <- as.factor(scores$Group)
group.cols <- as.numeric(group.factor)

plot(
  scores$PC1,
  scores$PC2,
  col = group.cols,
  pch = 19,
  xlab = paste0("PC1 (", pc1.var, "%)"),
  ylab = paste0("PC2 (", pc2.var, "%)"),
  main = "PCA of accessory orthogroup presence/absence"
)

legend(
  "topright",
  legend = levels(group.factor),
  col = seq_along(levels(group.factor)),
  pch = 19,
  bty = "n"
)

# Click points to show sample names.
identify(
  scores$PC1,
  scores$PC2,
  labels = scores$Sample
)

# To save the interactively labeled plot after identify(), run:
dev.copy(pdf, "accessory_pca_PC1_PC2_identified.pdf", width = 7, height = 6)
dev.off()
'''

    with open(instruction_file, "w") as f:
        f.write(text)

    return instruction_file


def main():
    args = parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)

    if not input_path.exists():
        sys.stderr.write(f"ERROR: Input file not found: {input_path}\n")
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    df, matrix = load_presence_absence(input_path)

    group_samples = None
    if args.group_samples is not None:
        group_file = Path(args.group_samples)
        if not group_file.exists():
            sys.stderr.write(f"ERROR: Group sample file not found: {group_file}\n")
            sys.exit(1)
        group_samples = read_group_samples(group_file)

    scores, loadings, variance_df, kept_features = run_pca(
        matrix,
        standardize=args.standardize
    )

    metadata = make_metadata(
        samples=scores.index.tolist(),
        group_samples=group_samples,
        group_name=args.group_name,
        nongroup_name=args.nongroup_name
    )

    scores_out = scores.reset_index().merge(metadata, on="Sample", how="left")
    matrix_out = matrix.reset_index()

    scores_out.to_csv(outdir / "accessory_pca_scores.tsv", sep="\t", index=False)
    loadings.reset_index().to_csv(outdir / "accessory_pca_loadings.tsv", sep="\t", index=False)
    variance_df.to_csv(outdir / "accessory_pca_variance.tsv", sep="\t", index=False)
    metadata.to_csv(outdir / "accessory_pca_metadata.tsv", sep="\t", index=False)
    matrix_out.to_csv(outdir / "accessory_presence_absence_matrix_for_R.tsv", sep="\t", index=False)

    kept_df = pd.DataFrame({"Orthogroup": kept_features})
    kept_df.to_csv(outdir / "accessory_pca_used_orthogroups.tsv", sep="\t", index=False)

    r_script = write_r_script(outdir)
    interactive_r_instructions = write_interactive_r_instructions(outdir)

    sys.stderr.write("Finished accessory PCA analysis.\n")
    sys.stderr.write(f"Input: {input_path}\n")
    sys.stderr.write(f"Output directory: {outdir}\n")
    sys.stderr.write("\n")

    sys.stderr.write("Main output files for R:\n")
    sys.stderr.write("  accessory_pca_scores.tsv\n")
    sys.stderr.write("  accessory_pca_loadings.tsv\n")
    sys.stderr.write("  accessory_pca_variance.tsv\n")
    sys.stderr.write("  accessory_pca_metadata.tsv\n")
    sys.stderr.write("  accessory_presence_absence_matrix_for_R.tsv\n")
    sys.stderr.write("  accessory_pca_used_orthogroups.tsv\n")
    sys.stderr.write("  plot_accessory_pca.R\n")
    sys.stderr.write("  interactive_R_identify_commands.txt\n")
    sys.stderr.write("\n")

    sys.stderr.write("To visualize PCA using the included R script, run:\n")
    sys.stderr.write(f"  cd {outdir}\n")
    sys.stderr.write("  Rscript plot_accessory_pca.R\n")
    sys.stderr.write("\n")

    sys.stderr.write("The R script will generate:\n")
    sys.stderr.write("  accessory_pca_PC1_PC2.pdf\n")
    sys.stderr.write("  accessory_pca_PC1_PC2.png\n")
    sys.stderr.write("\n")

    sys.stderr.write("To inspect points interactively in R using identify(), run:\n")
    sys.stderr.write(f"  cd {outdir}\n")
    sys.stderr.write("  R\n")
    sys.stderr.write("\n")
    sys.stderr.write("Then copy and paste the commands shown in:\n")
    sys.stderr.write("  interactive_R_identify_commands.txt\n")
    sys.stderr.write("\n")
    sys.stderr.write("In R, the identify() command lets you click points and display sample names.\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import sys
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from matplotlib import gridspec
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze shared accessory orthogroups from a presence/absence table.\n"
            "Input columns:\n"
            "  Orthogroup, n_strains_present, category, strain1, strain2, ...\n"
            "where strain columns are 0/1."
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
    return parser.parse_args()


def get_label_fontsize(n_samples: int) -> float:
    """
    Adaptive label font size for heatmaps and distance matrices.
    """
    if n_samples <= 10:
        return 12
    if n_samples <= 20:
        return 11
    if n_samples <= 40:
        return 9.5
    if n_samples <= 80:
        return 8
    if n_samples <= 150:
        return 7
    if n_samples <= 250:
        return 5.5
    return 4.5


def get_dendrogram_fontsize(n_samples: int) -> float:
    """
    Adaptive label font size for dendrograms.
    """
    if n_samples <= 10:
        return 12
    if n_samples <= 20:
        return 11
    if n_samples <= 40:
        return 9
    if n_samples <= 80:
        return 7.5
    if n_samples <= 150:
        return 6.5
    if n_samples <= 250:
        return 5
    return 4


def get_heatmap_figsize(n_samples: int, n_orthogroups: int):
    """
    Adaptive figure size for the presence/absence heatmap.
    """
    width = max(12, min(46, n_samples * 0.24))
    height = max(8, min(40, n_orthogroups * 0.012))
    return width, height


def get_matrix_figsize(n_samples: int):
    """
    Adaptive figure size for the Jaccard distance matrix.

    The matrix panel itself should be square.
    Extra horizontal space is added for the colorbar.
    """
    matrix_size = max(9, min(34, n_samples * 0.20))
    return matrix_size + 2.0, matrix_size


def get_dendrogram_figsize(n_samples: int):
    """
    Adaptive figure size for the dendrogram.
    """
    width = max(12, min(46, n_samples * 0.24))
    return width, 7


def get_combined_figsize(n_samples: int, n_orthogroups: int):
    """
    Adaptive figure size for the combined dendrogram and heatmap.
    """
    width = max(14, min(50, n_samples * 0.26))
    height = max(10, min(44, n_orthogroups * 0.012 + 5))
    return width, height


def save_png_and_emf(fig, png_path: Path, dpi: int = 300):
    """
    Save a matplotlib figure as PNG, SVG, and EMF if Inkscape is available.

    EMF export is optional. If Inkscape is not installed, PNG and SVG are still saved.
    """
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")

    svg_path = png_path.with_suffix(".svg")
    emf_path = png_path.with_suffix(".emf")

    fig.savefig(svg_path, bbox_inches="tight")

    try:
        subprocess.run(
            [
                "inkscape",
                str(svg_path),
                "--export-filename",
                str(emf_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except Exception as e:
        print(
            f"WARNING: EMF export failed for {png_path.name}: {e}",
            file=sys.stderr
        )

    plt.close(fig)


def load_presence_absence_table(input_path: Path):
    """
    Load a shared accessory presence/absence table.

    Expected columns:
      Orthogroup, n_strains_present, category, strain1, strain2, ...

    Strain columns must contain 0/1 values.
    """
    df = pd.read_csv(input_path, sep="\t")

    required_cols = {"Orthogroup", "n_strains_present", "category"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"Input table must contain columns: {sorted(required_cols)}"
        )

    strain_cols = [
        c for c in df.columns
        if c not in ["Orthogroup", "n_strains_present", "category", "Total"]
    ]

    if len(strain_cols) < 2:
        raise ValueError("Need at least two strain columns.")

    pa = df[strain_cols].copy()
    for c in strain_cols:
        pa[c] = pa[c].astype(int)

    return df, pa, strain_cols


def save_summary(df: pd.DataFrame, pa: pd.DataFrame, strain_cols, outdir: Path):
    """
    Save a simple summary table.
    """
    summary_rows = []
    summary_rows.append(["n_shared_accessory_orthogroups", len(df)])
    summary_rows.append(["n_strains", len(strain_cols)])

    per_strain = pa.sum(axis=0)
    for strain, n in per_strain.items():
        summary_rows.append([f"n_present_in_{strain}", int(n)])

    summary_df = pd.DataFrame(summary_rows, columns=["metric", "value"])
    summary_df.to_csv(outdir / "00_summary.tsv", sep="\t", index=False)


def save_distribution(df: pd.DataFrame, outdir: Path):
    """
    Save the distribution of the number of strains in which each orthogroup is present.
    """
    dist = (
        df["n_strains_present"]
        .value_counts()
        .sort_index()
        .rename_axis("n_strains_present")
        .reset_index(name="n_orthogroups")
    )

    dist.to_csv(outdir / "01_presence_distribution.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(
        dist["n_strains_present"],
        dist["n_orthogroups"],
        color="black"
    )
    ax.set_xlabel("Number of strains present")
    ax.set_ylabel("Number of shared accessory orthogroups")
    ax.set_title("Distribution of shared accessory orthogroups")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    save_png_and_emf(fig, outdir / "01_presence_distribution.png")


def compute_jaccard(pa: pd.DataFrame):
    """
    Compute pairwise Jaccard distances and similarities among strains.

    Rows of the input matrix are orthogroups and columns are strains.
    The matrix is transposed so that rows represent strains.
    """
    strain_matrix = pa.T.values
    dist_vec = pdist(strain_matrix, metric="jaccard")
    dist_mat = squareform(dist_vec)
    sim_mat = 1.0 - dist_mat
    return dist_mat, sim_mat


def compute_clustering(dist_mat, strain_cols):
    """
    Perform average-linkage hierarchical clustering using Jaccard distances.
    """
    dist_vec = squareform(dist_mat, checks=False)
    Z = linkage(dist_vec, method="average", optimal_ordering=True)

    d = dendrogram(Z, labels=strain_cols, no_plot=True)
    ordered_cols = d["ivl"]

    return Z, ordered_cols


def reorder_rows_for_heatmap(pa_ord: pd.DataFrame):
    """
    Reorder orthogroup rows for the heatmap.

    Rows are sorted by:
      1. Number of strains in which the orthogroup is present
      2. Binary presence/absence pattern
    """
    row_patterns = pa_ord.astype(str).agg("".join, axis=1)

    tmp = pd.DataFrame({
        "row_index": np.arange(pa_ord.shape[0]),
        "presence_sum": pa_ord.sum(axis=1).values,
        "pattern": row_patterns.values
    })

    tmp = tmp.sort_values(
        by=["presence_sum", "pattern"],
        ascending=[False, False],
        kind="mergesort"
    )

    return pa_ord.iloc[tmp["row_index"].values]


def save_distance_matrix(dist_mat, strain_cols, ordered_cols, outdir: Path):
    """
    Save Jaccard distance and similarity matrices.

    The plotted reordered distance matrix is drawn so that the matrix panel itself
    is square. The colorbar is placed outside the matrix panel.
    """
    n_samples = len(strain_cols)
    label_fontsize = get_label_fontsize(n_samples)

    dist_df = pd.DataFrame(dist_mat, index=strain_cols, columns=strain_cols)
    sim_df = 1.0 - dist_df

    dist_df.to_csv(outdir / "03_jaccard_distance_matrix.tsv", sep="\t")
    sim_df.to_csv(outdir / "03_jaccard_similarity_matrix.tsv", sep="\t")

    dist_df_ord = dist_df.loc[ordered_cols, ordered_cols]
    sim_df_ord = sim_df.loc[ordered_cols, ordered_cols]

    dist_df_ord.to_csv(outdir / "03_jaccard_distance_matrix_reordered.tsv", sep="\t")
    sim_df_ord.to_csv(outdir / "03_jaccard_similarity_matrix_reordered.tsv", sep="\t")

    fig = plt.figure(figsize=get_matrix_figsize(n_samples))
    ax = fig.add_subplot(111)

    im = ax.imshow(
        dist_df_ord.values,
        aspect="equal",
        interpolation="nearest",
        cmap="viridis"
    )

    ax.set_adjustable("box")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.15)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Jaccard distance")

    ax.set_xticks(range(len(ordered_cols)))
    ax.set_xticklabels(
        ordered_cols,
        rotation=90,
        fontsize=label_fontsize
    )
    ax.set_yticks(range(len(ordered_cols)))
    ax.set_yticklabels(
        ordered_cols,
        fontsize=label_fontsize
    )

    ax.set_title("Jaccard distance matrix (reordered by dendrogram)")

    plt.tight_layout()
    save_png_and_emf(fig, outdir / "03_jaccard_distance_matrix_reordered.png")

    return dist_df, sim_df, dist_df_ord, sim_df_ord


def save_clustering(Z, strain_cols, outdir: Path):
    """
    Save a hierarchical clustering dendrogram.
    """
    n_samples = len(strain_cols)
    dendro_fontsize = get_dendrogram_fontsize(n_samples)

    fig = plt.figure(figsize=get_dendrogram_figsize(n_samples))
    ax = fig.add_subplot(111)

    dendrogram(
        Z,
        labels=strain_cols,
        leaf_rotation=90,
        leaf_font_size=dendro_fontsize,
        ax=ax
    )

    ax.set_ylabel("Jaccard distance")
    ax.set_title("Hierarchical clustering based on shared accessory orthogroups")

    plt.tight_layout()
    save_png_and_emf(fig, outdir / "04_hierarchical_clustering.png")


def save_heatmap(pa: pd.DataFrame, ordered_cols, outdir: Path):
    """
    Save a reordered binary presence/absence heatmap.
    """
    n_samples = len(ordered_cols)

    pa_ord = pa[ordered_cols].copy()
    pa_ord = reorder_rows_for_heatmap(pa_ord)

    n_orthogroups = pa_ord.shape[0]
    label_fontsize = get_label_fontsize(n_samples)

    cmap = ListedColormap(["white", "black"])

    fig = plt.figure(figsize=get_heatmap_figsize(n_samples, n_orthogroups))
    ax = fig.add_subplot(111)

    ax.imshow(
        pa_ord.values,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1
    )

    ax.set_xticks(range(len(ordered_cols)))
    ax.set_xticklabels(
        ordered_cols,
        rotation=90,
        fontsize=label_fontsize
    )
    ax.set_yticks([])

    ax.set_xlabel("Strains")
    ax.set_ylabel("Shared accessory orthogroups")
    ax.set_title("Presence/absence heatmap of shared accessory orthogroups")

    plt.tight_layout()
    save_png_and_emf(fig, outdir / "02_presence_absence_heatmap_reordered.png")

    return pa_ord


def save_combined_dendrogram_heatmap(Z, pa_ord: pd.DataFrame, strain_cols, outdir: Path):
    """
    Save a combined dendrogram and presence/absence heatmap.
    """
    n_samples = len(strain_cols)
    n_orthogroups = pa_ord.shape[0]

    label_fontsize = get_label_fontsize(n_samples)
    dendro_fontsize = get_dendrogram_fontsize(n_samples)

    d = dendrogram(Z, labels=strain_cols, no_plot=True)
    dendro_order = d["ivl"]

    pa_for_plot = pa_ord[dendro_order]

    cmap = ListedColormap(["white", "black"])

    fig = plt.figure(figsize=get_combined_figsize(n_samples, n_orthogroups))
    gs = gridspec.GridSpec(
        nrows=2,
        ncols=1,
        height_ratios=[3, 12],
        hspace=0.05
    )

    ax_d = fig.add_subplot(gs[0, 0])
    dendrogram(
        Z,
        labels=strain_cols,
        leaf_rotation=90,
        leaf_font_size=dendro_fontsize,
        ax=ax_d
    )
    ax_d.set_ylabel("Jaccard distance")
    ax_d.set_title("Dendrogram and presence/absence heatmap")
    ax_d.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    ax_h = fig.add_subplot(gs[1, 0])
    ax_h.imshow(
        pa_for_plot.values,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1
    )

    ax_h.set_xticks(range(len(dendro_order)))
    ax_h.set_xticklabels(
        dendro_order,
        rotation=90,
        fontsize=label_fontsize
    )
    ax_h.set_yticks([])

    ax_h.set_xlabel("Strains")
    ax_h.set_ylabel("Shared accessory orthogroups")

    fig.subplots_adjust(hspace=0.05)
    save_png_and_emf(fig, outdir / "05_dendrogram_heatmap_combined.png")


def main():
    args = parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)

    if not input_path.exists():
        sys.stderr.write(f"ERROR: Input file not found: {input_path}\n")
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)

    df, pa, strain_cols = load_presence_absence_table(input_path)

    save_summary(df, pa, strain_cols, outdir)
    save_distribution(df, outdir)

    dist_mat, sim_mat = compute_jaccard(pa)
    Z, ordered_cols = compute_clustering(dist_mat, strain_cols)

    save_distance_matrix(
        dist_mat,
        strain_cols,
        ordered_cols,
        outdir
    )

    save_clustering(Z, strain_cols, outdir)

    pa_ord = save_heatmap(pa, ordered_cols, outdir)
    save_combined_dendrogram_heatmap(Z, pa_ord, strain_cols, outdir)

    sys.stderr.write("Finished.\n")
    sys.stderr.write(f"Input: {input_path}\n")
    sys.stderr.write(f"Output directory: {outdir}\n")


if __name__ == "__main__":
    main()

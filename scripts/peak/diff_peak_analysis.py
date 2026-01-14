#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from pybedtools import BedTool
import matplotlib.pyplot as plt
import seaborn as sns


# =================================================
# 1. Peak loading
# =================================================

def load_peaks(narrowpeak_file: str) -> pd.DataFrame:
    """
    Load narrowPeak file and generate unique peak_id
    """
    peaks = BedTool(narrowpeak_file)
    df = peaks.to_dataframe(names=[
        "chrom", "start", "end", "name", "score", "strand",
        "signal", "pval", "qval", "peak"
    ])

    df["peak_id"] = (
        df["chrom"] + ":" +
        df["start"].astype(str) + "-" +
        df["end"].astype(str)
    )
    return df


# =================================================
# 2. BAM quantification
# =================================================

def count_reads_in_peaks(
    bam_file: str,
    peaks_df: pd.DataFrame
) -> np.ndarray:
    """
    Count reads overlapping each peak
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    counts = []

    for _, row in peaks_df.iterrows():
        counts.append(
            bam.count(
                contig=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"])
            )
        )

    bam.close()
    return np.array(counts)


# =================================================
# 3. IP / Input normalization
# =================================================

def normalize_ip(
    ip_counts: np.ndarray,
    input_counts: np.ndarray | None = None,
    pseudocount: float = 1.0
) -> np.ndarray:
    """
    Normalize IP signal by Input if provided
    """
    if input_counts is None:
        return ip_counts + pseudocount

    return (ip_counts + pseudocount) / (input_counts + pseudocount)


# =================================================
# 4. log2 fold-change
# =================================================

def compute_log2fc(
    ko_signal: np.ndarray,
    wt_signal: np.ndarray
) -> np.ndarray:
    """
    Compute log2(KO / WT)
    """
    return np.log2(ko_signal / wt_signal)


# =================================================
# 5. Output tables
# =================================================

def export_tables(df: pd.DataFrame, outdir: Path):
    """
    Export result tables
    """
    outdir.mkdir(parents=True, exist_ok=True)

    df[[
        "peak_id", "chrom", "start", "end",
        "KO_IP", "WT_IP",
        "KO_norm", "WT_norm",
        "log2FC"
    ]].to_csv(
        outdir / "KO_vs_WT_peak_quantification.tsv",
        sep="\t",
        index=False
    )

    df[[
        "peak_id", "KO_IP", "WT_IP"
    ]].to_csv(
        outdir / "peak_counts_for_DESeq2.tsv",
        sep="\t",
        index=False
    )


# =================================================
# 6. Visualization
# =================================================

def plot_ma(df: pd.DataFrame, outdir: Path):
    A = np.log2(df["KO_norm"] + df["WT_norm"])
    M = df["log2FC"]

    plt.figure(figsize=(6, 5))
    plt.scatter(A, M, s=5, alpha=0.4)
    plt.axhline(0, color="black", lw=0.5)
    plt.axhline(1, color="red", linestyle="--")
    plt.axhline(-1, color="blue", linestyle="--")

    plt.xlabel("A = log2(KO + WT)")
    plt.ylabel("M = log2(KO / WT)")
    plt.title("MA plot (KO vs WT)")
    plt.tight_layout()
    plt.savefig(outdir / "MA_plot.png")
    plt.close()


def plot_volcano_like(df: pd.DataFrame, outdir: Path):
    plt.figure(figsize=(6, 5))
    plt.scatter(
        df["log2FC"],
        np.log10(df["KO_norm"] + df["WT_norm"]),
        s=5,
        alpha=0.4
    )

    plt.axvline(1, color="red", linestyle="--")
    plt.axvline(-1, color="blue", linestyle="--")

    plt.xlabel("log2FC (KO / WT)")
    plt.ylabel("log10 signal")
    plt.title("Peak enrichment shift")
    plt.tight_layout()
    plt.savefig(outdir / "Volcano_like_plot.png")
    plt.close()


def plot_heatmap(
    df: pd.DataFrame,
    outdir: Path,
    top_n: int = 2000
):
    df_top = df.sort_values("KO_norm", ascending=False).head(top_n)
    mat = np.log2(df_top[["WT_norm", "KO_norm"]])

    sns.clustermap(
        mat,
        cmap="RdBu_r",
        center=0,
        col_cluster=False,
        row_cluster=True,
        yticklabels=False,
        figsize=(4, 8)
    )

    plt.savefig(outdir / "Heatmap_KO_defined_peaks.png")
    plt.close()


# =================================================
# 7. Main pipeline
# =================================================

def main():
    parser = argparse.ArgumentParser(
        description="Differential peak analysis using KO-defined peaks"
    )
    parser.add_argument("--peaks", required=True, help="KO narrowPeak")
    parser.add_argument("--ko-ip", required=True, help="KO IP BAM")
    parser.add_argument("--wt-ip", required=True, help="WT IP BAM")
    parser.add_argument("--ko-input", help="KO Input BAM (optional)")
    parser.add_argument("--wt-input", help="WT Input BAM (optional)")
    parser.add_argument("-o", "--outdir", default="diff_peak_out")

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    # # Load peaks
    # peak_df = load_peaks(args.peaks)

    # # Count IP reads
    # peak_df["KO_IP"] = count_reads_in_peaks(args.ko_ip, peak_df)
    # peak_df["WT_IP"] = count_reads_in_peaks(args.wt_ip, peak_df)

    # # Input normalization
    # if args.ko_input and args.wt_input:
    #     peak_df["KO_Input"] = count_reads_in_peaks(args.ko_input, peak_df)
    #     peak_df["WT_Input"] = count_reads_in_peaks(args.wt_input, peak_df)

    #     peak_df["KO_norm"] = normalize_ip(
    #         peak_df["KO_IP"].values,
    #         peak_df["KO_Input"].values
    #     )
    #     peak_df["WT_norm"] = normalize_ip(
    #         peak_df["WT_IP"].values,
    #         peak_df["WT_Input"].values
    #     )
    # else:
    #     peak_df["KO_norm"] = normalize_ip(peak_df["KO_IP"].values)
    #     peak_df["WT_norm"] = normalize_ip(peak_df["WT_IP"].values)

    # # log2FC
    # peak_df["log2FC"] = compute_log2fc(
    #     peak_df["KO_norm"].values,
    #     peak_df["WT_norm"].values
    # )

    # Output
    # export_tables(peak_df, outdir)
    peak_df = pd.read_csv("/disk5/luosg/DIPseq20251215/output/Diff/KO_vs_WT_peak_quantification.tsv",sep="\t")
    # Visualization
    plot_ma(peak_df, outdir)
    plot_volcano_like(peak_df, outdir)
    plot_heatmap(peak_df, outdir)

    print("✔ Analysis finished")


if __name__ == "__main__":
    main()

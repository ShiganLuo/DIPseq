import pandas as pd
import numpy as np
import pysam
from pathlib import Path
import argparse

def read_narrowpeak(path: str) -> pd.DataFrame:
    cols = [
        "chrom", "start", "end", "name", "score",
        "strand", "signal", "pvalue", "qvalue", "summit"
    ]
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        usecols=["chrom", "start", "end"]
    )


############################################
# 2. Load genes + promoter (BED semantics)
############################################

def load_genes_from_gtf_fast(
    gtf_path: str,
    promoter_up: int = 2000,
    promoter_down: int = 2000
) -> pd.DataFrame:

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        names=["chrom", "feature", "start", "end", "strand", "attr"]
    )

    gtf = gtf.loc[gtf.feature == "gene"].copy()

    # GTF → BED (0-based, half-open)
    gtf["gene_start"] = gtf["start"] - 1
    gtf["gene_end"] = gtf["end"]

    gene_name = gtf["attr"].str.extract(r'gene_name "([^"]+)"')[0]
    gene_id = gtf["attr"].str.extract(r'gene_id "([^"]+)"')[0]
    gtf["gene"] = gene_name.fillna(gene_id)

    plus = gtf["strand"] == "+"
    minus = ~plus

    gtf["p_start"] = 0
    gtf["p_end"] = 0

    gtf.loc[plus, "p_start"] = (
        gtf.loc[plus, "gene_start"] - promoter_up
    ).clip(lower=0)
    gtf.loc[plus, "p_end"] = (
        gtf.loc[plus, "gene_start"] + promoter_down
    )

    gtf.loc[minus, "p_start"] = (
        gtf.loc[minus, "gene_end"] - promoter_down
    ).clip(lower=0)
    gtf.loc[minus, "p_end"] = (
        gtf.loc[minus, "gene_end"] + promoter_up
    )

    return gtf[[
        "chrom",
        "gene_start",
        "gene_end",
        "p_start",
        "p_end",
        "gene"
    ]]


############################################
# 3. Peak → gene mapping (NO interval merging)
############################################

def assign_peaks_to_genes_fast(
    peaks: pd.DataFrame,
    genes: pd.DataFrame
) -> pd.DataFrame:
    """
    Each row = one peak × one gene
    Peak boundaries are NEVER modified
    """

    out = []

    for chrom in peaks.chrom.unique():
        p = peaks[peaks.chrom == chrom]
        g = genes[genes.chrom == chrom]

        if p.empty or g.empty:
            continue

        for _, peak in p.iterrows():
            hit = g[
                (
                    (peak.start < g.gene_end) &
                    (peak.end > g.gene_start)
                ) |
                (
                    (peak.start < g.p_end) &
                    (peak.end > g.p_start)
                )
            ]

            for _, gene in hit.iterrows():
                region_type = (
                    "promoter"
                    if (peak.start < gene.p_end and peak.end > gene.p_start)
                    else "gene_body"
                )

                out.append({
                    "chrom": chrom,
                    "start": peak.start,
                    "end": peak.end,
                    "gene": gene.gene,
                    "region_type": region_type
                })

    return pd.DataFrame(out)


############################################
# 4. BAM coverage (peak-level)
############################################

def bam_peak_counts(
    bam_path: str,
    peak_gene_df: pd.DataFrame
) -> pd.DataFrame:

    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = []

    for _, r in peak_gene_df.iterrows():
        n = bam.count(
            contig=r.chrom,
            start=int(r.start),
            end=int(r.end)
        )
        counts.append(n)

    bam.close()
    peak_gene_df = peak_gene_df.copy()
    peak_gene_df["reads"] = counts
    return peak_gene_df


############################################
# 5. Core analysis
############################################

def alkbh1_6ma_recovery_analysis(
    ko_peak_file: str,
    wt_peak_file: str,
    ko_bam: str,
    wt_bam: str,
    gtf: str,
    min_log2fc: float = 1.0
) -> pd.DataFrame:

    ko_peaks = read_narrowpeak(ko_peak_file)
    wt_peaks = read_narrowpeak(wt_peak_file)
    genes = load_genes_from_gtf_fast(gtf)

    # KO-only peaks
    ko_only = ko_peaks.merge(
        wt_peaks,
        on=["chrom", "start", "end"],
        how="left",
        indicator=True
    )
    ko_only = ko_only[ko_only["_merge"] == "left_only"]
    ko_only = ko_only[["chrom", "start", "end"]]

    if ko_only.empty:
        raise RuntimeError("No KO-only peaks found")

    # Peak → gene (no merging)
    peak_gene = assign_peaks_to_genes_fast(ko_only, genes)

    if peak_gene.empty:
        raise RuntimeError("No KO-only peaks overlap genes/promoters")

    # BAM counts
    ko_cov = bam_peak_counts(ko_bam, peak_gene)
    wt_cov = bam_peak_counts(wt_bam, peak_gene)

    peak_gene["KO_reads"] = ko_cov["reads"].values
    peak_gene["WT_reads"] = wt_cov["reads"].values

    peak_gene["log2FC_KO_WT"] = np.log2(
        (peak_gene["KO_reads"] + 1) /
        (peak_gene["WT_reads"] + 1)
    )

    # Gene-level summary (DO NOT change peak coordinates)
    gene_summary = (
        peak_gene
        .groupby("gene")
        .agg(
            KO_reads=("KO_reads", "sum"),
            WT_reads=("WT_reads", "sum"),
            KO_only_peaks=("gene", "count"),
            log2FC_KO_WT=("log2FC_KO_WT", "mean")
        )
        .reset_index()
        .sort_values("log2FC_KO_WT", ascending=False)
    )

    return gene_summary[gene_summary.log2FC_KO_WT >= min_log2fc]


############################################
# 6. CLI
############################################

def main():
    parser = argparse.ArgumentParser(
        description="Alkbh1 KO 6mA recovery analysis (peak-safe)"
    )
    parser.add_argument("--ko-peaks", required=True)
    parser.add_argument("--wt-peaks", required=True)
    parser.add_argument("--ko-ip", required=True)
    parser.add_argument("--wt-ip", required=True)
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("-o", "--outdir", default="alkbh1_6ma_out")
    parser.add_argument("--min-log2fc", type=float, default=0.58)

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    result = alkbh1_6ma_recovery_analysis(
        ko_peak_file=args.ko_peaks,
        wt_peak_file=args.wt_peaks,
        ko_bam=args.ko_ip,
        wt_bam=args.wt_ip,
        gtf=args.gtf,
        min_log2fc=args.min_log2fc
    )

    result.to_csv(
        outdir / "Alkbh1_KO_6mA_recovery_genes.tsv",
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import logging
import pandas as pd
import pyranges as pr
import numpy as np
from typing import Union, Literal
from scipy.stats import binom, poisson

Number = Union[int, float]
ArrayLike = Union[pd.Series, Number]

# ============================
# logging
# ============================
def setup_logger(level: str):
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

# ============================
# fast GTF parsing
# ============================
def load_gene_gtf_fast(path: str) -> pr.PyRanges:
    logging.info("Loading gene GTF (gene only)")
    cols = ["Chromosome","Source","Feature","Start","End","Score","Strand","Frame","Attributes"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols, usecols=[0,2,3,4,6,8])
    df = df[df.Feature=="gene"]
    attrs = df.Attributes.str.extract(r'gene_id "([^"]+)".*gene_type "([^"]+)"')
    df["gene_id"] = attrs[0]
    df["gene_type"] = attrs[1]
    df = df.drop(columns=["Attributes","Feature"])
    return pr.PyRanges(df)

def load_line1_gtf_fast(path: str) -> pr.PyRanges:
    logging.info("Loading TE GTF (LINE1 only)")
    cols = ["Chromosome","Source","Feature","Start","End","Score","Strand","Frame","Attributes"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols, usecols=[0,2,3,4,6,8])
    df = df[df.Feature=="exon"]
    df["te_family"] = df.Attributes.str.extract(r'family_id "([^"]+)"')[0]
    df = df[df.te_family=="L1"].drop(columns=["Attributes","Feature","te_family"])
    return pr.PyRanges(df)

# ============================
# genomic utilities
# ============================
def make_promoter(genes: pr.PyRanges, upstream: int) -> pr.PyRanges:
    df = genes.df.copy()
    plus = df.Strand == "+"
    df.loc[plus,"End"] = df.loc[plus,"Start"]
    df.loc[plus,"Start"] -= upstream
    df.loc[~plus,"Start"] = df.loc[~plus,"End"]
    df.loc[~plus,"End"] += upstream
    df["Start"] = df["Start"].clip(lower=0)
    return pr.PyRanges(df)

def enrichment(
    line1_bp: ArrayLike,
    region_bp: ArrayLike,
    total_line1_bp: Number,
    genome_bp: Number
) -> ArrayLike:
    """
    Calculate vectorized enrichment scores for LINE1 elements in genomic regions.

    This function computes the enrichment of LINE1 coverage in a region relative to 
    its expected proportion in the genome.

    Parameters
    ----------
    line1_bp : int, float, or pd.Series
        Observed number of LINE1 base pairs in the region (scalar or vector).
    region_bp : int, float, or pd.Series
        Length of the region in base pairs (scalar or vector). Must match the shape of line1_bp.
    total_line1_bp : int or float
        Total LINE1 base pairs in the genome.
    genome_bp : int or float
        Total base pairs in the genome.

    Returns
    -------
    int, float, or pd.Series
        Enrichment score(s). If `line1_bp` and `region_bp` are Series, returns a Series
        of enrichment scores; otherwise returns a scalar. The enrichment score represents:
        
        enrichment = (line1_bp / total_line1_bp) / (region_bp / genome_bp)

    Notes
    -----
    - An enrichment score > 1 indicates that LINE1 is overrepresented in the region 
      compared to random expectation.
    - An enrichment score < 1 indicates underrepresentation.
    - The function handles division by zero, NaN, and infinite values:
        - Replaces inf or -inf with 0
        - Fills NaN with 0
    - Supports vectorized computation for pd.Series input for efficiency.

    Examples
    --------
    >>> enrichment(50, 1000, 1_000_000, 3_000_000)
    0.15

    >>> enrichment(pd.Series([50, 100]), pd.Series([1000, 2000]), 1_000_000, 3_000_000)
    0    0.15
    1    0.15
    dtype: float64
    """
    if total_line1_bp == 0:
        return 0.0

    ratio = (line1_bp / total_line1_bp) / (region_bp / genome_bp)

    if isinstance(ratio, pd.Series):
        return ratio.replace([np.inf, -np.inf], 0).fillna(0)

    if ratio == float("inf") or ratio != ratio:  # check for NaN
        return 0.0

    return ratio


# ============================
# vectorized gene-level p-value
# ============================
def compute_gene_pvalues(
    line1_bp: pd.Series,
    region_bp: pd.Series,
    genome_bp: int,
    method: Literal["binomial","poisson"]="binomial"
) -> pd.Series:
    """
    Compute per-gene p-values for LINE1 enrichment using statistical tests.

    Parameters
    ----------
    line1_bp : pd.Series
        Observed number of LINE1-overlapping base pairs for each gene.
        Index should match the genes.
    region_bp : pd.Series
        Total base pair length of each gene or region.
        Must have the same index as line1_bp.
    genome_bp : int
        Total number of base pairs in the genome.
    method : str, default "binomial"
        Statistical method to compute p-values:
        - "binomial": use binomial test.
        - "poisson": use Poisson test.

    Returns
    -------
    pd.Series
        Per-gene p-values indicating the significance of LINE1 enrichment.
        NaN is returned for genes with region_bp <= 0.

    Notes
    -----
    Binomial method ("binomial"):
        - Assumes LINE1 elements are randomly distributed across the genome.
        - Null hypothesis (H0): Each base in the gene has independent probability 
          p = region_bp / genome_bp to be covered by LINE1.
        - Observed count k = line1_bp.
        - Test: right-tailed binomial test P(X >= k | n=region_bp, p).
        - k>n is corrected by setting n = max(n, k) to avoid invalid binomial input.

    Poisson method ("poisson"):
        - Assumes the total LINE1 coverage follows a Poisson process.
        - Null hypothesis (H0): LINE1 occurrences are uniformly distributed in the genome.
        - Expected LINE1 coverage for a gene: λ = (region_bp / genome_bp) * sum(line1_bp)
        - Observed count k = line1_bp.
        - Test: right-tailed Poisson test P(X >= k | λ).
        - Suitable when LINE1 counts are sparse and genome is large.

    Raises
    ------
    ValueError
        If `method` is not "binomial" or "poisson".
    """
    mask = region_bp > 0
    pvals = pd.Series(np.nan, index=line1_bp.index)

    if method == "binomial":
        prob = region_bp[mask] / genome_bp
        k = line1_bp[mask].round().astype(int)
        n = region_bp[mask].round().astype(int)
        # Fix cases where k > n
        n = n.combine(k, lambda n_val, k_val: max(n_val, k_val))
        pvals[mask] = binom.sf(k - 1, n, prob)
    elif method == "poisson":
        lam = region_bp[mask] / genome_bp * line1_bp.sum()
        pvals[mask] = poisson.sf(line1_bp[mask] - 1, lam)
    else:
        raise ValueError("method must be 'binomial' or 'poisson'")
    
    return pvals


# ============================
# core analysis
# ============================
def run_line1_enrichment(
    gene_gtf: str,
    te_gtf: str,
    upstream: int = 2000,
    out: str = "LINE1_enrichment_by_gene_type.tsv",
    per_gene_out: str = "LINE1_enrichment_per_gene.tsv",
    log_level: str = "INFO",
    genome_bp: int = 3_088_286_401,
    pvalue_method: Literal["binomial","poisson"]="binomial"
):
    setup_logger(log_level)

    # ---------- load ----------
    genes = load_gene_gtf_fast(gene_gtf)
    line1 = load_line1_gtf_fast(te_gtf)
    promoters = make_promoter(genes, upstream)
    total_line1_bp = line1.lengths().sum()
    logging.info("Total LINE1 bp: %d", total_line1_bp)

    # ---------- calculate overlaps ----------
    gene_l1 = genes.intersect(line1)
    promoter_l1 = promoters.intersect(line1)

    def calc_overlap_bp(pr_obj: pr.PyRanges, value_name: str) -> pd.DataFrame:
        if len(pr_obj)==0:
            return pd.DataFrame(columns=["gene_id", value_name])
        df = pr_obj.df.copy()
        df["overlap_bp"] = df.End - df.Start
        return df.groupby("gene_id", as_index=False)["overlap_bp"].sum().rename(columns={"overlap_bp": value_name})

    gene_l1_bp = calc_overlap_bp(gene_l1, "line1_gene_bp")
    promoter_l1_bp = calc_overlap_bp(promoter_l1, "line1_promoter_bp")

    genes_df = genes.df.copy()
    genes_df["gene_bp"] = genes_df.End - genes_df.Start
    promoters_df = promoters.df.copy()
    promoters_df["promoter_bp"] = promoters_df.End - promoters_df.Start
    promoter_bp = promoters_df.groupby("gene_id", as_index=False)["promoter_bp"].sum()

    # ---------- merge ----------
    df_gene = (
        genes_df[["gene_id","gene_type","gene_bp"]]
        .merge(promoter_bp, on="gene_id", how="left")
        .merge(gene_l1_bp, on="gene_id", how="left")
        .merge(promoter_l1_bp, on="gene_id", how="left")
        .fillna(0)
    )
    df_gene["total_bp"] = df_gene["gene_bp"] + df_gene["promoter_bp"]
    df_gene["line1_total_bp"] = df_gene["line1_gene_bp"] + df_gene["line1_promoter_bp"]

    # ---------- enrichment ----------
    df_gene["gene_enrichment"] = enrichment(df_gene["line1_gene_bp"], df_gene["gene_bp"], total_line1_bp, genome_bp)
    df_gene["promoter_enrichment"] = enrichment(df_gene["line1_promoter_bp"], df_gene["promoter_bp"], total_line1_bp, genome_bp)
    df_gene["total_enrichment"] = enrichment(df_gene["line1_total_bp"], df_gene["total_bp"], total_line1_bp, genome_bp)

    # ---------- p-values ----------
    df_gene["gene_pvalue"] = compute_gene_pvalues(df_gene["line1_gene_bp"], df_gene["gene_bp"], genome_bp, pvalue_method)
    df_gene["promoter_pvalue"] = compute_gene_pvalues(df_gene["line1_promoter_bp"], df_gene["promoter_bp"], genome_bp, pvalue_method)
    df_gene["total_pvalue"] = compute_gene_pvalues(df_gene["line1_total_bp"], df_gene["total_bp"], genome_bp, pvalue_method)

    # ---------- save per-gene ----------
    df_gene.sort_values("total_enrichment", ascending=False).to_csv(per_gene_out, sep="\t", index=False)
    logging.info("Per-gene result written to %s", per_gene_out)

    # ---------- gene_type aggregation ----------
    df_type = df_gene.groupby("gene_type", as_index=False).agg({
        "gene_bp":"sum",
        "promoter_bp":"sum",
        "total_bp":"sum",
        "line1_gene_bp":"sum",
        "line1_promoter_bp":"sum",
        "line1_total_bp":"sum"
    })
    df_type["gene_enrichment"] = enrichment(df_type["line1_gene_bp"], df_type["gene_bp"], total_line1_bp, genome_bp)
    df_type["promoter_enrichment"] = enrichment(df_type["line1_promoter_bp"], df_type["promoter_bp"], total_line1_bp, genome_bp)
    df_type["total_enrichment"] = enrichment(df_type["line1_total_bp"], df_type["total_bp"], total_line1_bp, genome_bp)
        # ---------- p-values ----------
    df_type["gene_pvalue"] = compute_gene_pvalues(df_type["line1_gene_bp"], df_type["gene_bp"], genome_bp, pvalue_method)
    df_type["promoter_pvalue"] = compute_gene_pvalues(df_type["line1_promoter_bp"], df_type["promoter_bp"], genome_bp, pvalue_method)
    df_type["total_pvalue"] = compute_gene_pvalues(df_type["line1_total_bp"], df_type["total_bp"], genome_bp, pvalue_method)
    df_type.sort_values("total_enrichment", ascending=False).to_csv(out, sep="\t", index=False)
    logging.info("Gene-type result written to %s", out)


# ============================
# CLI
# ============================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="LINE1 enrichment with per-gene p-values")
    p.add_argument("--gene-gtf", required=True)
    p.add_argument("--te-gtf", required=True)
    p.add_argument("--upstream", type=int, default=2000)
    p.add_argument("--out", default="LINE1_enrichment_by_gene_type.tsv")
    p.add_argument("--per-gene-out", default="LINE1_enrichment_per_gene.tsv")
    p.add_argument("--log-level", default="INFO")
    p.add_argument("--gsize", type=int, default=3_088_286_401)
    p.add_argument("--pvalue-method", default="binomial", choices=["binomial","poisson"])
    return p

if __name__=="__main__":
    # args = build_argparser().parse_args()
    # run_line1_enrichment(
    #     gene_gtf=args.gene_gtf,
    #     te_gtf=args.te_gtf,
    #     upstream=args.upstream,
    #     out=args.out,
    #     per_gene_out=args.per_gene_out,
    #     log_level=args.log_level,
    #     genome_bp=args.gsize,
    #     pvalue_method=args.pvalue_method
    # )


    # DEBUG example
    run_line1_enrichment(
        gene_gtf="/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf",
        te_gtf="/disk5/luosg/Reference/GENCODE/human/GRCh38/GRCh38_GENCODE_rmsk_TE.gtf",
        out = "/disk5/luosg/DIPseq20251215/output/LINE1/LINE1_enrichment_by_gene_type.tsv",
        per_gene_out= "/disk5/luosg/DIPseq20251215/output/LINE1/LINE1_enrichment_per_gene.tsv"
    )

# te_analysis.py
import re
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


class TEPeak_overlapAnalyzer:
    """
    Analyzer for computing overlap between TE loci and peak regions
    """

    # ===========================
    # class-level cache
    # ===========================
    _te_gtf_cache: dict = {}

    def __init__(
        self,
        te_gtf: str,
        peak_file: str,
        condition: str,
        sta_level_key: str = 'subfamily',
        method: str = 'fast'
    ):
        """
        te_gtf: TE annotation GTF file path
        peak_file: Peak narrowPeak file path
        condition: Condition label (e.g., "Control", "Treatment")
        sta_level_key: TE annotation level ('subfamily', 'family', or 'class')
        method: 'fast' or 'slow' for overlap computation
        """
        self.te_gtf = te_gtf
        self.peak_file = peak_file
        self.condition = condition
        self.sta_level_key = sta_level_key
        self.method = method

    # ===========================
    # Loaders
    # ===========================

    def load_narrowpeak(self) -> pr.PyRanges:
        cols = [
            "Chromosome", "Start", "End", "name", "score",
            "strand", "signalValue", "pValue", "qValue", "summit"
        ]
        df = pd.read_csv(self.peak_file, sep="\t", header=None, names=cols)
        return pr.PyRanges(df[["Chromosome", "Start", "End"]])

    def load_te_gtf(self) -> pr.PyRanges:
        """
        Load TE annotations from GTF file (with cache)
        Chromosome, Start, End are required by PyRanges
        """
        cache_key = (self.te_gtf, self.sta_level_key)

        # -------- cache hit --------
        if cache_key in self._te_gtf_cache:
            return self._te_gtf_cache[cache_key]

        # -------- cache miss --------
        gtf = pd.read_csv(
            self.te_gtf,
            sep="\t",
            comment="#",
            header=None,
            names=[
                "Chromosome", "Source", "Feature",
                "Start", "End", "Score", "Strand",
                "Frame", "Attribute"
            ]
        )

        gtf = gtf[gtf["Feature"] == "exon"].copy()

        if self.sta_level_key == 'family':
            gtf['family'] = gtf['Attribute'].str.extract(r'family_id\s*"(.*?)"')
        elif self.sta_level_key == 'class':
            gtf['class'] = gtf['Attribute'].str.extract(r'class_id\s*"(.*?)"')
        else:
            gtf['subfamily'] = gtf['Attribute'].str.extract(r'gene_id\s*"(.*?)"')

        gtf['locus_id'] = gtf['Attribute'].str.extract(r'transcript_id\s*"(.*?)"')
        gtf['locus_length'] = gtf["End"] - gtf["Start"]

        te_df = gtf[
            ["Chromosome", "Start", "End",
             "locus_id", self.sta_level_key, "locus_length"]
        ].dropna()

        te = pr.PyRanges(te_df)

        # store cache
        self._te_gtf_cache[cache_key] = te
        return te

    # ===========================
    # Computation (slow)
    # ===========================

    def compute_union_overlap_bp(self, df: pd.DataFrame) -> int:
        """
        Vectorized version.
        df: overlapping peak segments of ONE locus
        """
        if df.empty:
            return 0

        locus_start = df["Start"].iloc[0]
        locus_end = df["End"].iloc[0]

        s = np.maximum(df["Start_b"].to_numpy(), locus_start)
        e = np.minimum(df["End_b"].to_numpy(), locus_end)

        valid = s < e
        if not valid.any():
            return 0

        s = s[valid]
        e = e[valid]

        order = np.argsort(s)
        s = s[order]
        e = e[order]

        total = 0
        cur_start = s[0]
        cur_end = e[0]

        for i in range(1, len(s)):
            if s[i] <= cur_end:
                cur_end = max(cur_end, e[i])
            else:
                total += cur_end - cur_start
                cur_start = s[i]
                cur_end = e[i]

        total += cur_end - cur_start
        return int(total)

    def compute_locus_overlap(
        self,
        te: pr.PyRanges,
        peaks: pr.PyRanges,
        condition: str,
    ) -> pd.DataFrame:
        sta_level_key = self.sta_level_key

        overlap = te.join(peaks, how="left")
        df = overlap.df

        results = []
        for (locus_id, sta_level_value, locus_length), sub in df.groupby(
            ["locus_id", sta_level_key, "locus_length"]
        ):
            overlap_bp = self.compute_union_overlap_bp(sub)

            results.append({
                "locus_id": locus_id,
                sta_level_key: sta_level_value,
                "locus_length": locus_length,
                "overlap_bp": overlap_bp,
                "overlap_fraction": overlap_bp / locus_length,
                "condition": condition
            })

        return pd.DataFrame(results)

    # ===========================
    # Computation (fast)
    # ===========================

    def compute_locus_overlap_fast(
        self,
        te: pr.PyRanges,
        peaks: pr.PyRanges,
        condition: str,
    ) -> pd.DataFrame:
        """
        Fast overlap computation using PyRanges built-ins
        """
        sta_level_key = self.sta_level_key

        intersections = te.intersect(peaks).merge(
            by=["locus_id", sta_level_key]
        )

        intersections.overlap_bp = intersections.End - intersections.Start

        res_df = (
            intersections.df
            .groupby(["locus_id", sta_level_key], as_index=False)["overlap_bp"]
            .sum()
        )

        te_base = te.df[
            ["locus_id", sta_level_key, "Start", "End"]
        ].copy()
        te_base["locus_length"] = te_base["End"] - te_base["Start"]

        final_df = (
            pd.merge(
                te_base,
                res_df,
                on=["locus_id", sta_level_key],
                how="left"
            )
            .fillna(0)
        )

        final_df["overlap_fraction"] = (
            final_df["overlap_bp"] / final_df["locus_length"]
        )
        final_df["condition"] = condition

        return final_df

    # ===========================
    # Entry
    # ===========================

    def run_analysis(self) -> pd.DataFrame:
        te = self.load_te_gtf()
        peaks = self.load_narrowpeak()

        if self.method == "fast":
            return self.compute_locus_overlap_fast(
                te, peaks, self.condition
            )
        else:
            return self.compute_locus_overlap(
                te, peaks, self.condition
            )



def order_by_subfamily_mean_length(
    df: pd.DataFrame,
    sta_level_key: str = "subfamily",
    return_stats: bool = False,
    output_prefix: str = None
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute subfamily order by mean locus length and apply it to DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least:
        - sta_level_key
        - locus_length
    sta_level_key : str
        TE annotation level (default: 'subfamily')
    return_stats : bool
        If True, also return subfamily statistics DataFrame

    Returns
    -------
    pd.DataFrame
        Input DataFrame with sta_level_key converted to ordered categorical
    (optional) pd.DataFrame
        Subfamily statistics used for ordering
    """

    # 1. compute ordering
    stats = (
        df
        .groupby(sta_level_key, as_index=False)
        .agg(mean_length=("locus_length", "mean"))
        .sort_values("mean_length")
    )

    order = stats[sta_level_key].tolist()

    # 2. apply ordering
    out = df.copy()
    out[sta_level_key] = pd.Categorical(
        out[sta_level_key],
        categories=order,
        ordered=True
    )

    if output_prefix is not None:
        stats.to_csv(
            f"{output_prefix}_te_{sta_level_key}_stats.tsv",
            sep="\t",
            index=False
        )
        out.to_csv(
            f"{output_prefix}_te_{sta_level_key}_locus_overlap_fraction.tsv",
            sep="\t",
            index=False
        )
    if return_stats:
        return out, stats

    return out



def plot_control_vs_treatment(
        locus_overlap: pd.DataFrame,
        stats: pd.DataFrame,
        smooth_window: int = 5,
        figsize_scale: float = 0.45,
        sta_level_key: str = 'subfamily',
        outfig: str = None
    ) -> plt.Figure:
        """
        Upper: mean locus length (smoothed)
        Lower: overlap fraction boxplot
        """
        # sns.set_style("whitegrid")

        order = stats[sta_level_key].tolist()
        n = len(order)

        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1,
            figsize=(max(8, figsize_scale * n), 8),
            sharex=True,
            gridspec_kw={"height_ratios": [1, 3]}
        )

        # ---- Upper panel: mean length ----
        y = stats["mean_length"]
        y_smooth = y.rolling(
            window=smooth_window,
            center=True,
            min_periods=1
        ).mean()

        ax_top.plot(
            range(n),
            y_smooth,
            color="black",
            linewidth=2
        )
        ax_top.set_ylabel("Mean TE length (bp)")
        ax_top.set_title(f"TE {sta_level_key} length & peak overlap")

        # ---- Lower panel: boxplot ----
        sns.boxplot(
            data=locus_overlap,
            x=sta_level_key,
            y="overlap_fraction",
            hue="condition",
            showfliers=False,
            ax=ax_bottom
        )

        ax_bottom.set_ylabel("Overlap fraction per locus")
        ax_bottom.set_xlabel("")
        ax_bottom.tick_params(axis="x", rotation=90)

        ax_bottom.legend(
            title="Condition",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False
        )

        plt.tight_layout()
        if outfig is not None:
            plt.savefig(outfig, dpi=300)
        return fig
# ---------------------------
# Main
# ---------------------------

def main():
    gtf = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf"
    control_peak = "/disk5/luosg/DIPseq20251215/output/macs3/Alkbh1_wild_peaks.narrowPeak"
    treatment_peak = "/disk5/luosg/DIPseq20251215/output/macs3/Alkbh1_KO_peaks.narrowPeak"

    ControlAnalyzer = TEPeak_overlapAnalyzer(
        te_gtf=gtf,
        peak_file=control_peak,
        condition="Control",
        sta_level_key="family",
        method="fast"
    )
    df_control = ControlAnalyzer.run_analysis()

    TreatmentAnalyzer = TEPeak_overlapAnalyzer(
        te_gtf=gtf,
        peak_file=treatment_peak,
        condition="Treatment",
        sta_level_key="family",
        method="fast"
    )
    df_treatment = TreatmentAnalyzer.run_analysis()

    df_control = df_control[df_control["overlap_fraction"] > 0]
    df_treatment = df_treatment[df_treatment["overlap_fraction"] > 0]
    
    locus_overlap = pd.concat(
        [df_control, df_treatment],
        ignore_index=True
    )

    locus_overlap,stats = order_by_subfamily_mean_length(
        locus_overlap,
        sta_level_key="family",
        return_stats=True
    )


    plot_control_vs_treatment(
        locus_overlap, stats,
        sta_level_key="family",
        outfig="te_family_control_vs_treatment_overlap.png"
    )


if __name__ == "__main__":
    main()

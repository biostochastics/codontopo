"""CodonSafe: Reproducible analysis pipeline.

Generates all analysis outputs from raw data. Run with:
    python3.11 -m codon_topo.analysis.codonsafe.run_analyses

Outputs are written to output/codonsafe/.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import (
    binomtest,
    mannwhitneyu,
    spearmanr,
    wilcoxon,
)

from codon_topo.analysis.codonsafe.classify import (
    build_reference_context,
    classify_swap_event,
)
from codon_topo.analysis.codonsafe.loaders.genbank_utils import (
    CodonChange,
    compare_genomes,
    load_genbank,
)
from codon_topo.analysis.codonsafe.models import (
    CodonSwapEvent,
    StudyId,
    SwapModel,
)
from codon_topo.analysis.codonsafe.normalize import normalize_codon

OUT = Path("output/codonsafe")
OUT.mkdir(parents=True, exist_ok=True)

AA1_TO_3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}

SYN57_RECODE = {
    "TCG": "AGC",
    "TCT": "AGC",
    "TCC": "AGC",
    "TCA": "AGT",
    "GCA": "GCT",
    "GCG": "GCC",
    "TAG": "TAA",
}


def _classify_changes(
    changes: list[CodonChange],
    study_id: StudyId,
    ctx,
) -> pd.DataFrame:
    """Classify a list of CodonChange objects and return a DataFrame."""
    rows = []
    for ch in changes:
        aa3 = AA1_TO_3.get(ch.source_aa, ch.source_aa)
        tgt_aa3 = AA1_TO_3.get(ch.target_aa, ch.target_aa)
        src_rna = normalize_codon(ch.source_codon_dna)
        tgt_rna = normalize_codon(ch.target_codon_dna)

        event = CodonSwapEvent(
            study=study_id,
            event_id=f"dev_{ch.gene}_{ch.codon_index_in_cds}",
            source_codon=src_rna,
            target_codon=tgt_rna,
            table_id=11,
        )
        model = (
            SwapModel.SYNONYMOUS_REPLACEMENT
            if ch.is_synonymous
            else SwapModel.REASSIGN_SOURCE_CODON
        )
        topo = classify_swap_event(
            event=event, ctx=ctx, swap_model=model, metrics=("grantham",)
        )

        rows.append(
            {
                "gene": ch.gene or "unknown",
                "locus_tag": ch.locus_tag,
                "codon_idx": ch.codon_index_in_cds,
                "genome_pos": ch.genome_pos_nt,
                "design_codon": ch.source_codon_dna,
                "final_codon": ch.target_codon_dna,
                "design_aa": aa3,
                "final_aa": tgt_aa3,
                "is_synonymous": ch.is_synonymous,
                "is_reversion": ch.target_codon_dna in SYN57_RECODE,
                "hamming": topo.hamming,
                "crosses_boundary": topo.crosses_component_boundary_eps1,
                "changes_components": topo.changes_components_eps1,
                "local_mismatch_design": topo.local_mismatch_source["grantham"],
                "local_mismatch_final": topo.local_mismatch_target["grantham"],
                "delta_local": (
                    topo.local_mismatch_target["grantham"]
                    - topo.local_mismatch_source["grantham"]
                ),
            }
        )
    return pd.DataFrame(rows)


def run_syn57_deviation_analysis() -> dict:
    """Analysis 1: Syn57 S9 deviation directionality."""
    print("=== ANALYSIS 1: Syn57 S9 Deviation Directionality ===")

    BASE = Path("data/codonsafe/robertson2025_syn57/Data_files")

    ctx = build_reference_context(table_id=11, encoding_id=0)

    print("  Loading genomes...")
    design = load_genbank(BASE / "Data_S2_vS33A7_genome_design.gb")
    final = load_genbank(BASE / "Data_S8_Syn57_final_genome.gb")

    print("  Comparing design vs final...")
    all_diffs = compare_genomes(design, final)
    print(f"  Total CDS differences: {len(all_diffs)}")

    rdf = _classify_changes(all_diffs, StudyId.ROBERTSON_2025_SYN57, ctx)
    rdf.to_csv(OUT / "syn57_s9_deviations_classified.csv", index=False)

    # === Sensitivity analysis across filtering cutoffs ===
    print("  Running sensitivity analysis across cutoffs...")
    sensitivity = []
    for cutoff in [3, 5, 10, 15, 20, 50, None]:
        if cutoff is not None:
            gene_counts = rdf["gene"].value_counts()
            keep = gene_counts[gene_counts <= cutoff].index
            subset = rdf[rdf["gene"].isin(keep)]
            label = f"≤{cutoff}"
        else:
            subset = rdf
            label = "all"

        n = len(subset)
        n_better = int((subset["delta_local"] < 0).sum())
        n_worse = int((subset["delta_local"] > 0).sum())
        n_same = int((subset["delta_local"] == 0).sum())
        n_test = n_better + n_worse

        binom_p = float("nan")
        wilcox_p = float("nan")
        prop = float("nan")
        ci_lo = float("nan")
        ci_hi = float("nan")

        if n_test > 0:
            result = binomtest(n_better, n_test, 0.5, alternative="greater")
            binom_p = float(result.pvalue)
            ci = result.proportion_ci(confidence_level=0.95)
            prop = n_better / n_test
            ci_lo = float(ci.low)
            ci_hi = float(ci.high)

        nonzero = subset[subset["delta_local"] != 0]["delta_local"]
        if len(nonzero) >= 5:
            _, wilcox_p = wilcoxon(nonzero, alternative="less")
            wilcox_p = float(wilcox_p)

        sensitivity.append(
            {
                "cutoff": label,
                "n_deviations": n,
                "n_genes_removed": len(rdf) - n,
                "n_better": n_better,
                "n_worse": n_worse,
                "n_same": n_same,
                "proportion_better": round(prop, 4),
                "ci_95_lo": round(ci_lo, 4),
                "ci_95_hi": round(ci_hi, 4),
                "binomial_p": round(binom_p, 4),
                "wilcoxon_p": round(wilcox_p, 4),
                "mean_delta_local": round(float(subset["delta_local"].mean()), 2),
                "mean_hamming": round(float(subset["hamming"].mean()), 2),
            }
        )

    sens_df = pd.DataFrame(sensitivity)
    sens_df.to_csv(OUT / "syn57_s9_sensitivity_table.csv", index=False)
    print(f"  Sensitivity table saved ({len(sensitivity)} rows)")

    # Primary result at cutoff=10
    genuine = rdf[rdf["gene"].map(rdf["gene"].value_counts()) <= 10].copy()
    genuine.to_csv(OUT / "syn57_s9_genuine_deviations.csv", index=False)

    # R-ready figure data
    genuine[
        [
            "gene",
            "design_codon",
            "final_codon",
            "design_aa",
            "final_aa",
            "is_synonymous",
            "is_reversion",
            "hamming",
            "delta_local",
            "local_mismatch_design",
            "local_mismatch_final",
        ]
    ].to_csv(OUT / "syn57_s9_for_ggplot.csv", index=False)

    # Compute summary stats for the primary cutoff
    n_better = int((genuine["delta_local"] < 0).sum())
    n_worse = int((genuine["delta_local"] > 0).sum())
    n_same = int((genuine["delta_local"] == 0).sum())
    n_test = n_better + n_worse
    binom_result = binomtest(n_better, n_test, 0.5, alternative="greater")
    ci = binom_result.proportion_ci(confidence_level=0.95)

    summary = {
        "total_deviations": len(rdf),
        "genuine_deviations": len(genuine),
        "genes_removed": int(len(rdf) - len(genuine)),
        "n_better": n_better,
        "n_worse": n_worse,
        "n_same": n_same,
        "proportion_better": round(n_better / max(n_test, 1), 4),
        "ci_95_lo": round(float(ci.low), 4),
        "ci_95_hi": round(float(ci.high), 4),
        "binomial_p": round(float(binom_result.pvalue), 4),
        "mean_delta_local": round(float(genuine["delta_local"].mean()), 2),
        "median_delta_local": round(float(genuine["delta_local"].median()), 2),
        "mean_hamming": round(float(genuine["hamming"].mean()), 2),
        "n_reversions": int(genuine["is_reversion"].sum()),
    }
    pd.DataFrame([summary]).to_csv(OUT / "syn57_s9_genuine_summary.csv", index=False)

    print(f"  Result: {n_better} better, {n_worse} worse, {n_same} same")
    print(
        f"  Binomial p={binom_result.pvalue:.4f}, 95% CI [{ci.low:.3f}, {ci.high:.3f}]"
    )

    return summary


def run_ostrov_case_control() -> dict:
    """Analysis 2: Ostrov segment case-control."""
    print("\n=== ANALYSIS 2: Ostrov Segment Case-Control ===")

    BASE = Path("data/codonsafe/ostrov2016")

    df = pd.read_excel(BASE / "table_s4.xlsx", header=None)
    headers = df.iloc[1].tolist()
    data = df.iloc[3:].copy()
    data.columns = [
        str(c).strip() if pd.notna(c) else f"col_{i}" for i, c in enumerate(headers)
    ]

    # Column aliases
    seg = "Seg"
    n_rec = "Total Number Of Recoded Codons"
    n_ess_rec = "Number Of Recoded Codons in Essential Genes"
    n_ess = "Number Of Essential Genes"
    dt_del = "Relative Doubling Time After Chromosomal Deletion (mean)"
    dt_int = "Relative Doubling Time After Segment Integration (mean)"
    reversions = (
        "Number of Codon Reversions found in E. coli Segment After Chromosomal Deletion"
    )
    comments = "Comments"

    for col in [n_rec, n_ess_rec, n_ess]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    for col in [dt_del, dt_int]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    data[reversions] = pd.to_numeric(data[reversions], errors="coerce")

    data["has_lethal"] = (
        data[comments].fillna("").str.contains("Lethal|lethal|dead", case=False)
    )
    data["has_issue"] = data["has_lethal"]
    data["has_reversions"] = data[reversions].fillna(0) > 0
    data["status"] = np.where(data["has_issue"], "Problem", "Normal")

    # Fitness correlation
    fitness = data.dropna(subset=[dt_del])
    n_fitness = len(fitness)
    rho_ess, p_ess = spearmanr(fitness[n_ess_rec], fitness[dt_del])
    rho_tot, p_tot = spearmanr(fitness[n_rec], fitness[dt_del])

    # Bonferroni correction for 3 tests
    p_bonf_ess = min(float(p_ess) * 3, 1.0)
    p_bonf_tot = min(float(p_tot) * 3, 1.0)

    # Case-control
    problem = data[data["has_issue"]]
    normal = data[~data["has_issue"]]
    n_problem = len(problem)
    n_normal = len(normal)

    mw_stat, mw_p = mannwhitneyu(
        problem[n_ess_rec].dropna(),
        normal[n_ess_rec].dropna(),
        alternative="two-sided",
    )
    mw_p_bonf = min(float(mw_p) * 3, 1.0)

    # Save R-ready data
    fig_cols = [
        seg,
        n_rec,
        n_ess_rec,
        n_ess,
        dt_del,
        dt_int,
        reversions,
        "has_issue",
        "has_lethal",
        "has_reversions",
        "status",
    ]
    fig_data = data[[c for c in fig_cols if c in data.columns]].copy()
    fig_data.columns = [
        "segment",
        "n_recoded",
        "n_essential_recoded",
        "n_essential",
        "doubling_time_deletion",
        "doubling_time_integration",
        "n_reversions",
        "has_issue",
        "has_lethal",
        "has_reversions",
        "status",
    ]
    fig_data.to_csv(OUT / "ostrov_segments_for_ggplot.csv", index=False)

    # Save full analysis table
    data.to_csv(OUT / "ostrov_segment_analysis.csv", index=False)

    summary = {
        "n_segments": len(data),
        "n_with_fitness": n_fitness,
        "n_problem": n_problem,
        "n_normal": n_normal,
        "rho_essential_recoded_vs_dt": round(float(rho_ess), 4),
        "p_essential_recoded_vs_dt": round(float(p_ess), 4),
        "p_bonferroni_essential": round(p_bonf_ess, 4),
        "rho_total_recoded_vs_dt": round(float(rho_tot), 4),
        "p_total_recoded_vs_dt": round(float(p_tot), 4),
        "p_bonferroni_total": round(p_bonf_tot, 4),
        "mw_essential_problem_mean": round(float(problem[n_ess_rec].mean()), 1),
        "mw_essential_normal_mean": round(float(normal[n_ess_rec].mean()), 1),
        "mw_p": round(float(mw_p), 4),
        "mw_p_bonferroni": round(mw_p_bonf, 4),
        "n_with_reversions": int(data["has_reversions"].sum()),
    }
    pd.DataFrame([summary]).to_csv(OUT / "ostrov_summary.csv", index=False)

    print(f"  N={n_fitness} segments with fitness data")
    print(
        f"  rho(essential_recoded, DT) = {rho_ess:.4f}, "
        f"p = {p_ess:.4f} (Bonf. = {p_bonf_ess:.4f})"
    )
    print(
        f"  Problem vs Normal: {problem[n_ess_rec].mean():.1f} vs "
        f"{normal[n_ess_rec].mean():.1f}, p = {mw_p:.4f} "
        f"(Bonf. = {mw_p_bonf:.4f})"
    )

    return summary


def main():
    """Run all CodonSafe analyses and save outputs."""
    print("CodonSafe Analysis Pipeline")
    print("=" * 60)

    syn57 = run_syn57_deviation_analysis()
    ostrov = run_ostrov_case_control()

    # Save combined statistics for R figure generation
    combined = {"syn57": syn57, "ostrov": ostrov}
    with open(OUT / "analysis_stats.json", "w") as f:
        json.dump(combined, f, indent=2, default=str)

    print(f"\n{'=' * 60}")
    print("All analyses complete. Outputs in output/codonsafe/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()

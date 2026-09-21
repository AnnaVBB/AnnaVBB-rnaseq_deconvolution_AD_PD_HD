configfile: "config/config.yaml"

DATASETS = config.get(
    "core_datasets",
    ["GSE53697", "GSE64810", "GSE68719"]
)


# ============================================================
# ALL
# ============================================================

rule all:
    input:
        # Auditoria
        expand(
            "results/qc/{dataset}_audit_report.txt",
            dataset=DATASETS
        ),

        # Metadata
        expand(
            "metadata/{dataset}_samples.csv",
            dataset=DATASETS
        ),

        # DEG
        expand(
            "results/01_core_original/{dataset}/{dataset}_deg_results.csv",
            dataset=DATASETS
        ),

        # PCA
        expand(
            "results/01_core_original/{dataset}/{dataset}_pca.png",
            dataset=DATASETS
        ),

        # Volcano
        expand(
            "results/01_core_original/{dataset}/{dataset}_volcano.png",
            dataset=DATASETS
        ),

        # Heatmap
        expand(
            "results/01_core_original/{dataset}/{dataset}_heatmap.png",
            dataset=DATASETS
        ),

        # GO enrichment
        expand(
            "results/01_core_original/{dataset}/{dataset}_pathways_GO_BP.csv",
            dataset=DATASETS
        ),

        # Intersection
        "results/01_core_original/intersection/venn_degs.png",
        "results/01_core_original/intersection/common_genes.csv",

        # Report
        "reports/01_core_original.html"


# ============================================================
# 00 — AUDIT
# ============================================================

rule audit_dataset:
    input:
        counts="data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        annot="data/raw/Human.GRCh38.p13.annot.tsv.gz",
        soft="data/raw/{dataset}_family.soft.gz"

    output:
        report="results/qc/{dataset}_audit_report.txt"

    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/00_auditoria.R \
        {input.counts} \
        {input.annot} \
        {input.soft} \
        {output.report}
        """


# ============================================================
# 01 — METADATA
# ============================================================

rule parse_metadata:
    input:
        soft="data/raw/{dataset}_family.soft.gz",
        counts="data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        audit="results/qc/{dataset}_audit_report.txt"

    output:
        meta="metadata/{dataset}_samples.csv"

    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/01_parse_metadata.R \
        {input.soft} \
        {input.counts} \
        {output.meta}
        """


# ============================================================
# 02b — DEG + PCA + VOLCANO + HEATMAP
# ============================================================

rule run_exploratory_analysis_original:
    input:
        counts="data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        meta="metadata/{dataset}_samples.csv"

    output:
        deg="results/01_core_original/{dataset}/{dataset}_deg_results.csv",
        pca="results/01_core_original/{dataset}/{dataset}_pca.png",
        volcano="results/01_core_original/{dataset}/{dataset}_volcano.png",
        heatmap="results/01_core_original/{dataset}/{dataset}_heatmap.png"

    shell:
        """
        Rscript scripts/02b_exploratory_deg.R \
        {input.counts} \
        {input.meta} \
        results/01_core_original/{wildcards.dataset}/{wildcards.dataset}
        """


# ============================================================
# 02c — GO ENRICHMENT
# ============================================================

rule run_pathway_enrichment:
    input:
        deg="results/01_core_original/{dataset}/{dataset}_deg_results.csv"

    output:
        go="results/01_core_original/{dataset}/{dataset}_pathways_GO_BP.csv"

    shell:
        """
        Rscript scripts/02c_pathway_enrichment.R \
        {input.deg} \
        {output.go}
        """


# ============================================================
# 02d — INTERSECTION
# ============================================================

rule run_intersection_analysis:
    input:
        degs=expand(
            "results/01_core_original/{dataset}/{dataset}_deg_results.csv",
            dataset=DATASETS
        )

    output:
        venn="results/01_core_original/intersection/venn_degs.png",
        common="results/01_core_original/intersection/common_genes.csv"

    shell:
        """
        Rscript scripts/02d_intersection_analysis.R
        """


# ============================================================
# REPORT
# ============================================================

rule generate_html_report_phase01:
    input:
        degs=expand(
            "results/01_core_original/{dataset}/{dataset}_deg_results.csv",
            dataset=DATASETS
        ),

        pcas=expand(
            "results/01_core_original/{dataset}/{dataset}_pca.png",
            dataset=DATASETS
        ),

        volcanos=expand(
            "results/01_core_original/{dataset}/{dataset}_volcano.png",
            dataset=DATASETS
        ),

        heatmaps=expand(
            "results/01_core_original/{dataset}/{dataset}_heatmap.png",
            dataset=DATASETS
        ),

        enrichment=expand(
            "results/01_core_original/{dataset}/{dataset}_pathways_GO_BP.csv",
            dataset=DATASETS
        ),

        venn="results/01_core_original/intersection/venn_degs.png",
        common="results/01_core_original/intersection/common_genes.csv",
        rmd="scripts/reports/01_core_original.Rmd"

    output:
        "reports/01_core_original.html"

    shell:
        """
        Rscript -e "rmarkdown::render(
          '{input.rmd}',
          output_dir = 'reports'
        )"
        """
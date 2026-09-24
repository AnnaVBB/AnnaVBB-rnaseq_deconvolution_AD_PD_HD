configfile: "config/config.yaml"

DATASETS = config.get(
    "core_datasets",
    ["GSE53697", "GSE64810", "GSE68719"]
)

# ============================================================
# ARQUIVOS ORIGINAIS DOS AUTORES — FASE 01
# ============================================================

AUTHOR_EXPRESSION = {
    "GSE53697": "data/raw/GSE53697/GSE53697_RNAseq_AD.txt.gz",
    "GSE64810": "data/raw/GSE64810/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz",
    "GSE68719": "data/raw/GSE68719/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz"
}

SOFT_FILES = {
    "GSE53697": "data/raw/GSE53697/GSE53697_family.soft.gz",
    "GSE64810": "data/raw/GSE64810/GSE64810_family.soft.gz",
    "GSE68719": "data/raw/GSE68719/GSE68719_family.soft.gz"
}

ANNOT_FILES = {
    "GSE53697": "data/raw/GSE53697/Human.GRCh38.p13.annot.tsv.gz",
    "GSE64810": "data/raw/GSE64810/Human.GRCh38.p13.annot.tsv.gz",
    "GSE68719": "data/raw/GSE68719/Human.GRCh38.p13.annot.tsv.gz"
}

# ============================================================
# FUNÇÕES AUXILIARES
# ============================================================
def author_expression(wildcards):
    return AUTHOR_EXPRESSION[wildcards.dataset]

def soft_file(wildcards):
    return SOFT_FILES[wildcards.dataset]

def annot_file(wildcards):
    return ANNOT_FILES[wildcards.dataset]

# ALL
rule all:
    input:
        # Auditoria
        expand(
            "results/qc/{dataset}_audit_report.txt",
            dataset=DATASETS
        ),

        # Metadata final
        expand(
            "metadata/{dataset}_metadata.csv",
            dataset=DATASETS
        ),

        # Matrizes preparadas
        expand(
            "data/processed/{dataset}_expression_prepared.csv",
            dataset=DATASETS
        ),

        # DEG
        expand(
            "results/01_core_original/{dataset}/deg/deg_all.csv",
            dataset=DATASETS
        ),

        expand(
            "results/01_core_original/{dataset}/deg/deg_summary.csv",
            dataset=DATASETS
        ),

        # PCA
        expand(
            "results/01_core_original/{dataset}/plots/pca.png",
            dataset=DATASETS
        ),

        # Volcano principal
        expand(
            "results/01_core_original/{dataset}/plots/volcano_primary.png",
            dataset=DATASETS
        ),

        # Volcano artigo
        expand(
            "results/01_core_original/{dataset}/plots/volcano_article.png",
            dataset=DATASETS
        ),

        # Heatmap variável
        expand(
            "results/01_core_original/{dataset}/plots/heatmap_top50_variable.png",
            dataset=DATASETS
        ),

        # Enriquecimento
        expand(
            "results/01_core_original/{dataset}/enrichment/enrichment_summary.csv",
            dataset=DATASETS
        ),

        # Integração
        "results/01_core_original/intersection/intersection_summary.csv",
        "results/01_core_original/intersection/integrated_AD_HD_PD.csv",
        "results/01_core_original/intersection/common_testable_universe.csv",
        # Plots da integração
        "results/01_core_original/intersection/plots/primary_deg_patterns.png",
        "results/01_core_original/intersection/plots/logFC_HD_vs_PD.png",
        "results/01_core_original/intersection/plots/venn_primary_AD_HD_PD.png",
        "results/01_core_original/intersection/plots/venn_article_AD_HD_PD.png",
        
        # Relatório final da Fase 01
        "reports/01_core_original.html"

# ============================================================
# 00 — AUDITORIA
# ============================================================

rule audit_dataset:
    input:
        expression=author_expression,
        soft=soft_file,
        annot=annot_file

    output:
        report="results/qc/{dataset}_audit_report.txt"

    shell:
        """
        mkdir -p results/qc

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/00_auditoria.R \
            {wildcards.dataset} \
            {input.expression} \
            {input.soft} \
            {input.annot} \
            {output.report}
        """


# ============================================================
# 01 — METADATA
# ============================================================

rule parse_metadata:
    input:
        soft=soft_file,
        expression=author_expression,
        audit="results/qc/{dataset}_audit_report.txt"

    output:
        meta="metadata/{dataset}_metadata.csv"

    shell:
        """
        mkdir -p metadata

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/01_parse_metadata.R \
            {wildcards.dataset} \
            {input.soft} \
            {input.expression} \
            {output.meta}
        """


# ============================================================
# 02a — PREPARAÇÃO AD
# ============================================================

rule prepare_AD:
    input:
        expression=AUTHOR_EXPRESSION["GSE53697"],
        annot=ANNOT_FILES["GSE53697"],
        meta="metadata/GSE53697_metadata.csv"

    output:
        prepared="data/processed/GSE53697_expression_prepared.csv"

    shell:
        """
        mkdir -p data/processed

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/02a_prepare_AD.R \
            {input.expression} \
            {input.annot} \
            {input.meta} \
            {output.prepared}
        """


# ============================================================
# 02b — PREPARAÇÃO HD
# ============================================================

rule prepare_HD:
    input:
        expression=AUTHOR_EXPRESSION["GSE64810"],
        annot=ANNOT_FILES["GSE64810"],
        meta="metadata/GSE64810_metadata.csv"

    output:
        prepared="data/processed/GSE64810_expression_prepared.csv"

    shell:
        """
        mkdir -p data/processed

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/02b_prepare_HD_PD.R \
            GSE64810 \
            {input.expression} \
            {input.annot} \
            {input.meta} \
            {output.prepared}
        """


# ============================================================
# 02b — PREPARAÇÃO PD
# ============================================================

rule prepare_PD:
    input:
        expression=AUTHOR_EXPRESSION["GSE68719"],
        annot=ANNOT_FILES["GSE68719"],
        meta="metadata/GSE68719_metadata.csv"

    output:
        prepared="data/processed/GSE68719_expression_prepared.csv"

    shell:
        """
        mkdir -p data/processed

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/02b_prepare_HD_PD.R \
            GSE68719 \
            {input.expression} \
            {input.annot} \
            {input.meta} \
            {output.prepared}
        """

# ============================================================
# 03 — QC + PCA + DEG
# ============================================================
rule exploratory_deg:
    input:
        expression="data/processed/{dataset}_expression_prepared.csv",
        meta="metadata/{dataset}_metadata.csv"

    output:
        deg_all="results/01_core_original/{dataset}/deg/deg_all.csv",
        deg_primary="results/01_core_original/{dataset}/deg/deg_significant_primary.csv",
        deg_article="results/01_core_original/{dataset}/deg/deg_significant_article.csv",
        deg_summary="results/01_core_original/{dataset}/deg/deg_summary.csv",

        pca_coordinates="results/01_core_original/{dataset}/plots/pca_coordinates.csv",
        pca="results/01_core_original/{dataset}/plots/pca.png",
        distribution="results/01_core_original/{dataset}/plots/expression_distribution.png",
        volcano_primary="results/01_core_original/{dataset}/plots/volcano_primary.png",
        volcano_article="results/01_core_original/{dataset}/plots/volcano_article.png",
        heatmap_variable="results/01_core_original/{dataset}/plots/heatmap_top50_variable.png"

    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/03_exploratory_deg.R \
            {wildcards.dataset} \
            {input.expression} \
            {input.meta} \
            results/01_core_original/{wildcards.dataset}
        """


# ============================================================
# 04 — GO BP + KEGG
# ============================================================

rule pathway_enrichment:
    input:
        deg="results/01_core_original/{dataset}/deg/deg_all.csv"

    output:
        summary="results/01_core_original/{dataset}/enrichment/enrichment_summary.csv"

    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/04_pathway_enrichment.R \
            {wildcards.dataset} \
            {input.deg} \
            results/01_core_original/{wildcards.dataset}/enrichment
        """


# ============================================================
# 05 — INTEGRAÇÃO TRANSDIAGNÓSTICA
# ============================================================

rule intersection_analysis:
    input:
        ad="results/01_core_original/GSE53697/deg/deg_all.csv",
        hd="results/01_core_original/GSE64810/deg/deg_all.csv",
        pd="results/01_core_original/GSE68719/deg/deg_all.csv"

    output:
        integrated="results/01_core_original/intersection/integrated_AD_HD_PD.csv",
        common="results/01_core_original/intersection/common_testable_universe.csv",
        summary="results/01_core_original/intersection/intersection_summary.csv",

        primary_patterns="results/01_core_original/intersection/primary_pattern_summary.csv",
        article_patterns="results/01_core_original/intersection/article_pattern_summary.csv",
        directional_patterns="results/01_core_original/intersection/directional_pattern_summary.csv",

        shared_hd_pd="results/01_core_original/intersection/shared_primary_HD_PD.csv",
        shared_same="results/01_core_original/intersection/shared_primary_HD_PD_same_direction.csv",
        shared_opposite="results/01_core_original/intersection/shared_primary_HD_PD_opposite_direction.csv",

        concordant="results/01_core_original/intersection/directionally_concordant_AD_HD_PD.csv",

        correlations="results/01_core_original/intersection/logFC_spearman_correlations.csv",

        patterns_plot="results/01_core_original/intersection/plots/primary_deg_patterns.png",
        scatter_plot="results/01_core_original/intersection/plots/logFC_HD_vs_PD.png",
        venn_primary="results/01_core_original/intersection/plots/venn_primary_AD_HD_PD.png",
        venn_article="results/01_core_original/intersection/plots/venn_article_AD_HD_PD.png"

    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false \
        Rscript scripts/05_intersection_analysis.R \
            {input.ad} \
            {input.hd} \
            {input.pd} \
            results/01_core_original/intersection
        """

# ============================================================
# 06 — RELATÓRIO DA FASE 01
# ============================================================
rule report_core_original:
    input:
        rmd="scripts/reports/01_core_original.Rmd",

        # Resultados DEG
        deg_ad="results/01_core_original/GSE53697/deg/deg_all.csv",
        deg_hd="results/01_core_original/GSE64810/deg/deg_all.csv",
        deg_pd="results/01_core_original/GSE68719/deg/deg_all.csv",

        # PCA
        pca_ad="results/01_core_original/GSE53697/plots/pca.png",
        pca_hd="results/01_core_original/GSE64810/plots/pca.png",
        pca_pd="results/01_core_original/GSE68719/plots/pca.png",

        # Volcano principal
        volcano_primary_ad="results/01_core_original/GSE53697/plots/volcano_primary.png",
        volcano_primary_hd="results/01_core_original/GSE64810/plots/volcano_primary.png",
        volcano_primary_pd="results/01_core_original/GSE68719/plots/volcano_primary.png",

        # Volcano comparativo
        volcano_article_ad="results/01_core_original/GSE53697/plots/volcano_article.png",
        volcano_article_hd="results/01_core_original/GSE64810/plots/volcano_article.png",
        volcano_article_pd="results/01_core_original/GSE68719/plots/volcano_article.png",

        # Heatmaps principais disponíveis
        heatmap_hd="results/01_core_original/GSE64810/plots/heatmap_top50_primary.png",
        heatmap_pd="results/01_core_original/GSE68719/plots/heatmap_top50_primary.png",

        # Enriquecimento
        enrichment_hd="results/01_core_original/GSE64810/enrichment/enrichment_summary.csv",
        enrichment_pd="results/01_core_original/GSE68719/enrichment/enrichment_summary.csv",

        # Integração
        intersection="results/01_core_original/intersection/intersection_summary.csv",
        integrated="results/01_core_original/intersection/integrated_AD_HD_PD.csv",
        common="results/01_core_original/intersection/common_testable_universe.csv",
        venn_primary="results/01_core_original/intersection/plots/venn_primary_AD_HD_PD.png",
        venn_article="results/01_core_original/intersection/plots/venn_article_AD_HD_PD.png",
        scatter="results/01_core_original/intersection/plots/logFC_HD_vs_PD.png"

    output:
        html="reports/01_core_original.html"

    shell:
        """
        mkdir -p reports

        RENV_CONFIG_AUTO_LOAD=false \
        Rscript -e "rmarkdown::render(
            'scripts/reports/01_core_original.Rmd',
            output_file = '01_core_original.html',
            output_dir = 'reports',
            knit_root_dir = getwd()
        )"
        """
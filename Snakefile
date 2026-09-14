configfile: "config/config.yaml"

DATASETS = config.get("core_datasets", ["GSE53697", "GSE64810", "GSE68719"])

rule all:
    input:
        expand("results/qc/{dataset}_audit_report.txt", dataset=DATASETS),
        expand("metadata/{dataset}_samples.csv", dataset=DATASETS),
        expand("results/01_core_original/{dataset}/{dataset}_deg_results.csv", dataset=DATASETS),
        "reports/01_core_original.html"

rule audit_dataset:
    input:
        counts = "data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        annot  = "data/raw/Human.GRCh38.p13.annot.tsv.gz",
        soft   = "data/raw/{dataset}_family.soft.gz"
    output:
        report = "results/qc/{dataset}_audit_report.txt"
    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false Rscript scripts/00_auditoria.R {input.counts} {input.annot} {input.soft} {output.report}
        """

rule parse_metadata:
    input:
        soft   = "data/raw/{dataset}_family.soft.gz",
        counts = "data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        audit  = "results/qc/{dataset}_audit_report.txt"
    output:
        meta   = "metadata/{dataset}_samples.csv"
    shell:
        """
        RENV_CONFIG_AUTO_LOAD=false Rscript scripts/01_parse_metadata.R {input.soft} {input.counts} {output.meta}
        """

rule run_exploratory_analysis_original:
    input:
        counts = "data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        meta   = "metadata/{dataset}_samples.csv"
    output:
        deg     = "results/01_core_original/{dataset}/{dataset}_deg_results.csv",
        pca     = "results/01_core_original/{dataset}/{dataset}_pca.png",
        volcano = "results/01_core_original/{dataset}/{dataset}_volcano.png",
        heatmap = "results/01_core_original/{dataset}/{dataset}_heatmap.png"
    shell:
        """
        Rscript scripts/02b_exploratory_deg.R \
            {input.counts} {input.meta} results/01_core_original/{wildcards.dataset}/{wildcards.dataset}
        """

# Regra explícita para a Análise de Interseção (Venn / UpSet)
rule run_intersection_analysis:
    input:
        degs = expand("results/01_core_original/{dataset}/{dataset}_deg_results.csv", dataset=DATASETS)
    output:
        venn = "results/01_core_original/intersection/venn_degs.png"
    shell:
        """
        Rscript scripts/03_intersection_analysis.R
        """

#Regra final para renderizar o relatório HTML completo
rule generate_html_report_phase01:
    input:
        "results/01_core_original/GSE53697/GSE53697_deg_results.csv",
        "results/01_core_original/GSE64810/GSE64810_deg_results.csv",
        "results/01_core_original/GSE68719/GSE68719_deg_results.csv",
        rmd = "scripts/reports/01_core_original.Rmd"
    output:
        "reports/01_core_original.html"
    shell:
        """
        Rscript -e "rmarkdown::render('{input.rmd}', output_dir = 'reports')"
        """
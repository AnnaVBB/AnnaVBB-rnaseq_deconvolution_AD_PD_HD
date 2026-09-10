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
        volcano = "results/01_core_original/{dataset}/{dataset}_volcano.pdf"
    shell:
        """
        Rscript scripts/02b_exploratory_deg.R \
            {input.counts} {input.meta} results/01_core_original/{wildcards.dataset}/{wildcards.dataset}
        """

rule generate_html_report_phase01:
    input:
        degs = expand("results/01_core_original/{dataset}/{dataset}_deg_results.csv", dataset=DATASETS),
        rmd  = "scripts/reports/01_core_original.Rmd"
    output:
        html = "reports/01_core_original.html"
    shell:
        """
        Rscript -e "rmarkdown::render('scripts/reports/01_core_original.Rmd')"
        mv scripts/reports/01_core_original.html reports/01_core_original.html
        """
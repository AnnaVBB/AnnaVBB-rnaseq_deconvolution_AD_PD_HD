configfile: "config/config.yaml"

DATASETS = config.get("core_datasets", ["GSE53697", "GSE64810", "GSE68719"])

rule all:
    input:
        expand("results/qc/{dataset}_audit_report.txt", dataset=DATASETS),
        expand("metadata/{dataset}_samples.csv", dataset=DATASETS)

rule audit_dataset:
    input:
        counts = "data/raw/{dataset}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        annot  = "data/raw/Human.GRCh38.p13.annot.tsv.gz",
        soft   = "data/raw/{dataset}_family.soft.gz"
    output:
        report = "results/qc/{dataset}_audit_report.txt"
    shell:
        """
        Rscript scripts/00_auditoria.R {input.counts} {input.annot} {input.soft} {output.report}
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
        Rscript scripts/01_parse_metadata.R {input.soft} {input.counts} {output.meta}
        """

# Parâmetros documentados nas nossas etapas 06–08 e teste de replicação 08e

## MuSiC2 (06)
- Referência Mathys 2019 pathology-free; 24 doadores; 35.140 núcleos; oito tipos Ex, In, Ast, Oli, Opc, Mic, End, Per. Sensibilidade ROSMAP sem três doadores sobrepostos: 21 doadores.
- Entrada bulk e núcleos: counts brutos, restritos à união top300 por tipo celular (2.222 genes ROSMAP e 2.241 genes nas outras três coortes, depois dos filtros).
- `MuSiC::music2_prop` com `markers=NULL` porque os dados já foram filtrados; `select.ct = c("Ex", "In", "Ast", "Oli", "Opc", "Mic", "End", "Per")`, `cell_size=NULL`, `expr_low=20`, `prop_r=.10`, `eps_c=.05`, `eps_r=.01`, `n_resample=20`, `sample_prop=.50`, `cutoff_expr=.05`, `cutoff_c=.05`, `cutoff_r=.01`, `maxiter=200`, `ct.cov=FALSE`, `centered=FALSE`, `normalize=FALSE`.
- Arquivos: `scripts/phase02_deconvolution/06_run_MuSiC2.R`, `results/deconvolution/MuSiC2/06_run_summary.tsv`, `parameters.tsv` por ajuste em `results/deconvolution/MuSiC2/<dataset>/<reference_version>/<marker_set>/` (conferir o manifest do 06).

## BisqueRNA (07)
- Mesmo Mathys main24, 35.140 núcleos, oito tipos e top300 (2222/2241 genes após filtros); entrada `ExpressionSet` raw counts.
- `BisqueRNA::ReferenceBasedDecomposition` com `markers=genes_top300`, `cell.types="cellType"`, `subject.names="SubjectName"`, `use.overlap=FALSE`, `old.cpm=FALSE`; biblioteca completa preservada por complemento técnico que é excluído da regressão.
- Saídas: `results/deconvolution/BisqueRNA/07_run_summary.tsv`, `07_api_audit.tsv`, `07_celltype_variability_all.tsv`; parâmetros detalhados por ajuste em `results/deconvolution/BisqueRNA/<dataset>/main24/top300/parameters.tsv` (conferir nome no arquivo 07).

## dtangle original (08)
- Mathys main24, oito tipos, **181 pseudobulks** doador × tipo; genes elegíveis restritos à união top300 da etapa 04: **2.222 ROSMAP / 2.241 demais**.
- Genes que apareciam como marcadores de tipos diferentes foram resolvidos com prioridade por menor rank e depois escore / ordem fixa de tipos; 2.400 atribuições originais, 2.242 genes únicos antes dos filtros.
- Referência: somar counts por doador/tipo em cada pseudobulk e dividir pela biblioteca completa desse pseudobulk antes de log2(CPM+1).
- Bulk: contar por GeneID e usar `library_size_full_input` da etapa 05 para CPM; transformar `log2(CPM+1)`.
- `dtangle::dtangle(Y=..., references=..., pure_samples=..., markers=lista_explicita, data_type="rna-seq")`, sem `n_markers=100`, sem normalização quantílica.
- `gamma` retornado pelo pacote no run: **0.943390169106717**. Resumo 08: `n_markers` corresponde ao total de marcadores **2.222/2.241 distribuídos entre os oito tipos**, não a 2.222 marcadores por tipo.
- Arquivos: `scripts/phase02_deconvolution/08_run_dtangle.R`, `results/deconvolution/dtangle/08_run_summary.tsv`, `08_api_audit.tsv`, `08_marker_overlap_resolution.tsv`, `08_markers_unique_assignment.tsv`.

## dtangle auditoria (08c) e validação (08d)
- 08c: 181 pseudobulks versus oito perfis por tipo, com/sem quantílica; preservou o conjunto top300 e os marcadores explicitamente atribuídos. `pooled` significa SOMAR counts de todos os núcleos de um tipo antes do CPM, **não** necessariamente fazer a média dos CPM individuais.
- 08d: validação LODO, quatro doadores retidos, nove cenários por doador; marcador top300 reclassificado somente nos doadores de treinamento como análise principal. Não equivale a executar `dtangle(n_markers=100)` com todos os genes comuns.

## Configuração fornecida pelo orientador — diferença em relação ao nosso 08
- Mesma referência Mathys e mesmos oito tipos.
- Bulk e referência `log2(CPM+1)`.
- Um perfil **médio por tipo**, em vez de 181 perfis doador/tipo.
- **TODOS os genes compartilhados** (não top300).
- `n_markers=100` por tipo; conferir se os marcadores são escolhidos automaticamente por `marker_method="ratio"`, padrão do pacote, ou por outro critério.
- Confirmar a definição exata de “média por tipo” (soma de counts antes do CPM, média dos CPM por núcleo ou média dos perfis por doador); confirmar ainda `gamma`, normalização quantílica, versão do dtangle e denominadores de CPM.

## Novo teste 08e
O script `08e_reproduzir_dtangle_orientador.R` aplica todos os genes comuns, seleção automática de 100 marcadores por tipo, escala log2(CPM+1), `data_type="rna-seq"`, `marker_method="ratio"` e duas construções alternativas do perfil médio por tipo; preserva todos os resultados anteriores.

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(stringr)
  library(dplyr)
})

#Recebe os argumentos passados para esse script pelo terminal 
#trailingOnly= retorna apenas os argumentos fornecidos, ignorando os argumentos internos do script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Uso: Rscript 00_auditoria.R <counts.tsv.gz> <annot.tsv.gz> <soft.gz> <out_report.txt>") #Interrompe o processamento caso hajam arquivos inexistentes ou com ausência de algum dos argumentos 
}

#counts, annotation, arquivo soft e arquivo de saída são necessários para este script
counts_file <- args[1]
annot_file  <- args[2]
soft_file   <- args[3]
output_txt  <- args[4]

#Criação de um vetor vazio para elaborar um "rascunho" do relatório
report_lines <- c() #c=combine --> cria elementos em um vetor
add_to_report <- function(...) {
  report_lines <<- c(report_lines, paste(...))
}

#Cabeçalho do relatório, registrando a data atual e os arquivos utilizados
add_to_report("==================================================================")
add_to_report("           RELATÓRIO DE AUDITORIA DE INTEGRIDADE DE DADOS         ")
add_to_report("==================================================================")
add_to_report("Data: ", as.character(Sys.time()))
add_to_report("Counts: ", counts_file)
add_to_report("Annot:  ", annot_file)
add_to_report("SOFT:   ", soft_file)
add_to_report("------------------------------------------------------------------\n")

# 1. Auditoria de Counts
add_to_report(">>> 1. AUDITORIA DA MATRIZ DE COUNTS")
counts_df <- fread(cmd = paste("zcat", counts_file), data.table = FALSE) #O arquivo é descomprimido, a função fread() lê os arquivos como uma tabela e retorna um data.frame
sample_cols <- colnames(counts_df)[-1] #Pega todas as colunas, menos a primeira 
gene_ids <- counts_df[, 1] #pega todas as linhas da primeira coluna onde estão os identificadores dos genes

add_to_report(" - Nome da coluna de ID: ", colnames(counts_df)[1]) #Verifica qual é a primeira coluna
add_to_report(" - Total de genes (linhas): ", nrow(counts_df)) # Cada linha representa um gene, logo a soma corresponde a quantidade de genes
add_to_report(" - Total de amostras (colunas): ", length(sample_cols)) #Possui as colunas depois de GeneID, seu tamanho corresponde ao número de amostras 
add_to_report(" - GeneIDs duplicados: ", sum(duplicated(gene_ids))) #Verifica quais GeneIDs aparecem mais de uma vez através da soma de TRUE

counts_mat <- as.matrix(counts_df[, -1]) #A primeira coluna é removida e o restante é transformado em matriz, deixando apenas os valores numéricos de expressão
add_to_report(" - Valores ausentes (NA): ", sum(is.na(counts_mat))) #Se algum valor estiver ausente (NA), ele é marcado como TRUE. Depois os valores TRUE são somados. Para uma matriz de counts é esperado que os valores estejam preenchidos
add_to_report(" - Contém valores negativos: ", any(counts_mat < 0, na.rm = TRUE)) #Verifica se há algum valor negativo, pos counts brutos não devem ser negativos

lib_sizes <- colSums(counts_mat, na.rm = TRUE) #Representa o tamanho da biblioteca de cada amostra, ou seja, o total de reads/counts atribuídos aos genes na amostra
add_to_report(" - Mediana do tamanho de biblioteca: ", format(median(lib_sizes), big.mark=".")) #Calcula a mediana dos tamanhos de biblioteca. Os números são formatados sem o .
add_to_report("\n------------------------------------------------------------------\n")

# 2. Auditoria de Anotação
add_to_report(">>> 2. AUDITORIA DA ANOTAÇÃO GÊNICA")
annot_df <- fread(cmd = paste("zcat", annot_file), data.table = FALSE) #Mesmo mecanismo usado para ler os counts, mas agora para a tabela de anotação dos genes 
annot_id_col <- colnames(annot_df)[1] 

add_to_report(" - Total de genes anotados: ", nrow(annot_df)) #quantos genes existem na anotação -atraves da contagem de linhas
add_to_report(" - GeneIDs duplicados na anotação: ", sum(duplicated(annot_df[[annot_id_col]]))) #soma os IDs duplicados na anotação

common_genes <- intersect(as.character(gene_ids), as.character(annot_df[[annot_id_col]])) #Encontra genes presentes nos dois arquivos, tanto na matriz de counts, quanto na anotação, selecionando a interseção (como texto para evitar diferenças de tipo)
pct_mapped <- (length(common_genes) / length(gene_ids)) * 100 #Percentual mapeado, quantos % dos GeneIDs da matriz possuem correspondência com os da anotação
add_to_report(" - Mapeamento com a matriz: ", length(common_genes), sprintf(" (%.2f%%)", pct_mapped))
add_to_report("\n------------------------------------------------------------------\n")

# 3. Auditoria do SOFT
add_to_report(">>> 3. AUDITORIA DOS METADADOS (SOFT FILE)")
gds <- getGEO(filename = soft_file) #A função getGEO() lê o arquivo soft e transforma suas informações em objetos que o R consegue manipular
gsm_list <- GSMList(gds) #extraindo a lista de objetos correspondentes às amostras GSM
soft_gsms <- names(gsm_list) #extrai os nomes das amostras

add_to_report(" - Amostras no SOFT: ", length(soft_gsms)) #descreve quantas amostras estão descritas no SOFT
add_to_report(" - GSMs pareadas (SOFT e Counts): ", length(intersect(sample_cols, soft_gsms))) # comparação dos GSMs na matriz com os do SOFT, idealmente devem conter o mesmo número

missing_in_counts <- setdiff(soft_gsms, sample_cols) #significa que os elementos que estão em A, mas não em B
add_to_report(" - GSMs no SOFT ausentes nos Counts: ", length(missing_in_counts))
if (length(missing_in_counts) > 0) { #Verifica se pelo menos uma GSM está faltando
  add_to_report("   [Ausentes]: ", paste(missing_in_counts, collapse = ", ")) #Se tiver algum GSM ausente, ele é relatado
}

writeLines(report_lines, con = output_txt)
cat("[SUCCESS] Relatório de auditoria salvo em: ", output_txt, "\n")

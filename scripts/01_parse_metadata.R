#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(readr)
})
#O arquivo soft contêm os metadados do estudo, portanto iremos identificar a doença, a idade, sexo, RIN, PMI e região afetada através deste script
#Recebe os argumentos passados para esse script pelo terminal 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 01_parse_metadata.R <soft.gz> <counts.tsv.gz> <out_samples.csv>")
}

soft_file   <- args[1]
counts_file <- args[2]
output_csv  <- args[3]

#Lê a primeira linha da matriz e remove a primeira coluna (GeneID), deixando as amostras GSM 
counts_header <- colnames(fread(counts_file, nrows = 1, data.table = FALSE))
valid_samples <- counts_header[-1]

#Lê o arquivo soft e extrai as amostras GSM do objeto (cada GSM contém seus respectivos metadados)
gds <- getGEO(filename = soft_file, getGPL = FALSE)
gsm_list <- GSMList(gds)

#Recebe as características da amostra (chars) e o que estamos procurando (key_pattern; como "age", "sex"). O grep() procura padrões de texto e devolve o texto encontrado ao invés da posição
extract_char_value <- function(chars, key_pattern) {
  matched <- grep(key_pattern, chars, value = TRUE, ignore.case = TRUE) 
  #Se nenhuma caracteristica for encontrada, o valor estiver ausente, retorna um NA 
  if (length(matched) == 0) return(NA_character_)
  #Procura a informação, separa pelo : pegando o que vem após dos :
  val <- str_split_fixed(matched[1], ":", 2)[, 2]
  return(str_trim(val))
}

#Transformar as diferentes formas de descrever a doença em uma classificação padronizada, por exemplo, transformar todos "Parkinson's disease", "Parkinson" e "PD" somente em "PD"
harmonize_diagnosis <- function(raw_diag, title_text) {
  combined <- tolower(paste(raw_diag, title_text)) #juntamos o que está em raw_diag e title_tex e transforma tudo em minúsculas
  #Identificando controles
  if (str_detect(combined, "control|normal|non-demented|^c_")) { #o ^c significa comomeça com c_, dessa forma podemos reconhecer amostras como c_01 ou c_control
    return("Control")
    #Se não for classificado como controle, vamos identificar a qual das doenças pertence
  } else if (str_detect(combined, "alzheimer|ad")) {
    return("AD")
  } else if (str_detect(combined, "huntington|hd")) {
    return("HD")
  } else if (str_detect(combined, "parkinson|pd|^p_")) {
    return("PD")
  } else {
    return("Other") #Caso nenhum dos casos se aplique, será salvo em "outros casos"
  }
}

#Criamos uma lista para depois podermos combinar todas as "tabelas"
meta_rows <- list()

for (gsm_id in names(gsm_list)) {
  if (!gsm_id %in% valid_samples) next #verifica se a GSM está dentro de valid_samples, "essa GSM não está na matriz de counts?", se sim o loop prossegue para a próxima GSM 
  
  gsm <- gsm_list[[gsm_id]] #gsm pega os metadados da amostra 
  meta <- Meta(gsm) #extrai os metadados do objeto GSM
  chars <- meta$characteristics_ch1 #separa as características da amostra 
  title_val <- meta$title #separa o título da amostra 
  
  #Dentro das características, procuramos alguma variável relacionada a doença/diagnóstico/status/condição. 
  #Depois informações relacionadas a idade, sexo, RIN (RNA Integrity Number), PMI (post-mortem interval), região do tecido
  raw_disease <- extract_char_value(chars, "disease|diagnosis|status|condition|subject state")
  age_val     <- extract_char_value(chars, "age")
  sex_val     <- extract_char_value(chars, "sex|gender")
  rin_val     <- extract_char_value(chars, "rin|rna integrity|integrity")
  pmi_val     <- extract_char_value(chars, "pmi|post-mortem|postmortem")
  region_val  <- extract_char_value(chars, "brain region|region|tissue|source")
  
  #Cria-se uma pequena tabela com uma linha para a GSM com as colunas sample_id, title, group, raw_diagnosis, age, sex, rin, pmi, region
  meta_rows[[gsm_id]] <- data.frame(
    sample_id     = gsm_id, #guarda o identificador original
    title         = title_val, #guarda o título original da amostra, permitindo rastreabilidade
    group         = harmonize_diagnosis(raw_disease, title_val), #harmonizamos as nomenclaturas das doenças
    raw_diagnosis = ifelse(is.na(raw_disease), title_val, raw_disease), #se o SOFT possui explicitamente o nome da doença (raw_disease) usamos essa nomenclatura, se não houver informações, usamos o título da amostra para preservar a informação e não deixar a coluna vazia
    age           = as.numeric(str_extract(age_val, "\\d+(\\.\\d+)?")), #transforma a idade no formato de texto em número
    sex           = ifelse(str_detect(tolower(sex_val), "^m|male"), "M", ifelse(str_detect(tolower(sex_val), "^f|female"), "F", NA_character_)), #padroniza a nomenclatura para sexo
    rin           = as.numeric(str_extract(rin_val, "\\d+(\\.\\d+)?")), #extrai o valor de RIN e transforma-o em valor numérico
    pmi           = as.numeric(str_extract(pmi_val, "\\d+(\\.\\d+)?")), #extrai o valor de PMI e transforma-o em valor numérico
    region        = region_val, 
    stringsAsFactors = FALSE #mantêm os textos como strings em vez de transformá-los automaticamente em fatores. 
  )
}

#No final temos uma tabela com colunas de metadados, onde cada linha contêm as informações de uma GSM. Cria uma versão padronizada do soft
sample_metadata <- do.call(rbind, meta_rows) #juntamos todas as amostras pelas linhas 
write_csv(sample_metadata, output_csv)
cat("[SUCCESS] Metadados extraídos e salvos em: ", output_csv, "\n")
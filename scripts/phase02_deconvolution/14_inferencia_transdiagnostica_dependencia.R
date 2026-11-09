#!/usr/bin/env Rscript
# =============================================================================
# 14 | INFERENCIA TRANSDIAGNOSTICA CORRIGIDA PELOS CONTROLES COMPARTILHADOS
# TCC AD / HD / PD, dtangle 08e principal; BisqueRNA e MuSiC2 de sensibilidade.
# GSE125583 AD principal; ROSMAP AD validacao anatomica.
#
# Objetivos COMPLEMENTARES (nao confundir):
# (1) GLS de efeito COMUM medio com covariancia HC3 entre efeitos HD/PD,
#     reconstruida a partir de individuos partilhados (influencias do OLS).
# (2) Teste de CONJUNCAO DIRECIONAL: evidencia do efeito com MESMO SINAL
#     nas 3 doencas; validado sob dependencia arbitraria dos p-valores.
#
# GLS nao e' um modelo de efeitos aleatorios; Q avalia heterogeneidade com
# somente 3 doencas, possui baixo poder e nao estabelece efeito universal.
# A inferencia mede contribuicao relativa de RNA, nao numero absoluto de celulas.
#
# Requer entradas criadas nos scripts 09, 10, 11 (v3) e 11b.
# Nenhuma matrix bruta sera reprocessada e nenhum arquivo anterior sobrescrito.
# Executar na raiz do repositorio.
# =============================================================================
options(stringsAsFactors = FALSE)
if (!requireNamespace('data.table', quietly = TRUE)) stop('Falta data.table')
if (!requireNamespace('ggplot2', quietly = TRUE)) stop('Falta ggplot2')
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- normalizePath(Sys.getenv('DECONV_ROOT', '.'), mustWork = TRUE)
P <- function(...) file.path(ROOT, ...)
OUT <- P('results', 'deconvolution', 'Integrated_Comparison', '14_joint')
PLOTS <- file.path(OUT, 'plots')
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
STYLE <- P('scripts', 'visualization', '00_plot_style.R')
if (!file.exists(STYLE)) stop('Modulo de estilo obrigatorio nao encontrado: ', STYLE)
source(STYLE, local = TRUE)

INPUT <- P('results','deconvolution','Integrated_Comparison','09',
           '09_proportions_long_harmonized.tsv')
PLAN <- P('data','processed','deconvolution','10_approved_model_plan.tsv')
META <- P('results','deconvolution','Integrated_Comparison','10_metadata')
OLD_INDIV <- P('results','deconvolution','Integrated_Comparison','11_adjusted',
               '11_adjusted_pp_HC3.tsv')
OLD_GROUP <- P('results','deconvolution','Integrated_Comparison','11b_grouped',
               '11b_grouped_adjusted_pp_HC3.tsv')
REQUIRED_FILES <- c(INPUT, PLAN, OLD_INDIV, OLD_GROUP)
not_found <- REQUIRED_FILES[!file.exists(REQUIRED_FILES)]
if (length(not_found)) stop('Faltam entradas obrigatorias:\n',paste(not_found,collapse='\n'))

DS <- c('GSE125583', 'ROSMAP', 'GSE64810', 'GSE68719')
METHODS <- c('dtangle', 'BisqueRNA', 'MuSiC2')
CELLS <- c('Ex','In','Ast','Oli','Opc','Mic','End','Per')
GROUPS <- list(Neuronios=c('Ex','In'),Oligodendroglia=c('Oli','Opc'),
               Vascular=c('End','Per'))
AD_SCENARIOS <- c(principal='GSE125583', validacao='ROSMAP')
DISEASE <- c(GSE125583='AD', ROSMAP='AD', GSE64810='HD', GSE68719='PD')
write_tsv <- function(x, name) {
  fwrite(x, file.path(OUT, name), sep='\t', na='NA', quote=FALSE)
}
fields <- function(s) {
  if (length(s)!=1L || is.na(s) || !nzchar(trimws(s))) return(character())
  trimws(strsplit(s, ',', fixed=TRUE)[[1L]])
}
nonblank <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))

D <- fread(INPUT, na.strings=c('NA','NaN',''))
req <- c('dataset','method','sample_id','cell_type','group','proportion')
if (!all(req %in% names(D))) stop('09 harmonizado sem colunas: ',paste(setdiff(req,names(D)),collapse=','))
if (anyDuplicated(D, by=c('dataset','method','sample_id','cell_type'))) stop('Duplicatas no 09')
if (!setequal(unique(D$dataset),DS) || !setequal(unique(D$method),METHODS) ||
    !setequal(unique(D$cell_type),CELLS)) stop('Conjuntos inesperados na entrada 09')
if (anyNA(D$proportion) || any(!is.finite(D$proportion)) ||
    any(D$proportion<0 | D$proportion>1)) stop('Fracoes invalidas')
cdf <- D[,.(n_cells=.N,sum_fraction=sum(proportion)), by=.(dataset,method,sample_id)]
if (any(cdf$n_cells!=8L) || any(abs(cdf$sum_fraction-1)>1e-6))
  stop('Proporcoes incompletas ou soma diferente de 1')

plan <- fread(PLAN, na.strings=c('NA',''))
reqplan <- c('dataset','approved','numeric_covariates','categorical_covariates',
             'cluster_column','justification','provenance')
if (!all(reqplan %in% names(plan)) || nrow(plan)!=4L ||
    !setequal(plan$dataset,DS) || anyDuplicated(plan,by='dataset'))
  stop('Plano de ajustes incompleto/invalido')
if (any(toupper(trimws(plan$approved))!='YES')) stop('Todas as quatro coortes precisam estar aprovadas')
if (any(nonblank(plan$cluster_column))) stop('Plano exige modelo clusterizado; script 14 espera OLS HC3')

OLD11 <- fread(OLD_INDIV, na.strings=c('NA','NaN',''))
OLD11B <- fread(OLD_GROUP, na.strings=c('NA','NaN',''))
if (!all(c('dataset','method','cell_type','scale','effect','robust_se','status') %in% names(OLD11)))
  stop('Tabela da etapa 11 inadequada')
if (!all(c('dataset','method','group_type','scale','effect','robust_se','status') %in% names(OLD11B)))
  stop('Tabela da etapa 11b inadequada')
OLD11 <- OLD11[scale=='percentage']
OLD11B <- OLD11B[scale=='percentage']
if (nrow(OLD11)!=4L*3L*8L || nrow(OLD11B)!=4L*3L*3L) stop('11/11b incompletos: reexecute etapas anteriores')

# OLS HC3: retorna componente de influencia por individuo, cuja soma dos
# quadrados e' exatamente a variancia HC3 do coeficiente de diagnostico.
fit_hc3_influence <- function(y, X, sample_id) {
  n <- nrow(X); p <- ncol(X)
  if (n<=p+5L || qr(X)$rank!=p) stop('Matriz nao estimavel')
  Q <- chol2inv(chol(crossprod(X)))
  beta <- as.numeric(Q %*% crossprod(X,y))
  names(beta) <- colnames(X)
  residual <- as.numeric(y-X%*%beta)
  lever <- rowSums((X%*%Q)*X)
  if (any(lever>=1-1e-8)) stop('Alavancagem inviabiliza HC3')
  idx <- which(colnames(X)=='case')
  if (length(idx)!=1L) stop('Coeficiente case ausente')
  infl <- as.numeric((X%*%Q)[,idx] * residual/(1-lever))
  names(infl) <- sample_id
  se <- sqrt(sum(infl^2))
  b <- unname(beta[idx])
  df <- n-p
  if (!is.finite(se) || se<=0 || !is.finite(b)) stop('SE/efeito nao estimavel')
  list(effect=b, se=se, df=df, p=2*pt(-abs(b/se),df),
       influence=infl, n=n, max_leverage=max(lever))
}

fit_cache <- list(); audit <- list(); matches <- list(); idx <- 0L
for (ds in DS) {
  pp <- plan[dataset==ds]
  if (!nonblank(pp$justification) || !nonblank(pp$provenance))
    stop('Plano sem justificativa/procedencia: ',ds)
  num <- fields(pp$numeric_covariates)
  fac <- fields(pp$categorical_covariates)
  cv <- c(num,fac)
  if (!length(cv) || anyDuplicated(cv)) stop('Covariaveis invalidas: ',ds)
  mf <- file.path(META, paste0('10_metadata_',ds,'_audited.tsv'))
  if (!file.exists(mf)) stop('Metadados ausentes: ',mf)
  md <- fread(mf, na.strings=c('NA','','NaN','NULL'))
  if (!all(c('sample_id','group',cv)%in%names(md)) || anyDuplicated(md,by='sample_id'))
    stop('Metadados mal formados: ',ds)
  md <- md[, c('sample_id','group',cv), with=FALSE]
  if (!setequal(md$sample_id,D[dataset==ds & method=='dtangle' & cell_type=='Ex',sample_id]))
    stop('Metadados e proporcoes nao correspondem: ',ds)
  if (!setequal(unique(md$group),c('Control',DISEASE[[ds]]))) stop('Grupos inconsistentes: ',ds)
  for (v in num) {
    if (!is.numeric(md[[v]]) || any(!is.finite(md[[v]][!is.na(md[[v]])])))
      stop('Covariavel numerica invalida: ',ds,'/',v)
  }
  for (v in fac) md[[v]] <- factor(trimws(as.character(md[[v]])))
  keep <- complete.cases(md)
  for (v in cv) keep <- keep & nonblank(md[[v]])
  m <- copy(md[keep])
  for (v in fac) m[[v]] <- droplevels(m[[v]])
  m[,case:=as.integer(group==DISEASE[[ds]])]
  if (nrow(m)<30L || sum(m$case==0)<10L || sum(m$case==1)<10L)
    stop('Amostras completas insuficientes: ',ds)
  X <- model.matrix(reformulate(c('case',cv)), data=as.data.frame(m))
  if (!'case' %in% colnames(X) || qr(X)$rank!=ncol(X) ||
      nrow(X)<max(30L,10L*ncol(X))) stop('Matriz singular/superajustada: ',ds)
  if (max(rowSums((X%*%chol2inv(chol(crossprod(X))))*X))>=1-1e-8)
    stop('Alavancagem HC3 inestimavel: ',ds)
  audit[[ds]] <- data.table(dataset=ds,n_used=nrow(m),n_control=sum(m$case==0),
    n_case=sum(m$case==1),n_parameters=ncol(X),n_rank=qr(X)$rank,
    covariates=paste(cv,collapse=';'),status='PASS_INPUT')

  for (mt in METHODS) {
    # mt nao coincide com nenhuma coluna, evitando captura de escopo no data.table.
    dd <- D[dataset==ds & method==mt]
    wide <- dcast(dd, sample_id~cell_type, value.var='proportion')
    if (!setequal(wide$sample_id,md$sample_id)) stop('ID divergente: ',ds,' ',mt)
    wide <- wide[match(m$sample_id,wide$sample_id)]
    if (!identical(wide$sample_id,m$sample_id)) stop('Alinhamento incorreto: ',ds)
    M <- as.matrix(wide[,..CELLS]);storage.mode(M)<-'double'
    if (anyNA(M) || any(abs(rowSums(M)-1)>1e-6)) stop('Composicao invalida: ',ds,'/',mt)
    Ys <- as.data.frame(100*M)
    for (grp in names(GROUPS)) {
      Ys[[grp]] <- 100*rowSums(M[,GROUPS[[grp]],drop=FALSE])
    }
    for (endpoint in c(CELLS,names(GROUPS))) {
      fam <- if(endpoint %in% CELLS) 'cell' else 'cell_group'
      y <- Ys[[endpoint]]
      if(endpoint %in% CELLS){
        zero <- mean(M[,endpoint]<=1e-12)
      } else {
        zero <- mean(y<=1e-10)
      }
      sdpp <- sd(y)
      state <- if(zero>=0.90 || sdpp<0.20) 'SKIP_COLLAPSED' else 'PASS'
      z <- NULL
      if(state=='PASS') {
        z <- tryCatch(fit_hc3_influence(y,X,m$sample_id),error=function(e)e)
        if (inherits(z,'error')) stop('14 falhou em ',ds,'/',mt,'/',endpoint,': ',conditionMessage(z))
      }
      old <- if(fam=='cell') {
        OLD11[dataset==ds & method==mt & cell_type==endpoint]
      } else {
        OLD11B[dataset==ds & method==mt & group_type==endpoint]
      }
      if (nrow(old)!=1L) stop('Resultado anterior nao unico: ',ds,'/',mt,'/',endpoint)
      if (identical(state=='PASS',old$status=='PASS')==FALSE)
        stop('14/11 divergiram no status de ',ds,'/',mt,'/',endpoint)
      if (state=='PASS') {
        tol <- 1e-5
        if (!is.finite(old$effect) || !is.finite(old$robust_se) ||
            abs(old$effect-z$effect)>tol || abs(old$robust_se-z$se)>tol)
          stop('14 nao reproduziu coeficiente/HC3 do script anterior: ',ds,'/',mt,'/',endpoint)
      }
      key <- paste(ds,mt,endpoint,sep='|')
      fit_cache[[key]] <- z
      idx <- idx+1L
      matches[[idx]] <- data.table(dataset=ds,method=mt,
        family=fam,endpoint=endpoint,status=state,
        n_used=nrow(m), effect=if(is.null(z))NA_real_ else z$effect,
        robust_se=if(is.null(z))NA_real_ else z$se,
        p_hc3=if(is.null(z))NA_real_ else z$p,
        zero_fraction=zero,sd_pp=sdpp,
        agrees_with_previous=if(is.null(z))NA else TRUE)
    }
  }
}
write_tsv(rbindlist(audit),'14_design_checks.tsv')
write_tsv(rbindlist(matches),'14_model_reproduction_checks.tsv')

# Shared controls must be identical sample IDs with group Control in HD/PD;
# no other subjects overlap between AD and HD/PD.
get_controls <- function(ds) D[dataset==ds & method=='dtangle' & cell_type=='Ex' & group=='Control',sample_id]
get_all <- function(ds) D[dataset==ds & method=='dtangle' & cell_type=='Ex',sample_id]
shared <- intersect(get_controls('GSE64810'),get_controls('GSE68719'))
other_overlap <- intersect(get_all('GSE64810'),get_all('GSE68719'))
if (length(shared)!=44L || !setequal(shared,other_overlap))
  stop('Controles HD/PD nao correspondem exatamente a 44 IDs unicos: audite o 10c')
for(ad in c('GSE125583','ROSMAP')) {
  if(length(intersect(get_all(ad),union(get_all('GSE64810'),get_all('GSE68719')))))
    stop('Sobreposicao inesperada entre AD e HD/PD: ',ad)
}

pool_gls <- function(b, V) {
  if(any(!is.finite(V)) || max(abs(V-t(V)))>1e-8 || any(eigen(V,symmetric=TRUE,only.values=TRUE)$values<=0))
    stop('Matriz de covariancia nao definida positiva')
  Q <- solve(V)
  ones <- rep(1,length(b))
  den <- as.numeric(crossprod(ones,Q%*%ones))
  theta <- as.numeric(crossprod(ones,Q%*%b)/den)
  se <- sqrt(1/den)
  hetero_Q <- as.numeric(crossprod(b-theta,Q%*%(b-theta)))
  c(effect=theta,robust_se=se,ci_low=theta-qnorm(.975)*se,
    ci_high=theta+qnorm(.975)*se,p_gls=2*pnorm(-abs(theta/se)),
    Q=hetero_Q,p_Q=pchisq(hetero_Q,df=length(b)-1,lower.tail=FALSE))
}

joint <- list(); component <- list(); covar <- list(); ix <- 0L
for (scenario in names(AD_SCENARIOS)) {
  ads <- AD_SCENARIOS[[scenario]]
  three <- c(AD=ads,HD='GSE64810',PD='GSE68719')
  for (method in METHODS) {
    for (endpoint in c(CELLS,names(GROUPS))) {
      fam <- if(endpoint %in% CELLS) 'cell' else 'cell_group'
      fit <- lapply(three, function(ds) fit_cache[[paste(ds,method,endpoint,sep='|')]])
      ok <- vapply(fit,function(z)is.list(z)&&!inherits(z,'error'),logical(1))
      if(!all(ok)) {
        joint[[length(joint)+1L]] <- data.table(scenario=scenario,AD_dataset=ads,
          method=method,family=fam,endpoint=endpoint,status='SKIP_INCOMPLETE_THREE_DISEASES',
          n_diseases=sum(ok),gls_effect_pp=NA_real_,gls_se_pp=NA_real_,
          gls_ci_low_pp=NA_real_,gls_ci_high_pp=NA_real_,p_gls=NA_real_,
          Q_heterogeneity=NA_real_,p_heterogeneity=NA_real_,
          shared_controls=44L,corr_HD_PD=NA_real_,same_sign=NA,
          direction=NA_character_,p_conjunction=NA_real_)
        next
      }
      b <- vapply(fit,`[[`,numeric(1),'effect')
      se <- vapply(fit,`[[`,numeric(1),'se')
      p_two <- vapply(fit,`[[`,numeric(1),'p')
      df <- vapply(fit,`[[`,numeric(1),'df')
      V <- diag(se^2,nrow=3,ncol=3)
      h <- fit[['HD']]$influence[shared]
      p <- fit[['PD']]$influence[shared]
      if(anyNA(h)||anyNA(p))stop('Influencias HD/PD nao mapeadas')
      cov_hp <- sum(h*p)
      V[2,3] <- V[3,2] <- cov_hp
      corr_hp <- cov_hp/(se[2]*se[3])
      if(!is.finite(corr_hp)||abs(corr_hp)>=1)
        stop('Correlacao HD/PD extrema/invalida: ',method,'/',endpoint)
      g <- pool_gls(b,V)
      tstats <- b/se
      # Para efeito no mesmo sentido em TODAS as 3 doencas: IUT direcional.
      # p(sign positivo) = max(p_one_sided_positive_i).
      # p(sign negativo) = max(p_one_sided_negative_i).
      # Direcao escolhida apos observar os dados -> fator Bonferroni de 2.
      p_plus <- max(pt(-tstats,df))
      p_minus <- max(pt(tstats,df))
      p_conj <- min(1,2*min(p_plus,p_minus))
      signs <- sign(b)
      same <- all(signs>0) || all(signs<0)
      direct <- if(all(signs>0))'increase_all' else if(all(signs<0))'decrease_all' else 'mixed'
      ix <- ix+1L
      covar[[ix]] <- data.table(scenario=scenario,method=method,
        family=fam,endpoint=endpoint,n_shared=length(shared),
        covariance_HD_PD=cov_hp,correlation_HD_PD=corr_hp,
        var_AD=V[1,1],var_HD=V[2,2],var_PD=V[3,3],
        min_eigenvalue=min(eigen(V,symmetric=TRUE,only.values=TRUE)$values))
      for(dis in names(three)){
        component[[length(component)+1L]] <- data.table(
          scenario=scenario,AD_dataset=ads,method=method,family=fam,
          endpoint=endpoint,disease=dis,dataset=three[[dis]],
          effect_pp=b[[dis]],se_pp=se[[dis]],p_hc3=p_two[[dis]],df_resid=df[[dis]])
      }
      joint[[length(joint)+1L]] <- data.table(
        scenario=scenario,AD_dataset=ads,method=method,family=fam,
        endpoint=endpoint,status='PASS_GLS_AND_CONJUNCTION',n_diseases=3L,
        gls_effect_pp=unname(g['effect']),gls_se_pp=unname(g['robust_se']),
        gls_ci_low_pp=unname(g['ci_low']),gls_ci_high_pp=unname(g['ci_high']),
        p_gls=unname(g['p_gls']),Q_heterogeneity=unname(g['Q']),
        p_heterogeneity=unname(g['p_Q']),shared_controls=length(shared),
        corr_HD_PD=corr_hp,same_sign=same,direction=direct,
        p_conjunction=p_conj,
        max_component_p=max(p_two),
        inferential_warning='GLS assumes common effect; 3 diseases; distinct covariate sets; low-abundance types cautious')
    }
  }
}

J <- rbindlist(joint,fill=TRUE)
J[,`:=`(fdr_gls=NA_real_,fdr_conjunction=NA_real_,
         fdr_by_gls=NA_real_,fdr_by_conjunction=NA_real_)]
J[status=='PASS_GLS_AND_CONJUNCTION',`:=`(
  fdr_gls=p.adjust(p_gls,method='BH'),
  fdr_conjunction=p.adjust(p_conjunction,method='BH'),
  fdr_by_gls=p.adjust(p_gls,method='BY'),
  fdr_by_conjunction=p.adjust(p_conjunction,method='BY')),
  by=.(scenario,method,family)]
# p.adjust chamado separadamente por 8 tipos e 3 agrupamentos, sem misturar testes.
write_tsv(J,'14_joint_inference_all.tsv')
write_tsv(J[method=='dtangle'],'14_dtangle_primary_and_ROSMAP.tsv')
write_tsv(rbindlist(component,fill=TRUE),'14_component_effects.tsv')
write_tsv(rbindlist(covar,fill=TRUE),'14_HD_PD_shared_control_covariance.tsv')

# Resumo de conclusao: p_gls significativo != p_conj significativo.
J[,classification:=fifelse(status!='PASS_GLS_AND_CONJUNCTION','NOT_TESTABLE',
  fifelse(fdr_conjunction<.05 & same_sign,'ALL_THREE_DIRECTIONAL_SUPPORT',
  fifelse(fdr_gls<.05,'AVERAGE_GLS_ONLY_OR_NO_CONJUNCTION','NO_FDR_SUPPORT')))]
write_tsv(J,'14_joint_inference_with_classification.tsv')

# Graficos padronizados para o metodo principal (cores constantes por celula).
for(sc in names(AD_SCENARIOS)) {
  for(fa in c('cell','cell_group')) {
    z <- J[scenario==sc & method=='dtangle' & family==fa &
             status=='PASS_GLS_AND_CONJUNCTION']
    if(!nrow(z))next
    z[,endpoint:=factor(endpoint,levels=if(fa=='cell')CELLS else names(GROUPS))]
    z[,conclusion:=fifelse(fdr_conjunction<.05,'FDR conjuncao < 0,05','Nao significativo')]
    pal <- if(fa=='cell')'cell' else 'cell_group'
    pp <- ggplot(z,aes(x=gls_effect_pp,y=endpoint,color=endpoint)) +
      geom_vline(xintercept=0,linetype='dashed',color='#777777') +
      geom_segment(aes(x=gls_ci_low_pp,xend=gls_ci_high_pp,
                       y=endpoint,yend=endpoint),linewidth=.65) +
      geom_point(aes(shape=conclusion),size=2.7) +
      tcc_scale_color(pal) +
      scale_shape_manual(values=c('FDR conjuncao < 0,05'=16,'Nao significativo'=1)) +
      labs(x='Efeito medio conjunto (pontos percentuais)',
           y=if(fa=='cell')'Tipo celular' else 'Compartimento celular',
           title=paste('Inferencia conjunta AD–HD–PD |',sc),
           subtitle='GLS: covariancia HD–PD por controles compartilhados; conjuncao testada separadamente',
           color=if(fa=='cell')'Tipo celular' else 'Compartimento',shape='Teste de conjuncao') +
      tcc_theme()
    tcc_assert_axes(pp)
    tcc_save_plot(pp,file.path(PLOTS,paste0('14_',sc,'_',fa,'_GLS_300dpi.png')),
                  width_mm=185,height_mm=125)
  }
}

# Inventario de entradas, reproducibilidade e metodologia.
paths <- c(INPUT,PLAN,OLD_INDIV,OLD_GROUP,
           P('scripts','phase02_deconvolution','14_inferencia_transdiagnostica_dependencia.R'),
           vapply(DS,function(ds) file.path(META,paste0('10_metadata_',ds,'_audited.tsv')),''),STYLE)
write_tsv(data.table(path=paths,md5=as.character(tools::md5sum(paths))),
          '14_input_manifest_md5.tsv')
writeLines(c(
 '# Inferencia transdiagnostica - metodo',
 '',
 '- AD principal GSE125583; AD validacao ROSMAP, avaliados em cenarios SEPARADOS.',
 '- HD GSE64810 e PD GSE68719 compartilham 44 controles (sample_id exato).',
 '- Para cada coorte, metodo e tipo, o OLS percentual foi reajustado com as mesmas covariaveis do script 11.',
 '- A matriz V (3x3) usa variancias HC3 nas diagonais e covariancia HD-PD calculada como soma do produto das influencias individuais HC3 dos controles compartilhados.',
 '- Estimativa conjunta e teste bilateral GLS usam V completa. Esse eh um modelo de efeito COMUM, nao efeitos aleatorios.',
 '- A heterogeneidade Q (2 gl) foi testada como diagnostico exploratorio com baixo poder (apenas 3 doencas).',
 '- Teste de conjuncao bidirecional: 2 * min(max(p_i^+), max(p_i^-)), limitado a 1; tres efeitos na mesma direcao. Valido sob dependencia arbitraria se p individuais validos.',
 '- FDR-BH e FDR-BY por metodo, cenario AD e familia de desfechos: 8 tipos individuais ou 3 agrupamentos.',
 '- p_GLS significativo apenas sugere efeito medio diferente de zero; nao basta para concluir efeito em TODAS as doencas.',
 '- Balanço de covariaveis e anatomia diferem entre coortes; efeitos podem ser heterogeneos.',
 '- MuSiC2: In/Opc degenerados, nao resgatados pela agregacao; baixa abundancia de End/Per e Mic exige cautela.',
 '- Sem correcao da incerteza das fracoes de entrada (deconvolucao).',
 '- Conferir que o arquivo 09 harmonizado realmente utiliza dtangle 08e pooled_counts antes da interpretacao final.',
 '- Graficos padronizados do script 14 usam o modulo 00_plot_style.R, resolucao >=300 dpi, eixos nomeados, cores constantes por tipo celular.',
 ''),file.path(OUT,'14_METHODS_AND_LIMITATIONS.md'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'14_sessionInfo.txt'))
cat('14 CONCLUIDO | entradas reproduzidas, covariancia HD/PD calculada; resultados em: ',OUT,'\n',sep='')

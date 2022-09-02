##################################################################### #
##' Função para criar tabela de modelos aditivos
##' 
##' Autores: Mikael e Matheus
##' Data: 03/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----
source("R/fct_aux/func_create_table.R")
source("R/fct_aux/func_check_return.R")
source("R/fct_aux/func_testes_chi2_fisher.R")

# 1 - Função ----
fct_qui2_fisher <- function(df_full, path_tables, debug){
  # 2 - Qui-quadrado e teste de fisher ----
  
  ## 2.1 - Teste 1 (Presença) ----
  ## Teste 1 será rodado para agrupamento não possui (ausencia) VS possui (presenca) de mucosite bocal
  ## ausencia (PIORMB = 0) vs presença (PIORMB = 1, 2, 3, 4)
  ## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")
  
  ## DF agrupado para ausência (PIORMB = 0) ou presença (PIORMB = 1, 2, 3, 4)
  df_abs_or_pres <- df_full |> 
    dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD", "MU"))) |> 
    dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0))
  
  # df_abs_or_pres_aux <- df_abs_or_pres
  df_abs_or_pres[sapply(df_abs_or_pres, is.numeric)] <- lapply(df_abs_or_pres[sapply(df_abs_or_pres, is.numeric)], 
                                         as.factor)
  
  ### 2.1.1 Rodando testes chi-2 e de fisher (ausencia vs presença) ----
  list_testes <- fct_testes_chi2_fisher(df_abs_or_pres)
  
  if(list_testes[[1]]){
    ## Testes com e sem simulação de Monte Carlo
    df_chi2_a_p <- list_testes[[2]]
    df_fisher_a_p <- list_testes[[3]]
  }
  
  ### 2.1.2 - Avaliando modelos (dominante, recessivo, aditivo) ----
  df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])
  
  #### 2.1.2.1 - Gerando tabelas para modelos agrupados por ausência vs presença ----
  ## Adicionando Tabelas de contingência
  
  type_group = "pres"
  table <- NULL
  
  if(switch_write_table){
    return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_a_p, table,
                                            type_group, path_tables,
                                            use_MC, write_table,switch_teste)
    fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)
  }
  
  ## 2.2 - Teste 2 (Ulcerações) ----
  ## Teste 2 será rodado para agrupamento não ulcerados vs ulcerados
  ## não ulcerado (PIORMB = 0, 1) ou ulcerado (PIORMB = 2, 3, 4)
  ## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")
  
  ## DF agrupado para nao ulcerado (PIORMB = 0, 1) ou ulcerado (PIORMB = 2, 3, 4)
  df_ulc <- df_full |> 
    dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD", "MU"))) |> 
    dplyr::mutate(PIORMB = ifelse(PIORMB > 1, 1, 0))
  
  df_ulc[sapply(df_ulc, is.numeric)] <- lapply(df_ulc[sapply(df_ulc, is.numeric)], 
                                                               as.factor)
  
  ### 2.2.1 Rodando testes chi-2 e de fisher (não ulcerados vs ulcerados) ----
  list_testes <- fct_testes_chi2_fisher(df_ulc)
  
  if(list_testes[[1]]){
    ## Testes com e sem simulação de Monte Carlo
    df_chi2_u <- list_testes[[2]]
    df_fisher_u <- list_testes[[3]]
  }
  
  ### 2.2.2 - Avaliando modelos (dominante, recessivo, aditivo) ----
  # df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])
  
  #### 2.2.2.1. Gerando tabelas para modelo não ulcerados vs ulcerados ----
  ## Adicionando Tabelas de contingência
  
  type_group = "ulc"
  table <- NULL
  
  if(switch_write_table){
    return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_u, table,
                                            type_group, path_tables,
                                            use_MC, write_table, switch_teste)
    fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)
  }
  
  
  ## 2.3 - Teste 3 (Severidade) ----
  ## Teste 2 será rodado para agrupamento não severo vs severo
  ## não severo (PIORMB = 0, 1, 2) ou severo (PIORMB = 3, 4)
  ## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")
  
  ## DF agrupado para nao severo (PIORMB = 0, 1, 2) ou severo (PIORMB = 3, 4)
  df_sev <- df_full |> 
    dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD", "MU"))) |> 
    dplyr::mutate(PIORMB = ifelse(PIORMB > 2, 1, 0))
  
  df_sev[sapply(df_sev, is.numeric)] <- lapply(df_sev[sapply(df_sev, is.numeric)], 
                                               as.factor)
  
  ### 2.3.1 Rodando testes chi-2 e de fisher (não severo vs severo) ----
  list_testes <- fct_testes_chi2_fisher(df_sev)
  
  if(list_testes[[1]]){
    ## Testes com e sem simulação de Monte Carlo
    df_chi2_s <- list_testes[[2]]
    df_fisher_s <- list_testes[[3]]
  }
  
  
  ### 2.3.2 - Avaliando modelos (dominante, recessivo, aditivo) ----
  # df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])
  
  #### 2.3.2.1 - Gerando tabelas para modelo mb severo vs não severo ----
  ## Adicionando Tabelas de contingência
  
  type_group = "sev"
  table <- NULL
  
  if(switch_write_table){
    return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_s, table, 
                                            type_group, path_tables,
                                            use_MC, write_table, switch_teste)
    fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)
  }
}
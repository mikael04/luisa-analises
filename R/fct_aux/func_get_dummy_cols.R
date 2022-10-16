######################################################################### #
##' Script para gerar novo dataframe com dummies nas colunas que precisarem
##' 
##' Autores: Mikael e Matheus
##' Data: 16/10/2022
######################################################################### #

# 0 - Scripts e bibliotecas ----
source("R/fct_aux/func_testes_chi2_fisher.R")

# 1 - Funções ----
fct_get_dummy_cols <- function(df_model, debug){
  teste_interno <- F
  if(teste_interno){
    df_model <- df_ctx_pres_aus
  }
  dummy_cols <- NULL
  for(i in 2:ncol(df_model)){
    col_tested <- as.character(df_model[[i]])
    if(length(unique(col_tested)) >= 3){
      # names(df_model)[i]
      dummy_cols <- append(dummy_cols, names(df_model)[i])
    }
  }
  # dummy_cols <- fct_get_dummy_cols(df_ctx_pres_aus, debug)
  dummy_cols_all <- names(df_model)[2:length(df_model)]
  
  df_model_char <- df_model |>
    dplyr::mutate(across(everything(), as.character))
  
  df_model_fast_dummies <- fastDummies::dummy_cols(df_model_char, select_columns = dummy_cols_all) |>
    dplyr::select(PIORMB, all_of(contains(c("_1", "_2")))) |> 
    dplyr::mutate(across(everything(), as.factor))
  
  df_model_fast_dummies
  }

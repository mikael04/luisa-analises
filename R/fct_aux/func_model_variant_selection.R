######################################################################### #
##' Script para selecionar as variantes que serão utilizadas para cada
##' modelo (ctx, doxo, mtx, outros) e agrupamento (sev, ulc, pres/aus)
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################################### #

# 0 - Scripts e bibliotecas ----
source("R/fct_aux/func_testes_chi2_fisher.R")

# 1 - Funções ----

## CTX
fct_mod_var_ctx_ulc <- function(){
  c("ABCC6_rs12931472_MD", "GSTA1_rs1051775_MA", "HSP90AA1_rs10873531_MA", "HSP90AA1_rs4947_MA",
    "SLCO6A1_rs6884141_MR")
}

fct_mod_var_ctx_sev <- function(){
  c()
}

fct_mod_var_ctx_pres <- function(){
  c("ABCC2_rs2273697_MD", "ABCC6_rs9940825_MD", "GSTP1_rs1695_MA", "HSP90AA1_rs4947_MA",
    "HSP90AA1_rs8005905_MA", "SLC19A1_rs1051266_MR", "SLC19A1_rs12659_MR")
}

## DOXO
fct_mod_var_doxo_ulc <- function(){
  c("ABCC2_rs2273697_MD", "ABCC4_rs2274407_MA", "CYP2A7_rs4079366_MD", "GSTM1_rs1065411_MD",
    "GSTP1_rs4891_MR", "SLCO6A1_rs10055840_MA")
}

fct_mod_var_doxo_sev <- function(){
  c("ABCC1_rs2230671_MR", "ABCC4_rs1751034_MR", "ABCC6_rs12931472_MR", "ABCC6_rs2238472_MR",
    "ABCC6_rs8058694_MR", "ABCC6_rs8058696_MR", "CYP2A7_rs4142867_MD", "GSTM1_rs1056806_MD",
    "GSTM1_rs147668562_MD", "MTHFR_rs1801133_MR", "SLC19A1_chr2145537775_MD",
    "SLCO6A1_rs6884141_MD")
}

fct_mod_var_doxo_pres <- function(){
  c("ABCC1_rs35587_MD", "CYP2A7_rs4079366_MD", "CYP2A7_rs4142867_MD", "MTHFR_rs1801133_MR")
}

## MTX
fct_mod_var_mtx_ulc <- function(){
  c("ABCC1_rs35587_MR", "ABCC2_rs2273697_MD", "ABCC2_rs3740066_MR", "ABCC3_rs1051640_MD",
    "ABCG2_rs2231137_MD", "GSTM1_rs1056806_MR", "SLCO6A1_rs6884141_MR", "TPRA1_chr3127579846_MD")
}

fct_mod_var_mtx_sev <- function(){
  c()
}

fct_mod_var_mtx_pres <- function(){
  c("ABCC1_rs35587_MR", "ABCC2_rs2273697_MD", "ABCC2_rs3740066_MR", "ABCC3_rs1051640_MD" ,
    "ABCC3_rs2277624_MD", "ABCC4_rs2274405_MA", "ABCC4_rs2274406_MA", "ABCC6_rs12931472_MA",
    "CYP2A7_rs4142867_MD", "GSTM1_rs1056806_MR", "GSTM1_rs147668562_MD", "GSTP1_rs4891_MD", 
    "MTHFR_rs4846051_MD", "SLC19A1_rs12659_MD", "SLCO6A1_rs10055840_MD", "SLCO6A1_rs6884141_MA",
    "TPRA1_chr3127579846_MD"
  )
}

## Outros tratamentos
fct_mod_var_out_ulc <- function(){
  c("ABCB1_rs1128503_MA", "GSTA1_rs1051775_MD", "GSTP1_rs1695_MA")
}

fct_mod_var_out_sev <- function(){
  c("CYP2A7_rs117539170_MA", "GSTM1_rs1065411_MR")
}

fct_mod_var_out_pres <- function(){
  c("ABCB1_rs1128503_MR", "ABCC3_rs1051640_MA", "ABCC4_rs2274407_MD", "ABCC6_rs2238472_MD",
    "CYP3A7CYP3A51P_chr799713534_MR", "SLCO6A1_rs10055840_MR")
}

fct_mod_var_out_df <- function(df_protocol, vars_selected, agroup, p_value){
  teste_interno <- F
  if(teste_interno){
    df_protocol <- df_out
    vars_selected <- out_ulc
    agroup = 0
    p_value = 0.05
  }
  ### 1 - Variantes MU selecionadas por p-valor ----
  #### DF agrupado para nao severo (PIORMB = 0, 1, 2) ou severo (PIORMB = 3, 4)
  df_protocol_MU <- df_protocol |> 
    dplyr::select(PIORMB, dplyr::ends_with("MU")) |> 
    dplyr::mutate(PIORMB = ifelse(PIORMB > agroup, 1, 0))
  
  ####  Rodando testes chi-2 e de fisher (não severo vs severo)
  list_testes <- fct_testes_chi2_fisher(df_protocol_MU)
  
  if(list_testes[[1]]){
    ## Testes com e sem simulação de Monte Carlo
    df_chi2_MU <- list_testes[[2]]
  }
  
  #### Variáveis selecionadas MU
  vars_MU_selected <- df_chi2_MU |> 
    dplyr::filter(`p-value(Chi-2)-MC` < p_value) |> 
    dplyr::select(variant) |> 
    dplyr::pull()
  
  ### 2 - Seleção de variantes e agrupamento (presença ausencia nesse caso) ----
  ### Selcionando variáveis MU com p-valor definido e as variáveis recebidas
  df_protocol_agr <-  df_protocol |> 
    dplyr::select(PIORMB, dplyr::matches(vars_selected), dplyr::matches(vars_MU_selected)) |>  
    dplyr::mutate(PIORMB = ifelse(PIORMB > agroup, 1, 0))
  
  ## Alterando valores para factors
  df_protocol_agr <-  df_protocol_agr |> 
    dplyr::mutate(dplyr::across(where(~length(unique(.)) > 1),
                                factor,
                                levels = c(0, 1, 2)))
  df_protocol_agr
}

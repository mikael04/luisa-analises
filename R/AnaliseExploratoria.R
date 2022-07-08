######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_calc_perc.R")
source("R/fct_aux/func_create_table.R")
source("R/fct_aux/func_remove_columns.R")
source("R/fct_aux/func_tables_aux.R")
source("R/fct_aux/func_test_assump.R")
source("R/fct_aux/func_testes_chi2_fisher.R")
source("R/fct_aux/func_writ_organ_xlsx.R")

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav") |> 
  dplyr::filter(!is.na(PIORMB))

set.seed(42)

# 2 - Rodando testes ----

## 2.1 - Teste 1 (Presença) ----
## Primeiro será rodado para agrupamento não possui VS possui
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

df_abs_or_pres <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD"))) |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0))

## Rodando teste
list_testes <- fct_testes_chi2_fisher(df_abs_or_pres)

if(list_testes[[1]]){
  ## Testes sem simulação de Monte Carlo
  df_chi2 <- list_testes[[2]]
  df_fisher <- list_testes[[3]]
}

## Verificando pressupostos (n < 5 em alguma categoria; tabela com mais de uma coluna)
df_checks <- fct_test_assump(df_abs_or_pres)

df_chi2_fisher <- dplyr::inner_join(df_chi2, df_fisher, by = "variant")
df_chi2_fisher_checks <- dplyr::inner_join(df_chi2_fisher, df_checks, by = "variant")

p_value <- 0.05

df_chi2_fisher_sig <- df_chi2_fisher_checks |> 
  dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
                `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value)

# 3 - Avaliando modelos (dominante, recessivo, aditivo) ----
df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])

## 3.1 - Gerando tabelas para modelo presença vs ausência ----

type_group = "pres"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2, table, type_group, write_table, switch_teste)
if(return_create_table == 1){
  print("Geração de todas as tabelas para parâmetro presença vs ausência feitas com sucesso")
}
if(return_create_table > 1 & return_create_table <= 5){
  print("Geração de tabela individual feita com sucesso")
  print(paste0("Gerada apenas tabela para modelo ", table))
}

## 3.2 - Gerando tabelas para modelo não ulcerados vs ulcerados ----

type_group = "ulc"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2, table, type_group, write_table, switch_teste)
if(return_create_table == 1){
  print("Geração de todas as tabelas para parametro ulcerados vs não ulcerados feitas com sucesso")
}
if(return_create_table > 1 & return_create_table <= 5){
  print("Geração de tabela individual feita com sucesso")
  print(paste0("Gerada apenas tabela para modelo ", table))
}


## 3.3 - Gerando tabelas para modelo mb severo vs não severo ----

type_group = "sev"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2, table, type_group, write_table, switch_teste)
if(return_create_table == 1){
  print("Geração de todas as tabelas para parametro mb severo vs não severo feitas com sucesso")
}
if(return_create_table > 1 & return_create_table <= 5){
  print("Geração de tabela individual feita com sucesso")
  print(paste0("Gerada apenas tabela para modelo ", table))
}


#### Avaliando modelo
# summary()

#### Resultado com multicolinearidade, tratar multicolinearidade
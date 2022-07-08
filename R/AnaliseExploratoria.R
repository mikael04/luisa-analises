######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_calc_perc.R")
source("R/fct_aux/func_check_return.R")
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
## Teste 1 será rodado para agrupamento não possui (ausencia) VS possui (presenca) de mucosite bocal
## ausencia (PIORMB = 0) vs presença (PIORMB = 1, 2, 3, 4)
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para ausência (PIORMB = 0) ou presença (PIORMB = 1, 2, 3, 4)
df_abs_or_pres <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD"))) |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0))

### 2.1.1 Rodando testes chi-2 e de fisher (ausencia vs presença) ----
list_testes <- fct_testes_chi2_fisher(df_abs_or_pres)

if(list_testes[[1]]){
  ## Testes com e sem simulação de Monte Carlo
  df_chi2_a_p <- list_testes[[2]]
  df_fisher_a_p <- list_testes[[3]]
}

### 2.1.2 - Avaliando modelos (dominante, recessivo, aditivo) ----
df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])

## 2.1.2.1 - Gerando tabelas para modelos agrupados por ausência vs presença ----

type_group = "pres"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_a_p, table, type_group, write_table, switch_teste)
fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)

## 2.2 - Teste 2 (Ulcerações) ----
## Teste 2 será rodado para agrupamento não ulcerados vs ulcerados
## não ulcerado (PIORMB = 0, 1) ou ulcerado (PIORMB = 2, 3, 4)
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para nao ulcerado (PIORMB = 0, 1) ou ulcerado (PIORMB = 2, 3, 4)
df_ulc <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD"))) |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 1, 1, 0))

### 2.2.1 Rodando testes chi-2 e de fisher (não ulcerados vs ulcerados) ----
list_testes <- fct_testes_chi2_fisher(df_ulc)

if(list_testes[[1]]){
  ## Testes com e sem simulação de Monte Carlo
  df_chi2_u <- list_testes[[2]]
  df_fisher_u <- list_testes[[3]]
}

## Gerando tabelas para modelo não ulcerados vs ulcerados ----

type_group = "ulc"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_u, table, type_group, write_table, switch_teste)
fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)


## 2.3 - Teste 3 (Severidade) ----
## Teste 2 será rodado para agrupamento não severo vs severo
## não severo (PIORMB = 0, 1, 2) ou severo (PIORMB = 3, 4)
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para nao severo (PIORMB = 0, 1, 2) ou severo (PIORMB = 3, 4)
df_sev <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD"))) |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 2, 1, 0))

### 2.2.1 Rodando testes chi-2 e de fisher (não ulcerados vs ulcerados) ----
list_testes <- fct_testes_chi2_fisher(df_sev)

if(list_testes[[1]]){
  ## Testes com e sem simulação de Monte Carlo
  df_chi2_s <- list_testes[[2]]
  df_fisher_s <- list_testes[[3]]
}

## 3.3 - Gerando tabelas para modelo mb severo vs não severo ----

type_group = "sev"
table <- NULL
switch_teste <- T

return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_s, table, type_group, write_table, switch_teste)
fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)


#### Avaliando modelo
# summary()

#### Resultado com multicolinearidade, tratar multicolinearidade
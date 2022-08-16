######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 02/08/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_qui2_fisher.R")


## 0.1 Parâmetros globais ----

## P-value padrão
p_value <- 0.05
## Se usaremos Monte Carlo ou não no teste qui-quadrado
use_MC <- T

options(dplyr.summarise.inform = FALSE)
switch_teste <- F
switch_write_table <- T
switch_write_binom <- F
switch_overwrite_binom <- F


# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav") |> 
  dplyr::filter(!is.na(PIORMB))

set.seed(42)

# 2 - Qui-quadrado e teste de fisher ----

## 2.1 - DOXO ----
## Teste 1 será rodado para os agrupamentos: 
## 1 - não possui (ausencia) VS possui (presenca) de mucosite bocal
## 2 - não ulcerados vs ulcerados
## 3 - não severo vs severo
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para DOXO
df_doxo <- df_full |>
  dplyr::filter(Protocolo == 2)

path_tables <- "data-raw/tabelas_modelos_protocolos/doxo/"
if(!dir.exists(path_tables)){
  dir.create(path_tables)
}

tryCatch({
  message("Iniciando qui-2 e fisher para protocolo DOXO dominante")
  fct_qui2_fisher(df_doxo, path_tables, debug)
  },
  error=function(cond) {
    message("Erro ao executar a função")
    message("Mensagem original:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  finally={
    message("")
    message("Finalizando execução para protocolo DOXO dominante")
  }
)

## 2.2 - CTX ----
## Teste 2 será rodado para os agrupamentos: 
## 1 - não possui (ausencia) VS possui (presenca) de mucosite bocal
## 2 - não ulcerados vs ulcerados
## 3 - não severo vs severo
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para CTX
df_ctx <- df_full |>
  dplyr::filter(Protocolo == 6)

path_tables <- "data-raw/tabelas_modelos_protocolos/ctx/"
if(!dir.exists(path_tables)){
  dir.create(path_tables)
}

tryCatch({
    message("Iniciando qui-2 e fisher para protocolo CTX predomina")
    fct_qui2_fisher(df_ctx, path_tables, debug)
  },
  error=function(cond) {
    message("Erro ao executar a função")
    message("Mensagem original:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  finally={
    message("")
    message("Finalizando execução para protocolo CTX predomina")
  }
)

## 2.3 - MTX ----
## Teste 3 será rodado para os agrupamentos: 
## 1 - não possui (ausencia) VS possui (presenca) de mucosite bocal
## 2 - não ulcerados vs ulcerados
## 3 - não severo vs severo
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para MTX
df_mtx <- df_full |>
  dplyr::filter(Protocolo == 1)

path_tables <- "data-raw/tabelas_modelos_protocolos/mtx/"
if(!dir.exists(path_tables)){
  dir.create(path_tables)
}

tryCatch({
    message("Iniciando qui-2 e fisher para protocolo MTX altas doses")
    fct_qui2_fisher(df_mtx, path_tables, debug)
  },
  error=function(cond) {
    message("Erro ao executar a função")
    message("Mensagem original:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  finally={
    message("")
    message("Finalizando execução para protocolo MTX altas doses")
  }
)

## 2.4 - Outros agrupamentos ----
## Teste 4 será rodado para os agrupamentos: 
## 1 - não possui (ausencia) VS possui (presenca) de mucosite bocal
## 2 - não ulcerados vs ulcerados
## 3 - não severo vs severo
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para outros
df_others <- df_full |>
  dplyr::filter(Protocolo %in% c(3,4,5,7,8))

path_tables <- "data-raw/tabelas_modelos_protocolos/outros/"
if(!dir.exists(path_tables)){
  dir.create(path_tables)
}

tryCatch({
    message("Iniciando qui-2 e fisher para protocolo outros agrupamentos de protocolos")
    fct_qui2_fisher(df_others, path_tables, debug)
  },
  error=function(cond) {
    message("Erro ao executar a função")
    message("Mensagem original:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  finally={
    message("")
    message("Finalizando execução para protocolo outros agrupamentos de protocolos")
  }
)

# 
# 1 - MTX altas doses
# 2 - DOXO predomina
# 3 - MTX + CTX + DOXO
# 4 - MTX baixa doses
# 5 - MTX + CTX
# 6 - CTX predomiona
# 7 - ARA C etoposide, carbo, outros
# 8 - CTX + DOXO
# 
# Grupos:
# 1 - DOXO predomina (2)
# 2 - CTX predomina (6)
# 3 - MTX altas doses (1)
# 4 - Outros (3, 4, 5, 7, 8)
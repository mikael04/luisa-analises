######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 02/08/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_qui2_fisher.R")
source("R/fct_aux/func_model_variant_selection.R")
source("R/fct_aux/func_print_glm_xlsx.R")


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

## 2.1 - CTX ----
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

## 2.2 - DOXO ----
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
df_out <- df_full |>
  dplyr::filter(Protocolo %in% c(3,4,5,7,8))

path_tables <- "data-raw/tabelas_modelos_protocolos/outros/"
if(!dir.exists(path_tables)){
  dir.create(path_tables)
}

tryCatch({
    message("Iniciando qui-2 e fisher para protocolo outros agrupamentos de protocolos")
    fct_qui2_fisher(df_out, path_tables, debug)
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

# 3. Regressão binomial multivariada ----
print("Iniciando regressão binomial multivariada")
## 3.1 - CTX ----

### 3.1.1 - Presença ou ausência ----
#### 3.1.1.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
ctx_pres_aus <- fct_mod_var_ctx_pres()

## Agroup = 0 -> Presença ausência
agroup = 0

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(ctx_pres_aus)){
  df_ctx_pres_aus <- fct_mod_var_out_df(df_ctx, ctx_pres_aus, agroup, p_value)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  print(summary(ctx_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_ctx_pres_aus, family = "binomial")))
  
  # hnp::hnp(ctx_abs_or_pres_binom_mult, resid.type = "deviance")
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC3_rs11568591_MU + ABCC6_rs9940825_MD + HSP90AA1_rs4947_MA +
      HSP90AA1_rs8005905_MA + SLC19A1_rs1051266_MR,
    data = df_ctx_pres_aus, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC3_rs11568591_MU + ABCC6_rs9940825_MD + HSP90AA1_rs4947_MA,
    data = df_ctx_pres_aus, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_abs_or_pres_binom_mult, switch_write_table)
  }
    
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo CTX")
}

### 3.1.2 - Severidade ----
#### 3.1.2.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
ctx_sev <- fct_mod_var_ctx_sev()

## Agroup = 2 -> Severidade
agroup = 2

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(ctx_sev)){
  df_ctx_sev <- fct_mod_var_out_df(df_ctx, ctx_sev, agroup, p_value)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  summary(ctx_sev_binom_mult <- glm(PIORMB ~ ., data = df_ctx_sev, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_ctx_sev, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_ctx_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_sev_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo CTX")
}

### 3.1.3 - Ulcerações ----
#### 3.1.3.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
ctx_ulc <- fct_mod_var_ctx_ulc()

## Agroup = 1 -> Ulcerações
agroup = 1

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(ctx_ulc)){
  df_ctx_ulc <- fct_mod_var_out_df(df_ctx, ctx_ulc, agroup, p_value)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  print(summary(ctx_ulc_binom_mult <- glm(PIORMB ~ ., data = df_ctx_ulc, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_ulc_binom_mult <- glm(
    PIORMB ~ ABCC6_rs12931472_MD + ABCC6_rs72657698_MU +
      GSTA1_rs1051775_MA + HSP90AA1_rs4947_MA,
    data = df_ctx_ulc, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_ulc_binom_mult <- glm(
    PIORMB ~ GSTA1_rs1051775_MA + HSP90AA1_rs4947_MA,
    data = df_ctx_ulc, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_ulc_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento ulcerações, protocolo CTX")
}
## 3.2 - DOXO ----

### 3.2.1 - Presença ou ausência ----
#### 3.2.1.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
doxo_pres_aus <- fct_mod_var_doxo_pres()

## Agroup = 0 -> Presença ausência
agroup = 0

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(doxo_pres_aus)){
  df_doxo_pres_aus <- fct_mod_var_out_df(df_doxo, doxo_pres_aus, agroup, p_value)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  print(summary(doxo_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_doxo_pres_aus, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MD + ABCC3_chr1750683660_MU + 
      CYP2A7_rs4079366_MD + MTHFR_rs1801133_MR,
    data = df_doxo_pres_aus, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MD + CYP2A7_rs4079366_MD + MTHFR_rs1801133_MR,
    data = df_doxo_pres_aus, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_abs_or_pres_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo DOXO")
}

### 3.2.2 - Severidade ----
#### 3.2.2.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
doxo_sev <- fct_mod_var_doxo_sev()

## Agroup = 2 -> Severidade
agroup = 2

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(doxo_sev)){
  df_doxo_sev <- fct_mod_var_out_df(df_doxo, doxo_sev, agroup, p_value)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  summary(doxo_sev_binom_mult <- glm(PIORMB ~ ., data = df_doxo_sev, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_sev_binom_mult <- glm(
    PIORMB ~ ABCC4_rs1751034_MR + CYP2A7_rs4142867_MD + SLC31A1_chr9113258719_MU,
    data = df_doxo_sev, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_sev_binom_mult <- glm(
    PIORMB ~ ABCC4_rs1751034_MR + CYP2A7_rs4142867_MD + SLC31A1_chr9113258719_MU,
    data = df_doxo_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_sev_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo DOXO")
}

### 3.2.3 - Ulcerações ----
#### 3.2.3.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
doxo_ulc <- fct_mod_var_doxo_ulc()

## Agroup = 1 -> Ulcerações
agroup = 1

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(doxo_ulc)){
  df_doxo_ulc <- fct_mod_var_out_df(df_doxo, doxo_ulc, agroup, p_value)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  print(summary(doxo_ulc_binom_mult <- glm(PIORMB ~ ., data = df_doxo_ulc, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD + ABCC4_rs2274407_MA + GSTM1_rs1065411_MD + 
      GSTP1_rs4891_MR,
    data = df_doxo_ulc, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD + GSTM1_rs1065411_MD + GSTP1_rs4891_MR,
    data = df_doxo_ulc, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_sev_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento ulcerações, protocolo DOXO")
}

## 3.3 - MTX ----

### 3.3.1 - Presença ou ausência ----
#### 3.3.1.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
mtx_pres_aus <- fct_mod_var_mtx_pres()

## Agroup = 0 -> Presença ausência
agroup = 0

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(mtx_pres_aus)){
  df_mtx_pres_aus <- fct_mod_var_out_df(df_mtx, mtx_pres_aus, agroup, p_value)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  print(summary(mtx_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_mtx_pres_aus, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(mtx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MR + ABCC2_rs2273697_MD + 
      ABCC2_rs3740066_MR + ABCC3_rs1051640_MD + ABCC3_rs2277624_MD + 
      ABCC4_rs2274405_MA + ABCC4_rs2274406_MA + ABCC6_rs12931472_MA + 
      CYP2A7_rs4142867_MD + GSTM1_rs1056806_MR + GSTM1_rs147668562_MD + 
      GSTP1_rs4891_MD + MTHFR_rs4846051_MD + SLC19A1_rs12659_MD + 
      SLCO6A1_rs10055840_MD + SLCO6A1_rs6884141_MA + TPRA1_chr3127579846_MD + 
      ABCA3_rs1319979593_MU + ABCA3_rs149532_MU + ABCC2_rs1137968_MU + 
      ABCC2_rs17222723_MU + ABCC2_rs8187707_MU + ABCC2_rs8187710_MU + 
      ABCC3_chr1750683660_MU + ABCC3_rs11568591_MU + CCND1_rs1181031465_MU + 
      CYP2A6_chr1940848742_MU,
    data = df_mtx_pres_aus, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(mtx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MR + ABCC2_rs2273697_MD + 
      ABCC2_rs3740066_MR + ABCC3_rs1051640_MD + ABCC3_rs2277624_MD + 
      ABCC4_rs2274405_MA + ABCC4_rs2274406_MA + ABCC6_rs12931472_MA + 
      CYP2A7_rs4142867_MD + GSTM1_rs1056806_MR + GSTM1_rs147668562_MD + 
      GSTP1_rs4891_MD + SLC19A1_rs12659_MD + 
      SLCO6A1_rs10055840_MD + SLCO6A1_rs6884141_MA + TPRA1_chr3127579846_MD + 
      ABCA3_rs1319979593_MU + ABCA3_rs149532_MU
      ,
    data = df_mtx_pres_aus, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_abs_or_pres_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo MTX")
}

### 3.3.2 - Severidade ----
#### 3.3.2.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
mtx_sev <- fct_mod_var_mtx_sev()

## Agroup = 2 -> Severidade
agroup = 2

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(mtx_sev)){
  df_mtx_sev <- fct_mod_var_out_df(df_mtx, mtx_sev, agroup, p_value)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  summary(mtx_sev_binom_mult <- glm(PIORMB ~ ., data = df_mtx_sev, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(mtx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_mtx_sev, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(mtx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_mtx_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_sev_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo MTX")
}

### 3.3.3 - Ulcerações ----
#### 3.3.3.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
mtx_ulc <- fct_mod_var_mtx_ulc()

## Agroup = 1 -> Ulcerações
agroup = 1

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(mtx_ulc)){
  df_mtx_ulc <- fct_mod_var_out_df(df_mtx, mtx_ulc, agroup, p_value)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  print(summary(mtx_ulc_binom_mult <- glm(PIORMB ~ ., data = df_mtx_ulc, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(mtx_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD + GSTM1_rs1056806_MR + 
      SLCO6A1_rs6884141_MR + ABCC2_rs1137968_MU,
    data = df_mtx_ulc, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(mtx_ulc_binom_mult <- glm(
    PIORMB ~ SLCO6A1_rs6884141_MR,
    data = df_mtx_ulc, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_ulc_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento ulcerações, protocolo MTX")
}

## 3.4 - Outros ----

### 3.4.1 - Presença ou ausência ----
#### 3.4.1.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
out_pres_aus <- fct_mod_var_out_pres()

## Agroup = 0 -> Presença ausência
agroup = 0

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(out_pres_aus)){
  df_out_pres_aus <- fct_mod_var_out_df(df_out, out_pres_aus, agroup, p_value)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  print(summary(out_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_out_pres_aus, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCB1_rs1128503_MR + ABCC3_rs1051640_MA + ABCC4_rs2274407_MD + 
      ABCC6_rs2238472_MD + SLCO6A1_rs10055840_MR + ABCA3_rs1319979593_MU,
    data = df_out_pres_aus, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCB1_rs1128503_MR + ABCC3_rs1051640_MA + ABCC4_rs2274407_MD + 
      ABCA3_rs1319979593_MU,
    data = df_out_pres_aus, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_abs_or_pres_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo Outros")
}

### 3.4.2 - Severidade ----
#### 3.4.2.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
out_sev <- fct_mod_var_out_sev()

## Agroup = 2 -> Severidade
agroup = 2

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(out_sev)){
  df_out_sev <- fct_mod_var_out_df(df_out, out_sev, agroup, p_value)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  summary(out_sev_binom_mult <- glm(PIORMB ~ ., data = df_out_sev, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_sev_binom_mult <- glm(
    PIORMB ~ CYP2A7_rs117539170_MA + GSTM1_rs1065411_MR + MTHFR_rs2066470_MU,
    data = df_out_sev, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_sev_binom_mult <- glm(
    PIORMB ~ GSTM1_rs1065411_MR + MTHFR_rs2066470_MU,
    data = df_out_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_sev_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento presença ou ausência, protocolo Outros")
}

### 3.4.3 - Ulcerações ----
#### 3.4.3.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
out_ulc <- fct_mod_var_out_ulc()

## Agroup = 1 -> Ulcerações
agroup = 1

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(out_ulc)){
  df_out_ulc <- fct_mod_var_out_df(df_out, out_ulc, agroup, p_value)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  print(summary(out_ulc_binom_mult <- glm(PIORMB ~ ., data = df_out_ulc, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_ulc_binom_mult <- glm(
    PIORMB ~ CYP2A7_rs117539170_MA + GSTM1_rs1065411_MR + MTHFR_rs2066470_MU,
    data = df_out_sev, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_ulc_binom_mult <- glm(
    PIORMB ~ GSTM1_rs1065411_MR + MTHFR_rs2066470_MU,
    data = df_out_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_ulc_binom_mult, switch_write_table)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento ulcerações, protocolo Outros")
}

# 4. Criando tabelas com análise da regressão binomial múltipla ----
## 4.1 - CTX ----

### 4.1.1 - Presença ou ausência ----
#### 4.1.1.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
ctx_pres_aus <- fct_mod_var_ctx_pres()
df_

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
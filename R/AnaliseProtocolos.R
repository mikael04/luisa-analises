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
source("R/fct_aux/func_res_bin_mult.R")
source("R/fct_aux/func_get_dummy_cols.R")


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
switch_test_chi2_protocol <- F

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
if(switch_test_chi2_protocol){
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
}else{
  print("Não vai rodar o teste qui-2 para protocolo CTX")
}

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

if(switch_test_chi2_protocol){
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
}else{
  print("Não vai rodar o teste qui-2 para protocolo DOXO")
}

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

if(switch_test_chi2_protocol){
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
}else{
  print("Não vai rodar o teste qui-2 para protocolo MTX")
}

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

if(switch_test_chi2_protocol){
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
}else{
  print("Não vai rodar o teste qui-2 para protocolo 'outros protocolos'")
}

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
  df_ctx_pres_aus <- fct_mod_var_out_df(df_ctx, ctx_pres_aus, agroup, p_value = 0.05)
  
  df_ctx_pres_aus_fast_dummies <-  fct_get_dummy_cols(df_ctx_pres_aus, debug)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  print(summary(ctx_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_ctx_pres_aus_fast_dummies, family = "binomial")))
  
  # hnp::hnp(ctx_abs_or_pres_binom_mult, resid.type = "deviance")
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC6_rs9940825_MD_1 + HSP90AA1_rs4947_MA_1 + 
      HSP90AA1_rs8005905_MA_1 + SLC19A1_rs12659_MR_1 + ABCC3_rs11568591_MU_1,
    data = df_ctx_pres_aus_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC6_rs9940825_MD_1 + HSP90AA1_rs4947_MA_1 + 
      SLC19A1_rs12659_MR_1,
    data = df_ctx_pres_aus_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_abs_or_pres_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(ctx_abs_or_pres_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(ctx_abs_or_pres_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_ctx, vars_aux, p_value_vars, "ctx_abs_or_pres_binom_mult", 
                     switch_write_table = T)
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
  df_ctx_sev <- fct_mod_var_out_df(df_ctx, ctx_sev, agroup, p_value = 0.05)
  
  df_ctx_sev_fast_dummies <-  fct_get_dummy_cols(df_ctx_sev, debug)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  summary(ctx_sev_binom_mult <- glm(PIORMB ~ ., data = df_ctx_sev_fast_dummies, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_ctx_sev_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_ctx_sev, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_sev_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(ctx_sev_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(ctx_sev_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_ctx, vars_aux, p_value_vars, "ctx_sev_binom_mult", 
                     switch_write_table = T)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento severidade, protocolo CTX")
}

### 3.1.3 - Ulcerações ----
#### 3.1.3.1 - Variantes selecionadas pela cliente (MA, MD, MR) ----
ctx_ulc <- fct_mod_var_ctx_ulc()

## Agroup = 1 -> Ulcerações
agroup = 1

## Check para ver se existe alguma variante na seleção, senão existir, o modelo não será rodado
if(!is.null(ctx_ulc)){
  df_ctx_ulc <- fct_mod_var_out_df(df_ctx, ctx_ulc, agroup, p_value = 0.05)
  ## Gerando dummies
  df_ctx_ulc_fast_dummies <-  fct_get_dummy_cols(df_ctx_ulc, debug)
  #### 3.1.1.2 Modelo com variantes selecionadas ----
  print(summary(ctx_ulc_binom_mult <- glm(PIORMB ~ ., data = df_ctx_ulc_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(ctx_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(ctx_ulc_binom_mult <- glm(
    PIORMB ~ ABCC6_rs12931472_MD_1 + GSTA1_rs1051775_MA_1 + 
      HSP90AA1_rs10873531_MA_1 + SLCO6A1_rs6884141_MR_1,
    data = df_ctx_ulc_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(ctx_ulc_binom_mult <- glm(
    PIORMB ~ HSP90AA1_rs10873531_MA_1,
    data = df_ctx_ulc_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(ctx_ulc_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(ctx_ulc_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(ctx_ulc_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_ctx, vars_aux, p_value_vars, "ctx_ulc_binom_mult", 
                     switch_write_table = T)
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
  df_doxo_pres_aus <- fct_mod_var_out_df(df_doxo, doxo_pres_aus, agroup, p_value = 0.05)
  ## Gerando dummies
  df_doxo_pres_aus_fast_dummies <-  fct_get_dummy_cols(df_doxo_pres_aus, debug)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  print(summary(doxo_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_doxo_pres_aus_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MD_1 + CYP2A7_rs4079366_MD_1 + 
      MTHFR_rs1801133_MR_1,
    data = df_doxo_pres_aus_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MD_1 + CYP2A7_rs4079366_MD_1 + 
      MTHFR_rs1801133_MR_1,
    data = df_doxo_pres_aus_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_abs_or_pres_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(doxo_abs_or_pres_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(doxo_abs_or_pres_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_doxo, vars_aux, p_value_vars, "doxo_abs_or_pres_binom_mult", 
                     switch_write_table = T)
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
  df_doxo_sev <- fct_mod_var_out_df(df_doxo, doxo_sev, agroup, p_value = 0.05)
  ## Gerando dummies
  df_doxo_sev_fast_dummies <-  fct_get_dummy_cols(df_doxo_sev, debug)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  summary(doxo_sev_binom_mult <- glm(PIORMB ~ ., data = df_doxo_sev_fast_dummies, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_sev_binom_mult <- glm(
    PIORMB ~ ABCC4_rs1751034_MR_1 + CYP2A7_rs4142867_MD_1 + 
      SLC31A1_chr9113258719_MU_1,
    data = df_doxo_sev_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_sev_binom_mult <- glm(
    PIORMB ~ ABCC4_rs1751034_MR_1,
    data = df_doxo_sev_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_sev_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(doxo_sev_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(doxo_sev_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_doxo, vars_aux, p_value_vars, "doxo_sev_binom_mult", 
                     switch_write_table = T)
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
  df_doxo_ulc <- fct_mod_var_out_df(df_doxo, doxo_ulc, agroup, p_value = 0.05)
  ## Gerando dummies
  df_doxo_ulc_fast_dummies <-  fct_get_dummy_cols(df_doxo_ulc, debug)
  #### 3.2.1.2 Modelo com variantes selecionadas ----
  print(summary(doxo_ulc_binom_mult <- glm(PIORMB ~ ., data = df_doxo_ulc_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(doxo_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(doxo_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD_1 + GSTM1_rs1065411_MD_1 + 
      GSTP1_rs4891_MR_1 + ABCC4_rs2274407_MA_2,
    data = df_doxo_ulc_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(doxo_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD_1 + GSTM1_rs1065411_MD_1 + 
      GSTP1_rs4891_MR_1,
    data = df_doxo_ulc_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(doxo_ulc_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(doxo_ulc_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(doxo_ulc_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_doxo, vars_aux, p_value_vars, "doxo_ulc_binom_mult", 
                     switch_write_table = T)
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
  df_mtx_pres_aus <- fct_mod_var_out_df(df_mtx, mtx_pres_aus, agroup, p_value = 0.05)
  ## Gerando dummies
  df_mtx_pres_aus_fast_dummies <-  fct_get_dummy_cols(df_mtx_pres_aus, debug)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  print(summary(mtx_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_mtx_pres_aus_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step ANTEPENULTIMA
  print(summary(mtx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC1_rs35587_MR_1 + ABCC2_rs2273697_MD_1 + 
      ABCC2_rs3740066_MR_1 + ABCC3_rs1051640_MD_1 + ABCC3_rs2277624_MD_1 + 
      ABCC4_rs2274405_MA_1 + ABCC4_rs2274406_MA_1 + ABCC6_rs12931472_MA_1 + 
      CYP2A7_rs4142867_MD_1 + GSTM1_rs1056806_MR_1 + GSTM1_rs147668562_MD_1 + 
      GSTP1_rs4891_MD_1 + MTHFR_rs4846051_MD_1 + SLC19A1_rs12659_MD_1 + 
      SLCO6A1_rs10055840_MD_1 + SLCO6A1_rs6884141_MA_1 + ABCA3_rs149532_MU_1 + ABCC2_rs17222723_MU_1 + 
      ABCC2_rs8187707_MU_1 + ABCC2_rs8187710_MU_1 + ABCC3_chr1750683660_MU_1 + 
      ABCC3_rs11568591_MU_1 + CCND1_rs1181031465_MU_1
    ## removidos a pedido
    ## ABCA3_rs1319979593_MU_1 + TPRA1_chr3127579846_MD_1
    ,
    data = df_mtx_pres_aus_fast_dummies, family = "binomial")))
  
  # ## Modelo apenas com variantes significativas (a partir do modelo step) (Removendo apenas ABCA3)
  # 
  # print(summary(mtx_abs_or_pres_binom_mult <- glm(
  #   PIORMB ~ ABCC1_rs35587_MR_1 + TPRA1_chr3127579846_MD_1 +
  #     ABCA3_rs149532_MU_1 + ABCC3_chr1750683660_MU_1,
  #   data = df_mtx_pres_aus_fast_dummies, family = "binomial")))
  # with(summary(mtx_abs_or_pres_binom_mult), 1 - deviance/null.deviance)
  
  ## Modelo com variantes significativas (a partir do modelo step) (Removendo ABCA3 e TPRA1)
  
  print(summary(mtx_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD_1 + ABCC4_rs2274406_MA_1 +  
      GSTM1_rs1056806_MR_1 + ABCC2_rs17222723_MU_1,
    data = df_mtx_pres_aus_fast_dummies, family = "binomial")))
  
  # with(summary(mtx_abs_or_pres_binom_mult), 1 - deviance/null.deviance)
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_abs_or_pres_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(mtx_abs_or_pres_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(mtx_abs_or_pres_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_mtx, vars_aux, p_value_vars, "mtx_abs_or_pres_binom_mult", 
                     switch_write_table = T)
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
  df_mtx_sev <- fct_mod_var_out_df(df_mtx, mtx_sev, agroup, p_value = 0.05)
  ## Gerando dummies
  df_mtx_sev_fast_dummies <-  fct_get_dummy_cols(df_mtx_sev, debug)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  summary(mtx_sev_binom_mult <- glm(PIORMB ~ ., data = df_mtx_sev_fast_dummies, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(mtx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_mtx_sev_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(mtx_sev_binom_mult <- glm(
    PIORMB ~ .,
    data = df_mtx_sev_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_sev_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(mtx_sev_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(mtx_sev_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_mtx, vars_aux, p_value_vars, "mtx_sev_binom_mult", 
                     switch_write_table = T)
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
  df_mtx_ulc <- fct_mod_var_out_df(df_mtx, mtx_ulc, agroup, p_value = 0.05)
  ## Gerando dummies
  df_mtx_ulc_fast_dummies <-  fct_get_dummy_cols(df_mtx_ulc, debug)
  #### 3.3.1.2 Modelo com variantes selecionadas ----
  print(summary(mtx_ulc_binom_mult <- glm(PIORMB ~ ., data = df_mtx_ulc_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(mtx_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(mtx_ulc_binom_mult <- glm(
    PIORMB ~ ABCC2_rs2273697_MD_1 + ABCC3_rs1051640_MD_1 + 
      GSTM1_rs1056806_MR_1 + SLCO6A1_rs6884141_MR_1,
    data = df_mtx_ulc_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(mtx_ulc_binom_mult <- glm(
    PIORMB ~ ABCC3_rs1051640_MD_1 + SLCO6A1_rs6884141_MR_1,
    data = df_mtx_ulc_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(mtx_ulc_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(mtx_ulc_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(mtx_ulc_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_mtx, vars_aux, p_value_vars, "mtx_ulc_binom_mult", 
                     switch_write_table = T)
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
  df_out_pres_aus <- fct_mod_var_out_df(df_out, out_pres_aus, agroup, p_value = 0.05)
  ## Gerando dummies
  df_out_pres_aus_fast_dummies <-  fct_get_dummy_cols(df_out_pres_aus, debug)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  print(summary(out_abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_out_pres_aus_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_abs_or_pres_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCB1_rs1128503_MR_1 + ABCC3_rs1051640_MA_1 + 
      ABCC4_rs2274407_MD_1 + ABCC6_rs2238472_MD_1 + SLCO6A1_rs10055840_MR_1,
    data = df_out_pres_aus_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_abs_or_pres_binom_mult <- glm(
    PIORMB ~ ABCB1_rs1128503_MR_1 + ABCC3_rs1051640_MA_1,
    data = df_out_pres_aus_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_abs_or_pres_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(out_abs_or_pres_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(out_abs_or_pres_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_out, vars_aux, p_value_vars, "out_abs_or_pres_binom_mult", 
                     switch_write_table = T)
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
  df_out_sev <- fct_mod_var_out_df(df_out, out_sev, agroup, p_value = 0.05)
  ## Gerando dummies
  df_out_sev_fast_dummies <-  fct_get_dummy_cols(df_out_sev, debug)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  summary(out_sev_binom_mult <- glm(PIORMB ~ ., data = df_out_sev_fast_dummies, family = "binomial"))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_sev_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_sev_binom_mult <- glm(
    PIORMB ~ CYP2A7_rs117539170_MA_1 + GSTM1_rs1065411_MR_1,
    data = df_out_sev_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_sev_binom_mult <- glm(
    PIORMB ~ GSTM1_rs1065411_MR_1,
    data = df_out_sev_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_sev_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(out_sev_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(out_sev_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_out, vars_aux, p_value_vars, "out_sev_binom_mult", 
                     switch_write_table = T)
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
  df_out_ulc <- fct_mod_var_out_df(df_out, out_ulc, agroup, p_value = 0.05)
  ## Gerando dummies
  df_out_ulc_fast_dummies <-  fct_get_dummy_cols(df_out_ulc, debug)
  #### 3.4.1.2 Modelo com variantes selecionadas ----
  print(summary(out_ulc_binom_mult <- glm(PIORMB ~ ., data = df_out_ulc_fast_dummies, family = "binomial")))
  
  ## Não utilizaremos o step, já que as features já foram definidas
  print(step(out_ulc_binom_mult)) # Poderia criar uma função para remover uma variante por vez, mas
  # acho que resultaria nisso de qualquer forma.
  
  ## Modelo com variantes sugeridas pelo step
  print(summary(out_ulc_binom_mult <- glm(
    PIORMB ~ GSTP1_rs1695_MA_1 + ABCB1_rs1128503_MA_2 + 
      GSTP1_rs1695_MA_2,
    data = df_out_ulc_fast_dummies, family = "binomial")))
  
  ## Modelo apenas com variantes significativas (a partir do modelo step)
  print(summary(out_ulc_binom_mult <- glm(
    PIORMB ~ ABCB1_rs1128503_MA_2 + 
      GSTP1_rs1695_MA_2,
    data = df_out_ulc_fast_dummies, family = "binomial")))
  
  if(switch_write_table){
    fct_print_glm_xlsx(out_ulc_binom_mult, switch_write_table)
    
    ## por algum motivo deu pau passando por parâmetro, vou tentar jogar aqui
    vars_aux <- gsub(".{1}$", "", names(out_ulc_binom_mult$qr$qr[2,]))
    p_value_vars <- round(coef(summary(out_ulc_binom_mult))[,"Pr(>|z|)"], 4)
    
    fct_res_bin_mult(df_out, vars_aux, p_value_vars, "out_ulc_binom_mult", 
                     switch_write_table = T)
  }
}else{
  print("Modelo não será gerado, não existem variantes selecionadas para o agrupamento ulcerações, protocolo Outros")
}


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

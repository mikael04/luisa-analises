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
source("R/fct_aux/func_regression_poisson.R")
source("R/fct_aux/func_regression_binomial.R")
source("R/fct_aux/func_select_variant.R")

## 0.1 Parâmetros globais ----

p_value <- 0.05

options(dplyr.summarise.inform = FALSE)
switch_teste <- F
switch_write_table <- F


# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav") |> 
  dplyr::filter(!is.na(PIORMB))

set.seed(42)

# 2 - Qui-quadrado e teste de fisher ----

## 2.1 - Teste 1 (Presença) ----
## Teste 1 será rodado para agrupamento não possui (ausencia) VS possui (presenca) de mucosite bocal
## ausencia (PIORMB = 0) vs presença (PIORMB = 1, 2, 3, 4)
## Ajustaremos a base de dados para o novo agrupamento (utilizando a variavel "PIORMB")

## DF agrupado para ausência (PIORMB = 0) ou presença (PIORMB = 1, 2, 3, 4)
df_abs_or_pres <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD", "MU"))) |> 
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

#### 2.1.2.1 - Gerando tabelas para modelos agrupados por ausência vs presença ----
## Adicionando Tabelas de contingência

type_group = "pres"
table <- NULL

if(switch_write_table){
  return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_a_p, table, type_group, write_table, switch_teste)
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
  return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_u, table, type_group, write_table, switch_teste)
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
  return_create_table <- fct_create_table(df_modelos_genotipos, df_chi2_s, table, type_group, write_table, switch_teste)
  fct_check_return(return_create_table, "create_table", table, type_group, switch_teste)
}

# # 3 Regressão de Poisson Univariada ----
# ## 3.1 Presença vs ausência ----
# 
# abs_or_pres_sig <- dplyr::inner_join(df_chi2_a_p, df_fisher_a_p, by = "variant") |>
#   dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
#                   `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
#   dplyr::select(variant) |>
#   dplyr::pull()
# 
# df_abs_or_pres_pois <- df_abs_or_pres |>
#   dplyr::select(PIORMB, dplyr::matches(abs_or_pres_sig))
# 
# df_abs_or_pres_pois <- fct_regression_poison_uni(df_abs_or_pres_pois)
# 
# 
# ## 3.2 Ulcerados vs não ulcerados ----
# 
# ulc_sig <- dplyr::inner_join(df_chi2_u, df_fisher_u, by = "variant") |>
#   dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
#                   `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
#   dplyr::select(variant) |>
#   dplyr::pull()
# 
# df_ulc_pois <- df_full |>
#   dplyr::select(PIORMB, dplyr::matches(ulc_sig))
# 
# df_ulc_pois <- fct_regression_poison_uni(df_ulc_pois)
# 
# 
# ## 3.3 Severos vs não severos ----
# 
# sev_sig <- dplyr::inner_join(df_chi2_s, df_fisher_s, by = "variant") |>
#   dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
#                   `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
#   dplyr::select(variant) |>
#   dplyr::pull()
# 
# df_sev_pois <- df_sev |>
#   dplyr::select(PIORMB, dplyr::matches(sev_sig))
# 
# df_sev_pois <- fct_regression_poison_uni(df_sev_pois)
# 
# 
# # 4. Regressão Poisson Multivariada ----
# ## 4.1 Presença vs ausência ----
# 
# df_abs_or_pres_pois_uni <- df_abs_or_pres_pois |>
#   dplyr::filter(`p-value(poisson)` <= p_value) # Nenhuma variante significativa
# 
# ## 4.2 Ulcerados vs não ulcerados ----
# 
# ulc_sig_pois_uni <- df_ulc_pois |>
#   dplyr::filter(`p-value(poisson)` <= p_value) |> # 16 variantes significativas
#   dplyr::pull(1)
# 
# df_ulc_sig_pois_uni <- df_full |>
#   dplyr::select(PIORMB, dplyr::matches(ulc_sig_pois_uni))
# 
# summary(ulc_pois_mult <- glm(PIORMB ~ ., data = df_ulc_sig_pois_uni, family = "poisson"))
# step(ulc_pois_mult) # Poderia criar uma função para remover uma variante por vez, mas
#                     # Acho que resultaria nisso de qualquer forma.
# 
# summary(ulc_pois_mult <- glm(formula = PIORMB ~ ABCC2_rs2273697_MD + ABCC6_rs2856585_MA +
#                        GSTA1_rs1051775_MD + GSTP1_rs4891_MA + GSTP1_rs4891_MD +
#                        HSP90AA1_rs4947_MD, family = "poisson", data = df_ulc_sig_pois_uni))
# 
# # Nem todas variantes são significativas para 0.05, mas a remoção da variante não significativa
# # não parece melhorar o modelo. Contudo, como o objetivo é encontrar aquelas que são significativas
# # a 0.05 no modelo multivariado, podemos removê-la ainda.
# 
# hnp::hnp(ulc_pois_mult, resid.type = "deviance") # é desejado que os pontos (as observações)
# # estejam contidos dentro dos limites de confiança (que são gerados por simulações de MC).
# # Isso representa o quão bem ajustados os dados estão ao modelo assumido (nesse caso poisson).
# # É aceitável que até 5% das observações estejam fora dos limites, porém que mantenham-se
# # em torno do centro dos limites de confiança. Como há um comportamento horizontal (ao fim), pode
# # ser evidênciad de que a poisson não tem sentido para esses dados.
# 
# 
# ## 4.3 Severos vs não severos ----
# 
# df_sev_sig_pois_uni <- df_sev_pois |>
#   dplyr::filter(`p-value(poisson)` <= p_value) # Três variantes significativas
# 

# 5. Regressão binomial univariada ----
## 5.1 Presença vs ausência ----
abs_or_pres_sig <- dplyr::inner_join(df_chi2_a_p, df_fisher_a_p, by = "variant") |> 
  dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
                  `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
  dplyr::select(variant) |>
  dplyr::pull()

df_abs_or_pres_binom <- df_abs_or_pres |> 
  dplyr::select(PIORMB, dplyr::matches(abs_or_pres_sig))

df_abs_or_pres_binom <- fct_regression_binom_uni(df_abs_or_pres_binom)
df_abs_or_pres_binom$`p-value(binomial)` <- round(as.numeric(df_abs_or_pres_binom$`p-value(binomial)`), 6)

df_abs_or_pres_binom <- fct_select_variant(df_abs_or_pres_binom)

## 5.2 Não ulcerados vs ulcerados ----
ulc_sig <- dplyr::inner_join(df_chi2_u, df_fisher_u, by = "variant") |> 
  dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
                  `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
  dplyr::select(variant) |>
  dplyr::pull()

df_ulc_binom <- df_ulc |> 
  dplyr::select(PIORMB, dplyr::matches(ulc_sig))

df_ulc_binom <- fct_regression_binom_uni(df_ulc_binom)
df_ulc_binom$`p-value(binomial)` <- round(as.numeric(df_ulc_binom$`p-value(binomial)`), 6)

df_ulc_binom <- fct_select_variant(df_ulc_binom)


## 5.3 Não severo vs severo ----
sev_sig <- dplyr::inner_join(df_chi2_s, df_fisher_s, by = "variant") |> 
  dplyr::filter(`p-value(Chi-2)` < p_value | `p-value(Chi-2)-MC` < p_value |
                  `p-value(fisher)` < p_value | `p-value(fisher)-MC` < p_value) |>
  dplyr::select(variant) |>
  dplyr::pull()

df_sev_binom <- df_sev |> 
  dplyr::select(PIORMB, dplyr::matches(sev_sig))

df_sev_binom <- fct_regression_binom_uni(df_sev_binom)
df_sev_binom$`p-value(binomial)` <- round(as.numeric(df_sev_binom$`p-value(binomial)`), 6)

df_sev_binom <- fct_select_variant(df_sev_binom)

# 6. Regressão binomial multivariada ----
## 6.1 Presença vs ausência ----

abs_or_pres_sig_binom_uni <- df_abs_or_pres_binom |> 
  dplyr::filter(`p-value(binomial)` <= p_value) |> # 18 variantes significativas
  dplyr::pull(1)

df_abs_or_pres_sig_binom_uni <- df_abs_or_pres |> 
  dplyr::select(PIORMB, dplyr::matches(abs_or_pres_sig_binom_uni))

summary(abs_or_pres_binom_mult <- glm(PIORMB ~ ., data = df_abs_or_pres_sig_binom_uni, family = "binomial"))
step(abs_or_pres_binom_mult) # Poderia criar uma função para remover uma variante por vez, mas
# acho que resultaria nisso de qualquer forma.

## Modelo sugerido pelo step
summary(abs_or_pres_binom_mult  <- glm(formula = PIORMB ~ ABCA3_rs1319979593_MU + ABCC2_rs1137968_MU + 
                                         SLC19A1_rs12659_MD + ABCC2_rs3740066_MR + GSTP1_rs4891_MA + 
                                         ABCC4_chr1395164412_MA + ABCC1_rs35587_MA + SLC31A1_chr9113258719_MU,
                                       family = "binomial", data = df_abs_or_pres_sig_binom_uni))

# Nem todas variantes são significativas a 0.05 também, contudo a remoção de mais variantes
# pioraria a qualidade do ajuste. Contudo, como o interesse aqui é inferencial, poderia ser feito
# o procedimento de remoção de variantes, uma a uma, até que todas fossem significativas.


## Continuando remoção até chegarmos no modelo apenas com variantes significativas
summary(abs_or_pres_binom_mult  <- glm(formula = PIORMB ~ ABCA3_rs1319979593_MU + ABCC2_rs1137968_MU + 
                                         SLC19A1_rs12659_MD + GSTP1_rs4891_MA + ABCC1_rs35587_MA,
                                       family = "binomial", data = df_abs_or_pres_sig_binom_uni))

hnp::hnp(abs_or_pres_binom_mult, resid.type = "deviance")

# Agora, como esperado, a distribuição assumida para os dados é muito melhor (e sem dúvidas certa),
# visto que nem mesmo uma observação sai dos limites de confiança.

## 6.2 Não ulcerados vs ulcerados ----

ulc_sig_binom_uni <- df_ulc_binom |> 
  dplyr::filter(`p-value(binomial)` <= p_value) |> # 21 variantes significativas
  dplyr::pull(1)

df_ulc_sig_binom_uni <- df_ulc |> 
  dplyr::select(PIORMB, dplyr::matches(ulc_sig_binom_uni))

summary(ulc_binom_mult <- glm(PIORMB ~ ., data = df_ulc_sig_binom_uni, family = "binomial"))
step(ulc_binom_mult) # Poderia criar uma função para remover uma variante por vez, mas
# acho que resultaria nisso de qualquer forma.

## Modelo sugerido pelo step
summary(ulc_binom_mult  <- glm(formula = PIORMB ~ ABCC2_rs2273697_MA + GSTP1_rs1695_MD + 
                                 ABCC6_rs2856585_MA + ABCA3_rs1319979593_MU + SLCO6A1_rs6884141_MR + 
                                 HSP90AA1_rs4947_MD + GSTA1_rs1051775_MD + ABCC2_rs1137968_MU,
                               family = "binomial", data = df_ulc_sig_binom_uni))

# Nem todas variantes são significativas a 0.05 também, contudo a remoção de mais variantes
# pioraria a qualidade do ajuste. Contudo, como o interesse aqui é inferencial, poderia ser feito
# o procedimento de remoção de variantes, uma a uma, até que todas fossem significativas.


## Continuando remoção até chegarmos no modelo apenas com variantes significativas
summary(ulc_binom_mult  <- glm(formula = PIORMB ~ ABCC2_rs2273697_MA + GSTP1_rs1695_MD + 
                                 ABCC6_rs2856585_MA + HSP90AA1_rs4947_MD + GSTA1_rs1051775_MD + ABCC2_rs1137968_MU,
                               family = "binomial", data = df_ulc_sig_binom_uni))

hnp::hnp(ulc_binom_mult, resid.type = "deviance")

# Agora, como esperado, a distribuição assumida para os dados é muito melhor (e sem dúvidas certa),
# visto que nem mesmo uma observação sai dos limites de confiança.


## 6.3 Não severo vs severo ----

sev_sig_binom_uni <- df_sev_binom |> 
  dplyr::filter(`p-value(binomial)` <= p_value) |> # 8 variantes significativas
  dplyr::pull(1)

df_sev_sig_binom_uni <- df_sev |> 
  dplyr::select(PIORMB, dplyr::matches(sev_sig_binom_uni))

summary(sev_binom_mult <- glm(PIORMB ~ ., data = df_sev_sig_binom_uni, family = "binomial"))
step(sev_binom_mult) # Poderia criar uma função para remover uma variante por vez, mas
# acho que resultaria nisso de qualquer forma.

## Modelo sugerido pelo step
summary(sev_binom_mult  <- glm(formula = PIORMB ~ ABCA3_rs1319979593_MU + ABCC2_rs1137968_MU + 
                                         SLC19A1_rs12659_MD + ABCC2_rs3740066_MR + GSTP1_rs4891_MA + 
                                         ABCC4_chr1395164412_MA + ABCC1_rs35587_MA + SLC31A1_chr9113258719_MU,
                                       family = "binomial", data = df_sev_sig_binom_uni))

# Nem todas variantes são significativas a 0.05 também, contudo a remoção de mais variantes
# pioraria a qualidade do ajuste. Contudo, como o interesse aqui é inferencial, poderia ser feito
# o procedimento de remoção de variantes, uma a uma, até que todas fossem significativas.


## Continuando remoção até chegarmos no modelo apenas com variantes significativas
summary(sev_binom_mult  <- glm(formula = PIORMB ~ ABCA3_rs1319979593_MU + ABCC2_rs1137968_MU + 
                                         SLC19A1_rs12659_MD + GSTP1_rs4891_MA + ABCC1_rs35587_MA,
                                       family = "binomial", data = df_sev_sig_binom_uni))

hnp::hnp(sev_binom_mult, resid.type = "deviance")

# Agora, como esperado, a distribuição assumida para os dados é muito melhor (e sem dúvidas certa),
# visto que nem mesmo uma observação sai dos limites de confiança.
######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_calc_perc.R")
source("R/fct_aux/func_remove_columns.R")
source("R/fct_aux/func_tables.R")
source("R/fct_aux/func_testes_chi2_fisher.R")
source("R/fct_aux/func_writ_organ_xlsx.R")


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

# 3 - Avaliando modelos (dominante, recessivo, aditivo)

df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])
### 3.1 - Tabela para modelos aditivos ----
df_geno_aditivo <- df_modelos_genotipos |> 
  dplyr::select(dplyr::ends_with(c("PIORMB", "MO")))

df_tabela_aditivo <- fct_table_ad(df_geno_aditivo, "pres")

#### Adicionar p-value a tabela
df_chi2_adit <- fct_break_gene_variant_ends(df_chi2, "_MA")

df_chi2_adit$`p-value(Chi-2)` <- round(as.numeric(df_chi2_adit$`p-value(Chi-2)`), 4)
df_tabela_aditivo_p <- dplyr::inner_join(df_tabela_aditivo, df_chi2_adit, by = c("Gene", "variant"))

## Calculando percentual da coluna de frequencia
df_tabela_aditivo_p_perc <- fct_calc_perc(df_tabela_aditivo_p)

#### Criar tabela XLSX do df criado (modelo aditivo)
excel_name <- "tabela_modelo_aditivo"
fct_merge_cels(df_tabela_aditivo_p_perc, excel_name)


### 1.2.2 - Tabela para modelo recessivo ----
df_geno_recessivo <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MU"))

df_tabela_recessivo <- fct_table_rec_dom_un(df_geno_recessivo, "pres")

#### Adicionar p-value a tabela
df_chi2_rec <- fct_break_gene_variant_ends(df_chi2, "_MR")

df_chi2_rec$`p-value(Chi-2)` <- round(as.numeric(df_chi2_rec$`p-value(Chi-2)`), 4)
df_tabela_recessivo_p <- dplyr::inner_join(df_tabela_recessivo, df_chi2_rec, by = c("Gene", "variant"))

## Calculando percentual da coluna de frequencia
df_tabela_recessivo_p_perc <- fct_calc_perc(df_tabela_recessivo_p)

#### Criar tabela XLSX do df criado (modelo recessivo)
excel_name <- "tabela_modelo_recessivo"
fct_merge_cels(df_tabela_recessivo_p_perc, excel_name)


### 1.2.3 - Tabela para modelo dominante ----
df_geno_dominante <- fct_remove_columns(df_modelos_genotipos, c("MA", "MR", "MU"))

df_tabela_dominante <- fct_table_rec_dom_un(df_geno_dominante, "pres")

#### Adicionar p-value a tabela
df_chi2_dom <- fct_break_gene_variant_ends(df_chi2, "_MR")

df_chi2_dom$`p-value(Chi-2)` <- round(as.numeric(df_chi2_dom$`p-value(Chi-2)`), 4)
df_tabela_dominante_p <- dplyr::inner_join(df_tabela_dominante, df_chi2_rec, by = c("Gene", "variant"))

## Calculando percentual da coluna de frequencia
df_tabela_dominante_p_perc <- fct_calc_perc(df_tabela_dominante_p)

#### Criar tabela XLSX do df criado (modelo dominante)

excel_name <- "tabela_modelo_dominante"
fct_merge_cels(df_tabela_dominante_p_perc, excel_name)


#### Avaliando modelo
# summary()

#### Resultado com multicolinearidade, tratar multicolinearidade

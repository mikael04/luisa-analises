######################################################### #
##' Script para avaliar os testes chi-quadrado e de fisher
##' Avaliando com simulação de monte carlo e sem
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_remove_columns.R")
source("R/fct_aux/func_tables.R")
source("R/fct_aux/func_writ_organ_xlsx.R")
source("R/fct_aux/func_calc_perc.R")


# 1 - Lendo Base de dados ----
df_abs_or_pres <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav") |>
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD"))) |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0))

set.seed(42)

# 2 - Rodando testes ----

#### df_chi2
df_chi2 <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_chi2) <- c("variant", "p-value(Chi-2)", "p-value(Chi-2)-MC")

#### df_fisher
df_fisher <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_fisher) <- c("variant", "p-value(fisher)", "p-value(fisher)-MC")

colnames_df <- colnames(df_abs_or_pres)
## Rodando teste para todas as variantes
for(i in (2:ncol(df_abs_or_pres))){
  ## Chi-quadrado
  testes_chi2 <- chisq.test(table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i])), simulate.p.value = F)
  testes_chi2_mc <- chisq.test(table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i])), simulate.p.value = TRUE)
  df_chi2[i-1, ] = c(colnames_df[i], testes_chi2$p.value, testes_chi2_mc$p.value)
  ## Fisher
  tabela_fisher <- table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i]))
  if(ncol(tabela_fisher) > 1){
    testes_fisher <- fisher.test(tabela_fisher, simulate.p.value = F, workspace = 2e8)
    testes_fisher_mc <- fisher.test(tabela_fisher, simulate.p.value = TRUE)
    df_fisher[i-1, ] = c(colnames_df[i], testes_fisher$p.value, testes_fisher_mc$p.value)
  }else{
    df_fisher[i-1, ] = c(colnames_df[i], -1, -1)
  }
}

## Ordenando por p-value
df_chi2 <- df_chi2 |> 
  dplyr::arrange(variant)
# 2 - Investigando testes chi-quadrado e 


df_both_tests <- data.frame(df_chi2$variant, df_chi2$`p-value(Chi-2)`, df_chi2$`p-value(Chi-2)-MC`, df_fisher$`p-value(fisher)`, df_fisher$`p-value(fisher)-MC`)
colnames(df_both_tests) <- c("variant", "p-value(Chi-2)", "p-value(Chi-2)-MC", "p-value(fisher)", "p-value(fisher)-MC")
write.csv(df_both_tests, "data-raw/df_testes_chi2_fisher_abs_or_pres.csv")

## Avaliando variantes genéticas significativas
# max_var <-  as.data.frame(abs(as.numeric(df_both_tests$`p-value(Chi-2)`) - as.numeric(df_both_tests$`p-value(Chi-2)-MC`)))
# max(max_var)
# 
# count_sig_chi2 <- df_both_tests |> 
#   dplyr::filter(`p-value(Chi-2)` < 0.05)
#   # dplyr::filter(`p-value(Chi-2)` < 0.05) |> 
#   # dplyr::count() |> 
#   # dplyr::pull()
# 
# count_sig_chi2_mc <- df_both_tests |> 
#   dplyr::filter(`p-value(Chi-2)-MC` < 0.05)
# 
# count_sig_fisher <- df_both_tests |> 
#   dplyr::filter(`p-value(fisher)` < 0.05)
# 
# count_sig_fisher_mc <- df_both_tests |> 
#   dplyr::filter(`p-value(fisher)-MC` < 0.05)

## Investigando alguns testes e seus resultados individualmente
df_ABCB1_rs9282564 <- df_abs_or_pres |> 
  dplyr::select(PIORMB, ABCB1_rs9282564_MA)

table(df_ABCB1_rs9282564)
testes_chi2 <- chisq.test(table(df_ABCB1_rs9282564), simulate.p.value = T)

df_ABCC1_rs8187858 <- df_abs_or_pres |> 
  dplyr::select(PIORMB, ABCC1_rs8187858_MR)

table(df_ABCC1_rs8187858)
testes_chi2 <- chisq.test(table(df_ABCC1_rs8187858), simulate.p.value = T)

df_SLC19A1_rs12659 <- df_abs_or_pres |> 
  dplyr::select(PIORMB, SLC19A1_rs12659_MR)

set.seed(44)
table(df_SLC19A1_rs12659)
chisq.test(table(df_SLC19A1_rs12659), simulate.p.value = T)

df_ABCC4_chr1395164412 <- df_abs_or_pres |> 
  dplyr::select(PIORMB, ABCC4_chr1395164412_MR)

table(df_ABCC4_chr1395164412)
chisq.test(table(df_ABCC4_chr1395164412), simulate.p.value = T)

#### df tabela de contingencia (checagem de valores menor que 5)
df_checks <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_checks) <- c("variant", "Pass (<5 assump)", "ncol > 1")
colnames_df <- colnames(df_abs_or_pres)

for(i in (2:ncol(df_abs_or_pres))){
  tabela_conting <- table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i]))
  check_ncol <- F
  if(ncol(tabela_conting) > 1){
    check_ncol <- T
  }
  check_assump <- T
  for(j in 1:length(tabela_conting)){
    if(tabela_conting[j] < 5){
      check_assump <- F
    }
  }
  df_checks[i-1, ] = c(colnames_df[i], check_assump, check_ncol)
}

df_MTHFR_rs4846051 <- df_abs_or_pres |> 
  dplyr::select(PIORMB, MTHFR_rs4846051_MA)

table(df_MTHFR_rs4846051)

# tabela_fisher <- table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i]))
# testes_fisher <- fisher.test(tabela_fisher, simulate.p.value = F, workspace = 2e8)

# i <- 64
# testes_chi2 <- table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i]))

## Investigando diferença entre Presença2 e PIORMB

df_presenca_piormb <- df_full |>
  dplyr::select(ID, Presença2, PIORMB) |>
  dplyr::mutate(PIORMB_cat = ifelse(PIORMB > 0, 1, 0)) |>
  dplyr::select(ID, Presença2, PIORMB_cat, PIORMB)

# 
# count(unique(df_full$Paciente))
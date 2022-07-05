######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)
source("R/fct_aux/func_remove_columns.R")
source("R/fct_aux/func_tables.R")
source("R/fct_aux/func_writ_organ_xlsx.R")


# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav") |> 
  dplyr::filter(!is.na(PIORMB))

df <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD", "MU")))

### 1.0.1 - Rodando um qui-quadrado para testes ----
table_ABCC2_RS2273697 <- table(df$PIORMB, df$ABCC2_rs2273697_MA)
chisq.test(table_ABCC2_RS2273697)

### 1.0.2 - Rodando um Teste de Fisher para testes ----
fisher.test(table_ABCC2_RS2273697)
# fisher.test(table_ABCC2_RS2273697, simulate.p.value = T)

## 1.1 - Rodando qui-quadrado para todas as variantes ----

#### df_chi2
df_chi2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_chi2) <- c("variant", "p-value(Chi-2)")

#### df_fisher
df_fisher <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_fisher) <- c("variant", "p-value(fisher)")

colnames_df <- colnames(df)
#### Rodando teste para todas as variantes
for(i in (2:ncol(df))){
  ## Chi-quadrado
  testes_chi2 <- chisq.test(table(unlist(df[,1]), unlist(df[,i])), simulate.p.value = TRUE)
  df_chi2[i-1, ] = c(colnames_df[i], testes_chi2$p.value)
  ## Fisher
  tabela_fisher <- table(unlist(df[,1]), unlist(df[,i]))
  if(ncol(tabela_fisher) > 1){
    testes_fisher <- fisher.test(tabela_fisher, simulate.p.value = TRUE)
    df_fisher[i-1, ] = c(colnames_df[i], testes_fisher$p.value)
  }else{
    df_fisher[i-1, ] = c(colnames_df[i], -1)
  }
}

#### Ordenando por p-value
df_chi2 <- df_chi2 |> 
  dplyr::arrange(variant)

## 1.2 - Criando tabela Chi-2 ----

df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])

### 1.2.1 - Tabela para modelos aditivos ----
df_geno_aditivo <- df_modelos_genotipos |> 
  dplyr::select(dplyr::ends_with(c("PIORMB", "MO")))

df_tabela_aditivo <- fct_table_ad(df_geno_aditivo)

#### Adicionar p-value a tabela
df_chi2_adit <- fct_break_gene_variant_ends(df_chi2, "_MA")

df_chi2_adit$`p-value(Chi-2)` <- round(as.numeric(df_chi2_adit$`p-value(Chi-2)`), 4)
df_tabela_aditivo_p <- dplyr::inner_join(df_tabela_aditivo, df_chi2_adit, by = c("Gene", "variant"))

#### Criar tabela XLSX do df criado (modelo aditivo)
excel_name <- "tabela_modelo_aditivo"
fct_merge_cels(df_tabela_aditivo_p, excel_name)


### 1.2.2 - Tabela para modelo recessivo ----
df_geno_recessivo <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MU"))

df_tabela_recessivo <- fct_table_rec_dom_un(df_geno_recessivo)

#### Adicionar p-value a tabela
df_chi2_rec <- fct_break_gene_variant_ends(df_chi2, "_MR")

df_chi2_rec$`p-value(Chi-2)` <- round(as.numeric(df_chi2_rec$`p-value(Chi-2)`), 4)
df_tabela_recessivo_p <- dplyr::inner_join(df_tabela_recessivo, df_chi2_rec, by = c("Gene", "variant"))
#### Criar tabela XLSX do df criado (modelo recessivo)

excel_name <- "tabela_modelo_recessivo"
fct_merge_cels(df_tabela_recessivo_p, excel_name)


### 1.2.3 - Tabela para modelo dominante ----
df_geno_dominante <- fct_remove_columns(df_modelos_genotipos, c("MA", "MR", "MU"))

df_tabela_dominante <- fct_table_rec_dom_un(df_geno_dominante)

#### Adicionar p-value a tabela
df_chi2_dom <- fct_break_gene_variant_ends(df_chi2, "_MR")

df_chi2_dom$`p-value(Chi-2)` <- round(as.numeric(df_chi2_dom$`p-value(Chi-2)`), 4)
df_tabela_dominante_p <- dplyr::inner_join(df_tabela_dominante, df_chi2_rec, by = c("Gene", "variant"))

#### Criar tabela XLSX do df criado (modelo dominante)

excel_name <- "tabela_modelo_dominante"
fct_merge_cels(df_tabela_dominante_p, excel_name)


### 1.2.4 - Tabela para modelo unico ----
df_geno_unico <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MR"))

df_tabela_unico <- fct_table_rec_dom_un(df_geno_unico)

#### Adicionar p-value a tabela
df_chi2_un <- fct_break_gene_variant_ends(df_chi2, "_MU")

df_chi2_un$`p-value(Chi-2)` <- round(as.numeric(df_chi2_un$`p-value(Chi-2)`), 4)
df_tabela_unico_p <- dplyr::inner_join(df_tabela_unico, df_chi2_un, by = c("Gene", "variant"))

#### Criar tabela XLSX do df criado (modelo unico)

excel_name <- "tabela_modelo_unico"
fct_merge_cels(df_tabela_unico_p, excel_name)


# 2 - Criando modelo 


#### Avaliando modelo
# summary()

#### Resultado com multicolinearidade, tratar multicolinearidade

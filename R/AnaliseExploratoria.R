######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)


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
colnames(df_chi2) <- c("Variante", "p-value(Chi-2)")

#### df_fisher
df_fisher <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_fisher) <- c("Variante", "p-value(fisher)")

colnames_df <- colnames(df)
#### Rodando teste para todas as variantes
for(i in (2:160)){
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
  dplyr::arrange(Variante)

## 1.2 Criando tabela Chi-2

df_modelos_genotipos <- data.frame(df_full[, 32], df_full[, 57:324])

# df_rs17222723 <- data.frame(df_modelos_genotipos[, 1], df_modelos_genotipos[, 8:9])

df_aux_geno_aditivo <- df_modelos_genotipos |> 
  dplyr::select(!ends_with(c("MR", "MD", "MU")))

colnames_df_ad <- colnames(df_aux_geno_aditivo)
colnames_just_ad <- c("PIORMB")
for(i in 2:(length(colnames_df_ad)-1)){
  this_str <- substr(colnames_df_ad[i], 1, nchar(colnames_df_ad[i])-3)
  next_str <- substr(colnames_df_ad[i+1], 1, nchar(colnames_df_ad[i+1])-3)
  previous_str <- substr(colnames_df_ad[i-1], 1, nchar(colnames_df_ad[i-1])-3)
  if(this_str == previous_str){
    ## Mantém a variante
    colnames_just_ad <- c(colnames_just_ad, colnames_df_ad[i], colnames_df_ad[i-1])
  }else{
    ## Não adiciona a variante, pois não tem modelo aditivo
  }
  if(this_str == next_str){
    ## Mantém a variante
    colnames_just_ad <- c(colnames_just_ad, colnames_df_ad[i], colnames_df_ad[i+1])
  }else{
    ## Não adiciona a variante, pois não tem modelo aditivo
  }
}

colnames_just_ad <- unique(unlist(strsplit(colnames_just_ad, " ")))

df_geno_aditivo <- df_aux_geno_aditivo |> 
  dplyr::select(any_of(colnames_just_ad))

df_geno_aditivo_test <- df_geno_aditivo[,1:3] |> 
  dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0)) |> 
  dplyr::group_by(ABCC2_rs2273697_MO, PIORMB) |> 
  dplyr::mutate(count = n()) |> 
  dplyr::distinct(ABCC2_rs2273697_MO, PIORMB, .keep_all = T) |>
  dplyr::select(-ABCC2_rs2273697_MA) |> 
  dplyr::ungroup()

df_geno_aditivo_test <- df_geno_aditivo_test |> 
  dplyr::arrange(ABCC2_rs2273697_MO, PIORMB)


df_aux <- df_geno_aditivo_test |> 
  tidyr::pivot_wider(names_from = PIORMB, values_from = count)

# df_aux <- df_full |> 
#   dplyr::select(PIORMB, ABCC2_rs2273697_MO, ABCC2_rs2273697_MA, ABCC2_rs2273697_MR, ABCC2_rs2273697_MD)

## Lendo outra  tabela auxiliar

# df_modelos_genotipos <- readxl::read_xlsx("data-raw/modelos_genotipos.xlsx")
# 
# df_modelos_genotipos_transposed <- as.data.frame(t(df_modelos_genotipos)) %>% 
#   tibble::rownames_to_column()
# 
# df_modelos_genotipos_transposed[1,6] <- "ALELO REF. 2"
# colnames(df_modelos_genotipos_transposed) <- df_modelos_genotipos_transposed[1,]
# df_m_g_t <- df_modelos_genotipos_transposed[-1, ]
# 
# df_m_g_t <-df_m_g_t |> 
#   dplyr::filter(!is.na(MODELOS))

# 2 - Criando modelo 


#### Avaliando modelo
# summary()

#### Resultado com multicolinearidade, tratar multicolinearidade

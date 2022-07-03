######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)


# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav")

df <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD")))

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
table_ABCC2_RS2273697
colnames_df[2]

df_aux <- df_full |> 
  dplyr::select(PIORMB, ABCC2_rs2273697_MO, ABCC2_rs2273697_MA, ABCC2_rs2273697_MR, ABCC2_rs2273697_MD)

## Lendo outra  tabela auxiliar

df_modelos_genotipos <- readxl::read_xlsx("data-raw/modelos_genotipos.xlsx")

df_modelos_genotipos_transposed <- as.data.frame(t(df_modelos_genotipos)) %>% 
  tibble::rownames_to_column()

df_modelos_genotipos_transposed[1,6] <- "ALELO REF. 2"
colnames(df_modelos_genotipos_transposed) <- df_modelos_genotipos_transposed[1,]
df_m_g_t <- df_modelos_genotipos_transposed[-1, ]

df_m_g_t <-df_m_g_t |> 
  dplyr::filter(!is.na(MODELOS))

# 2 - Criando modelo 


#### Avaliando modelo
summary()

#### Resultado com multicolinearidade, tratar multicolinearidade

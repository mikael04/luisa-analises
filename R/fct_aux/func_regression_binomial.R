######################################################### #
##' Script para avaliar o retorno das funções usadas
##' 
##' Autores: Mikael e Matheus
##' Data: 09/07/2022
######################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_regression_binom_uni <- function(df){
  # df <- df_abs_or_pres
  ## Criando dfs base regressão binomial
  df_binom <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_binom) <- c("variant", "p-value(binomial)")
  colnames_df <- colnames(df)
  ## Rodando teste para todas as variantes
  for(i in (2:ncol(df))){
    ## Binomial
    df_binomial <- data.frame(df[,1], df[,i])
    binom.model <- summary(glm(PIORMB ~ ., df_binomial, family = "binomial"))
    
    if(nrow(binom.model$coefficients) > 1){
      ## Para um modelo com apenas uma variavel independente, 
      ## o coeficiente com p-valor dela será nesta posição
      df_binom[i-1, ] = c(colnames_df[i], binom.model$coefficients[2,4])
    }else{
      ## Caso seja NA o p-value, terá apenas uma linha e adicionaremos por aqui
      df_binom[i-1, ] = c(colnames_df[i], NA)
    }
  }
  df_binom
}

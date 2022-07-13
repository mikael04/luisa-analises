######################################################### #
##' Script para avaliar o retorno das funções usadas
##' 
##' Autores: Mikael e Matheus
##' Data: 08/07/2022
######################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_regression_poison_uni <- function(df){
  # df <- df_abs_or_pres
  ## Criando dfs base regressão de poisson
  ### df_poi
  df_poi <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df_poi) <- c("variant", "p-value(poisson)")
  colnames_df <- colnames(df)
  ## Rodando teste para todas as variantes
  for(i in (2:ncol(df))){
    ## Poisson
    df_poisson <- data.frame(df[,1], df[,i])
    poisson.model <- summary(glm(PIORMB ~ ., df_poisson, family = "poisson"))
    
    if(nrow(poisson.model$coefficients) > 1){
      ## Para um modelo com apenas uma variavel independente, o coeficiente com p-valor dela será nesta posição
      df_poi[i-1, ] = c(colnames_df[i], poisson.model$coefficients[2,4])
    }else{
      ## Caso seja NA o p-value, terá apenas uma linha e adicionaremos por aqui
      df_poi[i-1, ] = c(colnames_df[i], NA)
    }
  }
  df_poi
}

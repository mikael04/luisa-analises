##################################################################### #
##' Função para alterar valor de testes (chi2 e fisher) para modelo
##' usando simulação de Monte Carlo
##' 
##' Autores: Mikael e Matheus
##' Data: 18/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----

## Função que altera valor base de teste para valor de teste com simulação
## de monte carlo para o teste qui-quadrado

fct_use_MC_chi2 <- function(df_chi2, use_MC){
  if(use_MC){
    df_chi2$`p-value(Chi-2)` <- round(as.numeric(df_chi2$`p-value(Chi-2)-MC`), 6)
  }else{
    df_chi2$`p-value(Chi-2)` <- round(as.numeric(df_chi2$`p-value(Chi-2)`), 6)
  }
  df_chi2
}


## Função que altera valor base de teste para valor de teste com simulação
## de monte carlo para o teste de fisher

fct_use_MC_fisher <- function(df_chi2, use_MC){
  if(use_MC){
    df_chi2$`p-value(fisher)` <- round(as.numeric(df_chi2$`p-value(fisher)-MC`), 6)
  }else{
    df_chi2$`p-value(fisher)` <- round(as.numeric(df_chi2$`p-value(fisher)`), 6)
  }
  df_chi2
}
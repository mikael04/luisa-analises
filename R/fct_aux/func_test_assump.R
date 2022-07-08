################################################################# #
##' Script para avaliar pressupostos do chi-quadrado e de fisher
##' Testando n > 5 em todas as linhas e colunas
##' Testando se possui mais de uma coluna
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
################################################################# #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_test_assump <- function(df_abs_or_pres){
  #### df tabela de contingencia (checagem de valores menor que 5)
  df_checks <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_checks) <- c("variant", "Pass (<5 assump)", "ncol > 1")
  colnames_df <- colnames(df_abs_or_pres)
  
  for(i in (2:ncol(df_abs_or_pres))){
    tabela_conting <- table(unlist(df_abs_or_pres[,1]), unlist(df_abs_or_pres[,i]))
    ## Rodando teste para ter tabela de valores esperados
    teste_chi2 <- chisq.test(tabela_conting, simulate.p.value = TRUE, B = 10000)
    expected <- teste_chi2$expected
    
    check_assump <- F
    ## Testando para todos os valores da tabela de contingencia
    for(j in 1:length(expected)){
      if(expected[j] < 5){
        check_assump <- T
      }
    }
    ## Checando se possui mais de uma coluna
    check_ncol <- F
    if(ncol(tabela_conting) > 1){
      check_ncol <- T
    }
    ## adicionando ao df
    df_checks[i-1, ] = c(colnames_df[i], check_assump, check_ncol)
  }
  df_checks
}

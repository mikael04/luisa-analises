######################################################### #
##' Script para avaliar os testes chi-quadrado e de fisher
##' Avaliando com simulação de monte carlo e sem
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_testes_chi2_fisher <- function(df){
  set.seed(42)
    # df <- df_abs_or_pres
  ## Criando dfs base chi2 e de fisher
  ### df_chi2
  df_chi2 <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_chi2) <- c("variant", "p-value(Chi-2)", "p-value(Chi-2)-MC")
  
  ### df_fisher
  df_fisher <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_fisher) <- c("variant", "p-value(fisher)", "p-value(fisher)-MC")
  colnames_df <- colnames(df)
  
  ## Rodando teste para todas as variantes
  for(i in (2:ncol(df))){
    ## Chi-quadrado
    testes_chi2 <- chisq.test(table(unlist(df[,1]), unlist(df[,i])), simulate.p.value = F)
    testes_chi2_mc <- chisq.test(table(unlist(df[,1]), unlist(df[,i])), simulate.p.value = TRUE, B = 10000)
    df_chi2[i-1, ] = c(colnames_df[i], testes_chi2$p.value, testes_chi2_mc$p.value)
    ## Fisher
    tabela_fisher <- table(unlist(df[,1]), unlist(df[,i]))
    if(ncol(tabela_fisher) > 1){
      testes_fisher <- fisher.test(tabela_fisher, simulate.p.value = F, workspace = 2e8)
      testes_fisher_mc <- fisher.test(tabela_fisher, simulate.p.value = TRUE, B = 10000)
      df_fisher[i-1, ] = c(colnames_df[i], testes_fisher$p.value, testes_fisher_mc$p.value)
    }else{
      df_fisher[i-1, ] = c(colnames_df[i], -1, -1)
    }
  }
  
  ## Ordenando por p-value
  df_chi2_num <- df_chi2 |> 
    dplyr::arrange(variant) |> 
    dplyr::mutate(`p-value(Chi-2)` = round(as.numeric(`p-value(Chi-2)`), 6),
                  `p-value(Chi-2)-MC` = round(as.numeric(`p-value(Chi-2)-MC`), 6))
  
  df_fisher_num <- df_fisher |> 
    dplyr::arrange(variant) |> 
    dplyr::mutate(`p-value(fisher)` = round(as.numeric(`p-value(fisher)`), 6),
                  `p-value(fisher)-MC` = round(as.numeric(`p-value(fisher)-MC`), 6))
  
  return(list(1, df_chi2_num, df_fisher_num))
}

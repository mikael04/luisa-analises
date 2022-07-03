##################################################################### #
##' Função para criar tabela de modelos aditivos
##' 
##' Autores: Mikael
##' Data: 03/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----

fct_table_ad <- function(df_geno_aditivo){
  df_aux <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(df_aux) <- c("variant", "genotype", "0", "1")
  for(i in 2:ncol(df_geno_aditivo)){
    df_this_col <- data.frame(df_geno_aditivo[,1], df_geno_aditivo[,i])
    cols_df <- c("PIORMB", colnames(df_geno_aditivo[i]))
    colnames(df_this_col) <- cols_df
    
    df_geno_aditivo_test <- df_this_col |> 
      dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0)) |> 
      dplyr::group_by_at(cols_df) |> 
      dplyr::summarise(count = n()) |> 
      dplyr::ungroup()
    
    cols_arrange <- c(colnames(df_geno_aditivo[i]), "PIORMB")
    df_geno_aditivo_test <- df_geno_aditivo_test |> 
      dplyr::arrange_at(cols_arrange)
    
    df_ <- df_geno_aditivo_test |> 
      tidyr::pivot_wider(names_from = PIORMB, values_from = count)
    
    df_$variant <- cols_df[2]
    
    colnames(df_) <- c("genotype", "0", "1", "variant")
    df_ <- df_ |> 
      dplyr::select("variant", "genotype", "0", "1")
    
    df_aux <- rbind(df_aux, df_)
  }
  df_aux
}

##################################################################### #
##' Função para remover colunas de modelo que não deveriam ter ficado
##' 
##' Autores: Mikael
##' Data: 03/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----

fct_remove_columns <- function(df, cols_to_remove){
  df_aux_geno_aditivo <- df |> 
    dplyr::select(!ends_with(cols_to_remove))
  
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
  
  df_aux_geno_aditivo |> 
    dplyr::select(any_of(colnames_just_ad))
}

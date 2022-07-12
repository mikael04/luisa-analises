###################################################################### #
##' Função para selecionar as variantes dentre modelos, de acordo com
##' menor p-value
##' 
##' Autores: Mikael e Matheus
##' Data: 12/07/2022
###################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_select_variant <- function(df){
  teste <- F
  if(teste){
    df <- df_abs_or_pres_binom
  }
  df_filt <- fct_break_gene_variant_model(df)
  df_filt <- df_filt |> 
    dplyr::arrange(`p-value(binomial)`) |> 
    dplyr::distinct(Gene, Variant, .keep_all = T)
  
  df <- fct_unite_gene_variant(df_filt) |> 
    dplyr::select(variant, `p-value(binomial)`)
  df
}

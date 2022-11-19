################################################################# #
##' Script para adicionar intervalo de confiança a tabela de resultados
##' 
##' Autores: Mikael e Matheus
##' Data: 03/11/2022
################################################################# #

# 0 - Scripts e bibliotecas ----

# 1 - Funções ----
fct_confint <- function(df_conf, debug){
  df_conf_aux <- df_conf |> 
    tibble::rownames_to_column(var = "gene_variant_model") |> 
    dplyr::filter(gene_variant_model != "(Intercept)") |> 
    dplyr::mutate(gene_variant_model = gsub(".{3}$", "", gene_variant_model)) |> 
    dplyr::mutate(`2.5 %` = round(`2.5 %`, 4),
                  `97.5 %` = round(`97.5 %`, 4))
}

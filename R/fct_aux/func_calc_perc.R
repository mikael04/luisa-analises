##################################################################### #
##' Função para calcular porcentagem do genótipo
##' 
##' Autores: Mikael e Matheus
##' Data: 04/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_calc_perc <- function(df_tabela_aditivo_p){
  df_tabela_aditivo_p$n_tot <- df_tabela_aditivo_p$`0` + df_tabela_aditivo_p$`1`
  df_tabela_aditivo_p_perc <- df_tabela_aditivo_p |> 
    dplyr::mutate(`abs_mb(%)` = paste0(`0`, " ", "(", round(`0`/n_tot,3)*100 , ")")) |> 
    dplyr::mutate(`pres_mb(%)` = paste0(`1`, " ", "(", round(`1`/n_tot,3)*100 , ")")) |> 
    dplyr::select(Gene, variant, genotype, `abs_mb(%)`, `pres_mb(%)`, `p-value(Chi-2)`)
}

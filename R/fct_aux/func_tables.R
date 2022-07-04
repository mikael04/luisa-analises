##################################################################### #
##' Função para criar tabela de modelos aditivos
##' 
##' Autores: Mikael e Matheus
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
  fct_break_gene_variant(df_aux)
}


fct_break_gene_variant <- function(df_break){
  # df_break <- df_tabela_aditivo
  df_break <- df_break |> 
    dplyr::mutate(variant = substr(variant,1,nchar(variant)-3))
  df_break |> 
    tidyr::separate(variant, c("Gene", "Variant"), sep = "_")
  
}

## Função que retornará uma tabela com os nomes dos genótipos agrupados pelo valor do modelo recessivo
fct_table_rec_get_names <- function(df_geno_others_sel, cols_names){
  # cols_names <- cols_sel
  # df_geno_others_sel <- df_recessivo_names
  df_recessivo_test <- df_geno_others_sel |>
    dplyr::select(all_of(cols_names)) |> ## Cria o df de teste para um modelo
    dplyr::group_by(across(cols_names)) |> ## Determina o agrupamento
    dplyr::summarise() |> 
    dplyr::ungroup() |> 
    dplyr::group_by(across(cols_names[2])) |> ## Define outro agrupamento
    dplyr::mutate(new_geno = paste((!!as.symbol(cols_names[1])), collapse = " or ")) |>
    dplyr::select(all_of(c(cols_names[2], "new_geno"))) |> ## remove a primeira coluna, mantendo apenas a coluna com valores agrupados
    dplyr::distinct(new_geno) |> 
    dplyr::ungroup()
  
  colnames(df_recessivo_test) <- c(cols_names[2], cols_names[1])
  
  df_recessivo_test
}

## Função que recebe tabela com colunas genótipo (MO) e modelo recessivo (MR)
fct_table_rec_dom_un <- function(df_geno_others){
  df_aux <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(df_aux) <- c("variant", "genotype", "0", "1")
  for(i in seq(3, ncol(df_geno_others), 2)){
    df_this_col <- data.frame(df_geno_others[,1], df_geno_others[,i])
    cols_df <- c("PIORMB", colnames(df_geno_others[i]))
    colnames(df_this_col) <- cols_df
    
    df_recessivo_names <- data.frame(df_geno_others[,i-1], df_geno_others[,i])
    cols_sel <- c(colnames(df_geno_others[i-1]), colnames(df_geno_others[i]))
    colnames(df_recessivo_names) <- cols_sel
    
    df_recessivo_names_aux <- fct_table_rec_get_names(df_recessivo_names, cols_sel)
    
    df_geno_others_test <- df_this_col |> 
      dplyr::mutate(PIORMB = ifelse(PIORMB > 0, 1, 0)) |> 
      dplyr::group_by_at(cols_df) |> 
      dplyr::summarise(count = n()) |> 
      dplyr::ungroup()
    
    cols_arrange <- c(colnames(df_geno_others[i]), "PIORMB")
    df_geno_others_test <- df_geno_others_test |> 
      dplyr::arrange(across(all_of(cols_arrange)))
    
    df_ <- df_geno_others_test |> 
      tidyr::pivot_wider(names_from = PIORMB, values_from = count)
    
    df_$variant <- cols_df[2]
    df_ <- dplyr::inner_join(df_, df_recessivo_names_aux, by = cols_df[2])
    
    colnames(df_) <- c("old", "0", "1", "variant", "genotype")
    df_ <- df_ |> 
      dplyr::select("variant", "genotype", "0", "1")
    
    df_aux <- rbind(df_aux, df_)
  }
  fct_break_gene_variant(df_aux)
}

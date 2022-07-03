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

fct_break_gene_variant <- function(df){
  # df <- df_tabela_aditivo
  df <- df |> 
    dplyr::mutate(variant = substr(variant,1,nchar(variant)-3))
  df <- df |> 
    tidyr::separate(variant, c("Gene", "Variant"), sep = "_")
  
}

fct_merge_cels <- function(df){
  ## Criar planilha
  library(openxlsx)
  wb <- createWorkbook()
  ## Add a worksheet
  addWorksheet(wb, "Sheet 1")
  openxlsx::writeData(wb, "Sheet 1", df, startCol = 1, startRow = 1)
  df_aux <- df
  ## Primeiro ordenar df
  df_aux <- fct_merge_cels_pre_method(df_aux)
  ## depois fazer o merge de genes
  fct_merge_cels_gene(df, wb)
  ## e por último fazer o merge de variantes
  fct_merge_cels_variant(df, wb)
  
  ## Save workbook
  saveWorkbook(wb, "data-raw/mergeCellsExample.xlsx", overwrite = TRUE)
}

## Ordenando df para próximos métodos
fct_merge_cels_pre_method <- function(df){
  df |> 
    dplyr::arrange(Gene, Variant, genotype)
}

fct_merge_cels_gene <- function(df, wb){
  df_count_gene <- df |> 
    dplyr::select(Gene) |> 
    dplyr::group_by(Gene) |> 
    dplyr::summarise(n = n())
    dplyr::ungroup()
  
  ## Merge cells: gene
  rows_position_first <- 2
  for (i in 1:nrow(df_count_gene)){
    rows_count <- as.numeric(df_count_gene[i, 2])
    rows_position_last <- rows_position_first - 1 + rows_count
    mergeCells(wb, "Sheet 1", cols = 1, rows = (rows_position_first):rows_position_last)
    rows_position_first <- rows_position_last + 1
  }
}
fct_merge_cels_variant <- function(df, wb){
  df_count_variant <- df |> 
    dplyr::select(Variant) |> 
    dplyr::group_by(Variant) |> 
    dplyr::mutate(n = dplyr::n()) |> 
    dplyr::distinct(Variant, .keep_all = T) |> 
    dplyr::ungroup()
  
  ## Merge cells: gene
  rows_position_first <- 2
  for (i in 1:nrow(df_count_variant)){
    rows_count <- as.numeric(df_count_variant[i, 2])
    rows_position_last <- rows_position_first - 1 + rows_count
    mergeCells(wb, "Sheet 1", cols = 2, rows = (rows_position_first):rows_position_last)
    rows_position_first <- rows_position_last + 1
  }
}

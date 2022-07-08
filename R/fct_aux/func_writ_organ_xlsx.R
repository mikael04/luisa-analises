##################################################################### #
##' Função para organizar e escrever tabelas no formato excel
##' inicialmente apenas contagens, mas é provável que seja usado
##' para agrupar modelos e para adicionar p-value a tabela final
##' 
##' Autores: Mikael e Matheus
##' Data: 04/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_merge_cels <- function(df_merge, excel_name){
  # df_merge <- df_tabela_aditivo_p
  ## Criar planilha
  library(openxlsx)
  wb <- createWorkbook()
  ## Add a worksheet
  addWorksheet(wb, "Sheet 1")
  ## Primeiro ordenar df
  df_aux <- fct_merge_cels_pre_method(df_merge)
  ## Criando workbook
  openxlsx::writeData(wb, "Sheet 1", df_aux, startCol = 1, startRow = 1)
  ## depois fazer o merge de genes
  fct_merge_cels_gene(df_aux, wb)
  ## e por último fazer o merge de variantes
  fct_merge_cels_variant(df_aux, wb)
  
  ## Save workbook
  saveWorkbook(wb, paste0("data-raw/tabelas_modelos/", excel_name ,".xlsx"), overwrite = TRUE)
}

## Ordenando df para próximos métodos
fct_merge_cels_pre_method <- function(df_pre){
  df_pre |> 
    dplyr::arrange(Gene, variant, genotype)
}

## Agrupando linhas de genótipos iguais
fct_merge_cels_gene <- function(df_merge_gene, wb){
  # df_merge_gene <- df_aux
  df_count_gene <- df_merge_gene |> 
    dplyr::select(Gene) |> 
    dplyr::group_by(Gene) |> 
    dplyr::summarise(n = n()) |> 
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

## Agrupando linhas de variantes iguais
fct_merge_cels_variant <- function(df_merge_variant, wb){
  # df_merge_variant <- df_aux
  df_count_variant <- df_merge_variant |> 
    dplyr::select(variant) |> 
    dplyr::group_by(variant) |> 
    dplyr::mutate(n = dplyr::n()) |> 
    dplyr::distinct(variant, .keep_all = T) |> 
    dplyr::ungroup()
  
  ## Merge cells: variant
  rows_position_first <- 2
  for (i in 1:nrow(df_count_variant)){
    rows_count <- as.numeric(df_count_variant[i, 2])
    rows_position_last <- rows_position_first - 1 + rows_count
    mergeCells(wb, "Sheet 1", cols = 2, rows = (rows_position_first):rows_position_last)
    mergeCells(wb, "Sheet 1", cols = 6, rows = (rows_position_first):rows_position_last)
    rows_position_first <- rows_position_last + 1
  }
}

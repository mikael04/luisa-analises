##################################################################### #
##' Função para criar tabela de modelos generalizada
##' Vai ser usada para os três agrupamentos (pres, ulc e sev)
##' 
##' Autores: Mikael e Matheus
##' Data: 06/07/2022
##################################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----

## Função que recebe tabela com colunas genótipo (MO) e modelo aditivo
fct_create_table <- function(df_modelos_genotipos, df_chi2, table, type_group, write_table, teste){
  if(is.null(table)){
    ####################################################################### #
    ## Gerando para o aditivo
    if(teste){
      print("Gerando todas as tabelas")
      print("Iniciando geração de tabela dos modelos aditivos")
    }
    df_geno_aditivo <- df_modelos_genotipos |> 
      dplyr::select(dplyr::ends_with(c("PIORMB", "MO")))
    
    df_tabela_aditivo <- fct_table_ad(df_geno_aditivo, type_group)
    
    #### Adicionar p-value a tabela
    df_chi2_adit <- fct_break_gene_variant_ends(df_chi2, "_MA")
    
    df_chi2_adit$`p-value(Chi-2)` <- round(as.numeric(df_chi2_adit$`p-value(Chi-2)`), 4)
    df_tabela_aditivo_p <- dplyr::inner_join(df_tabela_aditivo, df_chi2_adit, by = c("Gene", "variant"))
    
    ## Calculando percentual da coluna de frequencia
    df_tabela_aditivo_p_perc <- fct_calc_perc(df_tabela_aditivo_p)
    
    #### Criar tabela XLSX do df criado (modelo aditivo)
    excel_name <- paste0("tabela_modelo_aditivo_", type_group)
    fct_merge_cels(df_tabela_aditivo_p_perc, excel_name)
    
    if(teste){
      print("Finalizando geração de tabela dos modelos aditivos")
    }
    ####################################################################### #
    ## Gerando para o recessivo
    if(teste){
      print("Iniciando geração de tabela dos modelos recessivos")
    }
    
    df_geno_recessivo <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MU"))
    
    df_tabela_recessivo <- fct_table_rec_dom_un(df_geno_recessivo, type_group)
    
    #### Adicionar p-value a tabela
    df_chi2_rec <- fct_break_gene_variant_ends(df_chi2, "_MR")
    
    df_chi2_rec$`p-value(Chi-2)` <- round(as.numeric(df_chi2_rec$`p-value(Chi-2)`), 4)
    df_tabela_recessivo_p <- dplyr::inner_join(df_tabela_recessivo, df_chi2_rec, by = c("Gene", "variant"))
    
    ## Calculando percentual da coluna de frequencia
    df_tabela_recessivo_p_perc <- fct_calc_perc(df_tabela_recessivo_p)
    
    #### Criar tabela XLSX do df criado (modelo recessivo)
    excel_name <- paste0("tabela_modelo_recessivo_", type_group)
    fct_merge_cels(df_tabela_recessivo_p_perc, excel_name)
    
    if(teste){
      print("Finalizando geração de tabela dos modelos recessivos")
    }
    
    ####################################################################### #
    ## Gerando para o dominante
    if(teste){
      print("Iniciando geração de tabela dos modelos dominantes")
    }
    
    df_geno_dominante <- fct_remove_columns(df_modelos_genotipos, c("MA", "MR", "MU"))
    
    df_tabela_dominante <- fct_table_rec_dom_un(df_geno_dominante, type_group)
    
    #### Adicionar p-value a tabela
    df_chi2_dom <- fct_break_gene_variant_ends(df_chi2, "_MR")
    
    df_chi2_dom$`p-value(Chi-2)` <- round(as.numeric(df_chi2_dom$`p-value(Chi-2)`), 4)
    df_tabela_dominante_p <- dplyr::inner_join(df_tabela_dominante, df_chi2_rec, by = c("Gene", "variant"))
    
    ## Calculando percentual da coluna de frequencia
    df_tabela_dominante_p_perc <- fct_calc_perc(df_tabela_dominante_p)
    
    #### Criar tabela XLSX do df criado (modelo dominante)
    
    excel_name <- paste0("tabela_modelo_dominante_", type_group)
    fct_merge_cels(df_tabela_dominante_p_perc, excel_name)
    
    if(teste){
      print("Finalizando geração de tabela dos modelos dominantes")
      print("Finalizando geração de todas as tabelas de uma vez só")
    }
    ## Retornando sucesso
    return(1)
    
  }else{
    if(table == "aditivo"){
      if(teste){
        print("Iniciando geração individual de tabela dos modelos aditivos")
      }
      df_geno_aditivo <- df_modelos_genotipos |> 
        dplyr::select(dplyr::ends_with(c("PIORMB", "MO")))
      
      df_tabela_aditivo <- fct_table_ad(df_geno_aditivo, type_group)
      
      #### Adicionar p-value a tabela
      df_chi2_adit <- fct_break_gene_variant_ends(df_chi2, "_MA")
      
      df_chi2_adit$`p-value(Chi-2)` <- round(as.numeric(df_chi2_adit$`p-value(Chi-2)`), 4)
      df_tabela_aditivo_p <- dplyr::inner_join(df_tabela_aditivo, df_chi2_adit, by = c("Gene", "variant"))
      
      ## Calculando percentual da coluna de frequencia
      df_tabela_aditivo_p_perc <- fct_calc_perc(df_tabela_aditivo_p)
      
      #### Criar tabela XLSX do df criado (modelo aditivo)
      excel_name <- paste0("tabela_modelo_aditivo_", type_group)
      fct_merge_cels(df_tabela_aditivo_p_perc, excel_name)
      
      if(teste){
        print("Finalizando geração de tabela dos modelos aditivos")
      }
      ## Retornando sucesso
      return(2)
    }
    if(table == "recessivo"){
      if(teste){
        print("Iniciando geração individual de tabela dos modelos recessivos")
      }
      df_geno_recessivo <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MU"))
      
      df_tabela_recessivo <- fct_table_rec_dom_un(df_geno_recessivo, type_group)
      
      #### Adicionar p-value a tabela
      df_chi2_rec <- fct_break_gene_variant_ends(df_chi2, "_MR")
      
      df_chi2_rec$`p-value(Chi-2)` <- round(as.numeric(df_chi2_rec$`p-value(Chi-2)`), 4)
      df_tabela_recessivo_p <- dplyr::inner_join(df_tabela_recessivo, df_chi2_rec, by = c("Gene", "variant"))
      
      ## Calculando percentual da coluna de frequencia
      df_tabela_recessivo_p_perc <- fct_calc_perc(df_tabela_recessivo_p)
      
      #### Criar tabela XLSX do df criado (modelo recessivo)
      excel_name <- paste0("tabela_modelo_recessivo_", type_group)
      fct_merge_cels(df_tabela_recessivo_p_perc, excel_name)
      
      if(teste){
        print("Finalizando geração de tabela dos modelos recessivos")
      }
      ## Retornando sucesso
      return(3)
    }
    
    if(table == "dominante"){
      if(teste){
        print("Iniciando geração individual de tabela dos modelos dominantes")
      }
      df_geno_dominante <- fct_remove_columns(df_modelos_genotipos, c("MA", "MR", "MU"))
      
      df_tabela_dominante <- fct_table_rec_dom_un(df_geno_dominante, type_group)
      
      #### Adicionar p-value a tabela
      df_chi2_dom <- fct_break_gene_variant_ends(df_chi2, "_MR")
      
      df_chi2_dom$`p-value(Chi-2)` <- round(as.numeric(df_chi2_dom$`p-value(Chi-2)`), 4)
      df_tabela_dominante_p <- dplyr::inner_join(df_tabela_dominante, df_chi2_rec, by = c("Gene", "variant"))
      
      ## Calculando percentual da coluna de frequencia
      df_tabela_dominante_p_perc <- fct_calc_perc(df_tabela_dominante_p)
      
      #### Criar tabela XLSX do df criado (modelo dominante)
      
      excel_name <- paste0("tabela_modelo_dominante_", type_group)
      fct_merge_cels(df_tabela_dominante_p_perc, excel_name)
      
      if(teste){
        print("Finalizando geração de tabela dos modelos dominantes")
      }
      ## Retornando sucesso
      return(4)
    }
    
    if(table == "unico"){
      if(teste){
        print("Iniciando geração individual de tabela dos modelos unicos")
      }
      df_geno_unico <- fct_remove_columns(df_modelos_genotipos, c("MA", "MD", "MR"))
      
      df_tabela_unico <- fct_table_rec_dom_un(df_geno_unico, type_group)
      
      #### Adicionar p-value a tabela
      df_chi2_un <- fct_break_gene_variant_ends(df_chi2, "_MU")
      
      df_chi2_un$`p-value(Chi-2)` <- round(as.numeric(df_chi2_un$`p-value(Chi-2)`), 4)
      df_tabela_unico_p <- dplyr::inner_join(df_tabela_unico, df_chi2_un, by = c("Gene", "variant"))
      
      ## Calculando percentual da coluna de frequencia
      df_tabela_unico_p_perc <- fct_calc_perc(df_tabela_unico_p)
      
      #### Criar tabela XLSX do df criado (modelo unico)
      
      excel_name <- paste0("tabela_modelo_unico_", type_group)
      fct_merge_cels(df_tabela_unico_p_perc, excel_name)
      
      if(teste){
        print("Finalizando geração de tabela dos modelos unicos")
      }
      ## Retornando sucesso
      return(5)
    }
    
  }
}

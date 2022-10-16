################################################################# #
##' Script para criar tabela com saídas da reg binomial múltipla
##' 1 - retornar apenas variantes do modelo selecionado
##' 2 - Escrever tabela de resultados da binomial completa
##' 
##' Autores: Mikael e Matheus
##' Data: 08/10/2022
################################################################# #

# 0 - Scripts e bibliotecas ----

# 1 - Funções ----
## 1.1 - vars_pvalue ----
## Recebe o modelo, separa o nome das variantes em um df com colunas
## gene, variante, genótipo, gene_variante_genótipo e p-valor
fct_vars_pvalue <- function(vars, p_value, debug){
  teste_interno <- F
  if(teste_interno){
    model <- ctx_abs_or_pres_binom_mult
  }
  ## Buscando dados do objeto do modelo
  # vars <- names(model$qr$qr[2,])
  # browser()
  # p_value <- round(coef(summary(model))[,"Pr(>|z|)"], 4)
  ## Removendo constante
  vars <- vars[-c(1)]
  p_value <- p_value[-c(1)]
  # browser()
  
  ## Criando tabela de variante separado por colunas, e com última coluna para futuro join
  vars_aux <- as.data.frame(stringr::str_split_fixed(vars, "_", n = 4), stringsAsFactors = F) |> 
    dplyr::rename(Gene = 1, Variant = 2, Model = 3, dummy = 4) |> 
    # tidyr::separate(Model_dum, into = c("Model", "dummy"), sep = 2) |> 
    dplyr::mutate(gene_variant_model = paste0(Gene, "_", Variant, "_", Model))
  
  cbind(vars_aux, p_value)
}

## 1.2 -Cria dataframe usado para o join (genótipo e dummie)  ----
## Recebe o dataframe separado apenas com variantes selecionadas
## retorna um df com os genes, variantes e seus diferentes genótipos
fct_list_var_geno_dummy <- function(df_all_vars, df_vars_pvalue, debug){
  teste_interno <- F
  if(teste_interno){
    debug <- F
    df_all_vars <- df_all_vars
    df_vars_pvalue <- df_vars_pvalue
  }
  # browser()
  vars <- unique(df_vars_pvalue$gene_variant_model)
  list_dfs <- NULL
  df_full <- NULL
  for(i in 1:length(vars)){
    var <- gsub(pattern = '.{1}$', "", vars[i])
    df_aux <- df_all_vars |> 
      dplyr::select(dplyr::contains(var)) |> 
      dplyr::distinct()
    df_aux <- df_aux |>
      dplyr::select(!!as.name(paste0(var, "O")), all_of(names(df_aux))) |> 
      dplyr::rename(genotype = 1, dummy = 2) |> 
      ## Ifelse pq quero renomear as colunas, porém algumas (com final "MA" e "MD" vem antes de MO,
      ## Outras ("MR" e MU") vem depois)
      # dplyr::rename(genotype = if_else(sub('.*(?=.$)', '', vars[i], perl=T) %in% (c("A", "D"))
      #                                  , 2, 1),
      #               dummy =  if_else(sub('.*(?=.$)', '', vars[i], perl=T) %in% (c("A", "D")), 1, 2)) |>
      dplyr::mutate(gene_variant_model = paste0(vars[i]),
                    dummy = as.character(dummy))
    
    df_auxx <- dplyr::inner_join(df_vars_pvalue, df_aux, by=c("gene_variant_model", "dummy"))
    df_dummy_base <- dplyr::anti_join(df_aux, df_vars_pvalue, by=c("gene_variant_model", "dummy")) |> 
      dplyr::rename(dummy_base = dummy, genotype_base = genotype)
    
    df_auxxx <- dplyr::left_join(df_auxx, df_dummy_base, by=c("gene_variant_model")) |> 
      dplyr::select(Gene, Variant, Model, dummy_base, genotype_base, dummy, genotype,
                    p_value)
    
    # ## i <- 3 problemático porque temos mais de dois níveis, lidar com ele
    # if(nrow(df_aux) != nrow(df_aux_t)){
    #   gene_variant_model <- c(paste0(var, "O"), rep(paste0(var, "A"), 2))
    # }
    df_full <- rbind(df_full, df_auxxx)
  }
  df_full
}


## 1.3 - Resultado da binomial múltipla ----
## Vai receber o df do protocolo selecionado, o modelo e o nome do modelo
## Fará primeiro o join da tabela já filtrada com os genótipos e suas dummies
## Para ter o gene base + a variação, e pode escrever isso na tabela de resultados
fct_res_bin_mult <- function(df, vars, p_value, model_name, switch_write_table){
  teste_interno <- F
  if(teste_interno){
    model <- ctx_abs_or_pres_binom_mult
    df <- df_ctx
    model_name <- "ctx_abs_or_pres_binom_mult"
    # model_
  }
  ## Recebendo tabela \
  # browser()
  
  df_vars_pvalue <- fct_vars_pvalue(vars, p_value, debug = F)
  rownames(df_vars_pvalue) <- NULL
  
  ## Separando nomes de coluna, coluna com nome do genótipo e com valor da dummy
  vars_model <- unique(df_vars_pvalue$gene_variant_model)
  vars_model_genotype <- gsub(pattern = '.{1}$', "O", vars_model)
  vars_selected <- stringr::str_sort(c(vars_model, vars_model_genotype))
  
  ## Criando dfs (um pouco simplificado) com variantes (genótipos e modelos usados)
  df_all_vars <- df |> 
    dplyr::select(all_of(vars_selected)) |> 
    dplyr::distinct(across(all_of(vars_selected)))
  
  ## Separando dfs por variante
  df_bin_mult <- fct_list_var_geno_dummy(df_all_vars, df_vars_pvalue, debug = F)
  
  # df_teste <- as.data.frame(t(df_all_vars)) |> 
  #   tibble::rownames_to_column("gene_variant_model")
  
  ## Criando nome de protocolo e agrupamento para nome do arquivo excel
  model_name_ <- gsub("_binom_mult", "", model_name)
  model_name_ <- stringr::str_split_fixed(model_name_, "_", n = 2)
  protocol <- model_name_[1]
  group <- model_name_[2]
  if(switch_write_table){
    if(!dir.exists("data-raw/tabela-resultados-bin/")){
      dir.create("data-raw/tabela-resultados-bin/")
    }
    path_write <- paste0("data-raw/tabela-resultados-bin/tabela_bin_mult_", 
                         protocol, "_", group, ".xlsx")
    writexl::write_xlsx(df_bin_mult, path_write)
  }
}


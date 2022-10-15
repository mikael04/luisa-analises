################################################################# #
##' Script para criar tabela com saídas da reg binomial múltipla
##' 1 - retornar apenas variantes do modelo selecionado
##' 2 - Escrever tabela de resultados da binomial completa
##' 
##' Autores: Mikael e Matheus
##' Data: 08/10/2022
################################################################# #

# 0 - Scripts e bibliotecas ----

# 1 - Funções
## 1.1 - vars_pvalue ----
## Recebe o modelo, separa o nome das variantes em um df com colunas
## gene, variante, genótipo, gene_variante_genótipo e p-valor
fct_vars_pvalue <- function(model, debug){
  teste_interno <- F
  if(teste_interno){
    model <- ctx_abs_or_pres_binom_mult
  }
  ## Buscando dados do objeto do modelo
  vars <- names(model$qr$qr[2,])
  p_value <- round(coef(summary(model))[,"Pr(>|z|)"], 4)
  ## Removendo constante
  vars <- vars[-c(1)]
  p_value <- p_value[-c(1)]
  
  ## Criando tabela de variante separado por colunas, e com última coluna para futuro join
  vars_aux <- as.data.frame(stringr::str_split_fixed(vars, "_", n = 3), stringsAsFactors = F) |> 
    dplyr::rename(Gene = 1, Variant = 2, Model_dum = 3) |> 
    tidyr::separate(Model_dum, into = c("Model", "dummy"), sep = 2) |> 
    dplyr::mutate(gene_variant_model = paste0(Gene, "_", Variant, "_", Model))
  
  cbind(vars_aux, p_value)
}

## 1.2 -Cria dataframe usado para o join (genótipo e dummie)  ----
## Recebe o dataframe separado apenas com variantes selecionadas

## 1.3 - Resultado da binomial múltipla ----
## Vai receber o df do protocolo selecionado, o modelo e o nome do modelo
## Fará primeiro o join da tabela já filtrada com os genótipos e suas dummies
## Para ter o gene base + a variação, e pode escrever isso na tabela de resultados
fct_res_bin_mult <- function(df, model, model_name, switch_write_table){
  if(teste_interno){
    model <- ctx_abs_or_pres_binom_mult
    df <- df_ctx
    model_name <- "ctx_abs_or_pres_binom_mult"
    model_
  }
  ## Recebendo tabela 
  df_vars_pvalue <- fct_vars_pvalue(model, debug = F) |> 
    `rownames<-`(NULL)
  
  ## Separando nomes de coluna, coluna com nome do genótipo e com valor da dummy
  vars_model <- unique(df_vars_pvalue$gene_variant_model)
  vars_model_genotype <- gsub(pattern = '.{1}$', "O", vars_model)
  vars_selected <- stringr::str_sort(c(vars_model, vars_model_genotype))
  
  ## Criando dfs (um pouco simplificado) com variantes (genótipos e modelos usados)
  df_all_vars <- df |> 
    dplyr::select(all_of(vars_selected)) |> 
    dplyr::distinct(across(all_of(vars_selected)))
  
  ## Separando dfs por variante
  num_vars <- length(unique(vars_model))
  
  df_teste <- as.data.frame(t(df_all_vars)) |> 
    tibble::rownames_to_column("gene_variant_model")
  
  ## Criando nome de protocolo e agrupamento para nome do arquivo excel
  model_name <- gsub("_binom_mult", "", model_name)
  model_name <- stringr::str_split_fixed(model_name, "_", n = 2)
  protocol <- model_name_[1]
  group <- model_name_[2]
  if(switch_write_table){
    if(!file.exists("data-raw/tabela-resultados-bin/")){
      file.create("data-raw/tabela-resultados-bin/")
    }
    path_write <- paste0("data-raw/tabela-resultados-bin/tabela_bin_mult_", 
                         protocol, agroup, ".xlsx")
    writexl::write_xlsx(tabela_resultados, path_write)
  }
}


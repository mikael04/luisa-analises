######################################################### #
##' Script para avaliar o retorno das funções usadas
##' 
##' Autores: Mikael e Matheus
##' Data: 08/07/2022
######################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_check_return <- function(return, function_name, table, type_group, switch_teste){
  if(function_name == "create_table"){
    if(return == 1){
      if(type_group == "pres"){
        print("Geração de todas as tabelas para parâmetro presença vs ausência feitas com sucesso")
      }
      if(type_group == "ulc"){
        print("Geração de todas as tabelas para parametro ulcerados vs não ulcerados feitas com sucesso")
      }
      if(type_group == "sev"){
        print("Geração de todas as tabelas para parametro mb severo vs não severo feitas com sucesso")
      }
      
    }
    if(return > 1 & return <= 5){
      print("Geração de tabela individual feita com sucesso")
      print(paste0("Gerada apenas tabela para modelo ", table))
    }
  }
}
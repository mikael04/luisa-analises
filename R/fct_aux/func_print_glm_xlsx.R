######################################################### #
##' Script para organizar e imprimir arquivo de resposta da modelagem glm
##' 
##' Autores: Mikael e Matheus
##' Data: 08/10/2022
######################################################### #

# 0 - Scripts e bibliotecas ----

# 1 - Função ----
fct_print_glm_xlsx <- function(model, switch_write_table){
  model_name <- deparse(substitute(model))
  vars <- names(model$qr$qr[2,])
  estimate <- round(coef(summary(model))[,"Estimate"], 4)
  z_value <- round(coef(summary(model))[,"z value"], 4)
  p_value <- round(coef(summary(model))[,"Pr(>|z|)"], 4)
  
  tabela <- data.frame(vars, estimate, z_value, p_value)
  
  if(switch_write_table){
    writexl::write_xlsx(tabela, paste0("data-raw/tabela-reg-bin/resultado_", model_name, ".xlsx"))
  }
}

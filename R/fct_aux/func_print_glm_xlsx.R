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
  vars <- gsub(".{1}$", "", vars)
  vars[1] <- "(Intercept)"
  estimate <- round(coef(summary(model))[,"Estimate"], 4)
  z_value <- round(coef(summary(model))[,"z value"], 4)
  odd.ratio <- round(exp(estimate), 2)
  p_value <- round(coef(summary(model))[,"Pr(>|z|)"], 4)
  df_confint <- as.data.frame(exp(confint(model))) |> 
    dplyr::mutate(`2.5 %` = round(`2.5 %`, 4),
                  `97.5 %` = round(`97.5 %`, 4))
  
  tabela <- data.frame(vars, estimate, odd.ratio, df_confint, z_value, p_value) |> 
    dplyr::rename(`2.5%` = `X2.5..`, `97.5%` = `X97.5..`)
  
  if(switch_write_table){
    writexl::write_xlsx(tabela, paste0("data-raw/tabela-reg-bin/resultado_", model_name, ".xlsx"))
  }
}

######################################################### #
##' Script de análise exploratória dos dados
##' 
##' Autores: Mikael e Matheus
##' Data: 23/06/2022
######################################################### #

# 0 - Scripts e bibliotecas ----
library(dplyr)


# 1 - Lendo Base de dados ----
df_full <- haven::read_sav("data-raw/banco lla e linfoma 16.05.sav")

df <- df_full |> 
  dplyr::select(PIORMB, dplyr::ends_with(c("MA", "MR", "MD")))

############################################################ #
## Aprendendo a usar Chi-quadrado
data("HairEyeColor")
df_teste <- as.data.frame(HairEyeColor) 

HairEyeColor[,,1]

HairEyeColor
chisq.test(HairEyeColor[,,1])

chisq.test(df_teste$Sex == 'Male')
df_teste
df_teste_male = df_teste[df_teste$Sex == 'Male',]
table(df_teste_male)
table()

table(df_teste_male$Hair, df_teste_male$Eye)

library(MASS)        
print(str(survey))

survey_df <- survey
# Create a data frame from the main data set.
stu_data = data.frame(survey$Smoke,survey$Exer)

# Create a contingency table with the needed variables.           
stu_data = table(survey$Smoke,survey$Exer)
############################################################ #

table_ABCC2_RS2273697 <- table(df$PIORMB, df$ABCC2_rs2273697_MA)
chisq.test(table_ABCC2_RS2273697)

#create data frame with 0 rows and 3 columns
df_chi2 <- data.frame(matrix(ncol = 2, nrow = 0))
#provide column names
colnames(df_chi2) <- c("Variante", "Chi-2")

for(i in range(2:160)){
  testes_chi2 <- chisq.test(table(unlist(df[,1]), unlist(df[,i])), simulate.p.value = TRUE)
  colnames_df[i]
  testes_chi2$p.value
  df_chi2[i-1, ] = c(colnames_df[i], testes_chi2$p.value)
  print(i)
}

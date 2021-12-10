# libraries
library(readr)
library(sjPlot)

# working directory

Nosema_Proteins <- read_csv("~/Documents/Random/Nosema_Proteins.csv")

tab_df(Nosema_Proteins, alternate.rows = TRUE,
       title = "Nosema_Proteins",
       file = "Nosema_Proteins.doc")

getwd()
setwd("/Users/jaredflater/Documents/Ederson/")
library(readxl)
Accession_numbers_for_Jared_and_Adina <- read_excel("~/Documents/Ederson/Accession numbers for Jared and Adina.xlsx")
View(Accession_numbers_for_Jared_and_Adina)

df <- Accession_numbers_for_Jared_and_Adina$Refseq

write.(df, "IDs.txt", sep="\t")


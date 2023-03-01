
install.packages("devtools")

devtools::install_local("refseqR.zip")



library(devtools)
library(refseqR)
library(tidyr)
library(dplyr)
library(readr)
library(geneHummus)

Blast1 <- read.csv("XP_Loc_unici_lista_Mads.csv")

IDs <- Blast1$list

#or you can use 

 IDs <- (Blast1 %>% pull("Accession"))
 IDs

 
 
 save_AAfasta_from_xps(IDs, "ilaria_tenta.txt")

 
 




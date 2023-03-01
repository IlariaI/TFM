library(devtools)
library(refseqR)
library(tidyr)
library(dplyr)
library(readr)
library(geneHummus)

Blast1 <- read.csv("XP_Mads_unici.csv")

IDs <- Blast1$list

#or you can use 

 IDs <- (Blast1 %>% pull("Accession"))
 IDs

 
 
 save_AAfasta_from_xps(IDs, "Alignment_multifasta.txt")

 
 




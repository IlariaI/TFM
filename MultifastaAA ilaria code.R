library(devtools)
library(refseqR)
library(tidyr)
library(dplyr)
library(readr)
library(geneHummus)
devtools::source_url("https://github.com/jdieramon/my_scripts/blob/master/mads_box/functions.R?raw=TRUE")


              #  Iterative BLAST 

# Define the 'hitFile' variable with a file path
hitFile <- "data/chickpea_hit_table.csv"

# Perform BLAST best homolog function
res <- blast_best_homolog(hitFile)

# Extract the 'subject' column f
ca_mads <- res %>% pull(subject)

# Remove duplicates from 'xps' and store the unique XP ids
xps <- unique(xps)

# Add the first XP id from 'res' to 'ca_MADSbox'
ca_MADSbox <- c(ca_mads, res$query[1])


#checks the e value and identity score of the obtained HITS

# Create a histogram
hist(res$evalue)

hist(res$identity)

# Filter 'res' to select rows with identity >= 70 
ids70 <- res %>% 
  filter(identity >=70) %>% 
  pull(subject)

# Check if the 'ids70' values are in 'xps' and if not display them
ids70 %in% xps

ids70[!ids70 %in% xps]
#Save the newly obtaned results and add them to the main list
round_2 <- ids70[!ids70 %in% xps]
              
#Once enstablished the optimal paramiter, will be possible to set them to automatically filter the result of the BLAST analysis 
blast_res <- blast_homologs(hitFile, ident = 55, cover = 90)


# Addin new HITS to the database
ca_MADSbox <- unique(c(xps, blast_res$subject))              
#repeat until no new results are retrived

# Get the mRNA sequence 


# Extract the 'list'
ca_mads <- ca_mads$list

# Get XM ids from XP ids using the 'getXM' function
xms <- sapply(xps, function(xp) getXM(xp), USE.NAMES = FALSE)

# Create a list containing the corresponding XM and XP
seqs <- list(unname(xms), xps)

# Save the 'seqs' object to a file named "sequences.rda"
save(seqs, file = "data/sequences.rda")

# Make a multi-fasta file from XM ids
save_CDSfasta_from_xms(xms, nameFile = "data/mads_cds")


              #Characterization table

# Define the list of genes we want to characterize
targets <-ca_mads

# Characterize the gene family table
char_w_iso <- characterizeTable(targets)

              # Sort the data set by chromosome and start coordinate
tdat %>% 
  arrange(Chr, chr_s)

# Create a data set based on 1 gene model / locus
tdat %>% 
  arrange(Chr, chr_s) %>% 
  distinct(LOC, .keep_all = TRUE)

# Get information on isoforms: number of loci with more than 1 protein
tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC, sort = TRUE) %>% 
  filter(n > 1)

# Count the number of isoforms for each LOC
nisof <- tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC) %>% 
  pull(n)

# Update 'tdat' with the count of isoforms
tdat <- tdat %>%  
  arrange(LOC, desc(AA)) %>% 
  distinct(LOC, .keep_all = TRUE) %>% 
  mutate(n_isof = nisof)

# Extract the LOC column from 'char_w_iso' and store it in 'loc'
loc <- Char_w_iso %>%  
    as_tibble() %>% 
    pull(LOC)

# Create a new table 'unique_char_tab' by grouping 'char_w_iso' by LOC and selecting the first entry in each group (the biggest isoform)
unique_char_tab <- char_w_iso %>%  
  as_tibble() %>% 
  group_by(LOC) %>%
  dplyr::slice(1) %>% 
  ungroup()
              

 
 




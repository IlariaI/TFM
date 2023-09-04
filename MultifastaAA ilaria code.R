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
tdat <- characterizeTable(targets)

# Sort the data set by chromosome and start coordinate
tdat %>% 
  arrange(Chr, chr_s)

# Create a data set based on the model 1 gene / locus
tdat %>% 
  arrange(Chr, chr_s) %>% 
  distinct(LOC, .keep_all = TRUE)

# number of loci with more than 1 protein
tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC, sort = TRUE) %>% 
  filter(n > 1)

# Count the number of isoforms for each LOC
nisof <- tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC) %>% 
  pull(n)

# Update with the count of isoforms
tdat <- tdat %>%  
  arrange(LOC, desc(AA)) %>% 
  distinct(LOC, .keep_all = TRUE) %>% 
  mutate(n_isof = nisof)

# Write the gene family table to a CSV file
write.csv(tdat, file = "data/char_table.csv")


#                    GRanges object

# Negative widths are not allowed by IRanges, so let's define the true 
# coordinates and build a new data frame, which is the input for the function 
# 'makeGRangesFromDataFrame'

BiocManager::install("GenomicRanges")
library(GenomicRanges)

# Create a data frame 'gr' by transforming and modifying columns from 'tdat'
gr <- tdat %>% 
  transmute(LOC, Chr, Strand, 
            Coordstart = ifelse(Strand == "-", chr_e, chr_s),
            Coordend = ifelse(Strand == "+", chr_e, chr_s))

# Add a new column 'bp_cut' to map the TSS coordinate for promoter analysis
# Customize 'bp_cut' based on specific LOC values
gr <- gr %>% 
  dplyr::mutate(bp_cut = 150) %>% 
  dplyr::mutate(bp_cut = ifelse(LOC == "LOC101493118", 14, bp_cut), 
                bp_cut = ifelse(LOC == "LOC101504656", 27, bp_cut), 
                bp_cut = ifelse(LOC == "LOC101509413", 3,  bp_cut), 
                bp_cut = ifelse(LOC == "LOC101510075", 39, bp_cut))

# Create a GRanges object from the modified 'gr' data frame
gr <- makeGRangesFromDataFrame(gr, start.field = "Coordstart", 
                               end.field = "Coordend", 
                               strand.field = "Strand", ignore.strand = FALSE, 
                               seqnames.field = "Chr", 
                               keep.extra.columns = TRUE)

# Add genome information to the GRanges object
genome(gr) = "ASM33114v1"

# Filter out unplaced sequences 
gr <- gr[seqnames(gr) != "Un"]


# Order the 'GRanges' object by chromosome and region location
sort(table(seqnames(gr)), decreasing = TRUE)
gr <- sort(gr)
              #     Promoter sequence 


library(Biostrings)


# Load mRNA sequences and convert them into a DNAStringSet
mads_cds = readDNAStringSet("data/mads_cds.fasta")

# Convert sequence names to LOC names
my_names = names(mads_cds)
my_names <- sapply(my_names, function(i) names2LOC(i), USE.NAMES = FALSE)
names(mads_cds) = my_names

# Load the Cicer arietinum genome from the 'BSgenome' package
library(BSgenome.Carietinum.NCBI.v1)
genome = BSgenome.Carietinum.NCBI.v1

# Create a GRangesList object containing TSS coordinates for each chromosome
gr.tss = GRangesList(sapply(paste0("Ca", 1:8), function(i) TSScoordinates(gr, mads_cds, i)))
gr.tss = unlist(gr.tss)

# Save objects including 'seqs', 'gr', and 'gr.tss'
save(seqs, gr, gr.tss, file = "res/sequences.rda")

# Extract promoters (1500 nt upstream and 9 nt downstream of TSS)
gr.promoters = promoters(gr.tss, upstream=1500, downstream=9)

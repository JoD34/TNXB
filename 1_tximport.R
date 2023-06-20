# Import libraries ----
library(AnnotationHub)
library(rhdf5)
library(stringr)
library(tidyverse)
library(tximport)

# Get .h5 paths ----

files.cwd <- list.files(pattern = "^SRR.*$")
paths <- file.path(files.cwd,"abundance.h5")

# Verify the existence of the paths
if (!all(file.exists(paths))) stop("All given paths aren't working!")

# Get column names for DESeq2
col.names <- c("TNXB_LF_1", "TNXB_LF_2", "TNXB_WT_1", "TNXB_WT_2")

# Annotation ----

# Retrieve db for annotation
hub <- AnnotationHub()
req <- query(hub, "EnsDb.Hsapiens.v103")
ensdb <- hub[[names(req)]]

# get annotations 
suppressWarnings({
  tx <- transcripts(ensdb, columns = c("tx_id", "gene_name")) %>% 
    as_tibble() %>% 
    dplyr::rename(target_id = tx_id) %>% 
    dplyr::select(target_id, gene_name)
})

# Import kallisto .h5 ----
txi.kallisto <- tximport(paths, 
                         type = "kallisto", 
                         tx2gene = tx, 
                         ignoreTxVersion = T)

colnames(txi.kallisto$abundance) <- col.names
colnames(txi.kallisto$counts) <- col.names
colnames(txi.kallisto$length) <- col.names

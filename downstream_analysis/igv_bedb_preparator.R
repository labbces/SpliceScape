library(dplyr)
library(tidyr)
library(RSQLite)
library(DBI)
library(stringr)

options(scipen = 999)

species_of_interest <- "Oryza sativa" 
output_bed_file <- paste0("all_events_", gsub(" ", "_", species_of_interest), ".bed")
db_path <- "/home/bia/LandscapeSplicingGrasses/data/grasses.db"

con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# 1. Carregamos os dados 
splicing_data <- tbl(con, "splicing_events") %>%
  filter(species == species_of_interest) %>%
  select(event_id, seqid, strand, full_coord, event_type, mean_psi_majiq) %>% 
  collect()

dbDisconnect(con)

# 2. Processamento
bed_data <- splicing_data %>%
  filter(!is.na(full_coord), full_coord != "") %>%
  
  mutate(
    is_RI = str_detect(event_type, "RI|intron") & !str_detect(full_coord, "-"),
    coords_extracted = if_else(
      is_RI,
      str_extract(full_coord, ":(\\d+),(\\d+)", group = c(1, 2)),
      str_extract(full_coord, ":(\\d+),(\\d+)-(\\d+),(\\d+)", group = c(1, 2, 3, 4))
    )
  ) %>%
  
  separate(full_coord, into = c("discard_seq", "raw_coords"), sep = ":", remove = FALSE) %>%
  
  mutate(
      c1 = as.numeric(str_match(raw_coords, "^(\\d+),")[,2]),
      c2 = as.numeric(str_match(raw_coords, ",(\\d+)[,-]")[,2]), 
      c3 = as.numeric(str_match(raw_coords, "-(\\d+),")[,2]),    
      c4 = as.numeric(str_match(raw_coords, ",(\\d+)$")[,2])     
  ) %>%
  
  transmute(
    chrom = seqid,
    chromStart = c1 - 1, 
    chromEnd = if_else(is_RI, c2, c4),
    name = event_id,
    
    # --- CORREÇÃO DO PSI AQUI ---
    # Transforma 0.5 (50%) em 500 para o padrão BED
    score = as.integer(coalesce(mean_psi_majiq, 0) * 1000),
    
    strand = strand,
    thickStart = chromStart,
    thickEnd = chromEnd,
    itemRgb = "0,0,0", 
    
    blockCount = if_else(is_RI, 1, 2),
    
    blockSizes = if_else(
      is_RI,
      as.character(c2 - c1 + 1),
      paste((c2 - c1 + 1), (c4 - c3 + 1), sep = ",")
    ),
    
    blockStarts = if_else(
      is_RI,
      "0",
      paste("0", (c3 - c1), sep = ",")
    )
  ) %>%
  filter(!is.na(chromStart), !is.na(chromEnd))

write.table(
  bed_data, 
  output_bed_file, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
)

print(paste("BED file com PSI Scores gerado:", output_bed_file))
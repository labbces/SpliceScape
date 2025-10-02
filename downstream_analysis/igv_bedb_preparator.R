library(dplyr)
library(tidyr)
library(RSQLite)
library(DBI)
library(stringr)


species_of_interest <- "Oryza sativa" #### update
output_bed_file <- paste0("all_events_", gsub(" ", "_", species_of_interest), ".bed")

db_path <- "/home/bia/LandscapeSplicingGrasses/data/grasses.db"
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

splicing_data <- tbl(con, "splicing_events") %>%
  filter(species == species_of_interest) %>%
  collect()

dbDisconnect(con)

# Format BED12 ---
bed_data <- splicing_data %>%
  filter(
    !is.na(upstream_exon_coord) & 
    !is.na(downstream_exon_coord) &
    str_detect(upstream_exon_coord, "^\\d+-\\d+$") &
    str_detect(downstream_exon_coord, "^\\d+-\\d+$")
  ) %>%
  
  separate(upstream_exon_coord, into = c("up_start", "up_end"), sep = "-", convert = TRUE) %>%
  separate(downstream_exon_coord, into = c("down_start", "down_end"), sep = "-", convert = TRUE) %>%
  
  transmute(
    chrom = seqid,
    chromStart = up_start - 1,
    chromEnd = down_end,
    name = event_id,
    score = 0,
    strand = strand,
    thickStart = up_start - 1,
    thickEnd = down_end,
    itemRgb = 0,
        
    # Config to deal with reteined and spliced introns
    blockCount = if_else(str_detect(event_id, "RI_intron"), 1, 2),
    
    blockSizes = if_else(
      str_detect(event_id, "RI_intron"),
      as.character(down_end - (up_start - 1)), 
      paste(up_end - up_start, down_end - down_start, sep = ",")
    ),
    
    blockStarts = if_else(
      str_detect(event_id, "RI_intron"),
      "0",
      paste(0, down_start - (up_start - 1) -1, sep = ",")
    )
  )

write.table(
  bed_data, 
  output_bed_file, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE
)

print(paste("BED file to", species_of_interest, "successfully created:", output_bed_file))
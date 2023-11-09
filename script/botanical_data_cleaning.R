# Botanical data prep / cleaning
# Jakob J. Assmann 14 August 2023

# Dependencies
library(tidyverse)
library(readxl)

# Load XLS spread sheet
botanical_obs <- read_xlsx("data/botanical_obs/2023-08-11_botanical_ops_final.xlsx") %>%
    filter(Cleaned == "Yes") %>%
    select(
        date = "Date",
        observer = "Observer",
        scribe = "Scribe",
        site = "Site",
        plot_id = "Plot_ID",
        species = "Full name"
    ) %>%
    mutate(plot_id = paste0("S", gsub("Site ", "", site), plot_id)) %>%
    filter(!grepl("NA*", species)) %>%
    filter(!is.na(species))

# Check unique species
distinct(botanical_obs, species) %>% print(n = 100)

# Remove all entries without a species id, but keeping cf's
n_removed <- c("Eudicot", "Monocot", "sp\\. [1-9]") %>%
    map(function(x) length(grepl(x, botanical_obs$species))) %>%
    unlist() %>%
    sum()
botanical_obs <- botanical_obs %>%
    filter(!grepl("Eudicot", species)) %>%
    filter(!grepl("Monocot", species)) %>%
    filter(!grepl( "sp\\. [1-9]", species))

# Check unique species again
distinct(botanical_obs, species) %>% print(n = 100)

# write out cleaned data file
write_csv(botanical_obs, "data/botanical_obs/presence_absence_cleaned.csv")

# Load point frame data
point_frame_obs <- read_xlsx("data/botanical_obs/2023-08-11_botanical_ops_final.xlsx", sheet = "Pointframing")
write_csv(botanical_obs, "data/botanical_obs/pointframing.csv")

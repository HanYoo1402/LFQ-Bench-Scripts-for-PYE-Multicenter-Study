###############################################################################
# Protein Rank Plot — Intensity-ranked proteins with log10 normalized values
#
# Description:
# Processes .tsv report files (DIA) or .txt report files (DDA), merges with
# design metadata, normalizes intensities within PYE, computes per-protein 
# mean normalized intensity per dilution, ranks proteins by intensity within
# each dilution, and plots protein rank vs. log10 normalized intensity.
#   
# Usage:
#   1. Update the folder_path, save_path, and design_file with your directories.
#   2. Run the script in R. A PDF figure will be saved in the output directory.
#
# Requirements:
#   - R (4.3.2) (used version for generating plots)
#   - Packages: data.table, ggplot2
###############################################################################

library(data.table)
library(ggplot2)

generate_protein_rank_plot <- function(folder_path, design_path, save_path) {
  # List only TSV files (DIA), for DDA data use pattern = "\\.txt" and correct column header
  files <- list.files(folder_path, pattern = "\\.tsv", full.names = TRUE)
  
  # Required columns from report files
  selected_columns <- c("Protein.Names", "PG.MaxLFQ", "Run")
  
  dfs <- vector("list", length(files))
  names(dfs) <- files
  for (file in files) {
    dfs[[file]] <- fread(file, select = selected_columns)
  }
  
  # Combine into one data.table
  combined_df <- rbindlist(dfs, use.names = TRUE, fill = TRUE)
  
  # Read design and merge on Run
  design <- fread(design_path)
  merged_data <- merge(combined_df, design, by = "Run")
  
  # Mean PG.MaxLFQ per protein/run/PYE/dilution
  merged_data <- merged_data[
    , .(PG.MaxLFQ = mean(PG.MaxLFQ, na.rm = TRUE)),
    by = .(Protein.Names, Run, PYE, Dilution)
  ]
  
  # Normalize within PYE
  merged_data[, Normalized_Intensity := PG.MaxLFQ / max(PG.MaxLFQ, na.rm = TRUE), by = PYE]
  
  # Mean normalized intensity per protein & dilution
  merged_data_unique <- merged_data[
    , .(Intensity = mean(Normalized_Intensity, na.rm = TRUE)),
    by = .(Protein.Names, Dilution)
  ]
  
  # Avoid log10(0) and ensure finite values
  merged_data_unique <- merged_data_unique[is.finite(Intensity) & Intensity > 0]
  
  # Rank proteins within each dilution (higher intensity → smaller rank)
  merged_data_unique[, Rank := frank(-Intensity, ties.method = "min"), by = Dilution]
  
  # Color palette mapping (example HUMAN)
  # c("#fd862a", "#d95f02", "#8d3e01") HUMAN
  # c("#2bdba6", "#1b9e77", "#105d46") ECOLI
  # c("#a5a2ce", "#7570b3", "#4f4a8c") YEAST
  pal <- c("PYE1" = "#8d3e01", "PYE3" = "#d95f02", "PYE9" = "#fd862a")
  
  # Plot (adjust sizes and labels as needed)
  p1 <- ggplot(merged_data_unique, aes(x = Rank, y = log10(Intensity), color = Dilution)) +
    geom_point(size = 2) +
    labs(x = "Protein Rank", y = "Log10 Normalized Intensity", title = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 40),
      axis.text  = element_text(size = 40),
      strip.text = element_text(size = 40),
      legend.position = "none"
    ) +
    scale_color_manual(values = pal) +
    ylim(NA, 0)  # normalized ∈ (0,1] ⇒ log10 ≤ 0
  
  # Save to PDF (adjust width/height as needed)
  ggsave(
    p1,
    path = save_path,
    filename = "output.pdf",
    width = 7440, height = 5262, units = "px"
  )
}

# ----------------------- Fill in your paths and run ---------------------------
folder_path <- "<FILL_IN_PATH>/data" 
design_path <- "<FILL_IN_PATH>/design.csv"  
save_path   <- "<FILL_IN_PATH>/output"

generate_protein_rank_plot(folder_path, design_path, save_path)

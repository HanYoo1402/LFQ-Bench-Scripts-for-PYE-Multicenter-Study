###############################################################################
# Spaghetti plot: dilution series per protein group (ECOLI)
#
# Description:
#   For each dataset (base name), read PYE1/PYE3/PYE9 DIA reports, merge with
#   the design table, apply FDR filters and species filter (ECOLI), collapse to
#   a single log10(PG.MaxLFQ) per Protein.Group × Dilution × Condition, enforce
#   monotonic decrease across dilutions, bin proteins into 10 deciles by PYE1
#   intensity, and plot all protein trajectories with decile medians overlaid.
#   Same script structure can be applied to DDA reports —
#   just replace input files and adjust the column headers
#   to match the DDA report format.
#
# Inputs (edit in “Paths” below):
#   • dir1/dir3/dir9: folders containing *_PYE1.tsv, *_PYE3.tsv, *_PYE9.tsv
#   • design_path: design table with columns {Run, Dilution, Condition, ...}
#   • save_path: output folder for PDFs and summary CSV
#
# Outputs:
#   • For each base name:  spaghetti_<base>.pdf
#   • Summary CSV:         protein_filter_summary.csv
#
# Requirements:
#   • R ≥ 4.3.2
#   • Packages: data.table, dplyr, ggplot2, stringr
#
# Notes:
#   • “Valid” proteins must show a non-increasing trend across the ordered
#     levels: PYE1_A, PYE1_B, PYE3_A, PYE3_B, PYE9_A, PYE9_B.
#   • Deciles are defined on PYE1_PYE_A intensities (one bin per 10%).
#   • X-axis tick labels reflect Dilution_Condition; x-scale is reversed.
###############################################################################

# Load libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

generate_protein_dilution_plot <- function(file_paths, design_path, save_path,
                                           plot_title, output_name, summary_list) {
  selected_columns <- c(
    "Protein.Names","Protein.Group","PG.MaxLFQ","Run",
    "Stripped.Sequence","Q.Value","PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value"
  )
  
  # Read and bind all files
  dfs <- lapply(file_paths, function(f) data.table::fread(f, select = selected_columns))
  combined_df <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)
  
  design <- data.table::fread(design_path)
  
  merged_data <- merge(combined_df, design, by = "Run") %>%
    dplyr::filter(
      Q.Value      <= 0.01,
      PG.Q.Value   <= 0.01,
      Lib.Q.Value  <= 0.01,
      Lib.PG.Q.Value <= 0.01,
      PG.MaxLFQ > 0,
      !stringr::str_detect(Protein.Names, "_CONTA"),
      stringr::str_detect(Protein.Names, "ECOLI")
    ) %>%
    data.table::setDT()
  
  original_protein_count <- data.table::uniqueN(merged_data$Protein.Group)
  
  # Collapse replicates to one log10 mean per protein/dilution/condition
  merged_data <- unique(merged_data, by = c("Protein.Group","PG.MaxLFQ","Run"))
  merged_data <- merged_data[
    , .(Mean_PG.MaxLFQ = mean(log10(PG.MaxLFQ), na.rm = TRUE)),
    by = .(Protein.Group, Dilution, Condition)
  ]
  
  # Ordered factor for dilution × condition
  merged_data[
    , Dilution_Condition := factor(
      paste(Dilution, Condition, sep = "_"),
      levels = c("PYE1_PYE_A","PYE1_PYE_B","PYE3_PYE_A","PYE3_PYE_B","PYE9_PYE_A","PYE9_PYE_B")
    )
  ]
  data.table::setorder(merged_data, Protein.Group, Dilution_Condition)
  
  # Keep proteins with monotonic non-increasing trend across the sequence
  check_trend <- function(vals) all(diff(vals) <= 0)
  valid_proteins <- merged_data[, .SD[check_trend(Mean_PG.MaxLFQ)], by = Protein.Group]
  
  remaining_protein_count <- data.table::uniqueN(valid_proteins$Protein.Group)
  removed_count <- original_protein_count - remaining_protein_count
  
  summary_list[[output_name]] <- data.table::data.table(
    Plot = output_name,
    Original_Proteins  = original_protein_count,
    Removed_Proteins   = removed_count,
    Remaining_Proteins = remaining_protein_count
  )
  
  if (remaining_protein_count == 0) {
    message("Skipping ", output_name, ": no valid proteins remain.")
    return(summary_list)
  }
  
  # Bin by PYE1_PYE_A deciles
  pye1_data <- valid_proteins[Condition == "PYE_A" & Dilution == "PYE1"]
  pye1_data <- pye1_data[order(-Mean_PG.MaxLFQ)]
  pye1_data[, Group := cut(seq_len(.N), breaks = 10, labels = FALSE)]
  valid_proteins <- merge(valid_proteins, pye1_data[, .(Protein.Group, Group)],
                          by = "Protein.Group", all.x = TRUE)
  
  # Bin-wise summaries (mean and median) per Dilution_Condition
  summary_data <- valid_proteins[!is.na(Group),
                                 .(Mean = mean(Mean_PG.MaxLFQ, na.rm = TRUE),
                                   Median = median(Mean_PG.MaxLFQ, na.rm = TRUE)),
                                 by = .(Group, Dilution_Condition)
  ]
  
  # X mapping for ordered dilution × condition (reversed scale)
  x_values <- c(
    "PYE1_PYE_A" = log10(8),   "PYE1_PYE_B" = log10(4),
    "PYE3_PYE_A" = log10(2.67),"PYE3_PYE_B" = log10(1.33),
    "PYE9_PYE_A" = log10(0.89),"PYE9_PYE_B" = log10(0.44)
  )
  valid_proteins[,  Dilution_Condition_Num := x_values[as.character(Dilution_Condition)]]
  summary_data[, Dilution_Condition_Num := x_values[as.character(Dilution_Condition)]]
  
  # Plot
  # White strip used to "skip" unused y-ranges (adjust for setups) and shrink whitespace 
  # so relevant intensity ranges are emphasized.
  p <- ggplot(valid_proteins,
              aes(x = Dilution_Condition_Num, y = Mean_PG.MaxLFQ, group = Protein.Group)) +
    geom_line(linewidth = 1, alpha = 0.05) +
    geom_point(size = 0.5, alpha = 0.1) +
    geom_point(data = summary_data,
               aes(x = Dilution_Condition_Num, y = Median),
               color = "red", size = 6, inherit.aes = FALSE) +
    geom_line(data = summary_data,
              aes(x = Dilution_Condition_Num, y = Median, group = Group, color = "Median"),
              alpha = 1, linewidth = 2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1.5,
             fill = "white", alpha = 1) +
    geom_hline(yintercept = c(1, 1.5), linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = c("Median" = "black")) +
    scale_x_continuous(breaks = x_values, labels = names(x_values), trans = "reverse") +
    scale_y_continuous(
      limits = c(0.5, 6.5),
      breaks = c(0.5, 1, 1.5, 2, 3, 4, 5, 6),
      labels = c("0", "0.5", "1.5", "2", "3", "4", "5", "6")
    ) +
    labs(title = plot_title, x = "Dilution & Condition", y = "Log10 Mean Protein Abundance") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_line(size = 0.5),
      panel.grid.major = element_line(size = 1),
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
      plot.title = element_text(size = 40, hjust = 0.5, color = "black"),
      axis.title = element_text(size = 40, color = "black"),
      axis.text.x = element_text(size = 40, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 40, color = "black"),
      legend.position = "none"
    )
  
  # Ensure output folder exists and save
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  ggsave(p, path = save_path, filename = paste0(output_name, ".pdf"),
         width = 5262, height = 5262, units = "px")
  
  return(summary_list)
}

# ======== Paths (replace <FILL_IN_PATH> with your directories) ================
dir1        <- "<FILL_IN_PATH>/DIA_only1"
dir3        <- "<FILL_IN_PATH>/DIA_only3"
dir9        <- "<FILL_IN_PATH>/DIA_only9"
design_path <- "<FILL_IN_PATH>/design.csv"
save_path   <- "<FILL_IN_PATH>/output"

# Ensure output directory exists
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# ======== Discover base names from PYE1 folder ================================
files1     <- list.files(dir1, pattern = "\\.tsv$", full.names = FALSE)
base_names <- unique(gsub("_PYE[139]\\.tsv$", "", files1))

# ======== Run per base name and collect summary ===============================
summary_list <- list()

for (base in base_names) {
  file_paths <- c(
    file.path(dir1, paste0(base, "_PYE1.tsv")),
    file.path(dir3, paste0(base, "_PYE3.tsv")),
    file.path(dir9, paste0(base, "_PYE9.tsv"))
  )
  
  plot_title  <- paste(base, "(ECOLI)")
  output_name <- paste0("spaghetti_", base)
  
  summary_list <- generate_protein_dilution_plot(
    file_paths, design_path, save_path, plot_title, output_name, summary_list
  )
}

# ======== Save run-level summary ==============================================
summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
fwrite(summary_dt, file.path(save_path, "protein_filter_summary.csv"))

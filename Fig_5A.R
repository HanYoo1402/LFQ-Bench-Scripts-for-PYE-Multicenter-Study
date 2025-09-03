###############################################################################
# DIA vs DDA PYE_A vs PYE_B — per-feature log2 ratio (violin plot)
#
# Purpose
#   For each report file, compute per-feature (e.g., Protein.Group / Protein IDs)
#   condition means (PYE_A, PYE_B) after filtering and contaminant removal,
#   derive log2 ratio A/B, and visualize distributions per Lab.Setup as
#   violin plots (separately for Ecoli, Human, Yeast). Horizontal dashed lines
#   mark expected ratios:
#       Human = 0,  Yeast = −log2(3),  Ecoli = 1.
#   This script processes DIA and DDA in parallel and combines them in one figure.
#
# Inputs
#   DIA:
#     - data_path_dia   : folder with DIA reports (.tsv)
#     - design_path_dia : design table with at least {Run, Condition, TechnicalReplicate, Lab.Setup}
#     - feature_id_dia  : e.g., "Protein.Group"
#     - quantity_dia    : e.g., "PG.MaxLFQ"
#     - replicate_count_filter_dia : min non-NA replicates per condition (default 3)
#   DDA:
#     - data_path_dda   : folder with DDA reports (e.g., MaxQuant .txt exported or pre-processed)
#     - design_path_dda : design with at least {Sample, Condition, TechnicalReplicate, Lab.Setup}
#     - feature_id_dda  : e.g., "Protein IDs"
#     - quantity_dda    : e.g., "LFQ_Intensity"
#     - replicate_count_filter_dda : min non-NA replicates per condition (default 3)
#
# Output
#   - combined_plot_labeled.pdf  (three stacked panels: Ecoli, Human, Yeast;
#                                 DIA and DDA alternating per Lab.Setup)
#
# Reproducibility / Methods
#   • DIA filtering: Q.Value, Lib.Q.Value, PG.Q.Value, Lib.PG.Q.Value ≤ 0.01;
#                    remove contaminants via "_CONTA" in Protein.Names;
#                    organism assignment via suffixes: _HUMAN, _YEAS8, _ECOLI.
#   • DDA filtering: remove contaminants/decoys via patterns in
#                    `Fasta headers` and `Protein IDs` (CON|REV).
#   • Per-feature means are computed within Condition (PYE_A, PYE_B);
#     features retained if Count_PYE_A ≥ ReplicateCountFilter and
#     Count_PYE_B ≥ ReplicateCountFilter.
#   • Response variable: Ratio = log2(Mean_PYE_A / Mean_PYE_B).
#   • Visualization: violin with median crossbar; dashed reference at expected ratio.
#   • Color scheme: distinct palettes for DIA vs DDA; legend suppressed.
#
###############################################################################

# Load necessary libraries
library(tidyverse)
library(data.table)
library(patchwork)

# Define the function for processing DIA data
PYE_function_DIA <- function(ReportPath, DesignPath, FeatureId, Quantity, ReplicateCountFilter = 3) {
  
  report <- fread(ReportPath)
  design <- fread(DesignPath)
  
  report_filtered <- report %>%
    filter(Q.Value <= 0.01) %>%
    filter(Lib.Q.Value <= 0.01) %>%
    filter(PG.Q.Value <= 0.01) %>%
    filter(Lib.PG.Q.Value <= 0.01) %>%
    filter(!str_detect(Protein.Names, pattern = "_CONTA")) %>%
    left_join(y = design, by = "Run") %>%
    mutate(Human = str_detect(Protein.Names, pattern = "_HUMAN")) %>%
    mutate(Yeast = str_detect(Protein.Names, pattern = "_YEAS8")) %>%
    mutate(Ecoli = str_detect(Protein.Names, pattern = "_ECOLI")) %>%
    mutate(Organism = case_when(
      Human & !Yeast & !Ecoli ~ "Human",
      !Human & Yeast & !Ecoli ~ "Yeast",
      !Human & !Yeast & Ecoli ~ "Ecoli",
      TRUE ~ "NA"
    )) %>%
    filter(Organism != "NA")
  
  report_processed <- report_filtered %>%
    select(Run, Condition, TechnicalReplicate, "FeatureId" = FeatureId, "Quantity" = Quantity, Organism, Lab.Setup) %>%
    distinct() %>%
    group_by(FeatureId, Condition, Organism, Lab.Setup) %>%
    summarize(Mean = mean(Quantity, na.rm = TRUE), Count = sum(!is.na(Quantity))) %>%
    ungroup() %>%
    pivot_wider(names_from = Condition, values_from = c(Mean, Count)) %>%
    mutate(Ratio = log2(Mean_PYE_A / Mean_PYE_B)) %>%
    drop_na(Ratio) %>%
    filter(Count_PYE_A >= ReplicateCountFilter & Count_PYE_B >= ReplicateCountFilter)
  
  return(report_processed)
}

# Define the function for processing DDA data
PYE_function_DDA <- function(ReportPath, DesignPath, FeatureId, Quantity, ReplicateCountFilter = 3) {
  
  report <- fread(ReportPath)
  design <- fread(DesignPath)
  
  report_filtered <- report %>%
    filter(!str_detect(`Fasta headers`, pattern = "_CONTA")) %>%
    filter(!str_detect(`Protein IDs`, pattern = "CON|REV")) %>%
    left_join(y = design, by = "Sample") %>%
    mutate(Human = str_detect(`Fasta headers`, pattern = "_HUMAN")) %>%
    mutate(Yeast = str_detect(`Fasta headers`, pattern = "_YEAS8")) %>%
    mutate(Ecoli = str_detect(`Fasta headers`, pattern = "_ECOLI")) %>%
    mutate(Organism = if_else(Human == TRUE & Yeast == FALSE & Ecoli == FALSE, true = "Human",
                              false = if_else(Human == FALSE & Yeast == TRUE & Ecoli == FALSE, true = "Yeast",
                                              false = if_else(Human == FALSE & Yeast == FALSE & Ecoli == TRUE, true = "Ecoli", false = "NA")))) %>%
    filter(LFQ_Intensity >= 1) %>%
    filter(Organism == "Human" | Organism == "Yeast" | Organism == "Ecoli") %>%
    mutate(ExpectedRatio = if_else(Organism == "Human", true = 0,
                                   false = if_else(Organism == "Yeast", true = -log2(3),
                                                   false = if_else(Organism == "Ecoli", true = 1, false = NaN))))
  
  report_processed <- report_filtered %>%
    select(Sample, Condition, TechnicalReplicate, "FeatureId" = FeatureId, "Quantity" = Quantity, Organism, Lab.Setup) %>%
    distinct() %>%
    group_by(FeatureId, Condition, Organism, Lab.Setup) %>%
    summarize(Mean = mean(Quantity, na.rm = TRUE), Count = sum(!is.na(Quantity))) %>%
    ungroup() %>%
    pivot_wider(names_from = Condition, values_from = c(Mean, Count)) %>%
    mutate(Ratio = log2(Mean_PYE_A / Mean_PYE_B)) %>%
    drop_na(Ratio) %>%
    filter(Count_PYE_A >= ReplicateCountFilter & Count_PYE_B >= ReplicateCountFilter) %>%
    filter(!is.infinite(Ratio))  # Filter out rows where Ratio is Inf
  
  return(report_processed)
}

# Paths for DIA
data_path_dia <- "<FILL_IN_PATH>/DIA"                
design_path_dia <- "<FILL_IN_PATH>/design.csv"   
feature_id_dia <- "Protein.Group"
quantity_dia <- "PG.MaxLFQ"
replicate_count_filter_dia <- 3

# Paths for DDA
data_path_dda <- "<FILL_IN_PATH>/DDA"        
design_path_dda <- "<FILL_IN_PATH>/design.csv"
feature_id_dda <- "Protein IDs"
quantity_dda <- "LFQ_Intensity"
replicate_count_filter_dda <- 3

# Collect data from all files in DIA folder
files_dia <- list.files(path = data_path_dia, full.names = TRUE, recursive = FALSE)
all_data_dia <- lapply(files_dia, PYE_function_DIA, DesignPath = design_path_dia, FeatureId = feature_id_dia, Quantity = quantity_dia, ReplicateCountFilter = replicate_count_filter_dia)

# Combine all DIA data into one data frame
combined_data_dia <- bind_rows(all_data_dia)

# Collect data from all files in DDA folder
files_dda <- list.files(path = data_path_dda, full.names = TRUE, recursive = FALSE)
all_data_dda <- lapply(files_dda, PYE_function_DDA, DesignPath = design_path_dda, FeatureId = feature_id_dda, Quantity = quantity_dda, ReplicateCountFilter = replicate_count_filter_dda)

# Combine all DDA data into one data frame
combined_data_dda <- bind_rows(all_data_dda)

# Combine both DIA and DDA data into one data frame
combined_data_dia$Lab.Setup <- paste("DIA_", combined_data_dia$Lab.Setup, sep = "")
combined_data_dda$Lab.Setup <- paste("DDA_", combined_data_dda$Lab.Setup, sep = "")
combined_data <- bind_rows(combined_data_dia, combined_data_dda)

# Custom color scale
custom_colors_DDA <- c("Human" = "#8d3e01", "Yeast" = "#4f4a8c", "Ecoli" = "#105d46")
custom_colors_DIA <- c("Human" = "#fd862a", "Yeast" = "#a5a2ce", "Ecoli" = "#2bdba6")
# List of organisms
organisms <- c("Ecoli", "Human", "Yeast")

# Expected ratios for organisms
expected_ratios <- c("Human" = 0, "Yeast" = -log2(3), "Ecoli" = 1)

# Desired order of Lab.Setups (adjusted with DIA and DDA setups)
desired_order <- c("DDA_D_ulti_ecl", "DIA_D_ulti_ecl", "DDA_H_ulti_ecl", "DIA_H_ulti_ecl", "DDA_A_evo_ex", "DIA_A_evo_ex",
                   "DDA_D_Vanq_ex", "DIA_D_Vanq_ex", "DDA_E_Vanq_ex", "DIA_E_Vanq_ex", "DDA_A_ulti_ex", "DIA_A_ulti_ex",
                   "DDA_G_ulti_ex", "DIA_G_ulti_ex", "DDA_H_ulti_ex", "DIA_H_ulti_ex", "DDA_B_ulti_HF", "DIA_B_ulti_HF",
                   "DDA_B_ulti_HFX", "DIA_B_ulti_HFX", "DDA_E_ulti_lumos", "DIA_E_ulti_lumos",
                   "DDA_C_nE_tTOF", "DIA_C_nE_tTOF", "DDA_G_nLC_tTOF",  "DIA_G_nLC_tTOF") 

# Update Lab.Setup to be a factor with the desired order
combined_data$Lab.Setup <- factor(combined_data$Lab.Setup, levels = desired_order)

# Initialize a list to store individual plots
plots <- list()

# Loop through each organism and create separate plots
for (org in organisms) {
  org_data <- combined_data %>%
    filter(Organism == org)
  
  expected_ratio <- expected_ratios[org]
  y_limits <- c(expected_ratio - 1, expected_ratio + 1)
  
  # Determine the color scale to use based on DIA or DDA
  org_data <- org_data %>%
    mutate(Color = ifelse(str_detect(Lab.Setup, "DIA_"), custom_colors_DIA[Organism], custom_colors_DDA[Organism]))
  
  violin_plot <- org_data %>%
    ggplot(aes(x = Lab.Setup, y = Ratio, fill = Color)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
    geom_hline(yintercept = expected_ratio, linetype = "dashed", color = "red") +
    labs(x = NULL, y = NULL, color = NULL) +
    scale_fill_identity() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size = 40),
          axis.title = element_text(size = 40),
          plot.title = element_text(size = 40, hjust = 0.5),
          axis.text.x = element_text(size = 40, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 40)) +
    coord_cartesian(ylim = y_limits) +
    ggtitle(org) +
    guides(color = "none", fill = "none") +
    geom_vline(xintercept = seq(2.5, length(desired_order), by = 2),
               linetype = "dashed", color = "black")
  
  plots[[org]] <- violin_plot
}

# Combine individual plots into one grid
combined_plot <- wrap_plots(plots[[1]], plots[[2]], plots[[3]], ncol = 1)

# Save the combined plot
output_path <- "<FILL_IN_PATH>/output"
combined_output_name <- "combined_violin_plot.pdf"
ggsave(filename = file.path(output_path, combined_output_name),
       plot = combined_plot, device = "pdf", bg = "white",
       width = 9920, height = 7016, units = "px")

###############################################################################
# Violin plots of log2(PYE_A / PYE_B) ratios per Lab.Setup and organism (DIA)
#
# Purpose
#   For each DIA report in `data_path`, compute per-feature (e.g., Protein.Group)
#   condition means (PYE_A, PYE_B) after 1% FDR filtering and contaminant removal.
#   Derive per-feature log2 ratios (Mean_PYE_A / Mean_PYE_B), require a minimum
#   number of non-NA replicates per condition, and visualize distributions as
#   violins for Ecoli, Human, and Yeast across Lab.Setups. Dashed horizontal
#   lines indicate expected ratios:
#       Human = 0,  Yeast = −log2(3),  Ecoli = 1.
#
# Inputs
#   - data_path     : folder with DIA reports (.tsv), each containing columns used below
#   - design_path   : design table with at least {Run, Condition, TechnicalReplicate, Lab.Setup}
#   - feature_id    : feature identifier column name (e.g., "Protein.Group")
#   - quantity      : intensity column name (e.g., "PG.MaxLFQ")
#   - replicate_count_filter : minimum non-NA replicates per condition (default 3)
#   - output_path   : folder to save the combined PDF
#
# Output
#   - combined_<feature_id>.pdf  : 3 stacked violin plots (Ecoli, Human, Yeast)
#
# Reproducibility notes
#   - 1% FDR thresholds: Q.Value, Lib.Q.Value, PG.Q.Value, Lib.PG.Q.Value ≤ 0.01
#   - Contaminants filtered via "_CONTA" in Protein.Names
#   - Organism assignment from Protein.Names suffix: _HUMAN, _YEAS8, _ECOLI
#   - Per-feature ratios computed from condition means (no random ops)
###############################################################################

library(tidyverse)
library(data.table)
library(stringr)
library(patchwork)

# Define the function for processing data (DIA)
PYE_function <- function(ReportPath, DesignPath, FeatureId, Quantity, ReplicateCountFilter = 3) {
  
  report <- fread(ReportPath)
  design <- fread(DesignPath)
  
  # --- Filtering + metadata merge + organism annotation ---
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
  
  # --- Per-feature condition means + ratio; replicate count gating ------------
  report_processed <- report_filtered %>%
    select(Run, Condition, TechnicalReplicate,
           "FeatureId" = FeatureId, "Quantity" = Quantity,
           Organism, Lab.Setup) %>%
    distinct() %>%
    group_by(FeatureId, Condition, Organism, Lab.Setup) %>%
    summarize(Mean = mean(Quantity, na.rm = TRUE),
              Count = sum(!is.na(Quantity)),
              .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = c(Mean, Count)) %>%
    mutate(Ratio = log2(Mean_PYE_A / Mean_PYE_B)) %>%
    drop_na(Ratio) %>%
    filter(Count_PYE_A >= ReplicateCountFilter & Count_PYE_B >= ReplicateCountFilter)
  
  return(report_processed)
}

# ------------------------------- Paths & params -------------------------------
data_path              <- "<FILL_IN_PATH>/DIA" 
design_path            <- "<FILL_IN_PATH>/design.csv"
output_path            <- "<FILL_IN_PATH>/output"     
feature_id             <- "Protein.Group"              
quantity               <- "PG.MaxLFQ"                      
replicate_count_filter <- 3                                   

# ------------------------------- Read & process -------------------------------
files <- list.files(path = data_path, full.names = TRUE, recursive = FALSE)
all_data <- lapply(files, PYE_function,
                   DesignPath = design_path,
                   FeatureId = feature_id,
                   Quantity = quantity,
                   ReplicateCountFilter = replicate_count_filter)

# Combine all reports
combined_data <- bind_rows(all_data)

# ------------------------------- Plot settings --------------------------------
custom_colors <- c("Human" = "#D95F02", "Yeast" = "#7570B3", "Ecoli" = "#1B9E77")
organisms     <- c("Ecoli", "Human", "Yeast")
expected_ratios <- c("Human" = 0, "Yeast" = -log2(3), "Ecoli" = 1)

# Lab.Setup display order (adjustable)
desired_order <- c("D_ulti_ecl", "H_ulti_ecl", "A_evo_ex", "I_nLC_ex", "D_Vanq_ex",
                   "E_Vanq_ex", "A_ulti_ex", "G_ulti_ex", "H_ulti_ex", "L_ulti_ex",
                   "L_ulti_ex_FAIMS", "B_ulti_HF", "B_ulti_HFX", "E_ulti_lumos",
                   "L_nAcqu_tTOF", "C_nE_tTOF", "L_nE_tTOF", "G_nLC_tTOF",
                   "J_nLC_tTOF", "K_Mclass_zTOF")
combined_data$Lab.Setup <- factor(combined_data$Lab.Setup, levels = desired_order)

# (Optional) 95th percentile trimming per organism
filtered_data <- combined_data %>%
  group_by(Organism) %>%
  mutate(Percentile_95 = quantile(Ratio, 0.95, na.rm = TRUE)) %>%
  filter(Ratio <= Percentile_95) %>%
  ungroup()

# ------------------------------- Build plots ----------------------------------
plots <- list()

for (org in organisms) {
  org_data <- combined_data %>% filter(Organism == org)
  expected_ratio <- expected_ratios[org]
  y_limits <- c(expected_ratio - 1, expected_ratio + 1)
  
  violin_plot <- org_data %>%
    ggplot(aes(x = Lab.Setup, y = Ratio, fill = Organism)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
    geom_hline(yintercept = expected_ratio, linetype = "dashed", color = "red") +
    labs(x = NULL, y = NULL, color = NULL) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 40),
      axis.title = element_text(size = 40),
      plot.title = element_text(size = 40, hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 40)
    ) +
    coord_cartesian(ylim = y_limits) +
    ggtitle(org) +
    guides(color = "none", fill = "none")
  
  plots[[org]] <- violin_plot
}

# Stack Ecoli, Human, Yeast vertically (1 column)
combined_plot <- wrap_plots(plots[[1]], plots[[2]], plots[[3]], ncol = 1)

# ------------------------------- Save figure ----------------------------------
combined_output_name <- paste0("combined_", feature_id, ".pdf")
ggsave(
  filename = file.path(output_path, combined_output_name),
  plot = combined_plot, device = "pdf", bg = "white",
  width = 7440, height = 7016, units = "px"
)

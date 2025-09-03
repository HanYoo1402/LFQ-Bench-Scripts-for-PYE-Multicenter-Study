###############################################################################
# PYE_A vs PYE_B log2 ratio — lower tertile violin plots (Ecoli / Human / Yeast)
#
# Description:
#   For each input report, the script filters identifications (1% FDR on
#   multiple levels), removes contaminants, merges design, assigns Organism
#   from protein names, and summarizes per FeatureId (e.g., Protein.Group):
#     • Mean intensity per Condition (PYE_A, PYE_B)
#     • Replicate counts per Condition
#     • log2 ratio: log2(Mean_PYE_A / Mean_PYE_B)
#   It then keeps only the lower tertile of ratios per Organism and
#   generates 1×3 violin plots (ECOLI, HUMAN, YEAST) with median bars and
#   expected reference lines (Human = 0, Yeast = −log2(3), Ecoli = 1).
#
# Inputs (edit these):
#   • data_path    : folder with DIA *.tsv reports
#   • design_path  : CSV with at least {Run, Condition, TechnicalReplicate, Lab.Setup}
#   • output_path  : output folder for PDF (created if missing)
#   • feature_id   : column name for the biological feature (e.g., "Protein.Group")
#   • quantity     : column name for intensity (e.g., "PG.MaxLFQ")
#
# Output:
#   • combined_<FeatureId>_lower_tertile.pdf
#
# Requirements:
#   • R ≥ 4.3.2
#   • Packages: tidyverse, data.table, stringr, patchwork
#
# Notes:
#   • ReplicateCountFilter enforces minimum non-NA replicates per condition.
#   • Organism is inferred via suffixes in Protein.Names: _HUMAN / _YEAS8 / _ECOLI.
###############################################################################

library(tidyverse)
library(data.table)
library(stringr)
library(patchwork)

# -------- Function ------------------------------------------------------------
PYE_function <- function(ReportPath, DesignPath, FeatureId, Quantity, ReplicateCountFilter = 3) {
  report <- fread(ReportPath)
  design <- fread(DesignPath)
  
  report_filtered <- report %>%
    filter(
      Q.Value        <= 0.01,
      Lib.Q.Value    <= 0.01,
      PG.Q.Value     <= 0.01,
      Lib.PG.Q.Value <= 0.01,
      !str_detect(Protein.Names, "_CONTA")
    ) %>%
    left_join(design, by = "Run") %>%
    mutate(
      Human = str_detect(Protein.Names, "_HUMAN"),
      Yeast = str_detect(Protein.Names, "_YEAS8"),
      Ecoli = str_detect(Protein.Names, "_ECOLI"),
      Organism = case_when(
        Human & !Yeast & !Ecoli ~ "Human",
        !Human & Yeast & !Ecoli ~ "Yeast",
        !Human & !Yeast & Ecoli ~ "Ecoli",
        TRUE ~ "NA"
      )
    ) %>%
    filter(Organism != "NA")
  
  # Robust column capture using tidy-eval so FeatureId/Quantity can be strings
  report_processed <- report_filtered %>%
    select(
      Run, Condition, TechnicalReplicate,
      FeatureId = all_of(FeatureId),
      Quantity  = all_of(Quantity),
      Organism, Lab.Setup
    ) %>%
    distinct() %>%
    group_by(FeatureId, Condition, Organism, Lab.Setup) %>%
    summarize(
      Mean  = mean(Quantity, na.rm = TRUE),
      Count = sum(!is.na(Quantity)),
      .groups = "drop_last"
    ) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = Condition, values_from = c(Mean, Count)) %>%
    mutate(Ratio = log2(Mean_PYE_A / Mean_PYE_B)) %>%
    drop_na(Ratio) %>%
    filter(Count_PYE_A >= ReplicateCountFilter, Count_PYE_B >= ReplicateCountFilter)
  
  return(report_processed)
}

# -------- Paths / Params (edit these) -----------------------------------------
data_path  <- "<FILL_IN_PATH>/reports"   # folder with *.tsv reports
design_path <- "<FILL_IN_PATH>/design.csv"
output_path <- "<FILL_IN_PATH>/output"

feature_id <- "Protein.Group"
quantity   <- "PG.MaxLFQ"
replicate_count_filter <- 3

# -------- Read all reports and process ----------------------------------------
files <- list.files(path = data_path, pattern = "\\.tsv$", full.names = TRUE, recursive = FALSE)
all_data <- lapply(
  files,
  PYE_function,
  DesignPath = design_path,
  FeatureId = feature_id,
  Quantity  = quantity,
  ReplicateCountFilter = replicate_count_filter
)
combined_data <- bind_rows(all_data)

# -------- Styling / Expected values -------------------------------------------
custom_colors <- c("Human" = "#D95F02", "Yeast" = "#7570B3", "Ecoli" = "#1B9E77")
organisms     <- c("Ecoli", "Human", "Yeast")
expected_ratios <- c("Human" = 0, "Yeast" = -log2(3), "Ecoli" = 1)

desired_order <- c(
  "D_ulti_ecl","H_ulti_ecl","A_evo_ex","I_nLC_ex","D_Vanq_ex","E_Vanq_ex",
  "A_ulti_ex","G_ulti_ex","H_ulti_ex","L_ulti_ex","L_ulti_ex_FAIMS","B_ulti_HF",
  "B_ulti_HFX","E_ulti_lumos","L_nAcqu_tTOF","C_nE_tTOF","L_nE_tTOF","G_nLC_tTOF",
  "J_nLC_tTOF","K_Mclass_zTOF"
)
combined_data$Lab.Setup <- factor(combined_data$Lab.Setup, levels = desired_order)

# Keep only the lower tertile per organism
filtered_data <- combined_data %>%
  group_by(Organism) %>%
  mutate(tertile_group = ntile(Ratio, 3)) %>%
  filter(tertile_group == 1) %>%
  ungroup()

# -------- Build one violin per organism ---------------------------------------
plots <- list()
for (org in organisms) {
  org_data <- filtered_data %>% filter(Organism == org)
  expected_ratio <- expected_ratios[[org]]
  y_limits <- c(expected_ratio - 1, expected_ratio + 1)
  
  plots[[org]] <-
    ggplot(org_data, aes(x = Lab.Setup, y = Ratio, fill = Organism)) +
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
      axis.text.y = element_text(size = 40),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = y_limits) +
    ggtitle(org)
}

# -------- Combine in a 1×3 layout (fixed order) -------------------------------
combined_plot <- wrap_plots(plots[organisms], ncol = 1)

# -------- Save ----------------------------------------------------------------
out_pdf <- paste0("combined_", feature_id, "_lower_tertile.pdf")
ggsave(filename = file.path(output_path, out_pdf),
       plot = combined_plot, device = "pdf", bg = "white",
       width = 7440, height = 7016, units = "px")

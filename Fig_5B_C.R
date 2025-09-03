###############################################################################
# PYE_A vs PYE_B log2 ratio bench plot (per report)
#
# Purpose
#   For each DIA report (.tsv), compute per-feature (e.g., Protein.Group)
#   condition means (PYE_A, PYE_B) after standard 1% FDR filtering and
#   contaminant removal, annotate organism from Protein.Names, normalize
#   intensity to the report-wise maximum, and plot:
#       y = log2(Mean_PYE_A / Mean_PYE_B)  vs  x = log2(Mean_PYE_A)
#   Dashed horizontal lines indicate expected ratios:
#       Human = 0,  Yeast = −log2(3),  Ecoli = 1.
#   Same script structure can be applied to DDA reports —
#   just replace input files and adjust the column headers
#   to match the DDA report format.
#
# Inputs
#   - ReportPath  : DIA report (.tsv)
#   - DesignPath  : design table with {Run, Condition, TechnicalReplicate, Lab.Setup}
#   - FeatureId   : feature identifier column (e.g., "Protein.Group")
#   - Quantity    : intensity column (e.g., "PG.MaxLFQ")
#   - ReplicateCountFilter : minimum non-NA replicates per condition (default 3)
#   - OutputPath  : folder to save the PDF
#
# Output
#   - <report_basename>_<FeatureId>_max.pdf (scatter plot)
#
# Reproducibility notes
#   - 1% FDR thresholds: Q.Value, Lib.Q.Value, PG.Q.Value, Lib.PG.Q.Value ≤ 0.01
#   - Contaminants removed via "_CONTA" pattern in Protein.Names
#   - Organism assignment from Protein.Names suffix: _HUMAN, _YEAS8, _ECOLI
#   - Features retained if Count_PYE_A ≥ ReplicateCountFilter and
#     Count_PYE_B ≥ ReplicateCountFilter
#   - Plot uses no random operations; results are deterministic given inputs
###############################################################################

library(tidyverse)
library(ggplot2)
library(data.table)

PYE_function <- function(ReportPath, DesignPath, FeatureId, Quantity, ReplicateCountFilter = 3, OutputPath) {
  
  # Read report file
  report <- fread(ReportPath)
  
  # Read design file
  design <- fread(DesignPath)
  
  # Merge report with design
  merged_data <- report %>%
    filter(Q.Value <= 0.01,
           Lib.Q.Value <= 0.01,
           PG.Q.Value <= 0.01,
           Lib.PG.Q.Value <= 0.01,
           !str_detect(Protein.Names, pattern = "_CONTA")) %>%
    left_join(design, by = "Run") %>%
    mutate(Human = str_detect(Protein.Names, pattern = "_HUMAN"),
           Yeast = str_detect(Protein.Names, pattern = "_YEAS8"),
           Ecoli = str_detect(Protein.Names, pattern = "_ECOLI"),
           Organism = case_when(
             Human & !Yeast & !Ecoli ~ "Human",
             !Human & Yeast & !Ecoli ~ "Yeast",
             !Human & !Yeast & Ecoli ~ "Ecoli",
             TRUE ~ "NA"
           ),
           ExpectedRatio = case_when(
             Organism == "Human" ~ 0,
             Organism == "Yeast" ~ -log2(3),
             Organism == "Ecoli" ~ 1,
             TRUE ~ NaN
           ))
  
  # Normalize intensity to the maximum intensity in the report (paper-ready note)
  merged_data <- merged_data %>%
    mutate(MedianIntensity = max(get(Quantity), na.rm = TRUE),
           NormalizedQuantity = get(Quantity) / MedianIntensity) %>%
    ungroup()
  
  # Process data
  report_processed <- merged_data %>%
    select(Run, Condition, TechnicalReplicate, "FeatureId" = FeatureId, "Quantity" = NormalizedQuantity, Organism, ExpectedRatio, Lab.Setup) %>%
    distinct() %>%
    group_by(FeatureId, Condition, Organism, ExpectedRatio) %>%
    summarize(Mean = mean(Quantity, na.rm = TRUE), Count = sum(!is.na(Quantity))) %>%
    ungroup() %>%
    pivot_wider(names_from = Condition, values_from = c(Mean, Count)) %>%
    mutate(Ratio = log2(Mean_PYE_A / Mean_PYE_B)) %>%
    drop_na(Ratio) %>%
    filter(Count_PYE_A >= ReplicateCountFilter & Count_PYE_B >= ReplicateCountFilter)
  
  # Create plot
  point_plot <- report_processed %>%
    ggplot(aes(x = log2(Mean_PYE_A), y = Ratio, color = Organism)) +  
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log2(3), linetype = "dashed") +
    geom_point(size = 2, alpha = 0.5, shape = 19) +
    theme_classic(base_size = 10, base_family = "sans") +
    scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
    ylim(-3, 3) +
    labs(title = unique(merged_data$Lab.Setup)) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 4)) + 
    theme(
      axis.text = element_text(size = 40),
      axis.title = element_blank(),
      plot.title = element_text(size = 40, hjust = 0.5)
    ) +
    guides(color = "none")
  
  # Save the plot
  OutputName <- str_c(str_c(str_remove(basename(ReportPath), pattern = ".tsv"), FeatureId, sep = "_"), sep = "_max")
  ggsave(filename = str_c(OutputName, ".pdf"), path = OutputPath, plot = point_plot, device = "pdf", width = 222, height = 105, units = "mm")
}

# Apply function to all report files
lapply(X = list.files(path = "<FILL_IN_PATH>/DIA_labs",
                      full.names = TRUE, recursive = FALSE),
       FUN = PYE_function,
       DesignPath = "<FILL_IN_PATH>/design.csv",
       FeatureId = "Protein.Group",
       Quantity = "PG.MaxLFQ",
       ReplicateCountFilter = 3,
       OutputPath = "<FILL_IN_PATH>/results")

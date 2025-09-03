###############################################################################
# CV boxplots across Lab.Setups (Condition A and B) — DIA, DDA PYE1
#
# Description:
#   Reads DIA .tsv reports, merges with design, applies FDR filters, removes
#   contaminants, computes per-protein CV (%) within each Lab.Setup for
#   Condition PYE_A and PYE_B, and plots CV distributions with custom whiskers:
#     ymin=5th, lower=25th, middle=50th, upper=75th, ymax=90th percentile.
#   Same script structure can be applied to DDA reports —
#   just replace input files and adjust the column headers
#   to match the DDA report format.
#
# Inputs (edit “Paths”):
#   • report_path : folder with DIA *.tsv reports or DDA *.txt reports (adjust in script)
#   • design_file : CSV with at least {Run, Lab.Setup, Condition}
#   • save_path   : output folder for PDFs (created if missing)
#
# Outputs:
#   • CV_DIA_PYE1A.pdf
#   • CV_DIA_PYE1B.pdf
#
# Requirements:
#   • R ≥ 4.3.2
#   • Packages: data.table, dplyr, stringr, ggplot2
#
# Notes:
#   • CV is computed per Protein.Group within each Lab.Setup:
#       CV% = 100 * sd(PG.MaxLFQ) / mean(PG.MaxLFQ)
#   • The y-axis guidelines (10% and 25%) are visual references only.
###############################################################################

# Load required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# -------- Paths (replace <FILL_IN_PATH> with your directories) ----------------
report_path <- "<FILL_IN_PATH>/DIA_PYE1"
design_file <- "<FILL_IN_PATH>/design.csv"
save_path   <- "<FILL_IN_PATH>/output"

norm <- "PG.MaxLFQ"  # column used for CV calculation (DDA: "LFQ_Intensity")

# -------- Read reports & design ----------------------------------------------
files  <- list.files(report_path, pattern = "\\.tsv$", full.names = TRUE)
design <- data.table::fread(design_file)

# Read each file and stack
all_data <- lapply(files, data.table::fread)
combined_data <- data.table::rbindlist(all_data, use.names = TRUE, fill = TRUE)

# Merge with design and filter
combined_data <- merge(combined_data, design, by = "Run") %>%
  filter(
    Q.Value        <= 0.01,
    PG.Q.Value     <= 0.01,
    Lib.Q.Value    <= 0.01,
    Lib.PG.Q.Value <= 0.01,
    PG.MaxLFQ > 0,
    !str_detect(Protein.Names, "_CONTA")
  ) %>%
  distinct(Protein.Group, Run, PG.MaxLFQ, Lab.Setup, Condition) %>%
  setDT()

# Subset by Condition
data_A <- combined_data[grepl("PYE_A", Condition)]
data_B <- combined_data[grepl("PYE_B", Condition)]

# -------- CV per Protein.Group × Lab.Setup -----------------------------------
data_A[, cvplatetype := 100 * stats::sd(get(norm), na.rm = TRUE) /
         mean(get(norm), na.rm = TRUE),
       by = .(Protein.Group, Lab.Setup)]
data_B[, cvplatetype := 100 * stats::sd(get(norm), na.rm = TRUE) /
         mean(get(norm), na.rm = TRUE),
       by = .(Protein.Group, Lab.Setup)]

# Lab.Setup display order
desired_order <- c(
  "D_ulti_ecl","H_ulti_ecl","A_evo_ex","I_nLC_ex","D_Vanq_ex","E_Vanq_ex",
  "A_ulti_ex","G_ulti_ex","H_ulti_ex","L_ulti_ex","L_ulti_ex_FAIMS","B_ulti_HF",
  "B_ulti_HFX","E_ulti_lumos","L_nAcqu_tTOF","C_nE_tTOF","L_nE_tTOF","G_nLC_tTOF",
  "J_nLC_tTOF","K_Mclass_zTOF"
)
data_A[, Lab.Setup := factor(Lab.Setup, levels = desired_order)]
data_B[, Lab.Setup := factor(Lab.Setup, levels = desired_order)]

# -------- Custom boxplot stats (5–90% whiskers) -------------------------------
custom_boxplot_stats <- function(x) {
  qs <- stats::quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  names(qs) <- c("ymin","lower","middle","upper","ymax")
  qs
}

# Precompute stats for segments
stats_A <- data_A[, as.list(custom_boxplot_stats(cvplatetype)), by = Lab.Setup]
stats_B <- data_B[, as.list(custom_boxplot_stats(cvplatetype)), by = Lab.Setup]

# -------- Plot A (PYE_A) ------------------------------------------------------
plot_A <- ggplot(data_A, aes(x = Lab.Setup, y = cvplatetype)) +
  # Draw box using custom 5/25/50/75/90 percentiles
  geom_boxplot(
    aes(ymin = ..lower.., lower = ..lower.., middle = ..middle..,
        upper = ..upper.., ymax = ..ymax..),
    stat = "summary", fun.data = custom_boxplot_stats, fill = "darksalmon"
  ) +
  # Explicit whisker/hinge lines from precomputed stats
  geom_segment(data = stats_A, aes(x = Lab.Setup, xend = Lab.Setup, y = ymin,  yend = ymax)) +
  geom_segment(data = stats_A, aes(x = Lab.Setup, xend = Lab.Setup, y = lower, yend = upper)) +
  geom_segment(data = stats_A, aes(x = as.numeric(Lab.Setup) - 0.25,
                                   xend = as.numeric(Lab.Setup) + 0.25,
                                   y = ymin, yend = ymin)) +
  geom_segment(data = stats_A, aes(x = as.numeric(Lab.Setup) - 0.25,
                                   xend = as.numeric(Lab.Setup) + 0.25,
                                   y = ymax, yend = ymax)) +
  labs(x = NULL, y = "Coefficient of Variation (%)", title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40),
    axis.title  = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title= element_text(size = 40),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_blank()
  ) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "darkorange1") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#669900") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 115))

# -------- Plot B (PYE_B) ------------------------------------------------------
plot_B <- ggplot(data_B, aes(x = Lab.Setup, y = cvplatetype)) +
  geom_boxplot(
    aes(ymin = ..lower.., lower = ..lower.., middle = ..middle..,
        upper = ..upper.., ymax = ..ymax..),
    stat = "summary", fun.data = custom_boxplot_stats, fill = "darksalmon"
  ) +
  geom_segment(data = stats_B, aes(x = Lab.Setup, xend = Lab.Setup, y = ymin,  yend = ymax)) +
  geom_segment(data = stats_B, aes(x = Lab.Setup, xend = Lab.Setup, y = lower, yend = upper)) +
  geom_segment(data = stats_B, aes(x = as.numeric(Lab.Setup) - 0.25,
                                   xend = as.numeric(Lab.Setup) + 0.25,
                                   y = ymin, yend = ymin)) +
  geom_segment(data = stats_B, aes(x = as.numeric(Lab.Setup) - 0.25,
                                   xend = as.numeric(Lab.Setup) + 0.25,
                                   y = ymax, yend = ymax)) +
  labs(x = NULL, y = "Coefficient of Variation (%)", title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40),
    axis.title  = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title= element_text(size = 40),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_blank()
  ) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "darkorange1") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#669900") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 115))

# -------- Save ----------------------------------------------------------------
ggsave(filename = file.path(save_path, "CV_DIA_PYE1A.pdf"),
       plot = plot_A, width = 7440, height = 5262, units = "px")

ggsave(filename = file.path(save_path, "CV_DIA_PYE1B.pdf"),
       plot = plot_B, width = 7440, height = 5262, units = "px")

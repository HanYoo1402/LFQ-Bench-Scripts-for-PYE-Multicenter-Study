###############################################################################
# Retention-time CV boxplots across Lab.Setups
#
# Description:
#   Reads DIA .tsv reports, merges with the design table, computes per-precursor
#   RT CV (%) within each Lab.Setup, and visualizes distributions using custom
#   boxplots with 5th/25th/50th/75th/90th percentiles (whiskers at 5% and 90%).
#   Same script structure can be applied to DDA reports —
#   just replace input files and adjust the column headers
#   to match the DDA report format.
#
# Inputs (edit “Paths”):
#   • report_path : folder with DIA *.tsv reports or DDA *.txt reports (adjust in script)
#   • design_file : CSV with at least {Run, Lab.Setup}
#   • save_path   : output folder for PDF (created if missing)
#
# Output:
#   • RT_CV.pdf
#   • RT_CV_STATS.csv (export of 5/25/50/75/90 percentiles per Lab.Setup)
#
# Requirements:
#   • R ≥ 4.3.2
#   • Packages: data.table, ggplot2
#
# Notes:
#   • CV is computed per Precursor.Id × Lab.Setup as:
#       RT_CV% = 100 * sd(RT) / mean(RT)
#   • Within-setup downsampling to ≤ 20,000 rows is applied for plotting speed
#     only; it does not affect the statistics exported to CSV.
###############################################################################

# Load required libraries
library(data.table)
library(ggplot2)

# -------- Paths (replace <FILL_IN_PATH> with your directories) ----------------
report_path <- "<FILL_IN_PATH>/reports"
design_file <- "<FILL_IN_PATH>/design.csv"
save_path   <- "<FILL_IN_PATH>/output"

# Column to use for grouping CV
norm <- "RT"  # retention time column

# -------- Read reports & design ----------------------------------------------
files  <- list.files(report_path, pattern = "\\.tsv$", full.names = TRUE)
design <- fread(design_file)

# Read each file and stack
all_data <- lapply(files, fread)
combined_data <- rbindlist(all_data, use.names = TRUE, fill = TRUE)

# Merge with design
combined_data <- merge(combined_data, design, by = "Run")

# -------- Compute RT CV per Precursor.Id × Lab.Setup --------------------------
combined_data[, RT_CV := sd(get(norm), na.rm = TRUE) / mean(get(norm), na.rm = TRUE) * 100,
              by = .(Precursor.Id, Lab.Setup)]

# Display order for Lab.Setup
desired_order <- c(
  "D_ulti_ecl","H_ulti_ecl","A_evo_ex","I_nLC_ex","D_Vanq_ex","E_Vanq_ex",
  "A_ulti_ex","G_ulti_ex","H_ulti_ex","L_ulti_ex","L_ulti_ex_FAIMS","B_ulti_HF",
  "B_ulti_HFX","E_ulti_lumos","L_nAcqu_tTOF","C_nE_tTOF","L_nE_tTOF","G_nLC_tTOF",
  "J_nLC_tTOF","K_Mclass_zTOF"
)
combined_data[, Lab.Setup := factor(Lab.Setup, levels = desired_order)]

# Optional: per-setup downsampling for plotting speed (keeps stats independent)
plot_data <- combined_data[, .SD[sample(.N, min(.N, 20000))], by = Lab.Setup]

# -------- Custom boxplot stats (5–90% whiskers) -------------------------------
custom_boxplot_stats <- function(x) {
  qs <- quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  names(qs) <- c("ymin","lower","middle","upper","ymax")
  qs
}

# Precompute stats for export and whisker segments (use full data, not sampled)
stats_RT <- combined_data[, as.list(custom_boxplot_stats(RT_CV)), by = Lab.Setup]
fwrite(stats_RT, file.path(save_path, "RT_CV_STATS_DIA.csv"))

# -------- Plot ----------------------------------------------------------------
p1 <- ggplot(plot_data, aes(x = Lab.Setup, y = RT_CV)) +
  # Draw box using custom 5/25/50/75/90 percentiles
  geom_boxplot(
    aes(ymin = ..lower.., lower = ..lower.., middle = ..middle..,
        upper = ..upper.., ymax = ..ymax..),
    stat = "summary", fun.data = custom_boxplot_stats, fill = "azure3"
  ) +
  # Explicit whisker/hinge lines from precomputed stats
  geom_segment(data = stats_RT, aes(x = Lab.Setup, xend = Lab.Setup, y = ymin,  yend = ymax)) +
  geom_segment(data = stats_RT, aes(x = Lab.Setup, xend = Lab.Setup, y = lower, yend = upper)) +
  geom_segment(data = stats_RT, aes(x = as.numeric(Lab.Setup) - 0.25,
                                    xend = as.numeric(Lab.Setup) + 0.25,
                                    y = ymin, yend = ymin)) +
  geom_segment(data = stats_RT, aes(x = as.numeric(Lab.Setup) - 0.25,
                                    xend = as.numeric(Lab.Setup) + 0.25,
                                    y = ymax, yend = ymax)) +
  labs(x = "Lab.Setup", y = "Retention time CV [%]", title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 40),
    axis.title  = element_text(size = 40),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  coord_cartesian(ylim = c(0, 3))  # same as ylim(0,3) but safer for stats

# -------- Save ----------------------------------------------------------------
ggsave(filename = file.path(save_path, "RT_CV_DIA_PYE1.pdf"),
       plot = p1, units = "px", width = 4960, height = 3508)

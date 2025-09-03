###############################################################################
# Scatter plots: PYE1 vs PYE3 (log10 normalized intensities)
#
# Description:
# Reads reports from specified PYE directories (e.g., PYE1, PYE3, PYE9),
# applies FDR filters, averages PG.MaxLFQ per protein within each PYE/Lab.Setup,
# normalizes within PYE, and reshapes data to compare PYE pairs
# (PYE1 vs PYE3, PYE1 vs PYE9, PYE3 vs PYE9). Produces two 1:1 scatter plots:
#   • HUMAN vs ECOLI proteins
#   • HUMAN vs YEAST proteins
# Each plot includes R² annotations from simple linear models.
#
# Upstream filtering (performed before this script):
#  Protein groups identified by ≥2 peptides were retained, e.g.:
#     group_by(Protein.Group) %>%
#       mutate(npep = n_distinct(Stripped.Sequence)) %>%
#       ungroup() %>%
#       filter(npep >= 2)
#
# Usage:
#   1. Update the folder_path, save_path, and design_file with your directories.
#   2. Run the script in R. A PDF figure will be saved in the output directory.
#
# Requirements:
#   - R (4.3.2) (used version for generating plots)
#   - Packages: data.table, ggplot2, dplyr, stringr, tidyr
###############################################################################

# Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

# -------- I/O helper ----------------------------------------------------------
# Read and row-bind all .tsv (DIA) files in a directory (.txt for DDA)
read_and_combine <- function(directory_path) {
  files <- list.files(path = directory_path, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(files)) stop("No .tsv files found in: ", directory_path)
  data_list <- lapply(files, fread)
  bind_rows(data_list)
}

# -------- Paths (example) -----------------------------------------------------
path1       <- "<FILL_IN_PATH>/PYE_Plots/PYE1"   # PYE1 folder
path3       <- "<FILL_IN_PATH>/PYE_Plots/PYE3"   # PYE3 folder
design      <- fread("<FILL_IN_PATH>/design.csv")  # design table
save_path   <- "<FILL_IN_PATH>/output"                # output folder for PDFs

# -------- Read data -----------------------------------------------------------
pye1_data <- read_and_combine(path1)
pye3_data <- read_and_combine(path3)

# -------- Processing pipeline -------------------------------------------------
process_data <- function(data, design, quantity) {
  merged_data <- data %>%
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
             Human & !Yeast & !Ecoli ~ "HUMAN",
             !Human & Yeast & !Ecoli ~ "YEAST",
             !Human & !Yeast & Ecoli ~ "ECOLI",
             TRUE ~ "Other"
           )) %>%
    select(Run, Protein.Names, PG.MaxLFQ, Organism, Lab.Setup, PYE, Dilution) %>%
    group_by(PYE, Lab.Setup, Protein.Names) %>%
    summarize(PG.MaxLFQ = mean(PG.MaxLFQ, na.rm = TRUE)) %>%
    mutate(NormalizedQuantity = log10(PG.MaxLFQ / max(PG.MaxLFQ, na.rm = TRUE))) %>%
    ungroup() %>%
    drop_na()
  return(merged_data)
}

# Process both datasets
pye1_processed <- process_data(pye1_data, design, "NormalizedQuantity")
setDT(pye1_processed)
pye1_processed <- dcast(pye1_processed, Protein.Names ~ Lab.Setup,
                        value.var = "NormalizedQuantity", fun.aggregate = mean)

pye3_processed <- process_data(pye3_data, design, "NormalizedQuantity")
setDT(pye3_processed)
pye3_processed <- dcast(pye3_processed, Protein.Names ~ Lab.Setup,
                        value.var = "NormalizedQuantity", fun.aggregate = mean)

# Merge and reshape to long format
combined_data <- inner_join(pye1_processed, pye3_processed, by = "Protein.Names",
                            suffix = c(".PYE1", ".PYE3"))

combined_data_long <- combined_data %>%
  select(Protein.Names, matches("\\.PYE1$"), matches("\\.PYE3$")) %>%
  pivot_longer(cols = -Protein.Names, names_to = "Condition", values_to = "NormalizedQuantity") %>%
  separate(Condition, into = c("Lab.Setup", "PYE"), sep = "\\.") %>%
  mutate(Organism = case_when(
    grepl("HUMAN", Protein.Names) ~ "HUMAN",
    grepl("ECOLI", Protein.Names) ~ "ECOLI",
    grepl("YEAS8", Protein.Names) ~ "YEAST",
    TRUE ~ "Other"
  )) %>%
  pivot_wider(names_from = PYE, values_from = NormalizedQuantity)

# -------- Styling -------------------------------------------------------------
custom_colors <- c("HUMAN" = "#D95F02", "YEAST" = "#7570B3",
                   "ECOLI" = "#1B9E77", "Other" = "grey50")

combined_data_long <- combined_data_long %>%
  mutate(Organism = factor(Organism, levels = c("HUMAN", "ECOLI", "YEAST", "Other"))) %>%
  arrange(Organism)

custom_breaks <- c(0, -1, -2, -3, -4, -5)
custom_labels <- c("100", "10", "1", "0.1", "0.01", "0.001")

# -------- Conservative rounding to avoid overstating fit ----------------------
round_down_if_5 <- function(x, digits = 2) {
  scaled <- x * 10^(digits + 1)          # shift to inspect 3rd decimal
  third_decimal <- floor(scaled) %% 10   # isolate the 3rd decimal digit
  
  if (third_decimal == 5) {
    # truncate (floor) instead of rounding
    out <- floor(x * 10^digits) / 10^digits
  } else {
    # normal rounding
    out <- round(x, digits)
  }
  return(out)
}

# === Plot 1 datasets & R^2 ====================================================
data_human_ecoli <- combined_data_long %>% filter(Organism %in% c("HUMAN", "ECOLI"))
data_ecoli <- data_human_ecoli %>% filter(Organism == "ECOLI")
data_human <- data_human_ecoli %>% filter(Organism == "HUMAN")

lm_model_ecoli <- lm(PYE1 ~ PYE3, data = data_ecoli)
r2_ecoli <- summary(lm_model_ecoli)$r.squared
r2_label_ecoli <- paste0("R² (ECOLI) = ",
                         format(round_down_if_5(r2_ecoli, 2), nsmall = 2))

lm_model_human <- lm(PYE1 ~ PYE3, data = data_human)
r2_human <- summary(lm_model_human)$r.squared
r2_label_human <- paste0("R² (HUMAN) = ", round(r2_human, 2))

# === Plot 2 datasets & R^2 ====================================================
data_human_yeast <- combined_data_long %>% filter(Organism %in% c("HUMAN", "YEAST"))
data_yeast <- data_human_yeast %>% filter(Organism == "YEAST")
data_human2 <- data_human_yeast %>% filter(Organism == "HUMAN")

lm_model_yeast <- lm(PYE1 ~ PYE3, data = data_yeast)
r2_yeast <- summary(lm_model_yeast)$r.squared
r2_label_yeast <- paste0("R² (YEAST) = ",
                         format(round_down_if_5(r2_yeast, 2), nsmall = 2))

lm_model_human2 <- lm(PYE1 ~ PYE3, data = data_human2)
r2_human2 <- summary(lm_model_human2)$r.squared
r2_label_human2 <- paste0("R² (HUMAN) = ", round(r2_human2, 2))

# -------- Plot: HUMAN & ECOLI (adjust labels as needed)------------------------
p_ecoli <- ggplot(data_human_ecoli, aes(x = PYE3, y = PYE1, color = Organism)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0,        linetype = "dashed",
              color = custom_colors["HUMAN"], linewidth = 1) +
  geom_abline(slope = 1, intercept = log10(3), linetype = "dashed", # log10(9) if you compare PYE1 vs PYE9
              color = custom_colors["ECOLI"], linewidth = 1) +
  annotate("text", x = -5.8, y = -0.6, label = r2_label_ecoli, size = 12, hjust = 0) +
  annotate("text", x = -5.8, y = -1.2, label = r2_label_human, size = 12, hjust = 0) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = custom_breaks, labels = custom_labels,
                     limits = c(-6, 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = custom_breaks, labels = custom_labels,
                     limits = c(-6, 0.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Log10 Normalized Intensity [%] (PYE3)",
       y = "Log10 Normalized Intensity [%] (PYE1)") +
  theme_bw() +
  theme(
    text = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    panel.grid = element_blank()
  )

# -------- Plot: HUMAN & YEAST (adjust labels as needed)------------------------
p_yeast <- ggplot(data_human_yeast, aes(x = PYE3, y = PYE1, color = Organism)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0,        linetype = "dashed",
              color = custom_colors["HUMAN"], linewidth = 1) +
  geom_abline(slope = 1, intercept = log10(3), linetype = "dashed", # log10(9) if you compare PYE1 vs PYE9
              color = custom_colors["YEAST"], linewidth = 1) +
  annotate("text", x = -5.8, y = -0.6, label = r2_label_yeast, size = 12, hjust = 0) +
  annotate("text", x = -5.8, y = -1.2, label = r2_label_human2, size = 12, hjust = 0) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = custom_breaks, labels = custom_labels,
                     limits = c(-6, 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = custom_breaks, labels = custom_labels,
                     limits = c(-6, 0.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Log10 Normalized Intensity [%] (PYE3)",
       y = "Log10 Normalized Intensity [%] (PYE1)") +
  theme_bw() +
  theme(
    text = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    panel.grid = element_blank()
  )

# -------- Save (adjust width/height as needed)---------------------------------
ggsave(p_ecoli,
       file = file.path(save_path, "ECOLI_output.pdf"),
       units = "px", width = 7440, height = 5262)

ggsave(p_yeast,
       file = file.path(save_path, "YEAST_output.pdf"),
       units = "px", width = 7440, height = 5262)

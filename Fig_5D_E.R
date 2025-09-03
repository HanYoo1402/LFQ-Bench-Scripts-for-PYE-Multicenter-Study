###############################################################################
# Missingness vs. abundance rank with sigmoid fit (PYE_B reference; PYE_A test)
#
# Purpose
#   From multiple DIA reports (.tsv), compute protein abundance (normalized
#   PG.MaxLFQ) in PYE_B for YEAST proteins, rank proteins by abundance, and
#   quantify their “missingness” in PYE_A (fraction of PYE_A runs where a
#   protein is absent). Fit a logistic (sigmoid) model to relate rank to
#   missingness and visualize observed points with the fitted curve.
#
# Key steps
#   1) Read all *.tsv reports; merge with design by Run.
#   2) Split by Condition (PYE_B, PYE_A).
#   3) PYE_B (reference abundance):
#        - Keep YEAST (“_YEAS8”) proteins.
#        - Normalize PG.MaxLFQ within Run by the run’s maximum.
#        - For each Protein.Group, compute mean of normalized intensities.
#        - Rank proteins by descending mean (ties = first).
#   4) PYE_A (test presence):
#        - For each Protein.Group, compute missingness = 1 − (#runs present / total PYE_A runs).
#   5) Fit nls logistic model:  missingness ~ 1 / (1 + exp(-(a + b * rank))).
#   6) Plot scatter (rank vs missingness) and overlay fitted sigmoid curve.
#
# Inputs (edit paths below)
#   - folder_path : directory with DIA *.tsv reports
#   - design_path : design table with at least {Run, Condition, Lab.Setup}
#   - save_path   : output directory for PDF
#
# Output
#   - Missingness_DIA_PYE1_sigmoid.pdf
#
# Reproducibility notes
#   - Rank is computed on PYE_B mean normalized PG.MaxLFQ per Protein.Group.
#   - Missingness in PYE_A uses total_runs_A = uniqueN(PYE_A runs).
#   - Logistic start values: a = 0, b = 0.01 (nls); change only if convergence fails.
#   - Y-axis is bounded to [0, 1]; no random operations used.
###############################################################################

# Load required libraries
library(data.table)
library(ggplot2)
library(scales)

# Define paths (replace <FILL_IN_PATH> with your directories)

# ⚠️ When using DDA instead of DIA:
#   • Point folder_path to the correct DDA report files
#   • Adjust column headers (e.g., intensity column name) as required
# ⚠️ When comparing DDA vs DIA:
#   • Load both sets separately and ensure consistent column selection/renaming
folder_path <- "<FILL_IN_PATH>/DIA_reports"
design_path <- "<FILL_IN_PATH>/design.csv"
save_path   <- "<FILL_IN_PATH>/output"

# List all .tsv files in the folder
files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

all_data <- list()

# Read each file and append data to the list
for (file in files) {
  data <- fread(file)
  all_data[[file]] <- data
}

# Combine all data into one data.table
combined_data <- rbindlist(all_data, fill = TRUE)
design_df <- fread(design_path)

# Merge the combined data with the design dataframe by "Run"
merged_df <- merge(combined_data, design_df, by = "Run")

# Filter for PYE_B and PYE_A conditions
PYE_B_prots <- merged_df[Condition == "PYE_B"]
PYE_A_prots <- merged_df[Condition == "PYE_A"]

# Function to normalize PG.MaxLFQ and calculate mean normalized abundance
process_condition <- function(data) {
  data <- data[grepl("_YEAS8", Protein.Names)]
  data[, normalized_PG_MaxLFQ := PG.MaxLFQ / max(PG.MaxLFQ, na.rm = TRUE), by = Run]
  data[, Mean_MaxLFQ := mean(normalized_PG_MaxLFQ, na.rm = TRUE), by = .(Protein.Group)]
  unique_protein_groups <- unique(data[, .(Mean_MaxLFQ, Protein.Group)])
  return(unique_protein_groups)
}

# Process PYE_B
unique_protein_groups_B <- process_condition(PYE_B_prots)

# Calculate ranks based on PYE_B
unique_protein_groups_B[, rank := frank(-Mean_MaxLFQ, ties.method = "first")]

# Function to calculate missingness of PYE_A proteins relative to PYE_B
calculate_missingness <- function(protein_group, condition_A, condition_B, total_runs_A) {
  runs_present_A <- uniqueN(condition_A[Protein.Group == protein_group, Run])
  missingness <- 1 - (runs_present_A / total_runs_A)
  return(missingness)
}

# Total runs for PYE_A
total_runs_A <- uniqueN(PYE_A_prots$Run)

# Add missingness to unique_protein_groups_B
unique_protein_groups_B[, missingness_A := calculate_missingness(Protein.Group, PYE_A_prots, PYE_B_prots, total_runs_A), by = Protein.Group]

# Fit a sigmoid model (logistic function)
sigmoid_model <- nls(
  missingness_A ~ 1 / (1 + exp(-(a + b * rank))),
  data = unique_protein_groups_B,
  start = list(a = 0, b = 0.01)  
)

# Add fitted values to the dataset
unique_protein_groups_B[, fitted_missingness := predict(sigmoid_model, newdata = .SD)]

# Create the plot
p_sigmoid <- ggplot(unique_protein_groups_B, aes(x = rank, y = missingness_A)) +
  geom_point(shape = 1, color = "black") +
  geom_line(aes(y = fitted_missingness), color = "red", linewidth = 1.5) +
  labs(x = "Rank", y = "Missingness in PYE_A") +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 40),
    axis.text  = element_text(size = 40)
  ) +
  coord_cartesian(ylim = c(0, 1))

# Save the plot
ggsave(filename = "Missingness_DIA_PYE1_sigmoid.pdf",
       plot = p_sigmoid, path = save_path,
       width = 7440, height = 5262, units = "px")

###############################################################################################
# Protein Group Completeness & Gradient Length Plot
#
# Description:
# This script processes .tsv report files, calculates the completeness of
# protein groups across runs (i.e., in how many runs each protein group was
# identified), merges the results with experimental design metadata, and
# generates a bar plot showing completeness categories together with the
# gradient-length–normalized counts.
#
# Usage:
#   1. Update the folder_path, save_path, and design_file with your directories.
#   2. Run the script in R. A PDF figure will be saved in the output directory.
#
# Requirements:
#   - R (4.3.2) (used version for generating plots)
#   - Packages: data.table, ggplot2, dplyr, stringr
#
###############################################################################################

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# ---------------------------------------------------------------------------------------------
# Define user-specific paths (replace <FILL_IN_PATH> with actual locations)
# ---------------------------------------------------------------------------------------------
folder_path <- "<FILL_IN_PATH>/data"       # folder containing .tsv files
save_path   <- "<FILL_IN_PATH>/output"          # where to save the plots
design_file <- "<FILL_IN_PATH>/design.csv"      # design file (CSV)

# ---------------------------------------------------------------------------------------------
# Initialize storage list
# ---------------------------------------------------------------------------------------------
result_list <- list()

# ---------------------------------------------------------------------------------------------
# Loop through each result file (for DDA data use pattern = "\\.txt" and correct column header)
# ---------------------------------------------------------------------------------------------
for (file_name in list.files(folder_path, pattern = "\\.tsv", full.names = TRUE)) {
  
  # Read the file
  df <- fread(file_name)
  
  # Extract sample and protein identifiers
  sample_names   <- sort(unique(df$Run))
  unique_proteins <- unique(df$Protein.Group)
  
  # Create presence/absence matrix
  protein_presence <- matrix(
    0, nrow = length(unique_proteins), ncol = length(sample_names),
    dimnames = list(unique_proteins, sample_names)
  )
  
  for (i in 1:length(unique_proteins)) {
    protein_group <- unique_proteins[i]
    samples_with_protein <- df$Run[df$Protein.Group == protein_group]
    protein_presence[protein_group, samples_with_protein] <- 1
  }
  
  protein_presence_df <- as.data.frame(protein_presence)
  protein_presence_df$Appeared_In <- rowSums(protein_presence)
  
  # Count how many proteins appear in X samples
  count_appearances <- table(protein_presence_df$Appeared_In)
  max_samples <- max(as.numeric(names(count_appearances)))
  
  counts_unique <- data.frame(
    No.Samples = 1:max_samples,
    Count = ifelse(
      1:max_samples %in% names(count_appearances),
      count_appearances[as.character(1:max_samples)], 0
    )
  )
  
  # Function to assign completeness class
  get_completeness <- function(count) {
    if (count == max_samples) {
      return("complete")
    } else if (count >= max_samples / 2) {
      return("shared at least 50%")
    } else if (count == 1) {
      return("unique")
    } else {
      return("sparse")
    }
  }
  
  counts_unique$Completeness <- sapply(counts_unique$No.Samples, get_completeness)
  counts_unique$Origin <- file_name
  
  result_list[[file_name]] <- counts_unique
}

# ---------------------------------------------------------------------------------------------
# Combine results across all files
# ---------------------------------------------------------------------------------------------
final_df <- do.call(rbind, result_list)
sum_counts <- aggregate(Count ~ Completeness + Origin, data = final_df, sum)
sum_counts$Origin <- basename(sum_counts$Origin)

# ---------------------------------------------------------------------------------------------
# Merge with design metadata
# ---------------------------------------------------------------------------------------------
design_df <- fread(design_file)
merged_df <- merge(sum_counts, design_df, by.x = "Origin", by.y = "Origin", all.x = TRUE)

# Order factors
merged_df$Completeness <- factor(
  merged_df$Completeness,
  levels = c("unique", "sparse", "shared at least 50%", "complete")
)
merged_df$Lab.Setup <- factor(
  merged_df$Lab.Setup,
  levels = rev(unique(design_df$Lab.Setup))
)

# Custom display order (edit as needed)
desired_order <- rev(c(
  "D_ulti_ecl", "H_ulti_ecl", "A_evo_ex", "I_nLC_ex", "D_Vanq_ex",
  "E_Vanq_ex", "A_ulti_ex", "G_ulti_ex", "H_ulti_ex", "L_ulti_ex",
  "L_ulti_ex_FAIMS", "B_ulti_HF", "B_ulti_HFX", "E_ulti_lumos", "L_nAcqu_tTOF",
  "C_nE_tTOF", "L_nE_tTOF", "G_nLC_tTOF", "J_nLC_tTOF", "K_Mclass_zTOF"
))

# Function to format labels
format_labels <- function(labels) {
  sapply(labels, function(label) {
    parts <- strsplit(as.character(label), "_")[[1]]
    if (length(parts) >= 2) {
      regular_part <- str_c(parts[1:(length(parts) - 1)], collapse = "_")
      bold_part <- parts[length(parts)]
      bquote(.(regular_part) * "_" * bold(.(bold_part)))
    } else {
      label
    }
  })
}
formatted_labels <- format_labels(levels(merged_df$Lab.Setup))

# ---------------------------------------------------------------------------------------------
# Compute normalized counts per gradient length
# ---------------------------------------------------------------------------------------------
sum_counts_per_lab <- aggregate(Count ~ Lab.Setup, data = merged_df, sum)
merged_df <- merge(merged_df, sum_counts_per_lab, by = "Lab.Setup", suffixes = c("", ".sum"))
merged_df$Count_per_GL <- merged_df$Count.sum / merged_df$Gradient.Length * 25

# ---------------------------------------------------------------------------------------------
# Plot: stacked bar (Completeness) + normalized points/line
# ---------------------------------------------------------------------------------------------
p1 <- ggplot(merged_df, aes(x = Count, y = Lab.Setup, fill = Completeness)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_point(aes(x = Count_per_GL, y = Lab.Setup), color = "red", show.legend = FALSE, size = 4) +
  geom_line(aes(x = Count_per_GL, y = Lab.Setup), color = "red", group = 1, linewidth = 2) +
  geom_text(
    data = subset(merged_df, Completeness == "complete"),
    aes(x = Count - 750, label = Count),
    size = 40 * 0.352777778, color = "white", vjust = 0.5
  ) +
  labs(x = "No. identified Protein Groups", y = "") +
  scale_fill_manual(
    values = c(
      "unique" = "red", "sparse" = "orange",
      "shared at least 50%" = "darkgrey", "complete" = "black"
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    legend.position = "none",
    axis.title.x.top = element_blank(),
    axis.title.x.bottom = element_blank()
  ) +
  ggtitle("No. Protein Group Hits & per Gradient Length") +
  scale_x_continuous(
    limits = c(0, 5000), breaks = seq(0, 5000, by = 1000),
    labels = c("0", "1k", "2k", "3k", "4k", "5k"),
    sec.axis = sec_axis(~ . / 25, name = "PG Hits/GL")
  ) +
  scale_y_discrete(labels = formatted_labels)

# ---------------------------------------------------------------------------------------------
# Save to PDF (adjust width/height as needed)
# ---------------------------------------------------------------------------------------------
ggsave(
  p1, path = save_path, filename = "plot.pdf",
  units = "px", width = 2480, height = 5262
)

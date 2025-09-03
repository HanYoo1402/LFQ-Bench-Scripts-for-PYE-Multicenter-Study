###############################################################################
# Bar plot: identified protein groups per Lab.Setup and Dilution (SPECIES)
#
# Description:
#   Reads report files from a folder, merges with a design table, removes
#   contaminants/decoys, counts unique Protein.Group per (PYE, Lab.Setup, Dilution),
#   and plots counts by Lab.Setup with Dilution as grouped bars. Counts for PYE1
#   are annotated on the bars. Same script structure can be applied to DDA reports —
#   just replace input files and adjust the column headers
#   to match the DDA report format.
#
# Requirements:
#   • R >= 4.3.2
#   • Packages: data.table, ggplot2
#
# Notes:
#   • Adjust the color palette and x-axis limits to match each SPECIES.
#   • The legend is hidden by design. The secondary x-axis is a fixed duplicate
#     (for aesthetics; to mirror the axis styling of Figure 2A).
###############################################################################

# Load libraries
library(data.table)
library(ggplot2)

# -------- Paths (replace <FILL_IN_PATH> with your directories) ----------------
folder_path <- "<FILL_IN_PATH>/xxx_Species"            # input .tsv folder
design_path <- "<FILL_IN_PATH>/design"                 # design file
save_path   <- "<FILL_IN_PATH>/output"                 # output folder

# Ensure output directory exists
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# -------- Read input tables ---------------------------------------------------
# Read all .tsv files in the folder
files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
dfs <- lapply(files, fread)
combined_df <- rbindlist(dfs, use.names = TRUE, fill = TRUE)

# Remove contaminants/decoys
combined_df <- combined_df[!grepl("CON|REV", Protein.Group)]

# Read design and merge
design <- fread(design_path)
dt <- merge(combined_df, design, by = "Run")

# -------- Count unique protein groups ----------------------------------------
n2 <- dt[, .(Count = uniqueN(Protein.Group)), by = .(PYE, Lab.Setup, Dilution)]

# Lab.Setup display order (edit as needed)
desired_order <- rev(c(
  "D_ulti_ecl", "H_ulti_ecl", "A_evo_ex", "I_nLC_ex", "D_Vanq_ex",
  "E_Vanq_ex", "A_ulti_ex", "G_ulti_ex", "H_ulti_ex", "L_ulti_ex",
  "L_ulti_ex_FAIMS", "B_ulti_HF", "B_ulti_HFX", "E_ulti_lumos", "L_nAcqu_tTOF",
  "C_nE_tTOF", "L_nE_tTOF", "G_nLC_tTOF", "J_nLC_tTOF", "K_Mclass_zTOF"
))
n2[, Lab.Setup := factor(Lab.Setup, levels = desired_order)]

# Manual colors for Dilution (choose one palette per species)
# HUMAN: colors <- c("#fd862a", "#d95f02", "#8d3e01")
# ECOLI: colors <- c("#2bdba6", "#1b9e77", "#105d46")
# YEAST: colors <- c("#a5a2ce", "#7570b3", "#4f4a8c")
colors <- c("#fd862a", "#d95f02", "#8d3e01")  # <-- set for SPECIES

# Subset for label annotations (PYE1 only)
n2_pye1 <- n2[Dilution == "PYE1"]

# -------- Plot ----------------------------------------------------------------
p1 <- ggplot(n2, aes(
  x = Count, y = Lab.Setup,
  fill = factor(Dilution, levels = rev(levels(factor(Dilution))))
)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(
    data = n2_pye1, aes(label = Count),
    hjust = -0.2, vjust = 0.2, size = 40 * 0.352777778
  ) +
  scale_fill_manual(values = colors) +
  labs(x = "No. of identified proteins", y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 40),
    axis.text.y = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  # Adjust limits/breaks per SPECIES if necessary
  scale_x_continuous(
    limits = c(0, 600),
    breaks = seq(0, 600, by = 300),
    labels = c("0", "300", "600"),
    sec.axis = sec_axis(~ ., name = "PG Hits/GL")
  )

# -------- Save ----------------------------------------------------------------
ggsave(
  filename = file.path(save_path, "DIA_SPECIES_UNLABELED_counts.pdf"),
  plot = p1, units = "px", width = 2480, height = 5262
)

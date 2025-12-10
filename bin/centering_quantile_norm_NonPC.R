#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(plotly)
library(ggplot2)
library(ggpointdensity)
library(ggpmisc)
library(ggridges)
library(gridExtra)
library(viridis)
library(preprocessCore)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)

#### I. Read in data ####
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 7) {
  stop("Usage: Rscript insufficient arguments.", call. = FALSE)
}

input_file <- args[1]
plot_file <- args[2]
ori_cr_bg <- args[3]
#corr_cr_bg <- args[4]
qnorm_cr_bg <- args[4]
qnorm_cr_file <- args[5]
drop_out_plot <- args[6]
species <- args[7]
#drop_out_filter <- args[8]

con.cr.dt <- fread(input_file)
colnames(con.cr.dt) <- c("chrom", "start", "end", "con_nb", "uncon_nb", "cr", "CG")
setDT(con.cr.dt)

# Get C-centered contexts based on reference genome of specified species
if (species == "mm10") {
	genome <- BSgenome.Mmusculus.UCSC.mm10
} else if (species == "hg38") {
	genome <- BSgenome.Hsapiens.UCSC.hg38
}
positions <- GRanges(seqnames = con.cr.dt$chrom, ranges = IRanges(start = con.cr.dt$start + 1, width = 1), strand = "+")
extended_positions <- resize(positions, width = 9, fix = "center")  # 9bp in total（+/- 4bp）
con.cr.dt$seq <- as.character(getSeq(genome, extended_positions))
con.cr.dt$fillin <- 1

# Convert to reverse complementary sequences for those centered around G bases
con.cr.dt[, seq_rc := fifelse(
  substr(seq, 5, 5) == "G",
  stringi::stri_reverse(chartr("ATCG", "TAGC", seq)),
  seq
)]
con.cr.dt[, seq_21 := substr(seq_rc, 3, 6)] # Extract the -2/+1 bp (can be specialized for certain needs)

# Filter out rows where the -2/+1 context contains "N" or lower cases
con.cr.dt <- con.cr.dt[!grepl("[a-zNn]", seq_21)]

#### II. Quality control and normalization ####
start.time <- Sys.time()

# Calculate M-values (log odds) with pseudo-count = 0.5
con.cr.dt[, m_value := log2((con_nb + 0.5) / (uncon_nb + 0.5))]

# Filter out contexts with high proportion of dropouts
DROPOUT_THRESHOLD <- 0.15
context_stats <- con.cr.dt[, .(dropout_prop = sum(con_nb == 0) / .N), by = .(seq_21)]
contexts_to_keep <- context_stats[dropout_prop <= DROPOUT_THRESHOLD, seq_21]
dt_filtered <- con.cr.dt[seq_21 %in% contexts_to_keep]

# Calculate the original mode for each distribution and select reference context
original_mode_estimates <- dt_filtered[con_nb > 0, {
  dens <- density(m_value, bw = "nrd0")
  .(original_mode = dens$x[which.max(dens$y)])
}, by = .(seq_21)]

reference_context_info <- original_mode_estimates[which.min(abs(original_mode))]
reference_context_name <- reference_context_info$seq_21

# Centralization and quantile normalization
dt_filtered <- merge(dt_filtered, original_mode_estimates, by = "seq_21")
dt_filtered[, m_value_centered := m_value - original_mode]

target_distribution <- sort(dt_filtered[seq_21 == reference_context_name, m_value_centered])
dt_filtered[, m_value_final := {
  p_rank <- frankv(m_value_centered, ties.method = "average") / .N
  quantile(target_distribution, probs = p_rank, type = 7, names = FALSE)
}, by = .(seq_21)]

dt_filtered[, cr_final := (2^m_value_final) / (1 + 2^m_value_final)]

# Save the results
write.table(con.cr.dt %>% select(chrom, start, end, cr), ori_cr_bg, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
#write.table(con.cr.dt %>% select(chrom, start, end, cr_corr), corr_cr_bg, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dt_filtered %>% select(chrom, start, end, cr_final), qnorm_cr_bg, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dt_filtered %>% select(chrom, start, end, con_nb, uncon_nb, seq, cr, fillin, cr_final), qnorm_cr_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#### III. Visualizaiton ####
theme_set(theme_minimal())
con.cr.subset.dt <- con.cr.dt %>% sample_n(min(1000000, nrow(con.cr.dt)))

p1 <- ggplot(con.cr.subset.dt, aes(x = cr, y = seq_21, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "C.R.", option = "C") +
    labs(title = "Conversion ratio across different contexts", 
         x = "Conversion ratio", y = "Contexts") +
    coord_cartesian(xlim = c(0,1))

p2 <- ggplot(con.cr.subset.dt, aes(x = m_value, y = seq_21, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "M-value", option = "C") +
    labs(title = "M-values across different contexts", 
         x = "M-values", y = "Contexts")

#p3 <- ggplot(con.cr.subset.dt, aes(x = cr_corr, y = seq_21, fill = stat(x))) +
#    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#    scale_fill_viridis_c(name = "C.R.", option = "C") +
#    labs(title = "Control-norm CR across different contexts", 
#         x = "Conversion ratio", y = "Contexts") +
#    coord_cartesian(xlim = c(0,1))

p4 <- ggplot(dt_filtered, aes(x = cr_final, y = seq_21, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "C.R.", option = "C") +
    labs(title = "Quantile-norm CR across different contexts", 
         x = "Conversion ratio", y = "Contexts") +
    coord_cartesian(xlim = c(0,1))

# Save the figures 
combined_plot <- grid.arrange(p1, p2, p4, nrow = 1)
ggsave(plot_file, combined_plot, width = 18, height = 10)
cat("Visualization saved to:", plot_file, "\n")

# Draw dropout figures
p2_dropout_rates <- ggplot(context_stats, aes(x = reorder(seq_21, dropout_prop), y = dropout_prop)) +
  geom_col(aes(fill = seq_21 %in% contexts_to_keep)) +
  geom_hline(yintercept = DROPOUT_THRESHOLD, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "firebrick", "TRUE" = "steelblue"),
                    labels = c("Filtered Out", "Kept")) +
  labs(
    title = "Quality Control by Dropout Rate",
    subtitle = paste("Contexts with dropout rates above", DROPOUT_THRESHOLD * 100, "% are removed"),
    x = "Sequence Context (-2/+1)",
    y = "Proportion of Dropout Sites (No converted reads)",
    fill = "Status"
  ) +
  coord_flip() + # Flip for better readability of context names
  theme_bw()

ggsave(drop_out_plot, plot = p2_dropout_rates, height = 8, width = 12)

# Print the run time
end.time <- Sys.time()
cat("Total processing time:", format(end.time - start.time), "\n")
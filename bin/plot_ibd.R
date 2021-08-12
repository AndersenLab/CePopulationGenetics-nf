library(tidyverse)
load("haplotype_plot_df.Rda")

#========#
#  Plot  #
#========#

#==============================#
# Plot ~ Swept haplotypes only #
#==============================#

strain_labels <- plot_df %>%
                    dplyr::ungroup() %>%
                    dplyr::select(plotpoint, isotype) %>%
                    dplyr::distinct() %>%
                    dplyr::arrange(plotpoint)

ggplot(plot_df,
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
           fill = swept_haplotype)) +
geom_rect() +
scale_fill_manual(values = c("Gray", "Red")) +
scale_y_continuous(breaks = strain_labels$plotpoint,
                   labels = strain_labels$isotype,
                   expand = c(0, 0)) +
xlab("Position (Mb)") +
theme_bw() +
facet_grid(.~chromosome, scales="free", space="free") +
theme(legend.position="none")

ggsave(paste("max_haplotype_genome_wide.pdf"),
       width = 14,
       height = 28)

#===============#
# Sweep summary #
#===============#

# Plot swept by chromosome & sorted by isotype

is_overlapping <- function(start_1, end_1, start_2, end_2) {
    if (
        (dplyr::between(start_1, start_2, end_2)) |
        (dplyr::between(end_1, start_2, end_2)) |
        (dplyr::between(start_2, start_1, end_1)) |
        (dplyr::between(end_2, start_1, end_1))
      ) {
      return(TRUE)
    }
  FALSE
}

sweep_summary <- plot_df %>%
  dplyr::select(chromosome, isotype, max_haplotype_shared, is_swept) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(is_swept = ifelse(is.na(is_swept), F, is_swept)) %>%
  dplyr::arrange(chromosome, isotype)

suffix <- function(x) {
  paste0(x,"_hapshare")
}

hap_share <- sweep_summary %>%
    dplyr::select(-is_swept) %>%
    tidyr::spread(chromosome, max_haplotype_shared) %>%
    dplyr::rename_at(.vars=vars(-isotype), funs(suffix))

sweep_summary %>%
  dplyr::select(-max_haplotype_shared) %>%
  tidyr::spread(chromosome, is_swept) %>%
  dplyr::mutate(II = F, III = F) %>%
  dplyr::left_join(hap_share) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(swept_chroms = sum(I, IV, V, X)) %>%
  readr::write_tsv("sweep_summary.tsv")

chrom_plots <- lapply(c("I", "II", "III", "IV", "V", "X"), function(x) {
ranked_by_sharing <- plot_df %>%
    dplyr::ungroup() %>%
    dplyr::filter(chromosome == x) %>%
    dplyr::group_by(isotype) %>%
    dplyr::select(isotype, max_haplotype_shared) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(max_haplotype_shared)) %>%
    dplyr::mutate(plotpoint=row_number()) %>%
    dplyr::select(-max_haplotype_shared) %>%
    dplyr::arrange(plotpoint)

strain_labels <- ranked_by_sharing$isotype

plot_df %>%
  dplyr::filter(chromosome == x) %>%
  dplyr::select(-plotpoint) %>%
  dplyr::left_join(ranked_by_sharing) %>%
  ggplot(.,
         aes(xmin = start/1E6, xmax = stop/1E6,
             ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
             fill = is_swept)) +
  geom_rect() +
  scale_fill_manual(values = c("TRUE"="Red", "FALSE"="Gray")) +
  scale_y_continuous(breaks = 1:length(strain_labels),
                     labels = strain_labels,
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none")
})

cowplot::plot_grid(plotlist = chrom_plots, ncol=6)

ggsave("max_haplotype_sorted_genome_wide.pdf",
       width = 32,
       height = 28)

#===================================#
# Distribution of sharing by strain #
#===================================#

ggplot(plot_df) +
  geom_histogram(aes(x = hap_length)) +
  scale_x_log10(labels = scales::comma, limits = c(1, 1E6)) +
  labs(x = "Haplotype Length", y = "Count")

ggsave("haplotype_length.pdf", height = 10, width = 10)


#=======================#
# Normal haplotype plot #
#=======================#

mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
mcolor <- mcolor_grp$color
names(mcolor) <- mcolor_grp$haplotype

strain_labels <- plot_df %>%
                    dplyr::select(isotype, plotpoint)

ggplot(plot_df,
       aes(xmin = start/1E6, xmax = stop/1E6,
           ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
           fill = haplotype)) +
  geom_rect() +
  scale_fill_manual(values = mcolor) +
  scale_y_continuous(breaks = strain_labels$plotpoint,
                     labels = strain_labels$isotype,
                     expand = c(0, 0)) +
  xlab("Position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme(legend.position="none")

ggsave("haplotype.pdf", height = 48, width = 24)

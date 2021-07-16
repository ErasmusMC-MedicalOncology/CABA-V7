# Author:                      Job van Riet
# Date:                        29-06-2021
# Function:                    Generate an overview of the CNV-analysis of the QIASeq panel.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(patchwork)

# Helper theme.
source('R/misc_themes.R')


# Load data ---------------------------------------------------------------

data.Pat1 <- readxl::read_xlsx('Misc./CNV_Regionlevel_Track-L6571_I20-1113-63_S63_L001_R1_001 (paired), L4929_....xlsx', sheet = 2, trim_ws = T, guess_max = 99999)
data.Pat2 <- readxl::read_xlsx('Misc./CNV_Genelevel_track-L6563_I20-1113-59_S59_L001_R1_001 (paired), L4921_I2....xlsx', sheet = 2, trim_ws = T, guess_max = 99999)
data.Pat3 <- readxl::read_xlsx('Misc./Mapped_UMI_reads_tumor-L6503I20-1101-10_S10_L001_R1_001 (Region CNVs)_20....xlsx', sheet = 3, trim_ws = T, guess_max = 99999)

data.Pat1$Sample <- 'L-6571 (No 8q amp.)'
data.Pat2$Sample <- 'L-6563 (No 8q amp.)'
data.Pat3$Sample <- 'L-6503 (With 8q. amp.'

data.CNV <- dplyr::bind_rows(data.Pat1, data.Pat2, data.Pat3)

# Clean-up.
data.CNV <- data.CNV %>%
    tidyr::separate(col = `Region (joined targets)`, into = c('start', 'end'), sep = '\\..') %>%
    dplyr::mutate(
        Chromosome = gsub('\\..*', '', Chromosome),
        start = as.numeric(start),
        end = as.numeric(end),
        region = sprintf('%s - %s', start, end),
        ) %>%
    dplyr::group_by(region, Sample, Chromosome) %>%
    dplyr::summarise(
        start = min(start),
        end = max(end),
        `Regional fold-change` = median(`Regional fold-change`),
        SYMBOL = ifelse(length(unique(unlist(strsplit(Name, ';')))) <= 4, paste0(unique(unlist(strsplit(Name, ';'))), collapse = ', '), '>4 genes'),
        SYMBOL = sprintf('<strong>%s</strong><br><sub>%s - %s</sub>', SYMBOL, start, end),
        consequence = ifelse('Loss' %in% `Regional consequence`, 'Loss', NA),
        consequence = ifelse('Gain' %in% `Regional consequence`, 'Gain', consequence),
        label = ifelse(any(`Regional p-value` <= 0.05), gtools::stars.pval(min(`Regional p-value`)), NA)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        `Regional fold-change` = ifelse(`Regional fold-change` < -5, -5.1, `Regional fold-change`),
        `Regional fold-change` = ifelse(`Regional fold-change` > 5, 5.1, `Regional fold-change`),
        Chromosome = factor(Chromosome, levels = gtools::mixedsort(unique(Chromosome)))
    ) %>%
    dplyr::distinct()

plotCNV <- function(data){
    ggplot2::ggplot(data, aes(x = SYMBOL, y = `Regional fold-change`,  fill = consequence, label = label)) +
        ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .7) +
        ggplot2::geom_hline(yintercept = 1.4, color = '#ED1F24', lty = 11) +
        ggplot2::geom_hline(yintercept = 0, color = 'black', lty = 'solid') +
        ggplot2::geom_hline(yintercept = -1.4, color = '#115E81', lty = 11) +
        ggplot2::scale_y_continuous(limits = c(-5.5, 5.5), breaks = c(-5, -2.5, -1.4, 0, 1.4, 2.5, 5), expand = c(0,0)) +
        ggplot2::geom_text(aes(y = `Regional fold-change` + .75 * sign(`Regional fold-change`)), size = 4, fontface = 'bold') +
        ggplot2::labs(x = 'Targeted Genes', y = 'Regional Fold Change') +
        ggplot2::scale_fill_manual(values = c('Gain' = '#ED1F24', 'Loss' = '#115E81'), na.value = 'white', guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
        ggplot2::facet_grid(Sample ~ Chromosome, scales = 'free', space = 'free_x') +
        theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = .5, size = 7))
}

plots <- lapply(split(data.CNV, data.CNV$Sample), plotCNV)

plots$`L-6503 (With 8q. amp.` + plots$`L-6563 (No 8q amp.)` + plots$`L-6571 (No 8q amp.)` +
    patchwork::plot_layout(nrow = 3, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))

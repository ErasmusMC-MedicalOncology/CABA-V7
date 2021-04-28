# Author:                      Job van Riet
# Date:                        03-03-2021
# Function:                    Generate supplementary figures on the CABA-V7 cohort.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Helper theme.
source('R/misc_themes.R')


# Import data -------------------------------------------------------------

overviewPatients <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 1)

# Add number of patients per group.
overviewPatients <- overviewPatients %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::mutate(
        `AR-V7 (Baseline)` = ifelse(`AR-V7 (Baseline)` == 'Pos.', 'AR-V7<sup>Pos.</sup>', ifelse(`AR-V7 (Baseline)` == 'Neg.', 'AR-V7<sup>Neg.</sup>', 'AR-V7<sup>Und.</sup>')),
        `AR-V7 (Baseline) with n` = sprintf('%s<br>(<i>n</i> = %s)', `AR-V7 (Baseline)`, dplyr::n_distinct(`Subject Number`))
    ) %>%
    dplyr::ungroup()

# Number of patients per Z-score grouping.
overviewPatients <- overviewPatients %>%
    dplyr::group_by(`Genome-wide status (Baseline)`) %>%
    dplyr::mutate(
        `Genome-wide status (Baseline)` = sprintf('%s<br>(<i>n</i> = %s)', `Genome-wide status (Baseline)`, dplyr::n_distinct(`Subject Number`))
    ) %>%
    dplyr::ungroup()


# Import Fast-Seq ---------------------------------------------------------

# Import T1.
dataFastSeq.T1 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 6) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'L-code', value.name = 'value.T1') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T1 = as.numeric(value.T1),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `L-code`))

# Import T2.
dataFastSeq.T2 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 7) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'Follow up L-code', value.name = 'value.T2') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T2 = as.numeric(value.T2),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `Follow up L-code`))

# Merge T1 and T2
dataFastSeq.Combined <- dataFastSeq.T1 %>%
    dplyr::left_join(dataFastSeq.T2, by = c('Chromosomal arm', 'Subject Number')) %>%
    dplyr::mutate(
        chr = gsub('p|q', '', `Chromosomal arm`),
        `Chromosomal arm` = factor(`Chromosomal arm`, levels = gtools::mixedsort(unique(`Chromosomal arm`))),
        value.deltaT2vsT1 = value.T2 - value.T1
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `CTC Count (Baseline – 7.5mL)`, `AR-V7 (Baseline)`, `AR-V7 Conversion`, `Genome-wide status (Baseline)`, `Genome-Wide Z Score (T2)`, `Inclusion (Treated with Caba)`), by = 'Subject Number') %>%

    # Add # of samples per group.
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::mutate(`AR-V7 (Baseline)` = sprintf('%s (<i>n</i> = %s)', `AR-V7 (Baseline)`, dplyr::n_distinct(`Subject Number`))) %>%
    dplyr::ungroup() %>%

    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(`AR-V7 Conversion` = sprintf('%s (<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', ifelse(`AR-V7 Conversion` == 'Pos.', 'AR-V7<sup>Pos.</sup>', 'Und.')), dplyr::n_distinct(`Subject Number`[!is.na(`Follow up L-code`)]))) %>%
    dplyr::ungroup() %>%

    # Sort.
    dplyr::mutate(`AR-V7 (Baseline)` = factor(`AR-V7 (Baseline)`, levels = c('AR-V7<sup>Pos.</sup> (<i>n</i> = 38)', 'AR-V7<sup>Neg.</sup> (<i>n</i> = 59)', 'AR-V7<sup>Und.</sup> (<i>n</i> = 34)')))


# Plot templates ----------------------------------------------------------

halfhalfplots <- function(){
    list(
        gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = T, show.legend = F),
        gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black'),
        ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5)
    )
}

# List to contain the various figures.
plots <- list()


# Main. Fig. 2 - Correlation AR-V7 vs. CTC --------------------------------

## AR-V7 status vs. CTC counts  ----

# Do statistical tests.
stat.test <- overviewPatients %>%
    dplyr::mutate(g = `AR-V7 (Baseline) with n`, value = `CTC Count (Baseline)`) %>%
    rstatix::pairwise_wilcox_test(value ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))



plots$ARv7vsCTC <- overviewPatients %>%
    dplyr::group_by(`AR-V7 (Baseline) with n`) %>%
    dplyr::mutate(medianCTC = median(`CTC Count (Baseline)`)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., aes(x = reorder(`AR-V7 (Baseline) with n`, -medianCTC), y = `CTC Count (Baseline)`, fill = `AR-V7 (Baseline) with n`, label = medianCTC, group = `AR-V7 (Baseline) with n`)) +

    # Add half-half plots with median labels.
    halfhalfplots() +

    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test, y.position = 9, step.increase = .02, tip.length = .01) +
    ggplot2::scale_fill_manual(values = c('#648FFF', '#FE6100', '#4D4D4D'), guide = F) +
    ggplot2::labs(x = 'AR-V7 determination (Baseline)', y = 'CTC Count (Baseline)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(-.1, 15000), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000)) +
    theme_Job


## Correlation between Genome-Wide Z Score and CTC count ----

plots$ZScorevsCTCCorr <- overviewPatients %>%
    dplyr::mutate(
        baselineAR = `AR-V7 (Baseline)`,
        groupStatus = ifelse(baselineAR == 'AR-V7<sup>Pos.</sup>' & (!is.na(`AR-V7 Conversion`) & `AR-V7 Conversion` == 'Neg.'), 'AR-V7<sup>Conv.</sup>', baselineAR)
    ) %>%
    dplyr::filter(`Genome-Wide Z Score (Baseline)` != '.') %>%
    dplyr::mutate(`Genome-Wide Z Score (Baseline)` = as.numeric(`Genome-Wide Z Score (Baseline)`)) %>%
    ggplot2::ggplot(., aes(x = `Genome-Wide Z Score (Baseline)`, y = `CTC Count (Baseline – 7.5mL)`, fill = groupStatus, group = '1')) +
    ggpmisc::stat_fit_glance(method = "cor.test", label.y = "top", method.args = list(formula = ~ x + y, method = "spearman", exact = FALSE), mapping = aes(label = sprintf('rho~"="~%.3f~~italic(P)~"="~%.2g',stat(estimate), stat(p.value))),parse = TRUE) +     ggplot2::geom_point(shape = 21) +
    ggplot2::geom_vline(aes(xintercept = 5), lty = '11', color = '#D00103') +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Pos.</sup>' = '#FE6100', 'AR-V7<sup>Neg.</sup>' = '#648FFF', 'AR-V7<sup>Und.</sup>' = '#4D4D4D', 'AR-V7<sup>Conv.</sup>' = '#00A94D'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 750), limits = c(-.5, 800), expand = c(0,0)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000), limits = c(-.5, 6000), expand = c(0,0)) +
    ggplot2::geom_smooth(color = 'red', method = 'lm') +
    ggplot2::labs(y = 'CTC Count (Baseline – 7.5mL)', x = 'Aneuploidy score (Baseline)<br>(i.e. Genome-wide Z-Score)') +
    theme_Job


## Fast-Seq genome-wide Z-Score (High / Low) vs. AR Status  ----

stat.test <- overviewPatients %>%
    dplyr::filter(!grepl('\\.', `Genome-wide status (Baseline)`)) %>%
    dplyr::mutate(`Genome-Wide Z Score (Baseline)` = as.numeric(`Genome-Wide Z Score (Baseline)`)) %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::mutate(`AR-V7 (Baseline)` = sprintf('%s<br>(<i>n</i> = %s)', gsub('\\(.*', '', `AR-V7 (Baseline)`), dplyr::n_distinct(`Subject Number`))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(g = `AR-V7 (Baseline)`, value = `Genome-Wide Z Score (Baseline)`) %>%
    rstatix::wilcox_test(value ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

plots$ZvsARStatus <- overviewPatients %>%
    dplyr::filter(!grepl('\\.', `Genome-wide status (Baseline)`)) %>%
    dplyr::mutate(`Genome-Wide Z Score (Baseline)` = as.numeric(`Genome-Wide Z Score (Baseline)`)) %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::mutate(`AR-V7 (Baseline)` = sprintf('%s<br>(<i>n</i> = %s)', gsub('\\(.*', '', `AR-V7 (Baseline)`), dplyr::n_distinct(`Subject Number`))) %>%
    dplyr::mutate(medianZ = round(median(`Genome-Wide Z Score (Baseline)`), 2)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., aes(x = reorder(`AR-V7 (Baseline)`, -medianZ), y = `Genome-Wide Z Score (Baseline)`, fill = `AR-V7 (Baseline)`, label = medianZ, group = `AR-V7 (Baseline)`)) +

    # Add half-half plots with median labels.
    halfhalfplots() +

    ggplot2::labs(x = 'AR-V7 determination (Baseline)<br>(Only samples with mFAST-SeqS)', y = 'Aneuploidy score (Baseline)<br>(i.e. Genome-wide Z-Score)') +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test, y.position = 7, step.increase = .02, tip.length = .01) +
    ggplot2::scale_fill_manual(values = c('#648FFF', '#FE6100', '#4D4D4D'), guide = F) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(-.1, 1500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000)) +
    theme_Job


## Combine plots  ----

layout <- "AB
CC"

plots$ARv7vsCTC + plots$ZvsARStatus + plots$ZScorevsCTCCorr +
    patchwork::plot_layout(design = layout, guides = 'keep') +
    patchwork::plot_annotation(tag_levels = 'a')


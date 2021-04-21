# Author:                      Job van Riet
# Date:                        03-03-2021
# Function:                    Generate supplementary figures on the CABA-V7 cohort.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Helper theme.
theme_Job <- ggplot2::theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
    axis.text.x = ggtext::element_markdown(),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
    legend.text = ggtext::element_markdown()
)


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
    ggpmisc::stat_fit_glance(method = "cor.test",label.y = "top", method.args = list(formula = ~ x + y, method = "spearman", exact = FALSE), mapping = aes(label = sprintf('rho~"="~%.3f~~italic(P)~"="~%.2g',stat(estimate), stat(p.value))),parse = TRUE) +     ggplot2::geom_point(shape = 21) +
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


# Sup. Fig. 1 - mFAST-SeqS ------------------------------------------------

## Z-scores per arm per AR-V7 status (T1) ----
stat.test.T1 <- dataFastSeq.Combined %>%
    dplyr::filter(!`Chromosomal arm` %in% c('13p', '14p', '15p', '22p')) %>%
    dplyr::filter(!grepl('Und.', `AR-V7 (Baseline)`)) %>%
    dplyr::mutate(g = `AR-V7 (Baseline)`) %>%
    dplyr::group_by(`Chromosomal arm`) %>%
    rstatix::wilcox_test(value.T1 ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns')) %>%
    dplyr::filter(p.adj.signif != 'ns')

dataFastSeq.Combined %>%
    dplyr::filter(`Chromosomal arm` %in% c('8q', '13q', '18p', 'Xq')) %>%
    dplyr::group_by(`Chromosomal arm`, `AR-V7 (Baseline)`) %>%
    dplyr::mutate(median = round(median(value.T1), 1)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., aes(x = `AR-V7 (Baseline)`, y = value.T1, fill = `AR-V7 (Baseline)`, color = `AR-V7 (Baseline)`, label = median)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    gghalves::geom_half_boxplot(outlier.shape = NA, color = 'black') +
    gghalves::geom_half_point_panel(size = .1, alpha = .5) +
    ggplot2::scale_fill_manual(values = c('#FE6100', '#648FFF', '#4D4D4D'), guide = F) +
    ggplot2::scale_color_manual(values = c('#FE6100', '#648FFF', '#4D4D4D'), guide = F) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150), limits = c(-50, 200)) +
    ggplot2::labs(x = 'Chromosomal Arm', y = 'Aneuploidy score / Z-score<br>(Baseline (T1))') +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test.T1, y.position = 5,  tip.length = .01) +
    ggplot2::stat_summary(fun = min, colour='black', geom='text', size = 3, show.legend = FALSE, vjust = .5, angle = 90, hjust = 1.5) +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(.~`Chromosomal arm`, nrow = 2)


plots$ArmsvsARv7T1 <- dataFastSeq.Combined %>%
    ggplot2::ggplot(., aes(x = `AR-V7 (Baseline)`, y = value.T1, fill = `AR-V7 (Baseline)`, color = `AR-V7 (Baseline)`)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    gghalves::geom_half_boxplot(outlier.shape = NA, color = 'black') +
    gghalves::geom_half_point_panel(size = .1, alpha = .5) +
    ggplot2::scale_fill_manual(values = c('#FE6100', '#648FFF', '#4D4D4D'), guide = F) +
    ggplot2::scale_color_manual(values = c('#FE6100', '#648FFF', '#4D4D4D'), guide = F) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150), limits = c(-50, 200)) +
    ggplot2::labs(x = 'Chromosomal Arm', y = 'Aneuploidy score / Z-score<br>(Baseline (T1))') +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test.T1, y.position = 4,  tip.length = .01) +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(.~`Chromosomal arm`, nrow = 2)

## Z-scores per arm per AR-V7 status (T2).
stat.test.T2 <- dataFastSeq.Combined %>%
    dplyr::filter(!is.na(`Follow up L-code`), !grepl('NA', `AR-V7 Conversion`)) %>%
    dplyr::filter(`Chromosomal arm` %in% c('8q', '13q')) %>%
    dplyr::mutate(g = `AR-V7 Conversion`) %>%
    dplyr::group_by(`Chromosomal arm`) %>%
    rstatix::wilcox_test(value.T2 ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns')) %>%
    dplyr::filter(p.adj.signif != 'ns')

plots$ArmsvsARv7T2 <- dataFastSeq.Combined %>%
    dplyr::filter(!is.na(`Follow up L-code`), !grepl('NA', `AR-V7 Conversion`)) %>%
    ggplot2::ggplot(., aes(x = `AR-V7 Conversion`, y = value.T2, fill = `AR-V7 Conversion`, color = `AR-V7 Conversion`)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    gghalves::geom_half_boxplot(outlier.shape = NA, color = 'black') +
    gghalves::geom_half_point_panel(size = .1, alpha = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150), limits = c(-50, 200)) +
    ggplot2::labs(x = 'Chromosomal Arm', y = 'Aneuploidy score / Z-score<br>(post-cabazitaxel (T2))') +
    #ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test.T2, y.position = 4.9,  tip.length = .01) +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(.~`Chromosomal arm`, nrow = 2)


## Delta Z-scores per arm between Lost / Retained AR-V7 (T2 - T1)  ----
stat.test.pairedZ <- dataFastSeq.Combined %>%
    dplyr::filter(!is.na(`Follow up L-code`), !grepl('NA', `AR-V7 Conversion`)) %>%
    dplyr::filter(`Chromosomal arm` %in% c('8q', '13q')) %>%
    dplyr::mutate(g = `AR-V7 Conversion`) %>%
    dplyr::group_by(`Chromosomal arm`) %>%
    rstatix::wilcox_test(value.deltaT2vsT1 ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

plots$ArmsvsARv7T2vsT1 <- dataFastSeq.Combined %>%
    dplyr::filter(!is.na(`Follow up L-code`), !grepl('NA', `AR-V7 Conversion`)) %>%
    ggplot2::ggplot(., aes(x = `AR-V7 Conversion`, y = value.deltaT2vsT1, fill = `AR-V7 Conversion`, color = `AR-V7 Conversion`)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    gghalves::geom_half_boxplot(outlier.shape = NA, color = 'black') +
    gghalves::geom_half_point_panel(size = .1, alpha = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-15, -7.5, -2.5, 0, 2.5, 7.5, 15), limits = c(-17.5, 17.5)) +
    ggplot2::labs(x = 'Chromosomal Arm', y = 'Aneuploidy score / Z-score<br>(pre- vs. post-cabazitaxel (T2 - T1))') +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test.pairedZ, y.position = 3,  tip.length = .01) +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(.~`Chromosomal arm`, nrow = 2)

## Combine Fast-Seq plots  ----
plots$ArmsvsARv7T1 + plots$ArmsvsARv7T2 + plots$ArmsvsARv7T2vsT1 +
    patchwork::plot_layout(guides = 'keep', nrow = 3) +
    patchwork::plot_annotation(tag_levels = 'a')


# Sup. Fig. 2 - Converters ------------------------------------------------

plots.SupFig2 <- list()

## Compare T1 vs. T2 within (non-)converters ----

dataFastSeq.Paired <- dataFastSeq.Combined %>%
    dplyr::filter(!is.na(`Follow up L-code`), !grepl('NA', `AR-V7 Conversion`)) %>%
    dplyr::filter(`Chromosomal arm` %in% c('8q', '13q')) %>%
    dplyr::select(`Chromosomal arm`, `Subject Number`, `AR-V7 Conversion`, value.T1, value.T2) %>%
    reshape2::melt(id.vars = c('Chromosomal arm', 'Subject Number', 'AR-V7 Conversion')) %>%
    dplyr::mutate(
        variable = ifelse(variable == 'value.T1', 'Pre-cabazitaxel (T1)', 'Post-cabazitaxel (T2)'),
        variable = factor(variable, levels = c('Pre-cabazitaxel (T1)', 'Post-cabazitaxel (T2)'))
    ) %>%
    dplyr::group_by(variable, `Chromosomal arm`, `AR-V7 Conversion`) %>%
    dplyr::mutate(median = round(median(value), 1)) %>%
    dplyr::ungroup()

stat.test.pairedZ <- dataFastSeq.Paired %>%
    dplyr::group_by(`Chromosomal arm`, `AR-V7 Conversion`) %>%
    rstatix::wilcox_test(value ~ variable, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = T) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

plots.SupFig2$sigArms <- dataFastSeq.Paired %>%
    ggplot2::ggplot(., aes(x = variable, y = value, fill = `AR-V7 Conversion`, color = `AR-V7 Conversion`, label = median)) +
    ggplot2::geom_hline(yintercept = 0, color = 'red', lty = 'dashed') +
    gghalves::geom_half_boxplot(outlier.shape = NA, center = T, side = 'l', alpha = .5, color = 'black') +
    ggplot2::geom_point(size = 1.5, alpha = 1, position = ggplot2::position_nudge(x = .25)) +
    ggplot2::geom_line(aes(group = `Subject Number`), position =  ggplot2::position_nudge(x = .25)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150), limits = c(-50, 50)) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = F) +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-3, angle = 90, hjust = .5) +
    ggplot2::labs(x = NULL, y = 'Aneuploidy score / Z-score') +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(`Chromosomal arm` ~ `AR-V7 Conversion`, nrow = 1)


## Compare CTC-counts vs. converter status ----

dataConv.CTC <- overviewPatients %>%
    dplyr::filter(`Inclusion (Treated with Caba)` == 'Yes', !is.na(`AR-V7 Conversion`)) %>%
    dplyr::select(`CTC Count (Baseline – 7.5mL)`, `AR-V7 Conversion`, `Subject Number`) %>%
    dplyr::distinct() %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(
        g = sprintf('%s<br>(<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', 'AR-V7<sup>Pos.</sup>'), dplyr::n_distinct(`Subject Number`)),
        v = `CTC Count (Baseline – 7.5mL)`,
        median = round(median(v), 1)
    ) %>%
    dplyr::ungroup()

stat.test.convCTC <- dataConv.CTC %>%
    rstatix::wilcox_test(v ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

plots.SupFig2$convCTC <- dataConv.CTC %>%
    ggplot2::ggplot(., aes(x = reorder(g, -median), y = `CTC Count (Baseline – 7.5mL)`, fill = g, label = median)) +
    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +

    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 10)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 10)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::labs(x = 'AR-V7 determination<br>(post-cabazitaxel; T2)', y = 'CTC Count (Baseline)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(0, 2500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500)) +
    theme_Job

## CTC at treatment 4----

dataConv.CTC.T2 <- overviewPatients %>%
    dplyr::filter(`Inclusion (Treated with Caba)` == 'Yes', !is.na(`AR-V7 Conversion`)) %>%
    dplyr::select(`CTC – Treatment 4`, `AR-V7 Conversion`, `Subject Number`) %>%
    dplyr::filter(!is.na(`CTC – Treatment 4`)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(
        g = sprintf('%s<br>(<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', 'AR-V7<sup>Pos.</sup>'), dplyr::n_distinct(`Subject Number`)),
        v = `CTC – Treatment 4`,
        median = round(median(v, na.rm = T), 1)
    ) %>%
    dplyr::ungroup()

stat.test.convCTC.T2 <- dataConv.CTC.T2 %>%
    rstatix::wilcox_test(v ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

plots.SupFig2$convCTC.T2 <- dataConv.CTC.T2 %>%
    ggplot2::ggplot(., aes(x = reorder(g, -median), y = `CTC – Treatment 4`, fill = g, label = median)) +
    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +

    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 9)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 9)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::labs(x = 'AR-V7 determination<br>(post-cabazitaxel; T2)', y = 'CTC Count (T2)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(0, 2500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500)) +
    theme_Job

## Combine Converter plots -----
plots.SupFig2$sigArms + plots.SupFig2$convCTC + plots.SupFig2$convCTC.T2 +
    patchwork::plot_layout(guides = 'keep', nrow = 1, widths = c(1, .25, .25)) +
    patchwork::plot_annotation(tag_levels = 'a')

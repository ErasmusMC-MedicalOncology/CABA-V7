# Author:                      Job van Riet
# Date:                        28-04-2021
# Function:                    Determine mFASTSeqS deviations between T1 and T2.

# Import libraries --------------------------------------------------------

library(dplyr)
library(ggplot2)

# Helper theme.
source('R/misc_themes.R')

# Import data -------------------------------------------------------------

## Import patient metadata ----

overviewPatients <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 1)

## Import mFAST-Seq - T1 ----
dataFastSeq.T1 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 6) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'L-code', value.name = 'value.T1') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T1 = as.numeric(value.T1),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `L-code`))

## Import mFAST-Seq - T2 ----
dataFastSeq.T2 <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 7) %>%
    reshape2::melt(id = 'Chromosomal arm', variable.name = 'Follow up L-code', value.name = 'value.T2') %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        value.T2 = as.numeric(value.T2),
        `Chromosomal arm` = gsub('chr', '', `Chromosomal arm`)
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `Follow up L-code`))

## Merge mFAST-Seq T1 and T2 ----
dataFastSeq.Combined <- dataFastSeq.T1 %>%
    dplyr::inner_join(dataFastSeq.T2, by = c('Chromosomal arm', 'Subject Number')) %>%
    dplyr::mutate(
        chr = gsub('p|q', '', `Chromosomal arm`),
        `Chromosomal arm` = factor(`Chromosomal arm`, levels = gtools::mixedsort(unique(`Chromosomal arm`))),
        value.deltaT2vsT1 = value.T2 - value.T1
    ) %>%
    dplyr::inner_join(overviewPatients %>% dplyr::select(`Subject Number`, `AR-V7 (Baseline)`, `AR-V7 Conversion`, `Genome-Wide Z Score (Baseline)`, `Genome-Wide Z Score (T2)`), by = 'Subject Number') %>%

    # Remove samples with no first and second AR-V7 determination.
    dplyr::filter(!is.na(`AR-V7 Conversion`), !is.na(`AR-V7 (Baseline)`)) %>%

    # Add nr. of patients per conversion group.
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(`AR-V7 Conversion` = sprintf('%s (<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', ifelse(`AR-V7 Conversion` == 'Pos.', 'AR-V7<sup>Pos.</sup>', 'Und.')), dplyr::n_distinct(`Subject Number`))) %>%
    dplyr::ungroup() %>%

    # Sort on Pos and then Conv patients.
    dplyr::mutate(`AR-V7 Conversion` = factor(`AR-V7 Conversion`, levels = c('AR-V7<sup>Pos.</sup> (<i>n</i> = 8)', 'AR-V7<sup>Conv.</sup> (<i>n</i> = 8)'))) %>%

    # Remove acrocentric chromosomes which we cannot measure
    dplyr::filter(!is.na(value.deltaT2vsT1)) %>%

    # Fix numerics.
    dplyr::mutate(`Genome-Wide Z Score (Baseline)` = as.numeric(`Genome-Wide Z Score (Baseline)`))


# Sup. Fig. 2 - mFAST-SeqS on T1 ------------------------------------------------

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


# Sup. Fig 3 - Differences between T1 and T2 ----------------------

dataFastSeq.Combined %>%
    dplyr::select(`Chromosomal arm`, `Subject Number`, `AR-V7 Conversion`, value.T1, value.T2) %>%
    dplyr::mutate(
        value.T2.arrow = ifelse(value.T2 > value.T1, value.T2 - 1, value.T2 + 1),
        value.T2.arrow.col = ifelse(value.T2 >= value.T1, 'ZΔ (T2 - T1): ≥1', 'ZΔ (T2 - T1): ≤1'),
        deltaZ = abs(value.T2 - value.T1),
        `Subject Number` = factor(`Subject Number`, levels = dataFastSeq.Combined %>% dplyr::distinct(`Subject Number`, `AR-V7 Conversion`) %>% dplyr::arrange(`AR-V7 Conversion`) %>% dplyr::pull(`Subject Number`))
    ) %>%

    ggplot2::ggplot(., aes(x = `Chromosomal arm`, xend = `Chromosomal arm`, group = `Subject Number`)) +
    # Points and arrows.
    ggplot2::geom_point(aes(y = `value.T2`), fill = '#3E3EF7', size = 1.25, shape = 23) +
    ggplot2::geom_segment(data = . %>% dplyr::filter(deltaZ >= 1), aes(y = value.T1, yend = value.T2.arrow, color = value.T2.arrow.col), linejoin = 'mitre', size = .5, arrow = ggplot2::arrow(type = 'closed', length = ggplot2::unit(0.1, "cm"))) +
    ggplot2::geom_point(aes(y = `value.T1`), fill = '#B9B9F7', size = 1.25, shape = 21) +
    # Guidelines.
    ggplot2::geom_hline(yintercept = 0, lwd = .2, color = 'red', lty = 'dashed') +
    ggplot2::geom_hline(yintercept = 5, lwd = .2, color = 'orange', lty = 'dashed') +
    ggplot2::geom_hline(yintercept = -5, lwd = .2, color = 'orange', lty = 'dashed') +
    # Themes.
    ggplot2::scale_color_manual(values = c('ZΔ (T2 - T1): ≥1' = '#4FBE4F', 'ZΔ (T2 - T1): ≤1' = '#E14055'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::scale_y_continuous(limits = c(-26, 50), trans = scales::pseudo_log_trans(), breaks = c(-150, -100, -50, -25, -10, -5, -2.5, 0, 2.5, 5, 10, 25, 50, 150)) +
    ggplot2::labs(x = 'Chromosomal Arm', y = 'Aneuploidy score / Z-score') +
    theme_Job + ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1), panel.grid.major.x = element_blank()) +
    ggplot2::facet_wrap(. ~ `Subject Number`, ncol = 2, dir = 'v')


# Genome-wide Z-Score per timepoint ---------------------------------------

statTests <- list()
plots.FastSeq <- list()

# Statistical tests. ----
statTests$genomeWideZ.T1 <- dataFastSeq.Combined %>%
    dplyr::distinct(`Subject Number`, `Genome-Wide Z Score (Baseline)`, `AR-V7 Conversion`) %>%
    dplyr::mutate(g = `AR-V7 Conversion`, v  = `Genome-Wide Z Score (Baseline)`) %>%
    rstatix::wilcox_test(v ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

statTests$genomeWideZ.T2 <- dataFastSeq.Combined %>%
    dplyr::distinct(`Subject Number`, `Genome-Wide Z Score (T2)`, `AR-V7 Conversion`) %>%
    dplyr::mutate(g = `AR-V7 Conversion`, v = `Genome-Wide Z Score (T2)`) %>%
    rstatix::wilcox_test(v ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

statTests$genomeWideZ.BetweenT1andT2 <- dataFastSeq.Combined %>%
    dplyr::distinct(`Subject Number`,  `Genome-Wide Z Score (Baseline)`, `Genome-Wide Z Score (T2)`, `AR-V7 Conversion`) %>%
    reshape2::melt(id.vars = c('Subject Number', 'AR-V7 Conversion')) %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    rstatix::wilcox_test(value ~ variable, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = T) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

## Between baseline Pos. vs. Conv. ----

plots.FastSeq$ZT1.Pos_Convdata <- dataFastSeq.Combined %>%
    dplyr::distinct(`Subject Number`, `Genome-Wide Z Score (Baseline)`, `AR-V7 Conversion`) %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(median = round(median(`Genome-Wide Z Score (Baseline)`), 1)) %>%
    dplyr::ungroup() %>%

    ggplot2::ggplot(., aes(x = `AR-V7 Conversion`, y = `Genome-Wide Z Score (Baseline)`, fill = `AR-V7 Conversion`, label = median)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), limits = c(-1, 130), breaks = c(-1, 0:5, 10, 20, 50, 100, 125)) +
    gghalves::geom_half_boxplot(outlier.shape = NA, alpha = .5, color = 'black') +
    gghalves::geom_half_point_panel(size = 2, position = ggbeeswarm::position_quasirandom(width = .2), shape = 21) +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::labs(x = 'AR-V7 determination<br>(Baseline; T1)', y = 'Genome-wide Z-score<br>Baseline (log<sub>10</sub>)') +
    theme_Job


## Between T2 Pos. vs. Conv. ----

plots.FastSeq$ZT2.Pos_Conv <- dataFastSeq.Combined %>%
    dplyr::distinct(`Subject Number`, `Genome-Wide Z Score (T2)`, `AR-V7 Conversion`) %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(median = round(median(`Genome-Wide Z Score (T2)`), 1)) %>%
    dplyr::ungroup() %>%

    ggplot2::ggplot(., aes(x = `AR-V7 Conversion`, y = `Genome-Wide Z Score (T2)`, fill = `AR-V7 Conversion`, label = median)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), limits = c(-1, 130), breaks = c(-1, 0:5, 10, 20, 50, 100, 125)) +
    gghalves::geom_half_boxplot(outlier.shape = NA, alpha = .5, color = 'black') +
    gghalves::geom_half_point_panel(size = 2, position = ggbeeswarm::position_quasirandom(width = .2), shape = 21) +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::labs(x = 'AR-V7 determination<br>(Post-cabazitaxel; T2)', y = 'Genome-wide Z-score<br>Post-cabazitaxel (log<sub>10</sub>)') +
    theme_Job

## Between T1 vs T2 for Pos. vs. Conv. ----

plots.FastSeq$BothT.Pos_Conv <- dataFastSeq.Combined %>%
    dplyr::select(`Subject Number`, `Baseline (T1)` = `Genome-Wide Z Score (Baseline)`, `Post-cabazitaxel (T2)` = `Genome-Wide Z Score (T2)`, `AR-V7 Conversion`) %>%
    dplyr::distinct() %>%
    reshape2::melt(id.vars = c('Subject Number', 'AR-V7 Conversion')) %>%
    dplyr::group_by(variable, `AR-V7 Conversion`) %>%
    dplyr::mutate(median = round(median(value), 1)) %>%
    dplyr::ungroup() %>%

    ggplot2::ggplot(., aes(x = variable, y = value, fill = `AR-V7 Conversion`, label = median)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), limits = c(-1, 130), breaks = c(-1, 0:5, 10, 20, 50, 100, 125)) +
    gghalves::geom_half_boxplot(outlier.shape = NA, center = T, side = 'l', alpha = .5, color = 'black') +
    ggplot2::geom_point(size = 1.5, alpha = 1, position = ggplot2::position_nudge(x = .25)) +
    ggplot2::geom_line(aes(group = `Subject Number`), position =  ggplot2::position_nudge(x = .25)) +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust = -4, angle = 90, hjust = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 8)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 8)' = '#FE6100'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::labs(x = 'Timepoints', y = 'Genome-wide Z-score (log<sub>10</sub>)') +
    theme_Job + facet_wrap(. ~ `AR-V7 Conversion`)


### Combine plots ----
plots.FastSeq$ZT1.Pos_Convdata + plots.FastSeq$ZT2.Pos_Conv + plots.FastSeq$BothT.Pos_Conv +
    patchwork::plot_layout(guides = 'keep', ncol = 3, widths = c(1, 1, 2)) +
    patchwork::plot_annotation(tag_levels = 'a')


## Compare CTC-counts vs. converter status ----

plots.CTCs <- list()

# Get CTCs of patients.
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

# Test CTC between Pos/Conv on Baseline
stat.test.convCTC <- dataConv.CTC %>%
    rstatix::wilcox_test(v ~ g, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

## Figure - CTC (Baseline) between Pos. / Conv. ----

plots.CTCs$convCTC <- dataConv.CTC %>%
    ggplot2::ggplot(., aes(x = reorder(g, -median), y = `CTC Count (Baseline – 7.5mL)`, fill = g, label = median)) +
    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', alpha = .5, outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +

    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 10)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 10)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::labs(x = 'AR-V7 determination<br>(post-cabazitaxel; T2)', y = 'CTC Count (Baseline)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(0, 2500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500)) +
    theme_Job


## Figure - CTC (T2) between Pos. / Conv. ----

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

plots.CTCs$convCTC.T2 <- dataConv.CTC.T2 %>%
    ggplot2::ggplot(., aes(x = reorder(g, -median), y = `CTC – Treatment 4`, fill = g, label = median)) +
    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', alpha = .5, outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +

    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 9)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::scale_color_manual(values = c('AR-V7<sup>Conv.</sup><br>(<i>n</i> = 9)' = '#00A94D', 'AR-V7<sup>Pos.</sup><br>(<i>n</i> = 9)' = '#FE6100'), guide = F) +
    ggplot2::labs(x = 'AR-V7 determination<br>(post-cabazitaxel; T2)', y = 'CTC Count (T2)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(0, 2500), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500)) +
    theme_Job

## Figure - CTC T1 -> T2 per Pos. / Conv. ----

dataConv.CTCBetweenT1AndT2 <- overviewPatients %>%
    dplyr::filter(`Inclusion (Treated with Caba)` == 'Yes', !is.na(`AR-V7 Conversion`), !is.na(`CTC – Treatment 4`)) %>%
    dplyr::select(`CTC Count (Baseline – 7.5mL)`, `CTC – Treatment 4`, `AR-V7 Conversion`, `Subject Number`) %>%
    dplyr::distinct() %>%
    dplyr::group_by(`AR-V7 Conversion`) %>%
    dplyr::mutate(
        g = sprintf('%s (<i>n</i> = %s)', ifelse(`AR-V7 Conversion` == 'Neg.', 'AR-V7<sup>Conv.</sup>', 'AR-V7<sup>Pos.</sup>'), dplyr::n_distinct(`Subject Number`)),
    ) %>%
    dplyr::ungroup()

dataConv.CTCBetweenT1AndT2 %>%
    reshape2::melt(id.vars = c('g', 'AR-V7 Conversion', 'Subject Number')) %>%
    dplyr::group_by(g) %>%
    rstatix::wilcox_test(value ~ variable, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = T) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))


dataConv.CTCBetweenT1AndT2 %>%
    reshape2::melt(id.vars = c('g', 'AR-V7 Conversion', 'Subject Number')) %>%
    dplyr::group_by(variable, g) %>%
    dplyr::mutate(median = round(median(value), 1)) %>%
    dplyr::ungroup() %>%

    ggplot2::ggplot(., aes(x = variable, y = value, fill = g, label = median)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), limits = c(-.2, 2500), breaks = c(0:5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000)) +
    gghalves::geom_half_boxplot(outlier.shape = NA, center = T, side = 'l', alpha = .5, color = 'black') +
    ggplot2::geom_point(size = 1.5, alpha = 1, position = ggplot2::position_nudge(x = .25)) +
    ggplot2::geom_line(aes(group = `Subject Number`), position =  ggplot2::position_nudge(x = .25)) +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust = -4, angle = 90, hjust = .5) +
    ggplot2::scale_fill_manual(values = c('AR-V7<sup>Conv.</sup> (<i>n</i> = 9)' = '#00A94D', 'AR-V7<sup>Pos.</sup> (<i>n</i> = 9)' = '#FE6100'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::labs(x = 'Timepoints', y = 'CTC Count (log<sub>10</sub>)') +
    theme_Job + facet_wrap(. ~ g)

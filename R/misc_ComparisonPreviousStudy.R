
# Compare against previous study ------------------------------------------

ctcDataPrevious <- readr::read_tsv('Misc./ctcPreviousStudies.txt')
ctcDataNew <- data.Patient$Overview %>% dplyr::select(`Subject Number`, `AR-V7 (Baseline)`, `CTC Count (Baseline – 7.5mL)`) %>% dplyr::mutate(Study = 'CABA-V7') %>% dplyr::filter(`Subject Number` != '42.2')

dataCombined <- dplyr::bind_rows(ctcDataNew, ctcDataPrevious)

plotCTC.All <- dataCombined %>%
    dplyr::filter(`AR-V7 (Baseline)` != 'Und.') %>%
    dplyr::filter(Study != 'PRELUDE') %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>% dplyr::mutate(median = median(`CTC Count (Baseline – 7.5mL)`)) %>% dplyr::ungroup() %>%
    ggplot(aes(x = `AR-V7 (Baseline)`, y = `CTC Count (Baseline – 7.5mL)`, fill = `AR-V7 (Baseline)`, label = median)) +
    geom_hline(aes(yintercept = 5), color = '#F78A02', lty = '11') +
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1) +
    scale_fill_manual(values = c('Pos.' = '#F61C54', 'Neg.' = '#0147B2', 'Und.' = '#E9B577'), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +
    scale_shape_manual(values = c('Und.' = 19, 'Pos.' = 24, 'Neg.' = 25), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000), limits = c(0, 10000)) +
    ggsignif::geom_signif(comparisons = list(c('Pos.', 'Neg.')), step_increase = .025, color = 'black', map_signif_level = F, test = 'wilcox.test', tip_length = .01) +
    theme_Job

plotCTC.AllMin5 <- dataCombined %>%
    dplyr::filter(`CTC Count (Baseline – 7.5mL)` >= 5, `AR-V7 (Baseline)` != 'Und.') %>%
    dplyr::filter(Study != 'PRELUDE') %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>% dplyr::mutate(median = median(`CTC Count (Baseline – 7.5mL)`)) %>% dplyr::ungroup() %>%
    ggplot(aes(x = `AR-V7 (Baseline)`, y = `CTC Count (Baseline – 7.5mL)`, fill = `AR-V7 (Baseline)`, color = Study, group = `AR-V7 (Baseline)`, label = median)) +
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1) +
    scale_fill_manual(values = c('Pos.' = '#F61C54', 'Neg.' = '#0147B2', 'Und.' = '#E9B577'), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +
    scale_shape_manual(values = c('Und.' = 19, 'Pos.' = 24, 'Neg.' = 25), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000),  limits = c(0, 10000)) +
    ggsignif::geom_signif(comparisons = list(c('Pos.', 'Neg.')), y_position = 9.1, step_increase = .025, color = 'black', map_signif_level = F, test = 'wilcox.test', tip_length = .01) +
    theme_Job

plot.CTCPerStudy <- dataCombined %>%
    dplyr::filter(`CTC Count (Baseline – 7.5mL)` >= 5, `AR-V7 (Baseline)` != 'Und.') %>%
    dplyr::group_by(`AR-V7 (Baseline)`, Study) %>% dplyr::mutate(median = median(`CTC Count (Baseline – 7.5mL)`)) %>% dplyr::ungroup() %>%
    ggplot(aes(x = `AR-V7 (Baseline)`, y = `CTC Count (Baseline – 7.5mL)`, fill = `AR-V7 (Baseline)`, label = median)) +
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1) +
    scale_fill_manual(values = c('Pos.' = '#F61C54', 'Neg.' = '#0147B2', 'Und.' = '#E9B577'), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-2, angle = 90, hjust = .5) +
    scale_shape_manual(values = c('Und.' = 19, 'Pos.' = 24, 'Neg.' = 25), guide = guide_legend(title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000),  limits = c(0, 10000)) +
    ggsignif::geom_signif(comparisons = list(c('Pos.', 'Neg.')), y_position = 9, step_increase = .025, color = 'black', map_signif_level = F, test = 'wilcox.test', tip_length = .01) +
    theme_Job +
    facet_grid(~ Study)

## Combine landscape tracks.
layout <- "AB
CC"

plotCTC.All + plotCTC.AllMin5 + plot.CTCPerStudy +
    patchwork::plot_layout(design = layout, nrow = 2, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a')


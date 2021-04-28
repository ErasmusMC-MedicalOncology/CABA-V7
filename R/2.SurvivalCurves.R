# Author:                      Job van Riet
# Date:                        03-03-2021
# Function:                    Generate Survival Curves for the CABA-V7 study.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(survminer)
library(extrafont)
library(patchwork)


# Functions ---------------------------------------------------------------

# Helper theme.
source('R/misc_themes.R')

# Generate survival plots with p-values and median OS.
plotSurvival <- function(fit, ylim, hr = NULL){

    # Generate survival plot.
    x <- survminer::ggsurvplot(
        fit = fit,
        pval = F,
        size = .5,
        break.time.by = 10,
        break.y.by = .2,
        palette = 'jco',
        risk.table = T,
        tables.height = .3,
        xlab = 'Time (in months)',
        axes.offset = F,
        ylim = c(0, 1.05),
        xlim = c(0, ylim), strata.labels = 'c',
        risk.table.col = 'strata', censor.shape = '+',
        fontsize = 3,
        risk.table.title = 'No. at risk',
        ggtheme = ggplot2::theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            legend.text = ggtext::element_markdown()
        )
    )

    # Add the log-rank p-value.
    p.logrank <- survminer::surv_pvalue(fit = fit, method = 'log-rank', data = data.Survival, test.for.trend = F)
    x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = 1, label = paste0('log-rank: ', p.logrank$pval.txt), size = 2.5)

    # Add HR (if two groups)
    if(!is.null(hr)){

        HR.CI <- round(summary(hr)$conf.int, 2)
        HR.p <- round(summary(hr)$waldtest[[3]], 2)
        HR.CI <- sprintf('HR (95% CI): %s (%s - %s)', HR.CI[[1]], HR.CI[[3]], HR.CI[[4]])
        x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = .9, label = HR.CI, size = 2.5)
    }

    # Add the median OS.
    medianOS <- x$data.survplot %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
            medianOS = round(median(time, na.rm = T), 2),
            label = sprintf('%s - %s mo', unique(strata), medianOS)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-medianOS)

    x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = .75, label = paste0('Median OS (Desc.):\n', paste(medianOS$label, collapse = '\n')), size = 2.5)

    return(x)
}


# Import data -------------------------------------------------------------

data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))


# Convert data ------------------------------------------------------------

data.Survival <- data.Patient$clinicalData %>%
    # Convert dates.
    dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        dateCensor = ifelse(!is.na(`Date: Death`), `Date: Death`, na.omit(c(`Date: Last follow-up`, `Date: End of study`, `Date: Pre-screening`))),
        dateCensor = as.Date(dateCensor, origin = '1970-01-01'),
        daysFromPreScreeningToEnd = dateCensor - `Date: Pre-screening`,
        monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12),
        Survival = ifelse(is.na(Survival), 0, Survival),
    ) %>%
    dplyr::inner_join(data.Patient$Overview) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        `Inclusion (Treated with Caba)` = ifelse(is.na(`Inclusion (Treated with Caba)`), 'No', `Inclusion (Treated with Caba)`)
    ) %>%

    # Combine
    data.frame()

plotFits <- list()


#  Survival - All inclusion ------------------------------------------------

fit.AllInClusion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline. + Inclusion..Treated.with.Caba., data = data.Survival)
names(fit.AllInClusion$strata) <-  c('AR-V7 Neg.', 'AR-V7 Pos. (Not treated with caba. as part of CABA-V7 trial)', 'AR-V7 Pos. (Treated with caba. as part of -V7 trial)', 'AR-V7 Und.')

plotFits$AllInclusion <- plotSurvival(fit.AllInClusion, ylim = 51)

# Hazard Ratio of multiple groups.
fit.AllInClusion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival)
plotFits$hazardAll <- survminer::ggforest(fit.AllInClusion.hr,  data = data.Survival, noDigits = 2)


# Survival - Caba-V7 (Pos.) vs. Neg. / Und. -------------------------------

fit.PosVsNeg <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival %>% dplyr::mutate(AR.V7..Baseline. = ifelse(AR.V7..Baseline. %in% c('Neg.', 'Und.'), 'Neg. / Und.', 'Pos.')))
names(fit.PosVsNeg$strata) <-  c('AR-V7 (Neg. / Und.)', 'AR-V7 Pos.')

# Calculate HR.
fit.PosVsNeg.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival %>% dplyr::mutate(AR.V7..Baseline. = ifelse(AR.V7..Baseline. %in% c('Neg.', 'Und.'), 'Neg. / Und.', 'Pos.')))

plotFits$fit.PosVsNeg <- plotSurvival(fit.PosVsNeg, ylim = 45, hr = fit.PosVsNeg.hr)


# Survival - Caba-treated (per AR-V7 conversion) --------------------------

fit.BetweenConversion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.Survival %>% dplyr::filter(Inclusion..Treated.with.Caba. == 'Yes'))
names(fit.BetweenConversion$strata) <-  c('AR-V7 Conversion (Pos. -> Neg.)', 'AR-V7 Retainment (Pos. -> Pos.)')

# Calculate HR.
fit.BetweenConversion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.Survival %>% dplyr::filter(Inclusion..Treated.with.Caba. == 'Yes'))

# Plot.
plotFits$fit.BetweenConversion <- plotSurvival(fit.BetweenConversion, ylim = 21, hr = fit.BetweenConversion.hr)


# Combine plots. ----------------------------------------------------------

layout <- "AB
C#
DE
FG
HI"

plotFits$AllInclusion$plot + plotFits$hazardAll +
    plotFits$AllInclusion$table +
    plotFits$fit.PosVsNeg$plot + plotFits$fit.BetweenConversion$plot +
    plotFits$fit.PosVsNeg$table + plotFits$fit.BetweenConversion$table +
    patchwork::plot_layout(design = layout, heights = c(1, .25, 1, .25), ncol = 2, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a')


# Calculate Summary -------------------------------------------------------

data.Survival %>%
    dplyr::filter(Inclusion..Treated.with.Caba. == 'Yes') %>%
    dplyr::mutate(
        totalTrialTime = Date..End.of.cabazitaxel.treatment - Date..Start.of.cabazitaxel.treatment
    ) %>%
    dplyr::summarise(
        median(daysFromPreScreeningToEnd, na.rm = T),
        median(monthsFromPreScreeningToEnd, na.rm = T),
        median(totalTrialTime, na.rm = T),
        median(Cycles.of.Cabazitaxel, na.rm = T)
    )



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



# Sup. Fig. X - Correlation to OS -----------------------------------------

plots.SupFigX <- list()

## OS vs. grouped CTC-count (Baseline) ----

data.OSvsCTC <- data.Survival %>%
    dplyr::distinct(Subject.Number, Survival, monthsFromPreScreeningToEnd, CTC.Count..Baseline...7.5mL.) %>%
    dplyr::filter(Survival == 1) %>%
    dplyr::mutate(
        g = ifelse(CTC.Count..Baseline...7.5mL. >= 5, 'CTC count ≥5', 'CTC count <5'),
        monthsFromPreScreeningToEnd = as.numeric(monthsFromPreScreeningToEnd),
    ) %>%
    dplyr::group_by(g) %>%
    dplyr::mutate(
        median = round(median(monthsFromPreScreeningToEnd), 1),
        g = sprintf('%s<br>(<i>n</i> = %s)', g, dplyr::n_distinct(Subject.Number))
    ) %>%
    dplyr::ungroup()

stat.data.OSvsCTC <- data.OSvsCTC %>%
    rstatix::wilcox_test(monthsFromPreScreeningToEnd ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))


plots.SupFigX$OSvsCTC <- ggplot2::ggplot(data.OSvsCTC, aes(x = g, y = monthsFromPreScreeningToEnd, label = median, fill = g)) +

    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-5, angle = 90, hjust = .5) +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.data.OSvsCTC, y.position = 33, step.increase = .02, tip.length = .01) +

    ggplot2::scale_fill_manual(values = c('CTC count <5<br>(<i>n</i> = 26)' = 'black', 'CTC count ≥5<br>(<i>n</i> = 70)' = 'grey50'), guide = F) +
    ggplot2::labs(x = 'CTC Count (Baseline)<br>On non-censured patients (<i>n</i> = 96)', y = 'Overall Survival (months)') +
    ggplot2::scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
    theme_Job

## OS vs. grouped Z-score (Baseline) ----

data.OSvsZ <- data.Survival %>%
    dplyr::distinct(Subject.Number, Survival, monthsFromPreScreeningToEnd, Genome.wide.status..Baseline., Genome.Wide.Z.Score..Baseline.) %>%
    dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.', Survival == 1) %>%
    dplyr::mutate(
        Genome.Wide.Z.Score..Baseline. = as.numeric(Genome.Wide.Z.Score..Baseline.),
        monthsFromPreScreeningToEnd = as.numeric(monthsFromPreScreeningToEnd)
    ) %>%
    dplyr::group_by(Genome.wide.status..Baseline.) %>%
    dplyr::mutate(
        median = round(median(monthsFromPreScreeningToEnd), 1),
        g = sprintf('%s<br>(<i>n</i> = %s)', Genome.wide.status..Baseline., dplyr::n_distinct(Subject.Number))
    ) %>%
    dplyr::ungroup()

stat.data.OSvsZ <- data.OSvsZ %>%
    rstatix::wilcox_test(monthsFromPreScreeningToEnd ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))


plots.SupFigX$OSvsZ <- ggplot2::ggplot(data.OSvsZ, aes(x = g, y = monthsFromPreScreeningToEnd, label = median, fill = g)) +

    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', size = 1.25, color = 'black') +
    ggplot2::stat_summary(fun = median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-5, angle = 90, hjust = .5) +
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.data.OSvsZ, y.position = 33, step.increase = .02, tip.length = .01) +

    ggplot2::scale_fill_manual(values = c('Genome-wide Z-score <5<br>(<i>n</i> = 37)' = 'black', 'Genome-wide Z-score ≥5<br>(<i>n</i> = 56)' = 'grey50'), guide = F) +
    ggplot2::labs(x = 'Genome-wide status (Baseline)<br>On non-censured patients with mFAST-SeqS (<i>n</i> = 93)', y = 'Overall Survival (months)') +
    ggplot2::scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
    theme_Job


## Combine plots ----
plots.SupFigX$OSvsCTC + plots.SupFigX$OSvsZ +
    patchwork::plot_layout(nrow = 1, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a')

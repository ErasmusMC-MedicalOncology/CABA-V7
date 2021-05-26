# Author:                      Job van Riet
# Date:                        05-05-2021
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
plotSurvival <- function(fit, ylim, data, hr = NULL){

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
    p.logrank <- survminer::surv_pvalue(fit = fit, method = 'log-rank', data = data, test.for.trend = F)
    x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = 1, label = paste0('log-rank: ', p.logrank$pval.txt), size = 2.5)

    # Add HR (if two groups)
    if(!is.null(hr)){

        HR.CI <- round(summary(hr)$conf.int, 2)
        HR.p <- round(summary(hr)$waldtest[[3]], 2)
        HR.CI <- sprintf('HR (.95%% CI): %s (%s - %s)', HR.CI[[1]], HR.CI[[3]], HR.CI[[4]])
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

fit.AllInClusion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival)
names(fit.AllInClusion$strata) <-  c('AR-V7 Neg.', 'AR-V7 Pos.', 'AR-V7 Und.')

plotFits$AllInclusion <- plotSurvival(fit.AllInClusion, data = data.Survival, ylim = 51)

# Hazard Ratio of multiple groups.
fit.AllInClusion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival)
plotFits$hazardAll <- survminer::ggforest(fit.AllInClusion.hr, data = data.Survival, noDigits = 2)


# Survival - Caba-V7 (Pos.) vs. Neg. / Und. -------------------------------

data.PosVsNeg <- data.Survival %>% dplyr::mutate(AR.V7..Baseline. = ifelse(AR.V7..Baseline. %in% c('Neg.', 'Und.'), 'Neg. / Und.', 'Pos.'))

fit.PosVsNeg <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.PosVsNeg)
names(fit.PosVsNeg$strata) <-  c('AR-V7 (Neg. / Und.)', 'AR-V7 Pos.')

# Calculate HR.
fit.PosVsNeg.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.PosVsNeg)

plotFits$fit.PosVsNeg <- plotSurvival(fit.PosVsNeg, ylim = 45, data = data.PosVsNeg, hr = fit.PosVsNeg.hr)


# Survival - Caba-treated (per AR-V7 conversion) --------------------------

data.BetweenConversion <- data.Survival %>% dplyr::filter(Inclusion..Treated.with.Caba. == 'Yes')
fit.BetweenConversion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.BetweenConversion)
names(fit.BetweenConversion$strata) <-  c('AR-V7 Conversion (Pos. -> Neg.)', 'AR-V7 Retainment (Pos. -> Pos.)')

# Calculate HR.
fit.BetweenConversion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.BetweenConversion)

# Plot.
plotFits$fit.BetweenConversion <- plotSurvival(fit.BetweenConversion, ylim = 21, data = data.BetweenConversion, hr = fit.BetweenConversion.hr)


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


# Sup. Fig. 2 - CTC and Z against OS -----------------------------------------

plots.SupFigX <- list()

## OS vs. grouped CTC-count (Baseline) ----

data.CTC <- data.Survival %>% dplyr::mutate(ctcCategory = ifelse(CTC.Count..Baseline...7.5mL. >= 5, 'CTC count ≥5', 'CTC count <5'))
fit.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategory, data = data.CTC)
names(fit.CTC$strata) <-  c('CTC count <5<br>(Baseline)', 'CTC count ≥5<br>(Baseline)')

# Calculate HR.
fit.CTC.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategory, data = data.CTC)

plotFits$fit.CTC.hr <- plotSurvival(fit.CTC, ylim = 45, data = data.CTC, hr = fit.CTC.hr)


## OS vs. grouped Z-score (Baseline) ----

data.Z <- data.Survival %>% dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.') %>% dplyr::mutate(zCategory = ifelse(Genome.Wide.Z.Score..Baseline. >= 5, 'Genome-wide Z-score ≥5', 'Genome-wide Z-score <5'))
fit.Z <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ zCategory, data = data.Z)
names(fit.Z$strata) <-  c('Aneuploidy score <5<br>(Baseline)', 'Aneuploidy ≥5<br>(Baseline)')

# Calculate HR.
fit.Z.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ zCategory, data = data.Z)

plotFits$fit.Z.hr <- plotSurvival(fit.Z, ylim = 45, data = data.Z, hr = fit.Z.hr)


## Combine plots ----
plotFits$fit.CTC.hr$plot + plotFits$fit.Z.hr$plot +
    plotFits$fit.CTC.hr$table + plotFits$fit.Z.hr$table +
    patchwork::plot_layout(heights = c(1, .25), nrow = 2, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a')


# CTC / Z T1 vs. T2 -------------------------------------------------------

data.Z <- data.Survival %>%
    dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.') %>%
    dplyr::mutate(
        zCategory.T1 = ifelse(Genome.Wide.Z.Score..Baseline. >= 5, 'Genome-wide Z-score ≥5', 'Genome-wide Z-score <5'),
        zCategory.T2 = ifelse(Genome.Wide.Z.Score..T2. >= 5, 'Genome-wide Z-score ≥5', 'Genome-wide Z-score <5'),
        convStatus  = ifelse(zCategory.T1 == 'Genome-wide Z-score ≥5' & zCategory.T2 == 'Genome-wide Z-score ≥5', 'T1 and T2 > 5', 'Other'),
        convStatus  = ifelse(zCategory.T1 == 'Genome-wide Z-score <5' & zCategory.T2 == 'Genome-wide Z-score ≥5', 'T1 < 5, T2 > 5', convStatus),
        convStatus  = ifelse(zCategory.T1 == 'Genome-wide Z-score <5' & zCategory.T2 == 'Genome-wide Z-score <5', 'T1 < 5, T2 < 5', convStatus),
        convStatus  = ifelse(zCategory.T1 == 'Genome-wide Z-score ≥5' & zCategory.T2 == 'Genome-wide Z-score <5', 'T1 > 5, T2 < 5', convStatus)
    )

fit.Z <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ convStatus, data = data.Z)
names(fit.Z$strata) <-  c('Aneuploidy score <5<br>(Baseline)', 'Aneuploidy ≥5<br>(Baseline)')

# Calculate HR.
fit.Z.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ convStatus, data = data.Z)

plotFits$fit.Z.hr <- plotSurvival(fit.Z, ylim = 22, data = data.Z, hr = fit.Z.hr)


# Author:                      Job van Riet
# Date:                        16-06-2021
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
plotSurvival <- function(fit, ylim, data, palette = 'jco', hr = NULL){

    # Generate survival plot.
    x <- survminer::ggsurvplot(
        fit = fit,
        pval = F,
        size = .8,
        break.time.by = 10,
        break.y.by = .2,
        palette = palette,
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
        monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12)
    ) %>%
    dplyr::inner_join(data.Patient$Overview) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        `Inclusion (Treated with Caba)` = ifelse(is.na(`Inclusion (Treated with Caba)`), 'No', `Inclusion (Treated with Caba)`),
        ctcCategoryOn3 = ifelse(`CTC Count (Baseline – 7.5mL)` >= 3, 'CTC Count (Baseline) ≥3', 'CTC Count (Baseline) <3'),
        ctcCategoryOn5 = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5, 'CTC Count (Baseline) ≥5', 'CTC Count (Baseline) <5'),
        ARV7.PosvsOther = ifelse(`AR-V7 (Baseline)` %in% c('Neg.', 'Und.'), 'Neg. / Und.', `AR-V7 (Baseline)`),
        WHO.Pooled = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1 - 2', `WHO/ECOG PS at registration`)
    ) %>%

    # Combine
    data.frame()

plotFits <- list()


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


# Survival Analysis (Cox regression) --------------------------------------

## Survival - AR-V7 (All included; n = 133) ----

fit.AllInClusion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival)
names(fit.AllInClusion$strata) <-  c('AR-V7 Neg.', 'AR-V7 Pos.', 'AR-V7 Und.')

plotFits$AllInclusion <- plotSurvival(fit.AllInClusion, data = data.Survival, ylim = 51, palette = c('#648FFF', '#FE6100', '#4D4D4D'))

# Hazard Ratio of multiple groups.
fit.AllInClusion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.Survival)
plotFits$hazardAll <- survminer::ggforest(fit.AllInClusion.hr, data = data.Survival, noDigits = 2)


## Survival - Caba-treated (per AR-V7 conversion) ----

data.BetweenConversion <- data.Survival %>% dplyr::filter(Inclusion..Treated.with.Caba. == 'Yes')
fit.BetweenConversion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.BetweenConversion)
names(fit.BetweenConversion$strata) <-  c('AR-V7 Conversion (Pos. -> Neg.)', 'AR-V7 Retainment (Pos. -> Pos.)')

# Calculate HR.
fit.BetweenConversion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.BetweenConversion)
plotFits$fit.BetweenConversion <- plotSurvival(fit.BetweenConversion, ylim = 21, data = data.BetweenConversion, hr = fit.BetweenConversion.hr, palette = c('#00A94D', '#FE6100'))


## Survival - CTC (All included; n = 133) ----

fit.CTC3 <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategoryOn3, data = data.Survival)
names(fit.CTC3$strata) <-  c('CTC Count <3<br>(Baseline)', 'CTC Count ≥3<br>(Baseline)')

# Calculate HR.
fit.CTC3.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategoryOn3, data = data.Survival)
plotFits$fit.CTC3 <- plotSurvival(fit.CTC3, ylim = 45, data = data.Survival, hr = fit.CTC3.hr, palette = c('#f23005', '#ffbe73'))


fit.CTC5 <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategoryOn5, data = data.Survival)
names(fit.CTC3$strata) <-  c('CTC Count <5<br>(Baseline)', 'CTC Count ≥5<br>(Baseline)')

# Calculate HR.
fit.CTC5.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ctcCategoryOn5, data = data.Survival)
plotFits$fit.CTC5 <- plotSurvival(fit.CTC5, ylim = 45, data = data.Survival, hr = fit.CTC5.hr, palette = c('#f23005', '#ffbe73'))

## Survival - Genome-wide Z (All included; n = 127) ----

fit.GenomeWideZ <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = data.Survival %>% dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.'))
names(fit.GenomeWideZ$strata) <-  c('Aneuploidy score <5<br>(Baseline)', 'Aneuploidy score ≥5<br>(Baseline)')

# Calculate HR.
fit.GenomeWideZ.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = data.Survival %>% dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.'))
plotFits$fit.GenomeWideZ <- plotSurvival(fit.GenomeWideZ, ylim = 45, data = data.Survival %>% dplyr::filter(Genome.Wide.Z.Score..Baseline. != '.'), hr = fit.GenomeWideZ.hr, palette = c('#2d3e50', '#1bbc9b'))


## Survival - WHO (All included; n = 127) ----

fit.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.Pooled, data = data.Survival)
names(fit.WHO$strata) <-  c('WHO Status: 0', 'WHO Status: 1-2')

# Calculate HR.
fit.WHO.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.Pooled, data = data.Survival)
plotFits$fit.WHO <- plotSurvival(fit.WHO, ylim = 45, data = data.Survival, hr = fit.WHO.hr, palette = c('#ED468B', '#5A86C5'))


# Multivariate Cox-regression ---------------------------------------------

## Determine relevant factors (p <= 0.1; n = 127 complete cases) ----

data.AIC <- data.Survival %>%
    dplyr::select(
        monthsFromPreScreeningToEnd,
        Survival,
        ARV7.PosvsOther,
        Max..VAF,
        Total.Gleason,
        Age.at.registration,
        Nr..of.coding.mutations,
        Genome.wide.status..Baseline.,
        ctcCategoryOn3,
        WHO.Pooled
    ) %>%
    dplyr::filter(!is.na(Survival), Genome.wide.status..Baseline. != '.', !is.na(Max..VAF))

fit.AIC <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ARV7.PosvsOther + Max..VAF + Total.Gleason + Age.at.registration + Nr..of.coding.mutations + Genome.wide.status..Baseline. + ctcCategoryOn3 + WHO.Pooled, data = data.AIC, ties = 'breslow')

# Determine relevant factors using backward selection.
MASS::stepAIC(fit.AIC, direction = 'backward')


# Multivariate Cox Regression on relevant factors -------------------------

data.MultiCox <- data.Survival %>%
    dplyr::select(
        monthsFromPreScreeningToEnd,
        Survival,
        ARV7.PosvsOther,
        Genome.wide.status..Baseline.,
        ctcCategoryOn3,
        WHO.Pooled
    ) %>%
    dplyr::filter(!is.na(Survival), Genome.wide.status..Baseline. != '.')

fit.MultiCox <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ ARV7.PosvsOther + Genome.wide.status..Baseline. + ctcCategoryOn3 + WHO.Pooled, data = data.MultiCox, ties = 'breslow')
plotFits$fit.MultiCox <- survminer::ggforest(fit.MultiCox, data = data.MultiCox, noDigits = 2)


# Generate Main Figure ----------------------------------------------------

layout <- "
ACD
BCE
"

plotFits$AllInclusion$plot +
    plotFits$AllInclusion$table +
    plotFits$fit.MultiCox +
    plotFits$fit.BetweenConversion$plot +
    plotFits$fit.BetweenConversion$table +
    patchwork::plot_layout(design = layout, heights = c(1, .25), widths = c(1, 1.25, 1), ncol = 3, guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))


# Generate Suppl. Fig. ----------------------------------------------------


layout <- "
ABC
DEF"

plotFits$fit.CTC$plot +
    plotFits$fit.GenomeWideZ$plot +
    plotFits$fit.WHO$plot +
    plotFits$fit.CTC$table +
    plotFits$fit.GenomeWideZ$table +
    plotFits$fit.WHO$table +
    patchwork::plot_layout(design = layout, heights = c(1, .25), guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))

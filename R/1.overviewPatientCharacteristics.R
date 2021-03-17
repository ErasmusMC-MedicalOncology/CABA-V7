# Author:                      Job van Riet
# Date:                        03-03-2021
# Function:                    Determine CABA-V7 patient characteristics.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(tableone)
library(patchwork)

# Helper theme.
theme_Job <- theme(
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
    strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
    panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
    panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    strip.background = ggplot2::element_rect(colour = 'black', fill = 'white')
)

# Import data -------------------------------------------------------------

data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$PSA <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'PSA Measurements', col_types = c('text', 'date', 'numeric'))
data.Patient$priorChemo <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Prior Treatment - Chemo', col_types = c('text', 'text', 'date', 'date', 'numeric', 'text'))
data.Patient$priorHormonal <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Prior Treatment - Hormonal', col_types = c('text', 'text', 'date', 'date', 'text'))


# Table 1 - Patient Characteristics ---------------------------------------

# Add AR-V7 status (on baseline).
data.Patient$clinicalData <- data.Patient$clinicalData %>%
    # Add AR-V7 status (on baseline).
    dplyr::left_join(data.Patient$Overview %>% dplyr::distinct(`Subject Number`, `AR-V7 (Baseline)`)) %>%

    # Determine which patients received prior local therapy
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        'Local therapy' = base::ifelse('Yes' %in% c(`Radical prostactomy`, `Bilateral orchiectomy`, `Prior internal RTX`, `Prior external radiotherapy prostate`), 'Yes', 'No')
    ) %>%
    dplyr::ungroup()

# Generate Table 1.
tableOne <- tableone::CreateTableOne(
    vars = c(
        'Age at registration',
        'WHO/ECOG PS at registration',
        'Total Gleason',
        'Local therapy',
        'Prior hormonal therapy',
        'Prior chemotherapy',
        'Prior RTX metastases',
        'Prior other treatment',
        'M-stage primary diagnosis',
        'PSA at primary diagnosis [ug/L]'
    ),
    strata = c("AR-V7 (Baseline)"),
    data = data.Patient$clinicalData,
    smd = F,
    test = F,
    addOverall = T
)

# Print Table 1.
print(
    tableOne,
    test = F,
    nonnormal = c('Age at registration', 'Total Gleason', 'M-stage primary diagnosis', 'PSA at primary diagnosis [ug/L]'),
    exact = c('n'),
    smd = F,
    minMax = F,
    missing = T,
    explain = T,
    showAllLevels = T,
    quote = F,
    noSpaces = F
) %>% write.table(row.names = F, sep = '\t', quote = F, file = 'asd.txt')


# Determine the number of hormonal treatments. ----------------------------

rbind(
    data.Patient$clinicalData %>%
        dplyr::select(`Subject Number`, `Radical prostactomy`, `Bilateral orchiectomy`) %>%
        reshape2::melt(id.vars = 'Subject Number', variable.name = 'Type of prior hormonal treatment', value.name = 'isADT'),

    data.Patient$priorHormonal %>%
        dplyr::select(`Subject Number`, `Type of prior hormonal treatment`) %>%
        dplyr::mutate(isADT = ifelse(grepl('Degralex|Lucrin|Zoladex', `Type of prior hormonal treatment`), 'Yes', 'No'))
) %>%
    dplyr::left_join(data.Patient$Overview %>% dplyr::distinct(`Subject Number`, `AR-V7 (Baseline)`) %>% dplyr::group_by(`AR-V7 (Baseline)`) %>% dplyr::mutate(totalInSub = dplyr::n_distinct(`Subject Number`)) %>% dplyr::ungroup()) %>%
    dplyr::group_by(isADT, `AR-V7 (Baseline)`, totalInSub) %>%
    dplyr::summarise(totalYes = dplyr::n_distinct(`Subject Number`), label = sprintf('%s (%s%%)', totalYes, round(totalYes / totalInSub, 3) * 100)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(id_cols = 'isADT', names_from = 'AR-V7 (Baseline)', values_from = 'label')


# Determine CTC response --------------------------------------------------

# CTC response is positive if CTC-count <5 in follow-up CTC count during caba-treatment.
responseCTC <- data.Patient$Overview %>%
    dplyr::filter(
        `Inclusion (Treated with Caba)` == 'Yes',
        (!is.na(`CTC – Treatment 3`) | !is.na(`CTC – Treatment 4`) | !is.na(`CTC – Treatment 5`))
    ) %>%
    dplyr::select(`Subject Number`, `CTC Count (Baseline – 7.5mL)`, `CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`) %>%
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(

        # Check if any CTC-measurement is below five CTC after caba-treatment.
        ctcResponse = ifelse(min(`CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`, na.rm = T) < 5, 'Formal Response (# CTC after any caba-treatment < 5)', 'No Response'),

        # Check if any CTC-measurement after caba-treatment is below 50% of the baseline CTC count.
        ctcDecline = ifelse(any((na.omit(c(`CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`)) / `CTC Count (Baseline – 7.5mL)`) <= .5), 'Formal Decline (50% drop of baseline CTC (7.5 mL) in any caba-treatment)', 'No Decline')

    ) %>%
    dplyr::ungroup()


# Group patients per CTC and genome-wide Z-score --------------------------

groupsCategory <- data.Patient$Overview %>%
    dplyr::filter(`Genome-Wide Z Score (Baseline)` != '.') %>%
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        `Genome-Wide Z Score (Baseline)` = as.double(`Genome-Wide Z Score (Baseline)`),
        group = 'Other',
        group = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5 & `Genome-Wide Z Score (Baseline)` >= 5, 'Group A (≥5 CTC & ≥5 Genome-wide Z-score)', group),
        group = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5 & `Genome-Wide Z Score (Baseline)` < 5, 'Group B (≥5 CTC & <5 Genome-wide Z-score)', group),
        group = ifelse(`CTC Count (Baseline – 7.5mL)` < 5 & `Genome-Wide Z Score (Baseline)` < 5, 'Group C (<5 CTC & <5 Genome-wide Z-score)', group),
        group = ifelse(`CTC Count (Baseline – 7.5mL)` < 5 & `Genome-Wide Z Score (Baseline)` >= 5, 'Group D (<5 CTC & ≥5 Genome-wide Z-score)', group)
    ) %>%
    dplyr::ungroup()


# Combine data ------------------------------------------------------------

data.Patient$Overview %>%
    dplyr::left_join(responseCTC %>% dplyr::distinct(`Subject Number`, ctcResponse, ctcDecline)) %>%
    dplyr::left_join(groupsCategory %>% dplyr::distinct(`Subject Number`, group)) %>%
    dplyr::select(`Subject Number`, group, ctcResponse, ctcDecline)


# Figure 2 - Primary Endpoint ---------------------------------------------

# Convert PSA measurements.
data.PSA <- data.Patient$PSA %>%

    # Filter on caba-treated patients.
    dplyr::full_join(data.Patient$clinicalData, by = 'Subject Number') %>%
    dplyr::full_join(data.Patient$Overview, by = 'Subject Number') %>%
    dplyr::filter(`Inclusion (Treated with Caba)` == 'Yes') %>%

    # Convert dates.
    dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%

    # Center on start of caba. treatment.
    dplyr::mutate(
        dayPSA.Caba = `Date: PSA measurement` - `Date: Start of cabazitaxel treatment`,
        dayDeath.Caba = `Date: Death` - `Date: Start of cabazitaxel treatment`,

        # Correct for actual end of study.
        `Date: End of cabazitaxel cycles` = `Date: End of cabazitaxel cycles` + 21,

        duringCaba = ifelse(dayPSA.Caba > 0 & `Date: PSA measurement` <= `Date: End of cabazitaxel cycles`, 'Yes', 'No')
    ) %>%

    # Center on PSA prior to caba.
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        PSA.First = ifelse(any(dayPSA.Caba <= 0), tail(`PSA [ug/L]`[dayPSA.Caba <= 0], 1), -9999999999)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        PSA.Centered = (`PSA [ug/L]` / max(PSA.First)) - 1,
        dayEndOfCabaCycles = `Date: End of cabazitaxel cycles` - `Date: Start of cabazitaxel treatment`,
        maxPSA.Decrease = `PSA [ug/L]`[duringCaba == 'Yes'][which.min(abs(`PSA [ug/L]`[duringCaba == 'Yes']))],
        maxDeltaPSA.Decrease = ifelse(all(duringCaba == 'No'), NA, maxPSA.Decrease - max(PSA.First))
    ) %>%
    dplyr::ungroup() %>%

    # Threshold PSA increase.
    dplyr::mutate(PSA.Centered = ifelse(PSA.Centered > 1, 1, PSA.Centered)) %>%

    # Sort patients on total length of caba-treatment.
    dplyr::mutate(`Subject Number` = factor(`Subject Number`, levels = unique(`Subject Number`[order(-dayEndOfCabaCycles)])))


# Remove data prior to CABA-V7 (except the first PSA) and post-CABA-V7
data.PSA <- data.PSA %>%
    dplyr::filter(duringCaba == 'Yes' | (`PSA [ug/L]` == PSA.First & dayPSA.Caba <= 0))

# PSA per caba. treatment (per patient, centered on start of caba.).
plot.PSA <- data.PSA %>%
    ggplot2::ggplot(aes(x = dayPSA.Caba, y = PSA.Centered)) +

    # Initial pre-caba PSA (baseline).
    ggplot2::geom_hline(yintercept = 0, color = '#F6A18C', alpha = .5, lty = 11, size = .25) +

    # Up/low line of CABA-V7 cycles.
    ggplot2::geom_rect(aes(ymin = -1.5, ymax = 1.5, xmin = 0, xmax = dayEndOfCabaCycles), fill = NA, color = 'royalblue2', lty = '11', lwd = .05) +

    # PSA measurements.
    ggplot2::geom_line(size = .4, lty = 11, color = 'grey50') +
    ggplot2::geom_point(size = .75, color = 'black') +

    # Show which PSA-point used as the baseline.
    ggplot2::geom_point(data = data.PSA %>% dplyr::filter(`PSA [ug/L]` == PSA.First, dayPSA.Caba <= 0), size = .9, color = 'royalblue2') +

    # Show which PSA-point used as the max. diff decrease PSA.
    ggplot2::geom_point(data = data.PSA %>% dplyr::filter(`PSA [ug/L]` == maxPSA.Decrease, duringCaba == 'Yes'), size = .9, color = '#E87E0D') +

    # Show days of death after Caba treatment.
    # ggplot2::geom_point(data = data.PSA %>% dplyr::filter(!is.na(`Date: Death`)) %>% dplyr::distinct(`Subject Number`, dayDeath.Caba), aes(x = dayDeath.Caba, y = 0), color = 'black', shape='\u2020', size = 3) +

    # Scales.
    ggplot2::scale_y_continuous(breaks = c(-1, 0, 1), labels = c('-100%', '0%', '≥100%'), expand = c(.1, .1)) +
    ggplot2::scale_color_discrete(guide = F) +
    scale_x_continuous(breaks = seq(-25, 301, 25), expand = c(0,0), limits = c(-30, 310)) +

    # Labs.
    ggplot2::labs(x = 'Number of days on cabazitaxel treatment', y = 'ΔPSA (in %), µg/L<br><span style = "font-size:6pt">(Centered on PSA prior to Cabazitaxel treatment)</span>') +

    # Facet + themes
    facet_grid(`Subject Number` ~ ., scales = 'free_x', switch = 'y') +
    theme_Job +
    theme(
        axis.text.y = element_text(size = 4.5),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.spacing.y = ggplot2::unit(.1, "cm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted', color = '#D9D9D930'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(color = '#D9D9D9')
    )

# Show initial + delta PSA
plot.deltaPSA <- data.PSA %>%
    dplyr::select(`Subject Number`, 'Prior PSA (µg/L)' = PSA.First, 'Max. ΔPSA Decrease (%)' = maxDeltaPSA.Decrease) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        `Max. ΔPSA Decrease (%)` = sprintf('(%s%s%%)', ifelse(`Max. ΔPSA Decrease (%)` >= 0, '↥', '↧'), round((`Max. ΔPSA Decrease (%)` / `Prior PSA (µg/L)`) * 100)),
        `Subject Number` = factor(`Subject Number`, levels = rev(levels(data.PSA$`Subject Number`)))
    ) %>%
    reshape2::melt(id.vars = 'Subject Number') %>%

    ggplot2::ggplot(aes(x = variable, y = `Subject Number`, label = value, fill = value,)) +

    # Add heatmap.
    ggplot2::geom_tile(width = 2, height = 1, colour = 'white', fill = 'white', lwd = .25, na.rm = T) +
    ggplot2::geom_text(size = 2.5, hjust = 0) +

    # Labs
    ggplot2::labs(x = NULL, y = NULL) +

    # Scales.
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +

    theme_Job +
    theme(
        text = element_text(family = 'Helvetica'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank()
    )

# Show CTC / PSA response per patient.
plot.Response <- data.PSA %>%
    dplyr::distinct(`Subject Number`, `AR-V7 Conversion`, `Response PSA`, `Response CTC`, `Response CTC-Decline`, `Genome-wide status (Baseline)`) %>%
    reshape2::melt(id.vars = 'Subject Number') %>%
    dplyr::mutate(
        `Subject Number` = factor(`Subject Number`, levels = rev(levels(data.PSA$`Subject Number`))),
        variable = factor(variable, levels = c('Response PSA', 'Response CTC', 'Response CTC-Decline', 'AR-V7 Conversion', 'Genome-wide status (Baseline)'))
    ) %>%

    ggplot2::ggplot(aes(x = variable, y = `Subject Number`, fill = value)) +

    # Add heatmap.
    ggplot2::geom_tile(width = .75, height = .75, colour = 'grey25', lwd = .25, na.rm = T) +

    # Labs
    ggplot2::labs(x = NULL, y = NULL) +

    # Scales.
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_manual(values = c('Formal Decline (50% drop of baseline CTC (7.5 mL) in any caba-treatment)' = '#659A37', 'Formal Response (# CTC after any caba-treatment < 5)' = '#659A37', 'Formal Response (Two consecutive PSA decreases ≤ 50% from baseline PSA)' = '#659A37', 'Genome-wide Z-score <5' = '#80b3ff', 'Genome-wide Z-score ≥5' = '#ff5555', 'Pos.' = '#1752A0', 'Neg.' = '#F5740F', 'No Decline' = 'white', 'No Response' = 'white', '.' = 'grey75'), na.value = 'grey75', guide = guide_legend(title = 'Response Type', title.position = 'top', title.hjust = .5, nrow = 4)) +

    theme_Job +
    theme(
        text = element_text(family = 'Helvetica'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank()
    )


# Combine plots -----------------------------------------------------------

## Combine landscape tracks.
layout <- "AAAB"

plot.PSA + plot.Response +
    patchwork::plot_layout(design = layout, heights = c(.4, .4, .4, 5, rep(.1, 8)), guides = 'collect') +
    patchwork::plot_annotation(tag_levels = 'a')


# Calculate CI ------------------------------------------------------------

# PSA response.
z <- binom.test(3, n = 25, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('PSA Response (3 of 25): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))

# CTC response.
z <- binom.test(3, n = 20, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('CTC Response (3 of 20): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))

# CTC decline.
z <- binom.test(10, n = 20, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('CTC Decline (10 of 20): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))

# AR-V7 conversion.
z <- binom.test(10, n = 19, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('AR-V7 Conversion (10 of 19): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))


# Comparison to TROPIC ----------------------------------------------------

# Test the PSA-response rates.
test <- data.frame(
    row.names = c('CABA-V7', 'TRPOC'),
    PSA = c(3, 129),
    noPSA = c(22, 200)
)

chisq.test(test, simulate.p.value = T)
fisher.test(test, simulate.p.value = T)


# Suppl. Figure - Swimmers plot -------------------------------------------

# Select patients included into CABA-V7 study.
data.Swimmers <- data.Patient$Overview %>%
    dplyr::full_join(data.Patient$clinicalData, by = 'Subject Number') %>%
    dplyr::filter(`Inclusion (Treated with Caba)` == 'Yes') %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date)

# Order samples.
data.Swimmers <- data.Swimmers %>%
    dplyr::mutate(`Subject Number` = factor(`Subject Number`, levels = data.Swimmers %>% dplyr::arrange(`Date: Start of cabazitaxel treatment`) %>% dplyr::pull(`Subject Number`) %>% rev)) %>%
    dplyr::arrange(`Subject Number`)

# Convert prior therapies for plotting.
data.Swimmers.PriorTherapy <- dplyr::bind_rows(
    data.Patient$priorChemo %>%
        dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
        dplyr::filter(`Subject Number` %in% data.Swimmers$`Subject Number`) %>%
        dplyr::select(`Subject Number`, startDate = `Date: Start prior chemo`, endDate = `Date: End prior chemo`, type = `Type of prior chemo`),
    data.Patient$priorHormonal %>%
        dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
        dplyr::filter(`Subject Number` %in% data.Swimmers$`Subject Number`) %>%
        dplyr::select(`Subject Number`, startDate = `Date: Start prior hormonal treatment`, endDate = `Date: End prior hormonal treatment`, type = `Type of prior hormonal treatment`),
) %>%
    dplyr::mutate(`Subject Number` = factor(`Subject Number`, levels = levels(data.Swimmers$`Subject Number`))) %>%
    dplyr::arrange(`Subject Number`)

# Determine the x-axis breaks.
breaks.date <- data.frame(date = zoo::as.Date(zoo::as.yearmon(seq.Date(from = as.Date('2002-01-01'), to = as.Date('2020-09-01'), by = '3 month'))))

# Plot Gannt/ Swimmers plot.
plot.Swimmer <- data.Swimmers %>%
    ggplot2::ggplot(aes(y = `Subject Number`, yend = `Subject Number`)) +

    # Background shades.
    ggplot2::geom_rect(data = breaks.date %>% dplyr::filter(row_number() %% 2 == 1), ggplot2::aes(xmin = date - 45, xmax = date + 45, ymin = -Inf, ymax = Inf), inherit.aes = F, size = 1, fill = '#D9D9D930', color = '#D9D9D930') +
    ggplot2::geom_segment(ggplot2::aes(x = as.Date('1980-01-01'), xend = as.Date('2080-01-01')), lty = 11, position = ggplot2::position_nudge(y = .5), inherit.aes = T, size = .1, color = 'grey50') +
    ggplot2::geom_segment(ggplot2::aes(x = as.Date('1980-01-01'), xend = as.Date('2080-01-01')), lty = 11, position = ggplot2::position_nudge(y = -.5), inherit.aes = T, size = .1, color = 'grey50') +

    # Line for death of patient.
    ggplot2::geom_linerange(data = data.Swimmers %>% dplyr::filter(!is.na(`Date: Death`)), aes(xmin = `Date: End of cabazitaxel treatment`, xmax = `Date: Death`), size = .25, alpha = .75, color = 'black', lty = '11') +

    # Total cabazitaxel treatment.
    ggplot2::geom_segment(aes(x = `Date: Start of cabazitaxel treatment`, xend = `Date: End of cabazitaxel treatment`), alpha = 0.5, size = 1, color = 'royalblue2') +
    ggplot2::geom_point(aes(x = `Date: Start of cabazitaxel treatment`), size = 1.5, color = 'navy') +
    ggplot2::geom_point(aes(x = `Date: End of cabazitaxel treatment`), size = 1.5, color = 'navy', fill = 'white', shape = 21) +

    # Prior treatments (Start to End).
    ggplot2::geom_linerange(data = data.Swimmers.PriorTherapy %>% dplyr::filter(!is.na(startDate), !is.na(endDate)), aes(xmin = startDate, xmax = endDate, color = type), size = .5, alpha = 0.5, lty = '11', position = ggplot2::position_dodge2(1)) +
    ggplot2::geom_point(data = data.Swimmers.PriorTherapy %>% dplyr::filter(!is.na(startDate), !is.na(endDate)), aes(x = startDate, color = type), size = 1, position = ggplot2::position_dodge2(1)) +
    ggplot2::geom_point(data = data.Swimmers.PriorTherapy %>% dplyr::filter(!is.na(startDate), !is.na(endDate)), aes(x = endDate, color = type), size = 1, fill = 'white', shape = 21, position = ggplot2::position_dodge2(1)) +

    # Prior treatments (Missing end-date).
    ggplot2::geom_point(data = data.Swimmers.PriorTherapy %>% dplyr::filter(is.na(endDate)), aes(x = startDate, color = type), size = 3, shape = '\u21F4', position = ggplot2::position_dodge2(1)) +

    # Death of patient.
    ggplot2::geom_point(data = data.Swimmers %>% dplyr::filter(!is.na(`Date: Death`)), aes(x = `Date: Death`), shape = '\u2020', size = 2) +

    # Color scale.
    ggplot2::scale_color_manual(values = c('Docetaxel' = '#912130', 'Other' = 'red', 'Abiraterone/Zytiga' = '#65A528', 'Bicalutamide' = '#00A2EE', 'Docetaxel' = '#d6c65c', 'Enzalutamide' = 'orange', 'Lucrin/Leuproreline' = '#019581', 'Zoladex/Goserelin' = 'purple', 'Degarelix' = 'black'), guide = guide_legend(title = 'Prior hormonal / chemo-therapy', title.position = 'top', title.hjust = .5, nrow = 2)) +

    # Axis.
    ggplot2::scale_x_date(name = '', breaks = breaks.date$date, date_labels = '%b - %Y', minor_breaks = NULL, expand = c(.02, 1)) +
    ggplot2::labs(x = 'Dates (months)', y = 'Patients included into CABA-V7 study<br><span style = "font-size:6pt">(Ordered by CABA-V7 inclusion data)</span>') +

    # Limit the range of dates without removing underlying date.
    tidyquant::coord_x_date(xlim = c('2014-01-01', '2020-09-01')) +

    # Theme
    theme_Job +

    # Show which patients are still alive
    theme(
        axis.text.y = element_text(color = c(ifelse(data.Swimmers %>% dplyr::arrange(`Subject Number`) %>% dplyr::pull(Survival) == 0, '#009500', 'grey5')), family = 'Helvetica'),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

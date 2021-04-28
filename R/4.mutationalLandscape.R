# Author:                      Job van Riet
# Date:                        03-03-2021
# Function:                    Generate an oncoplot of the (coding) QIASeq mutations.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(patchwork)

# Helper theme.
source('R/misc_themes.R')

# Import data -------------------------------------------------------------

# Patient characteristics.
data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))

# (Filtered) QiaSeq mutations.
dataMuts <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'QIASeq - Filtered')
dataMuts <- dataMuts %>% dplyr::filter(`L-code` %in% data.Patient$Overview$`L-code`)


# Clean-up data -----------------------------------------------------------

## Clean-up mutational consequences ----
dataMuts <- dataMuts %>% dplyr::mutate(
    SYMBOL = `[18]SYMBOL`,
    mutType = gsub('\\|.*', '', `[16]ann_annotation`),
    mutType = gsub('&.*', '', mutType),
    mutType = ifelse(grepl('splice', mutType), 'Splicing variant', gsub('_', ' ', Hmisc::capitalize(mutType))))


## Combine multiple mutations per sample. ----
dataOncoplot <- dataMuts %>%
    dplyr::group_by(`L-code`, SYMBOL, mutType) %>%
    dplyr::summarise(
        mutType = ifelse(n() > 1, 'Multiple mutations', mutType),
        isMutant = T
    ) %>%
    dplyr::ungroup()

## Complete the data. ----
dataOncoplot <- dataOncoplot %>%
    dplyr::mutate(`L-code` = factor(`L-code`, levels = unique(data.Patient$Overview %>% dplyr::filter(`L-code` != 'N/A') %>% dplyr::pull(`L-code`)))) %>%
    tidyr::complete(`L-code`, SYMBOL) %>%
    dplyr::mutate(
        isMutant = ifelse(is.na(isMutant), F, isMutant),
        mutType = ifelse(is.na(mutType), 'Neutral', mutType)
    )


## Add patient metadata ----
dataOncoplot <- dataOncoplot %>%
    dplyr::left_join(data.Patient$Overview %>% dplyr::select(`Subject Number`, `L-code`, `cfDNA yield (ng)`, `Response CTC`, `Response CTC-Decline`, `CTC Count (Baseline – 7.5mL)`, `AR-V7 (Baseline)`, `AR-V7 Conversion`, `Genome-wide status (Baseline)`, `Genome-wide status (T2)`, `Inclusion (Treated with Caba)`), by = 'L-code') %>%
    dplyr::left_join(data.Patient$clinicalData %>% dplyr::select(`Subject Number`, `Response PSA`), by = 'Subject Number')


## Determine total number of (coding) mutations and max. VAF ----
dataOncoplot <- dataOncoplot %>%
    dplyr::left_join(
        dataMuts %>%
            dplyr::group_by(`L-code`) %>%
            dplyr::summarise(
                totalMuts = dplyr::n_distinct(`[4]ID`),
                maxVAF = base::max(`[34]consensus-5:AF`, na.rm = T)
            )
    )

## Sort on mutually exclusiveness. ----
memoData <- reshape2::dcast(dataOncoplot, SYMBOL~`L-code`, value.var = 'isMutant', fun.aggregate = sum)
rownames(memoData) <- memoData$SYMBOL; memoData$SYMBOL <- NULL
memoData[is.na(memoData)] <- 0
memoData <- R2CCBC::memoSort(memoData)

dataOncoplot$SYMBOL <- factor(dataOncoplot$SYMBOL, levels = rev(rownames(memoData)))
dataOncoplot$`L-code` <- factor(dataOncoplot$`L-code`, levels = colnames(memoData))


# Calc. Mut. Excl. Genes. -------------------------------------------------

dataFisher <- dataOncoplot %>%
    dplyr::group_by(SYMBOL, `AR-V7 (Baseline)`) %>%
    dplyr::summarise(
        totalInGroup = dplyr::n_distinct(`L-code`),
        totalMut = sum(isMutant),
        totalNoMut = totalInGroup - totalMut
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(`AR-V7 (Baseline)` != 'Und.')

# Perform Fisher's Exact Test
dataFisher.Results <- do.call(rbind, lapply(unique(dataFisher$SYMBOL), function(gene){

    geneData <- dataFisher %>%
        dplyr::filter(SYMBOL == gene) %>%
        dplyr::summarise(
            SYMBOL = unique(SYMBOL),

            Pos.withMut = totalMut[`AR-V7 (Baseline)` == 'Pos.'],
            Pos.withoutMut = totalNoMut[`AR-V7 (Baseline)` == 'Pos.'],

            Neg.withMut = sum(totalMut[`AR-V7 (Baseline)` == 'Neg.']),
            Neg.withoutMut = sum(totalNoMut[`AR-V7 (Baseline)` == 'Neg.'])

        )

    test <- data.frame(
        row.names = c('Pos.', 'Neg.'),
        mut = c(geneData$Pos.withMut, geneData$Neg.withMut),
        noMut = c(geneData$Pos.withoutMut, geneData$Neg.withoutMut)
    )

    geneData$Fisher.p <- fisher.test(test, hybrid = F, alternative = 'two.sided', simulate.p.value = T)$p.value
    geneData$chi.p <- chisq.test(test, simulate.p.value = T)$p.value

    return(geneData)

}))

# Correct for multiple testing.
dataFisher.Results$Fisher.p.adj <- stats::p.adjust(dataFisher.Results$Fisher.p, method = 'BH')
dataFisher.Results$chi.p.adj <- stats::p.adjust(dataFisher.Results$chi.p, method = 'BH')

# Check direction and effect size.
dataFisher.Results <- dataFisher.Results %>% dplyr::mutate(
    effectSize.Pos = round((Pos.withMut / (Pos.withMut + Pos.withoutMut)) * 100, 1),
    effectSize.Neg = round((Neg.withMut / (Neg.withMut + Neg.withoutMut)) * 100, 1),
    effectSize = sprintf('%s%% vs. %s%%', effectSize.Pos, effectSize.Neg),
    Direction = ifelse(effectSize.Pos > effectSize.Neg, 'Enriched', 'Depleted'),
    Direction = ifelse(Fisher.p.adj <= 0.05, Direction, 'N.s.')
)

# Order on oncoplot appearance.
dataFisher.Results <- dataFisher.Results %>% dplyr::mutate(SYMBOL = factor(SYMBOL, levels(dataOncoplot$SYMBOL)))


# Generate figure ---------------------------------------------------------

tracks.oncoplot <- list()

## Comparison AR-V7 Pos. vs. Neg. and Und. per mutation ----
tracks.oncoplot$frequency <- dataOncoplot %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::mutate(totalMut = sum(isMutant)) %>%
    dplyr::ungroup() %>%
    # Calc. rel. perc. of samples per AR-V7 status.
    dplyr::group_by(SYMBOL, totalMut, `AR-V7 (Baseline)`) %>%
    dplyr::summarise(
        totalWithMut = sum(isMutant),
        totalWithMut.Rel = totalWithMut / totalMut,
    ) %>%
    dplyr::mutate(`AR-V7 (Baseline)` = factor(`AR-V7 (Baseline)`, levels = rev(c('Pos.', 'Neg.', 'Und.')))) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(aes(x = SYMBOL, fill = `AR-V7 (Baseline)`, y = totalWithMut)) +
    ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .7) +
    ggplot2::scale_y_continuous(limits = c(0, 45.1), breaks = c(0, 10, 20, 30, 40, 45), expand = c(0,0)) +
    ggplot2::labs(x = NULL, y = '\\# Mutant samples') +
    ggplot2::scale_fill_manual(values = c('Pos.' = '#FE6100', 'Neg.' = '#648FFF', 'Und.' = '#4D4D4D'), na.value = 'white') +
    ggplot2::coord_flip() +
    theme_Job + theme(axis.text.x = ggplot2::element_text(size = 5, angle = 0, hjust = .5), axis.text.y = ggplot2::element_blank())


## Number of coding mutations ----
tracks.oncoplot$TMB <- dataOncoplot %>%
    dplyr::distinct(`L-code`, totalMuts) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = totalMuts)) +
    ggplot2::geom_bar(stat = 'identity', color = 'black', fill = '#F7EAC6', lwd = .33, width = .6) +
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 2, 4, 6, 8, 10), limits = c(0, 10.1)) +
    ggplot2::labs(y = 'Nr. of<br>coding mutation(s)') +
    themeTrack_Job

## Max. VAF ----

tracks.oncoplot$maxVAF <- dataOncoplot %>%
    dplyr::distinct(`L-code`, maxVAF) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = maxVAF)) +
    ggplot2::geom_bar(stat = 'identity', color = 'black', fill = '#95C7D9', lwd = .33, width = .6) +
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, .2, .4, .6, .8, 1), limits = c(0, 1.1)) +
    ggplot2::labs(y = 'Max. VAF of<br> coding mutation(s)') +
    themeTrack_Job

## CTC Count ----

tracks.oncoplot$countCTC <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `CTC Count (Baseline – 7.5mL)`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = `CTC Count (Baseline – 7.5mL)`)) +
    ggplot2::geom_bar(stat = 'identity', color = 'black', fill = '#E6CAE4', lwd = .33, width = .6) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 5, 10, 25, 50, 100, 500, 1000, 2500, 5000), limits = c(0, 5200)) +
    ggplot2::labs(y = 'CTC Count<br>(Baseline – 7.5mL; log<sub>10</sub>)') +
    themeTrack_Job

## cfDNA yield ---

tracks.oncoplot$cfDNA <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `cfDNA yield (ng)`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = `cfDNA yield (ng)`)) +
    ggplot2::geom_bar(stat = 'identity', color = 'black', fill = '#F75050', lwd = .33, width = .6) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 5, 10, 25, 50, 100, 500, 1000, 2500, 5000), limits = c(0, 1001)) +
    ggplot2::labs(y = 'cfDNA yield (ng)<br>(Baseline; log<sub>10</sub>)') +
    themeTrack_Job

## Oncoplot ----
tracks.oncoplot$oncoplot <- dataOncoplot %>% ggplot(aes(x = `L-code`, y = SYMBOL, fill = mutType)) +
    ggplot2::geom_tile(lwd = .2, width = .8, height = .8, color = 'grey95', na.rm = T) +
    # Colors of mutations.
    ggplot2::scale_fill_manual(values = colorMuts, drop = T) +
    ggplot2::labs(x = 'Samples (CABA-V7) processed with target-panel (QIASeq)<br>(cfDNA; <i>n</i> = 131)', y = 'Genes') +
    # Legend settings.
    ggplot2::guides( fill = guide_legend(title = 'Mutational Categories', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
    theme_Job +
    ggplot2::theme(
        axis.text.y = element_text(size = 7, family = 'Helvetica', face = 'bold'),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

## Genome-wide Z-score ----
tracks.oncoplot$genomeWideZT1 <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `Genome-wide status (Baseline)` ) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'Genome-wide status (Baseline)', fill = `Genome-wide status (Baseline)`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('Genome-wide Z-score <5' = 'black', 'Genome-wide Z-score ≥5' = 'grey70')) +
    themeAnno_Job

tracks.oncoplot$genomeWideZT2 <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `Genome-wide status (T2)`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'Genome-wide status (Baseline - T2)', fill = `Genome-wide status (T2)`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('Genome-wide Z-score <5' = 'black', 'Genome-wide Z-score ≥5' = 'grey70')) +
    themeAnno_Job

## Baseline characteristics ----

tracks.oncoplot$baseline <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `AR-V7 (Baseline)`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'AR-V7 (Baseline)', fill = `AR-V7 (Baseline)`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('Pos.' = '#FE6100', 'Neg.' = '#648FFF', 'Und.' = '#4D4D4D'), na.value = 'white') +
    themeAnno_Job

## Conversion status ----

tracks.oncoplot$Conversion <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `AR-V7 Conversion`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'AR-V7 Conversion', fill = `AR-V7 Conversion`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('Neg.' = 'black', 'Pos.' = 'grey70'), na.value = 'white') +
    themeAnno_Job

## Responses ----

tracks.oncoplot$responseCTC <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `Response CTC`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'Response CTC', fill = `Response CTC`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('black', 'grey75'), na.value = 'white') +
    themeAnno_Job

tracks.oncoplot$responseCTCDecline <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `Response CTC-Decline`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'Response CTC-Decline', fill = `Response CTC-Decline`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('black', 'grey75'), na.value = 'white') +
    themeAnno_Job

tracks.oncoplot$responsePSA <- dataOncoplot %>%
    dplyr::distinct(`L-code`, `Response PSA`) %>%
    ggplot2::ggplot(., aes(x = `L-code`, y = 'Response PSA', fill = `Response PSA`)) +
    ggplot2::geom_tile(width = .8, colour = 'grey50', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = c('black', 'grey75'), na.value = 'white') +
    themeAnno_Job


# Combine plots -----------------------------------------------------------

## Combine landscape tracks.
layout <- "
A#
B#
C#
D#
EF
G#
H#
I#
J#
K#
L#
M#
"

tracks.oncoplot$TMB +
    tracks.oncoplot$maxVAF +
    tracks.oncoplot$countCTC +
    tracks.oncoplot$cfDNA +
    tracks.oncoplot$oncoplot + tracks.oncoplot$frequency +
    tracks.oncoplot$baseline +
    tracks.oncoplot$genomeWideZT1 +
    tracks.oncoplot$genomeWideZT2 +
    tracks.oncoplot$Conversion +
    tracks.oncoplot$responsePSA +
    tracks.oncoplot$responseCTC +
    tracks.oncoplot$responseCTCDecline +
    patchwork::plot_layout(design = layout, heights = c(.25, .25, .25, .25, 2, rep(.075, 7)), widths = c(1, .125), guides = 'collect') +
    patchwork::plot_annotation(tag_levels = 'a')
